"""
pipeline/gf_resolver.py
=======================
RYA-353 — resolve gf from the single canonical source at load time.

The canonical per-line gf table (`data/linelists/canonical_gf.csv`, built by
`scripts/migrate_gf_single_source.py`) is the one authoritative value per physical
line. Each reader that used to carry its own gf now overwrites it from here
**after read** (the two GES/VALD files live in the shared, read-only `ispec/` tree,
so we cannot edit them in place):

    #1 synth   atomic_lines.tsv      → apply_to_synth_array  (branching-preserved rescale)
    #3 EW      iSpec line regions    → apply_to_regions      (per-line overwrite)
    #2 aux     linelist_solar.csv df → apply_to_linelist_df  (per-line overwrite)

Granularity (RYA-350 sign-off Q2): canonical stores the per-line TOTAL gf. The EW /
region / linelist readers consume that total directly. The synth reader needs the
individual HFS components for Turbospectrum profile shape, so we **rescale** each
line's components by a single additive shift Δ = (canonical_total − current_total):
this forces Σ10^gf to the canonical total while leaving every component-to-component
ratio (the branching) untouched.

NO silent fallback (project rule): a line that does not resolve against the canonical
table RAISES (`GfResolutionError`). A 0-match is never defaulted.
"""
from __future__ import annotations

import functools
from pathlib import Path

import numpy as np
import pandas as pd

from pipeline.species import species_key, MOLECULE
# RYA-765/318: intake debug tracer — no-op unless a trace is active, never changes a
# returned value. Safe to leave here permanently (intake_debug invariants #1/#2).
from pipeline import intake_debug as dbg

# Match tolerances — must equal the canonical build (audit_gf_duplication).
WTOL = 0.02   # Å
EPTOL = 0.02  # eV
HFS_GAP = 0.10  # Å — components within one physical line (mirrors _aggregate)

_CANON = Path(__file__).resolve().parents[1] / 'data' / 'linelists' / 'canonical_gf.csv'


class GfResolutionError(RuntimeError):
    """A line could not be resolved against the canonical gf table."""


def _row_key(key_z, ion) -> tuple:
    """Reconstruct the RYA-345 canonical key from canonical_gf.csv columns."""
    s = str(key_z)
    if ion == '' or (isinstance(ion, float) and np.isnan(ion)):
        return (MOLECULE, s)
    return (int(float(key_z)), int(float(ion)))


@functools.lru_cache(maxsize=1)
def _index() -> dict:
    """key -> (wl_sorted: ndarray, ep: ndarray, loggf: ndarray). Cached per process."""
    if not _CANON.exists():
        raise FileNotFoundError(
            f"canonical gf table not found: {_CANON}\n"
            "Run: python3 scripts/migrate_gf_single_source.py --build")
    df = pd.read_csv(_CANON, low_memory=False)
    idx: dict = {}
    df['_key'] = [_row_key(z, i) for z, i in zip(df['key_z'], df['ion'])]
    for k, g in df.groupby('_key', sort=False):
        g = g.sort_values('wavelength_air_A')
        idx[k] = (g['wavelength_air_A'].to_numpy(float),
                  g['excitation_potential_eV'].to_numpy(float),
                  g['log_gf'].to_numpy(float))
    return idx


def resolve(key: tuple, wl: float, ep: float) -> float:
    """Canonical total log_gf for one physical line. Raises if unresolved."""
    entry = _index().get(key)
    if entry is None:
        raise GfResolutionError(f"species {key} absent from canonical gf table")
    wls, eps, gfs = entry
    near = np.where((np.abs(wls - wl) <= WTOL) & (np.abs(eps - ep) <= EPTOL))[0]
    if near.size == 0:
        raise GfResolutionError(
            f"no canonical gf for {key} at λ={wl:.4f} EP={ep:.4f} "
            f"(±{WTOL}Å/±{EPTOL}eV) — 0-match, not defaulting")
    return float(gfs[near[np.abs(wls[near] - wl).argmin()]])


# ── per-line consumers (#2, #3) ───────────────────────────────────────────────

def apply_to_regions(regions: np.ndarray) -> np.ndarray:
    """Overwrite an iSpec line-regions array's `loggf` from canonical (in place).
    Non-atomic region notes (molecular) resolve via their molecular key too."""
    names = regions.dtype.names or ()
    if 'loggf' not in names or 'lower_state_eV' not in names:
        return regions  # not a gf-bearing region array (nothing to resolve)
    n_set = 0
    for i in range(len(regions)):
        try:
            key = species_key(str(regions['note'][i]))
        except ValueError:
            continue  # un-keyable note; leave as-is (not a gf-bearing abundance line)
        regions['loggf'][i] = resolve(key, float(regions['wave_A'][i]),
                                      float(regions['lower_state_eV'][i]))
        n_set += 1
    return regions


def apply_to_linelist_df(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of a linelist_solar-style df with `log_gf` resolved from
    canonical. Requires element, ion, wavelength_air_A, excitation_potential_eV.
    Raises GfResolutionError on any unresolved line (no silent fallback).

    NOTE: per-row overwrite to the canonical TOTAL — correct only when each df row
    is a distinct physical line. Do NOT use on a store that splits a physical line
    into multiple component rows (e.g. the VALD3 HFS pairs in linelist_solar): each
    component would be set to the full total and the summed strength would double.
    For per-line CONSUMERS (a COG that reads one gf per line) use resolve_df_gf."""
    out = df.copy()
    gf = np.empty(len(out), dtype=float)
    for j, (el, ion, wl, ep) in enumerate(zip(
            out['element'], out['ion'], out['wavelength_air_A'],
            out['excitation_potential_eV'])):
        gf[j] = resolve(species_key(str(el), str(ion)), float(wl), float(ep))
    out['log_gf'] = gf
    return out


def resolve_df_gf(df: pd.DataFrame, *, element_col: str = 'element',
                  ion_col: str = 'ion', wl_col: str = 'wavelength_air_A',
                  ep_col: str = 'excitation_potential_eV', gf_col: str = 'log_gf',
                  keep_unresolved: bool = True) -> tuple[pd.DataFrame, dict]:
    """Resolve a per-line consumer's gf column to canonical, line by line.

    For a CONSUMER that reads one gf per spectral line (e.g. the COG A(Fe)
    diagnostic) — NOT for a store whose physical lines are split into HFS
    component rows. Each row resolves independently via `resolve(key, wl, ep)`:

      * resolved  → gf_col overwritten with the canonical value
      * unresolved (not in canonical): keep_unresolved=True keeps the row's
        existing gf (a CONSUMER may legitimately carry a line outside canonical);
        keep_unresolved=False raises (no silent fallback for an authoritative use).

    Returns (df_copy, stats) with stats = {n, n_resolved, n_unresolved, n_changed}.
    A `gf_canonical` boolean column flags which rows now hold the canonical value.
    """
    out = df.copy().reset_index(drop=True)
    gf = out[gf_col].to_numpy(dtype=float, copy=True)
    is_canon = np.zeros(len(out), dtype=bool)
    n_res = n_unres = n_changed = 0
    for j in range(len(out)):
        el, ion = out.at[j, element_col], out.at[j, ion_col]
        wl, ep = out.at[j, wl_col], out.at[j, ep_col]
        if pd.isna(wl) or pd.isna(ep):
            continue
        try:
            key = species_key(str(el), str(ion))
            c = resolve(key, float(wl), float(ep))
        except (GfResolutionError, ValueError) as exc:
            if not keep_unresolved:
                raise
            # RYA-765/318: `keep_unresolved` is this module's ONE tolerated fallback
            # (the docstring above says why a CONSUMER may legitimately carry a line
            # outside canonical). It is still a fallback: the row keeps a gf from an
            # unnamed source while every neighbour carries the canonical one, and only
            # the aggregate `n_unresolved` count records it. Named per line here.
            dbg.trace_fallback('gf_unresolved_kept',
                               f'{el} {ion} {float(wl):.3f}: {exc} -> keeping the '
                               f"consumer's own log_gf (not canonical)",
                               species=f'{el} {ion}', wavelength=float(wl),
                               excitation_potential_eV=float(ep))
            n_unres += 1
            continue
        n_res += 1
        if np.isfinite(gf[j]) and abs(gf[j] - c) > 1e-6:
            n_changed += 1
        gf[j] = c
        is_canon[j] = True
    out[gf_col] = gf
    out['gf_canonical'] = is_canon
    dbg.trace_check('gf_canonical_fraction', n_unres == 0, severity_if_fail='WARN',
                    detail=f'{n_res}/{len(out)} rows resolved to canonical gf '
                           f'({n_unres} kept their own)',
                    n=len(out), n_resolved=n_res, n_unresolved=n_unres,
                    n_changed=n_changed)
    return out, {'n': len(out), 'n_resolved': n_res,
                 'n_unresolved': n_unres, 'n_changed': n_changed}


def annotate_used_gf(df: pd.DataFrame, *, element_col: str = 'element',
                     ion_col: str = 'ion', wl_col: str = 'wavelength_air_A',
                     ep_col: str = 'excitation_potential_eV', gf_col: str = 'log_gf',
                     intake_source: str = 'linelist intake (VALD)') -> tuple:
    """Make a per-line table report THE gf THAT WAS USED, without losing the intake one.

    RYA-825. The defect this exists to close: a pool or inventory carried a `log_gf`
    copied from the linelist at intake, while the production synthesis path resolved a
    DIFFERENT value out of `canonical_gf.csv` before inverting (`_load_synth_resources`
    defaults `apply_canonical_gf=True`, which runs `apply_to_synth_array`). The column
    therefore described an input nobody used. RYA-799 read that column, concluded the Fe
    pool was "100 % Kurucz K14", and built a whole SCALE-MISMATCH finding on it; RYA-824
    re-inverted and found 18 of 29 lines were already on the canonical value. A column
    that does not describe the computation is not documentation, it is a trap.

    Emits four columns and changes no value that anything computes with:

      `log_gf`            the gf the inversion used — canonical where the resolver
                          reaches, otherwise the intake value unchanged.
      `log_gf_intake`     what the linelist said at intake. Never dropped (RYA-825).
      `gf_source_intake`  where that intake value came from.
      `gf_canonical`      True where the resolver supplied it. This is the column that
                          makes the other three readable: `log_gf == log_gf_intake` means
                          "they agree" only if you also know whether the resolver ran, and
                          outside canonical's 3780-9199.9 A range it does not run at all.

    NOTE ON MEANING. Where `gf_canonical` is True the value becomes the PHYSICAL LINE
    TOTAL that canonical stores, which is not always what a caller's intake column held —
    `scripts/line_accounting_rya709.py`, for instance, aggregated HFS components with
    `max`. That is a semantic change and it is the intended one: the total is what the
    synthesis integrates. The intake convention survives in `log_gf_intake`.
    """
    out, stats = resolve_df_gf(df, element_col=element_col, ion_col=ion_col,
                               wl_col=wl_col, ep_col=ep_col, gf_col=gf_col,
                               keep_unresolved=True)
    out['log_gf_intake'] = df[gf_col].to_numpy(dtype=float, copy=True)
    if 'gf_source_intake' not in out.columns:
        out['gf_source_intake'] = intake_source
    return out, stats


def assert_gf_column_is_honest(df: pd.DataFrame, *, element_col: str = 'element',
                               ion_col: str = 'ion', wl_col: str = 'wavelength_air_A',
                               ep_col: str = 'excitation_potential_eV',
                               gf_col: str = 'log_gf', tol: float = 1e-6) -> dict:
    """The drift guard (RYA-825 spec item 3).

    Re-derives the resolver's answer for every row the table claims is canonical and
    checks the reported column still matches. Raises on the first disagreement, because
    a gf column that has drifted from the resolver is the exact condition that produced
    RYA-799's misdiagnosis, and it is invisible by inspection.
    """
    if 'gf_canonical' not in df.columns:
        raise GfResolutionError(
            f"{gf_col!r} carries no `gf_canonical` flag — cannot tell which rows claim to "
            f"be resolver-sourced, so the column cannot be checked at all (RYA-825)")
    bad, n_checked = [], 0
    for j in range(len(df)):
        if not bool(df.iloc[j]['gf_canonical']):
            continue
        row = df.iloc[j]
        key = species_key(str(row[element_col]), str(row[ion_col]))
        expected = resolve(key, float(row[wl_col]), float(row[ep_col]))
        n_checked += 1
        if abs(float(row[gf_col]) - expected) > tol:
            bad.append((row[element_col], row[ion_col], float(row[wl_col]),
                        float(row[gf_col]), expected))
    if bad:
        head = "; ".join(f"{e} {i} {w:.4f}: reported {g:.4f} vs resolver {x:.4f}"
                         for e, i, w, g, x in bad[:5])
        raise GfResolutionError(
            f"{len(bad)} of {n_checked} rows flagged gf_canonical report a log_gf that "
            f"is NOT the resolver's value — the column has drifted from the gf the "
            f"inversion uses (RYA-825). {head}")

    # THE OTHER HALF, and the half the first version was missing. Checking only the rows
    # that CLAIM to be canonical validates the positive claims and is blind to coverage
    # GROWTH: when canonical_gf gains lines, rows correctly flagged False go stale and
    # nothing notices. That is not hypothetical — RYA-822 extended canonical_gf blueward
    # the same day RYA-825 landed, and 3,413 accounting rows became resolvable while the
    # committed table still reported them as outside coverage. A one-sided guard let a
    # freshly-corrected table go stale within hours.
    n_should = 0
    for j in range(len(df)):
        if bool(df.iloc[j]['gf_canonical']):
            continue
        row = df.iloc[j]
        try:
            resolve(species_key(str(row[element_col]), str(row[ion_col])),
                    float(row[wl_col]), float(row[ep_col]))
        except (GfResolutionError, ValueError):
            continue
        n_should += 1
    if n_should:
        raise GfResolutionError(
            f"{n_should} rows are flagged gf_canonical=False but DO resolve against the "
            f"current canonical gf table — the table has gained coverage since this "
            f"column was built and it is now stale (RYA-825). Regenerate it.")
    return {'n_checked': n_checked, 'n_bad': 0, 'n_newly_resolvable': 0}


# ── synth consumer (#1) — branching-preserved rescale ─────────────────────────

EP_CLUSTER_TOL = 0.005  # eV — HFS components share a lower level (identical EP)


#: Per-isotope totals of one physical line must agree within this for the source to be
#: RYA-684 form (A) — isotope-coded with the full gf on each isotope. Components print
#: 3-4 decimals and a cluster carries up to ~25 of them, so summation rounding alone
#: reaches ~0.01 dex; this is that, doubled. Measured maximum among the real form-(A)
#: clusters in the GES v6 delivery: 0.0097 dex. The one cluster that fails it (Li I 6707)
#: fails by 0.301.
ISOTOPE_FORM_A_SPREAD_DEX = 0.02


def cluster_physical_lines(keys, wls, eps):
    """Group row indices into physical lines — the ONE canonical clustering, used by
    both the build (audit_gf_duplication._aggregate) and the synth reroute so their
    clusters/centroids are identical.

    Two-level and order-robust: per key, group by EP (HFS components share the lower
    level → near-identical EP, within EP_CLUSTER_TOL), then within an EP group split on
    wl-gap > HFS_GAP. EP-first is essential — lines that coincide in wavelength but
    differ in EP (e.g. Sc I 5350.321 at EP 1.969 and 5.413) are DISTINCT transitions
    and must not merge; a wl-sorted scan splits such ties inconsistently."""
    by_key: dict = {}
    for i in range(len(keys)):
        by_key.setdefault(keys[i], []).append(i)
    clusters = []
    for idxs in by_key.values():
        idxs.sort(key=lambda i: (eps[i], wls[i]))
        gstart = 0
        for p in range(1, len(idxs) + 1):
            if p == len(idxs) or eps[idxs[p]] - eps[idxs[gstart]] > EP_CLUSTER_TOL:
                group = sorted(idxs[gstart:p], key=lambda i: wls[i])
                wstart = 0
                for q in range(1, len(group) + 1):
                    if q == len(group) or wls[group[q]] - wls[group[q - 1]] > HFS_GAP:
                        clusters.append(group[wstart:q])
                        wstart = q
                gstart = p
    return clusters


def physical_total(loggf, isotopes=None) -> float:
    """The TOTAL log gf of one physical line, from its components. Isotope-aware.

    🔴 RYA-1075. ``cluster_physical_lines`` is isotope-blind by construction — it groups
    on species + EP + wavelength gap, and two isotopes of the same transition share the
    lower level and sit milli-Angstroms apart, so they land in ONE cluster. The GES v6
    HFS/ISO list is RYA-684 *form (A)*: isotope-coded, and every isotope's component set
    carries the FULL oscillator strength, because the engine applies ``isotopfrac``
    afterwards and Sum(f_i) = 1. So ``log10(sum 10**gf)`` over such a cluster is
    ``n_isotopes x`` the physical gf — a count, with no dependence on the abundances.

    That is how 54 ``canonical_gf`` rows came to be published log10(7) and log10(2) too
    high (Nd II, Sm II, Eu II, Ba II, Cu I). Demonstrated on Eu II 6645: 11 components,
    isotope 151 summing to +0.1199 and isotope 153 to +0.1198 — each the whole
    transition — against a naive all-component sum of +0.4208.

    ⚠️ NOT RYA-684's defect and NOT its correction term. RYA-684 is an ENGINE-side double
    application of ``isotopfrac`` on already-folded gf, offset ``-log10(sum f_i^2)``.
    This is a CONSUMER-side aggregation over a correctly-formed source, offset
    ``log10(n)``. La II separates them: +0.0008 against +0.3010.

    Where the per-isotope sums do NOT agree, the source is not form (A) and the count
    rule is undefined — the naive sum is returned unchanged. Li I 6707 is that case.
    """
    g = np.asarray(loggf, float)
    if isotopes is None:
        return float(np.log10(np.sum(10.0 ** g)))
    a = np.asarray(isotopes)
    present = sorted({int(x) for x in a if int(x) != 0})
    if len(present) < 2:
        return float(np.log10(np.sum(10.0 ** g)))
    per = [float(np.log10(np.sum(10.0 ** g[a.astype(int) == i]))) for i in present]
    if max(per) - min(per) > ISOTOPE_FORM_A_SPREAD_DEX:
        return float(np.log10(np.sum(10.0 ** g)))     # not form (A); nothing to divide out
    return float(np.mean(per))


def apply_to_synth_array(arr: np.ndarray, *, report: dict | None = None) -> np.ndarray:
    """Overwrite the GES synth linelist's `loggf` (in place) so each physical line's
    TOTAL matches canonical, preserving HFS branching via a single additive shift.

    Raises on any line absent from canonical — no SILENT fallback.

    🔴 PASS `report` TO MAKE IT PARTIAL AND LOUD — RYA-1045. COVERAGE IS PER LINE, NOT
    PER BAND. The caller used to decide whether canonical applied by comparing the band's
    wavelength bound against the table's extent, which can only produce FALSE NEGATIVES:
    six Mn I hyperfine components 0.09 A past the red edge switched laboratory gf off for
    4,364 Fe I lines and collapsed IR Fe I to n=4 at rung 1. A hard range is also the
    wrong shape going forward — instruments span different ranges, and a bound tuned to
    one of them is wrong for the next.

    `resolve` is the authoritative test and it is already loud ("0-match, not
    defaulting"). So: attempt every line, and when one genuinely has no canonical entry,
    leave it AS DELIVERED and record it in `report` for the caller's provenance. Without
    `report` the strict raise stands, because a line quietly left on the list's own gf
    inside a product claiming canonical single-sourcing is a lie (RYA-833) — the caller
    must ASK for partial coverage and is then obliged to state it.

    `report` is filled with `applied`, `uncovered` (species/wavelength/reason) and
    `n_uncovered`. It is a caller-supplied dict rather than an attribute on `arr`: numpy
    structured arrays reject attribute assignment, and my first attempt at this lost the
    accounting silently through an `except AttributeError`.
    """
    n = len(arr)
    keys = [species_key(arr['element'][i], arr['ion'][i], arr['molecule'][i])
            for i in range(n)]
    wls = arr['wave_A'].astype(float)
    eps = arr['lower_state_eV'].astype(float)
    gf = arr['loggf'].astype(float)

    uncovered: list[dict] = []
    n_applied = 0
    #  RYA-1075: isotope column when the delivery carries one. Without it `cur_total`
    #  is the isotope-inflated sum, and against a de-inflated canonical total the shift
    #  would scale every component down by n. See `physical_total`.
    isos = (arr['spectrum_synthe_isotope']
            if 'spectrum_synthe_isotope' in (arr.dtype.names or ()) else None)
    for cl in cluster_physical_lines(keys, wls, eps):
        comp_gf = gf[cl]
        cur_total = physical_total(comp_gf, None if isos is None else isos[cl])
        w = 10.0 ** comp_gf
        centroid = float(np.sum(wls[cl] * w) / w.sum())
        ep_mean = float(np.mean(eps[cl]))
        try:
            canon_total = resolve(keys[cl[0]], centroid, ep_mean)
        except GfResolutionError as exc:
            if report is None:
                raise                       # strict: the caller did not ask for partial
            uncovered.append({"species": str(keys[cl[0]]), "wave_A": centroid,
                              "n_components": len(cl), "reason": str(exc)[:160]})
            continue                        # LEFT AS DELIVERED, and named
        shift = canon_total - cur_total
        for i in cl:
            arr['loggf'][i] = gf[i] + shift
        n_applied += 1
    if report is not None:
        report["applied"] = n_applied
        report["uncovered"] = uncovered
        report["n_uncovered"] = len(uncovered)
    return arr
