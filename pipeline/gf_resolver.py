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
        except (GfResolutionError, ValueError):
            if not keep_unresolved:
                raise
            n_unres += 1
            continue
        n_res += 1
        if np.isfinite(gf[j]) and abs(gf[j] - c) > 1e-6:
            n_changed += 1
        gf[j] = c
        is_canon[j] = True
    out[gf_col] = gf
    out['gf_canonical'] = is_canon
    return out, {'n': len(out), 'n_resolved': n_res,
                 'n_unresolved': n_unres, 'n_changed': n_changed}


# ── synth consumer (#1) — branching-preserved rescale ─────────────────────────

EP_CLUSTER_TOL = 0.005  # eV — HFS components share a lower level (identical EP)


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


def apply_to_synth_array(arr: np.ndarray) -> np.ndarray:
    """Overwrite the GES synth linelist's `loggf` (in place) so each physical line's
    TOTAL matches canonical, preserving HFS branching via a single additive shift.
    Raises on any line absent from canonical (no silent fallback)."""
    n = len(arr)
    keys = [species_key(arr['element'][i], arr['ion'][i], arr['molecule'][i])
            for i in range(n)]
    wls = arr['wave_A'].astype(float)
    eps = arr['lower_state_eV'].astype(float)
    gf = arr['loggf'].astype(float)

    for cl in cluster_physical_lines(keys, wls, eps):
        comp_gf = gf[cl]
        cur_total = float(np.log10(np.sum(10.0 ** comp_gf)))
        w = 10.0 ** comp_gf
        centroid = float(np.sum(wls[cl] * w) / w.sum())
        ep_mean = float(np.mean(eps[cl]))
        canon_total = resolve(keys[cl[0]], centroid, ep_mean)
        shift = canon_total - cur_total
        for i in cl:
            arr['loggf'][i] = gf[i] + shift
    return arr
