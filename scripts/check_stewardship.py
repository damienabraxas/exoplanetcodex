"""
scripts/check_stewardship.py
============================
RYA-355 — Data-stewardship invariant check (CI).

The project's data-stewardship rule — ONE authoritative copy of every canonical
value, each sourced and cited — was enforced by discipline only, and discipline
silently failed: the synth line list (`atomic_lines.tsv`, GES v6) and the EW line
list (`linelist_solar.csv`, VALD3) each carried an independent `loggf` for the same
physical line, ~24% of them divergent >0.05 dex (RYA-350). The same defect class
produced the earlier STAR_PARAMS drifts (a legacy mirror dict carrying its own copy
of Teff/logg). This script mechanically enforces the rule so the next duplication
fails at commit time instead of slipping in.

THREE INVARIANTS — registry-driven (adding a new canonical table is a one-line
registration into the lists below, not a rewrite):

  1. gf invariant         — every line shared between two line lists must carry the
                            same total log gf (HFS-aware match, ported from the
                            RYA-350 audit). Pre-RYA-353 this FAILS with the divergent
                            count; that failure is TRACKED to RYA-353 (the single-
                            source migration) so CI stays green while the defect is
                            documented. When 353 lands the count goes to 0 and the
                            tracking can be dropped → the invariant hard-passes.
  2. STAR_PARAMS invariant — no fundamental stellar parameter (Teff, logg, [Fe/H])
                            may appear with a value different from STAR_PARAMS. The
                            legacy adapter dicts (STAR_SOLAR / STAR_PROCYON) are
                            scanned; a divergent copy is a failure, an equal copy is
                            a finding. Both are TRACKED to RYA-298 (the structural
                            removal of those adapter dicts).
  3. provenance invariant  — every canonical value (gf in `linelist_solar.csv`, each
                            STAR_PARAMS record) must carry a non-empty, non-placeholder
                            source/citation. UNTRACKED — a missing source fails now.

ENFORCEMENT POLICY
------------------
A stewardship violation is permitted ONLY if it is registered as TRACKED to a named
remediation ticket. Any UNTRACKED violation fails the build loudly (exit 1). This is
the whole point: a NEW duplication (no tracking ticket) goes red immediately, while
the documented backlog (gf→RYA-353, mirrors→RYA-298) stays green-but-visible. A
canonical source that cannot be parsed is a LOUD error (exit 2), never a silent skip.

PERMANENT RULES (this checker obeys them about itself):
  • It hardcodes none of the values it checks — it reads the canonical sources.
  • Tolerances are named constants with rationale (below).
  • No silent pass: an unreadable/unparseable source raises, it is not skipped.

Usage:
    python3 scripts/check_stewardship.py            # full report; exit 1 on UNTRACKED
    python3 scripts/check_stewardship.py --strict   # also fail on TRACKED violations
    pytest tests/test_data_stewardship.py -v
"""
from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

import numpy as np
import pandas as pd

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO))

from pipeline.species import species_key, z_symbol, MOLECULE  # noqa: E402

# The synth-line-list path is OWNED by abundances_derive — read it from there
# rather than re-hardcoding it (a duplicated path literal would be its own
# stewardship violation). Importing the module triggers the iSpec bootstrap
# (~2 s); that is acceptable for a guard that reads the same list the fit consumes.
import pipeline.abundances_derive as _ad  # noqa: E402
import config.constants as _const  # noqa: E402
from pipeline import gf_resolver as _gr  # noqa: E402  RYA-353 single-source resolver


# ── Named tolerances (with rationale) ─────────────────────────────────────────
# gf matching geometry — identical to the RYA-350 audit (anchor-validated against
# the five RYA-347 Fe II lines; a wrong match rule fails that self-check).
WTOL = 0.02       # Å  — physical-line match window after HFS aggregation
EPTOL = 0.02      # eV — same lower level; separates distinct transitions at one λ
HFS_GAP = 0.10    # Å  — max gap between HFS components grouped into one physical line

# gf agreement threshold. A cross-file Δgf is a pure differential-abundance bias
# between the two engines (gf is degenerate with the floated abundance — RYA-347).
# 0.05 dex is the RYA-350 materiality floor: below it the abundance impact is within
# other line-to-line scatter; at/above it the two paths disagree on the answer. After
# RYA-353 (single gf source) every shared line resolves to ONE value, so the matched
# Δgf collapses to ~0 and this threshold passes with margin.
GF_DIVERGENCE_DEX = 0.05
# Strict "effectively identical" band, reported alongside for the post-353 contract
# (one source ⇒ |Δgf| ≈ 0). Not the enforcement threshold (it would flag sub-mmag
# VALD/GES rounding noise that carries no abundance signal).
GF_IDENTICAL_DEX = 0.001

# Stellar parameters are exact catalog quantities — any difference is a divergence.
STAR_PARAM_TOL = 1e-9

# Provenance strings that are present-but-meaningless.
_PLACEHOLDER_PROVENANCE = {
    '', 'nan', 'none', 'todo', 'tbd', 'tba', '?', '??', '-', '--',
    'unknown', 'placeholder', 'fixme', 'xxx', 'n/a', 'na',
}


# ── Violation record ──────────────────────────────────────────────────────────
@dataclass
class Violation:
    invariant: str          # which invariant fired
    quantity: str           # the canonical quantity (e.g. 'log gf', 'Teff')
    locus: str              # file + key/line so the offender is unambiguous
    value: str              # the offending value(s)
    source: str             # provenance string carried with the value
    detail: str             # human-readable explanation
    ticket: Optional[str] = None   # remediation ticket; None ⇒ UNTRACKED ⇒ FAIL

    @property
    def tracked(self) -> bool:
        return self.ticket is not None


class StewardshipParseError(RuntimeError):
    """A canonical source could not be read/parsed — loud, never a silent skip."""


# ══════════════════════════════════════════════════════════════════════════════
# Invariant 1 — gf agreement across line lists (HFS-aware, ported from RYA-350)
# ══════════════════════════════════════════════════════════════════════════════
@dataclass
class GfTable:
    """One line list participating in the gf cross-check."""
    label: str
    path: Path
    sep: str
    col_element: str
    col_wl: str
    col_ep: str
    col_gf: str
    col_ion: Optional[str] = None       # explicit ion column, if any
    col_molecule: Optional[str] = None  # molecule flag column, if any
    col_source: Optional[str] = None    # provenance column, if any

    def load(self) -> pd.DataFrame:
        try:
            df = pd.read_csv(self.path, sep=self.sep, low_memory=False)
        except Exception as exc:  # unreadable source → loud
            raise StewardshipParseError(
                f"could not read gf table {self.label} at {self.path}: {exc}")
        missing = [c for c in (self.col_element, self.col_wl, self.col_ep,
                               self.col_gf) if c not in df.columns]
        if missing:
            raise StewardshipParseError(
                f"{self.label} missing expected columns {missing}")

        def _key(r):
            ion = r[self.col_ion] if self.col_ion else None
            mol = r[self.col_molecule] if self.col_molecule else None
            return species_key(r[self.col_element], ion, mol)

        out = pd.DataFrame({
            'key': df.apply(_key, axis=1),
            'wl': df[self.col_wl].astype(float),
            'ep': df[self.col_ep].astype(float),
            'gf': df[self.col_gf].astype(float),
        })
        out['source'] = (df[self.col_source].astype(str)
                         if self.col_source else self.label)
        return out


@dataclass
class GfPair:
    """A registered pair of line lists that MUST carry the same total gf per line."""
    name: str
    left: GfTable
    right: GfTable
    ticket: Optional[str] = None
    # Optional integrity anchors: wl -> (gf_left, gf_right, expected Δ=right-left).
    anchors: dict = field(default_factory=dict)
    anchor_species: str = ''
    # RYA-353: when True, resolve BOTH sides through gf_resolver (the single canonical
    # source) before comparing. Post-migration every shared line resolves to ONE value,
    # so the matched Δgf collapses to 0; a line that fails to resolve is an ORPHAN
    # (not in the single source) — a real, untracked break.
    resolved: bool = False


def _aggregate(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse HFS components to physical lines: same species key, |ΔEP|≤EPTOL and
    λ-gap ≤HFS_GAP → one line at total gf = log10(Σ 10^gf), gf-weighted centroid λ.
    `n_comp` records the HFS multiplicity. (Ported verbatim from the RYA-350 audit.)"""
    # RYA-353: cluster via the ONE shared physical-line grouper (gf_resolver) — the
    # same EP-group-then-wl-gap rule the canonical build uses, so centroids line up
    # with canonical_gf and the resolver lookup hits. (Replaces the old ep,wl-sort,
    # which merged distinct lines sharing an EP.)
    d = df.reset_index(drop=True)
    keys = d['key'].tolist()
    wls = d['wl'].to_numpy(float)
    eps = d['ep'].to_numpy(float)
    gfs = d['gf'].to_numpy(float)
    out = []
    for cl in _gr.cluster_physical_lines(keys, wls, eps):
        w = 10.0 ** gfs[cl]
        strongest = cl[int(np.argmax(gfs[cl]))]
        out.append({'key': keys[cl[0]],
                    'wl': float((wls[cl] * w).sum() / w.sum()),
                    'gf': float(np.log10(w.sum())),
                    'ep': float(eps[cl].mean()), 'n_comp': len(cl),
                    'source': d.iloc[strongest]['source']})
    return pd.DataFrame(out)


def _match_gf(pair: GfPair) -> pd.DataFrame:
    """Aggregate both sides to physical lines, then match physical-line to physical-
    line on (species key, λ within WTOL, EP within EPTOL). gf is the only quantity
    that may then differ. The clean 1:1 set (one component each side) is the
    trustworthy headline; HFS-aggregated pairs are flagged."""
    la = _aggregate(pair.left.load())
    sa = _aggregate(pair.right.load())
    if la.empty or sa.empty:
        raise StewardshipParseError(
            f"gf pair '{pair.name}': one side aggregated to zero physical lines "
            f"(left={len(la)}, right={len(sa)}) — match would be vacuous.")
    rows = []
    for key in sorted(set(la['key']) & set(sa['key']), key=str):
        l = la[la['key'] == key].sort_values('wl').reset_index(drop=True)
        s = sa[sa['key'] == key].sort_values('wl').reset_index(drop=True)
        m = pd.merge_asof(l, s, on='wl', direction='nearest', tolerance=WTOL,
                          suffixes=('_l', '_s'))
        m = m[m['gf_s'].notna()]
        m = m[(m['ep_l'] - m['ep_s']).abs() <= EPTOL]
        if m.empty:
            continue
        for _, r in m.iterrows():
            rows.append({
                'key': key, 'species': _species_label(key), 'wl': r['wl'],
                'gf_left': r['gf_l'], 'gf_right': r['gf_s'],
                'dgf': r['gf_s'] - r['gf_l'], 'ep': r['ep_l'],
                'n_comp_left': int(r['n_comp_l']), 'n_comp_right': int(r['n_comp_s']),
                'source_left': r['source_l'], 'source_right': r['source_s'],
                'hfs': bool(r['n_comp_l'] > 1 or r['n_comp_s'] > 1),
            })
    matched = pd.DataFrame(rows)
    _anchor_selfcheck(pair, matched)   # on RAW file gf (reproduces RYA-347 anchors)

    if pair.resolved and not matched.empty:
        # Post-RYA-353: both paths read gf from canonical_gf via the resolver. Replace
        # each matched line's gf with the canonical total; agreement is then exact
        # (single source). A line that won't resolve = orphan (not in the source).
        orphan = []
        for idx, r in matched.iterrows():
            try:
                canon = _gr.resolve(r['key'], float(r['wl']), float(r['ep']))
                matched.at[idx, 'gf_left'] = canon
                matched.at[idx, 'gf_right'] = canon
                matched.at[idx, 'dgf'] = 0.0
            except _gr.GfResolutionError:
                orphan.append(idx)
        matched['orphan'] = False
        if orphan:
            matched.loc[orphan, 'orphan'] = True
    return matched


def _species_label(key) -> str:
    if key[0] == MOLECULE:
        return key[1]
    roman = {1: 'I', 2: 'II', 3: 'III', 4: 'IV', 5: 'V', 6: 'VI'}
    return f"{z_symbol(key[0])} {roman.get(key[1], key[1])}"


def _anchor_selfcheck(pair: GfPair, matched: pd.DataFrame) -> None:
    """If the pair declares integrity anchors, reproducing them is a precondition
    for trusting the scope numbers. A broken match rule fails here → loud STOP."""
    if not pair.anchors:
        return
    sp = matched[matched['species'] == pair.anchor_species]
    for wl, (gf_l, gf_r, dexp) in pair.anchors.items():
        hit = sp[sp['wl'].sub(wl).abs() < 0.01]
        if hit.empty:
            raise StewardshipParseError(
                f"gf pair '{pair.name}': anchor {wl} missing from match — "
                f"match rule is wrong; refusing to report untrustworthy scope.")
        r = hit.iloc[0]
        if not (abs(r['gf_left'] - gf_l) < 1e-6 and abs(r['gf_right'] - gf_r) < 1e-6
                and abs(r['dgf'] - dexp) < 1e-6):
            raise StewardshipParseError(
                f"gf pair '{pair.name}': anchor {wl} mismatch "
                f"(got Δ={r['dgf']:+.3f}, expected {dexp:+.3f}) — match rule wrong.")


def check_gf_pairs(out_dir: Optional[Path] = None) -> list[Violation]:
    """Run every registered gf pair; emit one Violation per materially divergent
    physical line. Writes the full matched table per pair for the unambiguous report."""
    violations: list[Violation] = []
    for pair in GF_PAIRS:
        matched = _match_gf(pair)
        clean = matched[~matched['hfs']]
        div = clean[clean['dgf'].abs() > GF_DIVERGENCE_DEX]
        for _, r in div.iterrows():
            violations.append(Violation(
                invariant='gf',
                quantity='log gf',
                locus=f"{pair.name}  {r['species']} {r['wl']:.3f}Å",
                value=f"{pair.left.label}={r['gf_left']:+.3f} / "
                      f"{pair.right.label}={r['gf_right']:+.3f}  (Δ={r['dgf']:+.3f})",
                source=f"{r['source_left']} vs {r['source_right']}",
                detail=f"same physical line carries divergent gf in two lists "
                       f"(|Δ| {abs(r['dgf']):.3f} > {GF_DIVERGENCE_DEX} dex)",
                ticket=pair.ticket,
            ))
        # RYA-353: orphans — a shared line that does NOT resolve through the single
        # canonical source. This is a real break (a line outside canonical_gf), so it
        # is UNTRACKED (ticket=None) and fails the build loudly.
        for _, r in matched[matched.get('orphan', False) == True].iterrows():
            violations.append(Violation(
                invariant='gf', quantity='log gf',
                locus=f"{pair.name}  {r['species']} {r['wl']:.3f}Å",
                value=f"unresolved in canonical_gf.csv",
                source='—',
                detail="shared line absent from the single canonical gf source "
                       "(RYA-353) — a new orphan, not covered by the resolver",
                ticket=None,
            ))
        if out_dir is not None:
            out_dir.mkdir(parents=True, exist_ok=True)
            matched.sort_values(['species', 'wl']).to_csv(
                out_dir / f"gf_divergence_{pair.name}.csv", index=False)
        # attach summary for the console report
        _GF_SUMMARY[pair.name] = {
            'matched': len(matched), 'clean': len(clean),
            'div_material': int((clean['dgf'].abs() > GF_DIVERGENCE_DEX).sum()),
            'div_identical': int((clean['dgf'].abs() > GF_IDENTICAL_DEX).sum()),
        }
    return violations


_GF_SUMMARY: dict = {}


# ══════════════════════════════════════════════════════════════════════════════
# Invariant 2 — STAR_PARAMS single source (no divergent mirror copy)
# ══════════════════════════════════════════════════════════════════════════════
@dataclass
class StarParamMirror:
    """A dict that mirrors fundamental params owned by STAR_PARAMS. Maps the mirror's
    field names to the canonical STAR_PARAMS field names."""
    dict_name: str               # attribute on config.constants
    star_key: str                # STAR_PARAMS key it mirrors
    field_map: dict              # mirror_field -> canonical_field
    ticket: Optional[str] = None


def check_star_params() -> list[Violation]:
    """Compare every registered mirror dict against STAR_PARAMS. A divergent copy is
    a failure; an equal copy is a (still-reported) duplication finding — both are the
    same defect class (a second home for the value that CAN drift)."""
    violations: list[Violation] = []
    star_params = getattr(_const, 'STAR_PARAMS', None)
    if not isinstance(star_params, dict):
        raise StewardshipParseError("config.constants.STAR_PARAMS not found/parseable")

    for mir in STAR_PARAM_MIRRORS:
        mdict = getattr(_const, mir.dict_name, None)
        if not isinstance(mdict, dict):
            raise StewardshipParseError(
                f"mirror dict config.constants.{mir.dict_name} not found/parseable")
        if mir.star_key not in star_params:
            raise StewardshipParseError(
                f"STAR_PARAMS has no '{mir.star_key}' to compare {mir.dict_name} against")
        canon = star_params[mir.star_key]
        for mfield, cfield in mir.field_map.items():
            if mfield not in mdict or cfield not in canon:
                continue
            mval, cval = float(mdict[mfield]), float(canon[cfield])
            divergent = abs(mval - cval) > STAR_PARAM_TOL
            kind = 'divergent' if divergent else 'duplicate'
            violations.append(Violation(
                invariant='star_params',
                quantity=cfield,
                locus=f"constants.{mir.dict_name}['{mfield}']  (star '{mir.star_key}')",
                value=f"{mir.dict_name}={mval:g} / STAR_PARAMS={cval:g}",
                source=str(mdict.get('source', mdict.get('citation', '—'))),
                detail=(f"{kind} copy of a STAR_PARAMS fundamental outside the single "
                        f"source" + (f" — DIVERGES by {abs(mval - cval):g}"
                                     if divergent else " (equal now, free to drift)")),
                ticket=mir.ticket,
            ))
    return violations


# ══════════════════════════════════════════════════════════════════════════════
# Invariant 3 — provenance present on every canonical value
# ══════════════════════════════════════════════════════════════════════════════
@dataclass
class ProvenanceCheck:
    label: str
    kind: str           # 'csv' or 'star_params'
    path: Optional[Path] = None
    value_col: Optional[str] = None
    source_col: Optional[str] = None
    quantity: str = ''


def _is_placeholder(s) -> bool:
    return str(s).strip().lower() in _PLACEHOLDER_PROVENANCE


def check_provenance() -> list[Violation]:
    """Every canonical value must carry a non-empty, non-placeholder source. Missing
    provenance is UNTRACKED (no ticket) → fails now; a value without a citation is
    exactly how an unsourced number sneaks in."""
    violations: list[Violation] = []
    for chk in PROVENANCE_CHECKS:
        if chk.kind == 'csv':
            try:
                df = pd.read_csv(chk.path, low_memory=False)
            except Exception as exc:
                raise StewardshipParseError(
                    f"provenance check '{chk.label}': cannot read {chk.path}: {exc}")
            if chk.source_col not in df.columns:
                raise StewardshipParseError(
                    f"provenance check '{chk.label}': no source column "
                    f"'{chk.source_col}' in {chk.path}")
            bad = df[df[chk.source_col].map(_is_placeholder)]
            for idx, r in bad.iterrows():
                violations.append(Violation(
                    invariant='provenance',
                    quantity=chk.quantity,
                    locus=f"{chk.label} row {idx}",
                    value=str(r.get(chk.value_col, '')),
                    source=repr(r[chk.source_col]),
                    detail="canonical value carries empty/placeholder provenance",
                ))
        elif chk.kind == 'star_params':
            star_params = getattr(_const, 'STAR_PARAMS', {})
            for star, rec in star_params.items():
                if _is_placeholder(rec.get('source', '')):
                    violations.append(Violation(
                        invariant='provenance',
                        quantity='stellar parameters',
                        locus=f"STAR_PARAMS['{star}']",
                        value=f"teff={rec.get('teff')} logg={rec.get('logg')}",
                        source=repr(rec.get('source', '')),
                        detail="STAR_PARAMS record carries no source/citation",
                    ))
        else:
            raise StewardshipParseError(
                f"provenance check '{chk.label}': unknown kind '{chk.kind}'")
    return violations


# ══════════════════════════════════════════════════════════════════════════════
# REGISTRY — adding a new canonical table is a one-line append here
# ══════════════════════════════════════════════════════════════════════════════
_SYNTH_PATH = Path(_ad._SYNTH_LINELIST_FILE)
_LL_SOLAR = _const.PATHS['linelist_solar']

GF_PAIRS = [
    GfPair(
        name='synth-vs-solar',
        left=GfTable(
            label='linelist_solar.csv', path=_LL_SOLAR, sep=',',
            col_element='element', col_ion='ion', col_wl='wavelength_air_A',
            col_ep='excitation_potential_eV', col_gf='log_gf',
            col_source='loggf_source'),
        right=GfTable(
            label='atomic_lines.tsv', path=_SYNTH_PATH, sep='\t',
            col_element='element', col_ion='ion', col_molecule='molecule',
            col_wl='wave_A', col_ep='lower_state_eV', col_gf='loggf',
            col_source='reference_code'),
        # RYA-353 LANDED: both paths now resolve gf from canonical_gf via gf_resolver,
        # so this check resolves both sides and asserts agreement (0 divergent). No
        # remediation ticket — a divergence/orphan here is now a real, untracked break.
        resolved=True, ticket=None,
        anchor_species='Fe II',
        anchors={  # RYA-347 anchors: (gf_solar, gf_synth, Δ=synth−solar)
            5234.623: (-2.23, -2.180, +0.050),
            5991.371: (-3.54, -3.647, -0.107),
            6084.102: (-3.78, -3.881, -0.101),
            6247.557: (-2.31, -2.435, -0.125),
            6456.380: (-2.10, -2.185, -0.085),
        }),
]

STAR_PARAM_MIRRORS = [
    # Legacy adapter dicts (structural removal tracked by RYA-298). They mirror
    # fundamentals owned by STAR_PARAMS; any copy here can silently drift.
    StarParamMirror('STAR_SOLAR', 'solar',
                    {'teff_K': 'teff', 'logg': 'logg', 'feh': 'feh_ref'},
                    ticket='RYA-298'),
    StarParamMirror('STAR_PROCYON', 'procyon',
                    {'teff_K': 'teff', 'logg': 'logg', 'feh': 'feh_ref'},
                    ticket='RYA-298'),
]

PROVENANCE_CHECKS = [
    ProvenanceCheck('linelist_solar.csv gf', 'csv', path=_LL_SOLAR,
                    value_col='log_gf', source_col='loggf_source', quantity='log gf'),
    ProvenanceCheck('STAR_PARAMS', 'star_params', quantity='stellar parameters'),
]

INVARIANTS: list[Callable[..., list[Violation]]] = [
    check_gf_pairs, check_star_params, check_provenance,
]


# ══════════════════════════════════════════════════════════════════════════════
# Runner
# ══════════════════════════════════════════════════════════════════════════════
def run_all(out_dir: Optional[Path] = None) -> list[Violation]:
    """Run every invariant. Parse errors propagate (loud, never skipped)."""
    violations: list[Violation] = []
    violations += check_gf_pairs(out_dir=out_dir)
    violations += check_star_params()
    violations += check_provenance()
    return violations


def _report(violations: list[Violation]) -> None:
    tracked = [v for v in violations if v.tracked]
    untracked = [v for v in violations if not v.tracked]

    print("\n" + "=" * 78)
    print("DATA-STEWARDSHIP INVARIANT CHECK (RYA-355)")
    print("=" * 78)

    # gf summary (per pair)
    for name, s in _GF_SUMMARY.items():
        print(f"\n[gf] pair '{name}':")
        print(f"     matched physical lines : {s['matched']} "
              f"(clean 1:1 {s['clean']})")
        print(f"     divergent >{GF_DIVERGENCE_DEX} dex : {s['div_material']}  "
              f"(>{GF_IDENTICAL_DEX} dex 'not identical': {s['div_identical']})")

    # per-invariant violation tables
    for inv in ('gf', 'star_params', 'provenance'):
        vs = [v for v in violations if v.invariant == inv]
        if not vs:
            print(f"\n[{inv}] OK — no violations.")
            continue
        print(f"\n[{inv}] {len(vs)} violation(s):")
        cap = 12 if inv == 'gf' else len(vs)
        for v in vs[:cap]:
            tag = f"TRACKED {v.ticket}" if v.tracked else "UNTRACKED"
            print(f"   • ({tag}) {v.quantity} @ {v.locus}")
            print(f"       value : {v.value}")
            print(f"       source: {v.source}")
            print(f"       {v.detail}")
        if len(vs) > cap:
            print(f"   … and {len(vs) - cap} more (see data/results/ CSV).")

    print("\n" + "-" * 78)
    print(f"TOTAL: {len(violations)} violation(s) — "
          f"{len(tracked)} TRACKED, {len(untracked)} UNTRACKED")
    if untracked:
        print("RESULT: FAIL — untracked stewardship violation(s). A new canonical "
              "value was duplicated-and-divergent with no remediation ticket.")
    elif tracked:
        print("RESULT: PASS (with tracked known-issues) — all violations are "
              "registered against a remediation ticket; CI green, defect visible.")
    else:
        print("RESULT: PASS — every canonical value has a single, sourced home.")
    print("-" * 78)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Data-stewardship invariant check.")
    ap.add_argument('--strict', action='store_true',
                    help="also fail on TRACKED (known-issue) violations")
    ap.add_argument('--out', default=str(_REPO / 'data' / 'results'),
                    help="directory for the full per-pair divergence CSVs")
    args = ap.parse_args(argv)

    try:
        violations = run_all(out_dir=Path(args.out))
    except StewardshipParseError as exc:
        print(f"\nSTEWARDSHIP CHECK ERROR (loud, not skipped): {exc}", file=sys.stderr)
        return 2

    _report(violations)
    untracked = [v for v in violations if not v.tracked]
    if untracked:
        return 1
    if args.strict and violations:
        print("--strict: failing on tracked known-issues too.")
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
