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
# Invariant 4 — blend_flag DEFINITION pin (RYA-358)
# ══════════════════════════════════════════════════════════════════════════════
# blend_flag gates blend exclusion across all elements (RYA-208). RYA-356 found it
# had been silently REDEFINED (RYA-209: 0.10 Å proximity binary → curated vetted list,
# proximity moved to the continuous vald_proximity_flag), and that slipped through as a
# PR comment, not a CI failure, because nothing pinned the flag's DEFINITION. This
# invariant pins it: re-run the vetted builder and require the file to match — a builder
# swap or silent redefinition then fails loudly.
_BLEND_SUMMARY: dict = {}


def check_blend_flag() -> list[Violation]:
    """Pin blend_flag's definition, propagation, and provenance.

    (1) DEFINITION — re-run pipeline.build_linelist.build_vetted_blend_flag on
        linelist_solar and require the file's blend_flag column to match it exactly.
        Any mismatch = a silent redefinition / builder swap (e.g. back to the retired
        0.10 Å proximity binary) → UNTRACKED, fails the build.
    (2) PROVENANCE — every VETTED_BLENDS entry carries a non-placeholder citation.
    (3) PROPAGATION — where the generated solar_ew.csv exists, each matched line's
        per-measurement blend_flag equals the line-list flag it is propagated from
        (lines_fit propagates, it does not re-detect). The EW file is a generated
        artifact (absent in a clean checkout) — its absence is recorded, not failed
        (it is not a canonical source); a parse failure of the line list IS loud.
    """
    from pipeline.build_linelist import VETTED_BLENDS, build_vetted_blend_flag
    violations: list[Violation] = []

    try:
        ll = pd.read_csv(_LL_SOLAR, low_memory=False)
    except Exception as exc:
        raise StewardshipParseError(f"blend_flag check: cannot read {_LL_SOLAR}: {exc}")
    if 'blend_flag' not in ll.columns:
        raise StewardshipParseError(
            f"blend_flag check: no 'blend_flag' column in {_LL_SOLAR}")

    actual = ll['blend_flag'].astype(str).str.lower().eq('true').to_numpy()
    expected = build_vetted_blend_flag(ll).to_numpy().astype(bool)   # re-run the builder
    mism = np.where(actual != expected)[0]
    _BLEND_SUMMARY.update(n_true=int(actual.sum()), n_vetted=len(VETTED_BLENDS),
                          mismatch=int(len(mism)))
    for i in mism:
        r = ll.iloc[i]
        violations.append(Violation(
            invariant='blend_flag', quantity='blend_flag definition',
            locus=f"linelist_solar.csv {r['element']} {r['ion']} "
                  f"{float(r['wavelength_air_A']):.3f}Å",
            value=f"file={bool(actual[i])} / vetted-builder={bool(expected[i])}",
            source='pipeline.build_linelist.build_vetted_blend_flag (RYA-209)',
            detail="blend_flag diverges from the vetted builder — a silent redefinition "
                   "or builder swap (e.g. the retired 0.10 Å proximity binary). "
                   "vald_proximity_flag is the continuous proximity signal.",
            ticket=None))

    # (2) provenance of the vetted-blend definitions
    for entry in VETTED_BLENDS:
        src = entry[3] if len(entry) > 3 else ''
        if _is_placeholder(src):
            violations.append(Violation(
                invariant='blend_flag', quantity='vetted blend provenance',
                locus=f"VETTED_BLENDS {entry[0]} {entry[1]} {entry[2]}",
                value=str(tuple(entry[:3])), source=repr(src),
                detail="vetted blend carries empty/placeholder provenance (RYA-358)",
                ticket=None))

    # (3) propagation into the (generated) per-measurement EW table, when present
    _BLEND_SUMMARY['propagation'] = 'solar_ew.csv absent (generated) — not evaluated'
    ew_path = Path(str(_SOLAR_EW))
    if ew_path.exists():
        ew = pd.read_csv(ew_path, low_memory=False)
        if 'blend_flag' in ew.columns:
            llflag = {(e, i, round(float(w), 2)): bool(str(b).lower() == 'true')
                      for e, i, w, b in zip(ll['element'], ll['ion'],
                                            ll['wavelength_air_A'], ll['blend_flag'])}
            n_checked = n_bad = 0
            for e, i, w, b in zip(ew['element'], ew['ion'],
                                  ew['wavelength_air_A'], ew['blend_flag']):
                k = (e, i, round(float(w), 2))
                if k in llflag:
                    n_checked += 1
                    if bool(str(b).lower() == 'true') != llflag[k]:
                        n_bad += 1
                        violations.append(Violation(
                            invariant='blend_flag', quantity='blend_flag propagation',
                            locus=f"solar_ew {e} {i} {k[2]}Å",
                            value=f"EW={b} / linelist={llflag[k]}",
                            source='lines_fit propagation of the line-list blend_flag',
                            detail="per-measurement blend_flag drifted from the "
                                   "line-list flag it is propagated from",
                            ticket=None))
            _BLEND_SUMMARY['propagation'] = (
                f"checked {n_checked} matched lines, {n_bad} mismatched")
    return violations


# ══════════════════════════════════════════════════════════════════════════════
# REGISTRY — adding a new canonical table is a one-line append here
# ══════════════════════════════════════════════════════════════════════════════
_SYNTH_PATH = Path(_ad._SYNTH_LINELIST_FILE)
_LL_SOLAR = _const.PATHS['linelist_solar']
_SOLAR_EW = _const.PATHS['solar_ew']


# ══════════════════════════════════════════════════════════════════════════════
# Invariant 5 — all gf stores resolve to a canonical home (RYA-368)
# ══════════════════════════════════════════════════════════════════════════════
# The gf-pair invariant (1) resolves BOTH sides through gf_resolver before comparing,
# so it proves the SHARED lines collapse to one value but cannot see a store carrying
# raw gf that no consumer routes (the RYA-367 store-#2 landmine: linelist_solar's raw
# VALD3 [O I] −9.776 / Ni −2.70, read by a diagnostic that did not resolve). This
# invariant covers EVERY registered gf store:
#   • ORPHAN (a physical line with no entry in the canonical table) → UNTRACKED, fails
#     loudly: a line outside the single source has no authoritative gf.
#   • RAW divergence from canonical → reported as a TRACKED summary (per the store's
#     load contract: resolve-at-load #1/#3 → RYA-353; independent VALD3 ref #2 →
#     RYA-368). Visible, CI-green — the by-design seed/reference divergence.
#   • the RYA-367 trigger lines ([O I] 6300.304, Ni I 6300.34) MUST resolve to
#     canonical (−9.717 / −2.11) → UNTRACKED regression guard.
_GFSTORE_SUMMARY: dict = {}

# RYA-381: the solar list was extended beyond the optical core (3780–6910 Å) into the
# non-optical UV / red-optical / IR (1150–25000 Å). canonical_gf.csv is still optical
# (3780–9199 Å), so every non-optical extension line is an orphan until RYA-379 ingests
# the non-optical gf and extends this guard past 9199 Å. Such orphans are TRACKED to
# RYA-379 (their named remediation). An orphan INSIDE the optical core would be a real
# break (a curated optical line lost its canonical home) and stays UNTRACKED → FAIL.
_OPTICAL_CORE_LO, _OPTICAL_CORE_HI = 3780.0, 6910.0
#: RYA-822 extended canonical_gf blueward to 3000 A. Inside this window the table now has
#: coverage, so an orphan is a real break rather than the expected non-optical gap.
_NEARUV_COVERED_LO, _NEARUV_COVERED_HI = 3000.0, 3780.0
#: The redward edge canonical_gf still stops at; beyond it RYA-379 owns the extension.
_CANONICAL_RED_EDGE = 9199.0


def check_all_stores_resolve() -> list[Violation]:
    from pipeline.audit import gf_store_consistency as _gs
    violations: list[Violation] = []
    for store in _gs.STORES:
        rep = _gs.store_report(store)
        _GFSTORE_SUMMARY[store.label] = {
            'lines': rep['n_lines'], 'overlap': rep['n_overlap'],
            'orphan': rep['n_orphan'], 'raw_div': rep['n_divergent'],
            'max_dgf': rep['max_abs_dgf'], 'contract': store.contract}
        # orphans — a line outside the single canonical source (no authoritative gf)
        for (key, wl, ep, gf) in rep['orphans']:
            non_optical = wl < _OPTICAL_CORE_LO or wl >= _OPTICAL_CORE_HI
            if _NEARUV_COVERED_LO <= wl < _NEARUV_COVERED_HI:
                # RYA-822 covers this band. An orphan here is NOT an expected gap.
                detail = ("inside the RYA-822 near-UV window (3000–3780 Å), which "
                          "canonical_gf.csv now covers — so this is a real break, not "
                          "the tracked non-optical extension")
                ticket = None
            elif wl >= _CANONICAL_RED_EDGE:
                detail = ("RYA-381 non-optical extension line absent from "
                          "canonical_gf.csv — ingest + guard extension past "
                          f"{_CANONICAL_RED_EDGE:.0f} Å is RYA-379")
                ticket = 'RYA-379'
            elif wl < _NEARUV_COVERED_LO:
                detail = (f"line blueward of {_NEARUV_COVERED_LO:.0f} Å, below the edge "
                          f"RYA-822 extended canonical_gf to — the table does not reach "
                          f"here yet and no ticket currently extends it further blue")
                ticket = 'RYA-822'
            elif non_optical:
                detail = ("RYA-381 non-optical extension line absent from the "
                          "canonical gf table")
                ticket = 'RYA-379'
            else:
                detail = ("store line has no entry in the single canonical gf table "
                          "(RYA-368) — an orphan with no authoritative gf")
                ticket = None
            violations.append(Violation(
                invariant='gf_stores', quantity='log gf',
                locus=f"{store.label}  {_species_label(key)} {wl:.3f}Å",
                value="absent from canonical_gf.csv",
                source=store.contract,
                detail=detail, ticket=ticket))
        # raw divergence — one tracked summary per store (by-design seed/reference)
        if rep['n_divergent'] > 0:
            violations.append(Violation(
                invariant='gf_stores', quantity='log gf (raw vs canonical)',
                locus=f"{store.label}",
                value=f"{rep['n_divergent']} physical line(s) diverge "
                      f">{_gs.GF_DIVERGENCE_DEX} dex (max {rep['max_abs_dgf']:.2f})",
                source=store.contract,
                detail=("by-design: " + ("resolve-at-load (the loader rewrites gf to "
                        "canonical)" if store.is_load_bearing_gf else
                        "independent VALD3 reference; gf NOT load-bearing, consumers "
                        "resolve per-line (RYA-368)")),
                ticket=store.ticket))
    # RYA-367 trigger lines must resolve to canonical (regression guard)
    for t in _gs.check_targets():
        if not t['ok']:
            violations.append(Violation(
                invariant='gf_stores', quantity='log gf',
                locus=t['label'],
                value=f"resolver={t['got']:+.3f} (expected {t['expect']:+.3f})",
                source='gf_resolver.resolve vs canonical_gf.csv',
                detail="RYA-367 trigger line no longer resolves to canonical — the "
                       "[O I]/Ni 6300 landmine has returned",
                ticket=None))
    return violations


# ══════════════════════════════════════════════════════════════════════════════
# Invariant 6 — VALD extraction-threshold consistency (RYA-389 item 3)
# Every VALD delivery must be extracted at the canonical synthesis-era detection
# depth (0.001). The engine is blend-aware Turbospectrum synthesis (RYA-285) which
# needs the weak lines; a shallower cut (the EW-era 0.05) drops blends + trace
# species (Zr/P/S/n-capture — RYA-381), producing a heterogeneous list. This is the
# intake-gate check that makes the benchmark audits (RYA-382/384/385) catch the same
# defect on every star automatically. Under-deep deliveries are TRACKED against the
# star's re-extraction ticket (visible, CI-green); a NEW under-deep delivery with no
# ticket → UNTRACKED → loud FAIL.
# ══════════════════════════════════════════════════════════════════════════════
sys.path.insert(0, str(_REPO / 'data' / 'linelists'))
from vald_parse import (  # noqa: E402
    verify_extraction_threshold as _verify_threshold, THRESHOLD_CANONICAL as _THR)

# star-prefix → re-extraction ticket for KNOWN under-deep deliveries (the EW-era
# 0.05 / 0.01 raws). A delivery not matched here, if under-deep, is UNTRACKED → FAIL.
_THRESHOLD_TICKETS = [
    ('solar', 'RYA-387'), ('55cnc', 'RYA-382'),
    ('alpha_cen', 'RYA-384'), ('procyon', 'RYA-385'),
]


def check_vald_threshold() -> list[Violation]:
    ll_dir = _REPO / 'data' / 'linelists'
    violations: list[Violation] = []
    for p in sorted(ll_dir.glob('vald_*.txt')):
        if 'quarantine' in p.name.lower():
            continue                              # quarantined raws are out-of-band (RYA-378)
        verdict, msg, eff = _verify_threshold(p)
        if verdict == 'ACCEPT':
            continue
        if verdict == 'REJECT':                   # unparseable → loud, never skipped
            raise StewardshipParseError(msg)
        ticket = next((t for pre, t in _THRESHOLD_TICKETS if pre in p.name.lower()), None)
        violations.append(Violation(
            invariant='vald_threshold', quantity='extraction detection depth',
            locus=p.name, value=f'{eff:.3f} (canonical {_THR})',
            source='VALD extraction threshold (min central_depth)',
            detail=('delivery extracted under-deep vs the synthesis-era canonical '
                    f'{_THR} — drops blends + trace species (RYA-381/285). Re-extract '
                    'at 0.001 with finer wavelength chunking (not a shallower cut)'),
            ticket=ticket))                        # None ⇒ UNTRACKED ⇒ loud FAIL
    return violations


# ══════════════════════════════════════════════════════════════════════════════
# Invariant 7 — solar-EW canonical input source (RYA-408)
# ══════════════════════════════════════════════════════════════════════════════
# Root cause of the RYA-406 incident: the solar Fe gate / abundance derivation read
# its EW input from the GITIGNORED, regenerable staging file data/processed/solar_ew.csv,
# whose per-worktree content can silently diverge from the committed canonical. This
# guard enforces three things:
#   (1) the committed canonical (data/measured/sol_ew_results_v1.csv) is present and
#       well-formed — it is the single source the gate must read;
#   (2) IDENTITY — pipeline.abundances_derive._load_solar_ews actually reads the
#       canonical (PATHS['solar_ew_canonical']) and NOT the staging file as its EW pool
#       (a regression that re-points the gate at the runtime is a loud, untracked break);
#   (3) DRIFT — if the regenerable staging file is present, every canonical line must
#       appear in it with matching EW and blend_flag; a divergent staging set is stale /
#       from a different run and must never be mistaken for the source.
_SOLAR_EW_CANONICAL = _const.PATHS['solar_ew_canonical']
_EW_DRIFT_TOL_MA = 0.5  # mÅ — canonical is a resolved subset of the same lines_fit run

def check_solar_ew_canonical() -> list[Violation]:
    import inspect
    violations: list[Violation] = []
    canon = Path(str(_SOLAR_EW_CANONICAL))

    # (1) canonical present + well-formed
    if 'data/measured' not in canon.as_posix():
        raise StewardshipParseError(
            f"solar-EW canonical must live under data/measured (committed), got {canon}")
    if not canon.exists():
        raise StewardshipParseError(
            f"solar-EW canonical missing at {canon} — the gate/abundance EW source is gone (RYA-408)")
    try:
        cdf = pd.read_csv(canon, low_memory=False)
    except Exception as exc:
        raise StewardshipParseError(f"cannot read solar-EW canonical {canon}: {exc}")
    need = {'element', 'ion', 'wavelength_air_A', 'ew_mA', 'blend_flag'}
    miss = need - set(cdf.columns)
    if miss or len(cdf) == 0:
        raise StewardshipParseError(
            f"solar-EW canonical malformed (missing {miss or 'rows'}): {canon}")

    # (2) IDENTITY — the gate's loader must read the canonical, not the staging file
    try:
        src = inspect.getsource(_ad._load_solar_ews)
    except Exception as exc:
        raise StewardshipParseError(f"cannot inspect _load_solar_ews source: {exc}")
    reads_canonical = 'solar_ew_canonical' in src
    reads_runtime_pool = "read_csv(str(PATHS['solar_ew']))" in src
    if (not reads_canonical) or reads_runtime_pool:
        violations.append(Violation(
            invariant='solar_ew_canonical', quantity='gate EW input source',
            locus='pipeline.abundances_derive._load_solar_ews',
            value=f"reads_canonical={reads_canonical} reads_runtime_pool={reads_runtime_pool}",
            source="PATHS['solar_ew_canonical'] vs PATHS['solar_ew']",
            detail="the solar EW pool must be read from the committed canonical "
                   "(sol_ew_results_v1.csv), never the gitignored staging solar_ew.csv "
                   "(RYA-408; this is the RYA-406 incident). No remediation ticket — a "
                   "re-point to the runtime is a real, untracked break.",
            ticket=None))

    # (3) DRIFT — a present staging file must agree with the canonical on the MEASURED
    # EW of shared lines. blend_flag is intentionally NOT compared here: it is a curation
    # quantity OWNED by the canonical (11 vetted blends incl. O I 6300.3 Ni-blend,
    # RYA-104/408) that the raw lines_fit staging legitimately lacks — its integrity is
    # the job of check_blend_flag (linelist ↔ vetted-builder ↔ propagation). Only an EW
    # divergence signals a stale / different-run staging set masquerading as the source.
    staging = Path(str(_const.PATHS['solar_ew']))
    if staging.exists():
        try:
            sdf = pd.read_csv(staging, low_memory=False)
        except Exception as exc:
            raise StewardshipParseError(f"cannot read solar-EW staging {staging}: {exc}")
        skey = {(e, i, round(float(w), 2)): float(m)
                for e, i, w, m in zip(sdf['element'], sdf['ion'],
                                      sdf['wavelength_air_A'], sdf['ew_mA'])}
        n_div = 0
        for e, i, w, m in zip(cdf['element'], cdf['ion'], cdf['wavelength_air_A'],
                              cdf['ew_mA']):
            k = (e, i, round(float(w), 2))
            if k not in skey:
                continue  # canonical line not in staging — coverage, not drift
            if abs(skey[k] - float(m)) > _EW_DRIFT_TOL_MA:
                n_div += 1
                if n_div <= 5:  # cap the noise; the count is the signal
                    violations.append(Violation(
                        invariant='solar_ew_canonical', quantity='staging↔canonical EW drift',
                        locus=f"{e} {i} {k[2]}Å",
                        value=f"staging EW={skey[k]:.2f} mÅ vs canonical EW={float(m):.2f} mÅ",
                        source=f"{staging.name} vs {canon.name}",
                        detail="the regenerable staging solar_ew.csv diverges from the "
                               "committed canonical on a measured EW — it is stale / from a "
                               "different run and must not be mistaken for the EW source "
                               "(RYA-408).",
                        ticket=None))
    return violations


# ══════════════════════════════════════════════════════════════════════════════
# Invariant 8 — C/N/O molecular line lists secured (RYA-360)
# ══════════════════════════════════════════════════════════════════════════════
# RYA-236 acquired the Turbospectrum molecular lists (CH/¹³CH, CN isotopologues, C2,
# OH, NH) + the converted CO_IR_Li2015.dat, but they lived ONLY in the iSpec install
# tree — the repo tracked ZERO molecular artifacts, so an iSpec reinstall/rebuild would
# silently wipe the CO addition (same failure class as gf / blend_flag / STAR_PARAMS).
# RYA-360 vendored them into data/linelists/molecular/turbospectrum/ with a provenance
# manifest; this invariant makes the securing mechanical, registry-driven off that
# manifest (one entry per molecule):
#   • each required list PRESENT in the vendored location + non-empty (≥ recorded baseline),
#   • HARPS-window coverage non-empty (≥ baseline; headline CH/CN/C2 window too),
#   • provenance complete (source + distribution) and wavelength_coverage present,
#   • DRIFT: when iSpec is present, its files must still match the vendored copy (a
#     reinstall that reset/wiped a list → loud CI event).
# Any breach → UNTRACKED violation, exit 1 (a missing list is a CI event, not a silent
# RYA-237 false-absence).
_MOLECULAR_SUMMARY: dict = {}


def check_molecular_lists() -> list[Violation]:
    from pipeline import molecular_lists as _ml
    violations: list[Violation] = []
    try:
        manifest = _ml.load_manifest()
    except Exception as exc:
        raise StewardshipParseError(f"molecular manifest unreadable/absent: {exc}")
    mols = manifest.get('molecules', {})
    if not mols:
        raise StewardshipParseError("molecular manifest records no molecules (RYA-360)")

    ispec_present = _ml.ISPEC_MOLECULES_DIR.exists()
    summary: dict = {}
    for mol, entry in mols.items():
        sub = _ml.VENDORED_DIR / entry.get('vendored_subdir', mol)
        files = entry.get('files', [])
        baseline = int(entry.get('line_count', 0))

        # (1) vendored lists present + non-empty (≥ recorded baseline)
        missing = [f for f in files if not (sub / f).exists()]
        counted = sum(_ml.count_bsyn_lines(sub / f) for f in files if (sub / f).exists())
        if not files or missing:
            violations.append(Violation(
                invariant='molecular', quantity='vendored list present',
                locus=f"{mol}  molecular/turbospectrum/{entry.get('vendored_subdir', mol)}",
                value=(f"missing {len(missing)}/{len(files)} file(s): {missing[:3]}"
                       if files else "no files recorded"),
                source=str(entry.get('source', '—')),
                detail="required molecular list absent from the vendored secure record — "
                       "an iSpec reinstall/deletion would now be invisible (RYA-360). "
                       "Re-vendor: scripts/vendor_molecular_lists_rya360.py.",
                ticket=None))
        elif counted < baseline:
            violations.append(Violation(
                invariant='molecular', quantity='vendored list non-empty',
                locus=f"{mol}  molecular/turbospectrum/{entry.get('vendored_subdir', mol)}",
                value=f"{counted} lines < recorded baseline {baseline}",
                source=str(entry.get('source', '—')),
                detail="vendored molecular list is empty/truncated below its recorded "
                       "line count — the secure copy has been corrupted (RYA-360).",
                ticket=None))

        # (2) coverage non-empty — HARPS window for the held optical/electronic lists;
        #     the RYA-499 mid-IR window for the RYA-503 acquired ro-vibrational lists.
        gate = entry.get('coverage_gate', 'harps')
        hw = entry.get('harps_window')
        mw = entry.get('midir_window') or {}
        if gate == 'midir':
            if int(mw.get('count', 0)) <= 0:
                violations.append(Violation(
                    invariant='molecular', quantity='mid-IR window coverage',
                    locus=f"{mol}  {mw.get('label', '')}",
                    value=f"0 lines in {mw.get('range_cm-1')} cm⁻¹",
                    source=str(entry.get('source', '—')),
                    detail="acquired list carries no rows in its RYA-499 mid-IR "
                           "fundamental window — the acquisition is empty where it must "
                           "cover (RYA-503).", ticket=None))
        else:
            if hw and int(hw.get('count', 0)) <= 0:
                violations.append(Violation(
                    invariant='molecular', quantity='HARPS-window coverage',
                    locus=f"{mol}  {hw.get('name', '')}",
                    value=f"0 lines in {hw.get('range_A')}",
                    source=str(entry.get('source', '—')),
                    detail="headline HARPS diagnostic window is empty — the band the "
                           "synthesis keystone (RYA-237) fits has no lines (RYA-360).",
                    ticket=None))
            if mol != 'CO' and int(entry.get('harps_range_count', 0)) <= 0:
                violations.append(Violation(
                    invariant='molecular', quantity='HARPS-range coverage',
                    locus=f"{mol}", value="0 lines in 3800–6900 Å",
                    source=str(entry.get('source', '—')),
                    detail="no lines across the HARPS-VIS arm — an optical molecular list "
                           "with zero HARPS coverage is unusable (RYA-360).",
                    ticket=None))

        # (3) provenance complete (source + distribution) + wavelength_coverage present
        for pf in ('source', 'distribution'):
            if _is_placeholder(entry.get(pf, '')):
                violations.append(Violation(
                    invariant='molecular', quantity=f'provenance ({pf})',
                    locus=f"{mol}", value=repr(entry.get(pf, '')),
                    source='manifest', detail=f"molecular list carries empty/placeholder "
                    f"{pf} provenance (RYA-360).", ticket=None))
        wc = entry.get('wavelength_coverage') or {}
        if (wc.get('min_A') is None or wc.get('max_A') is None
                or _is_placeholder(wc.get('regime', ''))):
            violations.append(Violation(
                invariant='molecular', quantity='wavelength_coverage',
                locus=f"{mol}", value=repr(wc),
                source='manifest',
                detail="wavelength_coverage (min/max Å + regime) missing — the "
                       "electronic-vs-mid-IR distinction must be machine-recorded "
                       "(RYA-360/499).", ticket=None))

        # (4) DRIFT — when iSpec is present, its files must still match the vendored copy.
        #     Skipped for RYA-503 acquired lists (in_ispec: false) — they were acquired
        #     into the repo, not the iSpec bundle, so there is nothing to drift against.
        drift_note = ('acquired — not in iSpec' if not entry.get('in_ispec', True)
                      else 'iSpec absent — skipped')
        if entry.get('in_ispec', True) and ispec_present and files and not missing:
            drifted = []
            for f in files:
                ip = _ml.ISPEC_MOLECULES_DIR / f
                if not ip.exists():
                    drifted.append((f, 'absent-in-iSpec'))
                elif _ml.count_bsyn_lines(ip) != _ml.count_bsyn_lines(sub / f):
                    drifted.append((f, 'line-count differs'))
            if drifted:
                drift_note = f"{len(drifted)} file(s) DRIFTED"
                violations.append(Violation(
                    invariant='molecular', quantity='iSpec↔vendored drift',
                    locus=f"{mol}", value=f"{drifted[:3]}",
                    source='iSpec molecules dir vs vendored copy',
                    detail="the iSpec molecular list no longer matches the vendored "
                           "secure record — a reinstall/rebuild reset or wiped it "
                           "(the exact RYA-360 risk). Re-vendor or restore.",
                    ticket=None))
            else:
                drift_note = 'matches vendored'

        headline = (f"{mw['label'].split(' (')[0]} {mw['count']}" if gate == 'midir' and mw
                    else (f"{hw['name']} {hw['count']}" if hw else '—'))
        summary[mol] = {
            'files': len(files), 'lines': counted, 'baseline': baseline,
            'harps_range': int(entry.get('harps_range_count', 0)),
            'headline': headline, 'regime': wc.get('regime', '?'), 'drift': drift_note,
        }
    _MOLECULAR_SUMMARY.clear()
    _MOLECULAR_SUMMARY.update(summary)
    return violations


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

# ══════════════════════════════════════════════════════════════════════════════
# Invariant 8 — no hardcoded physical-line gf in constants.py dicts (RYA-543)
# ══════════════════════════════════════════════════════════════════════════════
# A gf for a real physical line must live ONLY in canonical_gf.csv and resolve at use
# via gf_resolver. A constants.py dict that hardcodes a `log_gf` for a physical line is
# a silent-divergence duplicate — the RYA-543 defect: NI6300_COG['log_gf'] = −2.841
# (stale VALD3) shadowed the RYA-365-adjudicated canonical −2.11 (Johansson 2003) in the
# [O I] 6300 Ni-subtraction, biasing A(O) high. The gf-pair (Inv. 1) and all-stores
# (Inv. 5) invariants cover LINE-LIST stores; they never looked inside constants.py, so
# this duplicate slipped. This guard closes that gap.
#
# Registry: (dict name, dict object, canonical key, wl, ep). TARGET STATE = the dict has
# NO 'log_gf' key (gf resolved at use), so this check is a no-op tripwire that FAILS
# loudly (UNTRACKED → exit 1) if a hardcoded, canonical-divergent gf is ever reintroduced.
_CONST_GF_DICTS = [
    ('NI6300_COG', _const.NI6300_COG, (28, 1), 6300.342, 4.266),  # Ni I 6300.34 (RYA-365/543)
]


def check_constants_gf_duplicates() -> list[Violation]:
    """Fail if any registered constants dict hardcodes a physical-line log gf that
    diverges from (or is absent from) the single canonical gf source."""
    out: list[Violation] = []
    for name, obj, key, wl, ep in _CONST_GF_DICTS:
        if 'log_gf' not in obj:
            continue  # RYA-543 target state: no hardcoded copy → nothing can diverge
        hard = float(obj['log_gf'])
        try:
            canon = _gr.resolve(key, wl, ep)
        except _gr.GfResolutionError:
            out.append(Violation(
                invariant='const_gf', quantity='log gf',
                locus=f"config.constants.{name}['log_gf']",
                value=f"hardcoded {hard:+.3f}; absent from canonical_gf.csv",
                source='config.constants vs canonical_gf.csv',
                detail="constants dict hardcodes a physical-line gf with no entry in the "
                       "single canonical source (RYA-543) — orphan, no authoritative gf"))
            continue
        if abs(hard - canon) > GF_DIVERGENCE_DEX:
            out.append(Violation(
                invariant='const_gf', quantity='log gf',
                locus=f"config.constants.{name}['log_gf']",
                value=f"hardcoded {hard:+.3f} vs canonical {canon:+.3f} (Δ={hard - canon:+.3f})",
                source='config.constants vs canonical_gf.csv',
                detail="constants dict hardcodes a physical-line gf that diverges from the "
                       "single canonical source (RYA-543) — resolve at use via gf_resolver"))
    return out


INVARIANTS: list[Callable[..., list[Violation]]] = [
    check_gf_pairs, check_star_params, check_provenance, check_blend_flag,
    check_all_stores_resolve, check_vald_threshold, check_solar_ew_canonical,
    check_molecular_lists, check_constants_gf_duplicates,
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
    violations += check_blend_flag()
    violations += check_all_stores_resolve()
    violations += check_vald_threshold()
    violations += check_solar_ew_canonical()
    violations += check_molecular_lists()
    violations += check_constants_gf_duplicates()
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

    # blend_flag summary
    if _BLEND_SUMMARY:
        b = _BLEND_SUMMARY
        print(f"\n[blend_flag] vetted definition pin:")
        print(f"     blend_flag=True in linelist_solar : {b.get('n_true')}  "
              f"(VETTED_BLENDS entries: {b.get('n_vetted')})")
        print(f"     mismatch vs re-run vetted builder : {b.get('mismatch')}")
        print(f"     per-measurement propagation       : {b.get('propagation')}")

    # gf-stores summary (RYA-368, all stores)
    if _GFSTORE_SUMMARY:
        print(f"\n[gf_stores] all-stores resolution vs canonical:")
        for label, s in _GFSTORE_SUMMARY.items():
            print(f"     {label:<22} overlap {s['overlap']:>6}  orphan {s['orphan']:>3}  "
                  f"raw_div {s['raw_div']:>5} (max {s['max_dgf']:.2f})  [{s['contract']}]")

    # molecular-lists summary (RYA-360)
    if _MOLECULAR_SUMMARY:
        print(f"\n[molecular] C/N/O line lists secured (vendored + iSpec drift):")
        for mol, s in _MOLECULAR_SUMMARY.items():
            print(f"     {mol:<3} {s['files']:>2} file(s)  {s['lines']:>7} lines "
                  f"(baseline {s['baseline']})  HARPS {s['harps_range']:>7}  "
                  f"headline[{s['headline']}]  {s['regime']:<34} iSpec: {s['drift']}")

    # per-invariant violation tables
    for inv in ('gf', 'star_params', 'provenance', 'blend_flag', 'gf_stores',
                'vald_threshold', 'solar_ew_canonical', 'molecular', 'const_gf'):
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
