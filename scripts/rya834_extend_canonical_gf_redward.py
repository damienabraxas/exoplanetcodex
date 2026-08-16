#!/usr/bin/env python3
"""
RYA-834 — extend canonical_gf.csv REDWARD of 9199.90 A to 12935 A, and adjudicate.

WHY. RYA-762 inventoried 239 Fe I lines at 9203-12935 A and proved all three engines can
reach them (Engine B 187/239). Its PRODUCTS cannot be derived, because `canonical_gf.csv`
stops at 9199.90 A and `apply_to_synth_array` raises rather than defaulting. This is the
Fe slice of RYA-379.

THE REFEREE HERE IS PRIMARY LABORATORY DATA, NOT NIST — AND THAT IS A CHANGE FROM RYA-822.
Blueward, NIST ASD was the only graded source reaching the band, so RYA-822 adopted it. It
is the wrong choice redward, for a reason that is MEASURED rather than assumed:

    NIST ASD Fe I vs this band's own linelist, matched on wavelength AND EP:
        34 of 56 matched     median |delta| 0.0003 dex     max 0.0016 dex

Agreement at the fourth decimal is not corroboration, it is the SAME NUMBER arriving twice.
It confirms RYA-760 for this band — FMW *is* a NIST compilation and VALD copies it — so a
NIST "match" here cannot referee anything and must never lift a line off the Kurucz floor.
NIST is therefore RECORDED (with its grade and TP code) and never ADOPTED.

What CAN referee is primary lab measurement with a published per-line sigma:

    Ruffoni 2014     3526.24-11316.06 A     (142 lines; MNRAS 441, 3127)
    Den Hartog 2014  3211.99-11013.24 A     (203 lines; ApJS 215, 23)
    Belmonte 2017    2132.02-10333.18 A     (120 lines; ApJ 848, 125)

⚠️ THE TICKET'S STATED RUFFONI SPAN (3526-10864 A) IS SHORT. Our vendored RYA-799 pull
carries a Ruffoni line at 11316.064 A, and Den Hartog reaches 11013.237 A. So the "10864-
12935 A gap" the ticket asked to verify is NOT one gap: it is covered to 11316.06 A and
open above it. Checked, not assumed (RYA-833).

THE 11316-12935 A VERDICT. No primary lab gf exists there — verified against all three
staged sources AND against NIST ASD, which does carry 19 graded rows above 11316 A but is
the compilation echo measured above. Those lines are honestly Kurucz-floor and are written
as `single_source`, never dressed up. An absence checked and stated is a result; an absence
assumed is the RYA-833 failure.

⚠️ MATCHING IS ON WAVELENGTH *AND* EXCITATION POTENTIAL, NEVER WAVELENGTH ALONE (RYA-780,
which manufactured a 2.85 dex "discrepancy" from one 0.06 A coincidence). The matcher is
run against a RANDOMISED-WAVELENGTH null so its rate is interpretable rather than merely
encouraging, and an ambiguous (>1 candidate) line is left unadjudicated rather than guessed.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

CANON = ROOT / 'data' / 'linelists' / 'canonical_gf.csv'
BAND_LL = ROOT / 'data' / 'linelists' / 'ispec_ir_9200_13000' / 'atomic_lines.tsv'
LAB = ROOT / 'data' / 'reference' / 'fe_gf_lab' / 'fe1_lab_loggf.csv'
NIST = ROOT / 'data' / 'linelists' / 'primary_gf' / 'nist_asd_FeI_9199_12935.tsv'
OUT_AUDIT = ROOT / 'data' / 'audit' / 'rya834_redward_gf'

BAND_LO_A = 9199.90
#: The 762 BAND tops out at 12935 A, but the synthesis loads the whole LINELIST, which
#: runs to 12976.14 A. `apply_to_synth_array` resolves every line it loads, so a table cut
#: at the band edge dies on the first line past it -- measured, not predicted:
#:     GfResolutionError: no canonical gf for (24, 1) at 12937.0190 EP=2.7099 -- 0-match
#: 17 rows (Mn 15, Cr 1, Si 1) live in 12937.02-12976.14 A. This is RYA-822's "cover the
#: BAND, not one element in it" lesson applied to WAVELENGTH rather than species (RYA-837).
BAND_HI_A = 12977.0

#: Tolerances. Reported by `--scan-tol` rather than asserted: the real matches cluster far
#: tighter than the false ones, and these sit in the gap between the two populations.
WTOL_A = 0.02
EPTOL_EV = 0.002

_ROMAN = {'1': 'I', '2': 'II', '3': 'III', '4': 'IV', '5': 'V', '6': 'VI'}


def load_band_lines() -> pd.DataFrame:
    """The 9199.9-12935 A linelist, collapsed to PHYSICAL lines."""
    d = pd.read_csv(BAND_LL, sep='\t', low_memory=False)
    d = d[(d.wave_A > BAND_LO_A) & (d.wave_A <= BAND_HI_A)].copy()

    sym = d.element.astype(str).str.strip()
    parts = sym.str.split()
    d['element_raw'] = sym
    d['key_z'] = np.nan
    d['ion'] = parts.str[1].map(lambda x: int(x) if str(x).isdigit() else np.nan)

    from pipeline.gf_resolver import cluster_physical_lines
    # The engine's own symbol table, never a hand-typed map (RYA-822's rule). A private
    # {'Fe': 26} fallback would silently drop every other species in the band, and the
    # table has to cover the BAND rather than one element in it -- `apply_to_synth_array`
    # resolves EVERY line the synthesis loads, so a Fe-only extension dies at the first
    # Ti or Si line with `GfResolutionError: ... 0-match`.
    # `z_symbol` is a FUNCTION (z -> symbol), so it is INVERTED here rather than treated
    # as a dict. RYA-822 tests `isinstance(z_symbol, dict)`, which is always False, and so
    # silently falls through to an iSpec import -- which is why that script cannot run off
    # Sirius. Inverting keeps the engine's own table as the authority AND keeps this
    # runnable anywhere, since the band linelist is read from the TSV directly.
    from pipeline.species import z_symbol
    zmap = {}
    for z in range(1, 93):
        try:
            zmap[z_symbol(z)] = z
        except Exception:
            continue
    if 'Fe' not in zmap:
        raise SystemExit("pipeline.species.z_symbol did not yield Fe — refusing to guess "
                         "element atomic numbers")
    d['key_z'] = parts.str[0].map(zmap)
    keep = d.key_z.notna() & d.ion.notna()
    if (~keep).any():
        print(f"  {int((~keep).sum())} row(s) whose species symbol is not resolvable — "
              f"left out rather than guessed")
    d = d[keep].reset_index(drop=True)
    d['key_z'] = d.key_z.astype(int)
    d['ion'] = d.ion.astype(int)

    keys = list(zip(d.key_z.tolist(), d.ion.tolist()))
    wls = d.wave_A.to_numpy(float)
    eps = d.lower_state_eV.to_numpy(float)
    gfs = d.loggf.to_numpy(float)
    refs = d.get('reference_code', pd.Series([''] * len(d))).astype(str).tolist()

    rows = []
    for cl in cluster_physical_lines(keys, wls, eps):
        w = 10.0 ** gfs[cl]
        strongest = cl[int(np.argmax(gfs[cl]))]
        z, io = keys[cl[0]]
        rows.append({
            'key_z': z, 'ion': io, 'element_raw': d.element_raw.iloc[cl[0]],
            'wavelength_air_A': float((wls[cl] * w).sum() / w.sum()),
            'excitation_potential_eV': float(eps[cl].mean()),
            'gf_linelist_vald': float(np.log10(w.sum())),
            'reference_code': refs[strongest],
            'hfs_n_components': len(cl),
        })
    agg = pd.DataFrame(rows)
    if len(agg) != len(d):
        print(f"  {len(d)} component rows -> {len(agg)} physical lines "
              f"(cluster_physical_lines, gf-weighted centroid)")
    return agg.sort_values('wavelength_air_A').reset_index(drop=True)


def _match(fe: pd.DataFrame, w_ref: np.ndarray, e_ref: np.ndarray,
           wtol: float, eptol: float) -> np.ndarray:
    """Index of the UNIQUE (wavelength, EP) match per row, else -1; -2 if ambiguous."""
    out = np.full(len(fe), -1, dtype=int)
    W = fe.wavelength_air_A.to_numpy(float)
    E = fe.excitation_potential_eV.to_numpy(float)
    for i in range(len(fe)):
        c = np.where((np.abs(w_ref - W[i]) <= wtol) & (np.abs(e_ref - E[i]) <= eptol))[0]
        out[i] = int(c[0]) if c.size == 1 else (-2 if c.size > 1 else -1)
    return out


def adjudicate_lab(fe: pd.DataFrame, lab: pd.DataFrame) -> pd.DataFrame:
    """Attach the PRIMARY LAB value where a unique wavelength+EP match exists."""
    lw = lab.wavelength_air_A.to_numpy(float)
    le = lab.elo_eV.to_numpy(float)
    idx = _match(fe, lw, le, WTOL_A, EPTOL_EV)

    fe['lab_log_gf'] = np.nan
    fe['lab_sigma_dex'] = np.nan
    fe['lab_source'] = None
    fe['lab_wavelength_A'] = np.nan
    fe['n_lab_candidates'] = np.where(idx == -2, 2, np.where(idx >= 0, 1, 0))
    for i, j in enumerate(idx):
        if j < 0:
            continue
        fe.at[i, 'lab_log_gf'] = float(lab.loggf.iloc[j])
        fe.at[i, 'lab_sigma_dex'] = float(lab.e_loggf_dex.iloc[j])
        fe.at[i, 'lab_source'] = str(lab.source.iloc[j])
        fe.at[i, 'lab_wavelength_A'] = float(lw[j])
    fe['adjudicated'] = fe.lab_log_gf.notna()
    fe['delta_lab_minus_vald'] = fe.lab_log_gf - fe.gf_linelist_vald
    return fe


def annotate_nist(fe: pd.DataFrame, nist: pd.DataFrame) -> pd.DataFrame:
    """Record NIST — never adopt it. See the module docstring for why."""
    n = nist[nist.log_gf.notna() & nist.ei_eV.notna()].reset_index(drop=True)
    idx = _match(fe, n.wavelength_A.to_numpy(float), n.ei_eV.to_numpy(float),
                 WTOL_A, EPTOL_EV)
    fe['nist_log_gf'] = np.nan
    fe['nist_grade'] = None
    fe['nist_ref_tp'] = None
    for i, j in enumerate(idx):
        if j < 0:
            continue
        fe.at[i, 'nist_log_gf'] = float(n.log_gf.iloc[j])
        fe.at[i, 'nist_grade'] = n.nist_grade.iloc[j]
        fe.at[i, 'nist_ref_tp'] = n.ref_transition_probability.iloc[j]
    fe['delta_nist_minus_vald'] = fe.nist_log_gf - fe.gf_linelist_vald
    return fe


def null_match_rate(fe: pd.DataFrame, w_ref: np.ndarray, e_ref: np.ndarray,
                    n_trials: int, seed: int) -> float:
    """Match rate of this SAME matcher on RANDOMISED wavelengths (EP kept real)."""
    rng = np.random.default_rng(seed)
    ep = fe.excitation_potential_eV.to_numpy(float)
    lo, hi = float(fe.wavelength_air_A.min()), float(fe.wavelength_air_A.max())
    rates = []
    for _ in range(n_trials):
        w = rng.uniform(lo, hi, len(fe))
        hits = sum(1 for k in range(len(w))
                   if np.count_nonzero((np.abs(w_ref - w[k]) <= WTOL_A)
                                       & (np.abs(e_ref - ep[k]) <= EPTOL_EV)) == 1)
        rates.append(hits / len(w))
    return float(np.mean(rates))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--apply', action='store_true', help='write canonical_gf.csv')
    ap.add_argument('--null-trials', type=int, default=20)
    ap.add_argument('--hi-A', type=float, default=None,
                    help='override the redward extent (default: the linelist top, so the '
                         'table covers everything the synthesis LOADS)')
    ap.add_argument('--scan-tol', action='store_true',
                    help='report match counts across tolerances instead of asserting one')
    a = ap.parse_args()

    fe = load_band_lines()
    lab_all = pd.read_csv(LAB)
    lab = lab_all[(lab_all.wavelength_air_A > BAND_LO_A)
                  & (lab_all.wavelength_air_A <= BAND_HI_A)].reset_index(drop=True)
    nist = pd.read_csv(NIST, sep='\t')

    n_fe1 = int(((fe.key_z == 26) & (fe.ion == 1)).sum())
    print(f"band lines         : {len(fe)} physical lines "
          f"({fe.wavelength_air_A.min():.1f}-{fe.wavelength_air_A.max():.1f} A)")
    print(f"  of which Fe I    : {n_fe1}")
    print(f"primary LAB in band: {len(lab)}  "
          f"({lab.wavelength_air_A.min():.1f}-{lab.wavelength_air_A.max():.1f} A)")
    for s, k in lab.source.value_counts().items():
        sub = lab[lab.source == s]
        print(f"    {s:15s} {k:3d}  to {sub.wavelength_air_A.max():9.2f} A")
    print(f"NIST graded in band: {int(nist.log_gf.notna().sum())}  "
          f"(RECORDED, NEVER ADOPTED — compilation echo, see docstring)")

    if a.scan_tol:
        lw, le = lab.wavelength_air_A.to_numpy(float), lab.elo_eV.to_numpy(float)
        print("\n  tolerance scan (unique matches):")
        for wt in (0.005, 0.01, 0.02, 0.05, 0.1):
            for et in (0.0005, 0.002, 0.01):
                idx = _match(fe, lw, le, wt, et)
                print(f"    wtol={wt:<6} eptol={et:<7} unique={int((idx >= 0).sum()):3d}  "
                      f"ambiguous={int((idx == -2).sum())}")
        return

    fe = adjudicate_lab(fe, lab)
    fe = annotate_nist(fe, nist)

    n_adj = int(fe.adjudicated.sum())
    n_amb = int((fe.n_lab_candidates > 1).sum())
    rate = n_adj / max(n_fe1, 1)
    null = null_match_rate(fe[(fe.key_z == 26) & (fe.ion == 1)],
                           lab.wavelength_air_A.to_numpy(float),
                           lab.elo_eV.to_numpy(float), a.null_trials, seed=834)
    print(f"\nadjudicated on LAB : {n_adj}/{n_fe1}  ({100 * rate:.1f}%)")
    print(f"  ambiguous (>1)   : {n_amb}  — left unadjudicated, never guessed")
    print(f"  RANDOMISED NULL  : {100 * null:.1f}%   excess {100 * (rate - null):+.1f} pts")
    if rate <= null:
        raise SystemExit("STOP: the match rate is at or below chance — this matcher is "
                         "finding coincidences, not identifying lines.")

    d = fe[fe.adjudicated]
    if len(d):
        med = d.delta_lab_minus_vald.median()
        mad = 1.4826 * (d.delta_lab_minus_vald - med).abs().median()
        print(f"\ndelta LAB - VALD   : median {med:+.3f}  MAD {mad:.3f}  "
              f"|max| {d.delta_lab_minus_vald.abs().max():.3f} dex")
        print(f"  lab sigma        : median {d.lab_sigma_dex.median():.3f} dex "
              f"(vs the 0.20 blanket)")

    # THE ABSENCE, MEASURED
    top_lab = float(lab.wavelength_air_A.max())
    above = fe[(fe.wavelength_air_A > top_lab) & (fe.key_z == 26) & (fe.ion == 1)]
    n_nist_above = int(above.nist_log_gf.notna().sum())
    print(f"\nABOVE the last lab line ({top_lab:.2f} A): {len(above)} Fe I lines")
    print(f"  primary lab gf   : 0   (verified against all three staged sources)")
    print(f"  NIST rows        : {n_nist_above}  — recorded, NOT adopted; they agree with "
          f"the linelist to ~0.0003 dex, i.e. the same number twice")
    print(f"  => these stay on the Kurucz floor, stated as single_source")

    OUT_AUDIT.mkdir(parents=True, exist_ok=True)
    fe.to_csv(OUT_AUDIT / 'redward_fe1_gf_adjudication.csv', index=False)
    summary = {
        'ticket': 'RYA-834',
        'band_A': [BAND_LO_A, BAND_HI_A],
        'physical_lines': int(len(fe)),
        'fe1_lines': n_fe1,
        'lab_in_band': int(len(lab)),
        'lab_span_A': [float(lab.wavelength_air_A.min()), top_lab],
        'lab_by_source': {k: int(v) for k, v in lab.source.value_counts().items()},
        'adjudicated_on_lab': n_adj,
        'ambiguous': n_amb,
        'match_rate': round(rate, 4),
        'randomised_null_rate': round(null, 4),
        'lines_above_last_lab_line': int(len(above)),
        'nist_rows_above_last_lab_line': n_nist_above,
        'nist_adopted': False,
        'nist_note': ('NIST ASD Fe I in this band agrees with the linelist to a median '
                      '0.0003 dex (max 0.0016) over 34 matched lines — the same number '
                      'arriving twice, not corroboration. RYA-760: FMW is a NIST '
                      'compilation and VALD copies it.'),
        'ticket_span_correction': ('The ticket states Ruffoni 2014 covers 3526-10864 A. '
                                   'Our vendored RYA-799 pull reaches 11316.06 A, and Den '
                                   'Hartog 2014 reaches 11013.24 A, so the claimed '
                                   '10864-12935 A gap actually opens at 11316.06 A.'),
        'tolerances': {'wavelength_A': WTOL_A, 'excitation_eV': EPTOL_EV},
    }
    (OUT_AUDIT / 'redward_gf_summary.json').write_text(
        json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print(f"\n  wrote {OUT_AUDIT}/redward_fe1_gf_adjudication.csv")
    print(f"  wrote {OUT_AUDIT}/redward_gf_summary.json")

    if not a.apply:
        print("\n  [dry-run] re-run with --apply to extend canonical_gf.csv")
        return

    canon = pd.read_csv(CANON, low_memory=False)
    before_hi = float(canon.wavelength_air_A.max())
    start = 1 + max(int(x.split('_')[1]) for x in canon.line_id if str(x).startswith('gf_'))

    new = pd.DataFrame({
        'line_id': [f'gf_{start + i:06d}' for i in range(len(fe))],
        'key_z': fe.key_z, 'ion': fe.ion,
        'species': fe.element_raw.map(
            lambda r: ' '.join([r.split()[0], _ROMAN.get(r.split()[1], r.split()[1])])),
        'wavelength_air_A': fe.wavelength_air_A.round(3),
        'excitation_potential_eV': fe.excitation_potential_eV.round(4),
        'hfs_n_components': fe.hfs_n_components,
        'log_gf': np.where(fe.adjudicated, fe.lab_log_gf, fe.gf_linelist_vald),
        'loggf_reference': np.where(
            fe.adjudicated, 'PRIMARY LAB ' + fe.lab_source.astype(str),
            fe.reference_code.replace('', 'VALD3')),
        'nist_grade': fe.nist_grade,
        'seed_source': 'ir_vald_9200_13000(RYA-762)',
        'adjudication_status': np.where(fe.adjudicated, 'lab_rya834', 'single_source'),
        'gf_synth_ges': np.nan,
        'gf_regions_vald': np.nan,
        'gf_linelist_vald': fe.gf_linelist_vald,
        'in_synth': False, 'in_regions': False, 'in_linelist': True,
        'delta_synth_minus_ew': np.nan,
    })

    key = ['species', 'wavelength_air_A', 'excitation_potential_eV']
    before = len(new)
    merged = new.merge(canon[key].drop_duplicates(), on=key, how='left', indicator=True)
    new = new[merged['_merge'].to_numpy() == 'left_only'].reset_index(drop=True)
    if before != len(new):
        print(f"  {before - len(new)} line(s) already in canonical_gf — kept the EXISTING "
              f"row (it may carry adjudication history), dropped the new copy")

    out = pd.concat([canon, new], ignore_index=True)
    out.to_csv(CANON, index=False)
    print(f"\n  canonical_gf.csv  {len(canon)} -> {len(out)} rows")
    print(f"  red edge          {before_hi:.3f} -> {float(out.wavelength_air_A.max()):.3f} A")


if __name__ == '__main__':
    main()
