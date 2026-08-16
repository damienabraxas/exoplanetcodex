#!/usr/bin/env python3
"""
RYA-822 — extend canonical_gf.csv BLUEWARD of 3780.04 A to 3000 A for Fe I, and adjudicate.

WHY. RYA-759 measured the near-UV Fe I product at 7.487 with a line-to-line scatter of
0.354 dex whose cause correlates with NOTHING measurable (wavelength, chi2r, continuum,
depth, EP, log gf all |r| <= 0.21), and best/worst fit quartiles differ by only 0.08 dex.
The one thing left is PER-LINE gf error -- and `canonical_gf.csv` starts at 3780.04 A, so
this band has never had a single line adjudicated. Every gf in it is VALD3 as delivered.

WHAT ADJUDICATES IT. NIST ASD is the only graded gf source that reaches this band:

  * Ruffoni-2014      4423.84-10863.52 A  -> ZERO lines below 3780. Checked, not assumed.
  * Den Hartog-2014   7016.39- 9214.50 A  -> ZERO lines below 3780. Checked, not assumed.
  * Peterson & Kurucz 2022 (PK21)         -> 860 lines in band, but its `dgf` is a
    REFERENCE-FIT offset (log(gf)0 = log gf + dgf, "best matches observed spectra"), which
    is the RYA-161 tuning-firewall violation RYA-779 ratified against. Its CALCULATED
    log(gf) is ordinary Kurucz theory -- the same class as the K07/K10/K13 rows already in
    canonical_gf -- so it is admissible but can never outrank a graded measurement. PK21 is
    therefore used for LINE IDENTIFICATION, not adopted as gf here.

⚠️ MATCHING IS ON WAVELENGTH *AND* EXCITATION POTENTIAL, NEVER WAVELENGTH ALONE (RYA-780).
This band carries ~5.6 Fe I lines per Angstrom, so a +/-0.02 A wavelength coincidence
happens about 20% of the time BY CHANCE -- measured against a randomised null, not assumed.
Adopting a value on a wavelength coincidence would write the wrong gf onto the wrong line
and look perfectly clean doing it. A match must ALSO agree in EP, and must be UNIQUE:
if two NIST rows fit one line, that line is left unadjudicated rather than guessed.
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
NIST = ROOT / 'data' / 'linelists' / 'primary_gf' / 'nist_asd_FeI_3000_3780.tsv'
NEARUV_LL = ROOT / 'data' / 'linelists' / 'ispec_nearuv_3000_3780' / 'atomic_lines.tsv'
OUT_AUDIT = ROOT / 'data' / 'audit' / 'nearuv_gf_adjudication'

# ⚠️ THESE ARE NOT TUNED TO AN ANSWER — THEY ARE READ OFF A BIMODALITY IN THE MATCHES.
# At the loose tolerances first tried (0.02 A / 0.05 eV) the 668 matches split cleanly:
#
#   632 lines agreeing to |dgf| <= 0.05 :  |dwave| median 0.0003 A (max 0.0051)
#                                          |dEP|   median 0.00002 eV (max 0.00005)
#     9 lines disagreeing by |dgf| > 0.3 :  |dwave| p90 0.0112 A (max 0.0194)
#                                          |dEP|   p90 0.0304 eV (max 0.0361)
#
# A genuine identification reproduces NIST's own wavelength and energy to MACHINE
# precision. The nine wild disagreements sit against the tolerance WALLS in both
# coordinates, and there are exactly nine — which is what the randomised null (0.2% of
# 4,364) predicts for chance pairings. So they are not "VALD is wrong by 3.8 dex", they
# are the matcher pairing unrelated lines, and adopting them would have written a wrong
# gf onto a real line with full provenance attached.
WTOL_A = 0.006      # 10x the median real residual, still well inside the false-match gap
EPTOL_EV = 0.0005   # 10x the max real residual (5e-5 eV); the false matches are 600x that


def load_nearuv_lines(fe_only: bool = False) -> pd.DataFrame:
    """Every species in the near-UV list, collapsed to physical lines.

    ⚠️ SCOPE IS WIDER THAN THE TICKET ASKED, AND IT HAS TO BE. RYA-822 specifies
    "covering the near-UV Fe I pool", but `gf_resolver.apply_to_synth_array` resolves
    EVERY line the synthesis loads, not just the target element's. With Fe I alone
    extended, turning canonical single-sourcing on for this band dies at the first
    other species:

        GfResolutionError: no canonical gf for (24, 2) at 3118.6490 EP=2.4212 — 0-match

    So the table must cover the BAND, not one element in it, or the re-derivation the
    ticket asks for in step 4 cannot run at all. Non-Fe species are seeded from VALD as
    delivered and marked `single_source` — exactly the status 139,782 of the 145,886
    existing optical rows already carry. Only Fe I is adjudicated, because NIST was
    pulled for Fe I.
    """
    from pipeline.nearuv_synth import build_solar_context
    ctx = build_solar_context('Fe', 4e5, linelist_file=str(NEARUV_LL),
                              apply_canonical_gf=False)
    ll = ctx['linelist']
    el = np.array([str(x).strip() for x in ll['element']])
    m = (el == 'Fe 1') if fe_only else np.ones(len(el), bool)
    d = pd.DataFrame({
        'element_raw': el[m],
        'wavelength_air_A': np.asarray(ll['wave_A'], float)[m],
        'excitation_potential_eV': np.asarray(ll['lower_state_eV'], float)[m],
        'gf_linelist_vald': np.asarray(ll['loggf'], float)[m],
        'reference_code': [str(x).strip() for x in np.asarray(ll['reference_code'])[m]],
    })
    # 'Fe 1' -> (Z=26, ion=1). The engine's own symbol table, never a hand-typed map.
    from pipeline.species import z_symbol
    sym2z = {v: k for k, v in z_symbol.items()} if isinstance(z_symbol, dict) else None
    if sym2z is None:
        import ispec
        chem = ctx['chem_elements']
        sym2z = {str(r['symbol']).strip(): int(r['atomic_num']) for r in chem}
    zs, ions, keep = [], [], []
    for raw in d.element_raw:
        parts = str(raw).split()
        ok = len(parts) == 2 and parts[0] in sym2z and parts[1].isdigit()
        keep.append(ok)
        zs.append(sym2z.get(parts[0], -1) if ok else -1)
        ions.append(int(parts[1]) if ok else -1)
    d['key_z'], d['ion'] = zs, ions
    n_drop = int(len(d) - sum(keep))
    if n_drop:
        print(f"  {n_drop} row(s) whose species symbol is not in the engine's table — "
              f"left out rather than guessed")
    d = d[pd.Series(keep).values].reset_index(drop=True)
    # ⚠️ COLLAPSE TO PHYSICAL LINES WITH THE PROJECT'S OWN GROUPER, NOT WITH ROUNDING.
    # `pipeline.gf_resolver.cluster_physical_lines` is THE shared physical-line grouper
    # (RYA-353) — its own docstring says the audit, the canonical build and the synth
    # reroute must all agree through it. My first pass grouped by rounding (lambda, EP)
    # instead, which produced rows the resolver could not find: `apply_to_synth_array`
    # looks up a cluster's gf-WEIGHTED CENTROID and MEAN EP, and a rounded raw
    # wavelength is not that coordinate. It failed loudly —
    #     GfResolutionError: no canonical gf for (26, 1) at 3214.0366 EP=2.4496 — 0-match
    # — which is the only reason it was caught before the numbers were believed.
    from pipeline.gf_resolver import cluster_physical_lines

    keys = list(zip(d.key_z.tolist(), d.ion.tolist()))
    wls = d.wavelength_air_A.to_numpy(float)
    eps = d.excitation_potential_eV.to_numpy(float)
    gfs = d.gf_linelist_vald.to_numpy(float)
    refs = d.reference_code.tolist()
    els = d.element_raw.tolist()

    rows = []
    for cl in cluster_physical_lines(keys, wls, eps):
        w = 10.0 ** gfs[cl]
        strongest = cl[int(np.argmax(gfs[cl]))]
        z, io = keys[cl[0]]
        rows.append({
            'key_z': z, 'ion': io, 'element_raw': els[cl[0]],
            'wavelength_air_A': float((wls[cl] * w).sum() / w.sum()),
            'excitation_potential_eV': float(eps[cl].mean()),
            'gf_linelist_vald': float(np.log10(w.sum())),
            'reference_code': refs[strongest],
            'hfs_n_components': len(cl),
        })
    agg = pd.DataFrame(rows)
    if len(agg) != len(d):
        print(f"  {len(d)} VALD component rows -> {len(agg)} physical lines "
              f"(cluster_physical_lines, gf-weighted centroid)")
    return agg.sort_values('wavelength_air_A').reset_index(drop=True)


def adjudicate(fe: pd.DataFrame, nist: pd.DataFrame) -> pd.DataFrame:
    """Attach the NIST value where a UNIQUE wavelength+EP match exists."""
    n = nist[nist.log_gf.notna() & nist.ei_eV.notna()].reset_index(drop=True)
    nw = n.wavelength_A.to_numpy(float)
    ne = n.ei_eV.to_numpy(float)

    cols = ['nist_log_gf', 'nist_grade', 'nist_acc_pct', 'nist_ref_tp', 'nist_ref_line',
            'nist_wavelength_A', 'nist_ei_eV']
    for c in cols:
        fe[c] = np.nan if c not in ('nist_grade', 'nist_ref_tp', 'nist_ref_line') else None
    fe['n_nist_candidates'] = 0

    is_fe1 = (fe.key_z == 26) & (fe.ion == 1)
    for i, r in fe[is_fe1].iterrows():
        w, ep = float(r.wavelength_air_A), float(r.excitation_potential_eV)
        cand = np.where((np.abs(nw - w) <= WTOL_A) & (np.abs(ne - ep) <= EPTOL_EV))[0]
        fe.at[i, 'n_nist_candidates'] = int(cand.size)
        if cand.size != 1:
            continue                      # 0 = no match; >1 = ambiguous, never guess
        j = int(cand[0])
        fe.at[i, 'nist_log_gf'] = float(n.log_gf[j])
        fe.at[i, 'nist_grade'] = n.nist_grade[j]
        fe.at[i, 'nist_acc_pct'] = float(n.nist_acc_pct[j]) if pd.notna(n.nist_acc_pct[j]) else np.nan
        fe.at[i, 'nist_ref_tp'] = n.ref_transition_probability[j]
        fe.at[i, 'nist_ref_line'] = n.ref_line[j]
        fe.at[i, 'nist_wavelength_A'] = float(nw[j])
        fe.at[i, 'nist_ei_eV'] = float(ne[j])

    fe['adjudicated'] = fe.nist_log_gf.notna()
    fe['delta_nist_minus_vald'] = fe.nist_log_gf - fe.gf_linelist_vald
    return fe


def null_match_rate(fe: pd.DataFrame, nist: pd.DataFrame, n_trials: int, seed: int) -> float:
    """How often does this SAME matcher fire on randomised wavelengths?

    Without this the match rate is uninterpretable: dense bands match by accident. The
    EP coordinate is kept real and only the wavelength is randomised, which is the
    conservative choice -- it isolates exactly the coincidence being worried about.
    """
    n = nist[nist.log_gf.notna() & nist.ei_eV.notna()]
    nw, ne = n.wavelength_A.to_numpy(float), n.ei_eV.to_numpy(float)
    rng = np.random.default_rng(seed)
    ep = fe.excitation_potential_eV.to_numpy(float)
    lo, hi = float(fe.wavelength_air_A.min()), float(fe.wavelength_air_A.max())
    rates = []
    for _ in range(n_trials):
        w = rng.uniform(lo, hi, len(fe))
        hits = sum(1 for k in range(len(w))
                   if np.count_nonzero((np.abs(nw - w[k]) <= WTOL_A)
                                       & (np.abs(ne - ep[k]) <= EPTOL_EV)) == 1)
        rates.append(hits / len(w))
    return float(np.mean(rates))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--apply', action='store_true', help='write canonical_gf.csv')
    ap.add_argument('--null-trials', type=int, default=20)
    a = ap.parse_args()

    fe = load_nearuv_lines()
    nist = pd.read_csv(NIST, sep='\t')
    n_fe1 = int(((fe.key_z == 26) & (fe.ion == 1)).sum())
    print(f"near-UV lines      : {len(fe)} physical lines, {fe.element_raw.nunique()} "
          f"species  ({fe.wavelength_air_A.min():.1f}-{fe.wavelength_air_A.max():.1f} A)")
    print(f"  of which Fe I    : {n_fe1}")
    print(f"NIST graded rows   : {int(nist.log_gf.notna().sum())}")

    fe = adjudicate(fe, nist)
    n_adj = int(fe.adjudicated.sum())
    n_amb = int((fe.n_nist_candidates > 1).sum())
    rate = n_adj / max(n_fe1, 1)
    null = null_match_rate(fe[(fe.key_z == 26) & (fe.ion == 1)], nist,
                           a.null_trials, seed=822)
    print(f"\nadjudicated (Fe I) : {n_adj}/{n_fe1}  ({100 * rate:.1f}%)")
    print(f"  ambiguous (>1)   : {n_amb}  — left unadjudicated, never guessed")
    print(f"  RANDOMISED NULL  : {100 * null:.1f}%   excess {100 * (rate - null):+.1f} pts")
    if rate <= null:
        raise SystemExit("STOP: the match rate is at or below chance — this matcher is "
                         "not identifying lines, it is finding coincidences.")

    d = fe[fe.adjudicated]
    print(f"\ndelta NIST - VALD  : median {d.delta_nist_minus_vald.median():+.3f}  "
          f"MAD {1.4826 * (d.delta_nist_minus_vald - d.delta_nist_minus_vald.median()).abs().median():.3f}  "
          f"|max| {d.delta_nist_minus_vald.abs().max():.3f} dex")
    print(f"  grades           : {dict(d.nist_grade.value_counts())}")

    OUT_AUDIT.mkdir(parents=True, exist_ok=True)
    fe.to_csv(OUT_AUDIT / 'nearuv_fe1_gf_adjudication.csv', index=False)

    canon = pd.read_csv(CANON, low_memory=False)
    before_lo = float(canon.wavelength_air_A.min())
    existing = set(canon.line_id)
    start = 1 + max(int(x.split('_')[1]) for x in existing if str(x).startswith('gf_'))

    new = pd.DataFrame({
        'line_id': [f'gf_{start + i:06d}' for i in range(len(fe))],
        'key_z': fe.key_z, 'ion': fe.ion,
        'species': fe.element_raw.map(lambda r: ' '.join(
            [r.split()[0], {'1': 'I', '2': 'II', '3': 'III', '4': 'IV', '5': 'V',
                            '6': 'VI'}.get(r.split()[1], r.split()[1])])),
        'wavelength_air_A': fe.wavelength_air_A.round(3),
        'excitation_potential_eV': fe.excitation_potential_eV.round(4),
        'hfs_n_components': fe.hfs_n_components,
        'log_gf': np.where(fe.adjudicated, fe.nist_log_gf, fe.gf_linelist_vald),
        'loggf_reference': np.where(fe.adjudicated,
                                    'NIST ASD (astroquery) ' + fe.nist_ref_tp.astype(str),
                                    fe.reference_code.replace('', 'VALD3')),
        'nist_grade': fe.nist_grade,
        'seed_source': 'nearuv_vald(RYA-759)',
        'adjudication_status': np.where(fe.adjudicated, 'nist_rya822', 'single_source'),
        'gf_synth_ges': np.nan,            # GES v6 stops at 4200 A — never reaches here
        'gf_regions_vald': np.nan,         # the 42000-region file is optical
        'gf_linelist_vald': fe.gf_linelist_vald,
        'in_synth': False, 'in_regions': False, 'in_linelist': True,
        'delta_synth_minus_ew': np.nan,
    })

    # A near-UV line that ALREADY has a canonical row (the bands abut at 3780.04 A, and
    # rounding can put a line on either side) must not be added a second time. Drop the
    # new copy rather than the old: the existing row may carry adjudication history.
    key = ['species', 'wavelength_air_A', 'excitation_potential_eV']
    have = set(map(tuple, canon[key].round({'wavelength_air_A': 3,
                                            'excitation_potential_eV': 4}).values))
    is_new = ~pd.Series(list(map(tuple, new[key].values))).isin(have).values
    n_drop = int((~is_new).sum())
    if n_drop:
        print(f"  {n_drop} line(s) already had a canonical row — kept the existing one")
    new = new[is_new].reset_index(drop=True)
    out = pd.concat([canon, new], ignore_index=True)

    dup = out.duplicated(subset=['species', 'wavelength_air_A',
                                 'excitation_potential_eV'], keep=False)
    if int(dup.sum()) > int(canon.duplicated(subset=['species', 'wavelength_air_A',
                                                     'excitation_potential_eV'],
                                             keep=False).sum()):
        raise SystemExit(f"STOP: the extension introduces {int(dup.sum())} duplicate "
                         f"(species, lambda, EP) keys — canonical_gf is one row per "
                         f"physical line and a duplicate would make the resolver "
                         f"non-deterministic.")

    print(f"\ncanonical_gf rows  : {len(canon)} -> {len(out)}  (+{len(new)})")
    print(f"blue edge          : {before_lo:.3f} -> {out.wavelength_air_A.min():.3f} A")

    prov = {
        'ticket': 'RYA-822', 'n_added': int(len(new)),
        'band_A': [3000.0, 3780.0], 'species': 'Fe I',
        'n_adjudicated_nist': n_adj, 'n_ambiguous_left_alone': n_amb,
        'match_rate': rate, 'randomised_null_rate': null,
        'match_tolerances': {'wavelength_A': WTOL_A, 'excitation_potential_eV': EPTOL_EV},
        'delta_nist_minus_vald_median': float(d.delta_nist_minus_vald.median()),
        'delta_nist_minus_vald_absmax': float(d.delta_nist_minus_vald.abs().max()),
        'sources_checked_and_rejected': {
            'ruffoni2014': 'span 4423.84-10863.52 A — ZERO lines below 3780',
            'denhartog2014': 'span 7016.39-9214.50 A — ZERO lines below 3780',
            'PK21_dgf': 'reference-fit astrophysical offset — RYA-161 firewall, RYA-779 '
                        'ratified NOT a graded gf source; used for line ID only',
        },
        'caveat': 'NIST Fe I here is largely a COMPILATION; agreement proves only "no '
                  'transcription error", DISagreement is the informative outcome. The '
                  'per-line TP/Line reference codes are carried so the underlying '
                  'measurement stays visible.',
    }
    (OUT_AUDIT / 'nearuv_fe1_gf_adjudication.prov.json').write_text(json.dumps(prov, indent=2))

    if a.apply:
        out.to_csv(CANON, index=False)
        print(f"\nWROTE {CANON.relative_to(ROOT)}")
    else:
        print("\n[dry-run] re-run with --apply to write canonical_gf.csv")


if __name__ == '__main__':
    main()
