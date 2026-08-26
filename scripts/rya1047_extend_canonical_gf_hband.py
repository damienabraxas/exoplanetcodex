#!/usr/bin/env python3
"""
RYA-1047 — extend canonical_gf.csv from 12975.9 A into the H band, and adjudicate.

WHY. RYA-1046 ingested Ruffoni 2013 -- 28 primary-laboratory Fe I log gf at 14680-21386 A,
the only lab Fe I gf we hold above 11316 A -- and it graded NOTHING, because canonical_gf
stops at 12975.9 A and every one of those 28 lines is above it. 0/28 present; a 1745 A gap
with zero overlap. This is the third time a hard wavelength bound has stopped real work
(RYA-762's products, then 6 Mn I lines 0.09 A past the red edge disabling lab gf for 4,364
Fe I lines, now this), so the bound goes.

🔴 THE SOURCE IS A FILE MARKED `hfsoff_quarantine`, AND THAT IS ARGUED, NOT WAVED THROUGH.
`vald_55cnc_nir_5k30k_hfsoff_quarantine.txt` is a VALD extract-stellar pull for 55 Cnc,
5000-29995 A, taken with hyperfine structure OFF. Two things could disqualify it: that it
is star-specific, and that hfs is off. Both were MEASURED against the list we already ship
and trust for NIR products (`ispec_ir_9200_13000/atomic_lines.tsv`, itself RYA-762's VALD
extract), over their 9203-12976 A overlap:

    638 unambiguous 1-1 matches
    EP     : 637/638 EXACT
    log gf : 566/638 exact, 72 discordant

and every one of the 72 discordant lines is K I, Na I, Al I, Co I or Mn I -- odd-Z species
with nuclear spin. Same wavelength, IDENTICAL EP, different log gf: the hfs signature, one
physical line split into components that divide the total gf.

    Fe lines compared: 236        Fe lines discordant: 0

That is the physics, not a coincidence: 56Fe is 91.75% of iron and has nuclear spin I = 0,
so Fe I HAS no hyperfine structure to switch off. The quarantine flag is REAL and material
for those five species and physically vacuous for Fe.

⇒ THIS SCRIPT EXTENDS Fe ONLY. The other species stay out, and stay out for a stated
reason rather than by omission. Star-specificity is handled the same way: what 55 Cnc's
depth threshold changed is which lines were SELECTED, and we re-select per band anyway;
the atomic data underneath is VALD's and is byte-identical where the two lists overlap.

MATCHING IS ON WAVELENGTH *AND* EXCITATION POTENTIAL, NEVER WAVELENGTH ALONE (RYA-780,
which manufactured a 2.85 dex "discrepancy" out of one 0.06 A coincidence), it is run
against a RANDOMISED-WAVELENGTH null so the rate is interpretable rather than merely
encouraging, and an ambiguous (>1 candidate) line is left UNADJUDICATED rather than
guessed at. All three are RYA-834's design, reused deliberately.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

CANON = ROOT / 'data' / 'linelists' / 'canonical_gf.csv'
VALD = ROOT / 'data' / 'linelists' / 'vald_55cnc_nir_5k30k_hfsoff_quarantine.txt'
TRUSTED = ROOT / 'data' / 'linelists' / 'ispec_ir_9200_13000' / 'atomic_lines.tsv'
OUT_AUDIT = ROOT / 'data' / 'audit' / 'rya1047_hband_gf'

# The window: from canonical's existing red edge to just past Ruffoni's reddest line.
LO_A, HI_A = 12976.0, 21400.0

# RYA-834's tolerances, unchanged. A wavelength tolerance loose enough to admit a second
# candidate is what RYA-780 fell into; EP is what separates same-wavelength coincidences.
WTOL_A = 0.02
EPTOL_EV = 0.002

# The provenance check above, as constants the code re-derives rather than trusts.
HFS_SENSITIVE = ('K', 'Na', 'Al', 'Co', 'Mn')

_ROMAN = {'1': 'I', '2': 'II', '3': 'III'}
_DATA = re.compile(r"^'([A-Z][a-z]?) (\d+)',\s*([\d.]+),\s*(-?[\d.]+),\s*([\d.]+),"
                   r"\s*(-?[\d.]+),\s*([\d.]+),")
_GFREF = re.compile(r'gf:(\S+)')


def load_vald(path: Path) -> pd.DataFrame:
    """Parse VALD 'extract stellar' long format: a data line, two term lines, a reference
    line beginning `'_`. The gf provenance code lives on the reference line as `gf:CODE`."""
    # 🔴 EVERY DATA LINE PRODUCES A ROW, whatever follows it. An earlier version of this
    # parser only emitted a row once it had seen the `'_` reference line, and silently
    # dropped ~3000 lines whose trailer did not match. That is not a cosmetic loss: a
    # DROPPED NEIGHBOUR MAKES AN AMBIGUOUS PAIR LOOK LIKE A CLEAN 1-1 MATCH, which is
    # precisely the RYA-780 failure the wavelength+EP matcher exists to prevent -- and it
    # showed up as C and Mg appearing to disagree on log gf when neither has hyperfine
    # structure to explain it. The reference code is optional metadata; the line is not.
    rows, pending = [], None
    for raw in path.open(encoding='latin-1'):
        m = _DATA.match(raw)
        if m:
            if pending is not None:
                rows.append(pending)
            pending = dict(element=m.group(1), ion=int(m.group(2)),
                           wavelength_air_A=float(m.group(3)),
                           gf_linelist_vald=float(m.group(4)),
                           excitation_potential_eV=float(m.group(5)),
                           reference_code='')
            continue
        if pending is not None and raw.startswith("'_") and not pending['reference_code']:
            g = _GFREF.search(raw)
            pending['reference_code'] = g.group(1) if g else ''
    if pending is not None:
        rows.append(pending)
    return pd.DataFrame(rows)


def verify_provenance(v: pd.DataFrame) -> dict:
    """Re-derive the hfs argument in the docstring rather than asserting it.

    🔴 A claim recorded in a comment is not a check. If this ever stops holding -- a new
    pull, a different VALD version -- the script must refuse, not carry an obsolete
    justification forward in prose."""
    t = pd.read_csv(TRUSTED, sep='\t')
    t = t[(t.wave_A >= 9203.0) & (t.wave_A <= 12976.0)]
    vw = v.wavelength_air_A.to_numpy(float)
    recs = []
    for _, r in t.iterrows():
        c = np.where(np.abs(vw - float(r.wave_A)) < 0.01)[0]
        if c.size == 1:
            j = int(c[0])
            recs.append((str(r.element).split()[0],
                         abs(float(r.loggf) - float(v.gf_linelist_vald.iloc[j])),
                         abs(float(r.lower_state_eV) - float(v.excitation_potential_eV.iloc[j]))))
    d = pd.DataFrame(recs, columns=['elem', 'dgf', 'dep'])
    bad = d[d.dgf > 0.001]
    fe = d[d.elem == 'Fe']
    n_fe_bad = int((fe.dgf > 0.001).sum())

    print(f"  PROVENANCE CHECK vs the trusted NIR list, over their 9203-12976 A overlap")
    print(f"    {len(d)} unambiguous 1-1 matches")
    print(f"    EP exact      : {int((d.dep < 1e-6).sum())}/{len(d)}")
    print(f"    log gf exact  : {int((d.dgf < 1e-6).sum())}/{len(d)}   discordant {len(bad)}")
    if len(bad):
        print(f"    discordant by species: {bad.elem.value_counts().to_dict()}")
    print(f"    Fe compared   : {len(fe)}   Fe discordant: {n_fe_bad}")

    if n_fe_bad:
        raise SystemExit(
            f"REFUSING: {n_fe_bad} Fe line(s) disagree between this pull and the trusted "
            f"NIR list. The whole argument for using an hfs-off file is that Fe has no "
            f"hyperfine structure (56Fe, I=0) and is therefore untouched. That argument "
            f"just failed its own test -- do not extend canonical on it.")
    off = set(bad.elem) - set(HFS_SENSITIVE)
    if off:
        raise SystemExit(
            f"REFUSING: discordant species {sorted(off)} are outside the hfs-sensitive set "
            f"{HFS_SENSITIVE}. The discordance is then NOT explained by hfs and the "
            f"provenance argument does not hold.")
    return {'matches': len(d), 'ep_exact': int((d.dep < 1e-6).sum()),
            'loggf_exact': int((d.dgf < 1e-6).sum()), 'discordant': len(bad),
            'discordant_species': bad.elem.value_counts().to_dict(),
            'fe_compared': len(fe), 'fe_discordant': n_fe_bad}


def _match(fe: pd.DataFrame, w_ref: np.ndarray, e_ref: np.ndarray) -> np.ndarray:
    """Index of the UNIQUE (wavelength, EP) match per row, else -1; -2 if ambiguous."""
    out = np.full(len(fe), -1, dtype=int)
    W = fe.wavelength_air_A.to_numpy(float)
    E = fe.excitation_potential_eV.to_numpy(float)
    for i in range(len(fe)):
        c = np.where((np.abs(w_ref - W[i]) <= WTOL_A) & (np.abs(e_ref - E[i]) <= EPTOL_EV))[0]
        out[i] = int(c[0]) if c.size == 1 else (-2 if c.size > 1 else -1)
    return out


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
    ap.add_argument('--apply', action='store_true')
    ap.add_argument('--null-trials', type=int, default=20)
    a = ap.parse_args()

    from pipeline.gf_grades import CITATIONS, lab_lines

    OUT_AUDIT.mkdir(parents=True, exist_ok=True)
    print(f"RYA-1047 — canonical_gf H-band extension  {LO_A:.0f}-{HI_A:.0f} A\n")

    v = load_vald(VALD)
    print(f"  parsed {len(v)} lines from {VALD.name}  "
          f"({v.wavelength_air_A.min():.1f}-{v.wavelength_air_A.max():.1f} A)")
    prov = verify_provenance(v)

    fe = v[(v.element == 'Fe') & (v.wavelength_air_A > LO_A)
           & (v.wavelength_air_A <= HI_A)].reset_index(drop=True)
    fe['key_z'] = 26
    fe['hfs_n_components'] = 1        # 56Fe has I=0: one physical line, one component
    print(f"\n  Fe in {LO_A:.0f}-{HI_A:.0f} A: {len(fe)} lines "
          f"(Fe I {int((fe.ion == 1).sum())}, Fe II {int((fe.ion == 2).sum())})")

    lab = lab_lines('Fe I')
    lab = lab[(lab.wavelength_air_A > LO_A) & (lab.wavelength_air_A <= HI_A)].reset_index(drop=True)
    print(f"  primary lab Fe I in window: {len(lab)}  "
          f"({sorted(set(lab.source))})")

    fe1 = fe[fe.ion == 1].reset_index(drop=True)
    idx = _match(fe1, lab.wavelength_air_A.to_numpy(float), lab.elo_eV.to_numpy(float))
    fe1['lab_log_gf'] = np.nan
    fe1['lab_sigma_dex'] = np.nan
    fe1['lab_source'] = None
    fe1['n_lab_candidates'] = np.where(idx == -2, 2, np.where(idx >= 0, 1, 0))
    for i, j in enumerate(idx):
        if j < 0:
            continue
        fe1.at[i, 'lab_log_gf'] = float(lab.loggf.iloc[j])
        fe1.at[i, 'lab_sigma_dex'] = float(lab.e_loggf_dex.iloc[j])
        fe1.at[i, 'lab_source'] = str(lab.source.iloc[j])
    fe1['adjudicated'] = fe1.lab_log_gf.notna()
    fe1['delta_lab_minus_vald'] = fe1.lab_log_gf - fe1.gf_linelist_vald

    n_adj = int(fe1.adjudicated.sum())
    n_amb = int((fe1.n_lab_candidates == 2).sum())
    rate = n_adj / len(lab) if len(lab) else 0.0
    null = null_match_rate(fe1, lab.wavelength_air_A.to_numpy(float),
                           lab.elo_eV.to_numpy(float), a.null_trials, seed=1047)
    print(f"\n  ADJUDICATED      : {n_adj}/{len(lab)} lab lines matched ({100 * rate:.1f}%)")
    print(f"  ambiguous (left alone): {n_amb}")
    print(f"  RANDOMISED NULL  : {100 * null:.1f}%   excess {100 * (rate - null):+.1f} pts")
    if rate <= null:
        print("  🔴 match rate does NOT beat the null — the matching is not informative")
    if n_adj:
        d = fe1[fe1.adjudicated].delta_lab_minus_vald
        print(f"  lab - VALD       : median {d.median():+.3f} dex, "
              f"|max| {d.abs().max():.3f}, n={len(d)}")

    fe1.to_csv(OUT_AUDIT / 'hband_fe1_gf_adjudication.csv', index=False)
    summary = {
        'ticket': 'RYA-1047', 'window_A': [LO_A, HI_A],
        'source': VALD.name, 'provenance_check': prov,
        'fe_lines': int(len(fe)), 'fe1_lines': int(len(fe1)),
        'lab_lines_in_window': int(len(lab)),
        'adjudicated': n_adj, 'ambiguous_left_unadjudicated': n_amb,
        'match_rate': round(rate, 4), 'randomised_null_rate': round(null, 4),
        'null_trials': a.null_trials,
    }
    (OUT_AUDIT / 'hband_gf_summary.json').write_text(json.dumps(summary, indent=2) + '\n',
                                                     encoding='utf-8')
    print(f"\n  wrote {OUT_AUDIT}/hband_fe1_gf_adjudication.csv")
    print(f"  wrote {OUT_AUDIT}/hband_gf_summary.json")

    if not a.apply:
        print("\n  [dry-run] re-run with --apply to extend canonical_gf.csv")
        return

    canon = pd.read_csv(CANON, low_memory=False)
    before_hi = float(canon.wavelength_air_A.max())
    start = 1 + max(int(x.split('_')[1]) for x in canon.line_id if str(x).startswith('gf_'))

    allfe = pd.concat([fe1, fe[fe.ion != 1].reset_index(drop=True)], ignore_index=True)
    for c in ('adjudicated', 'lab_log_gf', 'lab_sigma_dex', 'lab_source'):
        if c not in allfe:
            allfe[c] = np.nan
    allfe['adjudicated'] = allfe.adjudicated.fillna(False).astype(bool)
    allfe = allfe.sort_values('wavelength_air_A').reset_index(drop=True)

    doi = {s: CITATIONS[s][1] for s in set(allfe.lab_source.dropna())}
    new = pd.DataFrame({
        'line_id': [f'gf_{start + i:06d}' for i in range(len(allfe))],
        'key_z': allfe.key_z, 'ion': allfe.ion.astype(float),
        'species': allfe.element + ' ' + allfe.ion.astype(str).map(_ROMAN),
        'wavelength_air_A': allfe.wavelength_air_A.round(3),
        'excitation_potential_eV': allfe.excitation_potential_eV.round(4),
        'hfs_n_components': allfe.hfs_n_components,
        'log_gf': np.where(allfe.adjudicated, allfe.lab_log_gf, allfe.gf_linelist_vald),
        'loggf_reference': np.where(
            allfe.adjudicated, 'PRIMARY LAB ' + allfe.lab_source.astype(str),
            allfe.reference_code.replace('', 'VALD3')),
        'nist_grade': None,
        'seed_source': 'hband_vald_12976_21400(RYA-1047)',
        'adjudication_status': np.where(allfe.adjudicated, 'lab_rya1047', 'single_source'),
        'gf_synth_ges': np.nan, 'gf_regions_vald': np.nan,
        'gf_linelist_vald': allfe.gf_linelist_vald,
        'in_synth': False, 'in_regions': False, 'in_linelist': True,
        'delta_synth_minus_ew': np.nan,
        'lab_source_tag': allfe.lab_source,
        'gf_sigma_dex': allfe.lab_sigma_dex,
        'gf_source_doi': allfe.lab_source.map(doi),
        'is_diagnostic': False,
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
