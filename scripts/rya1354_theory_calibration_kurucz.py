"""Theory-vs-lab calibration of Kurucz gf for EVERY species we hold primary lab gf for, plus Theoretical Codex
line candidates per species (RYA-1354 / RYA-1352). Generalises rya1354_theory_calibration_fe.py (same matching,
same quality classes, same gap rule); Fe I is re-run here so all species sit in one table.

Theory: data/reference/kurucz/gfZZII.pos.gz (gzip of kurucz.harvard.edu/atoms/ZZII/gfZZII.pos, 2026-11-18;
every row a Kurucz calculation, K07..K18) and data/reference/kurucz_fe1/gf2600.pos.gz.
Lab: data/reference/lab_gf_tables/lab_gf_lines.csv (VizieR + transcribed primaries), fe_gf_lab, si_gf_lab (DH23).

Matching: lab rows with both level energies in cm-1 -> both energies within 0.25 cm-1; lab rows with energies
in eV (4 decimals) -> lower level within 1.5 cm-1 AND Kurucz air wavelength within 0.03 A. Unique match only.

Status: CANDIDATE, HOLD. Nothing here writes a grade or sigma into a product.
"""
import csv
import json
import statistics as st
from collections import defaultdict
from pathlib import Path

import rya1354_theory_calibration_fe as F

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/audit/rya1354/theory_calibration_kurucz'
CM_PER_EV = F.CM_PER_EV
SPECIES = {'Si I': '1400', 'Ca I': '2000', 'Ti I': '2200', 'Ti II': '2201', 'V I': '2300', 'V II': '2301',
           'Cr I': '2400', 'Cr II': '2401', 'Mn I': '2500', 'Mn II': '2501', 'Fe I': '2600', 'Fe II': '2601',
           'Co I': '2700', 'Co II': '2701', 'Ni I': '2800'}
WL_TOL_A = 0.03
ELO_TOL_EV_ROWS = 1.5


def kurucz_path(code):
    return ROOT / ('data/reference/kurucz_fe1/gf2600.pos.gz' if code == '2600' else f'data/reference/kurucz/gf{code}.pos.gz')


def lab_rows(species):
    out = []
    for r in csv.DictReader(open(ROOT / 'data/reference/lab_gf_tables/lab_gf_lines.csv')):
        if r['species'] == species and r['loggf'] and r['loggf_sigma_dex']:
            out.append(dict(source=r['paper'], air_A=float(r['wavelength_air_A']), elo=float(r['ep_eV']) * CM_PER_EV, eup=None,
                            loggf=float(r['loggf']), sig=float(r['loggf_sigma_dex'])))
    files = {'Fe I': F.LAB, 'Fe II': ('data/reference/fe_gf_lab/fe2_lab_loggf_dh19.csv',),
             'Si I': ('data/reference/si_gf_lab/si1_lab_loggf_dh23.csv',)}.get(species, ())
    for p in files:
        for r in csv.DictReader(open(ROOT / p)):
            out.append(dict(source=r['source'], air_A=float(r['wavelength_air_A']), elo=float(r['elo_cm1']), eup=float(r['eup_cm1']),
                            loggf=float(r['loggf']), sig=float(r['e_loggf_dex'])))
    return out


def air_A(k):
    return k['wl_nm'] * 10


def match(idx, rows_by_elo, l):
    if l['eup'] is not None:
        return F.find(idx, l['elo'], l['eup'])
    return [k for k in rows_by_elo.get(round(l['elo'] / 10), []) + rows_by_elo.get(round(l['elo'] / 10) - 1, [])
            + rows_by_elo.get(round(l['elo'] / 10) + 1, [])
            if abs(k['elo'] - l['elo']) <= ELO_TOL_EV_ROWS and abs(air_A(k) - l['air_A']) <= WL_TOL_A]


def calibrate(species, krows, idx):
    by_elo = defaultdict(list)
    for k in krows:
        by_elo[round(k['elo'] / 10)].append(k)
    matched, n_lab, n_unmatched = [], 0, 0
    for l in lab_rows(species):
        n_lab += 1
        h = match(idx, by_elo, l)
        if len(h) != 1:
            n_unmatched += 1
            continue
        k = h[0]
        matched.append(dict(species=species, source=l['source'], air_A=l['air_A'], band=F.band(l['air_A']), elo=round(k['elo'], 3),
                            eup=round(k['eup'], 3), lab_loggf=l['loggf'], sig=l['sig'], k_loggf=k['loggf'], k_ref=k['ref'],
                            d=round(k['loggf'] - l['loggf'], 3), lower=k['lab_lo'], upper=k['lab_up'],
                            k_class=F.k_class(k['lab_lo'], k['lab_up'], k['loggf'])))
    return matched, dict(n_lab=n_lab, n_matched=len(matched), n_unmatched=n_unmatched)


def class_table(rows, min_n=8):
    out = {}
    for c in sorted({r['k_class'] for r in rows}):
        d = [r['d'] for r in rows if r['k_class'] == c]
        if len(d) >= min_n:
            out[c] = dict(n=len(d), bias=round(st.median(d), 3), sigma68=round(F.half_width(d, 0.68), 3),
                          sigma95=round(F.half_width(d, 0.95), 3),
                          frac_off_gt_0p5=round(sum(abs(x - st.median(d)) > 0.5 for x in d) / len(d), 2))
    return out


def candidates(species, krows, idx, lab, cls, own):
    """Same rule as Fe: VALD gf IS Kurucz's (gap), K class with a calibrated sigma68 <= 0.15 dex, 0.1 <= depth <= 0.8,
    not blend-flagged, ground-observable; not a primary-lab line."""
    sym, ion = species.split()
    lab_up = [(l['elo'], l['air_A']) for l in lab]
    out = []
    for r in csv.DictReader(open(ROOT / 'data/linelists/linelist_solar.csv')):
        if r['element'] != sym or r['ion'] != ion or r['blend_flag'] == 'True':
            continue
        depth = float(r['central_depth'] or 0)
        air = float(r['wavelength_air_A'])
        if depth < 0.05:
            continue
        elo = float(r['excitation_potential_eV']) * CM_PER_EV
        k = [x for x in krows_near(krows, elo, air)]
        if len(k) != 1:
            continue
        k = k[0]
        if any(abs(k['elo'] - a) < ELO_TOL_EV_ROWS and abs(air - b) < WL_TOL_A for a, b in lab_up):
            continue
        kc = F.k_class(k['lab_lo'], k['lab_up'], k['loggf'])
        cs = cls.get(kc)
        gap = abs(float(r['log_gf']) - k['loggf']) < 0.006
        good = bool(cs) and cs['sigma68'] <= 0.15
        tier = ('T0_VALD_HAS_NON_KURUCZ_GF' if not gap else
                ('T2_THEORETICAL_CODEX_CANDIDATE' if kc in own else 'T2P_CANDIDATE_POOLED_SIGMA_ONLY')
                if good and 0.1 <= depth <= 0.8 and air >= 3000 else
                'T3_GAP_LOWER_CLASS_OR_DEPTH' if cs and air >= 3000 else 'T4_GAP_UNCALIBRATED_OR_SPACE_UV')
        out.append(dict(species=species, air_A=round(air, 3), band=F.band(air), chi_eV=float(r['excitation_potential_eV']),
                        solar_depth=depth, vald_loggf=float(r['log_gf']), k_loggf=k['loggf'], k_class=kc,
                        k_loggf_bias_corrected=round(k['loggf'] - cs['bias'], 3) if cs else None,
                        sigma68_dex=cs['sigma68'] if cs else None, sigma95_dex=cs['sigma95'] if cs else None,
                        sigma_source=('species' if kc in own else 'pooled') if cs else None,
                        lower=k['lab_lo'], upper=k['lab_up'], tier=tier))
    return out


_NEAR = {}


def krows_near(krows, elo, air):
    key = id(krows)
    if key not in _NEAR:
        b = defaultdict(list)
        for k in krows:
            b[round(air_A(k))].append(k)
        _NEAR[key] = b
    b = _NEAR[key]
    return [k for d in (-1, 0, 1) for k in b.get(round(air) + d, [])
            if abs(k['elo'] - elo) < 1.0 and abs(air_A(k) - air) < 0.05]


def main():
    all_rows, per, cand = [], {}, []
    loaded = {}
    for sp, code in SPECIES.items():
        krows, idx = F.kurucz(kurucz_path(code))
        loaded[sp] = (krows, idx)
        rows, counts = calibrate(sp, krows, idx)
        all_rows += rows
        d = [r['d'] for r in rows]
        per[sp] = dict(counts, kurucz_lines=len(krows), kurucz_refs=sorted({k['ref'] for k in krows}),
                       **(dict(median=round(st.median(d), 3), robust_sigma=round(F.robust_sigma(d), 3), rms=round(F.rms(d), 3),
                               frac_within_0p1=round(sum(abs(x) <= 0.1 for x in d) / len(d), 2), classes=class_table(rows))
                          if len(d) >= 8 else {}))
    pooled = class_table(all_rows)
    for sp, (krows, idx) in loaded.items():
        # per-species class table where it has >= 8 lines in that class, else the pooled class
        cls = dict(pooled)
        cls.update(per[sp].get('classes', {}))
        cand += candidates(sp, krows, idx, lab_rows(sp), cls, set(per[sp].get('classes', {})))
    OUT.mkdir(parents=True, exist_ok=True)
    for name, rows in (('kurucz_vs_lab', all_rows), ('candidates', sorted(cand, key=lambda c: (c['tier'], c['species'], -c['solar_depth'])))):
        with open(OUT / f'{name}.csv', 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
            w.writeheader()
            w.writerows(rows)
    t2 = defaultdict(lambda: defaultdict(int))
    for c in cand:
        if c['tier'].startswith('T2_'):
            t2[c['species']][c['band']] += 1
    summary = dict(schema='rya1354.theory-calibration-kurucz.v1', status='CANDIDATE_HOLD_NOT_ADOPTED', per_species=per,
                   pooled_classes=pooled, n_lab_matched=len(all_rows),
                   candidates=dict(n=len(cand), by_tier={t: sum(c['tier'] == t for c in cand) for t in sorted({c['tier'] for c in cand})},
                                   t2_by_species_band={s: dict(sorted(v.items())) for s, v in sorted(t2.items())}))
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    return summary


if __name__ == '__main__':
    s = main()
    print(f"{'species':7} {'lab':>5} {'match':>5} {'med':>6} {'rob':>6} {'rms':>6} {'<0.1':>5}  best class")
    for sp, v in s['per_species'].items():
        best = v.get('classes', {}).get('dS0|strong', {})
        print(f"{sp:7} {v['n_lab']:5} {v['n_matched']:5} {v.get('median', ''):>6} {v.get('robust_sigma', ''):>6} {v.get('rms', ''):>6} "
              f"{v.get('frac_within_0p1', ''):>5}  {best}")
    print('POOLED', json.dumps(s['pooled_classes'], indent=1))
    print(json.dumps(s['candidates'], indent=1))
