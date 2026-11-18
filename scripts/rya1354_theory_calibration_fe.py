"""Calibrate THEORETICAL Fe I gf (Kurucz 2017 semi-empirical, K17) line by line against modern primary
laboratory gf, and rank Fe I lines with no lab gf as Theoretical Codex candidates (RYA-1354 / RYA-1352).
Mac-side; no solver.

Fe I is the opposite case to O I: hundreds of primary lab gf with per-line sigma exist (Ruffoni 2013/2014,
Den Hartog 2014, Belmonte 2017; data/reference/fe_gf_lab), so the theory error is MEASURED per line instead of
inferred from a handful of lifetimes.

Theory: data/reference/kurucz_fe1/gf2600.pos (kurucz.harvard.edu/atoms/2600, calculation of 22 Oct 2014,
reference code K17 on every row: lines between experimentally known levels, gf from Kurucz's least-squares
fitted (semi-empirical) eigenvectors). Matched to the lab by BOTH level energies.

Status: CANDIDATE, HOLD. Nothing here writes a grade or sigma into a product.
"""
import csv
import gzip
import json
import math
import re
import statistics as st
from collections import defaultdict
from pathlib import Path

from rya1354_theory_calibration_oi import fit_sigma_th, rms, robust_sigma

ROOT = Path(__file__).resolve().parents[1]
KURUCZ = ROOT / 'data/reference/kurucz_fe1/gf2600.pos.gz'   # gzip of kurucz.harvard.edu/atoms/2600/gf2600.pos
OUT = ROOT / 'data/audit/rya1354/theory_calibration_fe'
CM_PER_EV = 8065.54429
E_TOL = 0.25            # cm-1: lab papers use Nave 1994 energies, K17 the 2014 NIST revision
LAB = ('data/reference/fe_gf_lab/fe1_lab_loggf.csv', 'data/reference/fe_gf_lab/fe1_lab_loggf_ruffoni2013.csv')


def band(air_A):
    return ('UV<300' if air_A < 3000 else 'NUV300-400' if air_A < 4000 else 'VIS400-700' if air_A < 7000
            else 'NIR700-1000' if air_A < 10000 else 'IR>1um')


def strength(loggf):
    return 'loggf<-3' if loggf < -3 else '-3..-2' if loggf < -2 else '-2..-1' if loggf < -1 else '>=-1'


def eup_bin(e_cm):
    return 'Eup<40k' if e_cm < 40000 else '40-50k' if e_cm < 50000 else '50-55k' if e_cm < 55000 else 'Eup>=55k'


def multiplicity(label):
    m = re.search(r'(\d)[SPDFGHIK]\s*$', label.strip())
    return int(m.group(1)) if m else None


def k_class(lower, upper, k17_loggf):
    """Quality class for a K17 gf, from the two indicators the lab comparison shows matter (Kurucz gives no
    gauge or cancellation indicator): spin change between the LS labels (intercombination lines carry the
    mixing errors, up to 6 dex) and predicted strength."""
    a, b = multiplicity(lower), multiplicity(upper)
    spin = 'jK' if a is None or b is None else 'dS0' if a == b else 'dS'
    st_ = 'strong' if k17_loggf >= -1 else 'mid' if k17_loggf >= -2.5 else 'weak'
    return f'{spin}|{st_}'


def half_width(d, q):
    """Half-width of the central q-fraction of residuals about their median (tail-honest sigma)."""
    med = st.median(d)
    a = sorted(abs(x - med) for x in d)
    return a[min(len(a) - 1, int(math.ceil(q * len(a))) - 1)]


def kurucz():
    rows = []
    for l in gzip.open(KURUCZ, 'rt'):
        e1, e2 = float(l[24:36]), float(l[52:64])
        if l[98:101] != 'K17':
            raise ValueError(f'Non-K17 row in gf2600.pos: {l[:40]!r}')
        rows.append(dict(wl_nm=float(l[:11]), loggf=float(l[11:18]), elo=min(e1, e2), eup=max(e1, e2),
                         lab_lo=(l[42:52] if e1 <= e2 else l[70:80]).strip(), lab_up=(l[70:80] if e1 <= e2 else l[42:52]).strip()))
    idx = defaultdict(list)
    for r in rows:
        idx[(round(r['elo']), round(r['eup']))].append(r)
    return rows, idx


def find(idx, elo, eup, tol=E_TOL):
    hits = [r for di in (-1, 0, 1) for dk in (-1, 0, 1) for r in idx.get((round(elo) + di, round(eup) + dk), [])
            if abs(r['elo'] - elo) <= tol and abs(r['eup'] - eup) <= tol]
    return hits


def lab_lines():
    out = []
    for p in LAB:
        for r in csv.DictReader(open(ROOT / p)):
            out.append(dict(source=r['source'], air_A=float(r['wavelength_air_A']), elo=float(r['elo_cm1']), eup=float(r['eup_cm1']),
                            loggf=float(r['loggf']), sig=float(r['e_loggf_dex'])))
    return out


def summarise(rows):
    d = [r['d'] for r in rows]
    if len(d) < 3:
        return dict(n=len(d))
    med = st.median(d)
    return dict(n=len(d), median=round(med, 3), robust_sigma=round(robust_sigma(d), 3), rms=round(rms(d), 3),
                frac_within_0p1=round(sum(abs(x) <= 0.1 for x in d) / len(d), 2),
                extra_theory_sigma_about_median=fit_sigma_th([x - med for x in d], [r['sig'] for r in rows]),
                extra_theory_sigma_zero_bias=fit_sigma_th(d, [r['sig'] for r in rows]))


def calibrate(idx):
    matched, unmatched = [], []
    for l in lab_lines():
        h = find(idx, l['elo'], l['eup'])
        if len(h) != 1:
            unmatched.append(dict(l, n_hits=len(h)))
            continue
        k = h[0]
        matched.append(dict(source=l['source'], air_A=l['air_A'], band=band(l['air_A']), elo=l['elo'], eup=l['eup'],
                            lab_loggf=l['loggf'], sig=l['sig'], k17_loggf=k['loggf'], d=round(k['loggf'] - l['loggf'], 3),
                            strength=strength(l['loggf']), eup_bin=eup_bin(l['eup']), lower=k['lab_lo'], upper=k['lab_up'],
                            k_class=k_class(k['lab_lo'], k['lab_up'], k['loggf'])))
    res = dict(n_lab=len(matched) + len(unmatched), n_matched=len(matched), unmatched=unmatched, all=summarise(matched))
    for key in ('band', 'strength', 'eup_bin', 'source', 'k_class'):
        res[f'by_{key}'] = {g: summarise([r for r in matched if r[key] == g]) for g in sorted({r[key] for r in matched})}
    res['by_band_x_strength'] = {f'{b}|{s}': summarise([r for r in matched if r['band'] == b and r['strength'] == s])
                                 for b in sorted({r['band'] for r in matched}) for s in ('>=-1', '-2..-1', '-3..-2', 'loggf<-3')}
    # Calibrated sigma per K class: a K17 gf used in place of a lab gf is corrected by the class median bias;
    # sigma68 / sigma95 are the central 68% / 95% half-widths of the corrected residuals (heavy tails kept).
    res['class_sigma'] = {}
    for c in sorted({r['k_class'] for r in matched}):
        d = [r['d'] for r in matched if r['k_class'] == c]
        if len(d) >= 8:
            res['class_sigma'][c] = dict(n=len(d), bias=round(st.median(d), 3), sigma68=round(half_width(d, 0.68), 3),
                                         sigma95=round(half_width(d, 0.95), 3), frac_off_gt_0p5=round(sum(abs(x - st.median(d)) > 0.5 for x in d) / len(d), 2))
    res['class_sigma_ir_gt_1um'] = {c: summarise([r for r in matched if r['k_class'] == c and r['band'] == 'IR>1um'])
                                    for c in ('dS0|strong', 'jK|strong', 'dS0|mid', 'dS|strong')}
    return matched, res


def candidates(idx, cal, lab):
    """Fe I rows of our solar list with NO primary lab gf and no Asplund grading, measurable and unblended."""
    graded = {round(float(r['wavelength_air_A']), 1) for r in csv.DictReader(open(ROOT / 'data/reference/line_sets/fe_agss21_FeI.csv'))
              if r.get('wavelength_air_A')}
    lab_keys = [(l['elo'], l['eup']) for l in lab]
    out = []
    for r in csv.DictReader(open(ROOT / 'data/linelists/linelist_solar.csv')):
        if r['element'] != 'Fe' or r['ion'] not in ('1', 'I'):
            continue
        depth = float(r['central_depth'] or 0)
        if depth < 0.05 or r['blend_flag'] == 'True':
            continue
        air = float(r['wavelength_air_A'])
        elo = float(r['excitation_potential_eV']) * CM_PER_EV
        k = [x for x in find(idx, elo, elo + 1e8 / (air * 1.00027) if air > 2000 else elo + 1e8 / air, tol=60)
             if abs(x['elo'] - elo) < 1.0 and abs(x['wl_nm'] * 10 - air) < 0.05]
        if len(k) != 1:
            continue
        k = k[0]
        if any(abs(k['elo'] - a) < E_TOL and abs(k['eup'] - b) < E_TOL for a, b in lab_keys):
            continue                                   # has a primary lab gf: not a gap
        if round(air, 1) in graded:
            continue                                   # Asplund-graded line
        b = band(air)
        kc = k_class(k['lab_lo'], k['lab_up'], k['loggf'])
        cs = cal['class_sigma'].get(kc)
        out.append(dict(air_A=round(air, 3), band=b, chi_eV=float(r['excitation_potential_eV']), solar_depth=depth,
                        vald_loggf=float(r['log_gf']), k17_loggf=k['loggf'], d_vald_minus_k17=round(float(r['log_gf']) - k['loggf'], 3),
                        lower=k['lab_lo'], upper=k['lab_up'], k_class=kc,
                        k17_loggf_bias_corrected=round(k['loggf'] - cs['bias'], 3) if cs else None,
                        sigma68_dex=cs['sigma68'] if cs else None, sigma95_dex=cs['sigma95'] if cs else None,
                        vald_is_k17=abs(float(r['log_gf']) - k['loggf']) < 0.006, ground_observable=air >= 3000))
    for c in out:
        # GAP = VALD itself carries the Kurucz value (no lab gf anywhere VALD knows of).
        if not c['vald_is_k17']:
            c['tier'] = 'T0_VALD_HAS_NON_K17_GF'
        elif c['k_class'] in ('dS0|strong', 'jK|strong') and 0.1 <= c['solar_depth'] <= 0.8 and c['ground_observable']:
            c['tier'] = 'T2_THEORETICAL_CODEX_CANDIDATE'
        elif c['k_class'] in ('dS0|strong', 'jK|strong', 'dS0|mid') and c['ground_observable']:
            c['tier'] = 'T3_GAP_LOWER_CLASS_OR_DEPTH'
        else:
            c['tier'] = 'T4_GAP_POOR_CLASS_OR_SPACE_UV'
    return sorted(out, key=lambda c: (c['tier'], c['band'], -c['solar_depth']))


def main():
    _, idx = kurucz()
    matched, cal = calibrate(idx)
    lab = lab_lines()
    cand = candidates(idx, cal, lab)
    OUT.mkdir(parents=True, exist_ok=True)
    for name, rows in (('k17_vs_lab', matched), ('candidates', cand)):
        with open(OUT / f'{name}.csv', 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
            w.writeheader()
            w.writerows(rows)
    tiers = defaultdict(lambda: defaultdict(int))
    for c in cand:
        tiers[c['tier']][c['band']] += 1
    summary = dict(schema='rya1354.theory-calibration-fe.v1', status='CANDIDATE_HOLD_NOT_ADOPTED',
                   theory='Kurucz K17 Fe I (gf2600.pos, 22 Oct 2014 calculation)',
                   lab='Ruffoni 2013 (H band), Ruffoni 2014, Den Hartog 2014, Belmonte 2017 (primary, per-line sigma)',
                   calibration=cal, candidates=dict(n=len(cand), by_tier_and_band={t: dict(sorted(v.items())) for t, v in sorted(tiers.items())},
                                                    selection='solar depth >= 0.05, not blend-flagged, unique K17 match, no primary lab gf, not Asplund-graded; '
                                                              'T2 = VALD gf IS K17 (gap), K class dS0|strong or jK|strong, 0.1 <= depth <= 0.8, lambda >= 300 nm'))
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    return summary


if __name__ == '__main__':
    s = main()
    c = s['calibration']
    print('matched', c['n_matched'], 'of', c['n_lab'], 'unmatched', [(u['source'], u['air_A'], u['n_hits']) for u in c['unmatched']][:12])
    print('ALL', c['all'])
    print('class_sigma', json.dumps(c['class_sigma'], indent=1))
    for k in ('by_band', 'by_k_class'):
        print('--', k)
        for g, v in c[k].items():
            print(f'   {g:22s} {v}')
    print('IR class', json.dumps(c['class_sigma_ir_gt_1um']))
    print(json.dumps(s['candidates']['by_tier_and_band'], indent=1))
