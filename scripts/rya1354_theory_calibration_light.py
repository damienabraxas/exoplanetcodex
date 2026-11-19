"""Theory-vs-lab for the light elements with ORCH products (RYA-1354). Na/Mg/Al carry almost no Kurucz-gap lines
(VALD gf = Kurucz on 4 measurable Na I and 4 Mg I lines), because their line-list gf are NIST-compiled ab initio.
The question there is the quality of THAT theory:
  Mg I  Pehlivan Rhodin+2017 MCHF vs their lab BF x lifetime (data/reference/mg_gf_lab/), split by whether the
        lifetime scale is experimental (independent) or their own calculation (not independent).
  Al I  Papoulia+2019 MCDHF, Kurucz 1995 and TOPbase vs Burheim+2023 lab (data/reference/al_gf_lab/).
Status: CANDIDATE, HOLD."""
import csv
import json
import statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/audit/rya1354/theory_calibration_light'


def extra_sigma(d, s):
    if sum(v * v / (si * si) for v, si in zip(d, s)) <= len(d):
        return 0.0
    lo, hi = 0.0, 2.0
    for _ in range(60):
        m = (lo + hi) / 2
        lo, hi = (m, hi) if sum(v * v / (si * si + m * m) for v, si in zip(d, s)) > len(d) else (lo, m)
    return round(lo, 3)


def stats(d, s):
    md = st.median(d)
    return dict(n=len(d), median=round(md, 3), rms=round((sum(v * v for v in d) / len(d)) ** .5, 3),
                sigma68_about_median=round(sorted(abs(v - md) for v in d)[min(len(d) - 1, int(.68 * len(d)))], 3),
                extra_theory_sigma=extra_sigma(d, s))


def rows(p):
    return list(csv.DictReader(l for l in open(p) if not l.startswith('#')))


def main():
    mg = rows(ROOT / 'data/reference/mg_gf_lab/mg1_lab_loggf_pehlivan2017.csv')
    al = rows(ROOT / 'data/reference/al_gf_lab/al1_lab_loggf.csv')
    out = {'Mg I MCHF (Pehlivan Rhodin 2017)': {}, 'Al I (vs Burheim 2023 lab)': {}}
    for name, sel in (('all', mg), ('lab_lifetime_only', [r for r in mg if r['tau_theory'] == 'False']),
                      ('theory_lifetime_scale', [r for r in mg if r['tau_theory'] == 'True'])):
        out['Mg I MCHF (Pehlivan Rhodin 2017)'][name] = stats([float(r['loggf_calc']) - float(r['loggf_exp']) for r in sel],
                                                              [float(r['e_loggf_dex']) for r in sel])
    for col, lab in (('loggf_papoulia19', 'Papoulia 2019 MCDHF'), ('loggf_kurucz95', 'Kurucz 1995'), ('loggf_topbase00', 'TOPbase')):
        sel = [r for r in al if r[col]]
        out['Al I (vs Burheim 2023 lab)'][lab] = stats([float(r[col]) - float(r['loggf']) for r in sel], [float(r['e_loggf_dex']) for r in sel])
    summary = dict(schema='rya1354.theory-calibration-light.v1', status='CANDIDATE_HOLD_NOT_ADOPTED', results=out,
                   kurucz_gap_lines_measurable={'Na I': 4, 'Mg I': 4},
                   note='Na/Mg/Al line-list gf are NIST-compiled ab initio, not Kurucz: the Kurucz error model does not apply.')
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    return summary


if __name__ == '__main__':
    print(json.dumps(main(), indent=1))
