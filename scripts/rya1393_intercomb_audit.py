"""RYA-1393: spin-changing (intercombination) rows of linelist_solar.csv vs the observed Sun and our atomic model.

For each species with a two-valence model run (Mg I first), every line-list feature our model identifies as spin-changing
(lower and upper LS multiplicity differ) with line-list predicted depth >= MIN_DEPTH is measured in the observed solar
atlases (CRIRES+ Y/J/H/K, IAG, Kitt Peak; rya1354_candidate_spectrum_check). An absence is a hypothesis: every band
carries POSITIVE CONTROLS (spin-allowed lines of the same species, our model agreeing with the line list to <= 0.1 dex,
predicted depth >= 0.2), and a band whose controls fail cannot refute anything.

Verdict per feature (observed depth d_obs vs line-list predicted depth d_ll, model-predicted gf vs line-list gf):
  REFUTED            d_obs < 0.25 d_ll and d_obs < 0.05, clean window, band controls pass
  CONFIRMED          0.5 <= d_obs/d_ll <= 2 and the absorption is centred (|dv| <= 2 km/s)
  BLENDED_UNRESOLVED anything else (crowded / offset / asymmetric / partial)
  NO_ATLAS           no atlas coverage
Usage: rya1393_intercomb_audit.py MODEL_DIR   (e.g. data/audit/rya1354/mg1_model_v67)
"""
import csv
import json
import sys
from pathlib import Path

import rya1354_candidate_spectrum_check as C

ROOT = Path(__file__).resolve().parents[1]
MIN_DEPTH = 0.05


def band(w):
    for name, lo, hi in (('UV', 0, 3200), ('blue', 3200, 5001), ('VIS', 5001, 9800), ('Y', 9800, 11100), ('J', 11100, 13500),
                         ('H', 14500, 18200), ('K', 19500, 25000)):
        if lo <= w < hi:
            return name
    return 'other'


def main(model_dir):
    md = ROOT / model_dir
    fill = list(csv.DictReader(open(md / 'solar_fill.csv')))
    atl = C.load_atlases()

    def meas(w, pd):
        a = C.pick(atl, w)
        return C.measure(a, w, pd) if a else dict(verdict='NO_ATLAS', atlas=None)
    spin, ctrl = [], []
    def mult(label):
        return label.split(' ')[-1][0]
    for f in fill:
        if f['model_loggf'] in ('', None) or '{LSmix}' in f['lower'] + f['upper'] or f['lower'].startswith('?'):
            continue      # configuration-level (S not defined) and unassigned features cannot be spin-classified
        sc = mult(f['lower']) != mult(f['upper'])     # from the term labels: also covers J-tied features ('3D{1,2}')
        pd = float(f['max_depth'])
        d = float(f['d_model_minus_linelist'])
        if sc and pd >= MIN_DEPTH:
            spin.append(dict(f, band=band(float(f['air_A']))))
        elif not sc and pd >= 0.2 and abs(d) <= 0.1:
            ctrl.append(dict(f, band=band(float(f['air_A']))))
    # positive controls per band
    ctl_out = []
    for c in ctrl:
        m = meas(float(c['air_A']), float(c['max_depth']))
        ok = m.get('obs_depth') is not None and 0.5 <= m['obs_depth'] / float(c['max_depth']) <= 2 and abs(m.get('v_min_kms', 9)) <= 3
        ctl_out.append(dict(air_A=c['air_A'], lower=c['lower'], upper=c['upper'], band=c['band'], pred_depth=c['max_depth'],
                            obs_depth=m.get('obs_depth'), atlas=m.get('atlas'), dv=m.get('v_min_kms'), control_pass=ok))
    band_ok = {}
    for c in ctl_out:
        band_ok.setdefault(c['band'], []).append(c['control_pass'])
    band_pass = {b: (sum(v) >= max(1, len(v) // 2)) for b, v in band_ok.items()}
    out = []
    for s in spin:
        w, pd = float(s['air_A']), float(s['max_depth'])
        m = meas(w, pd)
        od = m.get('obs_depth')
        if od is None:
            v = 'NO_ATLAS'
        elif not band_pass.get(s['band'], False):
            v = 'NO_CONTROL'
        elif od < 0.25 * pd and od < 0.05:
            v = 'REFUTED'
        elif 0.5 <= od / pd <= 2 and abs(m.get('v_min_kms', 9)) <= 2:
            v = 'CONFIRMED'
        else:
            v = 'BLENDED_UNRESOLVED'
        out.append(dict(air_A=s['air_A'], lower=s['lower'], upper=s['upper'], band=s['band'], linelist_loggf=s['linelist_loggf'],
                        model_loggf=s['model_loggf'], d_model_minus_linelist=s['d_model_minus_linelist'], pred_depth_linelist=pd,
                        obs_depth=od, atlas=m.get('atlas'), dv_kms=m.get('v_min_kms'), fails=m.get('fails'), verdict=v))
    od_ = md / 'rya1393_intercomb_audit'
    od_.mkdir(exist_ok=True)
    for name, rows in (('features.csv', out), ('controls.csv', ctl_out)):
        with open(od_ / name, 'w', newline='') as fh:
            w_ = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator='\n')
            w_.writeheader()
            w_.writerows(rows)
    summ = dict(schema='rya1393.intercomb-audit.v1', model_dir=str(model_dir), min_depth=MIN_DEPTH, n_spin_features=len(out),
                verdicts={k: sum(o['verdict'] == k for o in out) for k in sorted({o['verdict'] for o in out})},
                by_band={b: {k: sum(o['verdict'] == k and o['band'] == b for o in out) for k in sorted({o['verdict'] for o in out})}
                         for b in sorted({o['band'] for o in out})},
                controls={b: dict(n=len(v), passed=sum(v), band_usable=band_pass[b]) for b, v in band_ok.items()},
                refuted_and_model_weaker=sum(o['verdict'] == 'REFUTED' and float(o['d_model_minus_linelist']) < -0.5 for o in out))
    (od_ / 'summary.json').write_text(json.dumps(summ, indent=1) + '\n')
    return summ


if __name__ == '__main__':
    print(json.dumps(main(sys.argv[1]), indent=1))
