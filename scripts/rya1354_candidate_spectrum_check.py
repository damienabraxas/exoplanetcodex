"""Check Theoretical Codex line candidates against the OBSERVED solar spectrum (RYA-1354). Mac-side.

The candidate tiers (rya1354_theory_calibration_{oi,kurucz}.py) used VALD-synthesis depths and blend flags. This
step looks at each candidate in a telluric-corrected (or telluric-free) normalised solar atlas:

  lambda < 5001 A          Kitt Peak flux atlas 1984 (Kurucz+), air, pseudo-residual flux (no telluric correction;
                           pipeline.telluric_policy bands flagged)
  5001-11083 A             IAG telluric-free atlas, Baker+2020 (vacuum wavenumber -> air)
  9800-10796 A             CRIRES+ Vesta Y, normalised (RYA-794) [preferred over IAG inside its range]
  11160-13490 / 14943-17585 / 19452-24846 A
                           CRIRES+ Vesta J / H / K, molecfit-corrected, rest frame (RYA-1214/1219/1232);
                           min_mtrans (molecfit transmission) flags residual telluric risk

Per line, inside +-3 km/s of the catalogue wavelength: observed depth (1 - min flux), velocity of the minimum, core
symmetry (flux difference at +-HWHM of the observed core, relative to the depth), the deepest OTHER minimum within
+-25 km/s, and the telluric flag. Verdict CLEAN when: depth >= 0.05, |v_min| <= 2 km/s, asymmetry <= 0.25,
no deeper-than-half neighbour inside 2-4 x HWHM, telluric OK. Observed/predicted depth is reported (depth_anomaly
outside 0.5-2) but is NOT a cleanliness test. Criteria validated on a positive control: the 40 Asplund 2021 Fe I
lines + O I 777 pass offset/asymmetry/telluric 43/43; 'close_neighbour' flags 3/40 (7%).

Status: CANDIDATE screening; not a measurement, not a grade.
"""
import csv
import glob
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline.telluric_policy import in_telluric_band  # noqa: E402
from pipeline.wavelength_util import vac_to_air  # noqa: E402

C_KMS = 299792.458
OUT = ROOT / 'data/audit/rya1354/candidate_spectrum_check'
KP_DIR = Path(os.environ.get('CODEX_KP_ATLAS', Path.home() / 'codex/spectra/Kitt Peak Flux Atlas'))
IAG_BAKER = Path(os.environ.get('CODEX_SRV_SOLAR_ROOT', Path.home() / 'codex/sirius_local_srv/solar_reference')) / 'iag_baker2020/iag_telfree_solaratlas.fits'
CRIRES = {'Y': ROOT / 'data/results/rya794/vesta_crires_plus_Y_9800_10796_normalized.csv',
          'J': ROOT / 'data/results/rya1214_crires_jk/solar_crires_plus_j_rya1219_rest.csv',
          'H': ROOT / 'data/results/rya1214_crires_jk/solar_crires_plus_h_rya1232_rest.csv',
          'K': ROOT / 'data/results/rya1214_crires_jk/solar_crires_plus_k_rya1219_rest.csv'}
CORE_KMS, NEIGH_KMS, MAX_DV, MAX_ASYM = 3.0, 25.0, 2.0, 0.25
MIN_MTRANS = 0.9


class Atlas:
    def __init__(self, name, lam, flux, mtrans=None, telluric_corrected=True):
        o = np.argsort(lam)
        self.name, self.lam, self.flux = name, lam[o], flux[o]
        self.mtrans = mtrans[o] if mtrans is not None else None
        self.telluric_corrected = telluric_corrected

    def covers(self, w, pad):
        return self.lam[0] + pad < w < self.lam[-1] - pad


def load_atlases():
    a = []
    for k, p in CRIRES.items():
        d = np.genfromtxt(p, delimiter=',', names=True)
        a.append(Atlas(f'CRIRES+_{k}', d['wavelength_air_A'], d['flux_normalized'],
                       d['min_mtrans'] if 'min_mtrans' in d.dtype.names else None))
    from astropy.io import fits
    t = fits.open(IAG_BAKER)[1].data
    good = (t['s'] > 0) & np.isfinite(t['s']) & (t['err'] < 0.05)
    a.append(Atlas('IAG_Baker2020', vac_to_air(1e8 / t['v'][good].astype(float)), t['s'][good].astype(float)))
    lam, flux = [], []
    for f in sorted(glob.glob(str(KP_DIR / 'lm*'))):
        d = np.loadtxt(f, usecols=(0, 1))
        lam.append(d[:, 0] * 10)
        flux.append(d[:, 1])
    lam, flux = np.concatenate(lam), np.concatenate(flux)
    a.append(Atlas('KP1984', lam[lam < 5001], flux[lam < 5001], telluric_corrected=False))
    return a          # order = preference


def pick(atlases, w):
    for a in atlases:
        if a.covers(w, w * NEIGH_KMS / C_KMS * 1.2):
            return a
    return None


def measure(a, w, pred_depth):
    vel = (a.lam - w) / w * C_KMS
    win = np.abs(vel) <= NEIGH_KMS
    v, f = vel[win], a.flux[win]
    core = np.abs(v) <= CORE_KMS
    if core.sum() < 3:
        return dict(atlas=a.name, verdict='NO_PIXELS')
    i = np.argmin(np.where(core, f, np.inf))
    fmin, vmin = f[i], v[i]
    depth = 1 - fmin
    # observed half-depth half-width: walk out from the minimum to flux 1 - depth/2
    half = 1 - depth / 2
    lo = i
    while lo > 0 and f[lo] < half:
        lo -= 1
    hi = i
    while hi < len(f) - 1 and f[hi] < half:
        hi += 1
    hwhm = max((v[hi] - v[lo]) / 2, 0.5)
    fl = np.interp(vmin - hwhm, v, f)
    fr = np.interp(vmin + hwhm, v, f)
    asym = abs(fl - fr) / depth if depth > 0 else np.inf
    # deepest other minimum outside the line core (+-2 HWHM) but within the window
    out = np.abs(v - vmin) > 2 * hwhm
    neigh = 1 - f[out].min() if out.any() else 0.0
    near = (np.abs(v - vmin) > 2 * hwhm) & (np.abs(v - vmin) <= 4 * hwhm)
    near_depth = 1 - f[near].min() if near.any() else 0.0
    if a.mtrans is not None:
        mt = float(a.mtrans[win][core].min())
        tell_ok = mt >= MIN_MTRANS
    else:
        mt = None
        tell_ok = a.telluric_corrected or not in_telluric_band(w)
    ratio = depth / pred_depth if pred_depth > 0 else np.inf
    # obs/pred depth is REPORTED, not a cleanliness test: the predicted depths are VALD LTE estimates with the
    # uncorrected theory gf (O I IR is strongly NLTE: 777 observed ~1.8x deeper; Kurucz gf run ~0.1 dex strong).
    fails = [n for n, ok in (('shallow', depth >= 0.05), ('offset', abs(vmin) <= MAX_DV), ('asymmetric', asym <= MAX_ASYM),
                             ('close_neighbour', near_depth <= 0.5 * depth),
                             ('telluric', tell_ok)) if not ok]
    return dict(atlas=a.name, obs_depth=round(float(depth), 4), v_min_kms=round(float(vmin), 2), hwhm_kms=round(float(hwhm), 2),
                asymmetry=round(float(asym), 3), obs_over_pred=round(float(ratio), 3),
                depth_anomaly=bool(not 0.5 <= ratio <= 2), near_neighbour_depth=round(float(near_depth), 4),
                deepest_in_window=round(float(neigh), 4), min_mtrans=mt, verdict='CLEAN' if not fails else 'FLAGGED',
                fails='+'.join(fails))


def candidate_list():
    rows = []
    for r in csv.DictReader(open(ROOT / 'data/audit/rya1354/theory_calibration_kurucz/candidates.csv')):
        if r['tier'].startswith('T2'):
            rows.append(dict(species=r['species'], air_A=float(r['air_A']), pred_depth=float(r['solar_depth']), tier=r['tier'],
                             theory_loggf=r['k_loggf_bias_corrected'], sigma68_dex=r['sigma68_dex'], theory='Kurucz (bias-corrected)',
                             transition=f"{r['lower']} - {r['upper']}"))
    for r in csv.DictReader(open(ROOT / 'data/audit/rya1354/theory_calibration_oi/candidates.csv')):
        if r['tier'][:2] in ('T1', 'T2') and float(r['lambda_air_nm']) >= 300:
            rows.append(dict(species='O I', air_A=round(float(r['lambda_air_nm']) * 10, 3), pred_depth=float(r['solar_depth']),
                             tier=r['tier'], theory_loggf=r['li_loggf_B'], sigma68_dex='0.064' if r['li_class'] in ('A', 'B') else '',
                             theory='Li 2023 MCDHF', transition=r['transition']))
    return rows


def features(cands, merge_kms=6.0):
    """Components of one species within merge_kms of each other are ONE observable feature (multiplets such as
    O I 615.8/844.6/926 nm, fine-structure clusters): checked once, at the component predicted deepest, against
    the combined predicted depth 1 - prod(1 - d_i)."""
    out = []
    for sp in sorted({c['species'] for c in cands}):
        cs = sorted((c for c in cands if c['species'] == sp), key=lambda c: c['air_A'])
        group = [cs[0]]
        for c in cs[1:] + [None]:
            if c is not None and (c['air_A'] - group[-1]['air_A']) / c['air_A'] * C_KMS <= merge_kms:
                group.append(c)
                continue
            lead = max(group, key=lambda g: g['pred_depth'])
            out.append({**lead, 'n_components': len(group), 'components_A': ' '.join(f"{g['air_A']:.3f}" for g in group),
                        'pred_depth': round(1 - float(np.prod([1 - g['pred_depth'] for g in group])), 4)})
            if c is not None:
                group = [c]
    return out


def line_set_index():
    """(species, wavelength) of every line in every Codex line set (data/reference/line_sets/*.csv, not the
    *_graded derivatives): a candidate already in a set is a KNOWN abundance line, not a new one."""
    idx = []
    for p in sorted((ROOT / 'data/reference/line_sets').glob('*.csv')):
        if '_graded' in p.name or p.name == 'REGISTRY.csv':
            continue
        for r in csv.DictReader(open(p)):
            if r.get('species') and r.get('wavelength_air_A'):
                idx.append((r['species'], float(r['wavelength_air_A']), p.stem))
    return idx


def main():
    atlases = load_atlases()
    sets = line_set_index()
    out = []
    for c in features(candidate_list()):
        a = pick(atlases, c['air_A'])
        m = measure(a, c['air_A'], c['pred_depth']) if a else dict(atlas=None, verdict='NO_ATLAS_COVERAGE')
        known = sorted({n for sp, w, n in sets if sp == c['species'] and abs(w - c['air_A']) <= 0.05})
        out.append({**c, **m, 'in_codex_line_sets': ' '.join(known), 'new_line': not known})
    OUT.mkdir(parents=True, exist_ok=True)
    keys = list(dict.fromkeys(k for r in out for k in r))
    with open(OUT / 'checked.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys, lineterminator='\n')
        w.writeheader()
        w.writerows(sorted(out, key=lambda r: (r['verdict'] != 'CLEAN', r['species'], r['air_A'])))
    from collections import Counter, defaultdict
    by = defaultdict(Counter)
    for r in out:
        b = 'IR>1um' if r['air_A'] >= 10000 else 'NIR700-1000' if r['air_A'] >= 7000 else 'VIS400-700' if r['air_A'] >= 4000 else 'NUV300-400'
        by[f"{r['species']}|{b}"][r['verdict']] += 1
    fails = Counter(x for r in out if r.get('fails') for x in r['fails'].split('+'))
    clean = [r for r in out if r['verdict'] == 'CLEAN']
    summary = dict(schema='rya1354.candidate-spectrum-check.v1', status='SCREENING_NOT_A_MEASUREMENT', n=len(out),
                   verdicts=dict(Counter(r['verdict'] for r in out)),
                   clean_new_lines=sum(r['new_line'] for r in clean), clean_already_in_line_sets=sum(not r['new_line'] for r in clean),
                   clean_new_by_species={k: v for k, v in sorted(Counter(r['species'] for r in clean if r['new_line']).items())},
                   clean_new_ir_ge_7000A={k: v for k, v in sorted(Counter(r['species'] for r in clean if r['new_line'] and r['air_A'] >= 7000).items())}, fail_reasons=dict(fails.most_common()),
                   by_species_band={k: dict(v) for k, v in sorted(by.items())},
                   atlases={a.name: [round(float(a.lam[0]), 1), round(float(a.lam[-1]), 1)] for a in atlases},
                   criteria=dict(core_kms=CORE_KMS, window_kms=NEIGH_KMS, max_dv_kms=MAX_DV, max_asym=MAX_ASYM, min_mtrans=MIN_MTRANS,
                                 depth_anomaly_outside=[0.5, 2], min_depth=0.05, close_neighbour='<= 0.5 x depth within 2-4 HWHM'))
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    return summary


if __name__ == '__main__':
    print(json.dumps(main(), indent=1))
