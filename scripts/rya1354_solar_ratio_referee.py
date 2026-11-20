"""RYA-1354 v69 (A1): the observed Sun as an abundance-free referee of RELATIVE gf (same-lower-level ratios).

Weak-line limit: W/lambda ~ gf * lambda * N_lower(J) / kappa_c(lambda), with N_lower(J)/g_J ~ exp(-E_J/kT). For lines that share
a lower TERM (same configuration + LS term; J partners differ by < 100 cm-1, Boltzmann applied at T = 5000 K), the abundance, the
term population and the atmosphere cancel in the ratio. kappa_c(lambda) does not cancel, so a group is restricted to
|dlambda|/lambda <= 10 %.

Per group (>= 2 lines) and per gf source X: r_i = [log(W_i/lambda_i) - log(lambda_i) + theta*E_i] - loggf_i(X); centring r over the
group removes the unknown constant (abundance x population x atmosphere). The scatter of centred r across all groups is the
source's score; lower = its RELATIVE gf matches the Sun better.

Null (declared): lines whose groups have every source agreeing (all |dgf| <= 0.03 after centring) measure the method's own noise
(EW, blends, kappa_c, saturation). Discrimination is claimed only for lines where the sources' centred gf differ by > 2 x null.
Saturation: primary = weak limit with log(W/lambda) <= -5.0; sensitivity = Milne-Eddington Gaussian-core correction and the
strict cut log(W/lambda) <= -5.3. REFEREE ONLY: no gf is changed.

A1b (v69, POST-HOC, disclosed on RYA-1354 before its run): the declared weak-limit run has no power (null 0.25 dex) because
solar Mg I lines are saturated far below log(W/lambda) = -5 (LTE MARCS COG slope 0.67 at 37 mA for 6318, 0.43 at 80 mA for 7659)
and Gaussian EWs miss the damping wings. Mode 'cog': EW from a Voigt fit; each line's EW is inverted through its own LTE COG
(iSpec/SPECTRUM, MARCS 5772/4.44/0.0, vmic 1.0; SPECTRUM default damping where the line list has none) to an astrophysical
loggf; lines with COG slope < 0.3 are dropped; groups are then centred as above. Single-line features only (J-tied excluded).
Usage: rya1354_solar_ratio_referee.py MODEL_DIR [--mode=cog]
"""
import os
import csv
import json
import math
import statistics as st
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

import rya1354_candidate_spectrum_check as C

ROOT = Path(__file__).resolve().parents[1]
THETA = 5040.0 / 5000.0          # eV^-1 (log10 Boltzmann at 5000 K)
CM_TO_EV = 1.239842e-4
WEAK, WEAK_STRICT = -5.0, -5.3
DLAM = 0.10
AGREE, D_MAX = 0.03, 0.90
SOURCES = ('model', 'linelist', 'abinitio', 'tuned', 'hybrid', 'pure_LS', 'lab')


def gauss_ew(a, w):
    """Gaussian + constant continuum fit to the observed line; EW in A (continuum-normalised)."""
    vel = (a.lam - w) / w * C.C_KMS
    m0 = C.measure(a, w, 0.1)
    if 'obs_depth' not in m0:
        return None, m0
    v0, h = m0['v_min_kms'], max(m0['hwhm_kms'], 1.0)
    win = np.abs(vel - v0) <= 2.5 * h
    v, f = vel[win], a.flux[win]
    if len(v) < 6:
        return None, m0
    try:
        p, _ = curve_fit(lambda x, c, d, x0, s: c - d * np.exp(-0.5 * ((x - x0) / s) ** 2), v, f,
                         p0=(1.0, m0['obs_depth'], v0, h / 1.1774), maxfev=4000)
    except RuntimeError:
        return None, m0
    c, d, x0, s = p
    if not (0 < d < c and 0.5 < abs(s) < 15 and abs(x0 - v0) < 2):
        return None, m0
    ew = d / c * abs(s) * math.sqrt(2 * math.pi) * w / C.C_KMS
    return dict(ew_A=ew, depth_fit=d / c, sigma_kms=abs(s), v_fit=x0), m0


def me_saturation(depth):
    """log10(W_linear / W) for a Gaussian-opacity Milne-Eddington line of central depth `depth` (d_max = D_MAX)."""
    if depth >= D_MAX * 0.98:
        return None
    eta = depth / (D_MAX - depth)
    x = np.linspace(-6, 6, 4001)
    phi = np.exp(-x ** 2)
    w = np.trapezoid(D_MAX * eta * phi / (1 + eta * phi), x)
    return math.log10(D_MAX * eta * math.sqrt(math.pi) / w)


ISP = Path('/Users/ryanschmitt/codex/engines/ispec_src')
ISPEC_LL = [ROOT / 'data/linelists/ispec_ir_9200_13000/atomic_lines.tsv', ROOT / 'data/linelists/ispec_h_15007_17494/atomic_lines.tsv',
            ROOT / 'data/linelists/ispec_k_19452_24846/atomic_lines.tsv',
            ISP / 'input/linelists/transitions/GESv6_atom_hfs_iso.420_920nm/atomic_lines.tsv']
COG_OFF = np.array([-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5])
MIN_COG_SLOPE = 0.3
VMIC = 1.0


def voigt_ew(a, w, half=25.0):
    """Voigt + linear continuum fit within +-half km/s; returns EW (A), Gaussian sigma, Lorentz gamma (km/s), rms."""
    from scipy.special import voigt_profile
    v = (a.lam - w) / w * C.C_KMS
    m = np.abs(v) <= half
    v, f = v[m], a.flux[m]
    if len(v) < 10:
        return None
    i = int(np.argmin(np.where(np.abs(v) < C.CORE_KMS, f, 9)))
    mod = lambda x, c, s1, A, x0, sg, gm: (c + s1 * x) - A * voigt_profile(x - x0, abs(sg), abs(gm))
    try:
        p, _ = curve_fit(mod, v, f, p0=(1, 0, (1 - f[i]) * 4, v[i], 2.5, 0.5), maxfev=20000)
    except RuntimeError:
        return None
    c, s1, A, x0, sg, gm = p
    return dict(ew_A=A / c * w / C.C_KMS, sigma_kms=abs(sg), gamma_kms=abs(gm), v_fit=x0,
                rms=float(np.sqrt(np.mean((f - mod(v, *p)) ** 2))))


def solar_cog(waves, vmic=VMIC):
    """LTE COG (EW in A at loggf_ispec + COG_OFF) for each wavelength found as 'Mg 1' in the iSpec line lists."""
    import ispec
    rows = []
    for f in ISPEC_LL:
        ll = ispec.read_atomic_linelist(str(f), wave_base=0, wave_top=1e9)
        rows.append(ll[ll['element'] == 'Mg 1'])
    rows = np.hstack(rows)
    idx = {}
    for w in waves:
        d = np.abs(rows['wave_A'] - w)
        i = int(np.argmin(d))
        if d[i] < 0.05:
            idx[w] = i
    ll = rows[list(idx.values())].copy()
    ll['spectrum_support'] = 'T'      # the codex IR lists are Turbospectrum-flagged; SPECTRUM handles them (default damping)
    base = ll['loggf'].copy()
    pack = ispec.load_modeled_layers_pack(str(ISP / 'input/atmospheres/MARCS.GES'))
    atm = ispec.interpolate_atmosphere_layers(pack, {'teff': 5772., 'logg': 4.44, 'MH': 0., 'alpha': 0.}, code='spectrum')
    iso = ispec.read_isotope_data(str(ISP / 'input/isotopes/SPECTRUM.lst'))
    ab = ispec.read_solar_abundances(str(ISP / 'input/abundances/Asplund.2009/stdatom.dat'))
    tmp = '/private/tmp/claude-501/ispec_rya1354'
    os.makedirs(tmp, exist_ok=True)
    ew = np.zeros((len(ll), len(COG_OFF)))
    for k, o in enumerate(COG_OFF):
        ll['loggf'] = base + o
        ew[:, k] = ispec.calculate_theoretical_ew_and_depth(atm, 5772., 4.44, 0., 0., ll, iso, ab, microturbulence_vel=vmic,
                                                            verbose=0, tmp_dir=tmp)['theoretical_ew'] / 1e3   # mA -> A
    return {w: dict(base=float(base[j]), ew=ew[j]) for j, w in enumerate(idx)}


def cog_invert(c, ew_obs):
    """Astrophysical loggf (at the model abundance) and the local COG slope; None outside the curve or where it is flat."""
    e = np.log10(np.maximum(c['ew'], 1e-9))
    if not np.all(np.diff(e) > 0) or not e[0] <= math.log10(ew_obs) <= e[-1]:
        return None, None
    x = float(np.interp(math.log10(ew_obs), e, COG_OFF))
    k = min(int(np.searchsorted(COG_OFF, x)), len(COG_OFF) - 1)
    slope = (e[k] - e[k - 1]) / (COG_OFF[k] - COG_OFF[k - 1])
    return c['base'] + x, float(slope)


def term(label):
    cfg, t = label.split(' ')
    return cfg + ' ' + t.rstrip('0123456789')


def main(model_dir, mode='weak'):
    md = ROOT / model_dir
    fill = [f for f in csv.DictReader(open(md / 'solar_fill.csv')) if f['status'] == 'graded' and float(f['max_depth']) >= 0.01]
    lines = {(r['lower'], r['upper']): r for r in csv.DictReader(open(md / 'lines.csv'))}
    elow = {r['lower']: float(r['E_lower']) for r in lines.values()}
    atl = C.load_atlases()
    meas = []
    for f in fill:
        lo, up = f['lower'], f['upper']
        if '{LSmix}' in lo + up or lo not in elow:
            continue
        w = float(f['air_A'])
        a = C.pick(atl, w)
        if a is None or a.name == 'KP1984' and w < 3200:
            continue
        g, m0 = gauss_ew(a, w)
        fails = set(filter(None, (m0.get('fails') or '').split('+')))
        if g is None or fails - {'shallow'}:
            continue
        r = lines.get((lo, up), {})
        gf = dict(model=float(f['model_loggf']), linelist=float(f['linelist_loggf']))
        for s in ('abinitio', 'tuned', 'hybrid'):
            if r.get(f'loggf_so_{s}'):
                gf[s] = float(r[f'loggf_so_{s}'])
        if r.get('loggf_pure_LS'):
            gf['pure_LS'] = float(r['loggf_pure_LS'])
        if r.get('lab_loggf'):
            gf['lab'] = float(r['lab_loggf'])
        if mode == 'cog':
            if '{' in up:
                continue
            vg = voigt_ew(a, w)
            if vg is None or vg['sigma_kms'] > 5.0 or vg['rms'] > 0.01 or abs(vg['v_fit']) > C.MAX_DV:
                continue
            g = dict(g, ew_A=vg['ew_A'])
        lw = math.log10(g['ew_A'] / w)
        meas.append(dict(air_A=w, lower=lo, upper=up, term=term(lo), E_lower_cm=elow[lo], atlas=a.name, ew_mA=round(1e3 * g['ew_A'], 2),
                         log_w_over_l=round(lw, 3), depth_fit=round(g['depth_fit'], 4), sat_me=me_saturation(g['depth_fit']),
                         obs_term=lw - math.log10(w) + THETA * elow[lo] * CM_TO_EV, gf=gf))
    if mode == 'cog':
        cg = solar_cog([m['air_A'] for m in meas])
        keep = []
        for m in meas:
            if m['air_A'] not in cg:
                continue
            ga, sl = cog_invert(cg[m['air_A']], m['ew_mA'] / 1e3)
            if ga is None or sl < MIN_COG_SLOPE:
                continue
            # obs_term in the 'gf' unit: astrophysical loggf (abundance, population, kappa_c all in the COG)
            m.update(obs_term=ga, cog_slope=round(sl, 3), astro_loggf=round(ga, 3))
            keep.append(m)
        meas = keep
    # groups: same lower term, clustered in wavelength so that |dl|/l <= DLAM from the cluster's first line
    by = {}
    for m in sorted(meas, key=lambda x: x['air_A']):
        by.setdefault(m['term'], []).append(m)
    groups = []
    for t, ms in by.items():
        cur = []
        for m in ms:
            if cur and (m['air_A'] - cur[0]['air_A']) / cur[0]['air_A'] > DLAM:
                groups.append(cur)
                cur = []
            cur.append(m)
        groups.append(cur)

    def score(cut, sat):
        res = {s: [] for s in SOURCES}
        spread, rows = [], []
        for gi, g in enumerate(groups):
            g = [m for m in g if (mode == 'cog' or m['log_w_over_l'] <= cut) and (not sat or m['sat_me'] is not None)]
            if len(g) < 2:
                continue
            obs = [m['obs_term'] + (m['sat_me'] if sat else 0) for m in g]
            for s in SOURCES:
                if not all(s in m['gf'] for m in g):
                    continue
                r = [o - m['gf'][s] for o, m in zip(obs, g)]
                mu = st.mean(r)
                res[s] += [x - mu for x in r]
            # source disagreement per line (centred gf), over sources present for the whole group
            have = [s for s in SOURCES if all(s in m['gf'] for m in g)]
            cg = {s: [m['gf'][s] - st.mean(x['gf'][s] for x in g) for m in g] for s in have}
            for i, m in enumerate(g):
                sp_ = max(cg[s][i] for s in have) - min(cg[s][i] for s in have)
                rows.append(dict(group=gi, m=m, obs=obs[i] - st.mean(obs), spread=sp_,
                                 cgf={s: cg[s][i] for s in have}, all_agree=None))
                spread.append(sp_)
        # null: groups where every line agrees across sources
        gagree = {}
        for x in rows:
            gagree.setdefault(x['group'], []).append(x['spread'] <= AGREE)
        null_r = [x['obs'] - x['cgf']['model'] for x in rows if all(gagree[x['group']])]
        null = math.sqrt(st.mean([v * v for v in null_r])) if len(null_r) >= 3 else None
        disc = [x for x in rows if null is not None and x['spread'] > 2 * null]
        rms = {s: (round(math.sqrt(st.mean([v * v for v in r])), 4), len(r)) for s, r in res.items() if len(r) >= 2}
        dscore = {}
        for s in SOURCES:
            v = [x['obs'] - x['cgf'][s] for x in disc if s in x['cgf']]
            if v:
                dscore[s] = dict(n=len(v), rms=round(math.sqrt(st.mean([u * u for u in v])), 4))
        return dict(n_groups=len({x['group'] for x in rows}), n_lines=len(rows), rms_by_source=rms,
                    null_rms=round(null, 4) if null else None, n_null_lines=len(null_r), n_discriminating=len(disc),
                    discriminating_rms_by_source=dscore), rows, disc

    out = md / ('rya1354_solar_ratio_referee' + ('' if mode == 'weak' else '_' + mode))
    out.mkdir(exist_ok=True)
    summ = dict(schema='rya1354.solar-ratio-referee.v1', model_dir=str(model_dir), n_measured=len(meas), weak_cut=WEAK, dlam=DLAM,
                theta=THETA, me_dmax=D_MAX)
    runs = (('primary_weak', WEAK, False), ('strict_weak', WEAK_STRICT, False), ('me_saturation', WEAK, True)) if mode == 'weak' \
        else (('primary_weak', None, False),)
    summ['mode'] = mode
    for name, cut, sat in runs:
        s, rows, disc = score(cut, sat)
        summ[name] = s
        if name == 'primary_weak':
            with open(out / 'lines.csv', 'w', newline='') as fh:
                cols = ['group', 'air_A', 'lower', 'upper', 'atlas', 'ew_mA', 'log_w_over_l', 'depth_fit', 'obs_centred', 'spread', 'discriminating'] + \
                       [f'cgf_{s}' for s in SOURCES]
                wr = csv.DictWriter(fh, fieldnames=cols, lineterminator='\n')
                wr.writeheader()
                ids = {id(x) for x in disc}
                for x in rows:
                    m = x['m']
                    wr.writerow(dict(group=x['group'], air_A=m['air_A'], lower=m['lower'], upper=m['upper'], atlas=m['atlas'], ew_mA=m['ew_mA'],
                                     log_w_over_l=m['log_w_over_l'], depth_fit=m['depth_fit'], obs_centred=round(x['obs'], 4),
                                     spread=round(x['spread'], 4), discriminating=id(x) in ids,
                                     **{f'cgf_{s}': round(x['cgf'][s], 4) if s in x['cgf'] else '' for s in SOURCES}))
    (out / 'summary.json').write_text(json.dumps(summ, indent=1) + '\n')
    return summ


if __name__ == '__main__':
    md_ = [a for a in sys.argv[1:] if not a.startswith('--')][0]
    mode_ = next((a.split('=')[1] for a in sys.argv[1:] if a.startswith('--mode=')), 'weak')
    print(json.dumps(main(md_, mode_), indent=1))
