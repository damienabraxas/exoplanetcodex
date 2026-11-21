#!/usr/bin/env python3
"""RYA-1354 v70c -- solar strong-line WING test of van der Waals broadening (report only).

Declared on RYA-1354 (comment bdf5d990, 2026-11-21) before any run.

Method (Anstee & O'Mara 1995; Bergemann 2017): gf fixed to lab, A(X) fixed to AAG21, so the
only free line parameter is d = Delta log gamma6 relative to ABO, applied by scaling the
packed ABO sigma with alpha kept (wing depth ~ gf N gamma6, so A and gf errors propagate 1:1
into d and are carried as error terms). One analytic continuum scalar per window.

Fit mask (fixed in the declaration):
  * core excluded where the ABO synthetic flux < 0.6
  * points where the target line's own depth (full - blend-only synthesis) < 0.01 dropped
  * points where the blend-only synthesis (targets removed) is deeper than 0.02 dropped
  * one 3-MAD residual clip at d = 0

Error: sigma(d) = chi2 curvature (chi2 scaled to reduced 1) (+) dd/dA*0.03 (+) dd/dgf*sigma_gf
(+) |MARCS - ATLAS9|. dd/dA from a refit at A+0.05; dd/dgf = dd/dA (same N*gf product).

Runs single-process, niced, on the Mac. Usage:
  PYTHONPATH=~/codex/engines/ispec_src CODEX_KP_ATLAS=... ~/codex/venv/bin/python \
      scripts/rya1354_vdw_wings.py [--bench] [--only NAME ...] [--atm MARCS.GES|ATLAS9.Castelli ...]
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

OUT = ROOT / 'data/audit/rya1354/mg1_vdw_wings_v70c'
D_GRID = np.round(np.arange(-0.6, 0.81, 0.1), 2)
DA = 0.05
SIG_A = 0.03
CORE_FLUX, OWN_DEPTH, BLEND_DEPTH, CLIP = 0.6, 0.01, 0.02, 3.0
# attempt 2 (POST HOC after attempt-1 gate failure, declared on RYA-1354 reply c3ecd398 before rerun):
OWN_OVER_BLEND, MIN_PIX = 3.0, 30

# AAG21 (Asplund, Amarsi & Grevesse 2021, Table 2)
A_REF = {'Mg': 7.55, 'Na': 6.22, 'Ca': 6.30}
CODE = {'Mg': 12, 'Na': 11, 'Ca': 20}

# name, element, window (A, air), target line centres (A, air), group (multiplet sharing one ABO
# sigma), role, sigma_gf (dex, lab; Pehlivan 2017 exp 0.02-0.04, NIST AA 0.01, NIST A 0.03)
WINDOWS = [
    dict(name='NaD', el='Na', win=(5884.0, 5901.0), lines=[5889.951, 5895.924], group='Na 3s-3p',
         role='control', sig_gf=0.01, atlas='IAG'),
    dict(name='Ca6102', el='Ca', win=(6100.0, 6105.5), lines=[6102.723], group='Ca 4s4p-4s5s',
         role='control', sig_gf=0.03, atlas='IAG'),
    dict(name='Ca6122', el='Ca', win=(6119.5, 6125.0), lines=[6122.217], group='Ca 4s4p-4s5s',
         role='control', sig_gf=0.03, atlas='IAG'),
    dict(name='Ca6162', el='Ca', win=(6159.5, 6165.0), lines=[6162.173], group='Ca 4s4p-4s5s',
         role='control', sig_gf=0.03, atlas='IAG'),
    dict(name='Mg5167', el='Mg', win=(5164.5, 5170.0), lines=[5167.322], group='Mg b',
         role='T1', sig_gf=0.04, atlas='IAG'),
    dict(name='Mg5172', el='Mg', win=(5169.0, 5177.0), lines=[5172.684], group='Mg b',
         role='T1', sig_gf=0.04, atlas='IAG'),
    dict(name='Mg5183', el='Mg', win=(5178.0, 5189.0), lines=[5183.604], group='Mg b',
         role='T1', sig_gf=0.04, atlas='IAG'),
    dict(name='Mg5528', el='Mg', win=(5526.0, 5531.0), lines=[5528.405], group='Mg 5528',
         role='T1', sig_gf=0.02, atlas='IAG'),
    dict(name='Mg8806', el='Mg', win=(8802.5, 8811.0), lines=[8806.757], group='Mg 8806',
         role='T1', sig_gf=0.03, atlas='IAG'),
    dict(name='Mg4703', el='Mg', win=(4701.0, 4705.0), lines=[4702.991], group='Mg 4703',
         role='T3', sig_gf=0.05, atlas='KP1984'),
    # T4: H band. The IR list carries no vdW (TS fdamp 2.5 = Unsoeld x 2.5); targets get OUR v70 ABO
    # (sigma.alpha) and lab gf (Pehlivan 2017 exp, mg1_lab_loggf_pehlivan2017.csv) where it exists,
    # else Bergemann 2017 Table 2 (PR17 theo, list value kept, sigma 0.05).
    dict(name='Mg15024', el='Mg', win=(15018.0, 15031.0), lines=[15024.997], group='Mg 4s-4p',
         role='T4', sig_gf=0.03, atlas='CRIRES+_H', ir=True, fdamp_set=954.227, gf_set={15024.997: 0.334}),
    dict(name='Mg15040', el='Mg', win=(15035.0, 15044.0), lines=[15040.246], group='Mg 4s-4p',
         role='T4', sig_gf=0.04, atlas='CRIRES+_H', ir=True, fdamp_set=954.227, gf_set={15040.246: 0.113}),
    dict(name='Mg15047', el='Mg', win=(15044.0, 15052.0), lines=[15047.714], group='Mg 4s-4p',
         role='T4', sig_gf=0.04, atlas='CRIRES+_H', ir=True, fdamp_set=954.227, gf_set={15047.714: -0.364}),
    dict(name='Mg15740', el='Mg', win=(15735.0, 15745.0), lines=[15740.706], group='Mg 4p-4d',
         role='T4', sig_gf=0.05, atlas='CRIRES+_H', ir=True, fdamp_set=1640.280),
    dict(name='Mg15748', el='Mg', win=(15744.5, 15754.0), lines=[15748.886, 15748.988], group='Mg 4p-4d',
         role='T4', sig_gf=0.05, atlas='CRIRES+_H', ir=True, fdamp_set=1640.280),
    dict(name='Mg15765', el='Mg', win=(15760.0, 15771.0), lines=[15765.747, 15765.839], group='Mg 4p-4d',
         role='T4', sig_gf=0.05, atlas='CRIRES+_H', ir=True, fdamp_set=1640.280),
]
IR_LIST = ROOT / 'data/linelists/ispec_h_15007_17494/atomic_lines.tsv'
MIN_MTRANS = 0.9      # existing telluric-pixel rule (rya1354_candidate_spectrum_check.MIN_MTRANS)


def load_atlas(kind):
    from rya1354_candidate_spectrum_check import load_atlases
    global _ATL
    if '_ATL' not in globals():
        _ATL = {a.name: a for a in load_atlases()}
    return _ATL['IAG_Baker2020' if kind == 'IAG' else kind]


class Synth:
    def __init__(self, atm_grid):
        import ispec
        from pipeline import abundances_derive as ad
        from config.constants import get_star_params
        self.ispec, self.ad = ispec, ad
        sp = get_star_params('solar')
        self.sp = sp
        self.atm = ad._load_atmosphere(sp['teff'], sp['logg'], sp['feh_ref'], sp['xi'], model_grid=atm_grid)
        self.ll, self.iso, _ = ad._load_synth_resources()
        self.ll_ir = None
        self.sab = ispec.read_solar_abundances(ad._ISPEC_SOLAR_ABUND_FILE)
        self.tmp = os.environ.get('V70C_TMP', '/private/tmp/claude-501/v70c_ispec')

    def window_ll(self, w):
        lo, hi = w['win']
        if w.get('ir') and self.ll_ir is None:
            self.ll_ir, _, _ = self.ad._load_synth_resources(linelist_file=str(IR_LIST), apply_canonical_gf=False)
        src = self.ll_ir if w.get('ir') else self.ll
        ll = src[(src['wave_A'] > lo - 3) & (src['wave_A'] < hi + 3)].copy()
        tgt = np.zeros(len(ll), bool)
        for c in w['lines']:
            tgt |= (np.abs(ll['wave_A'] - c) < 0.006) & (np.char.startswith(ll['element'].astype(str), w['el'] + ' 1'))
        if tgt.sum() < len(w['lines']):
            raise SystemExit(f"{w['name']}: only {tgt.sum()} of {len(w['lines'])} target lines found in the list")
        for i in np.where(tgt)[0]:
            if 'fdamp_set' in w:
                ll['turbospectrum_fdamp'][i] = w['fdamp_set']
            for c, g in w.get('gf_set', {}).items():
                if abs(ll['wave_A'][i] - c) < 0.006:
                    ll['loggf'][i] = g
        return ll, tgt

    def run(self, w, ll, A, wave_A):
        return self.ad._synth_flux_at_abund(
            wave_A / 10.0, self.atm, self.sp['teff'], self.sp['logg'], self.sp['feh_ref'], self.sp['xi'],
            ll, self.iso, self.sab, w['el'], CODE[w['el']], A, R=100000 if w.get('ir') else 700000,
            macroturbulence=self.sp['vmac'], vsini=self.sp['vsini'], tmp_dir=self.tmp)


def scale_fdamp(ll, tgt, d):
    out = ll.copy()
    for i in np.where(tgt)[0]:
        f = float(out['turbospectrum_fdamp'][i])
        if f <= 1.0:
            raise SystemExit(f"target at {out['wave_A'][i]} has no packed ABO (fdamp {f})")
        sig, alpha = np.floor(f), f - np.floor(f)
        out['turbospectrum_fdamp'][i] = round(sig * 10 ** d) + alpha
    return out


def fit_window(S, w, wave_A, obs, log, good=None):
    ll, tgt = S.window_ll(w)
    A0 = A_REF[w['el']]
    fd = [float(x) for x in ll['turbospectrum_fdamp'][tgt]]
    gf = [float(x) for x in ll['loggf'][tgt]]
    full0 = S.run(w, scale_fdamp(ll, tgt, 0.0), A0, wave_A)
    blend = S.run(w, ll[~tgt], A0, wave_A)
    own = blend - full0
    mask = (full0 >= CORE_FLUX) & (own >= OWN_DEPTH) & ((1 - blend) <= BLEND_DEPTH) & np.isfinite(obs)
    if good is not None:
        mask &= good
    mask &= own >= OWN_OVER_BLEND * (1 - blend)
    c0 = np.sum(obs[mask] * full0[mask]) / np.sum(full0[mask] ** 2)
    r = obs[mask] - c0 * full0[mask]
    mad = 1.4826 * np.median(np.abs(r - np.median(r)))
    keep = np.zeros_like(mask)
    keep[np.where(mask)[0][np.abs(r - np.median(r)) <= CLIP * mad]] = True
    log(f"  {w['name']}: {len(ll)} lines, targets fdamp {fd} gf {gf}; mask {mask.sum()} -> {keep.sum()} after clip")
    if keep.sum() < MIN_PIX:
        return dict(name=w['name'], el=w['el'], group=w['group'], role=w['role'], verdict='UNMEASURABLE',
                    n_mask=int(keep.sum()))

    def chi2_curve(A):
        out = []
        for d in D_GRID:
            f = S.run(w, scale_fdamp(ll, tgt, d), A, wave_A)
            c = np.sum(obs[keep] * f[keep]) / np.sum(f[keep] ** 2)
            out.append(float(np.sum((obs[keep] - c * f[keep]) ** 2)))
        return np.array(out)

    def best(ch):
        i = int(np.clip(np.argmin(ch), 1, len(D_GRID) - 2))
        x, y = D_GRID[i - 1:i + 2], ch[i - 1:i + 2]
        a, b, _ = np.polyfit(x, y, 2)
        edge = np.argmin(ch) in (0, len(D_GRID) - 1)
        return (-b / (2 * a) if a > 0 else float(D_GRID[np.argmin(ch)])), a, edge

    ch = chi2_curve(A0)
    d0, a2, edge = best(ch)
    dof = max(keep.sum() - 2, 1)
    s2 = ch.min() / dof                      # residual variance -> reduced chi2 = 1
    sig_curv = float(np.sqrt(s2 / a2)) if a2 > 0 else float('nan')
    chA = chi2_curve(A0 + DA)
    dA, _, _ = best(chA)
    dd_dA = (dA - d0) / DA
    return dict(name=w['name'], el=w['el'], group=w['group'], role=w['role'], atlas=w['atlas'],
                lines=w['lines'], fdamp_abo=fd, loggf=gf, A=A0, n_fit=int(keep.sum()),
                d_grid=D_GRID.tolist(), chi2=ch.tolist(), chi2_Aplus=chA.tolist(),
                d_fit=float(d0), at_grid_edge=bool(edge), sig_curv=sig_curv, dd_dA=float(dd_dA),
                sig_A=abs(dd_dA) * SIG_A, sig_gf=abs(dd_dA) * w['sig_gf'], rms=float(np.sqrt(s2)),
                cont_scale=float(c0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bench', action='store_true')
    ap.add_argument('--only', nargs='*')
    ap.add_argument('--atm', nargs='*', default=['MARCS.GES', 'ATLAS9.Castelli'])
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    logf = open(OUT / ('bench.log' if args.bench else 'run.log'), 'a')

    def log(s):
        print(s, flush=True)
        logf.write(s + '\n'); logf.flush()

    wins = [w for w in WINDOWS if not args.only or w['name'] in args.only]
    for atm in args.atm:
        S = Synth(atm)
        log(f"== atmosphere {atm}  {time.strftime('%H:%M:%S')}")
        for w in wins:
            a = load_atlas(w['atlas'])
            lo, hi = w['win']
            sel = (a.lam >= lo) & (a.lam <= hi)
            wave_A, obs = a.lam[sel], a.flux[sel]
            good = (a.mtrans[sel] >= MIN_MTRANS) if a.mtrans is not None else None
            if args.bench:
                ll, tgt = S.window_ll(w)
                t = time.time(); S.run(w, ll, A_REF[w['el']], wave_A)
                log(f"  bench {w['name']}: {len(wave_A)} px, {len(ll)} lines, {time.time()-t:.1f} s/synth")
                continue
            t = time.time()
            res = fit_window(S, w, wave_A, obs, log, good)
            res['atm'] = atm
            res['seconds'] = round(time.time() - t, 1)
            log(f"  {w['name']} [{atm}]: d_fit {res.get('d_fit', float('nan')):+.3f} "
                f"curv {res.get('sig_curv', float('nan')):.3f} dd/dA {res.get('dd_dA', float('nan')):+.2f} "
                f"n {res.get('n_fit')} ({res['seconds']} s)")
            (OUT / f"fit_{w['name']}_{atm.split('.')[0]}.json").write_text(json.dumps(res, indent=1))
        if args.bench:
            break


if __name__ == '__main__':
    main()
