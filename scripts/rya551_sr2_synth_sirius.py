#!/usr/bin/env python3
"""
RYA-551 — Sr II synthesis abundance measurement (RUNS ON SIRIUS).

Completes RYA-421 Step 1: measure A(Sr II) for 4077/4161/4215/4305 from the
observed solar spectrum via Turbospectrum LTE flux synthesis WITH blend handling
(the full VALD near-UV block in-window), NOT a naive Gaussian EW. Reuses the
RYA-534 babsma/bsyn wiring. Post-hoc NLTE delta from the committed
Sr_Bergemann2012_INSPECT grid (4077/4215 only; 4161/4305 = NLTE_unavailable).

Engine B (RYA-525) = TS LTE synthesis (this fit). Engine A (1D-NLTE) = Engine B
+ INSPECT delta for the covered doublet.

Two solar "arms":
  - HARPS  (primary, per RYA-421): the reflected-solar normalized spectrum.
  - IAG    (cross-check): the disk-integrated FTS solar flux atlas (Reiners 2016).
    ESPRESSO solar is NOT staged anywhere; IAG is the honest cross-source check.

gf/EP for the Sr II lines are the CANONICAL single-source values (linelist_solar);
the VALD-list gf is OVERWRITTEN so the measured A is on our gf scale.

Validate-don't-tune: report what the fit gives; never adjust toward Asplund or to
force Sr I = Sr II agreement.
"""
import json
import os
import re
import subprocess
import sys

import numpy as np

EXE   = "/mnt/codex-data/engines/Turbospectrum_NLTE/exec-gf"
MARCS = ("/mnt/codex-data/grids/model_atmospheres/marcs_standard_comp/marcs_standard_comp/"
         "p5750_g+4.5_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod")
VALD_DIR = "/mnt/codex-data/engines/TSFitPy/input_files/linelists/linelist_vald"
W = "/mnt/codex-data/codex/_solar_tmp/rya551_work"
HARPS = "/mnt/codex-data/codex/_solar_tmp/solar_normalized_harps.csv"
IAG   = "/mnt/codex-data/solar_reference/iag_reiners2016/spvis.dat.gz"
INSPECT = "/mnt/codex-data/codex/_solar_tmp/Sr_Bergemann2012_INSPECT.csv"

A_SUN = 2.83          # Asplund 2021 reference point (NOT a tuning target)
XI    = 1.0           # solar microturbulence (matches the t01 MARCS node)
VSINI = 1.8           # solar vsini (STAR_PARAMS)
CLIGHT = 299792.458

# Canonical Sr II lines: (wave_air, EP_eV, log_gf) from linelist_solar (single source).
LINES = {
    4077.709: dict(ep=0.0000, loggf=0.143,  vald='vald-3800-4300-for-grid.list',
                    fit_hw=1.6, feh_node_delta='4077.709'),
    4161.792: dict(ep=2.9403, loggf=-0.327, vald='vald-3800-4300-for-grid.list',
                    fit_hw=0.9, feh_node_delta=None),
    4215.519: dict(ep=0.0000, loggf=-0.173, vald='vald-3800-4300-for-grid.list',
                    fit_hw=1.6, feh_node_delta='4215.519'),
    4305.443: dict(ep=3.0397, loggf=-0.041, vald='vald-4300-4800-for-grid.list',
                    fit_hw=0.9, feh_node_delta=None),
}
LMARGIN = 8.0                              # synth half-window (A) for babsma/bsyn
A_GRID = np.round(np.arange(2.10, 3.36, 0.05), 3)   # A(Sr) trial grid
A_LO, A_HI = float(A_GRID.min()), float(A_GRID.max())


def sh(cmd, **kw):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True, **kw)


def inspect_delta(tag):
    """Solar-node (teff 5772, logg 4.4, feh 0) NLTE delta from the committed grid."""
    if tag is None:
        return None
    import csv
    with open(INSPECT) as f:
        for row in csv.DictReader(f):
            if (abs(float(row['wave_A']) - float(tag)) < 0.01 and
                    int(float(row['teff_K'])) == 5772 and abs(float(row['logg']) - 4.4) < 0.01
                    and abs(float(row['feh'])) < 1e-6):
                return float(row['delta_nlte'])
    return None


def patch_linelist(center, vald_name):
    """Copy the VALD near-UV block, overwrite the target Sr II line's log_gf with the
    canonical value. bsyn selects lines by LAMBDA window; blends stay verbatim."""
    src = os.path.join(VALD_DIR, vald_name)
    dst = f"{W}/ll_{center:.0f}.list"
    loggf = LINES[center]['loggf']
    out = []
    with open(src) as f:
        for ln in f:
            p = ln.split()
            if len(p) >= 3:
                try:
                    wl = float(p[0])
                except ValueError:
                    out.append(ln); continue
                if abs(wl - center) < 0.02 and "Sr II" in ln:
                    # columns: wave EP loggf ... ; replace the 3rd token in place
                    m = re.match(r"^(\s*\S+\s+\S+\s+)(\S+)(\s+.*)$", ln.rstrip("\n"))
                    if m:
                        newln = f"{m.group(1)}{loggf:.3f}{m.group(3)}\n"
                        out.append(newln)
                        continue
            out.append(ln)
    with open(dst, "w") as f:
        f.writelines(out)
    return dst


def babsma(center, opac, lmin, lmax):
    ctl = (f"'LAMBDA_MIN:'  '{lmin}'\n'LAMBDA_MAX:'  '{lmax}'\n'LAMBDA_STEP:' '0.005'\n"
           f"'MODELINPUT:' '{MARCS}'\n'MARCS-FILE:' '.true.'\n'MODELOPAC:' '{opac}'\n"
           f"'METALLICITY:'    '0.00'\n'ALPHA/Fe   :'    '0.00'\n'HELIUM     :'    '0.00'\n"
           f"'R-PROCESS  :'    '0.00'\n'S-PROCESS  :'    '0.00'\n'XIFIX:' 'T'\n{XI}\n")
    r = subprocess.run([f"{EXE}/babsma_lu"], input=ctl, capture_output=True, text=True, cwd=W)
    if not os.path.exists(opac):
        raise SystemExit(f"babsma failed @ {center}:\n{r.stdout[-1500:]}\n{r.stderr[-800:]}")


def bsyn(center, a_x, linelist, opac, lmin, lmax, tag):
    res = f"{W}/s_{center:.0f}_{tag}.spec"
    ctl = (f"'NLTE :'          '.false.'\n'LAMBDA_MIN:'     '{lmin}'\n'LAMBDA_MAX:'     '{lmax}'\n"
           f"'LAMBDA_STEP:'    '0.005'\n'INTENSITY/FLUX:' 'Flux'\n'ABFIND        :' '.false.'\n"
           f"'MODELOPAC:' '{opac}'\n'RESULTFILE :' '{res}'\n'METALLICITY:'    '0.00'\n"
           f"'INDIVIDUAL ABUNDANCES:'   '1'\n38  {a_x:.3f}\n'ISOTOPES : ' '0'\n"
           f"'NFILES   :' '1'\n{linelist}\n'SPHERICAL:'  'F'\n  30\n  300.00\n  15\n  1.30\n")
    r = subprocess.run([f"{EXE}/bsyn_lu"], input=ctl, capture_output=True, text=True, cwd=W)
    if not os.path.exists(res):
        raise SystemExit(f"bsyn failed @ {center} A={a_x}:\n{r.stdout[-1500:]}")
    d = np.loadtxt(res)
    return d[:, 0], d[:, 1]


def rot_kernel(dv, vsini, eps=0.6):
    x = dv / vsini
    m = np.abs(x) < 1.0
    k = np.zeros_like(dv)
    c1 = 2 * (1 - eps) / (np.pi * vsini * (1 - eps / 3.0))
    c2 = 0.5 * eps / (vsini * (1 - eps / 3.0))
    k[m] = c1 * np.sqrt(1 - x[m] ** 2) + c2 * (1 - x[m] ** 2)
    return k


def broaden(wl, fl, vsini, gsig_kms):
    """Rotational (vsini) then Gaussian (gsig_kms, absorbing vmac+instrumental)."""
    step = wl[1] - wl[0]
    cen = np.median(wl)
    dv = (wl - cen) / cen * CLIGHT
    out = fl
    if vsini > 0.1:
        kv = np.arange(-vsini * 1.2, vsini * 1.2 + step, step) / cen * CLIGHT
        rk = rot_kernel(kv, vsini)
        if rk.sum() > 0:
            rk /= rk.sum()
            out = np.convolve(1 - out, rk, mode='same')
            out = 1 - out
    if gsig_kms > 0.05:
        gsig_A = gsig_kms / CLIGHT * cen
        n = int(np.ceil(4 * gsig_A / step))
        gx = np.arange(-n, n + 1) * step
        gk = np.exp(-0.5 * (gx / gsig_A) ** 2); gk /= gk.sum()
        out = np.convolve(1 - out, gk, mode='same'); out = 1 - out
    return out


def local_renorm(w, f, center, hw):
    """Divide observed flux by a linear continuum fit to the top-percentile points
    in the [center-hw-1, center+hw+1] window edges."""
    win = (w > center - hw - 1.2) & (w < center + hw + 1.2)
    ww, ff = w[win], f[win]
    if len(ww) < 20:
        return w, f, win
    edge = (ww < center - hw * 0.7) | (ww > center + hw * 0.7)
    xe, ye = ww[edge], ff[edge]
    if len(xe) < 6:
        return ww, ff, win
    thr = np.percentile(ye, 70)
    keep = ye >= thr
    c = np.polyfit(xe[keep], ye[keep], 1)
    cont = np.polyval(c, ww)
    return ww, ff / np.clip(cont, 1e-3, None), win


def load_harps():
    import csv
    w, f = [], []
    with open(HARPS) as fh:
        r = csv.reader(fh); next(r)
        for row in r:
            w.append(float(row[0])); f.append(float(row[3]))
    return np.array(w), np.array(f)


def load_iag(repo="/mnt/codex-data/codex/repo"):
    """IAG Reiners+2016 visible FTS solar flux atlas: cols = vacuum wavenumber(cm^-1),
    normalized flux. Convert vacuum wavelength -> air via the pipeline's SINGLE-SOURCE
    converter (pipeline.wavelength_util.vac_to_air, RYA-264/501) — never re-derive the
    refractive-index formula locally."""
    import gzip
    sys.path.insert(0, repo)
    from pipeline.wavelength_util import vac_to_air
    wn, fl = [], []
    with gzip.open(IAG, 'rt') as fh:
        for ln in fh:
            p = ln.split()
            if len(p) < 2:
                continue
            try:
                wn.append(float(p[0])); fl.append(float(p[1]))
            except ValueError:
                continue
    wn = np.array(wn); fl = np.array(fl)
    lam_air = vac_to_air(1e8 / wn)         # 1e8/wavenumber = vacuum Angstrom
    idx = np.argsort(lam_air)
    return lam_air[idx], fl[idx]


def fit_line(center, obs_w, obs_f, synth):
    """synth = dict A -> (wl, fl_norm). Fit A by min chi2 over the fit window, with a
    joint Gaussian broadening nuisance (grid) + fixed rotational vsini."""
    hw = LINES[center]['fit_hw']
    ww, ff, _ = local_renorm(obs_w, obs_f, center, hw)
    fitm = (ww > center - hw) & (ww < center + hw)
    if fitm.sum() < 10:
        return None
    xo, yo = ww[fitm], ff[fitm]
    gsig_grid = np.arange(1.5, 7.0, 0.5)    # km/s combined macro+instrumental
    best = dict(chi2=1e30)
    for a, (sw, sf) in synth.items():
        for gs in gsig_grid:
            sb = broaden(sw, sf, VSINI, gs)
            ys = np.interp(xo, sw, sb)
            chi2 = float(np.sum((yo - ys) ** 2)) / max(len(xo) - 2, 1) / (0.01 ** 2)
            if chi2 < best['chi2']:
                best = dict(chi2=chi2, A=float(a), gsig=float(gs), npix=int(len(xo)))
    # parabolic refine in A around the min (fixed gsig)
    gs = best['gsig']
    As = np.array(sorted(synth))
    chis = []
    for a in As:
        sw, sf = synth[a]
        sb = broaden(sw, sf, VSINI, gs)
        ys = np.interp(xo, sw, sb)
        chis.append(np.sum((yo - np.interp(xo, sw, sb)) ** 2))
    chis = np.array(chis)
    k = int(np.argmin(chis))
    A_ref = float(As[k])
    if 0 < k < len(As) - 1:
        d = chis[k + 1] - 2 * chis[k] + chis[k - 1]
        if d > 0:
            A_ref = float(As[k] - 0.5 * (chis[k + 1] - chis[k - 1]) / d * (As[1] - As[0]))
    best['A'] = A_ref
    best['red_chi2'] = best['chi2']
    # sensitivity: change in synthetic core EW (mA) per dex A near the fit (dEW/dA).
    # A weak / blend-dominated line barely responds -> low sensitivity -> unreliable.
    def _core_ew(a):
        aa = min(A_HI, max(A_LO, a))
        k = int(np.argmin(np.abs(As - aa)))
        sw, sf = synth[float(As[k])]
        sb = broaden(sw, sf, VSINI, gs)
        m = (sw > center - 0.4) & (sw < center + 0.4)
        return float(np.trapz(1 - sb[m], sw[m]) * 1000.0)
    best['dEW_dA'] = round(abs(_core_ew(A_ref + 0.15) - _core_ew(A_ref - 0.15)) / 0.30, 1)
    best['railed'] = bool(A_ref <= A_LO + 0.03 or A_ref >= A_HI - 0.03)
    return best


def main():
    # RYA-567: Sirius-only heavy-compute leg (TS Turbospectrum NLTE + MARCS). Refuse to
    # run it off Sirius — loud-fail, never against local-Mac copies of engines/grids.
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from config.constants import assert_on_sirius
    assert_on_sirius("RYA-551 Sr II synthesis", require_subdirs=("engines", "grids"))
    os.makedirs(W, exist_ok=True)
    if not os.path.lexists(f"{W}/DATA"):
        os.symlink("/mnt/codex-data/engines/Turbospectrum_NLTE/DATA", f"{W}/DATA")
    arms = {'harps': load_harps()}
    try:
        arms['iag'] = load_iag()
        print(f"  IAG atlas: {len(arms['iag'][0])} pts, "
              f"{arms['iag'][0].min():.1f}-{arms['iag'][0].max():.1f} A")
    except Exception as e:
        print(f"  IAG load FAILED ({e}); HARPS-only")

    results = {}
    for center, cfg in LINES.items():
        print(f"\n=== Sr II {center} (EP {cfg['ep']}, canonical loggf {cfg['loggf']}) ===")
        ll = patch_linelist(center, cfg['vald'])
        lmin, lmax = center - LMARGIN, center + LMARGIN
        opac = f"{W}/opac_{center:.0f}"
        babsma(center, opac, lmin, lmax)
        synth = {}
        for a in A_GRID:
            synth[float(a)] = bsyn(center, a, ll, opac, lmin, lmax, f"{a:.2f}")
        delta = inspect_delta(cfg['feh_node_delta'])
        rec = dict(wave=center, ep=cfg['ep'], loggf=cfg['loggf'],
                   nlte_delta=delta,
                   nlte_status=('applied' if delta is not None else 'NLTE_unavailable'))
        for arm, (ow, of) in arms.items():
            fit = fit_line(center, ow, of, synth)
            if fit is None:
                rec[arm] = dict(status='no_coverage')
                print(f"  {arm:6s}: no coverage")
                continue
            a_lte = fit['A']
            a_nlte = a_lte + delta if delta is not None else a_lte
            rec[arm] = dict(A_LTE=round(a_lte, 3), A_NLTE=round(a_nlte, 3),
                            gsig_kms=round(fit['gsig'], 2),
                            red_chi2=round(fit['red_chi2'], 2), npix=fit['npix'],
                            dEW_dA_mA_dex=fit['dEW_dA'], railed=fit['railed'],
                            reliable=bool((not fit['railed']) and fit['dEW_dA'] >= 40.0))
            print(f"  {arm:6s}: A_LTE={a_lte:.3f}  A_NLTE={a_nlte:.3f}  "
                  f"(delta={delta}, gsig={fit['gsig']:.1f} km/s, rchi2={fit['red_chi2']:.1f}, "
                  f"dEW/dA={fit['dEW_dA']} mA/dex, railed={fit['railed']})")
        results[str(center)] = rec

    with open(f"{W}/rya551_sr2_synth_result.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {W}/rya551_sr2_synth_result.json")


if __name__ == '__main__':
    main()
