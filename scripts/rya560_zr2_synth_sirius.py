#!/usr/bin/env python3
"""
RYA-560 — Zirconium (Zr II) synthesis abundance measurement (RUNS ON SIRIUS).

Wires the synthesis route for Zr II (RYA-524 sweep flag, class MERGED-NOT-WIRED):
the strong optical Zr II lines EW-saturate (COG), so a naive Gaussian EW is culled;
we measure A(Zr) instead by Turbospectrum LTE flux synthesis WITH blend handling
(the full VALD block in-window). Adapts scripts/rya551_sr2_synth_sirius.py verbatim
for babsma/bsyn/broaden/local_renorm/fit_line/load_harps/load_iag.

MAJORITY-ION / LTE-ROBUST — NO NLTE GRID (registry problem_children 279/458, the
Sr II / V II precedent): Zr II is the dominant ionisation stage of Zr in the solar
photosphere, so Zr II line formation is close to LTE and there is no departure grid
to apply. We therefore set nlte_delta=None, nlte_status='NLTE_unavailable_LTE_robust'
and report A_LTE as THE value (not a floor awaiting a correction).

Step 0 (gf provenance) is ALREADY SATISFIED in the committed single source:
data/linelists/canonical_gf.csv carries the Zr II optical lines with cited gf —
LNAJ = Ljung, Nilsson, Asplund & Johansson 2006, A&A 456, 1181 (the FERRUM-project
experimental Zr II oscillator strengths, which themselves derive the solar Zr
abundance) for 4208.98/4258.04/4442.99/5372.47, and the Gaia-ESO (CC/CCout) codes
for 4629.08/5350.09. This script does NOT re-add or fabricate any gf; it READS the
canonical values and asserts the hardcoded copy matches the single source (loud-fail
on drift — SSOT invariant, no silent fallback). The ticket's "no Zr II in the line
list" premise predates the RYA-492 canonical_gf regeneration (commit 3e9c98b).

Two solar arms (as in 551): HARPS reflected-solar (primary) + IAG FTS atlas (cross).

Validate-don't-tune: report what the fit gives; never adjust toward Asplund 2021
A(Zr)=2.59 (which is only the grid-bracketing reference point, NOT a target).
"""
import csv
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
W = "/mnt/codex-data/codex/_solar_tmp/rya560_work"
HARPS = "/mnt/codex-data/codex/_solar_tmp/solar_normalized_harps.csv"
IAG   = "/mnt/codex-data/solar_reference/iag_reiners2016/spvis.dat.gz"

Z_ZR  = 40            # Zr atomic number (bsyn INDIVIDUAL ABUNDANCES element)
XI    = 1.0           # solar microturbulence (matches the t01 MARCS node)
VSINI = 1.8           # solar vsini (STAR_PARAMS)
CLIGHT = 299792.458

# Canonical Zr II lines: (wave_air, EP_eV, log_gf) — values are the CANONICAL
# single-source gf from data/linelists/canonical_gf.csv (line_id + reference cited);
# verify_canonical_gf() asserts these match the file at runtime (SSOT, loud-fail).
LINES = {
    4208.980: dict(ep=0.713, loggf=-0.51, gf_id='gf_145794', gf_ref='LNAJ (Ljung+2006)',
                   vald='vald-3800-4300-for-grid.list', fit_hw=1.6),
    4258.041: dict(ep=0.559, loggf=-1.20, gf_id='gf_145801', gf_ref='LNAJ (Ljung+2006)',
                   vald='vald-3800-4300-for-grid.list', fit_hw=1.2),
    4442.992: dict(ep=1.486, loggf=-0.42, gf_id='gf_145833', gf_ref='LNAJ (Ljung+2006)',
                   vald='vald-4300-4800-for-grid.list', fit_hw=1.4),
    4629.079: dict(ep=2.490, loggf=-0.59, gf_id='gf_145848', gf_ref='CCout (Gaia-ESO)',
                   vald='vald-4300-4800-for-grid.list', fit_hw=0.9),
    5350.089: dict(ep=1.827, loggf=-1.24, gf_id='gf_145861', gf_ref='CC (Gaia-ESO)',
                   vald='vald-5300-5800-for-grid.list', fit_hw=0.9),
    5372.466: dict(ep=2.409, loggf=-0.81, gf_id='gf_145863', gf_ref='LNAJ (Ljung+2006, RYA-354)',
                   vald='vald-5300-5800-for-grid.list', fit_hw=0.9),
}
LMARGIN = 8.0                                       # synth half-window (A) for babsma/bsyn
A_GRID = np.round(np.arange(1.90, 3.26, 0.05), 3)   # A(Zr) trial grid (brackets Asplund 2.59)
A_LO, A_HI = float(A_GRID.min()), float(A_GRID.max())

RELIABLE_DEWDA = 40.0     # mA/dex core-EW sensitivity floor (as RYA-551)


def verify_canonical_gf(root):
    """SSOT invariant: assert every LINES gf equals the committed canonical_gf.csv
    value (Step-0 provenance gate). Loud-fail on drift or missing line — never run
    the synthesis on an unsourced/fabricated gf, never silently fall back."""
    path = os.path.join(root, "data", "linelists", "canonical_gf.csv")
    if not os.path.exists(path):
        raise SystemExit(f"RYA-560 SSOT: canonical_gf.csv not found at {path}")
    have = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            if row["species"] == "Zr II":
                have[round(float(row["wavelength_air_A"]), 3)] = row
    for center, cfg in LINES.items():
        key = round(center, 3)
        row = have.get(key)
        if row is None:                      # try a small tolerance match
            for w, r in have.items():
                if abs(w - center) < 0.01:
                    row = r; break
        if row is None:
            raise SystemExit(f"RYA-560 SSOT: no Zr II {center} in canonical_gf.csv — "
                             f"cannot source gf (Step-0 gate). Fix the single source, do not hardcode.")
        cg = float(row["log_gf"])
        if abs(cg - cfg["loggf"]) > 0.005:
            raise SystemExit(f"RYA-560 SSOT DRIFT: Zr II {center} canonical log_gf={cg} "
                             f"!= script {cfg['loggf']} ({row['line_id']}). Reconcile to the single source.")
        print(f"  gf OK  Zr II {center}: loggf={cg:+.3f}  ({row['line_id']}, {cfg['gf_ref']})")


def patch_linelist(center, vald_name):
    """Copy the VALD block for the window and overwrite the target Zr II line's log_gf
    with the canonical value so the engine synthesises on OUR single-source gf scale.
    bsyn selects lines by LAMBDA window; blends stay verbatim (COG blend handling)."""
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
                if abs(wl - center) < 0.02 and "Zr II" in ln:
                    m = re.match(r"^(\s*\S+\s+\S+\s+)(\S+)(\s+.*)$", ln.rstrip("\n"))
                    if m:
                        out.append(f"{m.group(1)}{loggf:.3f}{m.group(3)}\n")
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
           f"'INDIVIDUAL ABUNDANCES:'   '1'\n{Z_ZR}  {a_x:.3f}\n'ISOTOPES : ' '0'\n"
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
    w, f = [], []
    with open(HARPS) as fh:
        r = csv.reader(fh); next(r)
        for row in r:
            w.append(float(row[0])); f.append(float(row[3]))
    return np.array(w), np.array(f)


def load_iag(repo="/mnt/codex-data/codex/repo"):
    """IAG Reiners+2016 visible FTS solar flux atlas: cols = vacuum wavenumber(cm^-1),
    normalized flux. Vacuum->air via the pipeline SINGLE-SOURCE converter."""
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
    lam_air = vac_to_air(1e8 / wn)
    idx = np.argsort(lam_air)
    return lam_air[idx], fl[idx]


def fit_line(center, obs_w, obs_f, synth):
    hw = LINES[center]['fit_hw']
    ww, ff, _ = local_renorm(obs_w, obs_f, center, hw)
    fitm = (ww > center - hw) & (ww < center + hw)
    if fitm.sum() < 10:
        return None
    xo, yo = ww[fitm], ff[fitm]
    gsig_grid = np.arange(1.5, 7.0, 0.5)
    best = dict(chi2=1e30)
    for a, (sw, sf) in synth.items():
        for gs in gsig_grid:
            sb = broaden(sw, sf, VSINI, gs)
            ys = np.interp(xo, sw, sb)
            chi2 = float(np.sum((yo - ys) ** 2)) / max(len(xo) - 2, 1) / (0.01 ** 2)
            if chi2 < best['chi2']:
                best = dict(chi2=chi2, A=float(a), gsig=float(gs), npix=int(len(xo)))
    gs = best['gsig']
    As = np.array(sorted(synth))
    chis = []
    for a in As:
        sw, sf = synth[a]
        sb = broaden(sw, sf, VSINI, gs)
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

    def _core_ew(a):
        aa = min(A_HI, max(A_LO, a))
        kk = int(np.argmin(np.abs(As - aa)))
        sw, sf = synth[float(As[kk])]
        sb = broaden(sw, sf, VSINI, gs)
        m = (sw > center - 0.4) & (sw < center + 0.4)
        return float(np.trapz(1 - sb[m], sw[m]) * 1000.0)
    best['dEW_dA'] = round(abs(_core_ew(A_ref + 0.15) - _core_ew(A_ref - 0.15)) / 0.30, 1)
    best['railed'] = bool(A_ref <= A_LO + 0.03 or A_ref >= A_HI - 0.03)
    return best


def main():
    # RYA-567: Sirius-only heavy-compute leg (TS Turbospectrum + MARCS). Refuse to run
    # off Sirius — loud-fail, never against local-Mac copies of engines/grids/spectra.
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, root)
    from config.constants import assert_on_sirius, SOLAR_ASPLUND2021
    assert_on_sirius("RYA-560 Zr II synthesis", require_subdirs=("engines", "grids"))

    a_sun = float(SOLAR_ASPLUND2021['Zr'])        # 2.59 — reference point, NOT a tuning target
    print(f"RYA-560 Zr II synthesis — A(Zr)_sun reference (Asplund 2021) = {a_sun} "
          f"[grid-bracket only]\n")
    print("Step-0 SSOT check (gf from data/linelists/canonical_gf.csv):")
    verify_canonical_gf(root)

    os.makedirs(W, exist_ok=True)
    if not os.path.lexists(f"{W}/DATA"):
        os.symlink("/mnt/codex-data/engines/Turbospectrum_NLTE/DATA", f"{W}/DATA")
    arms = {'harps': load_harps()}
    try:
        arms['iag'] = load_iag()
        print(f"\n  IAG atlas: {len(arms['iag'][0])} pts, "
              f"{arms['iag'][0].min():.1f}-{arms['iag'][0].max():.1f} A")
    except Exception as e:
        print(f"\n  IAG load FAILED ({e}); HARPS-only")

    results = {'_meta': dict(element='Zr', ion='II', z=Z_ZR,
                             a_sun_ref=a_sun, a_sun_ref_source='Asplund2021',
                             nlte='NLTE_unavailable_LTE_robust',
                             nlte_rationale=('Zr II is the majority ion in the solar '
                                             'photosphere -> LTE-robust (registry 279/458, '
                                             'Sr II/V II precedent); no departure grid exists'),
                             reliable_dEW_dA_floor=RELIABLE_DEWDA)}
    for center, cfg in LINES.items():
        print(f"\n=== Zr II {center} (EP {cfg['ep']}, canonical loggf {cfg['loggf']}, "
              f"{cfg['gf_ref']}) ===")
        ll = patch_linelist(center, cfg['vald'])
        lmin, lmax = center - LMARGIN, center + LMARGIN
        opac = f"{W}/opac_{center:.0f}"
        babsma(center, opac, lmin, lmax)
        synth = {}
        for a in A_GRID:
            synth[float(a)] = bsyn(center, a, ll, opac, lmin, lmax, f"{a:.2f}")
        # NLTE: none — Zr II majority ion, LTE-robust. A_LTE IS the value.
        rec = dict(wave=center, ep=cfg['ep'], loggf=cfg['loggf'],
                   gf_id=cfg['gf_id'], gf_ref=cfg['gf_ref'],
                   nlte_delta=None, nlte_status='NLTE_unavailable_LTE_robust')
        for arm, (ow, of) in arms.items():
            fit = fit_line(center, ow, of, synth)
            if fit is None:
                rec[arm] = dict(status='no_coverage')
                print(f"  {arm:6s}: no coverage")
                continue
            a_lte = fit['A']
            reliable = bool((not fit['railed']) and fit['dEW_dA'] >= RELIABLE_DEWDA)
            rec[arm] = dict(A_LTE=round(a_lte, 3), A_NLTE=None,
                            gsig_kms=round(fit['gsig'], 2),
                            red_chi2=round(fit['red_chi2'], 2), npix=fit['npix'],
                            dEW_dA_mA_dex=fit['dEW_dA'], railed=fit['railed'],
                            reliable=reliable)
            print(f"  {arm:6s}: A_LTE={a_lte:.3f}  (LTE-robust, no NLTE; "
                  f"gsig={fit['gsig']:.1f} km/s, rchi2={fit['red_chi2']:.1f}, "
                  f"dEW/dA={fit['dEW_dA']} mA/dex, railed={fit['railed']}, reliable={reliable})")
        results[str(center)] = rec

    out_dir = os.path.join(root, "data", "results")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "zr2_synthesis_rya560.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

    # Summary: reliable HARPS lines -> the reportable value (else still-owed, loud).
    rel = [(c, results[str(c)]['harps']['A_LTE'])
           for c in LINES
           if isinstance(results[str(c)].get('harps'), dict)
           and results[str(c)]['harps'].get('reliable')]
    print("\n" + "=" * 62)
    if rel:
        vals = np.array([v for _, v in rel])
        print(f"RELIABLE Zr II lines (HARPS): "
              f"{', '.join(f'{c}={v:.3f}' for c, v in rel)}")
        print(f"A(Zr II) LTE = {vals.mean():.3f}  (n={len(vals)}, "
              f"scatter {vals.std(ddof=1) if len(vals) > 1 else 0.0:.3f})")
    else:
        print("NO Zr II line cleared the reliability floor "
              f"(dEW/dA >= {RELIABLE_DEWDA} mA/dex, not railed).")
        print("=> Zr remains MEASURABLE-OWED. Report reason; do NOT emit a weak/railed value.")
    print(f"\nWrote {out_path}")


if __name__ == '__main__':
    main()
