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
import argparse
import csv
import json
import os
import re
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from solar_profile_fit import (CLIGHT, RCHI2_REVIEW,  # noqa: E402,F401
                               RELIABLE_DEWDA, assess_reliability, broaden,
                               fit_profile, fit_profile_deblend, local_renorm,
                               measure_arm_rv, require_arm_rv)
from config.constants import codex_path  # RYA-810 path register

EXE   = str(codex_path('engines.ts_exec'))
MARCS = str(codex_path('grids.marcs_standard')
           / "p5750_g+4.5_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod")
VALD_DIR = str(codex_path('engines.vald_linelists'))
W = str(codex_path('work.solar_tmp') / 'rya560_work')
HARPS = str(codex_path('work.solar_harps_normalized'))
IAG   = str(codex_path('data.solar_iag_atlas'))

Z_ZR  = 40            # Zr atomic number (bsyn INDIVIDUAL ABUNDANCES element)
XI    = 1.0           # solar microturbulence (matches the t01 MARCS node)
VSINI = 1.8           # solar vsini (STAR_PARAMS)

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

# ---------------------------------------------------------------- RYA-585 deblend
# The three strong low-EP lines RYA-560 could actually measure. 4629/5350/5372 are
# excluded: RYA-560 found them railed/junk, and a blend model cannot rescue a line
# whose fit has no interior minimum.
DEBLEND_LINES = (4208.980, 4258.041, 4442.992)
DEBLEND_PRIMARY = 4208.980
# A_X passed to bsyn to suppress the target element and synthesise the blends alone.
BLEND_ONLY_ABUND = -5.0


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


def load_harps():
    w, f = [], []
    with open(HARPS) as fh:
        r = csv.reader(fh); next(r)
        for row in r:
            w.append(float(row[0])); f.append(float(row[3]))
    return np.array(w), np.array(f)


def load_iag(repo=str(codex_path('work.repo'))):
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
    """Thin adapter onto the SHARED fitter (scripts/solar_profile_fit.fit_profile).

    RYA-643: this harness previously carried its OWN copy of the fitter, which is how
    it kept the two defects RYA-592 fixed in ITS copy — (a) no rest-frame handling at
    all, so the fit compared the observed profile against rest-frame synthesis while
    the solar arms carry a real velocity offset, and (b) a broadening grid starting at
    1.5 km/s that the fit railed against. Both now come from the single shared source;
    dv is fitted as a nuisance and cross-checked against the abundance-blind
    measure_arm_rv()."""
    return fit_profile(center, obs_w, obs_f, synth,
                       hw=LINES[center]['fit_hw'], vsini=VSINI,
                       a_lo=A_LO, a_hi=A_HI)


def _setup(label):
    """Shared prologue for both run modes: Sirius gate, SSOT gf gate, arms, arm RVs."""
    # RYA-567: Sirius-only heavy-compute leg (TS Turbospectrum + MARCS). Refuse to run
    # off Sirius — loud-fail, never against local-Mac copies of engines/grids/spectra.
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, root)
    from config.constants import assert_on_sirius, SOLAR_ASPLUND2021
    assert_on_sirius(label, require_subdirs=("engines", "grids"))

    a_sun = float(SOLAR_ASPLUND2021['Zr'])        # 2.59 — reference point, NOT a tuning target
    print(f"{label} — A(Zr)_sun reference (Asplund 2021) = {a_sun} "
          f"[grid-bracket only]\n")
    print("Step-0 SSOT check (gf from data/linelists/canonical_gf.csv):")
    verify_canonical_gf(root)

    os.makedirs(W, exist_ok=True)
    if not os.path.lexists(f"{W}/DATA"):
        os.symlink(str(codex_path('engines.ts_data')), f"{W}/DATA")
    arms = {'harps': load_harps()}
    try:
        arms['iag'] = load_iag()
        print(f"\n  IAG atlas: {len(arms['iag'][0])} pts, "
              f"{arms['iag'][0].min():.1f}-{arms['iag'][0].max():.1f} A")
    except Exception as e:
        print(f"\n  IAG load FAILED ({e}); HARPS-only")

    # RYA-643: the rest-frame correction must be SOURCED from a measurement on
    # each arm being fitted — never hardcoded, never a silent zero. Loud-fails if
    # the clean check lines cannot deliver one.
    arm_rv = {}
    print('  residual velocity per arm (clean-line centroids, abundance-blind):')
    for _arm, (_ow, _of) in arms.items():
        _v, _n, _sd = require_arm_rv(_ow, _of, _arm)
        arm_rv[_arm] = dict(v_kms=round(_v, 3), n_lines=_n, scatter_kms=round(_sd, 3))
        print(f'    {_arm:6s}: {_v:+.3f} km/s (n={_n}, scatter {_sd:.3f})')
    return root, a_sun, arms, arm_rv


def main():
    root, a_sun, arms, arm_rv = _setup("RYA-560 Zr II synthesis")

    results = {'_meta': dict(element='Zr', ion='II', z=Z_ZR,
                             a_sun_ref=a_sun, a_sun_ref_source='Asplund2021',
                             nlte='NLTE_unavailable_LTE_robust',
                             nlte_rationale=('Zr II is the majority ion in the solar '
                                             'photosphere -> LTE-robust (registry 279/458, '
                                             'Sr II/V II precedent); no departure grid exists'),
                             reliable_dEW_dA_floor=RELIABLE_DEWDA,
                             red_chi2_review_trigger=RCHI2_REVIEW,
                             red_chi2_gates_reliable=False)}
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
            rel_f = assess_reliability(fit)
            rec[arm] = dict(dv_fitted_kms=round(fit['dv'], 2),
                            dv_measured_kms=arm_rv[arm]['v_kms'],
                            gsig_railed=fit['gsig_railed'],
                            A_LTE=round(a_lte, 3), A_NLTE=None,
                            gsig_kms=round(fit['gsig'], 2),
                            red_chi2=round(fit['red_chi2'], 2), npix=fit['npix'],
                            dEW_dA_mA_dex=fit['dEW_dA'], railed=fit['railed'],
                            **rel_f)
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


def main_deblend():
    """RYA-585 — in-window blend-fit rescue of the three strong Zr II lines.

    RYA-560 measured A(Zr II) ~ 2.4 with cross-arm agreement but could not emit it:
    every line sat below the dEW/dA >= 40 mA/dex floor (best 36.8) with red_chi2
    41-91 and ~0.2 dex line-to-line scatter. The elevated red_chi2 said the profile
    model was wrong, which is a fixable defect; this mode fixes it and then re-tests
    the floor, so the two failure modes stop being confounded.

    Runs the RYA-560 syntheses UNCHANGED and fits each line BOTH ways — the RYA-560
    `fit_profile` and the RYA-585 `fit_profile_deblend` — so the before/after is a
    controlled comparison, not a re-run. The extra synthesis per line is
    `blend_only`: the same window with Zr suppressed, which is what makes the
    target's own contribution (and hence the continuum and chi2 masks) knowable.

    DECISION GATE (RYA-679 ratified rule, via `assess_reliability`): a line is
    `reliable` only if dEW/dA >= RELIABLE_DEWDA AND not railed. red_chi2 is REPORTED
    and review-flagged above RCHI2_REVIEW but does NOT gate — RYA-679 measured that
    full-window red_chi2 tracks blend-list fidelity rather than the element, using
    this very line set (4208.98: 83.12 full-window vs 0.39 deblended, same data).
    Zr's outcome is unchanged either way: it fails on sensitivity, not on fit quality.
    If no line clears, Zr stays MEASURABLE-OWED and this file records the quantified
    reason. Validate-don't-tune: nothing here is fitted toward Asplund 2.59.
    """
    root, a_sun, arms, arm_rv = _setup("RYA-585 Zr II in-window blend-fit rescue")

    results = {'_meta': dict(element='Zr', ion='II', z=Z_ZR, ticket='RYA-585',
                             basis='RYA-560 syntheses, refitted with blends modelled in-window',
                             a_sun_ref=a_sun, a_sun_ref_source='Asplund2021',
                             nlte='NLTE_unavailable_LTE_robust',
                             reliable_dEW_dA_floor=RELIABLE_DEWDA,
                             red_chi2_review_trigger=RCHI2_REVIEW,
                             red_chi2_gates_reliable=False,
                             primary_line=DEBLEND_PRIMARY,
                             arm_rv=arm_rv)}
    for center in DEBLEND_LINES:
        cfg = LINES[center]
        print(f"\n=== Zr II {center} (EP {cfg['ep']}, canonical loggf {cfg['loggf']}, "
              f"{cfg['gf_ref']}) — deblend ===")
        ll = patch_linelist(center, cfg['vald'])
        lmin, lmax = center - LMARGIN, center + LMARGIN
        opac = f"{W}/opac_{center:.0f}"
        babsma(center, opac, lmin, lmax)
        synth = {float(a): bsyn(center, a, ll, opac, lmin, lmax, f"{a:.2f}")
                 for a in A_GRID}
        # blends alone: identical window/linelist, Zr suppressed.
        blend_only = bsyn(center, BLEND_ONLY_ABUND, ll, opac, lmin, lmax, 'blendonly')
        rec = dict(wave=center, ep=cfg['ep'], loggf=cfg['loggf'],
                   gf_id=cfg['gf_id'], gf_ref=cfg['gf_ref'],
                   nlte_delta=None, nlte_status='NLTE_unavailable_LTE_robust')
        for arm, (ow, of) in arms.items():
            old = fit_profile(center, ow, of, synth, hw=cfg['fit_hw'], vsini=VSINI,
                              a_lo=A_LO, a_hi=A_HI)
            new = fit_profile_deblend(center, ow, of, synth, blend_only,
                                      hw=cfg['fit_hw'], vsini=VSINI,
                                      a_lo=A_LO, a_hi=A_HI)
            if old is None or new is None:
                rec[arm] = dict(status='no_coverage')
                print(f"  {arm:6s}: no coverage")
                continue
            rel_new = assess_reliability(new)
            reliable = rel_new['reliable']
            rec[arm] = dict(
                before=dict(A_LTE=round(old['A'], 3),
                            red_chi2=round(old['red_chi2'], 2),
                            dEW_dA_mA_dex=old['dEW_dA'], railed=old['railed'],
                            **assess_reliability(old)),
                A_LTE=round(new['A'], 3), A_NLTE=None,
                red_chi2=round(new['red_chi2'], 2),
                dEW_dA_mA_dex=new['dEW_dA'], railed=new['railed'],
                **rel_new,
                gsig_kms=round(new['gsig'], 2), gsig_railed=new['gsig_railed'],
                dv_fitted_kms=round(new['dv'], 2),
                dv_measured_kms=arm_rv[arm]['v_kms'],
                npix_target=new['npix'], npix_window=new['npix_window'],
                cont_slope=round(new['cont1'], 5), cont_level=round(new['cont0'], 5),
                cont_scatter=round(new['cont_scatter'], 5),
                cont_npix=new['cont_npix'],
                target_EW_mA=new['target_EW_mA'], sat_index=new['sat_index'],
                target_core_depth_frac=new['target_core_depth_frac'])
            print(f"  {arm:6s} BEFORE: A={old['A']:.3f} rchi2={old['red_chi2']:7.2f} "
                  f"dEW/dA={old['dEW_dA']:6.1f} railed={old['railed']}")
            print(f"  {arm:6s} AFTER : A={new['A']:.3f} rchi2={new['red_chi2']:7.2f} "
                  f"dEW/dA={new['dEW_dA']:6.1f} railed={new['railed']} "
                  f"reliable={reliable}")
            print(f"  {arm:6s}         cont={new['cont0']:.4f}{new['cont1']:+.4f}x "
                  f"(sd {new['cont_scatter']:.4f}, n {new['cont_npix']}) "
                  f"npix={new['npix']}/{new['npix_window']} "
                  f"targetEW={new['target_EW_mA']:.1f} mA sat={new['sat_index']} "
                  f"depthfrac={new['target_core_depth_frac']} "
                  f"gsig={new['gsig']:.1f} dv={new['dv']:+.2f} "
                  f"(measured {arm_rv[arm]['v_kms']:+.3f})")
        results[str(center)] = rec

    out_path = os.path.join(root, "data", "results", "zr2_deblend_rya585.json")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    rel = [(c, results[str(c)]['harps']['A_LTE'])
           for c in DEBLEND_LINES
           if isinstance(results[str(c)].get('harps'), dict)
           and results[str(c)]['harps'].get('reliable')]
    print("\n" + "=" * 72)
    if rel:
        vals = np.array([v for _, v in rel])
        verdict = 'CLEARED'
        print(f"RELIABLE Zr II lines (HARPS, deblended): "
              f"{', '.join(f'{c}={v:.3f}' for c, v in rel)}")
        print(f"A(Zr II) LTE = {vals.mean():.3f}  (n={len(vals)}, "
              f"scatter {vals.std(ddof=1) if len(vals) > 1 else 0.0:.3f})")
    else:
        verdict = 'OWED'
        best = max((results[str(c)]['harps']['dEW_dA_mA_dex'] for c in DEBLEND_LINES
                    if isinstance(results[str(c)].get('harps'), dict)), default=0.0)
        worst_chi = max((results[str(c)]['harps']['red_chi2'] for c in DEBLEND_LINES
                         if isinstance(results[str(c)].get('harps'), dict)), default=0.0)
        print("NO Zr II line cleared BOTH RYA-585 gates.")
        print(f"  red_chi2 (reported, NOT a gate — RYA-679): worst {worst_chi:.2f} "
              f"vs review trigger {RCHI2_REVIEW} "
              f"— the blend/continuum systematic was real and is now modelled.")
        print(f"  dEW/dA   gate  : FAILED (best {best} < {RELIABLE_DEWDA} mA/dex) "
              f"— residual failure is intrinsic line sensitivity, not the blend model.")
        print("=> Zr remains MEASURABLE-OWED. This LINE SET is exhausted; the next "
              "lever is acquiring cleaner blue Zr II lines (RYA-458), not refitting "
              "these. Do NOT emit a sub-floor value.")
    results['_meta']['verdict'] = verdict
    results['_meta']['reliable_lines'] = [c for c, _ in rel]
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {out_path}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--deblend', action='store_true',
                    help='RYA-585: refit the three strong Zr II lines with blends '
                         'modelled in-window and a blend-pixel continuum, and re-test '
                         'the reliability floor. Writes zr2_deblend_rya585.json. '
                         'Without this flag the script runs the original RYA-560 '
                         'measurement unchanged.')
    args = ap.parse_args()
    main_deblend() if args.deblend else main()
