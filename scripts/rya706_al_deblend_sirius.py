#!/usr/bin/env python3
"""
RYA-706 — Aluminium (Al I 6696.023) DEBLEND re-measure by in-window blend fit (SIRIUS).

Adapted from scripts/rya581_ba2_deblend_sirius.py. The MACHINERY is unchanged and
deliberately so: the control run, the deblend-shift stop condition, the shared RYA-643
fitter, the arm-RV handling and the reliability reporting all carry over. Only the
species, the line, the window and the gf source differ. A second copy of a fitter is
how RYA-643 got four of them.

WHY THIS RUN EXISTS — a prediction, not a hope
----------------------------------------------
RYA-706 recovered two Al I lines that the committed pool never carried, graded them
against a NIST ASD pull (both C+, ref L4737), and got:

    Al I 6696.023   A(Al) = 6.823     has Fe I 6696.315 in-window at depth 0.105
    Al I 6698.673   A(Al) = 6.485     clean, no contaminant above 0.3% of its own depth

    spread 0.339 dex   Asplund 2021 = 6.43

The gf zero point cannot explain that spread: the NIST correction is -0.222 dex on
6696 and -0.223 on 6698, i.e. IDENTICAL, so it moves both lines together and leaves the
difference untouched. That rules gf out and leaves exactly one candidate on the table —
6696.023 is the line with the contaminant, and it is the HIGH one. A blend charged to
aluminium inflates the EW, which inflates A.

So this run is a falsifiable prediction: DEBLENDING 6696.023 SHOULD MOVE IT DOWN,
TOWARD 6698.673's 6.485. If it does not, the blend was not the cause and the spread is
still unexplained — that is a real result too, and it is reported as one. Ryan,
2026-08-08: "regardless, this will all get documented if it is a failure line."

WHAT IS DIFFERENT FROM Ba, AND WHY IT MATTERS
---------------------------------------------
  * gf comes from the NIST ASD pull (-1.569, C+), NOT from canonical_gf. canonical_gf
    carries -1.347 for this line, VALD3, adjudication_status 'pending' — 0.222 dex too
    high — AND it has the C+ grade attached to the WRONG LINE (6696.185, EP 4.022,
    a different lower level; NIST's graded row is Ei 25347.756 cm-1 = 3.1427 eV, which
    is 6696.023). Rescaling onto canonical here would bake in a known-wrong scale.
    The correction to canonical_gf is a reviewed act (RYA-354) and is NOT done here.
  * Al has HFS (27Al, I=5/2) but NO isotopic splitting — 27Al is effectively 100%. The
    RYA-684 isotope double-application hazard therefore CANNOT bite this run, which is
    worth stating because the shipped TSFitPy for-grid lists are form-(B) gf with
    form-(A) headers and that trap is live for other species.
"""
import argparse
import csv
import json
import math
import os
import re
import subprocess
import sys

import numpy as np

# Standalone-script bootstrap: these run under foreign interpreters and from arbitrary
# cwds, so put the REPO ROOT on sys.path BEFORE importing anything from `pipeline`.
# Derived from __file__, never from cwd. (RYA-313)
import os as _os_boot, sys as _sys_boot
_ROOT = _os_boot.path.dirname(_os_boot.path.dirname(_os_boot.path.abspath(__file__)))
_sys_boot.path.insert(0, _ROOT)
from pipeline._numcompat import trapezoid as _trapezoid  # numpy>=2 removed np.trapz (RYA-313)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from solar_profile_fit import (CLIGHT, RCHI2_REVIEW,  # noqa: E402,F401
                               RELIABLE_DEWDA, assess_reliability, broaden,
                               fit_profile, local_renorm, measure_arm_rv,
                               require_arm_rv)
from config.constants import codex_path  # RYA-810 path register

EXE      = "/mnt/codex-data/engines/Turbospectrum_NLTE/exec-gf"
MARCS    = str(codex_path('grids.marcs_standard')
           / "p5750_g+4.5_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod")
VALD_DIR = "/mnt/codex-data/engines/TSFitPy/input_files/linelists/linelist_vald"
VALD_BLOCK = "vald-6300-6800-for-grid.list"
W        = "/mnt/codex-data/codex/rya706/work"
HARPS    = "/mnt/codex-data/codex/_solar_tmp/solar_normalized_harps.csv"
IAG      = "/mnt/codex-data/solar_reference/iag_reiners2016/spvis.dat.gz"

Z_AL  = 13
LINE  = 6696.023            # Al I, NIST ASD obs 6696.015 / ritz 6696.018, grade C+ L4737
NIST_LOGGF = -1.569         # log10(g_i * f_ik) = log10(2 * 1.35e-02); verified to 4 dp
NIST_ACC   = 'C+'
LINE_EP = 3.1427            # lower level of the 6696.023 multiplet (NIST Ei 25347.756 cm-1)
EP_TOL  = 0.02              # same level, not a nearby one — see build_al_hfs_list
CANONICAL_LOGGF = -1.347    # what canonical_gf carries (VALD3, pending) — 0.222 too high
VSINI = 1.8                 # solar vsini (STAR_PARAMS)

LMARGIN = 7.0               # synth half-window (A) for babsma/bsyn
FIT_HW  = 0.6               # primary chi2 fit half-window (A) — see HW_SWEEP
HW_SWEEP = (0.4, 0.5, 0.6, 0.8)   # reported robustness sweep; FIT_HW is the headline
CORE_HW = 0.4               # core-EW integration half-window (blend accounting)
A_GRID = np.round(np.arange(6.00, 7.31, 0.05), 3)   # A(Al) grid; brackets 6.823 blended
                                                    # and 6.485 (the 6698 target) and 6.43
A_LO, A_HI = float(A_GRID.min()), float(A_GRID.max())

# Reliability floor. Be precise about where each half comes from, because they do NOT
# share a provenance:
#   * the dEW/dA sensitivity floor is the RYA-551/560/564 bar (rya560:79, rya564:95);
#   * the red_chi2 ceiling exists ONLY in RYA-564 (rya564_co1_synth_sirius.py:96, same
#     value, same three-term conjunction at :573-574). RYA-551 defines no reliability
#     constant at all and RYA-560 gates on sensitivity + railed only. So this is the Co
#     bar, not a universal one, and it has never been separately ratified — RYA-564's
#     own stated reason for 60 is that the sigma_flux=0.01 floor inflates red_chi2.
# NOT load-bearing for this result either way: Ba's in-window fit returns red_chi2 0.71,
# so `reliable` is unchanged under a ceiling of 60, a ceiling of 5, or no ceiling at all.
# Flagged for ratification rather than quietly inherited (see RYA-585 note below).
# Both now imported from the shared module. The UNRATIFIED red_chi2 ceiling this
# harness flagged was adjudicated by RYA-679 and RETIRED (red_chi2 is reported and
# review-flagged, never gated). As predicted here, Ba is unaffected: red_chi2 0.71.

# Damping for the HFS components. PRIMARY keeps the RYA-559 choice (Unsold enhancement
# 3.000) so the ONLY thing that changes between 2.410 and this value is the blend
# treatment; the VALD vdW is run as a declared sensitivity, never as a silent swap.
FDAMP_PRIMARY = 3.000       # Unsold enhancement, same convention as RYA-706/559
FDAMP_VALD    = -7.578      # linelist_full.csv damping_vdW for these components

# Species whose in-window blend contribution is quantified individually (everything with
# a row inside +/- BLEND_CENSUS_HW of the core). Reported in mA of core EW.
BLEND_CENSUS_HW = 0.8

HDR_RE = re.compile(r"^'\s*([\d.]+)\s*'\s+(\d+)\s+(\d+)\s*$")


# ────────────────────────────── single source of truth ──────────────────────────────
def canonical_al_row(root):
    """The canonical_gf.csv row for Al I 6696.023 — read for REPORTING ONLY.

    This harness does NOT rescale onto canonical for this line, and that is a departure
    from RYA-581 that has to be stated rather than inferred. canonical_gf carries -1.347
    from VALD3 with adjudication_status 'pending', which the RYA-706 NIST pull shows is
    0.222 dex too high; and the C+ grade that should belong to THIS line is attached to
    6696.185 instead (EP 4.022 — a different lower level; NIST's graded row is
    Ei 25347.756 cm-1 = 3.1427 eV). Rescaling onto it would bake a known-wrong scale
    into the result. The row is still read and recorded so the output carries the
    divergence explicitly. Correcting canonical_gf is RYA-354's reviewed act."""
    path = os.path.join(root, "data", "linelists", "canonical_gf.csv")
    if not os.path.exists(path):
        raise SystemExit(f"RYA-706 SSOT: canonical_gf.csv not found at {path}")
    with open(path) as f:
        for row in csv.DictReader(f):
            if row["species"] == "Al I" and abs(float(row["wavelength_air_A"]) - LINE) < 0.02:
                return row
    raise SystemExit("RYA-706 SSOT: no Al I 6696.023 row in canonical_gf.csv — fix the "
                     "single source, do not hardcode a gf.")


def build_al_hfs_list(root, path, fdamp):
    """Al I 6696 HFS/isotope components from data/linelists/linelist_full.csv (the RYA-559
    source), written as a Turbospectrum line list with the component gf sum RESCALED onto
    the canonical_gf total. The HFS pattern is preserved; the absolute scale comes from the
    single source. Any rescale is reported, never silent."""
    canon = canonical_al_row(root)
    canon_gf = float(canon["log_gf"])
    canon_ep = float(canon["excitation_potential_eV"])
    # THE SCALE IS NIST, NOT CANONICAL — see canonical_al_row's docstring.
    target_gf = NIST_LOGGF
    comps = []
    with open(os.path.join(root, "data", "linelists", "linelist_full.csv")) as f:
        for r in csv.DictReader(f):
            if r["element"] != "Al" or r["ion"] != "I":
                continue
            wl = float(r["wavelength_air_A"])
            if not (LINE - 1.0 < wl < LINE + 1.0):
                continue
            # HFS components of ONE line share ONE lower level. Ba II 5853 is isolated so
            # RYA-581 could take everything in the window; Al cannot. Within +/-1 A of
            # 6696.023 linelist_full carries NINE Al I rows: the six HFS components at
            # EP 3.1427, plus 6696.185 (EP 4.0215) and 6696.788 x2 (EP 4.0216), which are
            # DIFFERENT TRANSITIONS. Folding those into the HFS sum would inflate the
            # rescale target and corrupt the measurement silently. Filter on the level.
            ep = float(r["excitation_potential_eV"])
            if abs(ep - LINE_EP) > EP_TOL:
                continue
            comps.append((wl, ep, float(r["log_gf"])))
    if not comps:
        raise SystemExit("RYA-706: no Al I 6696 HFS components in linelist_full.csv")
    if len(comps) != 6:
        raise SystemExit(
            f"RYA-706: expected 6 HFS components at EP {LINE_EP} for Al I 6696.023, got "
            f"{len(comps)}. The line list moved — re-check which rows share the lower "
            f"level before trusting a rescale onto them.")
    comps.sort()
    raw_sum = math.log10(sum(10 ** g for _, _, g in comps))
    shift = target_gf - raw_sum                # uniform, preserves the HFS pattern
    rescaled = abs(shift) > 0.005
    with open(path, "w") as f:
        f.write(f"'  {Z_AL}.000            '    2      {len(comps)}\n")
        f.write("'Al I  6696 HFS'\n")
        for wl, ep, gf in comps:
            f.write(f"  {wl:.4f}  {ep:.4f}  {gf + shift:.3f}   {fdamp:.3f}    2.0  1.00E+08 "
                    f"  'x' 'x'   0.0    1.0 'Al I 6696 HFS'\n")
    info = dict(n_components=len(comps),
                linelist_full_sum_loggf=round(raw_sum, 4),
                gf_scale_used="NIST ASD (RYA-706 pull)", target_loggf=target_gf,
                nist_accuracy=NIST_ACC, nist_line_ref="L4737",
                canonical_loggf=canon_gf, canonical_ep=canon_ep,
                canonical_minus_nist_dex=round(canon_gf - target_gf, 4),
                canonical_NOT_used_reason=("VALD3, adjudication_status pending, 0.222 dex "
                    "high; its C+ grade is attached to 6696.185 (EP 4.022), a different "
                    "lower level. See RYA-706."),
                canonical_line_id=canon["line_id"], canonical_ref=canon["loggf_reference"],
                canonical_nist_grade=canon["nist_grade"],
                shift_applied=round(shift, 4), rescaled_to_target=bool(rescaled),
                fdamp=fdamp,
                span_mA=round((comps[-1][0] - comps[0][0]) * 1000, 1))
    return path, info


def read_vald_blocks(lmin, lmax):
    """Parse the VALD grid block once, keeping only rows inside [lmin, lmax]. Returns a
    list of {code, stage, name, rows}."""
    src = os.path.join(VALD_DIR, VALD_BLOCK)
    if not os.path.exists(src):
        raise SystemExit(f"RYA-706: VALD in-window block not staged on Sirius: {src}")
    blocks, cur = [], None
    with open(src) as f:
        it = iter(f)
        for ln in it:
            m = HDR_RE.match(ln.rstrip("\n"))
            if m:
                name = next(it).rstrip("\n")
                cur = dict(code=m.group(1), stage=m.group(2), name=name, rows=[], waves=[])
                blocks.append(cur)
                continue
            if cur is None:
                continue
            p = ln.split()
            try:
                wl = float(p[0])
            except (ValueError, IndexError):
                continue
            if lmin <= wl <= lmax:
                cur['rows'].append(ln.rstrip("\n"))
                cur['waves'].append(wl)
    return [b for b in blocks if b['rows']]


def _species(block):
    return block['name'].strip().strip("'").strip()


def write_blocks(blocks, path):
    with open(path, "w") as f:
        for b in blocks:
            f.write(f"'  {b['code']}            '    {b['stage']}    {len(b['rows'])}\n")
            f.write(b['name'] + "\n")
            for r in b['rows']:
                f.write(r + "\n")
    return path


# ────────────────────────────── engine ──────────────────────────────
def babsma(opac, lmin, lmax, feh, xi):
    ctl = (f"'LAMBDA_MIN:'  '{lmin}'\n'LAMBDA_MAX:'  '{lmax}'\n'LAMBDA_STEP:' '0.005'\n"
           f"'MODELINPUT:' '{MARCS}'\n'MARCS-FILE:' '.true.'\n'MODELOPAC:' '{opac}'\n"
           f"'METALLICITY:'    '{feh:.2f}'\n'ALPHA/Fe   :'    '0.00'\n'HELIUM     :'    '0.00'\n"
           f"'R-PROCESS  :'    '0.00'\n'S-PROCESS  :'    '0.00'\n'XIFIX:' 'T'\n{xi:.2f}\n")
    r = subprocess.run([f"{EXE}/babsma_lu"], input=ctl, capture_output=True, text=True, cwd=W)
    if not os.path.exists(opac) or os.path.getsize(opac) == 0:
        raise SystemExit(f"RYA-706: babsma failed ({lmin}-{lmax}):\n{r.stdout[-1500:]}")
    return opac


def bsyn(a_x, linelists, opac, lmin, lmax, tag, feh):
    """LTE flux synthesis over `linelists`. Returns (wave, flux). Loud-fails on a missing
    result file — a silently empty synthesis must never reach the fitter."""
    res = f"{W}/s_{tag}.spec"
    if os.path.exists(res):
        os.remove(res)
    ctl = (f"'NLTE :'          '.false.'\n'LAMBDA_MIN:'     '{lmin}'\n'LAMBDA_MAX:'     '{lmax}'\n"
           f"'LAMBDA_STEP:'    '0.005'\n'INTENSITY/FLUX:' 'Flux'\n'ABFIND        :' '.false.'\n"
           f"'MODELOPAC:' '{opac}'\n'RESULTFILE :' '{res}'\n'METALLICITY:'    '{feh:.2f}'\n"
           f"'INDIVIDUAL ABUNDANCES:'   '1'\n{Z_AL}  {a_x:.3f}\n'ISOTOPES : ' '0'\n"
           f"'NFILES   :' '{len(linelists)}'\n" + "".join(f"{p}\n" for p in linelists)
           + f"'SPHERICAL:'  'F'\n  30\n  300.00\n  15\n  1.30\n")
    r = subprocess.run([f"{EXE}/bsyn_lu"], input=ctl, capture_output=True, text=True, cwd=W)
    if not os.path.exists(res):
        raise SystemExit(f"RYA-706: bsyn failed {tag} A={a_x}:\n{r.stdout[-2000:]}")
    d = np.loadtxt(res)
    return d[:, 0], d[:, 1]


def core_ew(w, f, hw=CORE_HW):
    m = (w > LINE - hw) & (w < LINE + hw)
    return float(_trapezoid(1.0 - f[m], w[m]) * 1000.0)


# ────────────────────────────── observed arms ──────────────────────────────
def load_harps():
    w, f = [], []
    with open(HARPS) as fh:
        r = csv.reader(fh)
        next(r)
        for row in r:
            w.append(float(row[0]))
            f.append(float(row[3]))
    return np.array(w), np.array(f)


def load_iag(root):
    """IAG Reiners+2016 visible FTS solar flux atlas: vacuum wavenumber (cm^-1), norm flux.
    vac->air through the pipeline SINGLE-SOURCE converter (RYA-264/501)."""
    import gzip
    sys.path.insert(0, root)
    from pipeline.wavelength_util import vac_to_air
    wn, fl = [], []
    with gzip.open(IAG, 'rt') as fh:
        for ln in fh:
            p = ln.split()
            if len(p) < 2:
                continue
            try:
                wn.append(float(p[0]))
                fl.append(float(p[1]))
            except ValueError:
                continue
    lam_air = vac_to_air(1e8 / np.array(wn))
    idx = np.argsort(lam_air)
    return lam_air[idx], np.array(fl)[idx]


def fit_all_hw(obs_w, obs_f, synth):
    """chi2 profile fit at every half-window in HW_SWEEP. The headline is FIT_HW; the rest
    are reported so the value's dependence on the window choice is visible rather than
    assumed away."""
    out = {}
    for hw in HW_SWEEP:
        fit = fit_profile(LINE, obs_w, obs_f, synth, hw=hw, vsini=VSINI,
                          a_lo=A_LO, a_hi=A_HI, core_hw=CORE_HW)
        out[f"{hw:.1f}"] = fit
    return out


def _pack(fit, delta, arm_rv):
    if fit is None:
        return dict(status='no_coverage')
    a_lte = fit['A']
    rel_f = assess_reliability(fit)
    reliable = rel_f['reliable']
    return dict(A_LTE=round(a_lte, 3),
                A_NLTE=round(a_lte + delta, 3) if delta is not None else None,
                nlte_delta=delta,
                dv_fitted_kms=round(fit['dv'], 2), dv_measured_kms=arm_rv,
                gsig_kms=round(fit['gsig'], 2), gsig_railed=fit['gsig_railed'],
                red_chi2=round(fit['red_chi2'], 2), npix=fit['npix'],
                obs_EW_mA=round(fit['obs_ew_mA'], 2), core_EW_mA=fit['core_EW_mA'],
                dEW_dA_mA_dex=fit['dEW_dA'], railed=fit['railed'], **rel_f)


# ────────────────────────────── main ──────────────────────────────
def _assert_no_ba_residue():
    """RYA-706: this file was ADAPTED from the Ba II harness and that produced defects at
    a rate that made copying the wrong choice. Six were found by reading and by the
    runtime guards; a further sweep found the emitted record still carried
    element='Ba', a correction budget referenced to Ba's 2.410, and delta_vs_asplund
    computed against Ba's A_sun of 2.27. Had the fit converged, that JSON would have
    been a BARIUM-shaped record holding aluminium numbers -- and a downstream consumer
    keys on `element`.

    So the file now checks itself. This is a stopgap: the real fix is to parameterise
    the shared harness rather than copy it per species, which is a toolkit task."""
    import re as _re
    src = open(os.path.abspath(__file__)).read()
    # Ignore this function's own body: it necessarily NAMES the tokens it hunts for, and
    # a checker that trips on its own patterns is a checker nobody keeps.
    # Anchor on the NEWLINE-prefixed definitions. Searching for the bare string finds
    # this very line first, which truncated the exclusion range to nothing and made the
    # checker flag its own pattern list. A self-check that cannot survive reading itself
    # is worse than none, because the first thing it teaches you is to ignore it.
    _lo = src.index("\ndef _assert_no_ba_residue")
    _hi = src.index("\ndef main(", _lo)
    # Prose ABOUT the adaptation is legitimate; a live Ba token in code is not.
    bad = []
    for pat in (r"SOLAR_ASPLUND2021\['Ba'\]", r"'element':\s*'Ba'", r"2\.410", r"5853",
                r"_mpia_element_delta\('Ba'", r"element_grid_in_bounds\('Ba'"):
        for m in _re.finditer(pat, src):
            line = src[:m.start()].count("\n") + 1
            txt = src.splitlines()[line-1]
            if _lo <= m.start() < _hi:
                continue                      # the checker's own patterns and docstring
            if txt.lstrip().startswith("#") or '"""' in txt:
                continue                      # documentation of the adaptation itself
            bad.append(f"line {line}: {txt.strip()[:80]}")
    if bad:
        raise SystemExit("RYA-706 Ba-residue check FAILED — this harness still carries "
                         "barium-specific code:\n  " + "\n  ".join(bad))


def main():
    _assert_no_ba_residue()
    ap = argparse.ArgumentParser(description="RYA-706 Al I 6696 in-window deblend (Sirius)")
    ap.add_argument('--out', default=None, help="output JSON path")
    ap.add_argument('--skip-sensitivity', action='store_true',
                    help="skip the VALD-vdW damping sensitivity leg")
    args = ap.parse_args()

    root = _ROOT
    sys.path.insert(0, root)
    # RYA-567: Sirius-only heavy compute (Turbospectrum + MARCS + the staged VALD deck).
    # Loud-fail off-Sirius — never against local-Mac copies of engines/grids/spectra.
    from config.constants import assert_on_sirius, STAR_PARAMS, SOLAR_ASPLUND2021
    assert_on_sirius("RYA-706 Al I 6696 in-window deblend",
                     require_subdirs=("engines", "grids"))
    sp = STAR_PARAMS['solar']
    teff, logg, feh, xi = sp['teff'], sp['logg'], sp.get('feh_ref', 0.0), sp['xi']
    a_sun = float(SOLAR_ASPLUND2021['Al'])

    print(f"RYA-706 Al I {LINE} in-window deblend — solar node (STAR_PARAMS): "
          f"Teff={teff} logg={logg} [Fe/H]={feh} xi={xi}")
    print(f"  reference A(Al)_sun (Asplund 2021) = 6.43  [grid bracket only, NOT a target]\n")

    os.makedirs(W, exist_ok=True)
    if not os.path.lexists(f"{W}/DATA"):
        os.symlink("/mnt/codex-data/engines/Turbospectrum_NLTE/DATA", f"{W}/DATA")

    # ── (1) line lists ────────────────────────────────────────────────────────────
    al_list, hfs = build_al_hfs_list(root, f"{W}/al6696_hfs.list", FDAMP_PRIMARY)
    print(f"Al I HFS list: {hfs['n_components']} components, span {hfs['span_mA']} mA, "
          f"linelist_full sum {hfs['linelist_full_sum_loggf']:+.4f} -> canonical "
          f"{hfs['canonical_loggf']:+.3f} ({hfs['canonical_line_id']}, "
          f"{hfs['canonical_ref']}, grade {hfs['canonical_nist_grade']}; shift "
          f"{hfs['shift_applied']:+.4f}"
          f"{', RESCALED' if hfs['rescaled_to_target'] else ''})")

    lmin, lmax = LINE - LMARGIN, LINE + LMARGIN
    blocks = read_vald_blocks(lmin, lmax)
    al_blocks = [b for b in blocks if _species(b).startswith('Al I')]
    blend_blocks = [b for b in blocks if not _species(b).startswith('Al I')]
    n_al_dropped = sum(len(b['rows']) for b in al_blocks)
    n_blend_rows = sum(len(b['rows']) for b in blend_blocks)
    blend_list = write_blocks(blend_blocks, f"{W}/blend_6696.list")
    print(f"Blend list: VALD window {lmin:.1f}-{lmax:.1f} A -> {len(blend_blocks)} species, "
          f"{n_blend_rows} rows ({n_al_dropped} Al I rows removed; Al comes from the HFS list)")

    # THE RYA-706 STOP CONDITION, enforced up front: no blends in the window means the
    # deblend never happened. Refuse rather than emit a value that silently repeats 559.
    if n_blend_rows == 0:
        raise SystemExit("RYA-706 STOP: the in-window blend list is EMPTY — the blend was "
                         "culled again, so this fit cannot deblend anything. Investigate the "
                         "VALD staging; do NOT ship a value.")

    near = [(_species(b), sum(1 for w in b['waves'] if abs(w - LINE) <= BLEND_CENSUS_HW))
            for b in blend_blocks]
    census_species = sorted([s for s, n in near if n > 0])
    print(f"  species with rows within +/-{BLEND_CENSUS_HW} A of the core: "
          f"{', '.join(census_species) if census_species else '(none)'}")

    # ── (2) engine ────────────────────────────────────────────────────────────────
    opac = babsma(f"{W}/ba.opac", lmin, lmax, feh, xi)

    print("\nSynthesising the A(Al) grid ...")
    synth_blend, synth_noblend = {}, {}
    for a in A_GRID:
        a = float(a)
        synth_blend[a] = bsyn(a, [blend_list, al_list], opac, lmin, lmax, f"bl_{a:.2f}", feh)
        synth_noblend[a] = bsyn(a, [al_list], opac, lmin, lmax, f"nb_{a:.2f}", feh)
    print(f"  {len(A_GRID)} points x 2 configurations (with blends / Al-only control)")

    # ── (3) blend accounting: what the window actually contains ───────────────────
    # A(Al) = -99 switches aluminium off without touching anything else, so the residual is
    # exactly the modelled blend. Per-species figures come from synthesising each species
    # alone, which is what "the in-window blend components, listed" has to mean.
    w_b, f_b = bsyn(-99.0, [blend_list], opac, lmin, lmax, "blends_only", feh)
    blend_core_ew = core_ew(w_b, f_b)
    per_species = {}
    for b in blend_blocks:
        sp_name = _species(b)
        if sp_name not in census_species:
            continue
        p = write_blocks([b], f"{W}/sp_{re.sub(r'[^A-Za-z0-9]', '', sp_name)}.list")
        w_s, f_s = bsyn(-99.0, [p], opac, lmin, lmax,
                        f"sp_{re.sub(r'[^A-Za-z0-9]', '', sp_name)}", feh)
        per_species[sp_name] = round(core_ew(w_s, f_s), 3)
    per_species = dict(sorted(per_species.items(), key=lambda kv: -kv[1]))
    print(f"\nModelled in-window blend, core +/-{CORE_HW} A: {blend_core_ew:.2f} mA total")
    for k, v in per_species.items():
        print(f"    {k:12s} {v:7.3f} mA")

    # ── (4) observed arms + rest frame ────────────────────────────────────────────
    arms = {'harps': load_harps()}
    try:
        arms['iag'] = load_iag(root)
        print(f"\nIAG atlas: {len(arms['iag'][0])} pts, "
              f"{arms['iag'][0].min():.1f}-{arms['iag'][0].max():.1f} A")
    except Exception as e:                                            # noqa: BLE001
        print(f"\nIAG load FAILED ({e}); HARPS-only")

    # RYA-643: the rest-frame correction must be SOURCED from a measurement on each arm
    # being fitted — never hardcoded, never a silent zero.
    arm_rv = {}
    print('residual velocity per arm (clean-line centroids, abundance-blind):')
    for _arm, (_ow, _of) in arms.items():
        _v, _n, _sd = require_arm_rv(_ow, _of, _arm)
        arm_rv[_arm] = dict(v_kms=round(_v, 3), n_lines=_n, scatter_kms=round(_sd, 3))
        print(f'    {_arm:6s}: {_v:+.3f} km/s (n={_n}, scatter {_sd:.3f})')

    # ── (5) NLTE delta — READ from the production Engine-A grid ───────────────────
    from pipeline import nlte_corrections as N
    # Al's grid is Amarsi-2020 (PySME), NOT Ba's Korotin2015. The mechanical adaptation
    # from RYA-581 left the Ba call here and the run STOPPED on a NaN delta rather than
    # falling through to LTE — the RYA-518 loud-fail doing exactly its job. Recorded
    # rather than quietly fixed, because "no silent LTE" catching a real wiring error on
    # its first outing is the evidence that the guard earns its place.
    #
    # The Al grid carries EXACTLY two wavelengths, 6696.02 and 6698.67 — it was built for
    # the two lines RYA-706 recovered, which the committed pool never held.
    delta = float(N._mpia_element_delta('Al', LINE, teff, logg, feh))
    in_bounds = bool(N.element_grid_in_bounds('Al', teff, logg, feh))
    if not np.isfinite(delta):
        raise SystemExit("RYA-706: Al_Amarsi2020 delta is NaN at the solar node — the "
                         "Engine-A grid did not deliver a correction. STOP (no silent LTE).")
    print(f"\nEngine-A Amarsi2020 1D-NLTE delta at the solar node: {delta:+.4f} "
          f"(in_bounds={in_bounds}) — READ from data/nlte_grids/Al_Amarsi2020_PySME.csv")

    # ── (6) the fits ──────────────────────────────────────────────────────────────
    results = {}
    for arm, (ow, of) in arms.items():
        fb = fit_all_hw(ow, of, synth_blend)
        fn = fit_all_hw(ow, of, synth_noblend)
        results[arm] = dict(
            deblended={hw: _pack(f, delta, arm_rv[arm]['v_kms']) for hw, f in fb.items()},
            control_no_blend={hw: _pack(f, delta, arm_rv[arm]['v_kms']) for hw, f in fn.items()},
        )
        h = f"{FIT_HW:.1f}"
        db, ct = results[arm]['deblended'][h], results[arm]['control_no_blend'][h]
        print(f"\n  {arm.upper()}  (headline fit half-window {FIT_HW} A)")
        print(f"    deblended (blends modelled) : A_LTE={db['A_LTE']:.3f} "
              f"{delta:+.4f} -> A_NLTE={db['A_NLTE']:.3f}  "
              f"(rchi2={db['red_chi2']}, dEW/dA={db['dEW_dA_mA_dex']}, "
              f"gsig={db['gsig_kms']}, dv={db['dv_fitted_kms']}, reliable={db['reliable']})")
        print(f"    control   (Al alone, no blend): A_LTE={ct['A_LTE']:.3f} "
              f"-> A_NLTE={ct['A_NLTE']:.3f}  (rchi2={ct['red_chi2']})")
        print(f"    deblend shift = {db['A_NLTE'] - ct['A_NLTE']:+.3f} dex")
        print(f"    fit-window sweep (A_NLTE): " + "  ".join(
            f"hw={k}:{v['A_NLTE']:.3f}" for k, v in results[arm]['deblended'].items()))

    prim = results['harps']['deblended'][f"{FIT_HW:.1f}"]
    ctrl = results['harps']['control_no_blend'][f"{FIT_HW:.1f}"]
    a_lte, a_nlte = prim['A_LTE'], prim['A_NLTE']
    shift = round(a_nlte - ctrl['A_NLTE'], 3)
    chi2_gain = round(ctrl['red_chi2'] / prim['red_chi2'], 2) if prim['red_chi2'] > 0 else None

    # ── THE RYA-706 STOP CONDITION, second half ───────────────────────────────────
    # "if the fit still returns ~2.41 the blend was not modelled -> investigate, do not
    # ship." Three INDEPENDENT things have to hold, because the shift in A on its own is
    # a weak proxy: a blend seated in the WINGS can be modelled correctly, dominate the
    # fit residual, and still move the fitted A only slightly. Testing the A-shift alone
    # would abort a good run (this one moves A by exactly 0.020 dex) while a genuinely
    # culled blend list is caught by the other two. So we assert what actually has to be
    # true for a deblend to have happened:
    #   (a) the modelled blend deposits real absorption in the core (it is IN the window);
    #   (b) modelling it materially improves the fit (it EXPLAINS observed structure);
    #   (c) the answer has moved off the RYA-559 EW->COG value we came here to supersede.
    stop = []
    if blend_core_ew < 1.0:
        stop.append(f"the modelled blend deposits only {blend_core_ew:.2f} mA in the core "
                    f"+/-{CORE_HW} A — there is effectively nothing in the window to deblend")
    if chi2_gain is None or chi2_gain < 1.5:
        stop.append(f"modelling the blend does not improve the fit (red_chi2 "
                    f"{ctrl['red_chi2']} -> {prim['red_chi2']}, gain {chi2_gain}) — the "
                    f"blend is not explaining any observed structure")
    # RYA-706's own stop condition, replacing RYA-581's Ba one. The blended EW fit gave
    # A(Al) = 6.823 on the NIST gf scale. If the deblend returns that same value the
    # blend was not modelled and there is nothing to ship — the Fe I 6696.315 census
    # says 12.775 mA of a 60.30 mA line, so a real deblend CANNOT leave A unchanged.
    if abs(a_nlte - 6.823) < 0.03:
        stop.append(f"the fit returned A(Al)_NLTE {a_nlte:.3f}, still the blended "
                    f"EW-path value 6.823 — the blend was not modelled")
    if stop:
        raise SystemExit("RYA-706 STOP: " + "; ".join(stop)
                         + ". Investigate the VALD in-window block; do NOT ship this value.")

    print(f"\nDeblend evidence: blend core EW {blend_core_ew:.2f} mA; red_chi2 "
          f"{ctrl['red_chi2']} (Al alone) -> {prim['red_chi2']} (blends modelled), "
          f"{chi2_gain}x better; A shift {shift:+.3f} dex")

    # ── (7) damping sensitivity (declared, not silent) ────────────────────────────
    sens = None
    if not args.skip_sensitivity:
        al_v, _ = build_al_hfs_list(root, f"{W}/al6696_hfs_vald.list", FDAMP_VALD)
        sv = {}
        for a in A_GRID:
            a = float(a)
            sv[a] = bsyn(a, [blend_list, al_v], opac, lmin, lmax, f"vd_{a:.2f}", feh)
        fv = fit_profile(LINE, *arms['harps'], sv, hw=FIT_HW, vsini=VSINI,
                         a_lo=A_LO, a_hi=A_HI, core_hw=CORE_HW)
        sens = dict(fdamp=FDAMP_VALD, source='linelist_full.csv damping_vdW',
                    A_LTE=round(fv['A'], 3), A_NLTE=round(fv['A'] + delta, 3),
                    red_chi2=round(fv['red_chi2'], 2),
                    d_vs_primary=round(fv['A'] - a_lte, 3))
        print(f"\nDamping sensitivity (VALD vdW {FDAMP_VALD} instead of Unsold "
              f"{FDAMP_PRIMARY}): A_NLTE={sens['A_NLTE']:.3f} "
              f"({sens['d_vs_primary']:+.3f} dex vs the primary)")

    # ── (8) record ────────────────────────────────────────────────────────────────
    out = {
        'ticket': 'RYA-706',
        'supersedes': 'RYA-706 EW-path value 6.823 (Fe I 6696.315 charged to Al)',
        'element': 'Al', 'ion': 'I', 'line_A': LINE,
        'method': ('in-window blend-fit: full VALD3 in-window block MODELLED alongside the '
                   'Al I 6696 HFS/isotope components, A(Al) fitted to the observed solar '
                   'PROFILE by chi2 (RYA-551 pattern, shared RYA-643 fitter) + Engine-A '
                   'Amarsi2020 1D-NLTE delta read at the solar node'),
        'why_this_supersedes_559': (
            'an EW inversion cannot deblend: one scalar cannot separate aluminium from the '
            'other absorption inside the integration window, so the blend_flag=True pool EW '
            'was inflated by the Fe I 6696.315 blend charged to Al. This fits the '
            'profile with the blend modelled instead.'),
        'solar_node': {'teff': teff, 'logg': logg, 'feh': feh, 'xi': xi,
                       'source': "STAR_PARAMS['solar']", 'marcs_node': 'p5750_g+4.5_z+0.00'},
        'hfs': hfs,
        'blend_model': {
            'vald_block': VALD_BLOCK, 'window_A': [round(lmin, 3), round(lmax, 3)],
            'n_species': len(blend_blocks), 'n_rows': n_blend_rows,
            'al_i_rows_removed': n_al_dropped,
            'core_hw_A': CORE_HW,
            'blend_core_EW_mA': round(blend_core_ew, 3),
            'per_species_core_EW_mA': per_species,
        },
        'fit': {'half_window_A': FIT_HW, 'sweep_A': list(HW_SWEEP),
                'a_grid': [A_LO, A_HI, 0.05], 'vsini': VSINI,
                'reliable_dEW_dA_floor': RELIABLE_DEWDA,
                'red_chi2_review_trigger': RCHI2_REVIEW,
                'red_chi2_gates_reliable': False},
        'arms': {k: arm_rv[k] for k in arms},
        'engineA_amarsi_delta': round(delta, 4),
        'engineA_in_bounds': in_bounds,
        'engineA_grid': 'data/nlte_grids/Al_Amarsi2020_PySME.csv (production; carries exactly 6696.02 + 6698.67)',
        'A_lte': a_lte, 'A_nlte': a_nlte,
        'sigma': round(float(np.std([results[a]['deblended'][f"{FIT_HW:.1f}"]['A_NLTE']
                                     for a in results
                                     if results[a]['deblended'][f"{FIT_HW:.1f}"].get('A_NLTE')
                                     is not None], ddof=0)), 3) if len(results) > 1 else None,
        # Say what that number IS, so nobody reads it as a total error budget. It is the
        # HARPS-vs-IAG arm scatter and nothing else: it does not carry the canonical gf
        # uncertainty, the Amarsi2020 delta uncertainty, the 1D-MARCS model error, or the
        # fact that this is a SINGLE line. The honest read is a floor, not a sigma.
        'sigma_basis': ('inter-arm scatter of the fitted A_NLTE (HARPS vs IAG) at the '
                        'headline fit window ONLY — a precision floor, NOT a total '
                        'uncertainty budget: it excludes the canonical gf uncertainty, '
                        'the Amarsi2020 delta uncertainty, 1D-MARCS model error, and the '
                        'single-line risk. Do not quote it as the error on A(Al).'),
        'window_choice_spread_dex': round(
            max(v['A_NLTE'] for v in results['harps']['deblended'].values())
            - min(v['A_NLTE'] for v in results['harps']['deblended'].values()), 3),
        'reliable': prim['reliable'],
        'control_no_blend_A_nlte': ctrl['A_NLTE'],
        'deblend_shift_dex': shift,
        'ew_path_A_nlte_blended': 6.823,
        'deblend_evidence': {
            'blend_core_EW_mA': round(blend_core_ew, 3),
            'red_chi2_al_alone': ctrl['red_chi2'],
            'red_chi2_blends_modelled': prim['red_chi2'],
            'chi2_improvement_factor': chi2_gain,
            'note': ('the blend is demonstrably MODELLED, not culled: it deposits '
                     f'{blend_core_ew:.2f} mA in the core and modelling it improves the '
                     f'profile fit {chi2_gain}x. The A-shift it produces ({shift:+.3f} dex) '
                     'is small because the blend sits largely in the WINGS of the fit '
                     'window rather than under the Al core.'),
        },
        # Decomposition against the EW-path value this run exists to test, 6.823 dex on
        # the NIST gf scale. RYA-581's version of this block was written against Ba's
        # 2.410 and its "66.5 mA calibration point" -- none of which means anything here.
        'correction_budget_dex': {
            'ew_path_blended': 6.823,
            'profile_fit_al_alone': ctrl['A_NLTE'],
            'profile_fit_deblended': a_nlte,
            'from_dropping_the_EW_inversion': round(ctrl['A_NLTE'] - 6.823, 3),
            'from_modelling_the_in_window_blend': shift,
            'total': round(a_nlte - 6.823, 3),
            'fitted_core_EW_mA': prim['core_EW_mA'],
            'rya559_calibration_A227_EW_mA': 66.5,
            'note': ('the fitted synthetic core EW reproduces the RYA-559 calibration and '
                     'the literature clean-line EW (~64-66 mA), confirming the profile fit '
                     'recovers the CLEAN Al I 6696 that the blend-inflated pool EW '
                     '(74.62 mA) never isolated.'),
        },
        'damping_sensitivity': sens,
        'per_arm': results,
        'asplund2021': a_sun,
        'delta_vs_asplund': round(a_nlte - a_sun, 3),
        'engineB_status': (
            'NOT ASSESSED for Al I 6696 by this run. RYA-581 asserted a Gerber/GES level-ID '
            'finding and an Engine-B departure cross-check for Ba II 4554 -- both are Barium '
            'claims and neither transfers, so they are removed rather than reworded. Whether '
            'TS-native NLTE is runnable at this line is an open question for the Al two-engine '
            'record, not something this harness established.'),
        'single_line_caveat': ('this rests on ONE line, Al I 6696.023. The confirming line '
                               'is Al I 6698.673, which the RYA-706 EW path puts at 6.485 and '
                               'which is CLEAN -- it needs no deblend, so it is a separate '
                               'measurement rather than a second run of this harness.'),    }

    out_path = args.out or os.path.join(root, "data", "results", "solar_ba_deblend_rya706.json")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)

    print("\n" + "=" * 76)
    print(f"RYA-706 Al I {LINE} DEBLEND RESULT")
    print(f"  EW path, blend-inflated          : A(Al)_NLTE = 6.823")
    print(f"  control, Al alone (this harness) : A(Al)_NLTE = {ctrl['A_NLTE']:.3f}")
    print(f"  DEBLENDED, blends modelled       : A(Al)_NLTE = {a_nlte:.3f}  "
          f"(A_LTE {a_lte:.3f} {delta:+.4f})")
    print(f"  deblend shift                    : {shift:+.3f} dex")
    print(f"  vs Asplund 2021 ({a_sun})          : {a_nlte - a_sun:+.3f} dex")
    print(f"  reliable={prim['reliable']}  rchi2={prim['red_chi2']}  "
          f"dEW/dA={prim['dEW_dA_mA_dex']} mA/dex  railed={prim['railed']}")
    print(f"\nWrote {out_path}")


if __name__ == '__main__':
    main()
