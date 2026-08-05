#!/usr/bin/env python3
"""
RYA-564 — Cobalt (Co I) red-line HFS-resolved synthesis measurement (RUNS ON SIRIUS).

Replaces the untrusted blue-edge Co I 3845 value (A(Co)=6.128, +1.188 vs Asplund;
Kitt Peak SNR~24, chi2r~3100 — the same artifact class as the Sr I +2.13 the RYA-524
audit was built to catch) with a real measurement on clean RED Co I lines.

METHOD (three separated concerns — the RYA-559 Ba pattern, per-line):

  1. A_LTE — Turbospectrum LTE flux fit against the observed solar spectrum, with the
     full VALD in-window block supplying the BLENDS and the GES NLTE line list supplying
     the Co I HFS COMPONENTS (Co is hyperfine-split: a single-profile EW fit cannot reach
     it — the RYA-354/466 finding). Two arms: HARPS reflected-solar (primary) + IAG
     Reiners-2016 FTS atlas (cross-check). Machinery reused verbatim from
     scripts/rya551_sr2_synth_sirius.py (babsma/bsyn/broaden/local_renorm/fit_line).

  2. delta_NLTE — the PER-LINE 1D-NLTE correction READ from the Sirius-staged Gerber
     grid through the RYA-534-validated deck (scripts/ts_gerber_gate.py method: bsyn NLTE
     at the solar node vs an LTE curve-of-growth, delta = A_sun - A*). NEVER fitted. The
     deck RAISES if departures don't engage (the RYA-533 silent-departure=1 trap), so a
     line can never be silently reported as LTE.

  3. A_NLTE = A_LTE + delta_NLTE, per line.

WHY NOT 5000/5013 ALONE (the two lines RYA-534 gated): they were chosen for the NLTE
gate because BOTH their GES model-atom levels are identified — not because they are
measurable. At the solar Co abundance they are sub-mA features (EW 0.67 / 1.08 mA in the
RYA-534 gate output) and carry essentially no abundance sensitivity. They are measured
here anyway, as DIAGNOSTIC lines, so the record shows that by evidence rather than by
assertion. The ticket explicitly permits adding clean red Co I lines within reach; the
reportable set is the strong HFS-resolved red lines whose GES component sum reconciles
with the canonical_gf single source AND whose model-atom levels are identified.

SINGLE SOURCE: every log gf is asserted against data/linelists/canonical_gf.csv at
runtime (loud-fail on drift). Where the GES HFS component sum disagrees with the
canonical total, the components are rescaled UNIFORMLY onto the canonical total (the HFS
pattern is preserved, the scale comes from the single source) and the line is FLAGGED —
never silently synthesised on an off-SSOT gf.

VALIDATE-DON'T-TUNE: nothing is adjusted toward Asplund 2021 A(Co)=4.94. A red-line value
that still sits high is a finding, not a PASS. If no line clears the reliability floor,
Co reports NO VALUE (owed) — falling back to the 3845 artifact is forbidden.
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
from pipeline._numcompat import trapezoid as _trapezoid  # numpy>=2 removed np.trapz (RYA-313)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from solar_profile_fit import (CLIGHT, broaden,  # noqa: E402,F401
                               fit_profile, local_renorm, measure_arm_rv,
                               require_arm_rv)

EXE    = "/mnt/codex-data/engines/Turbospectrum_NLTE/exec-gf"
INTERP = "/mnt/codex-data/engines/Turbospectrum_NLTE/interpolator/interpol_modeles_nlte"
GES    = ("/mnt/codex-data/engines/Turbospectrum_NLTE/COM/linelists/"
          "nlte_ges_linelist_jmg17feb2022_I_II")
GT     = "/mnt/codex-data/grids/nlte/gerber_ts"
MARCS  = ("/mnt/codex-data/grids/model_atmospheres/marcs_standard_comp/marcs_standard_comp/"
          "p5750_g+4.5_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod")
VALD_DIR = "/mnt/codex-data/engines/TSFitPy/input_files/linelists/linelist_vald"
W      = "/mnt/codex-data/codex/_solar_tmp/rya564_work"
HARPS  = "/mnt/codex-data/codex/_solar_tmp/solar_normalized_harps.csv"
IAG    = "/mnt/codex-data/solar_reference/iag_reiners2016/spvis.dat.gz"

Z_CO  = 27
ATOM  = 'atom.co247qm'                              # RYA-534 Co model atom (Gerber 2023)
AUX   = 'auxData_CO_MARCS_Nov-14-2024.dat'
GRID  = 'NLTEgrid4TS_CO_MARCS_Nov-14-2024.bin'
TREF, LOGGREF, ZREF = "5750", "+4.50", "+0.00"      # nearest MARCS node (RYA-534 deck)

XI     = 1.0            # solar microturbulence (matches the t01 MARCS node)
VSINI  = 1.8            # solar vsini (STAR_PARAMS)

# RYA-534 published/banked Co anchor — used ONLY to CHECK the per-line deltas this run
# computes (cross-engine sanity), never to set them.
ANCHOR_DELTA = 0.10     # Bergemann+2010 MNRAS 401 1334 (RYA-534 TS-Gerber median +0.099)
ANCHOR_TOL   = 0.12     # the RYA-534 gate tolerance for Co

# Reliability floor — identical to RYA-551/560 so Co is judged on the same bar as Sr/Zr.
RELIABLE_DEWDA = 40.0   # mA/dex core-EW sensitivity
RELIABLE_RCHI2 = 60.0   # in-window fit quality ceiling (sigma_flux=0.01 floor inflates rchi2)

LMARGIN = 7.0                                   # synth half-window (A)
A_GRID  = np.round(np.arange(4.20, 5.76, 0.05), 3)   # A(Co) trial grid (brackets 4.94)
A_LO, A_HI = float(A_GRID.min()), float(A_GRID.max())
COG_OFFSETS = [-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3]  # LTE COG for the NLTE delta (RYA-534)

# Candidate red Co I lines. `loggf` = the CANONICAL total from data/linelists/canonical_gf.csv
# (verified at runtime); `role`:
#   primary    — gf adjudication 'agreed' + GES HFS sum reconciles + both GES levels identified
#   secondary  — measurable + levels identified, but gf is 'pending'/'single_source' or the GES
#                HFS sum needed rescaling onto the canonical total (flagged, not silent)
#   diagnostic — measured for the record only (the RYA-534 gate lines; sub-mA, no sensitivity)
LINES = {
    4813.476: dict(ep=3.216, loggf=+0.120, gf_ref='VALD3', adj='agreed',
                   vald='vald-4800-5300-for-grid.list', fit_hw=0.9, role='primary'),
    5212.688: dict(ep=3.514, loggf=-0.110, gf_ref='VALD3', adj='agreed',
                   vald='vald-4800-5300-for-grid.list', fit_hw=0.9, role='primary'),
    5352.040: dict(ep=3.576, loggf=+0.060, gf_ref='VALD3', adj='agreed',
                   vald='vald-5300-5800-for-grid.list', fit_hw=0.9, role='primary'),
    5647.234: dict(ep=2.280, loggf=-1.560, gf_ref='VALD3', adj='agreed',
                   vald='vald-5300-5800-for-grid.list', fit_hw=0.8, role='primary'),
    6454.994: dict(ep=3.632, loggf=-0.250, gf_ref='VALD3', adj='agreed',
                   vald='vald-6300-6800-for-grid.list', fit_hw=0.9, role='primary'),
    6632.439: dict(ep=2.280, loggf=-2.000, gf_ref='VALD3', adj='agreed',
                   vald='vald-6300-6800-for-grid.list', fit_hw=0.8, role='primary'),
    6771.034: dict(ep=1.883, loggf=-1.970, gf_ref='VALD3', adj='agreed',
                   vald='vald-6300-6800-for-grid.list', fit_hw=0.9, role='primary'),
    5301.041: dict(ep=1.710, loggf=-2.000, gf_ref='VALD3', adj='pending',
                   vald='vald-5300-5800-for-grid.list', fit_hw=0.9, role='secondary'),
    5331.453: dict(ep=1.785, loggf=-1.960, gf_ref='VALD3', adj='pending',
                   vald='vald-5300-5800-for-grid.list', fit_hw=0.9, role='secondary'),
    5483.354: dict(ep=1.710, loggf=-1.490, gf_ref='VALD3', adj='pending',
                   vald='vald-5300-5800-for-grid.list', fit_hw=1.0, role='secondary'),
    5915.550: dict(ep=2.137, loggf=-2.030, gf_ref='1999ApJS (Nitz+)', adj='pending',
                   vald='vald-5800-6300-for-grid.list', fit_hw=0.8, role='secondary'),
    6188.997: dict(ep=1.710, loggf=-2.450, gf_ref='1982ApJ (Cardon+)', adj='single_source',
                   vald='vald-5800-6300-for-grid.list', fit_hw=0.8, role='secondary'),
    5000.875: dict(ep=3.252, loggf=-1.965, gf_ref='K08', adj='single_source',
                   vald='vald-4800-5300-for-grid.list', fit_hw=0.6, role='diagnostic'),
    5013.324: dict(ep=3.568, loggf=-1.454, gf_ref='K08', adj='single_source',
                   vald='vald-4800-5300-for-grid.list', fit_hw=0.6, role='diagnostic'),
}

# EXCLUDED up front, with the reason recorded in the output (never silently dropped):
EXCLUDED = {
    5369.589: ('the GES NLTE row nearest 5369.589 is a DIFFERENT transition (EP 4.259, '
               'log gf -1.174 = the K08 5369.585 line) and its upper level is UNIDENTIFIED '
               "('0') — no HFS group and no computable departure. Co I 5369.589 (EP 1.740, "
               '4 HFS components) therefore has no NLTE route in this deck.'),
    5454.570: ('upper model-atom level UNIDENTIFIED in the GES NLTE list (index 0) — bsyn '
               'would silently fall back to departure=1 (the RYA-533 trap). LTE-only; '
               'excluded from an NLTE-corrected value rather than reported as if corrected.'),
}

HDR_RE = re.compile(r"^'\s*([\d.]+)\s*'\s+(\d+)\s+(\d+)\s*$")


# ────────────────────────────── single source of truth ──────────────────────────────
def verify_canonical_gf(root):
    """Assert every LINES log gf equals the committed canonical_gf.csv value. Loud-fail on
    drift or absence — never synthesise on an unsourced gf, never silently fall back."""
    path = os.path.join(root, "data", "linelists", "canonical_gf.csv")
    if not os.path.exists(path):
        raise SystemExit(f"RYA-564 SSOT: canonical_gf.csv not found at {path}")
    have = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            if row["species"] == "Co I":
                have[round(float(row["wavelength_air_A"]), 3)] = row
    out = {}
    for center, cfg in LINES.items():
        row = have.get(round(center, 3))
        if row is None:
            for w, r in have.items():
                if abs(w - center) < 0.02:
                    row = r
                    break
        if row is None:
            raise SystemExit(f"RYA-564 SSOT: no Co I {center} in canonical_gf.csv — fix the "
                             f"single source, do not hardcode.")
        cg = float(row["log_gf"])
        if abs(cg - cfg["loggf"]) > 0.005:
            raise SystemExit(f"RYA-564 SSOT DRIFT: Co I {center} canonical log_gf={cg} != script "
                             f"{cfg['loggf']} ({row['line_id']}). Reconcile to the single source.")
        cep = float(row["excitation_potential_eV"])
        if abs(cep - cfg["ep"]) > 0.01:
            raise SystemExit(f"RYA-564 SSOT DRIFT: Co I {center} canonical EP={cep} != script "
                             f"{cfg['ep']}.")
        out[center] = dict(line_id=row["line_id"], hfs_n=row["hfs_n_components"],
                           adjudication=row["adjudication_status"], ref=row["loggf_reference"])
        print(f"  gf OK  Co I {center:9.3f}: loggf={cg:+.3f}  EP={cep:.3f}  "
              f"({row['line_id']}, {row['loggf_reference']}, {row['adjudication_status']}, "
              f"hfs_n={row['hfs_n_components']})")
    return out


# ────────────────────────────── line lists ──────────────────────────────
def _ges_co_rows():
    """Parse the Co I block of the GES NLTE line list ONCE. Returns (header, [row-dicts])
    with the level indices parsed verbatim (they are what makes bsyn engage departures)."""
    with open(GES) as f:
        lines = f.readlines()
    blk = None
    for i, ln in enumerate(lines):
        m = HDR_RE.match(ln.rstrip("\n"))
        if m and int(float(m.group(1))) == Z_CO and int(m.group(2)) == 1:
            blk = i
            break
    if blk is None:
        raise SystemExit("RYA-564: no GES NLTE block for Co I (Z=27 stage=1)")
    hdr = lines[blk].rstrip("\n")
    rows = []
    for ln in lines[blk + 2:]:
        if HDR_RE.match(ln.rstrip("\n")):
            break
        p = ln.split()
        try:
            wl, ep, gf = float(p[0]), float(p[1]), float(p[2])
        except (ValueError, IndexError):
            continue
        m = re.search(r"'\s+(\d+)\s+(\d+)\s+'", ln)
        rows.append(dict(wl=wl, ep=ep, gf=gf,
                         lv=(int(m.group(1)), int(m.group(2))) if m else None,
                         raw=ln.rstrip("\n")))
    return hdr, rows


def co_hfs_group(center, cfg, ges_rows):
    """The HFS component group for `center`: the GES rows sharing the nearest row's
    model-atom level pair within 0.3 A. Asserts both levels are identified (else bsyn
    silently drops to LTE) and reconciles the component gf sum onto the canonical total.
    Returns (rows, info-dict). Rows are returned VERBATIM unless a uniform rescale onto
    the canonical total was required, which is flagged."""
    near = min(ges_rows, key=lambda r: abs(r['wl'] - center))
    if abs(near['wl'] - center) > 0.08:
        raise SystemExit(f"RYA-564: Co I {center} absent from the GES NLTE block "
                         f"(nearest {near['wl']:.3f})")
    lv = near['lv']
    if not lv or lv[0] == 0 or lv[1] == 0:
        raise SystemExit(f"RYA-564: Co I {center} has an UNIDENTIFIED model-atom level "
                         f"{lv} — bsyn would run silent-LTE (RYA-533 trap). RAISE.")
    grp = [r for r in ges_rows if r['lv'] == lv and abs(r['wl'] - center) < 0.30]
    ges_sum = math.log10(sum(10 ** r['gf'] for r in grp))
    d = ges_sum - cfg['loggf']
    rows, scaled = [r['raw'] for r in grp], False
    if abs(d) > 0.005:                       # enforce the single source, preserve the pattern
        scaled = True
        rows = []
        for r in grp:
            new_gf = r['gf'] - d
            # replace the 3rd whitespace field (log gf) in place, keeping the row verbatim
            m = re.match(r"^(\s*\S+\s+\S+\s+)(\S+)(.*)$", r['raw'])
            if not m:
                raise SystemExit(f"RYA-564: cannot rewrite gf on GES row: {r['raw'][:80]}")
            rows.append(f"{m.group(1)}{new_gf:+.3f}{m.group(3)}")
    info = dict(n_components=len(grp), levels=list(lv),
                ges_sum_loggf=round(ges_sum, 4), canonical_loggf=cfg['loggf'],
                ges_minus_canonical=round(d, 4), rescaled_to_canonical=scaled,
                span_mA=round((max(r['wl'] for r in grp) - min(r['wl'] for r in grp)) * 1000, 1))
    return rows, info


def write_co_list(center, rows, path):
    with open(path, "w") as f:
        f.write(f"'  {Z_CO}.000            '    1    {len(rows)}\n")
        f.write("'Co I'\n")
        for r in rows:
            f.write(r + "\n")
    return path


def write_blend_list(center, vald_name, path, lmin, lmax):
    """The VALD in-window block for the BLENDS, with Co I REMOVED (Co is supplied by the
    HFS-resolved GES list — keeping both would double-count the element being measured).
    Per-species line counts are rewritten to match what is actually emitted."""
    src = os.path.join(VALD_DIR, vald_name)
    blocks, cur = [], None
    with open(src) as f:
        it = iter(f)
        for ln in it:
            m = HDR_RE.match(ln.rstrip("\n"))
            if m:
                name = next(it).rstrip("\n")
                cur = dict(code=m.group(1), stage=m.group(2), name=name, rows=[])
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
    n_co_dropped = 0
    with open(path, "w") as f:
        for b in blocks:
            if not b['rows']:
                continue
            if b['name'].strip().strip("'").strip().startswith('Co I'):
                n_co_dropped += len(b['rows'])
                continue
            f.write(f"'  {b['code']}            '    {b['stage']}    {len(b['rows'])}\n")
            f.write(b['name'] + "\n")
            for r in b['rows']:
                f.write(r + "\n")
    return path, n_co_dropped


# ────────────────────────────── engine ──────────────────────────────
def make_departure():
    """Interpolate the Gerber Co departure grid onto the solar MARCS node (once).
    This is the RYA-534-validated staged grid — READ, never fitted."""
    dep = f"{W}/Testout/sun_Co_coef.dat"
    if os.path.exists(dep):
        return dep
    os.makedirs(f"{W}/Testout", exist_ok=True)
    binf, aux = f"{GT}/{GRID}", f"{GT}/{AUX}"
    for p in (binf, aux):
        if not os.path.exists(p):
            raise SystemExit(f"RYA-564: Gerber Co grid artifact missing on Sirius: {p}")
    with open(aux) as f:
        al = f.readlines()
    nhead = 0
    for ln in al:
        if ln.lstrip().startswith('#'):
            nhead += 1
        else:
            break
    nrows = len(al) - nhead
    stdin = "\n".join([f"'{MARCS}'"] * 8 + [
        f"'{W}/Testout/sun_Co.interpol'", f"'{W}/Testout/sun_Co.alt'", f"'{dep}'",
        f"'{binf}'", f"'{aux}'", str(nrows),
        TREF, LOGGREF, ZREF, f"{A_SUN:.2f}", ".false.", ".false.", "'none'", ""])
    r = subprocess.run([INTERP], input=stdin, capture_output=True, text=True, cwd=W)
    if not os.path.exists(dep):
        raise SystemExit(f"RYA-564: interpol_modeles_nlte failed:\n{r.stdout[-1500:]}\n{r.stderr[-500:]}")
    return dep


def write_nlteinfo(path):
    with open(path, "w") as f:
        f.write(f"# RYA-564 Co NLTE (Gerber TS-native deck, RYA-533/534)\n"
                f"# path for model atom files\n{GT}/\n"
                f"# path for departure files\n{W}/Testout/\n# atomic NLTE setup\n"
                f"{Z_CO}  'Co'  'nlte'  '{ATOM}'  'sun_Co_coef.dat'  'ascii'\n")
    return path


def babsma(opac, lmin, lmax, model=None, marcs_file=True):
    ctl = (f"'LAMBDA_MIN:'  '{lmin}'\n'LAMBDA_MAX:'  '{lmax}'\n'LAMBDA_STEP:' '0.005'\n"
           f"'MODELINPUT:' '{model or MARCS}'\n"
           f"'MARCS-FILE:' '{'.true.' if marcs_file else '.false.'}'\n"
           f"'MODELOPAC:' '{opac}'\n"
           f"'METALLICITY:'    '0.00'\n'ALPHA/Fe   :'    '0.00'\n'HELIUM     :'    '0.00'\n"
           f"'R-PROCESS  :'    '0.00'\n'S-PROCESS  :'    '0.00'\n'XIFIX:' 'T'\n{XI}\n")
    r = subprocess.run([f"{EXE}/babsma_lu"], input=ctl, capture_output=True, text=True, cwd=W)
    if not os.path.exists(opac):
        raise SystemExit(f"RYA-564: babsma failed ({lmin}-{lmax}):\n{r.stdout[-1500:]}\n{r.stderr[-500:]}")


def bsyn(a_x, linelists, opac, lmin, lmax, tag, nlte=False, nlteinfo=None):
    """Run bsyn over `linelists` (list of paths). Returns (wave, flux, engaged).
    For NLTE, `engaged` is False unless bsyn reports it actually read departures for the
    species AND did not fall back to departure=1 — the caller RAISES on a silent LTE."""
    res = f"{W}/s_{tag}.spec"
    if os.path.exists(res):
        os.remove(res)
    ctl = (f"'NLTE :'          '{'.true.' if nlte else '.false.'}'\n"
           + (f"'NLTEINFOFILE:'  '{nlteinfo}'\n" if nlte else "")
           + f"'LAMBDA_MIN:'     '{lmin}'\n'LAMBDA_MAX:'     '{lmax}'\n"
           f"'LAMBDA_STEP:'    '0.005'\n'INTENSITY/FLUX:' 'Flux'\n'ABFIND        :' '.false.'\n"
           f"'MODELOPAC:' '{opac}'\n'RESULTFILE :' '{res}'\n'METALLICITY:'    '0.00'\n"
           f"'INDIVIDUAL ABUNDANCES:'   '1'\n{Z_CO}  {a_x:.3f}\n'ISOTOPES : ' '0'\n"
           f"'NFILES   :' '{len(linelists)}'\n" + "".join(f"{p}\n" for p in linelists)
           + f"'SPHERICAL:'  'F'\n  30\n  300.00\n  15\n  1.30\n")
    r = subprocess.run([f"{EXE}/bsyn_lu"], input=ctl, capture_output=True, text=True, cwd=W)
    if not os.path.exists(res):
        raise SystemExit(f"RYA-564: bsyn failed {tag} A={a_x}:\n{r.stdout[-2000:]}")
    engaged = True
    if nlte:
        so = r.stdout
        engaged = ("read departure coefficients" in so
                   and "missing level identification" not in so
                   and "continuing with departure coefficient = 1. for that level" not in so)
    d = np.loadtxt(res)
    return d[:, 0], d[:, 1], engaged


# ────────────────────────────── observed arms + fitting ──────────────────────────────
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
    vac->air through the pipeline SINGLE-SOURCE converter (RYA-501)."""
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


def fit_line(center, cfg, obs_w, obs_f, synth):
    """Thin adapter onto the SHARED fitter (scripts/solar_profile_fit.fit_profile).

    RYA-643 lineage audit: this harness was a FOURTH copy of the RYA-551 fitter and
    carried both defects RYA-592 fixed in its own copy — no rest-frame handling, and a
    broadening grid the fit railed against at its 1.5 km/s floor. Now single-sourced."""
    return fit_profile(center, obs_w, obs_f, synth,
                       hw=cfg['fit_hw'], vsini=VSINI, a_lo=A_LO, a_hi=A_HI)


def nlte_delta(center, co_list, nlteinfo, model, opac_nlte, lmin, lmax):
    """PER-LINE 1D-NLTE correction, by the RYA-534-validated method: EW of the NLTE
    synthesis at the solar A, inverted through an LTE curve-of-growth on the SAME
    (Co-only) line list. delta = A_sun - A*. READ from the staged grid, never fitted.
    RAISES if departures do not engage — a line is never silently reported as LTE."""
    def _ew(w, f):
        m = (w > center - 1.2) & (w < center + 1.2)
        return float(_trapezoid(1.0 - f[m], w[m]) * 1000.0)

    w, f, engaged = bsyn(A_SUN, [co_list], opac_nlte, lmin, lmax,
                         f"{center:.0f}_nlte", nlte=True, nlteinfo=nlteinfo)
    if not engaged:
        raise SystemExit(f"RYA-564: Co I {center} — NLTE departures NOT engaged (silent LTE). "
                         f"STOP (RYA-533 trap).")
    ew_n = _ew(w, f)
    As, ews = [], []
    for off in COG_OFFSETS:
        a = A_SUN + off
        w2, f2, _ = bsyn(a, [co_list], opac_nlte, lmin, lmax, f"{center:.0f}_lte_{a:.2f}")
        As.append(a)
        ews.append(_ew(w2, f2))
    order = np.argsort(ews)
    a_star = float(np.interp(ew_n, np.array(ews)[order], np.array(As)[order]))
    return dict(delta=round(A_SUN - a_star, 4), ew_nlte_mA=round(ew_n, 3),
                a_star=round(a_star, 4), engaged=True,
                cog_EW_mA=[round(e, 3) for e in ews], cog_A=[round(a, 3) for a in As])


# ────────────────────────────── main ──────────────────────────────
A_SUN = None    # set from SOLAR_ASPLUND2021 in main() — never hardcoded


def main():
    global A_SUN
    ap = argparse.ArgumentParser(description="RYA-564 Co I red-line HFS synthesis (Sirius)")
    ap.add_argument('--only', nargs='*', type=float, help="restrict to these line centres")
    ap.add_argument('--out', default=None, help="output JSON path (default data/results/...)")
    args = ap.parse_args()

    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, root)
    # RYA-567: Sirius-only heavy compute (TS engine + Gerber NLTE grid + MARCS). Loud-fail
    # off-Sirius — never against local-Mac copies of engines/grids/spectra.
    from config.constants import assert_on_sirius, SOLAR_ASPLUND2021
    assert_on_sirius("RYA-564 Co I red-line synthesis", require_subdirs=("engines", "grids"))
    A_SUN = float(SOLAR_ASPLUND2021['Co'])

    print(f"RYA-564 Co I red-line HFS synthesis — reference A(Co)_sun (Asplund 2021) = {A_SUN} "
          f"[grid bracket + delta node only, NOT a tuning target]\n")
    print("Step-0 SSOT check (log gf from data/linelists/canonical_gf.csv):")
    canon = verify_canonical_gf(root)

    os.makedirs(W, exist_ok=True)
    if not os.path.lexists(f"{W}/DATA"):
        os.symlink("/mnt/codex-data/engines/Turbospectrum_NLTE/DATA", f"{W}/DATA")

    ges_hdr, ges_rows = _ges_co_rows()
    print(f"\nGES NLTE Co I block: {len(ges_rows)} rows ({ges_hdr.strip()})")

    dep = make_departure()
    nlteinfo = write_nlteinfo(f"{W}/species_Co.dat")
    model_nlte = f"{W}/Testout/sun_Co.interpol"     # the interpolator's MARCS+departure model
    print(f"Gerber Co departures interpolated onto the solar node -> {os.path.basename(dep)}")

    arms = {'harps': load_harps()}
    try:
        arms['iag'] = load_iag(root)
        print(f"IAG atlas: {len(arms['iag'][0])} pts, "
              f"{arms['iag'][0].min():.1f}-{arms['iag'][0].max():.1f} A")
    except Exception as e:
        print(f"IAG load FAILED ({e}); HARPS-only")

    # RYA-643: the rest-frame correction must be SOURCED from a measurement on
    # each arm being fitted — never hardcoded, never a silent zero. Loud-fails if
    # the clean check lines cannot deliver one.
    arm_rv = {}
    print('  residual velocity per arm (clean-line centroids, abundance-blind):')
    for _arm, (_ow, _of) in arms.items():
        _v, _n, _sd = require_arm_rv(_ow, _of, _arm)
        arm_rv[_arm] = dict(v_kms=round(_v, 3), n_lines=_n, scatter_kms=round(_sd, 3))
        print(f'    {_arm:6s}: {_v:+.3f} km/s (n={_n}, scatter {_sd:.3f})')

    centers = sorted(LINES) if not args.only else [c for c in sorted(LINES)
                                                   if any(abs(c - o) < 0.05 for o in args.only)]
    results = {'_meta': dict(
        ticket='RYA-564', element='Co', ion='I', z=Z_CO,
        a_sun_ref=A_SUN, a_sun_ref_source='Asplund2021 (bracket/delta node, not a target)',
        method=('HFS-resolved Turbospectrum LTE flux fit (GES NLTE line list supplies the Co I '
                'HFS components; VALD in-window block supplies the blends, Co I removed) + '
                'PER-LINE 1D-NLTE delta read from the RYA-534-validated Gerber TS-native grid'),
        nlte_grid=GRID, nlte_atom=ATOM, nlte_engine='Engine-B TS-Gerber (RYA-533/534 deck)',
        nlte_anchor=ANCHOR_DELTA, nlte_anchor_ref='Bergemann+2010 MNRAS 401 1334 (RYA-534 +0.099)',
        reliable_dEW_dA_floor=RELIABLE_DEWDA, reliable_red_chi2_ceiling=RELIABLE_RCHI2,
        arms=sorted(arms), vsini=VSINI, xi=XI, marcs_node=f"{TREF}/{LOGGREF}/{ZREF}",
        excluded_lines={str(k): v for k, v in EXCLUDED.items()},
        demotes=('Co I 3845.468 (Kitt Peak blue edge, SNR~24, chi2r~3100, A=6.128 / +1.188) '
                 '-> DIAGNOSTIC-ONLY; never the reported value'))}

    for center in centers:
        cfg = LINES[center]
        print(f"\n=== Co I {center:.3f} ({cfg['role']}, EP {cfg['ep']}, canonical loggf "
              f"{cfg['loggf']:+.3f}, {cfg['gf_ref']}, gf {cfg['adj']}) ===")
        rows, hfs = co_hfs_group(center, cfg, ges_rows)
        print(f"  HFS: {hfs['n_components']} components, span {hfs['span_mA']} mA, levels "
              f"{hfs['levels']}, GES sum {hfs['ges_sum_loggf']:+.3f} vs canonical "
              f"{cfg['loggf']:+.3f} (d={hfs['ges_minus_canonical']:+.3f}"
              f"{', RESCALED to canonical' if hfs['rescaled_to_canonical'] else ''})")
        co_list = write_co_list(center, rows, f"{W}/co_{center:.0f}.list")
        lmin, lmax = center - LMARGIN, center + LMARGIN
        blend_list, n_drop = write_blend_list(center, cfg['vald'],
                                              f"{W}/blend_{center:.0f}.list", lmin, lmax)
        print(f"  blends: VALD window {lmin:.1f}-{lmax:.1f} A, {n_drop} Co I rows removed "
              f"(Co comes from the HFS list)")

        # (2) per-line NLTE delta — Co-only list, on the departure-bearing model
        opac_n = f"{W}/opac_n_{center:.0f}"
        babsma(opac_n, lmin, lmax, model=model_nlte, marcs_file=False)
        dl = nlte_delta(center, co_list, nlteinfo, model_nlte, opac_n, lmin, lmax)
        print(f"  NLTE: departures ENGAGED; EW_NLTE={dl['ew_nlte_mA']:.3f} mA, A*={dl['a_star']:.4f}"
              f"  ->  delta = {dl['delta']:+.4f}  (anchor {ANCHOR_DELTA:+.3f} +/- {ANCHOR_TOL})")

        # (1) LTE flux fit with blends, over the A grid
        opac_l = f"{W}/opac_l_{center:.0f}"
        babsma(opac_l, lmin, lmax)
        synth = {}
        for a in A_GRID:
            w, f, _ = bsyn(float(a), [blend_list, co_list], opac_l, lmin, lmax,
                           f"{center:.0f}_fit_{a:.2f}")
            synth[float(a)] = (w, f)

        rec = dict(wave=center, ep=cfg['ep'], loggf=cfg['loggf'], gf_ref=cfg['gf_ref'],
                   gf_adjudication=cfg['adj'], role=cfg['role'], canonical=canon[center],
                   hfs=hfs, nlte=dl, blends_dropped_co_rows=n_drop)
        for arm, (ow, of) in arms.items():
            fit = fit_line(center, cfg, ow, of, synth)
            if fit is None:
                rec[arm] = dict(status='no_coverage')
                print(f"  {arm:6s}: no coverage")
                continue
            a_lte = fit['A']
            a_nlte = a_lte + dl['delta']
            reliable = bool((not fit['railed'])
                            and fit['dEW_dA'] >= RELIABLE_DEWDA
                            and fit['red_chi2'] <= RELIABLE_RCHI2)
            rec[arm] = dict(dv_fitted_kms=round(fit['dv'], 2),
                            dv_measured_kms=arm_rv[arm]['v_kms'],
                            gsig_railed=fit['gsig_railed'],
                            A_LTE=round(a_lte, 3), A_NLTE=round(a_nlte, 3),
                            nlte_delta=dl['delta'], gsig_kms=round(fit['gsig'], 2),
                            red_chi2=round(fit['red_chi2'], 2), npix=fit['npix'],
                            core_EW_mA=fit['core_EW_mA'], dEW_dA_mA_dex=fit['dEW_dA'],
                            railed=fit['railed'], reliable=reliable)
            print(f"  {arm:6s}: A_LTE={a_lte:.3f} {dl['delta']:+.4f} -> A_NLTE={a_nlte:.3f}  "
                  f"(EW~{fit['core_EW_mA']:.1f} mA, gsig={fit['gsig']:.1f} km/s, "
                  f"rchi2={fit['red_chi2']:.1f}, dEW/dA={fit['dEW_dA']} mA/dex, "
                  f"railed={fit['railed']}, reliable={reliable})")
        results[f"{center:.3f}"] = rec

    # ── summary: the reportable value comes from RELIABLE primary/secondary lines only ──
    def _reliable(arm, roles):
        return [(c, results[f"{c:.3f}"][arm]['A_NLTE'])
                for c in centers
                if LINES[c]['role'] in roles
                and isinstance(results[f"{c:.3f}"].get(arm), dict)
                and results[f"{c:.3f}"][arm].get('reliable')]

    def _stats(pairs):
        if not pairs:
            return None
        v = np.array([x for _, x in pairs])
        return dict(median=round(float(np.median(v)), 3), mean=round(float(v.mean()), 3),
                    scatter=round(float(v.std(ddof=1)) if len(v) > 1 else 0.0, 3),
                    n=len(v), lines={f"{c:.3f}": x for c, x in pairs})

    rel = _reliable('harps', ('primary', 'secondary'))
    deltas = [results[f"{c:.3f}"]['nlte']['delta'] for c in centers]
    med_delta = float(np.median(deltas)) if deltas else float('nan')
    anchor_ok = bool(abs(med_delta - ANCHOR_DELTA) <= ANCHOR_TOL)
    print("\n" + "=" * 72)
    print(f"NLTE cross-check: median per-line delta {med_delta:+.4f} vs the RYA-534/Bergemann "
          f"anchor {ANCHOR_DELTA:+.3f} +/- {ANCHOR_TOL}  ->  {'CONSISTENT' if anchor_ok else 'CHECK'}")
    summary = dict(median_nlte_delta=round(med_delta, 4), nlte_anchor_consistent=anchor_ok,
                   n_reliable=len(rel),
                   # gf-quality sensitivity: 'primary' = gf adjudication 'agreed' AND the GES
                   # HFS sum already equals the canonical total (no rescale). Reporting both
                   # sets makes the value's dependence on the 'pending'-gf lines visible.
                   primary_only=_stats(_reliable('harps', ('primary',))),
                   primary_plus_secondary=_stats(rel),
                   iag_crosscheck=_stats(_reliable('iag', ('primary', 'secondary'))))
    if rel:
        vals = np.array([v for _, v in rel])
        a_co = float(np.median(vals))
        summary.update(A_Co=round(a_co, 3), A_Co_mean=round(float(vals.mean()), 3),
                       scatter=round(float(vals.std(ddof=1)) if len(vals) > 1 else 0.0, 3),
                       reliable_lines={f"{c:.3f}": v for c, v in rel},
                       delta_vs_asplund=round(a_co - A_SUN, 3))
        print(f"RELIABLE Co I lines (HARPS): {', '.join(f'{c:.3f}={v:.3f}' for c, v in rel)}")
        print(f"A(Co) = {a_co:.3f}  (median of n={len(vals)}, scatter "
              f"{summary['scatter']:.3f}; {a_co - A_SUN:+.3f} vs Asplund {A_SUN})")
        po, ia = summary['primary_only'], summary['iag_crosscheck']
        if po:
            print(f"  gf-quality check — gf-'agreed' lines only: median {po['median']:.3f} "
                  f"(n={po['n']}, scatter {po['scatter']:.3f})")
        if ia:
            print(f"  arm cross-check  — IAG FTS atlas:          median {ia['median']:.3f} "
                  f"(n={ia['n']}, scatter {ia['scatter']:.3f})")
    else:
        summary.update(A_Co=None, reason=('no red Co I line cleared the reliability floor '
                                          f'(dEW/dA >= {RELIABLE_DEWDA} mA/dex, red_chi2 <= '
                                          f'{RELIABLE_RCHI2}, not railed)'))
        print(f"NO red Co I line cleared the reliability floor -> Co reports NO VALUE (owed).")
        print("Do NOT fall back to the blue-edge 3845 artifact (ticket CRITICAL).")
    results['_summary'] = summary

    out_path = args.out or os.path.join(root, "data", "results", "co_synthesis_rya564.json")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nWrote {out_path}")


if __name__ == '__main__':
    main()
