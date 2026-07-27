#!/usr/bin/env python3
"""
RYA-592 — Magnesium (Mg I) 5528 SECOND-LINE measurement (RUNS ON SIRIUS).

WHY: Mg is frozen CURATION-OWED-with-value 7.614 (+0.064) on the SINGLE line
5711.088. RYA-561's ratified gate-3 (strict) requires a REAL independent
confirmation, and Mg has none: Engine-A is suppressed by design and 5711 is the
only clean line in the measured pool, so Mg is single-line AND single-engine.
This delivers the second clean line — Mg I 5528.405 — by in-window blend-fit
synthesis, so any Mg promotion rests on real evidence, never a rule relaxation.

METHOD (the RYA-551 Sr II / RYA-560 Zr II in-window pattern, NOT an EW inversion):
synthesise the FULL window with every VALD3 blend component present (344 components
within +-1.6 A of 5528 — TiO/MgH/CC/CN plus Zr I 5528.410 sitting 0.005 A off the Mg
core), renormalise a local continuum, and fit A(Mg) to the OBSERVED solar profile by
chi2 over the fit window. 5711 is re-measured with the SAME harness so the
line-to-line concordance compares like with like (the committed 7.614 comes from a
different channel, synth-v2 + Gerber delta, and is reported alongside as context).

BOTH NLTE ENGINES, deltas READ from the committed grids (never fitted):
  Engine-A = Amarsi2020 PySME (nlte_Mg_scatt_pysme.grd) via pipeline.pysme_nlte —
             the machinery that BUILT data/nlte_grids/Mg_Amarsi2020_PySME.csv. That
             CSV only tabulates 5711/4730, so 5528 is DERIVED here from the same
             committed .grd by the same code path.
  Engine-B = TS-native Gerber 2023 (atom.mg86d + NLTEgrid4TS_Mg) via
             scripts/ts_gerber_gate — the RYA-534 deck.
CONTROL (validate-don't-tune): each engine must reproduce its OWN committed 5711
delta before its new 5528 number is trusted. Engine-A control = the committed CSV
solar node (-0.0343); Engine-B control = the RYA-534 gate record (-0.023). A control
miss is a hard STOP — no silent acceptance of an unreproduced grid.

GATE-3 DISCIPLINE (RYA-561, Ryan 2026-07-27, STRICT): gate 3 needs a REAL
cross-engine dCE, and the cross-engine delta of two NLTE ATOMS on the same profile
fit is NOT that — it is gate 1 under a second name, explicitly rejected. So the
per-engine A_NLTE below is reported as required, but the Engine-A/Engine-B ATOM
agreement is recorded as diagnostic ONLY and must never be fed to gate 3. Whether a
real Engine-A (EW->A) value is even admissible on 5528 is decided by the ratified
TWO_ENGINE clause-3 classifier (saturation knee), read from config.constants — not
by this script's opinion.

Stays 1D-NLTE. The 1D-NLTE-vs-3D-NLTE scale caveat (RYA-561) rides on the value;
the class-wide 3D-metals term is a separate post-Beta ticket. No 3D here.

Usage (Sirius):  /mnt/codex-data/venv_pysme/bin/python scripts/rya592_mg_5528_synth_sirius.py
"""
import csv
import json
import os
import re
import subprocess
import sys

import numpy as np

EXE = "/mnt/codex-data/engines/Turbospectrum_NLTE/exec-gf"
MARCS = ("/mnt/codex-data/grids/model_atmospheres/marcs_standard_comp/marcs_standard_comp/"
         "p5750_g+4.5_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod")
VALD_DIR = "/mnt/codex-data/engines/TSFitPy/input_files/linelists/linelist_vald"
GES_TSV = "/mnt/codex-data/linelists/ges_v6/atomic_lines.tsv"
W = "/mnt/codex-data/codex/rya592/work"
W_GERBER = "/mnt/codex-data/codex/rya592/gerber"
HARPS = "/mnt/codex-data/codex/_solar_tmp/solar_normalized_harps.csv"
IAG = "/mnt/codex-data/solar_reference/iag_reiners2016/spvis.dat.gz"

Z_MG = 12             # Mg atomic number (bsyn INDIVIDUAL ABUNDANCES element)
XI = 1.0              # solar microturbulence (matches the t01 MARCS node)
VSINI = 1.8           # solar vsini (STAR_PARAMS)
CLIGHT = 299792.458

# Canonical Mg I lines: (wave_air, EP_eV, log_gf) — the CANONICAL single-source gf
# from data/linelists/canonical_gf.csv; verify_canonical_gf() asserts these match the
# file at runtime (SSOT, loud-fail). fit_hw is set from the OBSERVED profile: 5528 is
# deep (core flux 0.21) with wings out to ~1.4 A; 5711 must stay narrow because a
# strong neighbour sits ~3 A blueward (flux 0.74 at 5708).
LINES = {
    5528.405: dict(ep=4.346, loggf=-0.498, gf_id='gf_087251',
                   gf_ref='NIST ASD v5.12 grade B+ (RYA-592 Step 0)',
                   vald='vald-5300-5800-for-grid.list', fit_hw=1.4, role='target'),
    5711.088: dict(ep=4.346, loggf=-1.724, gf_id='gf_087252',
                   gf_ref='NIST ASD v5.12 grade B (RYA-592 Step 0)',
                   vald='vald-5300-5800-for-grid.list', fit_hw=0.6,
                   role='concordance_reference'),
}
TARGET = 5528.405
REFLINE = 5711.088

LMARGIN = 8.0                                       # synth half-window (A) for babsma/bsyn
A_GRID = np.round(np.arange(7.00, 8.16, 0.05), 3)   # A(Mg) trial grid (brackets Asplund 7.55)
A_LO, A_HI = float(A_GRID.min()), float(A_GRID.max())

RELIABLE_DEWDA = 40.0     # mA/dex core-EW sensitivity floor (as RYA-551/560)

# Line-to-line concordance band. Same 0.10 dex "agree within ~0.1 dex" convention the
# ratified two-engine floor uses for cross-engine agreement (TWO_ENGINE
# ['cross_engine_mix_gate']) — reused here for LINE-to-LINE agreement so the number is
# not a fresh invention. Declared BEFORE the run (reference-blind), never tuned to the
# result.
CONCORDANCE_BAND = 0.10

# Rest-frame fit (RYA-309): the flux fit compares REST-frame synthetic flux to the
# observed profile, so a residual velocity offset misaligns the cores, inflates chi2
# and biases the broadening. The solar HARPS arm carries an uncorrected +0.8 km/s
# (measured independently below from 9 clean Fe/Si line centroids), so dv is fitted as
# a NUISANCE parameter and cross-checked against that measurement. dv is never tuned
# toward an abundance — it is checked against the centroids, which know nothing about A.
DV_GRID = np.round(np.arange(-1.5, 2.55, 0.10), 2)      # km/s
GSIG_GRID = np.round(np.arange(0.4, 8.01, 0.20), 2)     # km/s (vmac + instrumental)

# Clean, unblended solar lines used ONLY to measure the arm's residual velocity, as an
# independent check on the fitted dv. Air wavelengths; no abundance information is taken
# from them.
RV_CHECK_LINES = [5522.446, 5525.544, 5615.644, 5679.023, 5686.530,
                  5701.104, 5701.544, 5709.378, 5731.762]

# Committed per-engine 5711 solar deltas — the validate-don't-tune CONTROLS.
CTRL_TOL = 0.005          # dex; a re-derivation must land this close to its committed value
ENGINE_B_COMMITTED_5711 = -0.023   # RYA-534 Mg_gate_result.txt / Mg_gerber2023.prov.json


# ── Step 0: gf provenance (SSOT gate) ────────────────────────────────────────

def verify_canonical_gf(root):
    """SSOT invariant: assert every LINES gf equals the committed canonical_gf.csv
    value (Step-0 provenance gate). Loud-fail on drift or missing line — never run
    the synthesis on an unsourced/fabricated gf, never silently fall back."""
    path = os.path.join(root, "data", "linelists", "canonical_gf.csv")
    if not os.path.exists(path):
        raise SystemExit(f"RYA-592 SSOT: canonical_gf.csv not found at {path}")
    have = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            if row["species"] == "Mg I":
                have[round(float(row["wavelength_air_A"]), 3)] = row
    for center, cfg in LINES.items():
        row = have.get(round(center, 3))
        if row is None:
            for w, r in have.items():
                if abs(w - center) < 0.01:
                    row = r
                    break
        if row is None:
            raise SystemExit(f"RYA-592 SSOT: no Mg I {center} in canonical_gf.csv — cannot "
                             f"source gf (Step-0 gate). Fix the single source, do not hardcode.")
        cg = float(row["log_gf"])
        if abs(cg - cfg["loggf"]) > 0.005:
            raise SystemExit(f"RYA-592 SSOT DRIFT: Mg I {center} canonical log_gf={cg} != script "
                             f"{cfg['loggf']} ({row['line_id']}). Reconcile to the single source.")
        if row["line_id"] != cfg["gf_id"]:
            raise SystemExit(f"RYA-592 SSOT: Mg I {center} line_id {row['line_id']} != expected "
                             f"{cfg['gf_id']} — the canonical row moved; reconcile before running.")
        print(f"  gf OK  Mg I {center}: loggf={cg:+.3f}  grade={row['nist_grade']}  "
              f"({row['line_id']}, {row['loggf_reference'][:48]}...)")


# ── Engine-A (Amarsi2020 PySME) NLTE delta, derived from the committed .grd ──

def _ges_atomic(el, wl, ion='1', dw=0.25, ep_tol=0.03):
    """(loggf_total, EP, Eup, gamvw) for the feature at `wl` from the GES v6 synthesis
    line list — a pandas port of resource_clamped_grids_rya409._resolve_atomic, which
    reads the same file through iSpec. iSpec is NOT installed on Sirius and RYA-567
    pins this leg to Sirius, so the ispec path is unavailable here; the clustering
    logic (group components by lower-level energy, keep the dominant group, gf-sum
    within it) is reproduced verbatim so the derivation matches the committed grid."""
    import pandas as pd
    if not os.path.exists(GES_TSV):
        raise SystemExit(f"RYA-592: GES v6 synthesis line list not staged at {GES_TSV} — "
                         f"Engine-A atomic data has no single source. Stage it, do not guess.")
    ll = pd.read_csv(GES_TSV, sep='\t', low_memory=False)
    m = ll[(ll['element'] == f'{el} {ion}') & ((ll['wave_A'] - wl).abs() < dw)]
    if len(m) == 0:
        raise SystemExit(f"{el} {wl}: not in the GES linelist (+-{dw} A).")
    eps = m['lower_state_eV'].to_numpy(float)
    gfs = 10 ** m['loggf'].to_numpy(float)
    vws = m['waals'].to_numpy(float)
    order = np.argsort(eps)
    groups, cur = [], [order[0]]
    for k in order[1:]:
        if eps[k] - eps[cur[-1]] <= ep_tol:
            cur.append(k)
        else:
            groups.append(cur)
            cur = [k]
    groups.append(cur)
    g = max(groups, key=lambda idx: gfs[idx].sum())
    return (float(np.log10(gfs[g].sum())), float(np.median(eps[g])),
            float(np.median(eps[g])) + 12398.42 / wl,
            float(np.median(vws[g])) if float(np.median(vws[g])) > 0 else 0.0)


def engine_a_deltas(root):
    """Per-line Engine-A delta from the COMMITTED Amarsi-2020 .grd via PySME.

    The committed Mg_Amarsi2020_PySME.csv tabulates only 5711/4730, so 5528 does not
    exist there yet; it is derived here through pipeline.pysme_nlte — the same module,
    same solar node, same derivation options that produced the committed rows. 5711 is
    re-derived alongside as the CONTROL and must reproduce the committed value."""
    import pandas as pd
    from pipeline import pysme_nlte as P
    from config.constants import SOLAR_ASPLUND2021

    grid_csv = os.path.join(root, 'data', 'nlte_grids', 'Mg_Amarsi2020_PySME.csv')
    committed = pd.read_csv(grid_csv)
    ctl_rows = committed[(committed.wave_A == REFLINE) & (committed.teff_K == 5772)
                         & (committed.logg == 4.44) & (committed.feh == 0.0)]
    if ctl_rows.empty:
        raise SystemExit(f"RYA-592: no committed solar-node row for Mg {REFLINE} in "
                         f"{grid_csv} — no control to validate the re-derivation against.")
    ctl = float(ctl_rows.delta_nlte.iloc[0])

    lines = []
    for wl in sorted(LINES):
        loggf, ep, eup, vw = _ges_atomic('Mg', wl)
        if abs(loggf - LINES[wl]['loggf']) > 0.005:
            raise SystemExit(f"RYA-592 SSOT: GES synthesis loggf {loggf:+.3f} for Mg {wl} "
                             f"disagrees with canonical {LINES[wl]['loggf']:+.3f} — the "
                             f"Engine-A leg would run on a different gf scale than the fit. "
                             f"Reconcile the stores before deriving.")
        tl, tu, jl, ju = P.auto_labels('Mg', ep, eup)
        print(f"    {wl:.3f}: loggf={loggf:+.3f} EP={ep:.3f} Eup={eup:.3f} vdW={vw:.3f} "
              f"-> '{tl}' -> '{tu}' (J{jl:.0f}->{ju:.0f})")
        lines.append((wl, loggf, ep, jl, eup, ju, tl, tu, vw))

    P.NLTE_LINES['Mg'] = lines
    P._A_SUN.setdefault('Mg', SOLAR_ASPLUND2021['Mg'])
    res = P.nlte_delta('Mg', star={'teff': 5772, 'logg': 4.44, 'feh': 0.0, 'vmic': 1.0})
    per_line = {float(k): float(v) for k, v in res['per_line'].items()}

    resid = per_line[REFLINE] - ctl
    print(f"    CONTROL Engine-A {REFLINE}: re-derived {per_line[REFLINE]:+.4f} vs committed "
          f"{ctl:+.4f}  (residual {resid:+.4f}, tol {CTRL_TOL})")
    if abs(resid) > CTRL_TOL:
        raise SystemExit(
            f"RYA-592 STOP: Engine-A re-derivation does NOT reproduce the committed "
            f"Mg {REFLINE} solar delta ({per_line[REFLINE]:+.4f} vs {ctl:+.4f}). The 5528 "
            f"delta from this run is therefore not comparable to the committed 5711 chain. "
            f"Do NOT ship a value on an unreproduced grid.")
    return per_line, ctl


# ── Engine-B (TS-native Gerber 2023) NLTE delta, from the committed grid ─────

def engine_b_deltas(root):
    """Per-line Engine-B delta from the committed Gerber TS grid, via the RYA-534 deck
    (scripts/ts_gerber_gate) — only the wave list and work dir differ from that gate.
    The GES NLTE rows are used VERBATIM for the level identification (a plain VALD row
    makes bsyn silently drop to LTE, the RYA-533 trap) with only log_gf patched to the
    canonical single-source value. 5711 is the CONTROL against the RYA-534 record."""
    sys.path.insert(0, os.path.join(root, 'scripts'))
    os.environ.setdefault('TSGERBER_REPO_ROOT', root)
    import ts_gerber_gate as G
    G.W = W_GERBER

    cfg = dict(G.ELEMENTS['Mg'])
    waves = sorted(LINES)
    cfg['waves'] = waves

    os.makedirs(f"{G.W}/work", exist_ok=True)
    os.makedirs(f"{G.W}/Testout", exist_ok=True)
    dsl = f"{G.W}/work/DATA"
    if not os.path.lexists(dsl):
        os.symlink("/mnt/codex-data/engines/Turbospectrum_NLTE/DATA", dsl)

    hdr, rows = G.ges_lines(cfg['Z'], 0, waves)
    patched = []
    for w, r in zip(waves, rows):
        m = re.match(r"^(\s*\S+\s+\S+\s+)(\S+)(\s+.*)$", r)
        if not m:
            raise SystemExit(f"RYA-592: cannot parse GES NLTE row for Mg {w}: {r!r}")
        print(f"    GES NLTE row {w}: loggf {m.group(2)} -> canonical {LINES[w]['loggf']:+.3f} "
              f"(level IDs verbatim)")
        patched.append(f"{m.group(1)}{LINES[w]['loggf']:.3f}{m.group(3)}")

    linelist = f"{G.W}/work/Mg.lines"
    mh = re.match(r"^('\s*\d+\.\d\d\d\s*')\s+(\d+)\s+(\d+)", hdr)
    with open(linelist, 'w') as f:
        f.write(f"{mh.group(1)}    {mh.group(2)}      {len(patched)}\n")
        f.write("'Mg I'\n")
        for r in patched:
            f.write(r + "\n")

    dep = G.make_departure('Mg', cfg)
    nlteinfo = f"{G.W}/work/species_Mg.dat"
    with open(nlteinfo, 'w') as f:
        f.write(f"# RYA-592 Mg NLTE\n# path for model atom files\n{G.GT}/\n"
                f"# path for departure files\n{G.W}/Testout/\n# atomic NLTE setup\n"
                f"{cfg['Z']}  'Mg'  'nlte'  '{cfg['atom']}'  'sun_Mg_coef.dat'  'ascii'\n")

    lmin, lmax = min(waves) - G.LMARGIN, max(waves) + G.LMARGIN
    model = f"{G.W}/Testout/sun_Mg.interpol"
    opac = f"{G.W}/work/sun_Mg.opac"
    G.babsma('Mg', model, opac, lmin, lmax)
    res_n, engaged = G.bsyn('Mg', cfg, cfg['a_sun'], True, 'nlte', nlteinfo, linelist,
                            opac, lmin, lmax)
    if not engaged:
        raise SystemExit("RYA-592 STOP: Gerber NLTE departures did NOT engage (silent LTE) — "
                         "the Engine-B delta would be a fabricated zero. Do NOT ship.")
    print(f"    departures engaged: {engaged}  (departure file {os.path.basename(dep)})")

    offs = [-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3]
    A = np.array([cfg['a_sun'] + o for o in offs])
    cog = {c: [] for c in waves}
    for a in A:
        res_l, _ = G.bsyn('Mg', cfg, a, False, f"lte_{a:.2f}", nlteinfo, linelist,
                          opac, lmin, lmax)
        for c in waves:
            cog[c].append(G.ew(res_l, c))

    per_line = {}
    for c in waves:
        ewn = G.ew(res_n, c)
        ewl = np.array(cog[c])
        a_star = float(np.interp(ewn, ewl, A))
        per_line[c] = cfg['a_sun'] - a_star
        print(f"    Mg {c:.3f}: EW_NLTE={ewn:7.2f} mA  A*={a_star:.4f}  "
              f"delta={per_line[c]:+.4f}  (LTE COG {ewl.min():.1f}..{ewl.max():.1f} mA)")

    resid = per_line[REFLINE] - ENGINE_B_COMMITTED_5711
    print(f"    CONTROL Engine-B {REFLINE}: re-derived {per_line[REFLINE]:+.4f} vs RYA-534 "
          f"record {ENGINE_B_COMMITTED_5711:+.4f}  (residual {resid:+.4f}, tol {CTRL_TOL})")
    if abs(resid) > CTRL_TOL:
        raise SystemExit(
            f"RYA-592 STOP: Engine-B re-derivation does NOT reproduce the RYA-534 Mg "
            f"{REFLINE} delta ({per_line[REFLINE]:+.4f} vs {ENGINE_B_COMMITTED_5711:+.4f}). "
            f"Do NOT ship a value on an unreproduced grid.")
    return per_line


# ── LTE in-window blend-fit synthesis (the RYA-551/560 harness) ──────────────

def patch_linelist(center, vald_name):
    """Copy the VALD block for the window and overwrite the target Mg I line's log_gf
    with the canonical value so the engine synthesises on OUR single-source gf scale.
    bsyn selects lines by LAMBDA window; blends stay VERBATIM — that is the whole point
    of the in-window fit (344 components live within +-1.6 A of 5528)."""
    src = os.path.join(VALD_DIR, vald_name)
    dst = f"{W}/ll_{center:.0f}.list"
    loggf = LINES[center]['loggf']
    n_patched = 0
    out = []
    with open(src, errors='replace') as f:
        for ln in f:
            p = ln.split()
            if len(p) >= 3:
                try:
                    wl = float(p[0])
                except ValueError:
                    out.append(ln)
                    continue
                if abs(wl - center) < 0.02 and 'Mg I' in ln:
                    m = re.match(r"^(\s*\S+\s+\S+\s+)(\S+)(\s+.*)$", ln.rstrip("\n"))
                    if m:
                        out.append(f"{m.group(1)}{loggf:.3f}{m.group(3)}\n")
                        n_patched += 1
                        continue
            out.append(ln)
    if n_patched != 1:
        raise SystemExit(f"RYA-592: patched {n_patched} Mg I rows at {center} in {vald_name} "
                         f"(expected exactly 1) — the target line is missing or duplicated in "
                         f"the VALD window. Do NOT synthesise on an ambiguous target.")
    with open(dst, 'w') as f:
        f.writelines(out)
    return dst


def count_blends(linelist, center, hw):
    """How many NON-target components the synthesis actually carries in the fit window.
    The RYA-592 CRITICAL check: if the blend was culled, this collapses and the fit is
    an isolated-line EW in disguise."""
    n = 0
    with open(linelist, errors='replace') as f:
        for ln in f:
            p = ln.split()
            if len(p) < 3 or ln.startswith("'"):
                continue
            try:
                wl = float(p[0])
            except ValueError:
                continue
            if abs(wl - center) <= hw and abs(wl - center) > 0.001:
                n += 1
    return n


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
           f"'INDIVIDUAL ABUNDANCES:'   '1'\n{Z_MG}  {a_x:.3f}\n'ISOTOPES : ' '0'\n"
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
        gk = np.exp(-0.5 * (gx / gsig_A) ** 2)
        gk /= gk.sum()
        out = np.convolve(1 - out, gk, mode='same')
        out = 1 - out
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
        r = csv.reader(fh)
        next(r)
        for row in r:
            w.append(float(row[0]))
            f.append(float(row[3]))
    return np.array(w), np.array(f)


def load_iag(root):
    """IAG Reiners+2016 visible FTS solar flux atlas: cols = vacuum wavenumber(cm^-1),
    normalized flux. Vacuum->air via the pipeline SINGLE-SOURCE converter."""
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
    wn = np.array(wn)
    fl = np.array(fl)
    lam_air = vac_to_air(1e8 / wn)
    idx = np.argsort(lam_air)
    return lam_air[idx], fl[idx]


def measure_arm_rv(obs_w, obs_f):
    """Residual velocity of an arm from parabolic core centroids of clean, unblended
    solar lines (RYA-309 _measure_rv_kms in spirit). Abundance-blind: it is the
    independent check that the fitted nuisance dv is a real wavelength offset and not
    the fit absorbing a profile mismatch."""
    vs = []
    for lam0 in RV_CHECK_LINES:
        m = (obs_w > lam0 - 0.12) & (obs_w < lam0 + 0.12)
        if m.sum() < 8:
            continue
        x, y = obs_w[m], obs_f[m]
        k = int(np.argmin(y))
        sl = slice(max(k - 3, 0), min(k + 4, len(x)))
        if sl.stop - sl.start < 4:
            continue
        c = np.polyfit(x[sl], y[sl], 2)
        if c[0] <= 0:
            continue
        vs.append((-c[1] / (2 * c[0]) - lam0) / lam0 * CLIGHT)
    if not vs:
        return None, 0, 0.0
    return float(np.median(vs)), len(vs), float(np.std(vs))


def fit_line(center, obs_w, obs_f, synth):
    """chi2 profile fit of A(Mg) over the fit window (NOT an isolated-line EW inversion).

    Free parameters: A(Mg) [the measurement], plus two nuisance parameters — gsig (vmac
    + instrumental Gaussian) and dv (residual velocity, RYA-309: the fit is against
    REST-frame synthetic flux, and this solar arm carries ~+0.8 km/s). Returns the
    fitted A_LTE, the quality metrics the reliability gate needs, and the OBSERVED EW
    over the fit window (used for the ratified saturation-knee routing)."""
    hw = LINES[center]['fit_hw']
    ww, ff, _ = local_renorm(obs_w, obs_f, center, hw)
    pad = (ww > center - hw - 0.05) & (ww < center + hw + 0.05)
    if pad.sum() < 10:
        return None
    wp, fp = ww[pad], ff[pad]
    best = dict(chi2=1e30)
    for a, (sw, sf) in synth.items():
        for gs in GSIG_GRID:
            sb = broaden(sw, sf, VSINI, gs)
            for dv in DV_GRID:
                xo = wp / (1 + dv / CLIGHT)      # shift the OBSERVED to rest
                sel = (xo > center - hw) & (xo < center + hw)
                if sel.sum() < 10:
                    continue
                r = fp[sel] - np.interp(xo[sel], sw, sb)
                chi2 = float(np.sum(r ** 2)) / max(int(sel.sum()) - 3, 1) / (0.01 ** 2)
                if chi2 < best['chi2']:
                    best = dict(chi2=chi2, A=float(a), gsig=float(gs), dv=float(dv),
                                npix=int(sel.sum()))
    gs, dv = best['gsig'], best['dv']
    xo = wp / (1 + dv / CLIGHT)
    sel = (xo > center - hw) & (xo < center + hw)
    xo, yo = xo[sel], fp[sel]
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
    best['gsig_railed'] = bool(gs <= GSIG_GRID[0] + 1e-9 or gs >= GSIG_GRID[-1] - 1e-9)
    best['obs_ew_mA'] = float(np.trapz(1 - yo, xo) * 1000.0)

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


# ── main ────────────────────────────────────────────────────────────────────

def main():
    # RYA-567: Sirius-only heavy-compute leg (TS Turbospectrum + MARCS + the NLTE grids).
    # Refuse to run off Sirius — loud-fail, never against local-Mac copies.
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sys.path.insert(0, root)
    from config.constants import assert_on_sirius, SOLAR_ASPLUND2021, TWO_ENGINE
    assert_on_sirius("RYA-592 Mg I 5528 second-line synthesis",
                     require_subdirs=("engines", "grids"))

    a_sun = float(SOLAR_ASPLUND2021['Mg'])   # 7.55 — reference point, NOT a tuning target
    knee = float(TWO_ENGINE['saturation_knee_mA'])
    print(f"RYA-592 Mg I 5528 second line — A(Mg)_sun reference (Asplund 2021) = {a_sun} "
          f"[grid-bracket + reconciliation reference only, NEVER a fit target]\n")

    print("Step 0 — SSOT gf check (data/linelists/canonical_gf.csv):")
    verify_canonical_gf(root)

    print("\nStep 1a — Engine-A (Amarsi2020 PySME) NLTE deltas from the committed .grd:")
    delta_a, ctrl_a = engine_a_deltas(root)

    print("\nStep 1b — Engine-B (TS-native Gerber 2023) NLTE deltas from the committed grid:")
    delta_b = engine_b_deltas(root)

    os.makedirs(W, exist_ok=True)
    if not os.path.lexists(f"{W}/DATA"):
        os.symlink("/mnt/codex-data/engines/Turbospectrum_NLTE/DATA", f"{W}/DATA")

    arms = {'harps': load_harps()}
    try:
        arms['iag'] = load_iag(root)
        print(f"\n  IAG atlas: {len(arms['iag'][0])} pts, "
              f"{arms['iag'][0].min():.1f}-{arms['iag'][0].max():.1f} A")
    except Exception as e:
        print(f"\n  IAG load FAILED ({e}); HARPS-only")

    # Independent, abundance-blind residual-velocity measurement per arm (RYA-309).
    arm_rv = {}
    print("\n  residual velocity per arm (clean-line core centroids, abundance-blind):")
    for arm, (ow, of) in arms.items():
        v, n, sd = measure_arm_rv(ow, of)
        arm_rv[arm] = dict(v_kms=(round(v, 3) if v is not None else None), n_lines=n,
                           scatter_kms=round(sd, 3))
        print(f"    {arm:6s}: {v:+.3f} km/s (n={n}, scatter {sd:.3f})" if v is not None
              else f"    {arm:6s}: not measurable")

    results = {'_meta': dict(
        ticket='RYA-592', element='Mg', ion='I', z=Z_MG,
        a_sun_ref=a_sun, a_sun_ref_source='Asplund2021',
        method=('in-window blend-fit synthesis (chi2 over the fit window vs the observed '
                'solar profile, ALL VALD3 in-window components synthesised) — NOT an '
                'isolated-line EW inversion'),
        target_line=TARGET, concordance_reference_line=REFLINE,
        reliable_dEW_dA_floor=RELIABLE_DEWDA,
        saturation_knee_mA=knee,
        saturation_knee_source="config.constants.TWO_ENGINE['saturation_knee_mA'] (RYA-525 clause-3)",
        engine_a=dict(name='Mg_Amarsi2020_PySME', grid='nlte_Mg_scatt_pysme.grd',
                      control_line=REFLINE, control_committed=ctrl_a,
                      control_rederived=delta_a[REFLINE], control_tol=CTRL_TOL),
        engine_b=dict(name='Mg_gerber2023_TSnative', grid='NLTEgrid4TS_Mg_MARCS_Nov-13-2024.bin',
                      control_line=REFLINE, control_committed=ENGINE_B_COMMITTED_5711,
                      control_rederived=delta_b[REFLINE], control_tol=CTRL_TOL),
        gate3_note=('RYA-561 STRICT: the Engine-A/Engine-B agreement recorded here is an '
                    'ATOM-delta agreement on a shared profile fit — DIAGNOSTIC ONLY. It is '
                    'NOT a cross-engine dCE and must never satisfy gate 3.'),
        scale_caveat=('1D-NLTE value vs 3D-NLTE reference; un-applied 3D term folded into '
                      'the offset (RYA-561; class-wide 3D-metals fix is post-Beta).'))}

    for center, cfg in LINES.items():
        print(f"\n=== Mg I {center} ({cfg['role']}; EP {cfg['ep']}, canonical loggf "
              f"{cfg['loggf']}, {cfg['gf_ref']}) ===")
        ll = patch_linelist(center, cfg['vald'])
        n_blend = count_blends(ll, center, cfg['fit_hw'])
        print(f"  in-window blend components (fit window +-{cfg['fit_hw']} A): {n_blend}")
        if n_blend < 5:
            raise SystemExit(
                f"RYA-592 CRITICAL: only {n_blend} non-target components inside the {center} "
                f"fit window — the blend was CULLED, so this fit would be an isolated-line EW "
                f"in disguise (the exact defect this ticket exists to remove). Do NOT ship.")
        lmin, lmax = center - LMARGIN, center + LMARGIN
        opac = f"{W}/opac_{center:.0f}"
        babsma(center, opac, lmin, lmax)
        synth = {}
        for a in A_GRID:
            synth[float(a)] = bsyn(center, a, ll, opac, lmin, lmax, f"{a:.2f}")

        rec = dict(wave=center, role=cfg['role'], ep=cfg['ep'], loggf=cfg['loggf'],
                   gf_id=cfg['gf_id'], gf_ref=cfg['gf_ref'], fit_hw=cfg['fit_hw'],
                   n_blend_components_in_window=n_blend,
                   nlte_delta_engineA=round(delta_a[center], 4),
                   nlte_delta_engineB=round(delta_b[center], 4),
                   nlte_delta_atom_spread=round(abs(delta_a[center] - delta_b[center]), 4))
        for arm, (ow, of) in arms.items():
            fit = fit_line(center, ow, of, synth)
            if fit is None:
                rec[arm] = dict(status='no_coverage')
                print(f"  {arm:6s}: no coverage")
                continue
            a_lte = fit['A']
            reliable = bool((not fit['railed']) and fit['dEW_dA'] >= RELIABLE_DEWDA)
            rec[arm] = dict(
                A_LTE=round(a_lte, 3),
                A_NLTE_engineA=round(a_lte + delta_a[center], 3),
                A_NLTE_engineB=round(a_lte + delta_b[center], 3),
                gsig_kms=round(fit['gsig'], 2), gsig_railed=fit['gsig_railed'],
                dv_fitted_kms=round(fit['dv'], 2),
                dv_measured_kms=arm_rv.get(arm, {}).get('v_kms'),
                red_chi2=round(fit['red_chi2'], 2),
                npix=fit['npix'], obs_ew_mA=round(fit['obs_ew_mA'], 1),
                dEW_dA_mA_dex=fit['dEW_dA'], railed=fit['railed'], reliable=reliable,
                ew_route_admissible=bool(fit['obs_ew_mA'] <= knee))
            print(f"  {arm:6s}: A_LTE={a_lte:.3f}  A_NLTE(A)={a_lte+delta_a[center]:.3f}  "
                  f"A_NLTE(B)={a_lte+delta_b[center]:.3f}")
            print(f"          gsig={fit['gsig']:.2f} km/s (railed={fit['gsig_railed']}), "
                  f"dv_fit={fit['dv']:+.2f} vs dv_measured="
                  f"{arm_rv.get(arm, {}).get('v_kms')} km/s, rchi2={fit['red_chi2']:.2f}, "
                  f"obsEW={fit['obs_ew_mA']:.1f} mA, dEW/dA={fit['dEW_dA']} mA/dex, "
                  f"railed={fit['railed']}, reliable={reliable}")
        results[str(center)] = rec

    # ── Step 2: concordance + the ratified routing verdict ───────────────────
    print("\n" + "=" * 74)
    conc = {}
    for arm in arms:
        t = results[str(TARGET)].get(arm, {})
        r = results[str(REFLINE)].get(arm, {})
        if not (isinstance(t, dict) and isinstance(r, dict)
                and t.get('A_LTE') is not None and r.get('A_LTE') is not None):
            continue
        conc[arm] = {eng: round(t[f'A_NLTE_engine{eng}'] - r[f'A_NLTE_engine{eng}'], 3)
                     for eng in ('A', 'B')}
        conc[arm]['both_reliable'] = bool(t.get('reliable') and r.get('reliable'))
        print(f"CONCORDANCE {arm}: A({TARGET}) - A({REFLINE}) = "
              f"{conc[arm]['A']:+.3f} (Engine-A atom) / {conc[arm]['B']:+.3f} (Engine-B atom); "
              f"both lines reliable = {conc[arm]['both_reliable']}")
    results['_concordance'] = conc

    tgt = results[str(TARGET)].get('harps', {})
    ref = results[str(REFLINE)].get('harps', {})
    reliable = bool(isinstance(tgt, dict) and tgt.get('reliable'))
    ew_ok = bool(isinstance(tgt, dict) and tgt.get('ew_route_admissible'))
    # Concordance verdict: worst |delta| across arms x engine-atoms, vs the pre-declared
    # band. Both lines must themselves be reliable for the comparison to mean anything.
    spans = [abs(v) for arm in conc for k, v in conc[arm].items()
             if k in ('A', 'B') and conc[arm].get('both_reliable')]
    conc_worst = max(spans) if spans else None
    concordant = bool(conc_worst is not None and conc_worst <= CONCORDANCE_BAND)
    results['_verdict'] = dict(
        target_reliable=reliable,
        target_A_NLTE_engineB=tgt.get('A_NLTE_engineB'),
        reference_A_NLTE_engineB=ref.get('A_NLTE_engineB'),
        concordance_band=CONCORDANCE_BAND,
        concordance_worst_abs_dex=(round(conc_worst, 3) if conc_worst is not None else None),
        lines_concordant=concordant,
        engineA_ew_route_admissible=ew_ok,
        engineA_ew_route_reason=(
            f"observed EW {tgt.get('obs_ew_mA')} mA vs ratified saturation knee {knee} mA "
            f"(TWO_ENGINE clause-3: EW above the knee => SATURATED => HARD => Engine-B). "
            f"{'ADMISSIBLE' if ew_ok else 'NOT admissible — Engine-A (EW->A) cannot run on this line'}"),
        real_cross_engine_dCE_available=ew_ok,
        gate3_status=(
            'a REAL dCE is possible from this line' if ew_ok else
            'NO real dCE from this line (saturated => Engine-B by the ratified clause-3 '
            'classifier), so under the RYA-561 STRICT rule gate 3 FAILS. The separately-'
            'flagged >=2-concordant-line variant is NOT ratified AND would fail anyway: '
            f'the two lines disagree by {conc_worst:.3f} dex (band {CONCORDANCE_BAND}). '
            'Mg holds CURATION-OWED-with-value.'
            if conc_worst is not None and not concordant else
            'NO real dCE from this line; gate 3 FAILS under the RYA-561 STRICT rule. The '
            'lines DO concord, so the (unratified) >=2-concordant-line variant would be '
            'satisfied — flag for ratification, do NOT self-apply.'),
        disposition='CURATION-OWED-with-value (unchanged)')

    out_dir = os.path.join(root, 'data', 'results')
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, 'mg_5528_synthesis_rya592.json')
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)

    print("-" * 74)
    if reliable:
        print(f"Mg I {TARGET} CLEARED the reliability floor "
              f"(dEW/dA >= {RELIABLE_DEWDA} mA/dex, not railed) -> a real second line.")
    else:
        print(f"Mg I {TARGET} did NOT clear the reliability floor -> Mg stays owed, single-line. "
              f"Do NOT emit a weak/railed value.")
    if conc_worst is not None:
        print(f"CONCORDANCE VERDICT: worst |A({TARGET}) - A({REFLINE})| = {conc_worst:.3f} dex "
              f"vs band {CONCORDANCE_BAND} -> {'CONCORDANT' if concordant else 'DISCORDANT'}")
    print(f"Engine-A EW route on {TARGET}: {results['_verdict']['engineA_ew_route_reason']}")
    print(f"gate 3: {results['_verdict']['gate3_status']}")
    print(f"Mg disposition: {results['_verdict']['disposition']}")
    print(f"\nWrote {out_path}")


if __name__ == '__main__':
    main()
