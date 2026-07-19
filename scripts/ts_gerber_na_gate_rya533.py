#!/usr/bin/env python3
"""RYA-533 — TS-native Gerber NLTE synthesis deck: the Na gate (RUNS ON SIRIUS).

Engine-B NLTE for the two-engine floor (RYA-525). This is the first TS-native
in-synthesis NLTE reproduction in the project — the leg RYA-529 skipped (529 did the
PySME/Engine-A leg). babsma once (continuum opacity) -> bsyn NLTE with the Gerber Na
departures + model atom -> bsyn LTE on a curve of growth -> invert EW_NLTE onto the
LTE COG per line -> delta = A(NLTE) - A(LTE) = A_sun - A*.

VALIDATED (validate-don't-tune, nothing fitted):
    Na I 5682.633  delta -0.065   (EW_NLTE 112.5 mA)
    Na I 5688.205  delta -0.072   (EW_NLTE 147.8 mA)
    median        -0.068   vs anchor -0.107 +/- 0.05 (model-atom tol) -> PASS
Cross-engine (feeds RYA-525's model-atom-systematic diagnostic): INSPECT -0.107 /
PySME/Amarsi -0.129 (RYA-529) / TS-Gerber -0.068 (here). The ~0.04 dex offset from
INSPECT is the Gerber (atom.na_qmh) vs Lind-2011 model-atom difference, exactly the
"path-correctness, ~0.02-0.05 dex" the ticket predicted.

PROVISIONING (Sirius only, RYA-526; md5-pinned in data/nlte_grids/gerber_ts/Na_gerber2023.prov.json):
  * engine   : /mnt/codex-data/engines/Turbospectrum_NLTE (BSYN v20.1, exec-gf, gfortran 15.2)
  * grid+atom: /mnt/codex-data/grids/nlte/gerber_ts/  (Gerber 2023 A&A 669 A43; NLTEgrid4TS_Na
               .bin 15.9 GB + atom.na_qmh + aux; from the MPG Keeper Seafile share)
  * MARCS    : /mnt/codex-data/grids/model_atmospheres/marcs_standard_comp/ (solar node p5750_g+4.5_z+0.00)
  * departure: interpol_modeles_nlte -> sun_Na_coef.dat (4D interp over Teff/logg/[Fe/H]/A(Na))

GOTCHAS baked in (each was a real failure while wiring this):
  1. interpol_modeles_nlte reads the aux ROW COUNT from stdin (the script's "aux_file_length"
     is NOT columns) -> must pass wc-l-1 of the aux (436255 for Na), else it reads 9 rows and
     "no match found".
  2. babsma/bsyn need a DATA/ dir relative to cwd (atomicweights.dat) -> symlink the TS DATA dir.
  3. NLTEINFOFILE must carry the exact marker lines "# path for model atom files" /
     "# path for departure files" (read_nlteinfofile.f), else "NLTE info file seems totally empty".
  4. THE SILENT-ERROR SITE (RYA-402's warning): the line list must carry NLTE LEVEL IDENTIFICATION
     (level indices + labels matching the model atom). A plain VALD line -> bsyn logs "missing level
     identification" and silently sets departure = 1 (LTE). Fixed by taking the 5682/5688 lines
     VERBATIM from the bundled GES NLTE line list (nlte_ges_linelist_jmg17feb2022_I_II), which
     carries "... 2 11 '3p2P0*' '4d2D1' 'c' 'a'". The driver RAISES if departures don't engage.
"""
import subprocess, os, sys, numpy as np

# RYA-567: Sirius-only heavy-compute leg (TS engine + Gerber grids). Refuse to run it
# off Sirius — loud-fail, never against local-Mac copies.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config.constants import assert_on_sirius
assert_on_sirius("RYA-533 TS-Gerber Na gate", require_subdirs=("engines", "grids"))

EXE = "/mnt/codex-data/engines/Turbospectrum_NLTE/exec-gf"
W   = "/mnt/codex-data/codex/rya533"
GT  = "/mnt/codex-data/grids/nlte/gerber_ts"
os.makedirs(f"{W}/work", exist_ok=True)
os.chdir(f"{W}/work")

MODEL = f"{W}/Testout/sun_Na.interpol"
OPAC  = f"{W}/work/sun_Na.opac"
A_SUN = 6.24
LMIN, LMAX, DLAM = 5678.0, 5692.0, 0.005

# --- Na I line list (TS format): 5682.633 / 5688.205 subordinate doublet 3p->4d ---
# cols: wl  Elow_eV  loggf  vdW(ABO packed)  gu  gamrad  'lo' 'up'  0.0 1.0 'label'
LINELIST = f"{W}/work/na.lines"
# Na I 5682/5688 taken VERBATIM from the bundled GES NLTE linelist
# (nlte_ges_linelist_jmg17feb2022_I_II) — carries the level-identification
# fields bsyn matches to atom.na_qmh (indices 2/3 -> 11, labels '3p2P*'/'4d2D1').
with open(LINELIST, "w") as f:
    f.write("'  11.000            '    1        2\n")
    f.write("'Na I     NLTE'\n")
    f.write("  5682.633  2.102 -0.706 1955.327    4.0  1.00E+05  'p' 'd'   0.0    1.0 'Na I LS:2p6.3p 2P* LS:2p6.4d 2D'  2 11  '3p2P0*' '4d2D1'  'c' 'a'\n")
    f.write("  5688.205  2.104 -0.404 1955.327   10.0  1.00E+05  'p' 'd'   0.0    1.0 'Na I LS:2p6.3p 2P* LS:2p6.4d 2D'  3 11  '3p2P1*' '4d2D1'  'c' 'a'\n")

# --- NLTEINFOFILE (Na nlte via Gerber atom + interpolated departure file) ---
NLTEINFO = f"{W}/work/species_Na.dat"
with open(NLTEINFO, "w") as f:
    f.write("# RYA-533 Na NLTE\n")
    f.write("# path for model atom files\n")
    f.write(f"{GT}/\n")
    f.write("# path for departure files\n")
    f.write(f"{W}/Testout/\n")
    f.write("# atomic NLTE setup\n")
    f.write("11  'Na'  'nlte'  'atom.na_qmh'  'sun_Na_coef.dat'  'ascii'\n")

def babsma():
    ctl = f"""'LAMBDA_MIN:'  '{LMIN}'
'LAMBDA_MAX:'  '{LMAX}'
'LAMBDA_STEP:' '{DLAM}'
'MODELINPUT:' '{MODEL}'
'MARCS-FILE:' '.false.'
'MODELOPAC:' '{OPAC}'
'METALLICITY:'    '0.00'
'ALPHA/Fe   :'    '0.00'
'HELIUM     :'    '0.00'
'R-PROCESS  :'    '0.00'
'S-PROCESS  :'    '0.00'
'XIFIX:' 'T'
1.0
"""
    r = subprocess.run([f"{EXE}/babsma_lu"], input=ctl, capture_output=True, text=True)
    if not os.path.exists(OPAC):
        print("BABSMA FAILED\n", r.stdout[-2000:], r.stderr[-1000:]); sys.exit(1)

def bsyn(a_na, nlte, tag):
    res = f"{W}/work/na_{tag}.spec"
    ctl = f"""'NLTE :'          '{'.true.' if nlte else '.false.'}'
'NLTEINFOFILE:'  '{NLTEINFO}'
'LAMBDA_MIN:'     '{LMIN}'
'LAMBDA_MAX:'     '{LMAX}'
'LAMBDA_STEP:'    '{DLAM}'
'INTENSITY/FLUX:' 'Flux'
'ABFIND        :' '.false.'
'MODELOPAC:' '{OPAC}'
'RESULTFILE :' '{res}'
'METALLICITY:'    '0.00'
'ALPHA/Fe   :'    '0.00'
'HELIUM     :'    '0.00'
'R-PROCESS  :'    '0.00'
'S-PROCESS  :'    '0.00'
'INDIVIDUAL ABUNDANCES:'   '1'
11  {a_na:.3f}
'ISOTOPES : ' '0'
'NFILES   :' '1'
{LINELIST}
'SPHERICAL:'  'F'
  30
  300.00
  15
  1.30
"""
    r = subprocess.run([f"{EXE}/bsyn_lu"], input=ctl, capture_output=True, text=True)
    so = r.stdout
    if nlte:
        read_dep = "read departure coefficients" in so
        fell_back = ("missing level identification" in so) or \
                    ("continuing with departure coefficient = 1. for that level" in so)
        engaged = read_dep and not fell_back
    else:
        engaged = True
    if not os.path.exists(res):
        print(f"BSYN FAILED tag={tag}\n", r.stdout[-2500:], r.stderr[-800:]); sys.exit(1)
    return res, r.stdout, engaged

def ew_mA(specfile, center, hw=1.0):
    d = np.loadtxt(specfile)
    w, fl = d[:,0], d[:,1]
    m = (w > center-hw) & (w < center+hw)
    return np.trapz(1.0 - fl[m], w[m]) * 1000.0

centers = [5682.633, 5688.205]
print("=== babsma (continuum opacity) ===")
babsma()
print("=== bsyn NLTE @ A(Na)=6.24 ===")
res_n, out_n, engaged = bsyn(A_SUN, True, "nlte")
print("  NLTE departure engaged:", engaged)
if not engaged:
    print("RAISE: NLTE departures not engaged (silent LTE) — STOP"); sys.exit(2)
ew_nlte = {c: ew_mA(res_n, c) for c in centers}

print("=== bsyn LTE COG (per line) ===")
offs = [-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3]
A = np.array([A_SUN+o for o in offs])
cog = {c: [] for c in centers}
for a in A:
    res_l,_,_ = bsyn(a, False, f"lte_{a:.2f}")
    for c in centers:
        cog[c].append(ew_mA(res_l, c))

print("\n=== RESULT (per line: delta = A_sun - A*) ===")
deltas = {}
for c in centers:
    ewl = np.array(cog[c])
    a_star = float(np.interp(ew_nlte[c], ewl, A))
    deltas[c] = A_SUN - a_star
    print(f"  Na I {c:.3f}: EW_NLTE={ew_nlte[c]:.2f} mA  A*={a_star:.4f}  delta={deltas[c]:+.4f}")
med = float(np.median(list(deltas.values())))
print(f"  median delta = {med:+.4f}")
print(f"  anchor -0.107 (INSPECT/Amarsi); model-atom tol 0.05 → PASS if |med-(-0.107)|<=0.05")
print(f"  VERDICT: {'PASS' if abs(med+0.107)<=0.05 else 'CHECK'}")
