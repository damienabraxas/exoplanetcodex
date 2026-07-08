#!/usr/bin/env python3
"""
RYA-534 — element-general TS-native Gerber NLTE gate (RUNS ON SIRIUS).

Generalises the RYA-533 Na deck to any Family-A element. Per element:
  interpol_modeles_nlte (Gerber departures, 4D interp -> solar departure file)
  -> babsma (contopac) -> bsyn NLTE + bsyn LTE curve-of-growth
  -> per-line EW inversion -> delta = A(NLTE) - A(LTE) = A_sun - A*.
Reproduce the element's published solar NLTE anchor (validate-don't-tune).

Line lists are pulled from the bundled GES NLTE line list BY WAVELENGTH so the
level identification is verbatim (RYA-533's silent-departure=1 trap: a plain VALD
line makes bsyn drop to LTE with no error). The deck RAISES if departures don't engage.

Usage (on Sirius):  venv_pysme/bin/python ts_gerber_gate.py <ELEMENT>
All heavy artifacts live on Sirius only (RYA-526); grids md5-pinned in the prov JSON.
"""
import os
import re
import subprocess
import sys
import numpy as np

EXE = "/mnt/codex-data/engines/Turbospectrum_NLTE/exec-gf"
INTERP = "/mnt/codex-data/engines/Turbospectrum_NLTE/interpolator/interpol_modeles_nlte"
GES = "/mnt/codex-data/engines/Turbospectrum_NLTE/COM/linelists/nlte_ges_linelist_jmg17feb2022_I_II"
GT = "/mnt/codex-data/grids/nlte/gerber_ts"
MARCS = ("/mnt/codex-data/grids/model_atmospheres/marcs_standard_comp/marcs_standard_comp/"
         "p5750_g+4.5_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod")
W = "/mnt/codex-data/codex/rya534"
TREF, LOGGREF, ZREF = "5750", "+4.50", "+0.00"

# Asplund 2021 solar A(X); anchors = published solar 1D-NLTE corrections (sign+magnitude).
ELEMENTS = {
    'O':  dict(Z=8,  a_sun=8.69, atom='atom.o41f',   aux='auxData_O_MARCS_May-21-2021.txt',
               grid='NLTEgrid4TS_O_MARCS_May-21-2021.bin',
               waves=[7771.940, 7774.170, 7775.390], anchor=-0.134, tol=0.05,
               ref='O I 777, 1D-NLTE; anchor = our banked Amarsi-2019 nlte_cno 1D leg (table6) solar -0.134 (RYA-362 cross-check; the 3D leg is -0.169)'),
    'Mg': dict(Z=12, a_sun=7.55, atom='atom.mg86d',  aux='auxData_Mg_MARCS_Nov-13-2024.dat',
               grid='NLTEgrid4TS_Mg_MARCS_Nov-13-2024.bin',
               waves=[5711.088], anchor=-0.02, tol=0.05,
               ref='Mg I 5711, solar NLTE ~0 in FGK (our Mg_Amarsi2020_PySME -0.022)'),
    'Si': dict(Z=14, a_sun=7.51, atom='atom.si340',  aux='auxData_Si_MARCS_Feb-13-2022.dat',
               grid='NLTEgrid4TS_Si_MARCS_Feb-13-2022.bin',
               waves=[5772.146], anchor=-0.01, tol=0.05,
               ref='Si I 5772, solar NLTE ~0 (our Si_Amarsi2020_PySME -0.013)'),
    'Ca': dict(Z=20, a_sun=6.30, atom='atom.ca105b', aux='auxData_Ca_MARCS_Jun-02-2021.dat',
               grid='NLTEgrid4TS_Ca_MARCS_Jun-02-2021.bin',
               waves=[6122.217, 6162.173], anchor=+0.02, tol=0.05,
               ref='Ca I 6122/6162, solar NLTE small+ (our Ca_Mashonkina2017 +0.017; same atom lineage)'),
    'Ti': dict(Z=22, a_sun=4.90, atom='atom.ti503b', aux='auxData_TI_MARCS_Feb-21-2022.dat',
               grid='NLTEgrid4TS_TI_MARCS_Feb-21-2022.bin',
               waves=[5702.660, 5703.898], anchor=+0.05, tol=0.06,
               ref='Ti I 5702/5703 (both GES levels identified), solar NLTE small+ (Bergemann 2011 MPIA lineage)'),
    'Mn': dict(Z=25, a_sun=5.42, atom='atom.mn281kbc', aux='auxData_MN_MARCS_Mar-15-2023.dat',
               grid='NLTEgrid4TS_MN_MARCS_Mar-15-2023.bin',
               waves=[6013.510, 6021.800], anchor=+0.10, tol=0.06,
               ref='Mn I 6013/6021, solar NLTE +~0.1 (Bergemann&Gehren 2008 / our Mn_Bergemann_MPIA +0.107; same atom lineage)'),
    # Co/Ni line picks CORRECTED (RYA-534): the first-pass lines had an UNIDENTIFIED
    # upper level in the GES list ('none'/'x') -> bsyn silently set departure=1 (the deck
    # RAISED). These lines carry BOTH levels identified ('c'/'a'/'c').
    'Co': dict(Z=27, a_sun=4.94, atom='atom.co247qm', aux='auxData_CO_MARCS_Nov-14-2024.dat',
               grid='NLTEgrid4TS_CO_MARCS_Nov-14-2024.bin',
               waves=[5000.875, 5013.324], anchor=+0.10, tol=0.12,
               ref='Co I 5000/5013 (both GES levels identified), solar NLTE positive (Bergemann+2010 MNRAS 401 1334)'),
    'Ni': dict(Z=28, a_sun=6.20, atom='atom.ni538qm', aux='auxData_Ni_MARCS_Jan-21-2022.txt',
               grid='NLTEgrid4TS_Ni_MARCS_Jan-31-2022.bin',
               waves=[5018.282, 5035.972], anchor=+0.02, tol=0.05,
               ref='Ni I 5018/5035 (both GES levels identified), solar NLTE small+ (Bergemann+2021)'),
    # Sr/Ba: II (singly-ionised) — ges_lines ion=1 (GES stage 2). Line picks CORRECTED to
    # those actually present + identified in the GES Sr II / Ba II NLTE blocks.
    'Sr': dict(Z=38, a_sun=2.83, atom='atom.sr191',  aux='auxData_Sr_MARCS_Mar-15-2023.dat',
               grid='NLTEgrid4TS_Sr_MARCS_Mar-15-2023.bin',
               waves=[4215.519], anchor=-0.01, tol=0.06,
               ref='Sr II 4215 resonance (GES-identified; 4077 absent from GES block), near-LTE at solar (our Sr_Bergemann2012_INSPECT ~-0.005)'),
    'Ba': dict(Z=56, a_sun=2.27, atom='atom.ba111',  aux='auxData_Ba_MARCS_Jul-29-2023.dat',
               grid='NLTEgrid_Ba_MARCS_Jul-29-2023.bin',
               waves=[4554.029], anchor=-0.05, tol=0.10,
               ref='Ba II 4554 resonance (GES-identified; 5853/6141/6496 absent from GES block), solar NLTE moderate negative (Korotin+2015)'),
}

ION = {'Sr': 1, 'Ba': 1}   # singly-ionised diagnostics; the rest are neutral (I)

LMARGIN = 7.0   # window half-width (A) for babsma/bsyn around the lines
EW_HW = 1.2     # per-line EW integration half-width (A)


def sh(cmd, **kw):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True, **kw)


def ges_lines(Z, ion, waves, tol=0.02):
    """Extract the GES NLTE lines for element Z at ionisation `ion` (0=I, 1=II)
    nearest each requested wave, VERBATIM (keeps the level-identification fields).
    NOTE: the GES header encodes the species as `'  Z.000  '  STAGE  NLINES`, where
    STAGE = ion+1 is a SEPARATE middle field (1=neutral, 2=singly-ionised) — NOT the
    decimals of the species code. Returns (header_species, rows)."""
    stage = ion + 1
    rows = []
    with open(GES) as f:
        lines = f.readlines()
    hdr_re = re.compile(r"^'\s*(\d+)\.\d\d\d\s*'\s+(\d+)\s+(\d+)")
    blk_start = None
    for i, ln in enumerate(lines):
        m = hdr_re.match(ln)
        if m and int(m.group(1)) == Z and int(m.group(2)) == stage:
            blk_start = i
            break
    if blk_start is None:
        raise SystemExit(f"no GES NLTE block for Z={Z} stage={stage} (ion={ion})")
    hdr_species = lines[blk_start].rstrip("\n")
    # scan rows until the next species header
    body = []
    for ln in lines[blk_start + 2:]:
        if hdr_re.match(ln):
            break
        body.append(ln)
    for w in waves:
        best, bestd = None, 1e9
        for ln in body:
            try:
                wl = float(ln.split()[0])
            except (ValueError, IndexError):
                continue
            if abs(wl - w) < bestd:
                bestd, best = abs(wl - w), ln.rstrip("\n")
        if best is None or bestd > tol:
            raise SystemExit(f"Z={Z}: no GES line within {tol} A of {w} (nearest d={bestd:.3f})")
        rows.append(best)
    return hdr_species, rows


def make_departure(el, cfg):
    dep = f"{W}/Testout/sun_{el}_coef.dat"
    if os.path.exists(dep):
        return dep
    os.makedirs(f"{W}/Testout", exist_ok=True)
    binf = f"{GT}/{cfg['grid']}"
    aux = f"{GT}/{cfg['aux']}"
    if not os.path.exists(binf):
        raise SystemExit(f"grid .bin missing (provision first): {binf}")
    # aux row count = total lines minus leading comment/header lines
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
        f"'{W}/Testout/sun_{el}.interpol'", f"'{W}/Testout/sun_{el}.alt'", f"'{dep}'",
        f"'{binf}'", f"'{aux}'", str(nrows),
        TREF, LOGGREF, ZREF, f"{cfg['a_sun']:.2f}", ".false.", ".false.", "'none'", ""])
    r = subprocess.run([INTERP], input=stdin, capture_output=True, text=True,
                       cwd=f"{W}/work")
    if not os.path.exists(dep):
        raise SystemExit(f"interpolator failed for {el}:\n{r.stdout[-1500:]}\n{r.stderr[-500:]}")
    return dep


def babsma(el, model, opac, lmin, lmax):
    ctl = (f"'LAMBDA_MIN:'  '{lmin}'\n'LAMBDA_MAX:'  '{lmax}'\n'LAMBDA_STEP:' '0.005'\n"
           f"'MODELINPUT:' '{model}'\n'MARCS-FILE:' '.false.'\n'MODELOPAC:' '{opac}'\n"
           f"'METALLICITY:'    '0.00'\n'ALPHA/Fe   :'    '0.00'\n'HELIUM     :'    '0.00'\n"
           f"'R-PROCESS  :'    '0.00'\n'S-PROCESS  :'    '0.00'\n'XIFIX:' 'T'\n1.0\n")
    subprocess.run([f"{EXE}/babsma_lu"], input=ctl, capture_output=True, text=True, cwd=f"{W}/work")
    if not os.path.exists(opac):
        raise SystemExit(f"babsma failed for {el}")


def bsyn(el, cfg, a_x, nlte, tag, nlteinfo, linelist, opac, lmin, lmax):
    res = f"{W}/work/{el}_{tag}.spec"
    ctl = (f"'NLTE :'          '{'.true.' if nlte else '.false.'}'\n"
           f"'NLTEINFOFILE:'  '{nlteinfo}'\n"
           f"'LAMBDA_MIN:'     '{lmin}'\n'LAMBDA_MAX:'     '{lmax}'\n'LAMBDA_STEP:'    '0.005'\n"
           f"'INTENSITY/FLUX:' 'Flux'\n'ABFIND        :' '.false.'\n'MODELOPAC:' '{opac}'\n"
           f"'RESULTFILE :' '{res}'\n'METALLICITY:'    '0.00'\n"
           f"'INDIVIDUAL ABUNDANCES:'   '1'\n{cfg['Z']}  {a_x:.3f}\n'ISOTOPES : ' '0'\n"
           f"'NFILES   :' '1'\n{linelist}\n'SPHERICAL:'  'F'\n  30\n  300.00\n  15\n  1.30\n")
    r = subprocess.run([f"{EXE}/bsyn_lu"], input=ctl, capture_output=True, text=True, cwd=f"{W}/work")
    if nlte:
        so = r.stdout
        engaged = ("read departure coefficients" in so
                   and "missing level identification" not in so
                   and "continuing with departure coefficient = 1. for that level" not in so)
    else:
        engaged = True
    if not os.path.exists(res):
        raise SystemExit(f"bsyn failed {el} {tag}:\n{r.stdout[-1500:]}")
    return res, engaged


def ew(specfile, c, hw=EW_HW):
    d = np.loadtxt(specfile)
    w, fl = d[:, 0], d[:, 1]
    m = (w > c - hw) & (w < c + hw)
    return float(np.trapz(1.0 - fl[m], w[m]) * 1000.0)


def main():
    el = sys.argv[1]
    cfg = ELEMENTS[el]
    os.makedirs(f"{W}/work", exist_ok=True)
    # DATA/ symlink for babsma/bsyn
    dsl = f"{W}/work/DATA"
    if not os.path.lexists(dsl):
        os.symlink("/mnt/codex-data/engines/Turbospectrum_NLTE/DATA", dsl)

    ion = ION.get(el, 0)
    hdr, rows = ges_lines(cfg['Z'], ion, cfg['waves'])
    linelist = f"{W}/work/{el}.lines"
    m = re.match(r"^('\s*\d+\.\d\d\d\s*')\s+(\d+)\s+(\d+)", hdr)
    ion_roman = 'II' if ion else 'I'
    with open(linelist, "w") as f:
        f.write(f"{m.group(1)}    {m.group(2)}      {len(rows)}\n")
        f.write(f"'{el} {ion_roman}'\n")
        for r in rows:
            f.write(r + "\n")

    dep = make_departure(el, cfg)
    nlteinfo = f"{W}/work/species_{el}.dat"
    with open(nlteinfo, "w") as f:
        f.write(f"# RYA-534 {el} NLTE\n# path for model atom files\n{GT}/\n"
                f"# path for departure files\n{W}/Testout/\n# atomic NLTE setup\n"
                f"{cfg['Z']}  '{el}'  'nlte'  '{cfg['atom']}'  'sun_{el}_coef.dat'  'ascii'\n")

    lmin = min(cfg['waves']) - LMARGIN
    lmax = max(cfg['waves']) + LMARGIN
    model = f"{W}/Testout/sun_{el}.interpol"
    opac = f"{W}/work/sun_{el}.opac"
    print(f"=== {el}: babsma ==="); babsma(el, model, opac, lmin, lmax)
    print(f"=== {el}: bsyn NLTE @ A={cfg['a_sun']} ===")
    res_n, engaged = bsyn(el, cfg, cfg['a_sun'], True, "nlte", nlteinfo, linelist, opac, lmin, lmax)
    print("  NLTE departure engaged:", engaged)
    if not engaged:
        raise SystemExit("RAISE: NLTE departures not engaged (silent LTE) — STOP")
    ewn = {c: ew(res_n, c) for c in cfg['waves']}

    print(f"=== {el}: bsyn LTE COG ===")
    offs = [-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3]
    A = np.array([cfg['a_sun'] + o for o in offs])
    cog = {c: [] for c in cfg['waves']}
    for a in A:
        res_l, _ = bsyn(el, cfg, a, False, f"lte_{a:.2f}", nlteinfo, linelist, opac, lmin, lmax)
        for c in cfg['waves']:
            cog[c].append(ew(res_l, c))

    print(f"\n=== {el} RESULT (delta = A_sun - A*) ===")
    deltas = {}
    for c in cfg['waves']:
        ewl = np.array(cog[c])
        a_star = float(np.interp(ewn[c], ewl, A))
        deltas[c] = cfg['a_sun'] - a_star
        print(f"  {el} {c:.3f}: EW_NLTE={ewn[c]:.2f} mA  A*={a_star:.4f}  delta={deltas[c]:+.4f}")
    med = float(np.median(list(deltas.values())))
    passed = abs(med - cfg['anchor']) <= cfg['tol']
    print(f"  median delta = {med:+.4f}")
    print(f"  anchor {cfg['anchor']:+.3f} +/- {cfg['tol']} ({cfg['ref']})")
    print(f"  VERDICT: {'PASS' if passed else 'CHECK'}")


if __name__ == '__main__':
    main()
