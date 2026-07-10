#!/usr/bin/env python3
"""
RYA-542 — independent TSFitPy MARCS Bergemann-2011 Ti I solar NLTE synthesis.

DECISIVE discriminator deferred from RYA-535. Engine-B (our TS-Gerber deck,
scripts/ts_gerber_gate.py; MARCS + atom.ti503b) gives a solar Ti I NLTE
correction of +0.221; Engine-A (MPIA Bergemann-2011, MAFAGS-OS) gives +0.108.
Both use the SAME Bergemann-2011 atom (Gerber 2023 Table 1), so the +0.11 excess
is NOT a model-atom difference. It is either (a) a REAL MAFAGS-OS-vs-MARCS
atmosphere effect, or (b) a bug in OUR deck's Ti line handling.

This driver is an INDEPENDENT reading of the SAME physics:
  * MARCS standard-composition atmosphere (marcs_standard_comp) — same as the deck
  * atom.ti503b (Bergemann-2011) + the Gerber NLTEgrid4TS Ti departure grid
  * the VERBATIM GES Ti I line rows (loggf, chi, level IDs) — via the deck's ges_lines()
  * the compiled Turbospectrum_NLTE bsyn_lu/babsma_lu + interpol_modeles_nlte
but driven by TSFitPy's OWN Python layer (TurboSpectrum.make_babsma_bsyn_file,
departure-interpolation config, calculate_equivalent_width, root_scalar abundance
inversion) — i.e. TSFitPy's `generate_and_fit_atmosphere`. It shares ONLY the
compiled Fortran with the deck; the suspected "deck line handling" lives entirely in
the Python driver, which here is TSFitPy's, not ours.

  ~+0.22  -> MARCS genuinely gives the larger correction (driver-independent) =>
            REAL atmosphere effect; Ti's Engine-B value is legitimate, atmosphere-flagged.
  ~+0.11  -> agrees with MAFAGS-OS; the deck's +0.22 is a deck line-handling bug; Ti CHECK.

delta = A_NLTE - A_LTE per line (same convention as the deck's +0.221 and the ticket).
Solar node from STAR_PARAMS via the deck's _solar_node() (NO hardcode, RYA-298).
Validate-don't-tune: no fitting to a target. RUNS ON SIRIUS ONLY (grids Sirius-only, RYA-526).
"""
from __future__ import annotations
import argparse
import os
import pickle
import re
import sys

import numpy as np

# --- Sirius paths (grids/engines Sirius-only, RYA-526) -----------------------
TSFITPY = "/srv/codex/engines/TSFitPy"
DECK    = "/mnt/codex-data/codex/rya534"          # ts_gerber_gate: _solar_node + ges_lines (STAR_PARAMS SSOT)
RUN     = "/mnt/codex-data/codex/rya542/tsfitpy_run"
GT      = "/mnt/codex-data/grids/nlte/gerber_ts"  # atom.ti503b + auxData + the pulled NLTEgrid .bin
MARCS   = "/mnt/codex-data/grids/model_atmospheres/marcs_standard_comp/marcs_standard_comp/"
MARCS_LIST_SRC = f"{TSFITPY}/input_files/model_atmospheres/1D/model_atmosphere_list.txt"
MARCS_LIST     = f"{RUN}/marcs1d/model_atmosphere_list.txt"
EXE     = "/mnt/codex-data/engines/Turbospectrum_NLTE/exec-gf/"
INTERP  = "/mnt/codex-data/engines/Turbospectrum_NLTE/interpolator/"

BIN  = "NLTEgrid4TS_TI_MARCS_Feb-21-2022.bin"     # prov Ti_gerber2023.prov.json (RYA-540 cache)
AUX  = "auxData_TI_MARCS_Feb-21-2022.dat"
ATOM = "atom.ti503b"

# Ti diagnostic lines (deck gate; stronger, GES-identified, MPIA-grid overlap).
LINES = [5689.460, 5648.565, 5662.150]
Z_TI, ION_TI = 22, 0
LDELTA = 0.005
# Deck a_sun(Ti)=4.90 (Asplund 2021); TSFitPy solar A(Ti)=4.94 -> [Ti/Fe] offset to match
# the deck's line-strength regime (the delta itself is a differential, ~offset-insensitive).
DECK_ATI = 4.90
TSFITPY_SOLAR_ATI = 4.94

# Reference numbers to compare against (the whole point of this ticket).
DECK_DELTA   = {5689.460: +0.2658, 5648.565: +0.2073, 5662.150: +0.2206}  # Engine-B TS-Gerber (MARCS)
DECK_MEDIAN  = +0.2206
ENGINE_A_MAFAGS = +0.108                                                   # Engine-A MPIA Bergemann-2011 (MAFAGS-OS)


def _import_deck():
    if DECK not in sys.path:
        sys.path.insert(0, DECK)
    from ts_gerber_gate import _solar_node, ges_lines  # noqa: E402
    return _solar_node, ges_lines


def _import_tsfitpy():
    if TSFITPY not in sys.path:
        sys.path.insert(0, TSFITPY)
    from scripts.calculate_nlte_correction_line import generate_and_fit_atmosphere, AbusingClasses  # noqa
    from scripts.turbospectrum_class_nlte import TurboSpectrum  # noqa
    from scripts.synthetic_code_class import fetch_marcs_grid  # noqa
    return generate_and_fit_atmosphere, AbusingClasses, TurboSpectrum, fetch_marcs_grid


def _header_one(hdr_species: str) -> str:
    """Turn the GES block header ('  22.000  '  1  6803) into a 1-line-count header."""
    return re.sub(r"(\d+)\s*$", "1", hdr_species.rstrip("\n"))


def setup():
    """Stage the per-line minimal linelists (one verbatim GES Ti I line each, level IDs
    preserved) + the MARCS model list. Idempotent."""
    _, ges_lines = _import_deck()
    os.makedirs(f"{RUN}/marcs1d", exist_ok=True)
    os.makedirs(f"{RUN}/temp", exist_ok=True)
    # MARCS list: copy TSFitPy's shipped 1D standard-comp list (no p/.mod; fetch_marcs_grid
    # reconstructs pXXXX...mod, which exist in marcs_standard_comp).
    with open(MARCS_LIST_SRC) as f:
        lst = f.read()
    with open(MARCS_LIST, "w") as f:
        f.write(lst)
    hdr_species, rows = ges_lines(Z_TI, ION_TI, LINES)   # verbatim, SAME as the deck
    hdr1 = _header_one(hdr_species)
    row_by_wl = dict(zip(LINES, rows))
    for wl in LINES:
        d = f"{RUN}/ll_{wl:.3f}"
        os.makedirs(d, exist_ok=True)
        with open(f"{d}/ti_line.bsyn", "w") as f:
            f.write(hdr1 + "\n")
            f.write("'Ti I     NLTE'\n")
            f.write(row_by_wl[wl].rstrip("\n") + "\n")
    print(f"[setup] MARCS list -> {MARCS_LIST} ({len(lst.splitlines())} models)")
    print(f"[setup] linelists  -> {RUN}/ll_<wl>/ti_line.bsyn (verbatim GES rows):")
    for wl in LINES:
        print(f"           {wl:.3f}: {row_by_wl[wl].strip()[:80]}...")


def _build_abusingclasses(need_bin: bool):
    _, AbusingClasses, TurboSpectrum, fetch_marcs_grid = _import_tsfitpy()
    ac = AbusingClasses()
    ac.atmosphere_type = "1D"
    ac.global_temp_dir = f"{RUN}/temp/"
    ac.segment_file = None
    ac.spectral_code_path = EXE
    ac.interpol_path = INTERP
    ac.ldelta = LDELTA
    ac.output_folder = f"{RUN}/output/"
    ac.linemask_file = None
    # atmosphere (same grid for the NLTE and LTE=_1d branches; NLTE just adds departures)
    ac.model_atmosphere_grid_path = MARCS
    ac.model_atmosphere_grid_path_1d = MARCS
    ac.model_atmosphere_list = MARCS_LIST
    ac.model_atmosphere_list_1d = MARCS_LIST
    # NLTE data (Gerber grid via the RYA-540 persistent cache). NOTE: TSFitPy bsyn builds the
    # model-atom path by RAW string concat (model_atom_path + atom_file), so it MUST end in "/".
    ac.model_atom_path = GT + "/"
    ac.departure_file_path = GT
    ac.depart_bin_file_dict = {"Ti": BIN}
    ac.depart_aux_file_dict = {"Ti": AUX}
    ac.model_atom_file_dict = {"Ti": ATOM}
    if need_bin and not os.path.isfile(os.path.join(GT, BIN)):
        raise SystemExit(f"[run] Ti NLTE grid not present yet: {os.path.join(GT, BIN)} "
                         f"(RYA-540 pull still in flight?).")
    ac.aux_file_length_dict = {"Ti": int(len(np.loadtxt(os.path.join(GT, AUX), dtype="str")))}
    mt, ml, mm, mvk, mmod, mv = fetch_marcs_grid(MARCS_LIST, TurboSpectrum.marcs_parameters_to_ignore)
    ac.model_temperatures = mt; ac.model_logs = ml; ac.model_mets = mm
    ac.marcs_value_keys = mvk; ac.marcs_models = mmod; ac.marcs_values = mv
    ac.model_temperatures_1d = mt; ac.model_logs_1d = ml; ac.model_mets_1d = mm
    ac.marcs_value_keys_1d = mvk; ac.marcs_models_1d = mmod; ac.marcs_values_1d = mv
    os.makedirs(ac.global_temp_dir, exist_ok=True)
    return ac


def smoke_lte():
    """Validate the whole non-NLTE wiring (MARCS interpolation + babsma + bsyn + linelist)
    WITHOUT the NLTE grid. If EW_lte at A(Ti)=4.90 is sane (deck EWs 8.6/9.7/18.6 mA), the
    machinery is correct and only the NLTE grid read remains."""
    _import_tsfitpy()
    from scripts.calculate_nlte_correction_line import generate_atmosphere
    from scripts.auxiliary_functions import calculate_equivalent_width
    solar_node, _ = _import_deck()
    node = solar_node()
    teff, logg, met, xi = node["teff"], node["logg"], node["feh"], node["xi"]
    off = DECK_ATI - TSFITPY_SOLAR_ATI
    ac = _build_abusingclasses(need_bin=False)
    print(f"[smoke-lte] solar node teff={teff} logg={logg} feh={met} xi={xi}  [Ti/Fe]={off:+.3f} -> A(Ti)={DECK_ATI}")
    for wl in LINES:
        lldir = f"{RUN}/ll_{wl:.3f}/"
        w, fl = generate_atmosphere(ac, teff, logg, xi, met, wl - 7, wl + 7, LDELTA,
                                    lldir, "Ti", off, {}, False, verbose=False)
        if w is None:
            print(f"[smoke-lte] {wl:.3f}: LTE synthesis FAILED"); continue
        ew = calculate_equivalent_width(w, fl, wl - 5, wl + 5) * 1000.0
        print(f"[smoke-lte] {wl:.3f}: EW_LTE = {ew:.2f} mA   (deck EW_NLTE was 8.6/9.7/18.6)")


def debug_nlte():
    """Run ONE line's NLTE synthesis with verbose Fortran output + preserve the temp dir,
    to surface interpolator/bsyn errors (why EW_nlte came back as the 9999999 sentinel)."""
    import shutil
    _import_tsfitpy()
    from scripts.calculate_nlte_correction_line import generate_atmosphere
    solar_node, _ = _import_deck()
    node = solar_node()
    teff, logg, met, xi = node["teff"], node["logg"], node["feh"], node["xi"]
    off = DECK_ATI - TSFITPY_SOLAR_ATI
    ac = _build_abusingclasses(need_bin=True)
    # keep the temp dir so we can inspect the interpol config + coef.dat
    _orig_rmtree = shutil.rmtree
    shutil.rmtree = lambda *a, **k: None
    wl = LINES[0]
    print(f"[debug-nlte] {wl}: NLTE synth verbose (aux_len={ac.aux_file_length_dict['Ti']})")
    w, fl = generate_atmosphere(ac, teff, logg, xi, met, wl - 7, wl + 7, LDELTA,
                                f"{RUN}/ll_{wl:.3f}/", "Ti", off, {}, True, verbose=True)
    shutil.rmtree = _orig_rmtree
    print(f"[debug-nlte] result: wave is {'None' if w is None else 'OK (%d pts)' % len(w)}")
    print(f"[debug-nlte] temp dirs under {ac.global_temp_dir}:")
    import glob as _g
    for d in sorted(_g.glob(ac.global_temp_dir + "*")):
        print("   ", d)


def run():
    """Full independent NLTE correction per line: delta = A_NLTE - A_LTE."""
    gen_fit, _, _, _ = _import_tsfitpy()
    solar_node, _ = _import_deck()
    node = solar_node()
    teff, logg, met, xi = node["teff"], node["logg"], node["feh"], node["xi"]
    off = DECK_ATI - TSFITPY_SOLAR_ATI
    ac = _build_abusingclasses(need_bin=True)
    pkl = f"{RUN}/abusingclasses.pkl"
    with open(pkl, "wb") as f:
        pickle.dump(ac, f)
    print(f"=== RYA-542 independent TSFitPy MARCS Bergemann-2011 Ti I NLTE correction ===")
    print(f"  solar node (STAR_PARAMS, no hardcode): teff={teff} logg={logg} feh={met} xi={xi}")
    print(f"  A(Ti)={DECK_ATI} ([Ti/Fe]={off:+.3f}); grid={BIN}; atom={ATOM}; atmosphere=MARCS std-comp")
    deltas = {}
    for wl in LINES:
        lldir = f"{RUN}/ll_{wl:.3f}/"
        res = gen_fit(pkl, "solar", teff, logg, xi, met, wl - 2, wl + 2, LDELTA,
                      lldir, "Ti", off, {}, wl, verbose=False)
        # res = ["name\tteff\tlogg\tmet\tvmic\twave\tew_lte\tew_nlte\tnlte_abund\tnlte_correction"]
        parts = res[0].split("\t")
        ew_lte, ew_nlte, nlte_abund, corr = (float(parts[6]), float(parts[7]),
                                             float(parts[8]), float(parts[9]))
        deltas[wl] = corr
        print(f"  Ti {wl:.3f}: EW_LTE={ew_lte:.2f} EW_NLTE={ew_nlte:.2f} mA  "
              f"A_NLTE-A_LTE(delta)={corr:+.4f}   [deck {DECK_DELTA[wl]:+.4f}]")
    med = float(np.median(list(deltas.values())))
    print(f"\n  TSFitPy-MARCS median delta = {med:+.4f}")
    print(f"  deck (Engine-B TS-Gerber, MARCS) median = {DECK_MEDIAN:+.4f}")
    print(f"  Engine-A (MPIA Bergemann-2011, MAFAGS-OS) = {ENGINE_A_MAFAGS:+.4f}")
    d_atm = abs(med - DECK_MEDIAN)
    d_mafags = abs(med - ENGINE_A_MAFAGS)
    branch = ("ATMOSPHERE (~+0.22): TSFitPy(MARCS) reproduces the deck -> the larger Ti NLTE "
              "correction is REAL (MAFAGS-OS-vs-MARCS atmosphere), Ti value legit + atmosphere-flagged"
              if d_atm < d_mafags else
              "DECK (~+0.11): TSFitPy(MARCS) agrees with MAFAGS-OS, NOT the deck -> the deck's +0.22 "
              "is a Ti line-handling bug, Ti stays CHECK")
    print(f"\n  |median - deck +0.221| = {d_atm:.3f}   |median - MAFAGS +0.108| = {d_mafags:.3f}")
    print(f"  DECISION BRANCH: {branch}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="RYA-542 independent TSFitPy Ti NLTE correction")
    ap.add_argument("--setup", action="store_true", help="stage linelists + MARCS list")
    ap.add_argument("--smoke-lte", action="store_true", help="validate non-NLTE wiring (no grid needed)")
    ap.add_argument("--run", action="store_true", help="full NLTE correction per line")
    ap.add_argument("--debug-nlte", action="store_true", help="one-line verbose NLTE synth (debug)")
    a = ap.parse_args()
    if a.setup:
        setup()
    if a.smoke_lte:
        smoke_lte()
    if a.debug_nlte:
        debug_nlte()
    if a.run:
        run()
    if not (a.setup or a.smoke_lte or a.run or a.debug_nlte):
        ap.print_help()
