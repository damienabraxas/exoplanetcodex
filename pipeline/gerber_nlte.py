"""
pipeline/gerber_nlte.py
=======================
RYA-798 — the adapter that lets the PRODUCTION flux-fit path use the TS-native Gerber
departure decks. Engine-B NLTE had a validated deck (RYA-785) and no route to it.

WHAT WAS MISSING, precisely
---------------------------
iSpec's Turbospectrum wrapper has always accepted NLTE:
`ispec.synth.turbospectrum.generate_spectrum(..., nlte_departure_coefficients=...)` writes
the departure + `nlteinfo` files and drives `bsyn_lu`. Two things kept us out of it:

  * `_synth_flux_at_abund` — the ONE generator both the v1 EW path and the v2 flux fit call
    — had no NLTE parameter, so there was nowhere to hand them in;
  * iSpec's own interpolator reads `input/dep-grid/{El}_nlte_grid_data.h5`, and that store
    is EMPTY on Sirius. Our decks are TS-native (`atom.*` + `auxData_*.dat` + the `.bin`),
    consumed only by `scripts/ts_gerber_gate.py`'s direct `interpol_modeles_nlte` + `bsyn`.

So this module produces, from the TS-native deck, exactly the tuple iSpec builds from its
own HDF5 grids:

    {element_name: (departures, tau, ndep, nk, Z, abundance, atom_model_path)}

⚠️ THE DECK HAS NO ABUNDANCE AXIS — MEASURED, NOT ASSUMED
----------------------------------------------------------
`interpol_modeles_nlte` takes an abundance argument, so it looks as though departures can
follow a trial abundance. They cannot. Running the interpolator at A = 7.36 / 7.46 / 7.56
gives departure files that are **byte-identical except for line 8, the abundance stamp
itself** (444836 bytes each; 1 differing line of 132). The Gerber aux table says why:
`A(X)` is not an axis at all, it is exactly `7.50 + [Fe/H]` across all 15229 nodes. The
grid's axes are Teff / logg / [Fe/H] / [alpha/Fe] / vturb.

Two consequences, both load-bearing:

  1. **The departures are a property of the ATMOSPHERE NODE, not of the trial abundance.**
     A flux fit varies A(X) to minimise chi2; the departure coefficients cannot track it.
     We therefore hold them FIXED across the fit and let only the stamp follow the trial
     value — which is what the deck does anyway, not an approximation we introduced. It is
     still an approximation *of the physics* (more Fe means more line opacity means a
     different radiation field means different departures), and it is second-order over the
     ~0.3 dex a fit explores, so it is STATED on the product rather than buried here.
  2. **The interpolation can be cached.** One interpolator run per (element, node) serves
     every trial abundance, instead of one per chi2 evaluation.

⚠️ SOLAR NODE ONLY, TODAY
--------------------------
`ts_gerber_gate` passes ONE concrete MARCS model eight times as the eight interpolation
corners — the TS idiom for a degenerate box, i.e. "use this model, do not interpolate".
That is correct for a solar gate and is why the departure file is unchanged by the
[Fe/H] argument. Serving another star needs real corner selection, so `for_node()` refuses
any node it cannot prove it has the model for, rather than silently returning solar
departures for Procyon. The eight corner paths are recorded in the departure file's own
footer, so the check is against the file rather than against our intent.
"""
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

INTERP = "/mnt/codex-data/engines/Turbospectrum_NLTE/interpolator/interpol_modeles_nlte"
GT = "/mnt/codex-data/grids/nlte/gerber_ts"
PROV_DIR = ROOT / "data" / "nlte_grids" / "gerber_ts"
MARCS_SOLAR = ("/mnt/codex-data/grids/model_atmospheres/marcs_standard_comp/"
               "marcs_standard_comp/p5750_g+4.5_m0.0_t01_st_z+0.00_a+0.00_c+0.00_n+0.00_"
               "o+0.00_r+0.00_s+0.00.mod")

# element -> (Z, atom, aux, grid). Mirrors ts_gerber_gate.ELEMENTS but carries only what
# the adapter needs; the GATE remains the authority on whether a deck may be used at all.
DECKS = {
    "Fe": dict(Z=26, atom="atom.fe607a",
               aux="auxData_Fe_MARCS_May-07-2021.dat",
               grid="NLTEgrid4TS_Fe_MARCS_May-07-2021.bin"),
}


class GerberDeckError(RuntimeError):
    """Raised loudly. There is no degraded mode: silent LTE is the failure to avoid."""


def deck_abundance(element: str) -> float:
    """The abundance the deck's departures were computed at, from its provenance record.

    NEVER a fitted or reference value. Fe's grid was computed at A(Fe) = 7.46, not the 7.50
    printed in the `atom.fe607a` header, and bsyn STOPs on the mismatch.
    """
    p = PROV_DIR / f"{element}_gerber2023.prov.json"
    if not p.exists():
        raise GerberDeckError(
            f"no provenance record at {p} — an unregistered deck has not passed the "
            f"RYA-534/785 gate and must not be used as an Engine-B leg")
    d = json.loads(p.read_text())
    try:
        return float(d["deck_abundance"]["a_sun"])
    except (KeyError, TypeError, ValueError) as e:
        raise GerberDeckError(f"{p} carries no usable deck_abundance.a_sun: {e}") from e


def read_departure_file(path: str | os.PathLike) -> dict:
    """Parse `interpol_modeles_nlte`'s output.

    Layout, established by reading a real file rather than from documentation:

        0-7                 eight `#` comment lines (the interpolation weights)
        8                   abundance A(X) the file is STAMPED with
        9                   ndep   (depth points)
        10                  nk     (model-atom levels; 607 for atom.fe607a)
        11 .. 10+ndep       ndep log-tau values
        next ndep lines     nk departure coefficients each  -> (ndep, nk)
        final 8 lines       the eight MARCS corner model paths (provenance footer)
    """
    lines = Path(path).read_text(errors="replace").splitlines()
    if len(lines) < 12:
        raise GerberDeckError(f"{path}: too short to be a departure file")
    try:
        abundance = float(lines[8])
        ndep = int(lines[9])
        nk = int(lines[10])
    except ValueError as e:
        raise GerberDeckError(f"{path}: unreadable header ({e})") from e

    tau = np.array([float(x) for x in lines[11:11 + ndep]], dtype=float)
    if tau.size != ndep:
        raise GerberDeckError(f"{path}: expected {ndep} tau values, got {tau.size}")

    body = lines[11 + ndep:]
    rows = [r.split() for r in body if len(r.split()) == nk]
    if len(rows) != ndep:
        raise GerberDeckError(
            f"{path}: expected {ndep} departure rows of {nk} values, found {len(rows)}")
    dep = np.array(rows, dtype=float)

    corners = [r.strip() for r in body[len(rows):] if r.strip()]
    return dict(abundance=abundance, ndep=ndep, nk=nk, tau=tau,
                departures=dep, corners=corners)


def _interpolate(element: str, node: str, teff: float, logg: float, feh: float,
                 abundance: float, workdir: str) -> str:
    cfg = DECKS[element]
    binf, aux = f"{GT}/{cfg['grid']}", f"{GT}/{cfg['aux']}"
    for p in (INTERP, binf, aux):
        if not os.path.exists(p):
            raise GerberDeckError(f"missing Gerber asset (Sirius-only): {p}")
    with open(aux) as fh:
        nrows = sum(1 for ln in fh if not ln.lstrip().startswith("#"))

    os.makedirs(f"{workdir}/work", exist_ok=True)
    os.makedirs(f"{workdir}/Testout", exist_ok=True)
    dep = f"{workdir}/Testout/{node}_{element}_coef.dat"
    if os.path.exists(dep):
        os.remove(dep)          # never inherit a previous run's file (RYA-785 stale guard)

    stdin = "\n".join([f"'{MARCS_SOLAR}'"] * 8 + [
        f"'{workdir}/Testout/{node}_{element}.interpol'",
        f"'{workdir}/Testout/{node}_{element}.alt'", f"'{dep}'",
        f"'{binf}'", f"'{aux}'", str(nrows),
        f"{teff:.0f}", f"{logg:+.2f}", f"{feh:+.2f}", f"{abundance:.2f}",
        ".false.", ".false.", "'none'", ""])
    r = subprocess.run([INTERP], input=stdin, capture_output=True, text=True,
                       cwd=f"{workdir}/work")
    if not os.path.exists(dep):
        raise GerberDeckError(
            f"interpol_modeles_nlte produced no departure file for {element}:\n"
            f"{r.stdout[-1200:]}\n{r.stderr[-400:]}")
    return dep


_CACHE: dict[tuple, dict] = {}


def for_node(element: str, teff: float, logg: float, feh: float,
             node: str = "solar",
             workdir: str = "/tmp/rya798_gerber") -> dict:
    """Interpolated departures for one atmosphere node. Cached — see the module docstring:
    the departures do not depend on the abundance, so one run serves every trial."""
    if element not in DECKS:
        raise GerberDeckError(
            f"no TS-native Gerber deck registered for {element!r} in this adapter "
            f"(registered: {sorted(DECKS)}). Staging a deck is RYA-710; validating it is "
            f"the RYA-534/785 gate. Both must happen before it can be used here.")
    key = (element, round(teff, 1), round(logg, 3), round(feh, 3))
    if key in _CACHE:
        return _CACHE[key]

    # SOLAR NODE ONLY — see the docstring. Refuse rather than quietly hand back solar
    # departures for another star.
    if not (abs(teff - 5750.0) <= 60.0 and abs(logg - 4.50) <= 0.10
            and abs(feh - 0.0) <= 0.10):
        raise GerberDeckError(
            f"the Gerber departure path is solar-node-only today: it passes ONE MARCS "
            f"model ({os.path.basename(MARCS_SOLAR)}) eight times as the eight "
            f"interpolation corners, which is a degenerate box. Asked for "
            f"teff={teff:.0f} logg={logg:.2f} feh={feh:+.2f}. Serving another star needs "
            f"real corner selection — do that rather than accepting solar departures.")

    a_deck = deck_abundance(element)
    dep_path = _interpolate(element, node, teff, logg, feh, a_deck, workdir)
    parsed = read_departure_file(dep_path)

    # Verify against the file's OWN footer, not against our intent.
    if parsed["corners"] and not all(os.path.basename(MARCS_SOLAR) in c
                                     for c in parsed["corners"]):
        raise GerberDeckError(
            f"departure file corners are not the requested solar model: "
            f"{parsed['corners'][:2]}")
    parsed["atom_path"] = f"{GT}/{DECKS[element]['atom']}"
    parsed["Z"] = DECKS[element]["Z"]
    parsed["deck_abundance"] = a_deck
    _CACHE[key] = parsed
    return parsed


def as_ispec_tuple(parsed: dict, abundance: float) -> tuple:
    """Pack a parsed departure file into iSpec's tuple.

    ⚠️ ORIENTATION. `interpol_modeles_nlte` writes the block as (ndep, nk) — one line per
    depth, nk values along it — and iSpec's writer wants **(nk, ndep)**. It checks, and
    raises `Shape of depart_coeffs is (56, 607), but expected (607, 56)`, so this is a
    caught error rather than a silent one; the transpose lives here so exactly one place
    knows about the difference.
    """
    return (parsed["departures"].T, parsed["tau"], parsed["ndep"], parsed["nk"],
            parsed["Z"], float(abundance), parsed["atom_path"])


def coefficients(element: str, teff: float, logg: float, feh: float,
                 abundance: float, **kw) -> dict:
    """The dict `ispec.generate_spectrum(nlte_departure_coefficients=...)` expects.

    `abundance` is the value the SYNTHESIS will run at. bsyn checks it against the stamp in
    the departure file and STOPs on a mismatch, so the stamp must follow the trial value.
    The coefficients themselves are node-fixed (module docstring) — this is the stamp
    following the synthesis, not the physics following the abundance.
    """
    return {element: as_ispec_tuple(for_node(element, teff, logg, feh, **kw),
                                    abundance)}


def assert_depth_match(parsed: dict, atmosphere) -> None:
    """iSpec OVERWRITES the departure tau with the atmosphere's own (`atmosphere_layers[:, 7]`)
    to dodge Turbospectrum's "tau scales differ" error. That silently assumes the two have
    the same number of layers — if they do not, depth i of the departures is paired with a
    different physical depth of the atmosphere and the correction is quietly wrong.
    """
    n = int(len(atmosphere))
    if n != parsed["ndep"]:
        raise GerberDeckError(
            f"depth mismatch: the departure grid has {parsed['ndep']} layers, the model "
            f"atmosphere has {n}. iSpec overwrites the departure tau with the "
            f"atmosphere's, so these MUST match or the departures are applied at the "
            f"wrong depths.")


def assert_linelist_supports_nlte(linelist, Z: int, element: str,
                                  wave_lo_A: float | None = None,
                                  wave_hi_A: float | None = None) -> int:
    """iSpec fails SOFT — replicate its skip test and RAISE instead.

    `generate_spectrum` appends an element to `nlte_ignored` and continues **in LTE without
    raising** when no line of that element carries NLTE labels. That is the RYA-764 defect
    (2,644 in-band Fe I lines at `nlte_label_up='none'` running silently in LTE), and it is
    the single worst outcome available here: an LTE synthesis shipped under an NLTE label.
    iSpec does not report `nlte_available` back to the caller, so the check has to happen
    on this side of the call.
    """
    try:
        species = np.floor(np.asarray(linelist["turbospectrum_species"], dtype=float))
        sel = species == Z
        # ⚠️ WINDOW-LOCAL, NOT ELEMENT-GLOBAL. bsyn applies departures per LINE and falls
        # back to departure = 1 for any line whose levels are unidentified, so "Fe is
        # labelled somewhere in the list" says nothing about the lines actually being
        # synthesised here. Fe answers the global question yes 15706 times while every
        # line in the Fe II 6910-9199 window is unlabelled.
        if wave_lo_A is not None and wave_hi_A is not None:
            names = linelist.dtype.names or ()
            if "wave_A" in names:
                w = np.asarray(linelist["wave_A"], dtype=float)
            elif "wave_nm" in names:
                w = np.asarray(linelist["wave_nm"], dtype=float) * 10.0
            else:
                w = None
            if w is not None:
                sel = sel & (w >= wave_lo_A) & (w <= wave_hi_A)
        rows = linelist[sel]
        if len(rows) == 0:
            raise GerberDeckError(
                f"no {element} lines in this window's linelist — NLTE cannot engage")
        flagged = rows[rows["nlte"] == "T"]
        n = int(np.sum((flagged["nlte_label_low"] != "none")
                       | (flagged["nlte_label_up"] != "none")))
    except GerberDeckError:
        raise
    except Exception as e:
        raise GerberDeckError(f"cannot read NLTE labels from the linelist: {e}") from e
    if n == 0:
        where = (f" in {wave_lo_A:.1f}-{wave_hi_A:.1f} A"
                 if wave_lo_A is not None and wave_hi_A is not None else "")
        raise GerberDeckError(
            f"{element}: {len(rows)} lines{where} but NONE carry NLTE level labels. bsyn "
            f"sets departure = 1 for an unidentified level and iSpec drops an unlabelled "
            f"element into `nlte_ignored` — either way the synthesis runs in LTE WITHOUT "
            f"RAISING (RYA-534 Co/Ni, RYA-764). Refusing: an LTE spectrum under an NLTE "
            f"label is worse than no product. This is a LINE-LIST coverage gap, not a "
            f"deck failure — the deck is fine, these transitions have no level "
            f"identification.")
    return n
