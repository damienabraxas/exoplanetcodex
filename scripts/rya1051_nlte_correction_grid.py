#!/usr/bin/env python3
"""RYA-1051: compute OUR Fe I NLTE abundance corrections the way Amarsi's grid defines
them, across metallicity — so the two are the same quantity and can be differenced.

🔴 WHY THIS REPLACES THE PROXY. The first pass estimated the correction as
`-log10(b_lower)` sampled at one optical depth. That proxy has two fatal limits:

  * it is STRICTLY POSITIVE whenever b < 1, so it can never produce a negative
    correction -- yet Amarsi's solar value IS negative (-0.029). The proxy could
    explain our sign and was structurally incapable of explaining his.
  * its answer depended on the tau I sampled: the gradient ratio ran 1.96x to 7.73x
    and the solar difference +0.048 to +0.126 dex across log tau -0.5 .. -2.0. That
    is a chosen cut doing the work (RYA-161), not a measurement.

This computes the ACTUAL correction instead, by its definition:

    synthesise the line in LTE at A0
    find the abundance A at which the NLTE synthesis has the SAME equivalent width
    Delta_NLTE = A - A0

which is what an observer applies to an LTE abundance, and what Amarsi's grid tabulates.
Both the opacity term (b_lower) and the source-function term (b_upper/b_lower) are then
present automatically, because the radiative transfer is actually solved.

⚠️ vturb IS MATCHED TO THE ANCHOR, NOT TO THE SUN. Amarsi's grid has no vturb = 1.0 node
and the values extracted for comparison are at 1.50, so the syntheses run at 1.50. A
solar-vturb run would be a different quantity from the one it is being compared against.

⚠️ EW, NOT FLUX. The correction is defined on line STRENGTH, so it is measured as an
equivalent width and not as a flux-fit chi2 -- a flux fit would fold in the profile shape,
which is not what the anchor tabulates.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

TEFF, LOGG, VTURB = 5750.0, 4.5, 1.50
HALF_A = 0.62          # the band run's half-width (RYA-967), so the window is the same
DA = 0.15              # trial offset for the local dEW/dA slope


def ew_mA(wave_nm: np.ndarray, flux: np.ndarray) -> float:
    """Equivalent width in mA from a normalised flux array."""
    # 1 nm = 10 A = 1e4 mA. (Had this as 1e7 first: the reported EWs came out 1000x too
    # large -- 232653 mA for a single line -- while Delta was UNAFFECTED, because it is a
    # ratio (W_lte - W_nlte)/slope and any constant scale cancels exactly. A wrong number
    # that cannot change the answer is still a wrong number in the artifact.)
    return float(np.trapezoid(1.0 - flux, wave_nm) * 1e4)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lines", required=True,
                    help="CSV with wavelength_air_A (our measured pool)")
    ap.add_argument("--n-lines", type=int, default=12)
    ap.add_argument("--feh", type=float, nargs="+",
                    default=[0.0, -1.0, -2.0, -3.0])
    ap.add_argument("--out", default=str(ROOT / "data" / "results" / "rya1051" /
                                         "nlte_correction_grid.json"))
    a = ap.parse_args()

    import pandas as pd
    from pipeline import gerber_nlte as gnlte
    import ispec
    from pipeline.abundances_derive import (_load_atmosphere, _synth_flux_at_abund,
                                            _load_synth_resources,
                                            _ISPEC_SOLAR_ABUND_FILE)

    d = pd.read_csv(a.lines)
    # spread over excitation potential -- the axis over-ionisation acts on
    d = d.sort_values("ep_eV").reset_index(drop=True)
    idx = np.linspace(0, len(d) - 1, a.n_lines).astype(int)
    waves = d.wavelength_air_A.values[idx].astype(float)
    eps = d.ep_eV.values[idx].astype(float)
    print(f"{len(waves)} lines, EP {eps.min():.2f}-{eps.max():.2f} eV")

    linelist, isotopes, _chem = _load_synth_resources()
    solar_abund = ispec.read_solar_abundances(_ISPEC_SOLAR_ABUND_FILE)
    doc = {"ticket": "RYA-1051",
           "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
           "teff": TEFF, "logg": LOGG, "vturb": VTURB,
           "definition": "Delta = A(NLTE matching the LTE equivalent width) - A(LTE)",
           "vturb_note": "1.50 to match the Amarsi grid node used for comparison; the "
                         "grid has no vturb=1.0 node and the Sun is 1.0.",
           "nodes": {}}

    for feh in a.feh:
        A0 = 7.50 + feh
        print(f"\n[Fe/H] = {feh:+.2f}   A0 = {A0:.2f}")
        atm = _load_atmosphere(TEFF, LOGG, feh, VTURB, model_grid="MARCS.GES")
        # 🔴 read_deck_node, NOT for_node. `for_node` routes through the vendor
        # interpolator, which is SOLAR-NODE-ONLY and says so:
        #     GerberDeckError: the Gerber departure path is solar-node-only today: it
        #     passes ONE MARCS model eight times as the eight interpolation corners,
        #     which is a degenerate box. Asked for teff=5750 logg=4.50 feh=-1.00.
        # That guard is right and it stopped this grid at its second node. The direct
        # binary read has no such limit -- and it is the BETTER instrument here anyway:
        # the Gerber grid carries an EXACT node at each of [Fe/H] = 0, -1, -2, -3 with
        # A(Fe) = 7.50 + [Fe/H], so nothing is interpolated at all.
        #
        # ⚠️ read_deck_node returns only (departures, tau, ndep, nk, abundance, corners),
        # so the three fields `as_ispec_tuple` needs are supplied here from the deck
        # registry. This is the same gap RYA-1040 found in the <3D> direct-read branch,
        # where a missing 'Z' meant no <3D> deck could reach a synthesis at all.
        dep = gnlte.read_deck_node("Fe", TEFF, LOGG, feh, A0)
        _d = gnlte.DECKS["Fe"]
        dep["Z"] = _d["Z"]
        dep["atom_path"] = f"{gnlte.GT}/{_d['atom']}"
        dep["deck_abundance"] = float(dep.get("abundance", A0))
        gnlte.assert_depth_match(dep, atm)
        # 🔴 THE DECK TUPLE IS REBUILT PER TRIAL ABUNDANCE, NOT ONCE.
        # `as_ispec_tuple(dep, A)` stamps A into the tuple as "the abundance these
        # departures belong to", and bsyn REFUSES if that disagrees with the abundance it
        # is synthesising at:
        #     Bsyn: NLTE departure coeff calculated for abundance = 7.50 while it is 7.65
        # Building it once at A0 and reusing it for the A0+DA call is exactly that
        # mismatch, and it cost a smoke test. The production fit rebuilds it every
        # evaluation (`as_ispec_tuple(dep, float(a_x))`) for this reason.
        # ⚠️ The DEPARTURES do not change with A: Fe's deck has a single A(X) node, so
        # this is a relabelling, and RYA-1035 established that the deck's own abundance is
        # 7.50. Keeping bsyn's STOP is correct -- it is the check that caught my error.
        def _nd(A):
            return {"Fe": gnlte.as_ispec_tuple(dep, float(A))}
        print(f"  atmosphere {len(atm)} layers, deck ndep={dep['ndep']} nk={dep['nk']}")

        rows = []
        for w, ep in zip(waves, eps):
            wn = np.arange((w - HALF_A) / 10.0, (w + HALF_A) / 10.0, 0.00002)
            kw = dict(waveobs_nm=wn, atmosphere=atm, teff=TEFF, logg=LOGG, feh=feh,
                      vturb=VTURB, linelist=linelist, isotopes=isotopes,
                      solar_abund=solar_abund, element="Fe", atom_code=26)
            try:
                f_lte = _synth_flux_at_abund(**kw, trial_A=A0)
                f_n0 = _synth_flux_at_abund(**kw, trial_A=A0, nlte_departures=_nd(A0))
                f_n1 = _synth_flux_at_abund(**kw, trial_A=A0 + DA,
                                            nlte_departures=_nd(A0 + DA))
            except Exception as e:
                print(f"    {w:.3f}  FAILED {type(e).__name__}: {str(e)[:60]}")
                continue
            W_lte, W_n0, W_n1 = (ew_mA(wn, f) for f in (f_lte, f_n0, f_n1))
            slope = (W_n1 - W_n0) / DA          # mA per dex, local
            if abs(slope) < 1e-6 or W_lte <= 0 or W_n0 <= 0:
                print(f"    {w:.3f}  degenerate (W_lte={W_lte:.2f} slope={slope:.3g})")
                continue
            delta = (W_lte - W_n0) / slope      # linear solve for the matching abundance
            rows.append({"wave_A": float(w), "ep_eV": float(ep),
                         "ew_lte_mA": round(W_lte, 3), "ew_nlte_mA": round(W_n0, 3),
                         "dEW_dA": round(slope, 3), "delta_nlte": round(float(delta), 4)})
            print(f"    {w:.3f}  ep={ep:4.2f}  W_LTE={W_lte:7.2f}  W_NLTE={W_n0:7.2f}  "
                  f"Delta={delta:+.4f}")
        if rows:
            dv = np.array([r["delta_nlte"] for r in rows])
            doc["nodes"][f"{feh:+.2f}"] = {
                "A0": A0, "n": len(rows),
                "median": round(float(np.median(dv)), 4),
                "mean": round(float(np.mean(dv)), 4),
                "min": round(float(dv.min()), 4), "max": round(float(dv.max()), 4),
                "lines": rows}
            print(f"  => median Delta_NLTE = {np.median(dv):+.4f}  (n={len(rows)})")

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.out).write_text(json.dumps(doc, indent=2) + "\n")
    print(f"\nwrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
