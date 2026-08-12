#!/usr/bin/env python3
"""RYA-798 CONTROL — does the wired path reproduce the RYA-785 gate's own numbers?

The gate drives Turbospectrum directly (`interpol_modeles_nlte` + `bsyn`). The production
path goes through iSpec's wrapper. If the wiring is right, the SAME deck on the SAME lines
at the SAME abundance must give the same equivalent widths. This compares them line by
line, against numbers that were produced before this adapter existed.

Validate-don't-tune: nothing here is fitted. Both legs synthesise at the deck's own
A(Fe) = 7.46 and the EWs are integrated identically to the gate (+/-1.2 A).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, os.environ.get("ISPEC_DIR", "/srv/codex/engines/ispec_src"))

from pipeline._numcompat import trapezoid            # noqa: E402
from pipeline import gerber_nlte as gn               # noqa: E402

# The RYA-785 gate's own result, from its provenance record. Reference, never a target.
GATE_EW_NLTE = {5905.671: 63.86, 6597.559: 41.51, 6705.117: 57.28,
                6726.666: 42.26, 6828.591: 58.85}
A_DECK = 7.46
EW_HW = 1.2          # ts_gerber_gate.EW_HW
WSTEP_NM = 0.0002


def ew_mA(w_nm, flux, centre_A, hw=EW_HW):
    w_A = np.asarray(w_nm) * 10.0
    m = (w_A > centre_A - hw) & (w_A < centre_A + hw)
    return float(trapezoid(1.0 - np.asarray(flux)[m], w_A[m]) * 1000.0)


def main() -> None:
    from scripts.control_synthesis_handler import build_context
    from pipeline.abundances_derive import _synth_flux_at_abund
    from config.constants import STAR_PARAMS

    p = STAR_PARAMS["solar"]
    teff, logg = float(p["teff"]), float(p["logg"])
    feh = float(p.get("feh", p.get("feh_ref", 0.0)))
    ctx = build_context("Fe", "I", 500000.0)
    # THE DECK IS MARCS. The production default is ATLAS9.Castelli (72 layers, and only 10
    # columns so index 7 -- the column iSpec writes as the departure tau -- is not a tau at
    # all and reads as zeros). NLTEgrid4TS_Fe_MARCS was computed on MARCS, and iSpec ships
    # MARCS.GES: 56 layers, col7 = -4.9154..1.7709 against the departure grid's
    # -4.9177..1.7744. Same structure, same depth scale. Using ATLAS9 here would apply
    # MARCS departures to a different atmosphere at mismatched depths.
    from pipeline.abundances_derive import _load_atmosphere
    ctx["atmosphere"] = _load_atmosphere(teff, logg, feh, ctx["vturb"],
                                         model_grid="MARCS.GES")
    print(f"atmosphere: MARCS.GES, {len(ctx['atmosphere'])} layers")

    dep = gn.for_node("Fe", teff, logg, feh)
    print(f"deck: atom.fe607a  ndep={dep['ndep']} nk={dep['nk']}  "
          f"A_deck={dep['deck_abundance']}")
    print(f"node: teff={teff:.0f} logg={logg:.2f} feh={feh:+.2f}\n")

    print(f"{'line':>10} {'EW_LTE':>8} {'EW_NLTE':>8} {'delta_EW':>9} "
          f"{'gate EW_NLTE':>13} {'diff':>7}")
    rows = []
    for c in sorted(GATE_EW_NLTE):
        lo_nm, hi_nm = (c - 3.0) / 10.0, (c + 3.0) / 10.0
        sw = np.arange(lo_nm, hi_nm + WSTEP_NM * 0.5, WSTEP_NM)
        kw = dict(atmosphere=ctx["atmosphere"], teff=teff, logg=logg, feh=feh,
                  vturb=ctx["vturb"], linelist=ctx["linelist"], isotopes=ctx["isotopes"],
                  solar_abund=ctx["solar_abund"], element="Fe", atom_code=26,
                  R=0, macroturbulence=0.0, vsini=0.0)

        gn.assert_linelist_supports_nlte(ctx["linelist"], 26, "Fe")
        gn.assert_depth_match(dep, ctx["atmosphere"])
        nd = {"Fe": gn.as_ispec_tuple(dep, A_DECK)}

        f_lte = _synth_flux_at_abund(sw, trial_A=A_DECK, **kw)
        f_nlte = _synth_flux_at_abund(sw, trial_A=A_DECK, nlte_departures=nd, **kw)

        e_l, e_n = ew_mA(sw, f_lte, c), ew_mA(sw, f_nlte, c)
        g = GATE_EW_NLTE[c]
        rows.append((c, e_l, e_n, g))
        print(f"{c:10.3f} {e_l:8.2f} {e_n:8.2f} {e_n - e_l:+9.2f} {g:13.2f} "
              f"{e_n - g:+7.2f}")

    # ── the equivalence test that matters: reproduce the GATE'S DELTA, not its EWs ──
    # Absolute EWs cannot agree: the gate synthesises against Turbospectrum's own
    # nlte_ges_linelist_jmg17feb2022_I_II while production uses our GES list, so the blends
    # inside the +/-1.2 A window differ (6597.559 is +46 mA on that account alone). The
    # DELTA is what the product carries, and it is what has to reproduce.
    print()
    A = np.array([A_DECK + o for o in (-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3)])
    print(f"{'line':>10} {'A*':>8} {'delta':>8}   (gate median +0.0579)")
    deltas = []
    for c, e_l, e_n, g in rows:
        lo_nm, hi_nm = (c - 3.0) / 10.0, (c + 3.0) / 10.0
        sw = np.arange(lo_nm, hi_nm + WSTEP_NM * 0.5, WSTEP_NM)
        kw = dict(atmosphere=ctx["atmosphere"], teff=teff, logg=logg, feh=feh,
                  vturb=ctx["vturb"], linelist=ctx["linelist"], isotopes=ctx["isotopes"],
                  solar_abund=ctx["solar_abund"], element="Fe", atom_code=26,
                  R=0, macroturbulence=0.0, vsini=0.0)
        cog = [ew_mA(sw, _synth_flux_at_abund(sw, trial_A=float(a), **kw), c) for a in A]
        a_star = float(np.interp(e_n, cog, A))
        d_ = A_DECK - a_star
        deltas.append(d_)
        print(f"{c:10.3f} {a_star:8.4f} {d_:+8.4f}")
    dv = np.array(deltas)
    print()
    print(f"  n={len(dv)} median {np.median(dv):+.4f}  mean {np.mean(dv):+.4f}  "
          f"sd {dv.std(ddof=1):.4f}")
    print(f"  gate (RYA-785, 7 isolated lines): median +0.0579")
    print(f"  |this - gate| = {abs(np.median(dv) - 0.0579):.4f} dex")

    d = np.array([r[2] - r[3] for r in rows])
    same = np.array([abs(r[1] - r[2]) < 1e-6 for r in rows])
    print(f"\n  vs the gate: median {np.median(d):+.2f} mA, max |diff| {np.abs(d).max():.2f} mA")
    if same.any():
        print(f"  ⚠️  {int(same.sum())} line(s) gave IDENTICAL LTE and NLTE flux — "
              f"departures did NOT engage (silent LTE, RYA-764).")
    else:
        print(f"  NLTE differs from LTE on all {len(rows)} lines: departures engaged.")


if __name__ == "__main__":
    main()
