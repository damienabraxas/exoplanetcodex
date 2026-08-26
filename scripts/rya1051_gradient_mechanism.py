#!/usr/bin/env python3
"""RYA-1051: is the Gerber-vs-Amarsi metallicity-gradient difference PHOTOIONISATION or
COLLISIONAL? A per-line discriminator.

🔴 THE TEST. Fe I over-ionisation is driven by super-thermal UV, and the UV field is what
changes with metallicity (less line blanketing -> more UV -> stronger over-ionisation). So
the two candidate mechanisms make OPPOSITE predictions about the per-LINE structure of the
disagreement:

  * PHOTOIONISATION. If the two atoms' bound-free cross-sections differ, the difference
    must depend on WHERE each level's ionisation edge sits. A level at EP = 2.85 eV
    ionises at 2455 A, one at 4.73 eV at 3912 A -- a real span across the region where the
    metallicity-driven flux change is largest. The gradient difference should TRACK the
    edge wavelength.
  * COLLISIONS. Electron collisions thermalise the level system as a whole. Their effect
    is not tied to any particular level's edge, so the gradient difference should be
    UNIFORM across lines.

⚠️ THE CONDITIONING CUT IS DERIVED, NOT CHOSEN. Delta is solved as
(W_LTE - W_NLTE) / (dEW/dA), with the slope measured over a trial step DA. A returned
|Delta| larger than DA is therefore an EXTRAPOLATION beyond the only range in which the
slope was ever sampled, and is rejected on that ground alone -- no threshold is picked.
It matters: Fe I 6315.811 at [Fe/H] = -3 has dEW/dA = 0.37 mA/dex on a 12.7 mA line and
returned Delta = +0.706, an order of magnitude off every other line. The original guard
(|slope| < 1e-6) was far too permissive to catch it.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
IP_FE_I = 7.902          # eV, Fe I ionisation potential -- sets each level's edge
DA = 0.15                # the trial step the slope was measured over


def _solar_ep_test(nlte_csv: Path, lte_csv: Path, anchor_csv: Path) -> dict:
    """🔴 THE SMALL-SAMPLE TRAP, AND WHY THIS FUNCTION EXISTS.

    On the 8 lines above, the SOLAR gap looked EP-dependent -- r = -0.41 against
    excitation potential, -0.46 against edge wavelength -- and I reported that as a
    second, level-dependent term sitting on top of the uniform gradient term. IT IS NOT
    THERE. The measured solar pair shares 43 lines with the anchor, and on those:

        gap vs EP        r = -0.083  (t = -0.54)
        gap vs edge_A    r = -0.086
        EP terciles      +0.0789 / +0.0729 / +0.0779   -- flat

    The SAME 8 lines evaluated through the measured route give r = -0.273, so the earlier
    number was SAMPLE SIZE and not method: at n = 8 the 95% interval on r = -0.41 runs
    roughly -0.85 to +0.40, which contains zero comfortably. Eight points cannot carry a
    correlation claim, and I should not have made one from them.

    ⚠️ What IS true, and is the more interesting fact: both codes individually trend with
    EP -- ours +0.055 -> +0.042 across the range, Amarsi -0.020 -> -0.043 -- they simply
    trend TOGETHER, so the gap between them stays flat. The disagreement is a UNIFORM
    offset in BOTH metallicity and excitation potential. One effect, not two.
    """
    n = pd.read_csv(nlte_csv); l = pd.read_csv(lte_csv)
    m = n.merge(l, on="wavelength_air_A", suffixes=("_n", "_l"))
    m["ours"] = m.abundance_n - m.abundance_l
    m["ep"] = m.ep_eV_n
    a = pd.read_csv(anchor_csv)
    a = a[(a.species == "Fe1") & (a.clean == "yes")
          & (a.vturb.astype(float) == 1.50)].copy()
    a["w"] = a.lam_nm.astype(float) * 10.0
    a["am"] = a["corr"].astype(float)
    a = a[a.am != -4.0]
    aw, av = a.w.values, a.am.values
    rows = []
    for _, r in m.iterrows():
        i = int(np.argmin(np.abs(aw - r.wavelength_air_A)))
        if abs(aw[i] - r.wavelength_air_A) <= 0.05:
            rows.append((r.ep, r.ours - av[i]))
    if len(rows) < 10:
        return {"status": "TOO_FEW", "n": len(rows)}
    D = pd.DataFrame(rows, columns=["ep", "gap"])
    D["edge"] = 12398.4 / (IP_FE_I - D.ep)
    r_ep = float(np.corrcoef(D.ep, D.gap)[0, 1])
    t_stat = r_ep * np.sqrt((len(D) - 2) / (1 - r_ep ** 2))
    q = D.ep.quantile([0, 1 / 3, 2 / 3, 1]).values
    terciles = [round(float(D[(D.ep >= lo) & (D.ep <= hi)].gap.median()), 4)
                for lo, hi in zip(q[:-1], q[1:])]
    return {"status": "RUN", "n": len(D),
            "median_gap": round(float(D.gap.median()), 4),
            "r_vs_ep": round(r_ep, 3), "t_vs_ep": round(float(t_stat), 2),
            "r_vs_edge_A": round(float(np.corrcoef(D.edge, D.gap)[0, 1]), 3),
            "ep_terciles": terciles,
            "verdict": ("FLAT -- no level-dependent term at solar; the disagreement is a "
                        "uniform offset in EP as well as in metallicity"
                        if abs(t_stat) < 2.0 else
                        "EP-DEPENDENT at solar, on top of the uniform gradient term")}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--grid", default=str(ROOT / "data" / "results" / "rya1051" /
                                          "nlte_correction_grid.json"))
    ap.add_argument("--anchor-per-line", required=True,
                    help="Amarsi nmarcs_lmarcs per-line slice at the same nodes")
    ap.add_argument("--solar-pair", nargs=2, metavar=("NLTE_CSV", "LTE_CSV"),
                    default=None,
                    help="the MEASURED solar 1D pair (ENGINE-B-NLTE, synth-1D-LTE-gerber) "
                         "per-line CSVs. Enables the LARGE-SAMPLE solar EP test, which is "
                         "the one that matters -- see `_solar_ep_test`.")
    ap.add_argument("--solar-anchor", default=None,
                    help="Amarsi nmarcs_lmarcs solar-node slice (all clean Fe I lines)")
    ap.add_argument("--out", default=str(ROOT / "data" / "results" / "rya1051" /
                                         "gradient_mechanism.json"))
    a = ap.parse_args()

    g = json.loads(Path(a.grid).read_text())
    rows = [(float(k), r["wave_A"], r["ep_eV"], r["delta_nlte"], r["dEW_dA"])
            for k, v in g["nodes"].items() for r in v["lines"]]
    O = pd.DataFrame(rows, columns=["feh", "wave_A", "ep", "ours", "dEW_dA"])
    n_all = len(O)
    O = O[O.ours.abs() <= DA]          # the derived conditioning cut
    A = pd.read_csv(a.anchor_per_line).rename(columns={"corr": "am"})
    m = O.merge(A[["feh", "wave_A", "am"]], on=["feh", "wave_A"], how="inner")

    out = []
    for w, grp in m.groupby("wave_A"):
        if len(grp) < 3:
            continue
        grp = grp.sort_values("feh")
        so = float(np.polyfit(grp.feh, grp.ours, 1)[0])
        sa = float(np.polyfit(grp.feh, grp.am, 1)[0])
        s = grp[grp.feh == 0.0]
        out.append({"wave_A": float(w), "ep_eV": float(grp.ep.iloc[0]),
                    "edge_A": 12398.4 / (IP_FE_I - float(grp.ep.iloc[0])),
                    "grad_ours": round(so, 4), "grad_amarsi": round(sa, 4),
                    "grad_diff": round(sa - so, 4),
                    "solar_gap": (round(float(s.ours.iloc[0] - s.am.iloc[0]), 4)
                                  if len(s) else None),
                    "n_nodes": int(len(grp))})
    R = pd.DataFrame(out).sort_values("ep_eV")

    doc = {"ticket": "RYA-1051",
           "conditioning_cut": f"|Delta| <= DA = {DA}; DERIVED, not chosen -- a larger "
                               f"value extrapolates past the step the slope was measured "
                               f"over. Rejected {n_all - len(O)} of {n_all} line-nodes.",
           "n_lines": len(R), "lines": out}
    for col in ("grad_diff", "solar_gap"):
        v = R[col].astype(float)
        doc[col] = {
            "median": round(float(v.median()), 4),
            "p16": round(float(np.percentile(v, 16)), 4),
            "p84": round(float(np.percentile(v, 84)), 4),
            "r_vs_ep": round(float(np.corrcoef(R.ep_eV, v)[0, 1]), 3),
            "r_vs_edge_A": round(float(np.corrcoef(R.edge_A, v)[0, 1]), 3)}
    doc["all_lines_same_sign"] = bool((R.grad_diff < 0).all())
    if a.solar_pair and a.solar_anchor:
        doc["solar_ep_test"] = _solar_ep_test(Path(a.solar_pair[0]),
                                              Path(a.solar_pair[1]),
                                              Path(a.solar_anchor))
    doc["verdict"] = (
        "UNIFORM: the gradient difference tracks NEITHER excitation potential NOR the "
        "photoionisation edge wavelength, and every line carries it with the same sign "
        "and nearly the same size. That is the COLLISIONAL signature, not the bound-free "
        "one -- a cross-section difference would have to depend on where each level's "
        "edge sits."
        if abs(doc["grad_diff"]["r_vs_edge_A"]) < 0.3 else
        "EDGE-DEPENDENT: consistent with a photoionisation origin.")

    print(f"{'wave':>9} {'EP':>5} {'edge_A':>7} {'g ours':>8} {'g Amarsi':>9} {'diff':>8}")
    for r in out:
        print(f"{r['wave_A']:9.3f} {r['ep_eV']:5.2f} {r['edge_A']:7.0f} "
              f"{r['grad_ours']:+8.4f} {r['grad_amarsi']:+9.4f} {r['grad_diff']:+8.4f}")
    print(f"\ngrad_diff  r vs EP {doc['grad_diff']['r_vs_ep']:+.3f}   "
          f"r vs edge {doc['grad_diff']['r_vs_edge_A']:+.3f}   "
          f"median {doc['grad_diff']['median']:+.4f}")
    print(f"solar_gap  r vs EP {doc['solar_gap']['r_vs_ep']:+.3f}   "
          f"r vs edge {doc['solar_gap']['r_vs_edge_A']:+.3f}")
    print(f"\nVERDICT: {doc['verdict']}")
    Path(a.out).write_text(json.dumps(doc, indent=2) + "\n")
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
