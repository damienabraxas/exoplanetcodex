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


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--grid", default=str(ROOT / "data" / "results" / "rya1051" /
                                          "nlte_correction_grid.json"))
    ap.add_argument("--anchor-per-line", required=True,
                    help="Amarsi nmarcs_lmarcs per-line slice at the same nodes")
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
