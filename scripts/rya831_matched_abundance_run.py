#!/usr/bin/env python3
"""
RYA-831 (MAC) — is RYA-819's +0.033 Fe 3D-atmosphere term a DOMAIN-RAIL artifact?

RYA-819 got two answers for A(3D-NLTE) - A(1D-NLTE): +0.033 from running the Amarsi-2022
MLP on our lines, and about -0.013 from Amarsi's own reproduced solar row. Both say gold's
Magic-2013 -0.05 is wrong; they disagree on magnitude, and the magnitude is what the
RYA-819 re-freeze note needs.

THE SUSPECT. The MLP's A(Fe) axis runs 4.5 .. 7.5 (read from the vendored grid). Our
per-line A(1D-LTE) has median 7.568 — ABOVE the ceiling — because it carries the RYA-161
gf zero-point. Amarsi's sit near 7.467, mid-axis. A railed input returns what the network
does at its boundary, which is not a measurement of the Sun.

⚠️ THE AXIS IS PINNED WITH `_compute_aberr`, WHICH IS HOW RYA-817's OWN CONTROL DOES IT.
My first attempt called `aberr_for_line` in per-line mode and let each line converge. That
REFUSED 17 of Amarsi's own 41 control lines — their converged A(Fe;3N) lands at 7.503 to
7.620, over the ceiling — and the surviving 24 reproduced 7.4325 against a published 7.47.
The control caught it. 817 pins the axis at ASPLUND21_FE = 7.46 and queries the network
there, which is exactly the in-domain query this ticket asks for; reusing that primitive
rather than a second one of my own is the whole point.

⚠️ THE MATCHED ABUNDANCE IS COMPUTED, NOT CHOSEN. The sweep is reported whole, so the term
is visible at every axis value including our railed one. Nothing is adjusted to make it
agree with anything.
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.nlte_corrections import _compute_aberr                   # noqa: E402

DECOMP = ROOT / "data" / "results" / "rya819" / "rya819_decomposition_per_line.csv"
CONTROL_LINES = (ROOT / "data" / "reference" / "amarsi2022_training"
                 / "amarsi2022_solar_control_lines.csv")
OUT = ROOT / "data" / "results" / "rya831"

TEFF, LOGG, VMIC = 5772.0, 4.438, 1.0
ASPLUND21_FE = 7.46            # the axis RYA-817's control pins, and the paper's own scale
AMARSI_1DLTE = 7.4671          # Amarsi+2022 Table 6 solar A(Fe I) 1D LTE
PAPER = {"1D_LTE": 7.47, "3D_NLTE": 7.46}
CONTROL_TOL = 0.01
AXIS_CEILING = 7.5


def reactivation_control() -> dict:
    """Reproduce Amarsi+2022 Table 6's solar Fe I row on HIS lines, on THIS machine."""
    c = pd.read_csv(CONTROL_LINES)
    sub = c[(c.ion.astype(str).str.strip() == "I") & c.weak_line_rew_lt_m49]
    ab = np.array([_compute_aberr("I", r.elo_eV, r.eup_eV, r.loggf,
                                  TEFF, LOGG, ASPLUND21_FE, VMIC)
                   for r in sub.itertuples()], float)
    a1d = float(sub.a_1d_lte_ap2002.mean())
    a3n = a1d + float(np.nanmean(ab))
    ok = (abs(a1d - PAPER["1D_LTE"]) <= CONTROL_TOL
          and abs(a3n - PAPER["3D_NLTE"]) <= CONTROL_TOL)
    print(f"REACTIVATION CONTROL on this machine  (Mac, python "
          f"{sys.version_info.major}.{sys.version_info.minor}, n={len(sub)} weak Fe I)")
    print(f"  A(Fe I) 1D LTE      {a1d:.4f}   published {PAPER['1D_LTE']}")
    print(f"  A(Fe I) 3D non-LTE  {a3n:.4f}   published {PAPER['3D_NLTE']}")
    print(f"  net (3D-NLTE - 1D-LTE)  {a3n - a1d:+.4f}")
    print(f"  -> {'PASS' if ok else 'FAIL'}")
    if not ok:
        raise SystemExit("STOP: the reactivation does not reproduce Table 6 here. "
                         "Every number below would be a different engine's answer.")
    return {"n": int(len(sub)), "a_1dlte": a1d, "a_3dnlte": a3n, "net": a3n - a1d,
            "mean_correction": float(np.nanmean(ab)), "pass": bool(ok)}


def main() -> None:
    ctrl = reactivation_control()

    d = pd.read_csv(DECOMP)
    d = d[(d.in_domain == True) & d.delta_nlte.notna()].copy()         # noqa: E712
    n_railed = int((d.a_1dlte > AXIS_CEILING).sum())
    med_1dlte = float(d.a_1dlte.median())
    dnlte = float(d.delta_nlte.median())

    print(f"\nOUR LINES: {len(d)} in-domain Fe I VIS carrying an MPIA delta")
    print(f"  A(1D-LTE) median              {med_1dlte:.4f}")
    print(f"  ABOVE the {AXIS_CEILING} axis ceiling     {n_railed}"
          f"  ({100 * n_railed / len(d):.0f}%)  <- queried at the boundary")
    print(f"  measured gf zero-point offset {med_1dlte - AMARSI_1DLTE:+.4f} dex "
          f"(ours - Amarsi's)")
    print(f"  delta_NLTE (MPIA, median)     {dnlte:+.4f}")

    axes = [7.10, 7.20, 7.30, 7.40, ASPLUND21_FE, AMARSI_1DLTE, 7.50]
    print(f"\n3D-ATMOSPHERE TERM vs WHERE ON THE A(Fe) AXIS THE NETWORK IS ASKED")
    print(f"  {'axis A(Fe)':>11}{'net 3D-NLTE':>14}{'3D-atmosphere':>16}")
    curve = []
    for ax in axes:
        ab = np.array([_compute_aberr("I", r.elo_eV, r.eup_eV, r.loggf,
                                      TEFF, LOGG, float(ax), VMIC)
                       for r in d.itertuples()], float)
        net = float(np.nanmedian(ab))
        t3d = float(np.nanmedian(ab - d.delta_nlte.to_numpy(float)))
        tag = ("  <- RYA-817 control axis" if ax == ASPLUND21_FE else
               "  <- Amarsi's 1D-LTE" if ax == AMARSI_1DLTE else
               "  <- THE RAIL (where 819 asked)" if ax == AXIS_CEILING else "")
        print(f"  {ax:>11.4f}{net:>+14.4f}{t3d:>+16.4f}{tag}")
        curve.append({"axis_afe": float(ax), "net_3dnlte": net,
                      "term_3d_atmosphere": t3d})

    at_rail = next(c for c in curve if c["axis_afe"] == AXIS_CEILING)
    at_asp = next(c for c in curve if c["axis_afe"] == ASPLUND21_FE)
    amarsi_term = ctrl["net"] - dnlte

    print(f"\nVERDICT INPUTS")
    print(f"  term at the RAIL (7.50)           {at_rail['term_3d_atmosphere']:+.4f}"
          f"   <- RYA-819 reported +0.033 here")
    print(f"  term at Asplund/Amarsi (7.46)     {at_asp['term_3d_atmosphere']:+.4f}")
    print(f"  moved by de-railing               "
          f"{at_asp['term_3d_atmosphere'] - at_rail['term_3d_atmosphere']:+.4f} dex")
    print(f"  Amarsi's own in-domain estimate   {amarsi_term:+.4f}"
          f"   (his net {ctrl['net']:+.4f} minus delta_NLTE {dnlte:+.4f})")
    print(f"  Magic-2013 as gold applies it     -0.0500")
    print(f"  gold over-corrects by             "
          f"{at_asp['term_3d_atmosphere'] - (-0.05):+.4f} dex (at the matched axis)")

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "rya831_axis_sensitivity.json").write_text(json.dumps({
        "ticket": "RYA-831",
        "machine": f"Mac, python {sys.version.split()[0]}, sklearn present "
                   f"(absent from all four Sirius venvs)",
        "method": "pipeline.nlte_corrections._compute_aberr with the A(Fe;3N) axis "
                  "PINNED — the same primitive RYA-817's reactivation control uses",
        "reactivation_control": ctrl,
        "grid_afe_axis": [4.5, AXIS_CEILING],
        "n_lines": int(len(d)), "n_railed_above_ceiling": n_railed,
        "our_median_a1dlte": med_1dlte, "amarsi_a1dlte": AMARSI_1DLTE,
        "gf_zero_point_offset_dex": med_1dlte - AMARSI_1DLTE,
        "delta_nlte_median": dnlte,
        "curve": curve,
        "term_at_rail": at_rail["term_3d_atmosphere"],
        "term_at_matched_axis": at_asp["term_3d_atmosphere"],
        "amarsi_in_domain_term": amarsi_term,
        "magic_as_applied": -0.05,
    }, indent=2, default=float))
    pd.DataFrame(curve).to_csv(OUT / "rya831_axis_sensitivity.csv", index=False)
    print(f"\n[out] {OUT}/rya831_axis_sensitivity.json")


if __name__ == "__main__":
    main()
