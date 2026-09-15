#!/usr/bin/env python3
"""RYA-1214 — what `SOLAR_VIS_GATES` actually tests, and two ways it is not what it reads as.

`cno_synthesis.validate_solar` prints a PASS/FAIL table against Asplund 2021 and returns an
all-pass boolean. Running the re-measurement this ticket owes surfaced two problems with the
test itself. Both are properties of the gate, not of the numbers, and neither is visible
from the table it prints.

🔴 1. THE TOLERANCE INCLUDES THE PRODUCT'S OWN sigma_tot, SO IMPRECISION BUYS A PASS.

    within = abs(val - tgt) <= (terr + sig)          # sig = result.uncertainty[q]['tot']

Nitrogen this run is A(N) = 7.384 against a target of 7.83 — off by 0.446 dex, more than
SIX TIMES the gate's own 0.07 tolerance — and it PASSES, because sigma_tot is 1.135 and
1.135 + 0.07 = 1.205. A measurement that determined nothing clears a gate that a good
measurement 0.1 dex off would fail. The gate rewards a wide bar.

⚠️ AND THIS TICKET MADE IT WORSE BEFORE IT MADE IT VISIBLE. Retiring the
`np.clip(sigma_fit, 0.0, 1.0)` sentinel let nitrogen's sigma report its true value, which is
LARGER than the clip allowed (one iteration read 1.38). So the honest sigma widens the
tolerance and makes the gate easier to pass than the dishonest one did. Fixing the sentinel
without saying this would have quietly loosened an acceptance test.

🔴 2. IT COMPARES A 1D-LTE VALUE TO A 3D REFERENCE.

The targets are AGSS21's published abundances, which are 3D. `val` is
`result.abundances[q]`, the 1D-LTE fit, taken BEFORE the Phase-A cited-correction layer
that exists precisely to move 1D-LTE onto the cited 3D/NLTE scale. For [O I] 6300 that
layer applies Caffau et al. 2015's CO5BOLD -0.080 dex, so the gate is testing a number
0.080 dex away from the frame its target lives in — a mismatch built into the test.

Comparing like with like changes the answer for oxygen and for C/O, in the direction of
agreement, which is exactly why it must be stated rather than quietly adopted: 8.810 fails
a 0.05 gate on |delta| alone, 8.730 passes it.

WHAT THIS SCRIPT DOES. It recomputes the gate three ways on the emitted product — as coded,
with the sigma term removed, and on the Phase-A corrected values — and reports all three. It
changes no gate and no number; `SOLAR_VIS_GATES` is ratified science and moving it is Ryan's
call, not a side effect of a re-measurement.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CNOSYNTH = ROOT / "data" / "audit" / "cno_synthesis"
OUT = ROOT / "data" / "audit" / "rya1214_cno_products"

#: Quoted from `cno_synthesis.SOLAR_VIS_GATES` rather than re-typed — if that moves, this
#: audit must move with it, and a copy would let the two disagree (RYA-845).
def _gates() -> dict:
    src = (ROOT / "pipeline" / "cno_synthesis.py").read_text()
    i = src.index("SOLAR_VIS_GATES = {")
    j = src.index("}", i)
    return eval(src[i + len("SOLAR_VIS_GATES = "): j + 1])   # noqa: S307 - a literal dict


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    gates = _gates()
    prod = pd.read_csv(CNOSYNTH / "solar_vis_cno_product.csv").set_index("element")
    bands = pd.read_csv(CNOSYNTH / "solar_vis_cno_per_band.csv")
    prov = json.loads((CNOSYNTH / "solar_vis_cno_provenance.json").read_text())

    # The Phase-A cited corrections, per element, from the run's own provenance.
    # ⚠️ NESTED under `provenance`, not top level — `_write_product` puts
    # `phase_a_corrections` inside `result.provenance` (cno_synthesis.py:1304) while
    # `abundances` / `uncertainty` sit at the root. A top-level `.get` returns [] and the
    # audit silently reports "0 verdicts change", which is the answer you get when you
    # never looked (RYA-833).
    _prov_inner = prov.get("provenance", {}) or {}
    _corrections = (prov.get("phase_a_corrections")
                    or _prov_inner.get("phase_a_corrections") or [])
    if not _corrections:
        raise SystemExit(
            "no phase_a_corrections found in the provenance — refusing to report "
            "'0 verdicts change' from an empty list. Check the key path.")
    corr = {}
    for c in _corrections:
        if c.get("kind") == "cited_3d_anchor" and c.get("a_corr") is not None:
            corr[c["element"]] = (float(c["a_lte"]), float(c["a_corr"]),
                                  float(c["delta"]), c.get("source", ""))

    rows = []
    for q, (tgt, terr) in gates.items():
        if q == "C/O":
            a_lte = float(prod.loc["C", "A_X"])
            a_o = float(prod.loc["O", "A_X"])
            val = 10 ** (a_lte - a_o)
            sig = 0.0
            o_c = corr.get("O")
            val_corr = 10 ** (a_lte - o_c[1]) if o_c else None
        else:
            val = float(prod.loc[q, "A_X"])
            s = prod.loc[q, "sigma_tot"]
            sig = float(s) if pd.notna(s) else 0.0
            val_corr = corr[q][1] if q in corr else None
        rows.append({
            "qty": q, "target": tgt, "tol": terr,
            "ours_1D_LTE": round(val, 3), "sigma_tot": (round(sig, 3) if sig else None),
            "abs_delta": round(abs(val - tgt), 3),
            "as_coded_PASS": bool(abs(val - tgt) <= terr + sig + 1e-9),
            "tolerance_as_coded": round(terr + sig, 3),
            "sigma_free_PASS": bool(abs(val - tgt) <= terr + 1e-9),
            "ours_phase_a_corrected": (round(val_corr, 3) if val_corr is not None else None),
            "corrected_abs_delta": (round(abs(val_corr - tgt), 3)
                                    if val_corr is not None else None),
            "corrected_sigma_free_PASS": (bool(abs(val_corr - tgt) <= terr + 1e-9)
                                          if val_corr is not None else None),
        })
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "solar_vis_gate_audit.csv", index=False)

    pd.set_option("display.width", 220)
    print("=== RYA-1214 — SOLAR_VIS_GATES, recomputed three ways ===")
    print(d[["qty", "target", "tol", "ours_1D_LTE", "sigma_tot", "abs_delta",
             "tolerance_as_coded", "as_coded_PASS", "sigma_free_PASS",
             "ours_phase_a_corrected", "corrected_abs_delta",
             "corrected_sigma_free_PASS"]].to_string(index=False))

    bought = d[d.as_coded_PASS & ~d.sigma_free_PASS]
    print(f"\n🔴 gates that pass ONLY because sigma_tot widens the tolerance: {len(bought)}")
    for _, r in bought.iterrows():
        print(f"   {r.qty}: |delta| {r.abs_delta} vs tol {r.tol} — passes on "
              f"{r.tolerance_as_coded} because sigma_tot = {r.sigma_tot}")
    moved = d[d.ours_phase_a_corrected.notna() &
              (d.sigma_free_PASS != d.corrected_sigma_free_PASS)]
    print(f"\n🔴 verdicts that CHANGE once the 1D-LTE value is moved onto the "
          f"reference's own 3D frame: {len(moved)}")
    for _, r in moved.iterrows():
        print(f"   {r.qty}: 1D-LTE {r.ours_1D_LTE} (|d| {r.abs_delta}) -> corrected "
              f"{r.ours_phase_a_corrected} (|d| {r.corrected_abs_delta})  "
              f"{r.sigma_free_PASS} -> {r.corrected_sigma_free_PASS}")
    for el, (a, ac, dl, src) in sorted(corr.items()):
        print(f"\n   cited correction applied to {el}: {a:.3f} -> {ac:.3f} ({dl:+.3f}) "
              f"— {src[:80]}")

    (OUT / "solar_vis_gate_audit.prov.json").write_text(json.dumps({
        "ticket": "RYA-1214", "read_only": True,
        "what": "SOLAR_VIS_GATES recomputed as coded, sigma-free, and on the Phase-A "
                "cited-3D values",
        "defect_1": "the tolerance includes the product's own sigma_tot, so imprecision "
                    "buys a pass. A(N) is 0.446 dex off a 0.07 gate and passes on a "
                    "1.135 dex bar. ⚠️ Retiring the RYA-1214 sigma clip made this WORSE "
                    "before it made it visible: the honest sigma is larger than the clip "
                    "allowed, so it widens the tolerance further.",
        "defect_2": "it compares result.abundances (1D-LTE) against AGSS21 targets that "
                    "are 3D, BEFORE the Phase-A cited-correction layer that exists to "
                    "bridge exactly that gap. For [O I] 6300 the gap is Caffau 2015's "
                    "-0.080 dex, and the oxygen verdict turns on it.",
        "changed_nothing": "no gate and no number was altered. SOLAR_VIS_GATES is ratified "
                           "science; moving it is Ryan's call, not a side effect of a "
                           "re-measurement.",
        "rows": rows,
    }, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
