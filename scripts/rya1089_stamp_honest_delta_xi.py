#!/usr/bin/env python3
"""RYA-1089: derive the solar Type B budget from the CODE's delta_p, and flag the gate.

🔴 WHY THIS EXISTS. `data/audit/uncertainty/solar_uncertainty_rya158.json` had no
generator and no GENERATORS.yaml entry — it was hand-maintained. So when RYA-1093 re-derived
the solar microturbulence allowance, the code moved and the artifact did not:

⚠️ AND IT COULD NOT HAVE HAD AN ENTRY. `check_result_generators.py` scans `data/results/`
and `data/processed/` ONLY, and refuses an entry naming anything outside them (tried:
"names an artifact that is not tracked under data/results/, data/processed/"). So RYA-686's
manifest — the mechanism whose entire purpose is "no artifact without its generator" —
structurally cannot cover `data/audit/`, which is where the solar uncertainty budget lives.
That blind spot is what let this file go stale for a whole ticket cycle. Widening the scan
is its own ticket; until then `tests/test_honest_delta_xi_rya1089.py` is the guard, and it
asserts the artifact against the CODE rather than against a manifest.

    pipeline/uncertainty_stack.py   SOLAR_DELTA_XI_KMS = 0.2912   (RYA-1093)
    solar_uncertainty_rya158.json   delta_p_solar.vmic = 0.05     (RYA-158, retired)

Both were committed on main at the same time, and **nothing compared them**. The three
readers of the artifact — `rya1088_record_sigma_params.py`,
`tests/test_uncertainty_stack_rya158.py`, `tests/test_solar_calibration_gate.py` — all pass,
because none of them asks whether the file agrees with the module that defines the value.
A consumer reading the artifact got the retired 0.05; a consumer calling the code got 0.2912.

⚠️ THE STALE VALUE IS THE ONE THAT FLATTERS US. 0.05 km/s gives sigma_B_vmic = 0.0120 dex,
comfortably inside the 0.05 dex solar gate. 0.2912 gives 0.0699 dex, which FAILS it. An
artifact that silently disagrees with the code in the direction of passing a gate is the
RYA-161 hazard arriving by accident rather than by choice, and it is why this generator
reads delta_p from the code rather than carrying its own copy.

WHAT IS DERIVED AND WHAT IS READ:

  * `delta_p` is READ FROM THE CODE — `uncertainty_stack.params_and_deltas('solar')`. It is
    never written here, so this artifact cannot drift from the module again.
  * The per-element MEASUREMENTS are read from the existing artifact and passed through
    untouched: `n_lines`, `A_X_mean`, `raw_sigma`, `flags`, and the measured derivatives
    `dA_dTeff_per100K` / `dA_dvmic_per_kms`. A derivative is a property of the star and the
    line set; changing the STEP SIZE we probe it with does not change it.
  * `sigma_B_*`, `sigma_solar` and `sigma_reported` are RECOMPUTED from those two.

🔴 --control IS THE PROOF THE GENERATOR IS FAITHFUL. Run with the artifact's OWN delta_p it
must reproduce the committed file EXACTLY. Regenerating a file with a new input, without
first showing the generator reproduces the old output from the old input, cannot tell "the
input changed" from "my generator is wrong". Run it before trusting the re-stamp:

    python3 scripts/rya1089_stamp_honest_delta_xi.py --control     # must print REPRODUCED
    python3 scripts/rya1089_stamp_honest_delta_xi.py               # then re-stamp

⚠️ ONLY Fe I MOVES. Every other element carries `dA_dvmic_per_kms = NaN` — the derivative
was never measured for it — and so reports `sigma_B_vmic = 0.0`. That is RYA-907's
"unmeasured is not zero" defect sitting in this file already; it predates this ticket and is
NOT fixed here, but it is now stated in the artifact's own `caveats` rather than left for a
reader to infer from a zero.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import uncertainty_stack  # noqa: E402

ARTIFACT = ROOT / "data" / "audit" / "uncertainty" / "solar_uncertainty_rya158.json"

#: The solar Type A/B gate. Reported, NEVER met by choosing a smaller bar (RYA-161).
SOLAR_GATE_DEX = 0.05

#: Keys this generator recomputes. Everything else in a per-element row is a measurement and
#: is passed through byte-for-byte.
DERIVED = ("sigma_B_Teff", "sigma_B_logg", "sigma_B_vmic", "sigma_B_FeH",
           "sigma_solar", "sigma_reported")

CAVEATS = [
    ("delta_p_solar is READ FROM pipeline/uncertainty_stack.py at generation time, never "
     "typed here. RYA-1093 set the solar vmic allowance to 0.2912 km/s -- the "
     "method+selection spread |1.0 - 0.709| the 58-line Fe I pool cannot resolve -- "
     "superseding RYA-158's uncited 0.05."),
    ("THE xi TERM ALONE EXCEEDS THE 0.05 dex SOLAR GATE (0.0699 dex) and is reported "
     "anyway. RYA-1093 2E forbids adopting the 0.0588 formal error because it would pass "
     "(RYA-161). The gate is a FLAG here, not a blocker; the adopt-or-relax call is "
     "Ryan's."),
    ("THE DERIVATIVES ARE IN THREE STATES, AND ONLY ONE OF THEM IS A RESULT. Fe I is "
     "MEASURED (dA_dTeff 0.0665, dA_dvmic -0.24). Eight elements (Ba, C, Co, Cu, Mn, N, "
     "O, V) carry NaN -- UNMEASURED -- and their sigma_B = 0.0 means 'not measured', not "
     "'measured as zero' (RYA-907). S I is a THIRD state: it carries an EXACT 0.0 for "
     "BOTH derivatives on n=2 lines (flags n_lines_low). A zero Teff sensitivity is "
     "physically implausible for a real line set, so that 0.0 is more likely a solve "
     "that failed into a default than a measurement. All three states predate RYA-1089 "
     "and NONE is fixed here -- they are stated so no zero in this file is read as a "
     "result. The S I case is worth its own ticket."),
]


def _round(x: float, n: int = 4) -> float:
    return round(x + 0.0, n)


def recompute(doc: dict, delta_p: dict) -> dict:
    """A new doc with delta_p applied and every derived field recomputed."""
    out = json.loads(json.dumps(doc))          # deep copy, NaN preserved
    out["delta_p_solar"] = {"Teff": delta_p["teff_K"], "logg": delta_p["logg"],
                            "vmic": delta_p["vturb_kms"], "FeH": delta_p["feh"]}
    for r in out["per_element"]:
        d_teff = r.get("dA_dTeff_per100K")
        d_vmic = r.get("dA_dvmic_per_kms")
        # A NaN derivative is UNMEASURED. It contributes nothing rather than raising,
        # which is the artifact's existing behaviour -- see CAVEATS.
        sb_teff = (0.0 if d_teff is None or math.isnan(d_teff)
                   else abs(d_teff) * (out["delta_p_solar"]["Teff"] / 100.0))
        sb_vmic = (0.0 if d_vmic is None or math.isnan(d_vmic)
                   else abs(d_vmic) * out["delta_p_solar"]["vmic"])
        # logg and FeH are exact by definition for the Sun -> delta_p = 0 -> term = 0.
        r["sigma_B_Teff"] = _round(sb_teff)
        r["sigma_B_logg"] = 0.0
        r["sigma_B_vmic"] = _round(sb_vmic)
        r["sigma_B_FeH"] = 0.0
        # ⚠️ sigma_SE is STORED ROUNDED but must be used at FULL precision: the committed
        # file only reproduces from raw_sigma/sqrt(n). Using the rounded 0.0178 gives
        # 0.0215 where the file says 0.0214.
        se = r["raw_sigma"] / math.sqrt(r["n_lines"]) if r.get("n_lines") else 0.0
        sigma_params = math.sqrt(sb_teff ** 2 + sb_vmic ** 2)
        rep = math.sqrt(se ** 2 + sigma_params ** 2)
        r["sigma_solar"] = _round(rep)
        r["sigma_reported"] = _round(rep)
    out["caveats"] = CAVEATS
    return out


def _dump(doc: dict) -> str:
    return json.dumps(doc, indent=2)          # no trailing newline, matches the file


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--artifact", default=str(ARTIFACT))
    ap.add_argument("--control", action="store_true",
                    help="recompute with the artifact's OWN delta_p and assert it "
                         "reproduces the committed file; write nothing")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    path = Path(a.artifact)
    raw = path.read_text()
    doc = json.loads(raw)

    if a.control:
        dp = doc["delta_p_solar"]
        got = recompute(doc, {"teff_K": dp["Teff"], "logg": dp["logg"],
                              "vturb_kms": dp["vmic"], "feh": dp["FeH"]})
        # Compare only the DERIVED fields: `caveats` is new, and the control is a claim
        # about the arithmetic, not about the commentary.
        bad = []
        for old, new in zip(doc["per_element"], got["per_element"]):
            for k in DERIVED:
                if old.get(k) != new.get(k):
                    bad.append(f"{old['element']} {old['ion']} {k}: committed "
                               f"{old.get(k)!r} != recomputed {new.get(k)!r}")
        if bad:
            print("CONTROL FAILED — the generator does not reproduce the committed file:")
            for b in bad:
                print(f"  🔴 {b}")
            return 1
        n = len(doc["per_element"])
        print(f"CONTROL REPRODUCED: all {len(DERIVED)} derived fields on all {n} rows "
              f"match the committed artifact at its own delta_p "
              f"(vmic={doc['delta_p_solar']['vmic']}).")
        return 0

    _, deltas = uncertainty_stack.params_and_deltas("solar")
    print(f"delta_p READ FROM pipeline/uncertainty_stack.py: {deltas}")
    was = doc["delta_p_solar"]["vmic"]
    now = deltas["vturb_kms"]
    got = recompute(doc, deltas)

    fe = next(r for r in got["per_element"] if r["element"] == "Fe" and r["ion"] == "I")
    fe_old = next(r for r in doc["per_element"] if r["element"] == "Fe" and r["ion"] == "I")
    print(f"  delta_xi   {was} -> {now} km/s")
    print(f"  Fe I sigma_B_vmic   {fe_old['sigma_B_vmic']} -> {fe['sigma_B_vmic']} dex")
    print(f"  Fe I sigma_reported {fe_old['sigma_reported']} -> {fe['sigma_reported']} dex")

    v = fe["sigma_B_vmic"]
    verdict = "FAILS" if v > SOLAR_GATE_DEX else "inside"
    print(f"\n  GATE: the xi term ALONE is {v:.4f} dex vs the {SOLAR_GATE_DEX} dex solar "
          f"gate -> {verdict}")
    if v > SOLAR_GATE_DEX:
        print("  ⚠️ REPORTED, NOT TUNED. RYA-1093 2E forbids adopting the 0.0588 formal "
              "error because it would pass the gate (RYA-161). This is a FLAG on the "
              "product, never a block; the adopt-or-relax call is Ryan's.")

    if a.dry_run:
        print("\n(--dry-run: nothing written)")
        return 0
    path.write_text(_dump(got))
    print(f"\nwrote {path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
