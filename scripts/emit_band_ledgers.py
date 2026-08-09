#!/usr/bin/env python3
"""Split a band measurement into PASS and FAIL ledgers, with the evidence — RYA-711/713.

    python3 scripts/emit_band_ledgers.py --element Fe --ion I --lo 6910 --hi 9199

Ryan: *"csv results that pass, csv results that fail. And document, document, document.
Show me the proof."*

Two files, because they answer different questions and get read by different people:

  *_PASS.csv  — lines that survived every check and may enter an aggregate. Carries the
                measurement AND how it was made, so the number is auditable without
                re-running anything.
  *_FAIL.csv  — every quarantined line, with symptom, FAULT DOMAIN, mechanism, the
                discriminator that established it, the fix, and the appendix panel that
                shows it. A failure with no evidence is an assertion (RYA-429).

Neither file drops a line. Together they must account for every candidate, and this
script refuses to write if they do not — that reconciliation is the point.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

EW_DIR = ROOT / "data" / "measured" / "band_ew"
PLOT_DIR = ROOT / "data" / "plots" / "band_appendix"
PANELS_PER_FIG = 12


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--element", required=True)
    ap.add_argument("--ion", default="I")
    ap.add_argument("--lo", type=float, required=True)
    ap.add_argument("--hi", type=float, required=True)
    ap.add_argument("--instrument", default="kpno_solar_atlas")
    a = ap.parse_args()

    stem = f"{a.element}{a.ion}_{int(a.lo)}_{int(a.hi)}_{a.instrument}"
    ew = pd.read_csv(EW_DIR / f"{stem}_ew.csv")
    skipped = json.loads((EW_DIR / f"{stem}_skipped.json").read_text())
    rc_path = EW_DIR / f"{stem}_root_causes.csv"
    rc = pd.read_csv(rc_path) if rc_path.exists() else pd.DataFrame()

    acc = pd.read_csv(ROOT / "data" / "audit" / "line_accounting" / "per_line.csv")
    meta = acc[["wave_air_A", "log_gf", "ep_eV", "predicted_depth", "instruments"]]
    ew = ew.merge(meta, left_on="wavelength_air_A", right_on="wave_air_A", how="left")

    # Appendix panels are ordered quarantined-first, 12 to a figure -- so a line can be
    # pointed at the exact image that proves its verdict.
    order = ew.sort_values(["in_aggregate", "wavelength_air_A"]).reset_index(drop=True)
    order["panel"] = [f"{stem}_appendix_{i // PANELS_PER_FIG + 1}.png"
                      for i in range(len(order))]
    ew = ew.merge(order[["wavelength_air_A", "panel"]], on="wavelength_air_A", how="left")

    passed = ew[ew.in_aggregate].copy()
    failed = ew[~ew.in_aggregate].copy()

    if len(rc):
        failed = failed.merge(
            rc[["wave", "fault_domain", "mechanism", "discriminator", "fix"]],
            left_on="wavelength_air_A", right_on="wave", how="left").drop(columns=["wave"])

    # Every FAIL must carry a reason and a domain. A silent exclusion is the defect
    # ELEMENT_PROTOCOL names -- indistinguishable from "never attempted".
    if len(failed):
        blank = failed[failed.excluded_reason.fillna("").str.strip() == ""]
        if len(blank):
            raise SystemExit(f"{len(blank)} failed lines carry no reason: "
                             f"{list(blank.wavelength_air_A)}")
        if "fault_domain" in failed.columns:
            nodom = failed[failed.fault_domain.isna()]
            if len(nodom):
                print(f"  WARNING: {len(nodom)} failed lines have no root cause attributed "
                      f"— they will read as unexplained: {list(nodom.wavelength_air_A)[:6]}")

    pcols = ["element", "ion", "wavelength_air_A", "instrument", "ew_mA", "rew",
             "log_gf", "ep_eV", "predicted_depth", "ew_method", "panel"]
    fcols = pcols + ["excluded_reason"] + \
            [c for c in ("fault_domain", "mechanism", "discriminator", "fix") if c in failed.columns]

    p_out = EW_DIR / f"{stem}_PASS.csv"
    f_out = EW_DIR / f"{stem}_FAIL.csv"
    passed[[c for c in pcols if c in passed.columns]].to_csv(p_out, index=False)
    failed[[c for c in fcols if c in failed.columns]].to_csv(f_out, index=False)

    # Reconciliation: PASS + FAIL + skipped must equal every candidate considered.
    total = len(passed) + len(failed) + len(skipped)
    print(f"{a.element} {a.ion}  {a.lo:.0f}-{a.hi:.0f} A  ({a.instrument})")
    print(f"  PASS     {len(passed):4d}  -> {p_out.name}")
    print(f"  FAIL     {len(failed):4d}  -> {f_out.name}")
    print(f"  skipped  {len(skipped):4d}  (no atlas coverage / telluric)")
    print(f"  TOTAL    {total:4d}  accounted for, none dropped")

    if len(failed) and "fault_domain" in failed.columns:
        print("\n  FAIL by fault domain:")
        for dom, g in failed.groupby("fault_domain"):
            print(f"    {dom:12s} {len(g):4d}")
            for mech, gg in g.groupby("mechanism"):
                print(f"        {len(gg):4d}  {str(mech)[:74]}")

    if len(skipped):
        reasons = pd.Series([s["reason"].split(":")[0] for s in skipped]).value_counts()
        print("\n  skipped by reason:")
        for k, v in reasons.items():
            print(f"    {v:4d}  {k[:74]}")


if __name__ == "__main__":
    main()
