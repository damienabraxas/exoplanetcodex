#!/usr/bin/env python3
"""RYA-963 — register the telluric-corrected α Cen A CRIRES+ holding.

A corrected SIBLING row, following the `solar_harps_molecfit_corrected` precedent
(RYA-931): the uncorrected `alpha_cen_a_crires_plus` EXTRACTC IDPs stay the archival
base and are NOT modified, and the corrected products get their own holding_id.

`telluric_applied` is **derived, not asserted** — `pipeline.telluric_intake` reads it off
the products' own MTRANS extension, and refuses `applied` if that transmission is all
unity. So a frozen molecfit fit (every column at its prior, transmission ≡ 1.0) cannot
reach the registry as a correction even if this script were run on it.

    python3 scripts/rya963_register_corrected.py --products <dir>          # report
    python3 scripts/rya963_register_corrected.py --products <dir> --write  # write

⚠️ The registry is edited with `csv`, never round-tripped through pandas: pandas rewrites
quoting and line endings across every untouched row and buries the real change (RYA-786).
"""
from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import io
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pipeline.telluric_intake import from_many                      # noqa: E402

HOLDINGS = ROOT / "data" / "catalog" / "holdings_manifest_registry.csv"
HOLDING_ID = "alpha_cen_a_crires_plus_molecfit"
BASE_HOLDING = "alpha_cen_a_crires_plus"


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def build_notes(run: dict, products: list) -> str:
    g = run["gdas_gate"]
    slots = sorted(set(g["slots"].values()))
    frames = run["frames"]
    confirmed = [f for f in frames if f["star_id"] == "A"]
    gates = [f"{f['wlen_id']} {f['gate_before']:.4f}->{f['gate_after']:.4f}"
             for f in frames]
    pwv = [f"{f['wlen_id']} {f['h2o_col_mm']:.3f}" for f in frames
           if f.get("h2o_col_mm") is not None]
    return (
        f"RYA-963 telluric-corrected sibling of {BASE_HOLDING}; {BASE_HOLDING} itself is "
        f"UNMODIFIED. {len(frames)} CRIRES+ EXTRACTC IDPs from {g['night']} "
        f"(one per setting: {', '.join(f['wlen_id'] for f in frames)}), corrected with "
        f"ESO molecfit (molecfit_model on derived fit windows, then molecfit_calctrans "
        f"over every chip) against the REAL per-night Paranal GDAS profile "
        f"{'/'.join(slots)} - all {len(frames)} exposures fall in the same 3-hourly slot, "
        f"and no standard-atmosphere fallback is in play (the RYA-373 failure mode). "
        f"ONE atmosphere is fitted per exposure and evaluated over all chips: one "
        f"exposure saw one sky, so a per-chip independent fit would let the water column "
        f"disagree with itself inside a single 4-second exposure. "
        f"THE FIT IS PROVEN TO HAVE MOVED, not merely to have exited 0: molecfit reports "
        f"success and writes every product even when mpfit takes its status-4 exit with a "
        f"zero Jacobian, leaving each column at its prior with uncertainty exactly 0 and "
        f"transmission 1.0 everywhere - an uncorrected spectrum wearing a correction's "
        f"provenance. Two distinct causes produced exactly that here and both are now "
        f"asserted against: WLC_CONST is a FRACTION of the chip's half wavelength range, "
        f"so its -0.05 default is -3.4 A on a 134 A chip but -135 A on a whole 540 nm "
        f"setting handed over as one chip; and without MAP_REGIONS_TO_CHIP molecfit says "
        f"'Assuming that all regions are mapped to Chip 1' and puts every fit window on "
        f"the first chip. Fit windows are DERIVED per frame (the unsaturated, non-stellar "
        f"absorbed fraction per chip), not tabulated, and the stellar exclusion mask comes "
        f"from linelist_solar rather than hand-listed wavelengths. Molecules are split: "
        f"the MODEL carries every molecule with a band in the frame (or calctrans would "
        f"leave that absorber out entirely), while only those with a band inside a fit "
        f"window are FITTED. D1 residual gate, at telluric-dominated stellar-clean pixels, "
        f"before->after: {'; '.join(gates)}. Fitted PWV (mm): {'; '.join(pwv)}. "
        f"STAR ID: {len(confirmed)}/{len(frames)} frames confirmed alpha Cen A by the "
        f"RYA-423 rule run on the measured RV - which is only available AFTER correction, "
        f"since RYA-423's own CRIRES branch returns INDETERMINATE for want of a telluric-"
        f"corrected spectrum. telluric_applied is DERIVED by pipeline.telluric_intake from "
        f"each product's own non-unity MTRANS extension, not asserted here. "
        f"Evidence: data/audit/rya963_crires_telluric/; reproduce with "
        f"python3 -m pipeline.crires_stellar_telluric --set alpha_cen_a_crires."
    ).replace("\n", " ")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--products", required=True, help="directory of corrected products")
    ap.add_argument("--run-json", required=True, help="run_set record")
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()

    run = json.loads(Path(a.run_json).read_text())
    products = sorted(glob.glob(str(Path(a.products) / "*_telluric.fits")))
    if not products:
        raise SystemExit(f"no corrected products under {a.products}")

    # DERIVE the state from the products; never assert it.
    state, evidence = from_many(products)
    print(f"products: {len(products)}")
    for p in products:
        print(f"  {Path(p).name}  sha256={sha256(p)}")
    print(f"derived telluric_applied = {state}")
    for e in evidence[:2]:
        print(f"  evidence: {e}")
    if state != "applied":
        raise SystemExit(
            f"REFUSING to register: intake derives telluric_applied={state!r} from the "
            f"products themselves. A holding is not corrected because this script says "
            f"so — if the transmission is all unity the fit never moved.")

    manifest_dir = ROOT / "data" / "audit" / "rya963_crires_telluric"
    manifest = manifest_dir / "corrected_manifest.csv"
    manifest_dir.mkdir(parents=True, exist_ok=True)
    with open(manifest, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["product", "sha256", "base_frame", "wlen_id", "band", "date_obs",
                    "gdas_profile", "gdas_md5", "star_id", "star_id_contested",
                    "gate_before", "gate_after", "gate_passed"])
        by_name = {Path(f["product"]).name: f for f in run["frames"]}
        for p in products:
            f = by_name.get(Path(p).name)
            if f is None:
                continue
            # A star-id whose underlying rule is contested must not be recorded as a
            # bare letter. The products carry STARIDQ; the manifest -- which is what the
            # live tracker reads -- carried nothing, so a reader of the dashboard saw a
            # settled verdict where the header said "disputed".
            w.writerow([Path(p).name, sha256(p), f["frame"], f["wlen_id"], f["band"],
                        f["date_obs"], f["gdas"], f["gdas_md5"], f["star_id"],
                        bool(f.get("star_id_branch_contested")),
                        f"{f['gate_before']:.5f}", f"{f['gate_after']:.5f}",
                        f["gate_passed"]])
    print(f"[manifest] {manifest}")

    rows = list(csv.DictReader(open(HOLDINGS, newline="")))
    fields = list(rows[0].keys())
    row = {"holding_id": HOLDING_ID, "system_id": "alpha_cen_a",
           "instrument_id": "crires_plus",
           "manifest_path": str(manifest.relative_to(ROOT)),
           "evidence_state": "verified", "source_issue_ids": "RYA-963",
           "notes": build_notes(run, products), "telluric_applied": state}
    existing = [i for i, r in enumerate(rows) if r["holding_id"] == HOLDING_ID]
    if existing:
        rows[existing[0]] = row
        print(f"UPDATE row {HOLDING_ID}")
    else:
        rows.append(row)
        print(f"APPEND row {HOLDING_ID}")

    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=fields, lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    if a.write:
        HOLDINGS.write_text(buf.getvalue())
        print(f"[written] {HOLDINGS}")
    else:
        print("(report only; pass --write to edit the registry)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
