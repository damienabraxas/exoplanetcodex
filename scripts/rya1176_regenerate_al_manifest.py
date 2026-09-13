#!/usr/bin/env python3
"""Regenerate the RYA-1176 Al manifest with explicit provenance axes."""
from pathlib import Path
import hashlib, json
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/audit/rya1132_al_intake/al_line_manifest.csv"
OUT = ROOT / "data/audit/rya1176_al_manifest/al_line_manifest_v2.csv"
OUR_GRADED = {"alphys_I_6696.0150_0352", "alphys_I_6698.6730_0355", "alphys_I_11253.1890_0406", "alphys_I_13123.4160_0423", "alphys_I_13150.7530_0425"}
OUR_DEEP = {"alphys_I_2652.4750_0318", "alphys_I_2660.3860_0319", "alphys_II_2669.1550_0320", "alphys_I_3082.1530_0329", "alphys_I_3092.7100_0330", "alphys_I_3944.0060_0335", "alphys_I_3961.5200_0336"}

def main():
    d = pd.read_csv(SOURCE).fillna("")
    if len(d) != 505 or not d.canonical_line_id.is_unique:
        raise SystemExit("frozen 505-row denominator changed")
    d["gf_provenance"] = d["gf_source"]
    d["selection_state"] = d["measurement_suitability_status"]
    d["line_set"] = "our-all"
    d.loc[d.canonical_line_id.isin(OUR_GRADED), "line_set"] = "our-graded"
    d.loc[d.canonical_line_id.isin(OUR_DEEP), "line_set"] = "our-deep-graded"
    d["telluric_applied"] = "unknown"
    d["normalization_state"] = "unknown"
    d["observed_conditioning"] = "unknown"
    d["conditioning_source"] = "holdings_manifest_registry;product_observed_conditioning"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    d.to_csv(OUT, index=False)
    prov = {"ticket": "RYA-1176", "source": str(SOURCE.relative_to(ROOT)), "output": str(OUT.relative_to(ROOT)), "denominator": len(d), "line_set_counts": d.line_set.value_counts().to_dict(), "our_graded_ids": sorted(OUR_GRADED), "our_deep_graded_ids": sorted(OUR_DEEP), "conditioning_policy": "unknown_until_holding_join", "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest()}
    (OUT.parent / "provenance.json").write_text(json.dumps(prov, indent=2, sort_keys=True) + "\n")

if __name__ == "__main__":
    main()
