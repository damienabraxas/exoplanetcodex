#!/usr/bin/env python3
"""RYA-1218/RYA-587: audit Si uncertainty evidence without admitting products."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
COMPONENTS = (
    "measurement", "transition_data", "stellar.teff", "stellar.logg", "stellar.xi",
    "stellar.metallicity", "continuum", "profile_ew", "pseudo_continuum", "telluric",
    "holding_instrument", "nlte", "model_atmosphere", "hfs_isotopes", "blends",
    "molecular_coupling",
)


def audit(root: Path) -> dict:
    products = []
    for path in sorted(root.rglob("*_products.csv")):
        frame = pd.read_csv(path)
        if frame.empty or not {"A", "n_lines", "stat_dex"}.issubset(frame.columns):
            continue
        row = frame.iloc[0]
        line_path = path.with_name(path.name.replace("_products.csv", "_lines.csv"))
        if not line_path.exists():
            candidates = sorted(path.parent.glob(path.stem.replace("_products", "") + "*_lines.csv"))
            line_path = candidates[0] if candidates else line_path
        lines = pd.read_csv(line_path) if line_path.exists() else pd.DataFrame()
        accepted = lines[lines.get("in_aggregate", pd.Series(dtype=bool)).astype(str).str.lower().eq("true")]
        values = pd.to_numeric(accepted.get("abundance", pd.Series(dtype=float)), errors="coerce").dropna()
        scatter = float(np.std(values, ddof=1)) if len(values) > 1 else None
        sem = scatter / math.sqrt(len(values)) if scatter is not None else None
        products.append({
            "artifact": str(path.relative_to(ROOT)),
            "treatment": str(row.get("treatment", "")),
            "A": None if pd.isna(row.get("A")) else float(row["A"]),
            "n_lines": int(row.get("n_lines", 0)),
            "n_excluded": int(row.get("n_excluded", 0)),
            "reported_stat_dex": None if pd.isna(row.get("stat_dex")) else float(row["stat_dex"]),
            "reported_syst_dex": None if pd.isna(row.get("syst_dex")) else float(row["syst_dex"]),
            "line_scatter_dex": scatter,
            "scatter_sem_dex": sem,
            "measurement_state": "HOLD_estimator_or_independence" if len(values) > 1 else "HOLD_single_or_empty_pool",
            "accepted_line_count": int(len(values)),
        })
    return {
        "ticket": "RYA-1218",
        "uncertainty_contract": "RYA-587",
        "scope": "Si I production VIS and CRIRES+ synthesis products",
        "products": products,
        "components": [{"name": c, "state": "HOLD", "sigma_dex": None,
                        "reason": "Exact-pool evidence or covariance not supplied."}
                       for c in COMPONENTS],
        "admitted_products": 0,
        "disposition": "No total uncertainty is reported; numeric scatter/SE fields are diagnostics only.",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=ROOT / "data/results/rya1218")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    if args.out.exists():
        raise SystemExit(f"refusing to overwrite existing uncertainty audit: {args.out}")
    args.out.mkdir(parents=True)
    doc = audit(args.root)
    (args.out / "si_uncertainty_audit.json").write_text(json.dumps(doc, indent=2) + "\n")
    pd.DataFrame(doc["products"]).to_csv(args.out / "si_uncertainty_products.csv", index=False)
    print(json.dumps({"products": len(doc["products"]), "admitted_products": 0}, indent=2))


if __name__ == "__main__":
    main()
