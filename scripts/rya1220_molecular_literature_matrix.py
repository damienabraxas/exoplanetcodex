#!/usr/bin/env python3
"""
RYA-1220 Work Package A -- the molecular CNO literature matrix.

Rolls the ingested Amarsi et al. 2021 Table 2 (408 used lines) up per species, and
cross-matches our own live molecular indicator windows against the reference used-line
ranges.

The point of the rollup is the MODEL-FORM TERM. Table 2 publishes each line's abundance
under five model atmospheres, so 3D-minus-MARCS is a measurement, not an assumption --
and it comes out molecule-specific, +0.035 (CH) to -0.140 (CO). A single uniform
molecular 3D offset would be wrong by up to 0.18 dex.

The point of the cross-match is BAND IDENTITY. Two of our four molecular N diagnostics
sit on bands this analysis does not use at all, and they carry our worst pathologies.
Coverage is compared on declared windows; a wavelength overlap is not a transition match,
which is why the per-line rows stay CROSSMATCH_REVIEW.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

import pandas as pd

REF = "data/reference/amarsi2021_cno/derived/amarsi2021_cno_molecular_lines.csv"
INV = "data/reference/nitrogen_method_rya1220/indicator_inventory.json"
OUT = "data/reference/molecular_cno_literature_rya1220"
TOL = 0.02  # fractional slack when asking "is our window inside their used range"


def by_species(ref: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for mol, g in ref.groupby("species"):
        rows.append({
            "molecule": mol,
            "system": g.system.mode()[0],
            "n_lines": len(g),
            "element": g.element_parameter.mode()[0].replace("logeps", ""),
            "lambda_vac_nm_min": round(g.wavelength_vac_nm.min(), 1),
            "lambda_vac_nm_max": round(g.wavelength_vac_nm.max(), 1),
            "bands": ";".join(sorted(g.source_band.unique())),
            "A_3D": round(g.abundance_3d.median(), 3),
            "A_mean3D": round(g.abundance_mean_3d.median(), 3),
            "A_MARCS_1D": round(g.abundance_marcs.median(), 3),
            "A_ATMO": round(g.abundance_atmo.median(), 3),
            "A_HM74": round(g.abundance_hm74.median(), 3),
            "model_form_3D_minus_MARCS": round(g.abundance_3d.median() - g.abundance_marcs.median(), 3),
            "granulation_3D_minus_mean3D": round(g.abundance_3d.median() - g.abundance_mean_3d.median(), 3),
            "line_to_line_sd_3D": round(g.abundance_3d.std(), 3),
        })
    return pd.DataFrame(rows).sort_values("molecule").reset_index(drop=True)


def _molecule_of(key: str) -> str:
    for m in ("CN", "NH", "CO", "OH", "CH", "C2"):
        if m in key:
            return m
    return "?"


def cross_match(inv: dict, species: pd.DataFrame) -> pd.DataFrame:
    cov = {r.molecule: (r.lambda_vac_nm_min * 10, r.lambda_vac_nm_max * 10)
           for r in species.itertuples()}
    ours: dict[str, dict] = {}
    for it in inv["indicators"]:
        if it.get("kind") != "molecular_band":
            continue
        w = it["windows_air_A"]
        lo, hi = min(x[0] for x in w), max(x[1] for x in w)
        cur = ours.setdefault(it["key"], {"element": it["element"], "lo": lo, "hi": hi,
                                          "n_windows": len(w), "holdings": set()})
        cur["lo"] = min(cur["lo"], lo)
        cur["hi"] = max(cur["hi"], hi)
        cur["holdings"].add(it["holding"])

    rows = []
    for key, v in ours.items():
        mol = _molecule_of(key)
        lo, hi = cov.get(mol, (None, None))
        inside = bool(lo is not None and v["lo"] >= lo * (1 - TOL) and v["hi"] <= hi * (1 + TOL))
        rows.append({
            "our_key": key, "molecule": mol, "element": v["element"],
            "our_window_air_A": f"{v['lo']:.0f}-{v['hi']:.0f}",
            "n_windows": v["n_windows"],
            "amarsi2021_used_range_A": f"{lo:.0f}-{hi:.0f}" if lo else "n/a",
            "in_amarsi2021_used_set": "YES" if inside else "NO",
            "holdings": ";".join(sorted(v["holdings"])),
        })
    return pd.DataFrame(rows).sort_values(["molecule", "our_key"]).reset_index(drop=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", default=REF)
    ap.add_argument("--inventory", default=INV)
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args()

    ref_path = pathlib.Path(args.ref)
    if not ref_path.exists():
        print(f"reference table absent: {ref_path}", file=sys.stderr)
        return 2
    ref = pd.read_csv(ref_path)
    ref = ref[ref.use_status == "USED"]

    species = by_species(ref)
    match = cross_match(json.loads(pathlib.Path(args.inventory).read_text()), species)

    out = pathlib.Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    species.to_csv(out / "amarsi2021_molecular_reference_by_species.csv", index=False)
    match.to_csv(out / "our_diagnostics_vs_amarsi2021.csv", index=False)

    print("RYA-1220 molecular CNO literature matrix")
    print(f"  reference lines used: {len(ref)}\n")
    print(species.to_string(index=False))
    print()
    print(match.to_string(index=False))
    off = match[match.in_amarsi2021_used_set == "NO"]
    print(f"\n  diagnostics OUTSIDE the reference used-line set: {list(off.our_key)}")
    print(f"  -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
