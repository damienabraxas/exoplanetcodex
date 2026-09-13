#!/usr/bin/env python3
"""RYA-1214 — our CNO indicators against AGSS21 Table 3, indicator by indicator.

This is the campaign's answer sheet. AGSS21 publishes a per-INDICATOR abundance, not one
number per element, and the whole point of replicating it is to compare like with like: a
CH G-band value belongs beside AGSS21's CH row, not beside its C I row.

🔴 THE COMPARISON RULE, AND THE TRAP IT AVOIDS. AGSS21 publishes each indicator in TWO
columns, 3D-LTE and 3D-NLTE, and our 1D-LTE number is comparable to NEITHER without
saying which and why:

  * A FORBIDDEN line ([C I] 8727, [O I] 6300) is LTE-insensitive — AGSS21's own two
    columns are IDENTICAL for both — so 1D-LTE vs 3D-LTE isolates the ATMOSPHERE, and that
    is a real, interpretable difference.
  * A MOLECULAR band has no NLTE grid at all (`lte_molecular_band`), so AGSS21 publishes
    only a 3D-LTE value for it and the same comparison applies.
  * A PERMITTED atomic line is the one case where the two AGSS21 columns differ, and there
    our 1D-LTE is comparable to the 3D-LTE column while the 3D-NLTE column differs by
    physics we did not apply.

So every row below names WHICH AGSS21 column it is against and what the residual therefore
measures. Differencing against whichever column is closer would manufacture agreement
(RYA-161), and quoting one number per element would hide that the two nitrogen routes
disagree with each other by ~0.9 dex.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "audit" / "rya1214_cno_products"
AGSS21 = ROOT / "data" / "audit" / "rya1160_cno_nist_gf" / "agss21_table3_indicators.csv"
CNOSYNTH = ROOT / "data" / "audit" / "cno_synthesis"
BP = ROOT / "data" / "results" / "band_products"

#: Our diagnostic -> the AGSS21 Table 3 row it replicates, and the comparison basis.
#: `column` is which AGSS21 column the residual is taken against; `basis` says what the
#: residual therefore MEASURES. Declared, never inferred from which one is closer.
MAP = {
    "CH_Gband":  ("C", "CH (A-X)",   "lte_3D", "1D vs 3D atmosphere; molecular band, no NLTE grid either side"),
    "C2_Swan":   ("C", "C2 Swan",    "lte_3D", "1D vs 3D atmosphere; molecular band, no NLTE grid either side"),
    "CN_red":    ("N", "CN (dnu>=1)", "lte_3D", "1D vs 3D atmosphere; molecular band, no NLTE grid either side"),
    "OI_6300":   ("O", "[O i]",      "nlte_3D", "FORBIDDEN and LTE-insensitive — AGSS21's two columns are identical, so this isolates the atmosphere"),
    "CI_5052":   ("C", "C i",        "lte_3D", "permitted atomic; AGSS21's 3D-NLTE column differs by physics we did not apply"),
    "CI_5380":   ("C", "C i",        "lte_3D", "permitted atomic; AGSS21's 3D-NLTE column differs by physics we did not apply"),
    "NH_AX":     ("N", "NH (dnu=1)", "lte_3D", "near-UV molecular band — UPPER BOUND (near-UV opacity deficit, RYA-1204/1207)"),
    "OH_AX":     ("O", "OH (dnu=1)", "lte_3D", "near-UV molecular band — UPPER BOUND (near-UV opacity deficit, RYA-1204/1207)"),
}

#: The band-product route's per-line values for the two FORBIDDEN indicators, which
#: `cno_synthesis` measures jointly with their blend partner and the band route measures
#: as one line of a named set. Both are legitimate and they are NOT the same measurement,
#: so both are reported.
FORBIDDEN_LINES = {"[C I] 8727.14": (8727.139, "C", "[C i]"),
                   "[O I] 6300.30": (6300.304, "O", "[O i]")}


def _agss21() -> pd.DataFrame:
    d = pd.read_csv(AGSS21)
    d["key"] = d.element + "|" + d.indicator
    return d.set_index("key")


def _per_band() -> pd.DataFrame:
    frames = []
    for p in sorted(CNOSYNTH.glob("solar_*_cno_per_band.csv")):
        f = pd.read_csv(p)
        f["region"] = p.name.split("_")[1]
        frames.append(f)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _forbidden_from_band_products() -> list[dict]:
    rows = []
    for name, (wave, el, ind) in FORBIDDEN_LINES.items():
        for f in sorted(BP.glob(f"{el}I_*_lines.csv")):
            d = pd.read_csv(f)
            m = d[(d.wavelength_air_A - wave).abs() < 0.05]
            m = m[m.in_aggregate.astype(str).str.lower() == "true"]
            for _, r in m.iterrows():
                hold = f.name.split("_SYNTH")[0]
                for inst in ("kpno_solar_atlas", "harps", "iag_fts_solar_atlas"):
                    if f"_{inst}_" in hold:
                        hold = hold.split(f"_{inst}_", 1)[1]
                        break
                rows.append({"indicator_line": name, "element": el,
                             "agss21_indicator": ind, "holding": hold,
                             "A_X": float(r.abundance),
                             "red_chi2": round(float(r.red_chi2), 1)})
    return rows


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    ref = _agss21()
    pb = _per_band()
    rows = []
    for _, r in pb.iterrows():
        spec = MAP.get(r.key)
        if spec is None:
            continue
        el, ind, col, basis = spec
        k = f"{el}|{ind}"
        if k not in ref.index:
            rows.append({"our_key": r.key, "region": r.region, "element": el,
                         "agss21_indicator": ind, "verdict": "NO AGSS21 ROW"})
            continue
        a = ref.loc[k]
        pub = a[col]
        pub = float(pub) if pd.notna(pub) else None
        ours = float(r.A_X)
        rows.append({
            "our_key": r.key, "region": r.region, "element": el, "role": r.role,
            "agss21_indicator": ind,
            "agss21_column": col, "agss21_A": pub,
            "agss21_error_dex": float(a.error_dex),
            "agss21_line_count": int(a.agss21_line_count),
            "our_A": round(ours, 3),
            "our_sigma_fit": r.sigma_fit,
            "our_red_chi2": round(float(r.red_chi2), 1),
            "delta_ours_minus_agss21": (round(ours - pub, 3) if pub is not None else None),
            "comparison_basis": basis,
            "upper_bound_only": bool(r.region == "nearuv"),
        })
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "agss21_indicator_comparison.csv", index=False)

    fb = pd.DataFrame(_forbidden_from_band_products())
    if len(fb):
        fb["agss21_A"] = fb.agss21_indicator.map(
            lambda i: float(ref.loc[f"{fb.element.iloc[0]}|{i}", "nlte_3D"])
            if f"{fb.element.iloc[0]}|{i}" in ref.index else None)
        for idx, r in fb.iterrows():
            k = f"{r.element}|{r.agss21_indicator}"
            fb.at[idx, "agss21_A"] = (float(ref.loc[k, "nlte_3D"])
                                      if k in ref.index else None)
        fb["delta"] = (fb.A_X - fb.agss21_A).round(3)
        fb.to_csv(OUT / "agss21_forbidden_line_comparison.csv", index=False)

    pd.set_option("display.width", 240)
    print("=== RYA-1214 — our CNO indicators vs AGSS21 Table 3 ===")
    if len(d):
        print(d[["region", "our_key", "element", "role", "agss21_indicator",
                 "agss21_column", "agss21_A", "our_A", "delta_ours_minus_agss21",
                 "our_sigma_fit", "our_red_chi2", "upper_bound_only"]].to_string(index=False))
    print("\n=== the two FORBIDDEN indicators, per line, from the band-product route ===")
    if len(fb):
        print(fb[["indicator_line", "holding", "A_X", "agss21_A", "delta",
                  "red_chi2"]].to_string(index=False))

    (OUT / "agss21_indicator_comparison.prov.json").write_text(json.dumps({
        "ticket": "RYA-1214",
        "what": "per-INDICATOR comparison against AGSS21 Table 3",
        "rule": "every row names WHICH AGSS21 column it is against and what the residual "
                "measures. A forbidden line is LTE-insensitive (AGSS21's two columns are "
                "identical) so 1D-LTE vs 3D isolates the atmosphere; a molecular band has "
                "no NLTE grid either side; a permitted atomic line is compared to the "
                "3D-LTE column because the 3D-NLTE column differs by physics we did not "
                "apply. Differencing against whichever column is closer would manufacture "
                "agreement (RYA-161).",
        "n_rows": int(len(d)),
        "n_forbidden_line_rows": int(len(fb)),
        "upper_bound_rows": int(d.upper_bound_only.sum()) if len(d) else 0,
    }, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
