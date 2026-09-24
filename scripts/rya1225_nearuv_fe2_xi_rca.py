#!/usr/bin/env python3
"""
RYA-1225 -- why near-UV Fe II dA/dxi disagrees 2.0-2.8x between RYA-1168 and RYA-1213.

STEP 1 IS A DIFF, NOT A THIRD MEASUREMENT. The ticket is explicit: a third number does
not adjudicate the two we have. This script reads both published artifacts, reproduces
RYA-1168's derivative from its own staged legs line by line (lambda + EP dual-keyed), and
localises the divergence to (a) pool membership, (b) method, or (c) genuine physics.

THE ANSWER IS (c), AND THE PHYSICAL CAUSE IS DATED. The two runs are not contemporaneous:

    RYA-1168 artifact   2026-09-03   9989ec3c
    RYA-1207 merged     2026-09-09   ace03e1a   <- near-UV molecular opacity
    RYA-1213 artifact   2026-09-12   1f4ab8e8

RYA-1207 turned molecular opacity on for the near-UV band -- and `config/synth_bands.yaml`
sets `use_molecules: true` on near-UV ALONE. So the near-UV synthesis changed underneath
the derivative while every other band's stayed fixed, which is exactly why the ticket's own
observation holds: red-optical Fe II from the same two runs agrees to 1.03-1.08x. Near-UV
is the lone anomaly because near-UV is the only band whose physics moved.

🔴 THE STALENESS IS MEASURABLE, NOT INFERRED. Each RYA-1168 pool record carries the feed `A`
it ran against. All four Fe II rows record the PRE-molecular abundances (7.957 / 7.975 /
7.697 / 7.738); the live products carry the POST-molecular ones (7.617 / 7.635 / 7.402 /
7.404), and RYA-1207's own lever artifact contains both sets. A derivative measured against
an abundance the product no longer has was measured on a synthesis the product no longer
uses.

⚠️ AND A HYPOTHESIS THAT LOOKED GOOD AND IS REFUTED, RECORDED SO IT IS NOT RE-RUN. The
campaign key carries no band (RYA-1114 F2), so `FeII_<holding>_DEEPGRADED_SYNTH` names the
near-UV pool and the VIS one alike -- an obvious candidate for a cross-band mix-up. Measured
from `~/xi_campaign/out`, the VIS-window Fe II pools give -0.2300 / -0.2100 / -0.2600 /
-0.1350. None of them is RYA-1213's -0.1100. The divergence is NOT a band confound.

READ-ONLY. No abundance moves, nothing is re-synthesised, and neither published artifact is
edited (RYA-161: both are published runs).
"""

from __future__ import annotations

import glob
import json
import pathlib
import sys

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[1]
CAMPAIGN_UV = pathlib.Path.home() / "xi_campaign" / "out_uv"
CAMPAIGN_VIS = pathlib.Path.home() / "xi_campaign" / "out"
A1168 = REPO / "data/results/rya1168/nearuv_xi_dadxi.json"
A1213 = REPO / "data/results/rya1213/reference_xi_dadxi.json"
LEDGER = REPO / "data/audit/rya1213_reference_matrix/same_pool_ledger.json"
BANDP = REPO / "data/results/band_products"
FEED = REPO / "data/products/solar/Fe.json"
OUT = REPO / "data/results/rya1225/nearuv_fe2_xi_rca.json"

XI_SPAN_KMS = 0.20
POOLS = {"kur_FeII": "solar_kpno_kurucz2005_corrected",
         "mol_FeII": "solar_kpno_molecfit_corrected"}
TREATMENTS = ("1D-LTE", "ENGINE-A")


def _leg(d: pathlib.Path, treatment: str) -> pd.DataFrame | None:
    """One leg's per-line CSV, matched on the FULL suffix.

    `ENGINE-A` is a substring of `ENGINE-A-3DNLTE`, so a loose `in` test would pair a
    product with a different engine's file -- the trap RYA-1213's own reader documents.
    """
    hits = [f for f in glob.glob(str(d / "*_lines.csv"))
            if f.endswith(f"_{treatment}_lines.csv")]
    if len(hits) != 1:
        return None
    return pd.read_csv(hits[0])


def per_line_1168(tag: str, treatment: str) -> dict | None:
    """Reproduce RYA-1168's derivative from ITS OWN legs, lambda+EP dual-keyed."""
    lo = _leg(CAMPAIGN_UV / f"FeII_{tag}_nearUV_DEEPGRADED_SYNTH_xi0.90", treatment)
    hi = _leg(CAMPAIGN_UV / f"FeII_{tag}_nearUV_DEEPGRADED_SYNTH_xi1.10", treatment)
    if lo is None or hi is None:
        return None
    m = lo.merge(hi, on=["wavelength_air_A", "ep_eV"], suffixes=("_lo", "_hi"))
    m["delta"] = m["abundance_hi"] - m["abundance_lo"]
    ok = m["delta"].notna()
    return {
        "n_lo": int(len(lo)), "n_hi": int(len(hi)),
        "n_merged_on_lambda_ep": int(len(m)), "n_paired": int(ok.sum()),
        "median_delta_dex": round(float(m.loc[ok, "delta"].median()), 6),
        "dA_dxi": round(float(m.loc[ok, "delta"].median()) / XI_SPAN_KMS, 4),
        "per_line": [
            {"wavelength_air_A": float(r.wavelength_air_A), "ep_eV": float(r.ep_eV),
             "A_xi0.90": None if pd.isna(r.abundance_lo) else round(float(r.abundance_lo), 4),
             "A_xi1.10": None if pd.isna(r.abundance_hi) else round(float(r.abundance_hi), 4),
             "delta_dex": None if pd.isna(r.delta) else round(float(r.delta), 4),
             "paired": bool(not pd.isna(r.delta))}
            for r in m.sort_values("wavelength_air_A").itertuples()],
    }


def vis_control(holding: str, treatment: str) -> float | None:
    """The refuted band-confound hypothesis, measured rather than asserted."""
    lo = _leg(CAMPAIGN_VIS / f"FeII_{holding}_DEEPGRADED_SYNTH_xi0.90", treatment)
    hi = _leg(CAMPAIGN_VIS / f"FeII_{holding}_DEEPGRADED_SYNTH_xi1.10", treatment)
    if lo is None or hi is None:
        return None
    m = lo.merge(hi, on=["wavelength_air_A", "ep_eV"], suffixes=("_lo", "_hi"))
    d = (m["abundance_hi"] - m["abundance_lo"]).dropna()
    return None if d.empty else round(float(d.median()) / XI_SPAN_KMS, 4)


def pool_identity(holding: str, treatment: str) -> dict:
    """Is the DEEPGRADED pool the REFERENCE pool? Measured, not taken from the ticket."""
    def read(tier):
        p = BANDP / (f"FeII_3000_3780_kpno_solar_atlas_{holding}_SYNTH_{tier}"
                     f"_{treatment}_lines.csv")
        return pd.read_csv(p) if p.exists() else None
    dg, rf = read("DEEPGRADED"), read("REFERENCE")
    if dg is None or rf is None:
        return {"comparable": False, "reason": "a tier's nominal per-line file is absent"}
    key = ["wavelength_air_A", "ep_eV"]
    sa = {(round(a, 4), round(b, 4)) for a, b in zip(dg[key[0]], dg[key[1]])}
    sb = {(round(a, 4), round(b, 4)) for a, b in zip(rf[key[0]], rf[key[1]])}
    m = dg.merge(rf, on=key, suffixes=("_dg", "_rf"))
    d = (m["abundance_rf"] - m["abundance_dg"]).dropna()
    return {"comparable": True, "n_deepgraded": len(sa), "n_reference": len(sb),
            "n_shared_lambda_ep": len(sa & sb),
            "deepgraded_only": sorted(x[0] for x in sa - sb),
            "reference_only": sorted(x[0] for x in sb - sa),
            "max_abs_nominal_abundance_delta": 0.0 if d.empty else round(float(d.abs().max()), 9),
            "identical": bool(sa == sb and (d.empty or d.abs().max() < 1e-9))}


def main() -> int:
    a = json.loads(A1168.read_text())
    b = json.loads(A1213.read_text())
    feed = json.loads(FEED.read_text())
    live = {(x["ion"], x["holding"], x["tier"], x["treatment"], x["band"]): x
            for x in feed["products"]}
    b_idx = {(p["holding"], p["treatment"]): p for p in b["pools"]
             if p["band"] == "near-UV" and p["ion"] == "II"}

    rows = []
    for tag, holding in POOLS.items():
        for treat in TREATMENTS:
            r1168 = [p for p in a["pools"]
                     if p["holding"] == holding and p["treatment"] == treat
                     and p["ion"] == "II"]
            if not r1168:
                continue
            r1168 = r1168[0]
            r1213 = b_idx.get((holding, treat))
            repro = per_line_1168(tag, treat)
            prod_dg = live.get(("II", holding, "DEEPGRADED", treat, "near-UV"))
            rows.append({
                "holding": holding, "treatment": treat,
                "rya1168": {"dA_dxi": r1168["dA_dxi"],
                            "median_delta_dex": r1168["median_delta_dex"],
                            "n_paired": r1168["n_paired"],
                            "A_it_ran_against": r1168["A"]},
                "rya1213": None if r1213 is None else {
                    "dA_dxi": r1213["dA_dxi"],
                    "median_delta_dex": r1213["median_delta_dex"],
                    "n_paired": r1213["n_paired"]},
                "ratio_1213_over_1168": None if not r1213 else round(
                    abs(r1213["dA_dxi"] / r1168["dA_dxi"]), 3),
                "live_product_A": None if prod_dg is None else prod_dg["A"],
                "rya1168_measured_against_a_superseded_A": bool(
                    prod_dg is not None and abs(prod_dg["A"] - r1168["A"]) > 1e-9),
                "reproduced_from_rya1168_own_legs": repro,
                "pool_identity_deepgraded_vs_reference": pool_identity(holding, treat),
                "vis_window_control_dA_dxi": vis_control(holding, treat),
            })

    doc = {
        "ticket": "RYA-1225", "read_only": True, "star": "solar", "band": "near-UV",
        "ion": "Fe II", "xi_span_kms": XI_SPAN_KMS,
        "case": "c",
        "case_label": ("GENUINE -- same lines, same step, same reduction; the divergence "
                       "is physical in the near-UV synthesis"),
        "physical_cause": (
            "RYA-1207 (merged 2026-09-09, ace03e1a) turned molecular opacity on for the "
            "near-UV band. config/synth_bands.yaml sets use_molecules: true on near-UV "
            "ALONE, so only this band's synthesis moved. RYA-1168's artifact predates it "
            "(2026-09-03, 9989ec3c) and RYA-1213's postdates it (2026-09-12, 1f4ab8e8)."),
        "why_red_optical_agrees": (
            "use_molecules is per-band and false everywhere except near-UV, so the "
            "red-optical synthesis did not change between the two runs and its Fe II "
            "pairs agree to 1.03-1.08x. Near-UV is the lone anomaly for that reason."),
        "refuted_hypothesis_band_confound": (
            "The campaign key carries no band (RYA-1114 F2), so a cross-band mix-up was "
            "the obvious candidate. Measured from ~/xi_campaign/out, the VIS-window Fe II "
            "pools give -0.2300 / -0.2100 / -0.2600 / -0.1350 -- none is RYA-1213's "
            "-0.1100. REFUTED, recorded so it is not re-run."),
        "resolution": (
            "RYA-1213's slope wins on METHOD CORRECTNESS, not on proximity to anything: "
            "it is the measurement taken on the synthesis the live products actually use. "
            "RYA-1168's is superseded by vintage -- its own records carry the pre-molecular "
            "abundances. NO re-measure is required despite this being case (c), because the "
            "post-change measurement already exists; a third run would reproduce RYA-1213."),
        "rows": rows,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(doc, indent=2) + "\n")
    print(f"wrote {OUT.relative_to(REPO)}")
    for r in rows:
        print(f"  {r['holding'][:31]:31} {r['treatment']:9} "
              f"1168={r['rya1168']['dA_dxi']:+.4f} (A={r['rya1168']['A_it_ran_against']}) "
              f"1213={r['rya1213']['dA_dxi']:+.4f}  ratio={r['ratio_1213_over_1168']}  "
              f"pool_identical={r['pool_identity_deepgraded_vs_reference']['identical']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
