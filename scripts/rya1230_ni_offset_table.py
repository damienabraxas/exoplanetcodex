"""RYA-1230 -- the per-line decomposition table, from `rya1230_ni_offset_rca.py`'s ladder JSON.

Apportions each line's (S0 - 7.813) into window / molecules / continuum / model / gf and a
residual, using the LADDER order, and prints the same table as ASCII for the EOS. The
"alone" lever for each rung is carried beside it so order-dependence is visible.
Read-only: no synthesis is run here.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
A21 = 7.813
SHORT = {"solar_iag": "IAG", "solar_kpno_kurucz2005_corrected": "KP-k05",
         "solar_kpno_molecfit_corrected": "KP-mf"}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ladder", type=Path, nargs="+",
                    default=sorted((REPO / "data/output/rya1230").glob("ni_offset_ladder*.json")))
    ap.add_argument("--out", type=Path, default=REPO / "data/output/rya1230/ni_offset_table.csv")
    args = ap.parse_args()
    rows = [r for f in args.ladder for r in json.loads(f.read_text())["rows"]]
    out = []
    for r in rows:
        lad, alone = r["ladder"], r["alone"]
        rec = dict(
            holding=SHORT.get(r["holding"], r["holding"]), line_A=r["wavelength_air_A"],
            S0_prod=r["S0_PROD"]["A"], total_vs_A21=r["S0_PROD"]["A"] - A21,
            window=lad["window"], molecules=lad["molecules"], continuum=lad["continuum"],
            model_marcs=lad["marcs"], gf=r["gf_term"],
            final_S4=r["S4_MARCS"]["A"], residual=r["S4_MARCS"]["A"] - A21 - r["gf_term"],
            window_alone=alone["window"], molecules_alone=alone["molecules"],
            continuum_alone=alone["continuum"], marcs_alone=alone["marcs"],
            continuum_p95=r["continuum_p95"], continuum_model_ratio=r["continuum_model_ratio"],
            continuum_bracket_lo=r["continuum_bracket_dex"][0],
            continuum_bracket_hi=r["continuum_bracket_dex"][1],
            dA_per_0p1pct_cont=r["dA_per_plus0p1pct_continuum"],
            red_chi2_S0=r["S0_PROD"]["red_chi2"], red_chi2_S4=r["S4_MARCS"]["red_chi2"],
            model_molecular_fraction=r["blend_at_fit"]["molecular_fraction"],
            model_ni_fraction=r["blend_at_fit"]["ni_fraction"],
            model_feature_ew_mA=r["blend_at_fit"]["feature_ew_mA"],
            ni_ew_model_at_7813_mA=r["ni_ew_model_at_A21_mA"],
            ni_ew_model_at_fit_mA=r["ni_ew_model_at_fit_mA"],
            amarsi_ni_ew_I_mA=r["amarsi2020_ni_ew_intensity_mA"],
            amarsi_feature_ew_I_mA=r["amarsi2020_feature_ew_intensity_mA"],
            our_loggf=r["our_loggf"], amarsi_loggf=r["published_loggf_amarsi2020"])
        out.append(rec)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        for rec in out:
            w.writerow({k: (round(v, 4) if isinstance(v, float) else v) for k, v in rec.items()})
    hdr = ("hold    line      S0    tot   | window   mol    cont   MARCS   gf   | final  resid | "
           "cont(p95)  dA/0.1%  molfrac")
    print(hdr)
    print("-" * len(hdr))
    for r in out:
        print(f"{r['holding']:7s} {r['line_A']:8.2f} {r['S0_prod']:6.3f} {r['total_vs_A21']:+6.3f} | "
              f"{r['window']:+6.3f} {r['molecules']:+6.3f} {r['continuum']:+6.3f} "
              f"{r['model_marcs']:+6.3f} {r['gf']:+6.3f} | {r['final_S4']:6.3f} {r['residual']:+6.3f} | "
              f"{r['continuum_p95']:8.5f} {r['dA_per_0p1pct_cont']:+7.3f}  {r['model_molecular_fraction']:5.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
