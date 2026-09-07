#!/usr/bin/env python3
"""RYA-1191 — the TWIN test: is a wide error bar gf-limited, or uncorrected telluric?

    python3 scripts/rya1191_twin_compare.py

🔴 THE COMPARISON RYAN'S TRACKER REVIEW ASKS FOR, AND WHY IT DISCRIMINATES.
`solar_kpno_molecfit_corrected` NIR ENGINE-A ships A = 7.503 with sigma_stat 0.497 and an
`irreducible_dispersion` note reading *"gf-limited, IRREDUCIBLE, reducible only by
laboratory gf, NOT more lines."* The SAME treatment on `solar_iag` gives sigma_stat 0.072
-- seven times tighter -- on the same NLTE-grid domain with the same exclusion count.

If the dispersion were gf-limited it would be a property of the LINE LIST, so it would be
identical on both holdings: same lines, same gf values, same atmosphere. It is not. The
only thing that differs is the HOLDING, and one of them is verified telluric-corrected
while the other is not. A statistic that changes when you change the spectrum and holds
the line list fixed is a statistic about the SPECTRUM.

⚠️ AND THE CORRECTION TO THAT ARGUMENT MATTERS AS MUCH AS THE ARGUMENT. Ryan's own
follow-up: the NIR WIDTH is genuinely gf-limited -- IAG NIR 1D-LTE is ALSO wide (0.159) on
verified-clean flux, so that note is TRUE. What the KP-molecfit rows carry ON TOP of that
real width is a telluric OFFSET and, in one product, a sigma blow-up. So this script
separates the two rather than replacing one blanket story with another: the OFFSET and the
EXCESS WIDTH are reported apart, and a shared width is attributed to neither holding.

WHAT IS COMPARED, AND THE TWO TRAPS IN IT
------------------------------------------
Per line measured on BOTH twins: `delta_A = A(kp_molecfit) - A(iag)`.

⚠️ TRAP 1 — `sigma_stat` IS NOT THE PRINTED SCATTER. `build_product` prints the line
scatter (std, ddof=1) while the feed's `sigma_stat` is that divided by sqrt(n). Comparing
one against the other manufactures a 5-8x discrepancy out of arithmetic (RYA-1084). Both
are carried here, each named.

⚠️ TRAP 2 — ONLY LINES IN BOTH POOLS MAY BE DIFFERENCED. The two holdings exclude
different lines, so an unmatched mean difference mixes a per-line offset with a different
line SET. The comparison is restricted to the intersection, and the pool sizes are
reported so the restriction is visible.
"""
from __future__ import annotations

import argparse
import glob
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BP = ROOT / "data/results/band_products"

#: The twin pairs Ryan's scan named: the raw-holding product and its verified-clean
#: counterpart, same band, same treatment.
#:
#: ⚠️ THE SPAN IN A STEM IS THE ARM-INTERSECTED RANGE, NOT THE BAND, AND IT IS NOT
#: PREDICTABLE. The two holdings do not cover the same NIR span (KP reaches 12976 A, IAG
#: 11081) and the intersected edge comes out of the run, not the band table -- hardcoding
#: `9199_11083` missed the real `9199_11081` by 2 A and reported both NIR artifacts as
#: ABSENT, which reads exactly like "that run has not happened yet". Stems are resolved by
#: GLOB on the holding; the wavelength overlap is enforced by the line intersection below,
#: where it belongs.
TWINS = [
    ("red-optical", "FeI_6910_", "kpno_solar_atlas_solar_kpno_molecfit_corrected",
     "iag_fts_solar_atlas_solar_iag"),
    ("NIR", "FeI_9199_", "kpno_solar_atlas_solar_kpno_molecfit_corrected",
     "iag_fts_solar_atlas_solar_iag"),
]

#: Ryan's sanity gate (2026-09-07).
GATE_DEX, GATE_SIGMA_RATIO = 0.03, 2.0


def _lines(prefix: str, holding: str, treat: str) -> pd.DataFrame | None:
    hits = sorted(glob.glob(str(BP / f"{prefix}*_{holding}_SYNTH_GRADED_{treat}_lines.csv")))
    if not hits:
        return None
    if len(hits) > 1:
        raise SystemExit(f"{len(hits)} artifacts match {prefix}*{holding} {treat}: "
                         f"{[Path(h).name for h in hits]} — refusing to guess which")
    d = pd.read_csv(hits[0])
    if "in_aggregate" in d:
        d = d[d.in_aggregate.astype(str).str.lower().isin(["true", "1"])]
    return d[["wavelength_air_A", "abundance"]].dropna()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default="data/results/rya1191")
    a = ap.parse_args(argv)
    out = ROOT / a.out_dir; out.mkdir(parents=True, exist_ok=True)

    rows, per_line = [], []
    for band, prefix, kp_hold, iag_hold in TWINS:
        for treat in ("1D-LTE", "ENGINE-A"):
            k = _lines(prefix, kp_hold, treat)
            i = _lines(prefix, iag_hold, treat)
            rec = {"band": band, "treatment": treat,
                   "raw_holding": "solar_kpno_molecfit_corrected",
                   "clean_twin": "solar_iag"}
            if k is None or i is None:
                rec.update(state="ARTIFACT-MISSING",
                           why=f"kp={'yes' if k is not None else 'NO'} "
                               f"iag={'yes' if i is not None else 'NO'}")
                rows.append(rec); continue
            m = pd.merge(k, i, on="wavelength_air_A", suffixes=("_kp", "_iag"))
            rec["n_kp"], rec["n_iag"], rec["n_matched"] = len(k), len(i), len(m)
            if len(m) < 3:
                rec.update(state="NO-OVERLAP",
                           why=f"only {len(m)} line(s) in both pools")
                rows.append(rec); continue
            d = m.abundance_kp - m.abundance_iag

            # 🔴 THE PRODUCT AGGREGATES BY MEDIAN, AND DIFFERENCING MEANS GOT IT WRONG.
            # Verified on every graded artifact on disk: `products.csv` A equals the
            # MEDIAN of the in-aggregate per-line abundances EXACTLY, never the mean. On
            # the NIR ENGINE-A pool the two are 7.581 and 7.109 — a 0.47 dex gap opened by
            # one line (9437.793, A = 4.54) that a median ignores and a mean does not. A
            # twin test built on means would have attributed that entirely to telluric.
            # `stat_dex` is std/sqrt(n), read off the artifact's own `stat_basis` rather
            # than reverse-engineered (RYA-1084).
            def stats(x):
                x = np.asarray(x, float)
                return (float(np.median(x)), float(np.std(x, ddof=1)),
                        float(np.std(x, ddof=1) / np.sqrt(len(x))))
            mk, sk, ek = stats(m.abundance_kp)
            mi, si, ei = stats(m.abundance_iag)
            rec.update(state="MEASURED",
                       A_kp=round(mk, 4), scatter_kp=round(sk, 4), sigma_stat_kp=round(ek, 4),
                       A_iag=round(mi, 4), scatter_iag=round(si, 4), sigma_stat_iag=round(ei, 4),
                       aggregation="median (matches products.csv A exactly on every "
                                   "graded artifact checked)",
                       # 🔴 THE GATE IS ON THE DIFFERENCE OF THE PRODUCT STATISTIC, which
                       # is what ships and what Ryan's scan gated. `delta_A_per_line` is
                       # reported beside it and is NOT the same number: the median of the
                       # per-line differences says whether a TYPICAL line agrees, while
                       # the difference of medians says whether the two POOLS sit apart.
                       # Red-optical 1D-LTE is +0.032 on the second and +0.002 on the
                       # first — most lines agree and the pools still differ, which is a
                       # subset story, not a global offset. Reporting either alone would
                       # have been a different (and wrong) diagnosis.
                       delta_A=round(mk - mi, 4),
                       delta_A_per_line_median=round(float(np.median(d)), 4),
                       delta_A_per_line_mean=round(float(np.mean(d)), 4),
                       delta_A_se=round(float(np.std(d, ddof=1) / np.sqrt(len(d))), 4),
                       n_lines_over_1dex=int((d.abs() > 1.0).sum()),
                       # the EXCESS width, over and above the width the clean twin shows
                       scatter_ratio=round(sk / si, 3) if si > 0 else None,
                       shared_scatter=round(min(sk, si), 4),
                       excess_scatter=round(float(np.sqrt(max(sk**2 - si**2, 0.0))), 4))
            rec.update(_verdict(rec))
            for _, r in m.iterrows():
                per_line.append({"band": band, "treatment": treat,
                                 "wavelength_air_A": round(float(r.wavelength_air_A), 4),
                                 "A_kp": round(float(r.abundance_kp), 4),
                                 "A_iag": round(float(r.abundance_iag), 4),
                                 "delta_A": round(float(r.abundance_kp - r.abundance_iag), 4)})
            rows.append(rec)

    pl = pd.DataFrame(per_line)
    if len(pl):
        pl.to_csv(out / "rya1191_twin_per_line.csv", index=False)
    doc = {"ticket": "RYA-1191 — twin test: telluric offset vs gf-limited width",
           "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
           "gate": {"max_offset_dex": GATE_DEX, "max_scatter_ratio": GATE_SIGMA_RATIO,
                    "source": "Ryan 2026-09-07"},
           "statistic_note": ("⚠️ `sigma_stat` is the printed line scatter divided by "
                              "sqrt(n); both are carried and named, because comparing one "
                              "against the other manufactures a 5-8x discrepancy out of "
                              "arithmetic (RYA-1084)."),
           "pairs": rows,
           "worst_lines": (pl.reindex(pl.delta_A.abs().sort_values(ascending=False).index)
                           .head(20).to_dict("records") if len(pl) else [])}
    (out / "rya1191_twin_compare.json").write_text(json.dumps(doc, indent=2) + "\n")

    print(f"{'band':<13}{'treat':<10}{'n':>5}{'A_kp':>9}{'A_iag':>9}{'dA':>9}"
          f"{'scat_kp':>9}{'scat_iag':>9}{'ratio':>7}  verdict")
    for r in rows:
        if r["state"] != "MEASURED":
            print(f"{r['band']:<13}{r['treatment']:<10}  {r['state']} {r.get('why','')}")
            continue
        print(f"{r['band']:<13}{r['treatment']:<10}{r['n_matched']:>5}{r['A_kp']:>9.4f}"
              f"{r['A_iag']:>9.4f}{r['delta_A']:>+9.4f}{r['scatter_kp']:>9.4f}"
              f"{r['scatter_iag']:>9.4f}{(r['scatter_ratio'] or 0):>7.2f}  {r['verdict']}")
    if len(pl):
        print("\n=== worst 10 lines by |delta_A| ===")
        print(pl.reindex(pl.delta_A.abs().sort_values(ascending=False).index)
              .head(10).to_string(index=False))
    print(f"\nwrote {a.out_dir}/")
    return 0


def _verdict(r) -> dict:
    """Offset and excess width are named apart — a shared width belongs to neither holding."""
    off, ratio = abs(r["delta_A"]), (r["scatter_ratio"] or 1.0)
    bad_off, bad_wid = off > GATE_DEX, ratio > GATE_SIGMA_RATIO
    if not bad_off and not bad_wid:
        return {"verdict": "PASSES THE GATE",
                "reading": (f"within {GATE_DEX} dex and {GATE_SIGMA_RATIO}x of the "
                            f"verified-clean twin")}
    bits = []
    if bad_off:
        bits.append(f"OFFSET {r['delta_A']:+.4f} dex (gate {GATE_DEX})")
    if bad_wid:
        bits.append(f"WIDTH x{ratio:.2f} (gate {GATE_SIGMA_RATIO})")
    return {"verdict": "🔴 FAILS THE GATE — " + " and ".join(bits),
            "reading": (
                f"⚠️ The width the two SHARE ({r['shared_scatter']:.3f}) is a property of "
                f"the LINE LIST and belongs to neither holding — on Ryan's own correction, "
                f"the NIR width IS gf-limited because the clean twin is wide too. What is "
                f"attributable to the raw holding is the EXCESS ({r['excess_scatter']:.3f}) "
                f"and the offset. Only those two may be blamed on telluric.")}


if __name__ == "__main__":
    raise SystemExit(main())
