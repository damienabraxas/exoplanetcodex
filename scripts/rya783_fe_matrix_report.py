#!/usr/bin/env python3
"""RYA-783 — assemble and report the Fe product matrix.

    python3 scripts/rya783_fe_matrix_report.py

The deliverable of RYA-713's "Owed before Fe closes" item 4: 1D-LTE / Engine-A /
Engine-B-LTE / Engine-B-NLTE, each with its own value, sigma, line count and plot, per
(instrument x band x engine) — presented SIDE BY SIDE and NEVER combined (RYA-712). There
is no "primary": *"All Engines, LTE and NLTE, are products that get presented."*

WHAT THIS DOES NOT DO, deliberately
------------------------------------
* It does not average, weight, rank or arbitrate the engines. The per-engine spread is a
  RYA-525 cross-engine DIAGNOSTIC and never enters an error bar.
* It does not apply the optical FE_GATE [7.41, 7.51] to the near-UV or IR cells. Verified
  by inspection: no gate exists anywhere in the band-product path — the only gate there is
  fit quality (RYA-342). A frontier value of 7.55-7.64 is a RESULT, not a gate failure.
* It does not compare the IR to the optical anchor as a pass/fail. Per the RYA-777 clause
  ratified on this ticket, the IR validation reference is SAME-REGIME IR literature, and
  the IR-vs-optical difference is REPORTED, never flagged as a miss. The anchor appears as
  context on the plot only (Ryan: *"fine to show it on the matrix, but we should not hold
  back measured science because it doesn't fit a particular VIS value"*).
"""
from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_IN = ROOT / "data" / "results" / "band_products"
OUT_DIR = ROOT / "data" / "results" / "rya783"

# SINGLE-SOURCED from pipeline.band_products (RYA-836). This was a local hardcoded copy
# and it silently dropped every treatment added after it was written: the RYA-836 lab-gf
# cell reached the CSV but never the printed table or the plot, so the figure came out
# BYTE-IDENTICAL while the data underneath had gained a row. The cross-engine spread
# below does not use this list, so the report even contradicted itself — quoting a spread
# "over 2 engines" for a band whose table showed one.
#
# A products vocabulary with two definitions has the same shape as the gf column with two
# sources (RYA-353/825): the copy that is not the source drifts, and drifts quietly.
from pipeline.band_products import TREATMENTS as _CANON_TREATMENTS  # noqa: E402
TREATMENTS = list(_CANON_TREATMENTS)

# Context only — never a target, never a gate (RYA-161 / RYA-777).
OPTICAL_ANCHOR = 7.466          # RYA-553 banked, VIS — on the 3D-NLTE scale
FE_1D3D = 0.05                  # config.constants.FE_1D3D_SOLAR_OFFSET
IR_INREGIME = {
    "RYA-780 primary-only Fe IR": 7.5508,
}


def load(indir: Path) -> pd.DataFrame:
    rows = []
    for f in sorted(glob.glob(str(indir / "**" / "*_products.csv"), recursive=True)):
        d = pd.read_csv(f)
        d["_src"] = f
        rows.append(d)
    if not rows:
        raise SystemExit(f"no products under {indir} — run scripts/rya783_run_matrix.sh")
    return pd.concat(rows, ignore_index=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="indir", default=str(DEFAULT_IN))
    ap.add_argument("--no-plot", action="store_true")
    a = ap.parse_args()

    d = load(Path(a.indir))
    d = d[d.element == "Fe"].copy()
    # Each deck writes its own directory, so 1D-LTE and ENGINE-A (which do not depend on
    # the deck) appear once per deck with identical values. Dedupe on the product key.
    d = d.drop_duplicates(subset=["element", "ion", "band", "instrument", "treatment"],
                          keep="first").reset_index(drop=True)
    # Engine-B under two decks lands in two directories with the SAME treatment label
    # 'ENGINE-B'; the deck is what distinguishes them, so keep the last-writer-free split.
    d["treatment"] = d["treatment"].astype(str)

    print("=" * 100)
    print("RYA-783 — Fe PRODUCT MATRIX.  Separate products, never combined (RYA-712).")
    print("There is no 'primary': all engines are products that get presented.")
    print("=" * 100)

    for ion in ("I", "II"):
        sub = d[d.ion == ion]
        if not len(sub):
            continue
        print(f"\n### Fe {ion}")
        print(f"{'band':<14}{'treatment':<16}{'A(Fe)':>9}{'stat':>8}{'syst':>8}"
              f"{'n':>6}{'excl':>6}   dominant systematic")
        for band in sorted(sub.band.unique()):
            for t in TREATMENTS:
                r = sub[(sub.band == band) & (sub.treatment == t)]
                if not len(r):
                    continue
                for _, x in r.iterrows():
                    v = x.get("A")
                    v = float(v) if pd.notna(v) else float("nan")
                    print(f"{band:<14}{t:<16}{v:9.3f}"
                          f"{float(x.get('stat_dex', np.nan)):8.3f}"
                          f"{float(x.get('syst_dex', np.nan)):8.3f}"
                          f"{int(x.get('n_lines', 0)):6d}"
                          f"{int(x.get('n_excluded', 0)):6d}   "
                          f"{str(x.get('dominant', ''))[:34]}")

        # cross-engine spread — a DIAGNOSTIC (RYA-525), never an error bar
        print(f"\n  cross-engine spread (RYA-525 DIAGNOSTIC — never folded into a bar):")
        for band in sorted(sub.band.unique()):
            vals = sub[(sub.band == band)][["treatment", "A"]].dropna()
            if len(vals) < 2:
                continue
            v = vals.A.astype(float)
            print(f"    {band:<14} range {v.min():.3f}-{v.max():.3f}  "
                  f"spread {v.max() - v.min():.3f} dex over {len(v)} engines")

    # ── the comparison discipline this ticket binds ─────────────────────────
    print("\n" + "=" * 100)
    print("COMPARISON — IR/near-UV compare to SAME-REGIME references, NOT the optical anchor")
    print("=" * 100)
    print("  ⚠️ SCALE. Every cell above is a 1D product. Gold v3's Fe I 7.466 is on the")
    print("  3D-NLTE scale — the Magic-2013 1D->3D correction is APPLIED (RYA-553 ported")
    print(f"  7.516 -> 7.466, FE_1D3D_SOLAR_OFFSET = {FE_1D3D:.2f}). Comparing a 1D cell to")
    print("  7.466 mixes two scales, which is the RYA-669 ratchet: gold v3 carries the")
    print("  label '1D-NLTE (Fe I)' on a 3D value, the guard read the label, re-armed the")
    print("  correction, and A(Fe I) drifted to 7.416 with ALL NINE GATES GREEN. So the")
    print(f"  1D reference is {OPTICAL_ANCHOR + FE_1D3D:.3f}, not {OPTICAL_ANCHOR}.")
    print(f"    1D reference  A(Fe I) = {OPTICAL_ANCHOR + FE_1D3D:.3f}   (gold 3D {OPTICAL_ANCHOR} + {FE_1D3D:.2f})")
    for k, v in IR_INREGIME.items():
        print(f"  in-regime IR reference: {k} = {v}")
    ir = d[(d.ion == "I") & (d.band.astype(str).str.contains("red|IR", case=False))]
    for _, x in ir.iterrows():
        v = x.get("A")
        if pd.isna(v):
            continue
        v = float(v)
        line = f"  Fe I {x.band} {x.treatment:<14} {v:.3f}"
        for k, ref in IR_INREGIME.items():
            line += f"   vs {k}: {v - ref:+.3f}"
        line += f"   [vs 1D ref {v - (OPTICAL_ANCHOR + FE_1D3D):+.3f} — REPORTED, not a miss]"
        print(line)

    # RYA-832. This block used to read "NEAR-UV: NOT derivable by this route ... the
    # synthesis route is RYA-759, which is not merged. Reported as blocked-by-route."
    # Every clause of that is now false: 759 merged, the synthesis route is wired into
    # derive_band_products, and the near-UV is a first-class matrix cell above. Leaving
    # the old text would have been the "accounting doesn't describe reality" class this
    # project keeps paying for — a report describing a state the repo left behind.
    nuv = d[d.band.astype(str).str.lower().str.startswith("near")]
    if len(nuv):
        print("\n  NEAR-UV 3000-3780 A is DERIVED and first-class (RYA-832): synthesis-only")
        print("  by band policy — profile-fit and interval-integration are FORBIDDEN there")
        print("  (median line gap 0.146 A leaves no interval containing one profile and")
        print("  excluding its neighbours), and RYA-759 falsified profile fitting in band")
        print("  (901 candidates -> 0 measurable). The cell is 1D-LTE ONLY: UV Fe I is")
        print("  heavily over-ionised, so the missing NLTE correction is large and POSITIVE")
        print("  and a low value is expected physics, not a defect. Never coadded (RYA-712).")
    else:
        print("\n  NEAR-UV 3000-3780 A: no product present in this input directory.")
        print("  It is synthesis-only (band policy); derive it with")
        print("  `derive_band_products.py --lo 3000 --hi 3780`, which takes the RYA-832")
        print("  synthesis route automatically. Absent here means NOT RUN — not 'blocked'.")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / "fe_product_matrix.csv"
    d.to_csv(out, index=False)
    print(f"\nwrote {out.relative_to(ROOT)}  ({len(d)} product rows)")

    if not a.no_plot:
        _plot(d)


def _plot(d: pd.DataFrame) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"  (plot skipped: {e})")
        return
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)
    for ax, ion in zip(axes, ("I", "II")):
        sub = d[d.ion == ion].dropna(subset=["A"])
        if not len(sub):
            ax.set_visible(False)
            continue
        bands = sorted(sub.band.unique())
        for i, t in enumerate(TREATMENTS):
            r = sub[sub.treatment == t]
            if not len(r):
                continue
            x = [bands.index(b) + (i - 1.5) * 0.16 for b in r.band]
            y = r.A.astype(float)
            e = r.get("syst_dex", pd.Series([0] * len(r))).astype(float)
            ax.errorbar(x, y, yerr=e, fmt="o", capsize=3, label=t, ms=5)
            for xi, yi, ni in zip(x, y, r.n_lines.astype(int)):
                ax.annotate(f"n={ni}", (xi, yi), textcoords="offset points",
                            xytext=(0, 7), ha="center", fontsize=6)
        # anchor as CONTEXT ONLY — dashed, labelled, never a gate
        ax.axhline(OPTICAL_ANCHOR, ls="--", lw=0.9, color="0.45")
        ax.annotate("optical anchor 7.466 (VIS) — context only, not a target",
                    (0.02, OPTICAL_ANCHOR), xycoords=("axes fraction", "data"),
                    fontsize=7, color="0.35", va="bottom")
        ax.set_xticks(range(len(bands)))
        ax.set_xticklabels(bands, fontsize=8)
        ax.set_title(f"Fe {ion}", fontsize=11)
        ax.grid(alpha=0.25)
    axes[0].set_ylabel("A(Fe)")
    axes[0].legend(fontsize=8, loc="best")
    fig.suptitle("RYA-783 — Fe products per band x engine. Separate, never combined "
                 "(RYA-712). Bars are the STATED systematic budget.", fontsize=10)
    fig.tight_layout()
    p = OUT_DIR / "fe_product_matrix.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    print(f"wrote {p.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
