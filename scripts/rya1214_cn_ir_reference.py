#!/usr/bin/env python3
"""RYA-1214 — the 1D-to-1D comparison basis for the IR CN runs, and what the window cap drops.

Two things this ticket needs before an IR A(N) means anything.

🔴 1. THE COMPARISON MUST BE 1D-TO-1D, AND AMARSI PUBLISHES HIS OWN 1D ANSWER.
Our IR CN fit is 1D-LTE. AGSS21's adopted A(N) = 7.83 is 3D. Differencing them measures the
atmosphere and the abundance at once, which is the trap `rya1214_agss21_comparison.py`
documents for the optical indicators. But Amarsi et al. 2021 Table 2 carries FOUR abundance
columns per line, and `abundance_marcs` is the SAME lines reduced through a 1D MARCS model —
so a 1D-to-1D comparison exists and costs nothing but reading the right column.

On the lines our windows actually cover that column reads 7.940 (KP) / 7.941 (IAG), and
their own 1D(MARCS) -> 3D shift on those same lines is -0.081 dex. That -0.081 is the size
of the effect we are NOT applying, published by the author of the reference — which is why
comparing to 7.83 would have charged us for their 3D correction.

🔴 2. THE WINDOW CAP IS BINDING ON KITT PEAK AND NOT ON IAG, AND THE COUNTS DIVERGE.
`_CN_IR_MAX_WINDOWS = 12` exists for compute (every window is synthesised on every chi2
evaluation). I had described the two runs as "IAG 16 lines, Kitt Peak 54 lines, the whole
band". That is the count of telluric-clean lines in each holding's REACH, not the count
inside the windows the region fits:

  IAG  only 12 candidate windows exist below its 11083 A red edge, so the cap does not
       bind at all: 12 windows, all 16 reachable lines, a COMPLETE replication of the
       reachable part of AGSS21's line set.
  KP   48 candidate windows exist across the whole band; the cap keeps the 12 strongest by
       AGSS21's OWN published equivalent width and leaves 36 windows / 39 lines OUTSIDE.
       15 lines in. That is a SUBSET selected by published strength, not the band.

So IAG is the cleaner replication even though it reaches less of the band, and the KP number
must be read as a strength-ranked subset. Stating it the other way round would have implied
the wider-reaching holding was the better-sampled one.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "audit" / "rya1214_cno_products"
SRC = (ROOT / "data" / "reference" / "amarsi2021_cno" / "derived"
       / "amarsi2021_cno_molecular_lines.csv")

#: Amarsi's four abundance columns, and what each one IS. The comparison column is named
#: here once so the choice is a declaration and not whichever happened to be closer
#: (RYA-161).
COLUMNS = {
    "abundance_marcs":   ("1D MARCS", True,
                          "the same lines through a 1D model — OUR comparison basis"),
    "abundance_hm74":    ("1D HM74", False, "a different 1D atmosphere, carried for scale"),
    "abundance_3d":      ("3D", False, "their adopted value's frame; not ours"),
    "abundance_mean_3d": ("<3D>", False, "mean-stratified 3D; not ours"),
}


def _constants() -> tuple[float, int, float, float]:
    """Read the region's own constants from its source — a copy here could drift (RYA-845)."""
    src = (ROOT / "pipeline" / "cno_synthesis.py").read_text()
    pad = float(re.search(r"^_CN_IR_PAD_A = ([\d.]+)", src, re.M).group(1))
    cap = int(re.search(r"^_CN_IR_MAX_WINDOWS = (\d+)", src, re.M).group(1))
    lo, hi = (float(x) for x in
              re.search(r"H2O_LO, H2O_HI = ([\d.]+), ([\d.]+)", src).groups())
    return pad, cap, lo, hi


def _lines() -> pd.DataFrame:
    d = pd.read_csv(SRC)
    c = d[(d.element_parameter == "logepsN") & (d.species == "CN")].copy()
    lv = c.wavelength_vac_nm * 10.0
    s2 = (1e4 / lv) ** 2
    c["air_A"] = lv / (1 + 0.0000834254 + 0.02406147 / (130 - s2)
                       + 0.00015998 / (38.9 - s2))
    c["ew_A"] = c.equivalent_width_pm * 10.0
    return c.sort_values("air_A").reset_index(drop=True)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    pad, cap, h2o_lo, h2o_hi = _constants()
    c = _lines()

    rows = list(zip(c.air_A, c.ew_A))
    clusters, cur = [], [rows[0]]
    for w, ew in rows[1:]:
        if w - cur[-1][0] <= 2 * pad:
            cur.append((w, ew))
        else:
            clusters.append(cur)
            cur = [(w, ew)]
    clusters.append(cur)

    cand, h2o_win, h2o_lines = [], 0, 0
    for cl in clusters:
        lo = min(x[0] for x in cl) - pad
        hi = max(x[0] for x in cl) + pad
        if hi > h2o_lo and lo < h2o_hi:
            h2o_win += 1
            h2o_lines += len(cl)
            continue
        cand.append((round(lo, 2), round(hi, 2), max(x[1] for x in cl), len(cl)))

    kept = {
        "nir_cn_kp": (sorted(sorted(cand, key=lambda t: -t[2])[:cap]), cand),
        "nir_cn_iag": (sorted(sorted([x for x in cand if x[1] <= 11083.0],
                                     key=lambda t: -t[2])[:cap]),
                       [x for x in cand if x[1] <= 11083.0]),
    }

    out, sizing = [], {}
    print("=== RYA-1214 — the IR CN comparison basis, and what the window cap drops ===")
    print(f"{len(c)} published CN lines -> {len(clusters)} clusters at +/-{pad} A; "
          f"{h2o_win} window(s) / {h2o_lines} line(s) excluded as H2O {h2o_lo:.0f}-{h2o_hi:.0f}; "
          f"{len(cand)} candidate windows; cap = {cap}")
    for region, (keep, pool) in kept.items():
        w2 = [(a, b) for a, b, _, _ in keep]
        ins = c[c.air_A.apply(lambda a: any(lo <= a <= hi for lo, hi in w2))]
        binding = len(pool) > len(keep)
        sizing[region] = {
            "candidate_windows_reachable": len(pool),
            "windows_kept": len(keep),
            "windows_dropped_by_cap": len(pool) - len(keep),
            "lines_in_windows": int(len(ins)),
            "lines_left_outside_the_cap": int(sum(x[3] for x in pool) - len(ins)),
            "span_A": round(sum(b - a for a, b in w2), 1),
            "cap_is_binding": bool(binding),
            "replication_status": ("SUBSET selected by AGSS21's own published equivalent "
                                  "width — not the band" if binding else
                                  "COMPLETE for the part of the band this holding reaches"),
        }
        print(f"\n  {region}: {len(keep)} windows, {sizing[region]['span_A']} A, "
              f"{len(ins)} line(s) in "
              f"({'CAP BINDS — ' if binding else 'cap does not bind — '}"
              f"{sizing[region]['windows_dropped_by_cap']} window(s) / "
              f"{sizing[region]['lines_left_outside_the_cap']} line(s) outside)")
        print(f"    {sizing[region]['replication_status']}")
        for col, (label, is_basis, why) in COLUMNS.items():
            v = ins[col]
            mark = "<--" if is_basis else "   "
            print(f"    {mark} {label:10s} median {v.median():.3f}  sd {v.std(ddof=1):.3f}  "
                  f"[{v.min():.3f}, {v.max():.3f}]   {why}")
            out.append({"region": region, "n_windows": len(keep),
                        "n_lines": int(len(ins)), "amarsi_column": col,
                        "column_label": label, "is_comparison_basis": is_basis,
                        "median": round(float(v.median()), 3),
                        "sd": round(float(v.std(ddof=1)), 3),
                        "min": round(float(v.min()), 3),
                        "max": round(float(v.max()), 3), "what_it_is": why})
        shift = float(ins.abundance_3d.median() - ins.abundance_marcs.median())
        print(f"    their OWN 1D(MARCS)->3D shift on these exact lines: {shift:+.3f} dex "
              f"— the effect we are not applying, sized by the reference's own author")
        sizing[region]["their_1D_to_3D_shift_dex"] = round(shift, 3)
        sizing[region]["comparison_basis_A_N"] = round(float(ins.abundance_marcs.median()), 3)

    pd.DataFrame(out).to_csv(OUT / "amarsi2021_cn_ir_reference.csv", index=False)
    (OUT / "amarsi2021_cn_ir_reference.prov.json").write_text(json.dumps({
        "ticket": "RYA-1214",
        "source": "Amarsi et al. 2021, A&A 656, A113, Table 2 (the CN A-X (0-0) rows), via "
                  "data/reference/amarsi2021_cno/derived/",
        "why_1D_to_1D": "our fit is 1D-LTE and AGSS21's adopted 7.83 is 3D. Amarsi publishes "
                        "`abundance_marcs` — the same lines through a 1D model — so the "
                        "comparison can be like-for-like. Differencing against 7.83 would "
                        "charge us for a 3D correction we did not apply.",
        "cap_finding": "`_CN_IR_MAX_WINDOWS = 12` BINDS on Kitt Peak (48 candidates, 36 "
                       "dropped) and does NOT bind on IAG (12 candidates below its 11083 A "
                       "edge). So IAG is the COMPLETE replication of what it reaches and KP "
                       "is a strength-ranked SUBSET — the opposite of what 'KP reaches the "
                       "whole band' implies about sampling.",
        "n_published_lines": int(len(c)), "n_clusters": len(clusters),
        "pad_A": pad, "window_cap": cap,
        "h2o_excluded": {"band_A": [h2o_lo, h2o_hi], "windows": h2o_win,
                         "lines": h2o_lines,
                         "note": "excluded, not corrected — both holdings are "
                                 "telluric-corrected and could be argued into fitting "
                                 "there; 440 A of 2333 is not worth leaning on it"},
        "regions": sizing,
    }, indent=2) + "\n")
    print(f"\n  wrote amarsi2021_cn_ir_reference.csv / .prov.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
