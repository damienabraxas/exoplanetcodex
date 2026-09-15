#!/usr/bin/env python3
"""RYA-1214 — WHY optical N I cannot be measured: N I is a minority absorber in its own lines.

Ryan asked what the story with nitrogen is. This is the first of the three answers, and it is
the one that cannot be engineered around: in every one of AGSS21's optical/red-optical N I
fit windows, the N I line itself carries a few per cent of the absorption and iron, silicon,
titanium or CN carries most of it. A fit that minimises chi2 over such a window is fitting
the blend, and the abundance it returns for nitrogen is whatever makes the blend fit.

HOW DOMINANCE IS MEASURED, AND THE TWO TRAPS
--------------------------------------------
🔴 PER LINE, NEVER PER WINDOW OR PER POOL. RYA-1195's lesson is that a window-level or
pool-level statistic credits a neighbour's property to this line: `pool_label_coverage`
credited any labelled row within 0.05 A and claimed 12 where 11 was right, the label
belonging to a line 2.5 dex weaker. So the unit here is THE LINE, and the window is only the
integration domain.

🔴 SUMMED `central_depth`, NOT A COUNT. Counting neighbours makes 200 shallow lines look
worse than one saturated one. `linelist_solar.central_depth` is the catalogue's own predicted
depth, so summing it over the window and taking the target line's share is a share of
ABSORPTION, which is what a chi2 actually weights. This is the same discriminant RYA-1189
used to show the near-UV continuum shift was blend-driven.

⚠️ THIS IS A CATALOGUE MEASUREMENT, NOT A SYNTHESIS. It reads predicted depths, so it says
what fraction of the CATALOGUED absorption belongs to N I. Uncatalogued absorption would
make the fraction smaller, never larger, so the number is an UPPER BOUND on N I's share —
which is the conservative direction for the claim being made.

The window half-width is the fit's own, read from `cno_synthesis` rather than chosen here, so
the number describes the windows the products were actually measured in.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "audit" / "rya1214_cno_products"
LINELIST = ROOT / "data" / "linelists" / "linelist_solar.csv"
SETS = ROOT / "data" / "linelists" / "reference_sets"

#: AGSS21's `source_band` -> our band name, so each line is measured over the half-width its
#: OWN product was fitted with.
BAND_OF = {"RED_OPTICAL": "red-optical", "NIR": "NIR", "VIS": "VIS",
           "NEAR_UV": "near-UV", "H": "H"}
#: Swept so the finding is not an artifact of one width. The point of the sweep is that the
#: share is a strong function of it: at +/-0.25 A a red-optical N I line carries the MAJORITY
#: of its window. The fit does not use 0.25 A - `config.synth_bands` sets red-optical to
#: 1.10 A and NIR to 1.40 A - so the number that describes the PRODUCT is the one at the
#: band's own width, and quoting the narrow one would describe a fit nobody ran.
SWEEP_A = (0.25, 0.5, 0.62, 1.1, 1.4, 2.0)


def _half_widths() -> dict:
    """Per-band fit half-width, from `config.synth_bands` - the same SSOT the products used."""
    import sys
    sys.path.insert(0, str(ROOT))
    from config.synth_bands import SYNTH_BANDS
    return {b: float(c.half_width_A) for b, c in SYNTH_BANDS.items()}


def _species(df: pd.DataFrame) -> pd.Series:
    ion = df.ion.astype(str).str.strip()
    return df.element.astype(str).str.strip() + " " + ion


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    ll = pd.read_csv(LINELIST)
    ll = ll[ll.central_depth.notna() & (ll.central_depth > 0)].copy()
    ll["sp"] = _species(ll)

    hws = _half_widths()
    rows, sweep_rows = [], []
    # The set files are read from their OWN manifest, never from a guessed name: C and O are
    # stored as `_grid_` (Amarsi 2019 analysis grid) and only N as `_adopted_`, and a hardcoded
    # `_adopted_` name made C and O vanish from this census without a word.
    manifest = json.loads((SETS / "agss21_cno_lineset_rya1214.prov.json").read_text())
    by_el = {f["element"]: f for f in manifest["files"]}
    for el, target_sp in (("C", "C I"), ("N", "N I"), ("O", "O I")):
        if el not in by_el:
            raise KeyError(f"the AGSS21 set manifest has no {el} entry; refusing to report a "
                           f"census with an element silently absent")
        use_status = by_el[el]["use_status"]
        want = pd.read_csv(ROOT / by_el[el]["file"])
        # THE UNIT IS THE PUBLISHED FEATURE (line_label), not the component row. The O grid
        # lists 615.6nm as three components and 926.1nm as three with different gf; taken
        # row by row, 615.6nm was reported three times and 844.6nm read 7.9 / 21.0 / 15.5 %
        # depending only on which component the window happened to be centred on.
        for label, comps in want.groupby("line_label", sort=False):
            w = comps.iloc[0]
            comp_A = comps.wavelength_air_A.astype(float).to_numpy()
            lam = float(comp_A.mean())
            band = BAND_OF.get(str(w.source_band).strip().upper())
            own_hw = hws.get(band)
            if own_hw is None:
                print(f"  {el} {label}: source_band {w.source_band!r} maps to no "
                      f"band in config.synth_bands - SKIPPED, not defaulted")
                continue
            # A SET, not a concatenation: the band's own width (1.1 / 1.4 A) is already IN
            # the sweep, and appending it again emitted every line twice.
            for hw in sorted(set(SWEEP_A) | {own_hw}):
                near = ll[(ll.wavelength_air_A >= comp_A.min() - hw)
                          & (ll.wavelength_air_A <= comp_A.max() + hw)]
                if not len(near):
                    continue
                tot = float(near.central_depth.sum())
                # THE LINE, not the species in the window: the target row is the one within
                # the set's own match tolerance, so a second N I line 0.4 A away is a BLEND
                # partner here and is counted as one (RYA-1195).
                tol = float(w.match_tol_A)
                near_A = near.wavelength_air_A.to_numpy(float)
                own = (abs(near_A[:, None] - comp_A[None, :]) <= tol).any(axis=1)
                mine = near[(near.sp == target_sp).to_numpy() & own]
                share = float(mine.central_depth.sum()) / tot if tot else float("nan")
                by_sp = (near.groupby("sp").central_depth.sum() / tot).sort_values(
                    ascending=False)
                sweep_rows.append({"element": el, "line_label": label, "line_A": round(lam, 3),
                                   "half_width_A": hw,
                                   "target_share": round(share, 4),
                                   "n_blend_rows": int(len(near) - len(mine))})
                if hw != own_hw:
                    continue
                top = [(s, round(float(v), 3)) for s, v in by_sp.head(4).items()]
                rows.append({
                    "element": el, "species": target_sp, "set_use_status": use_status,
                    "line_A": round(lam, 3),
                    "line_label": label, "n_components": int(len(comps)), "band": band,
                    "half_width_A": hw,
                    "half_width_basis": "config.synth_bands[band].half_width_A - the "
                                        "product's own fit width, not a width chosen here",
                    "target_line_found": bool(len(mine)),
                    "target_share_of_catalogued_depth": round(share, 4),
                    "target_share_pct": round(100 * share, 1),
                    "n_catalogued_rows_in_window": int(len(near)),
                    "n_blend_rows": int(len(near) - len(mine)),
                    "dominant_species": (top[0][0] if top else ""),
                    "dominant_share_pct": (round(100 * top[0][1], 1) if top else None),
                    "top4": "; ".join(f"{s} {100*v:.0f}%" for s, v in top),
                })

    d = pd.DataFrame(rows)
    s = pd.DataFrame(sweep_rows)
    d.to_csv(OUT / "n_blend_dominance.csv", index=False)
    s.to_csv(OUT / "n_blend_dominance_sweep.csv", index=False)

    pd.set_option("display.width", 250)
    print("=== RYA-1214 - whose absorption is in the CNO indicator windows? ===")
    print("    each line over ITS OWN band's fit half-width from config.synth_bands: "
          + ", ".join(f"{b} {v}A" for b, v in sorted(hws.items())))
    print(d[["element", "line_label", "band", "target_share_pct",
             "n_blend_rows", "dominant_species", "dominant_share_pct",
             "top4"]].to_string(index=False))

    print("\n--- the claim, per element: the target line's share of its own window ---")
    summary = {}
    for el, g in d.groupby("element"):
        v = g.target_share_pct
        summary[el] = {"n_lines": int(len(g)), "min_pct": round(float(v.min()), 1),
                       "max_pct": round(float(v.max()), 1),
                       "median_pct": round(float(v.median()), 1)}
        print(f"  {el}: {len(g)} indicator line(s), share {v.min():.1f}-{v.max():.1f}% "
              f"(median {v.median():.1f}%)")

    print("\n--- swept over half-width. 🔴 THE SHARE IS A STRONG FUNCTION OF IT: at "
          "+/-0.25 A a red-optical N I line carries the MAJORITY of its window. The fit "
          "does not use 0.25 A. ---")
    piv = s.pivot_table(index="element", columns="half_width_A",
                        values="target_share", aggfunc="median")
    print((100 * piv).round(1).to_string())

    (OUT / "n_blend_dominance.prov.json").write_text(json.dumps({
        "ticket": "RYA-1214",
        "question": "Ryan: 'So what is the story with N?' — answer 1 of 3",
        "what": "the target indicator line's share of the CATALOGUED absorption in its own "
                "fit window, per line, summed central_depth",
        "unit_is_the_published_feature": "one row per AGSS21 line_label; a multi-component "
                            "feature's target is every target-species row within the set's "
                            "match tolerance of ANY of its components, and its window spans "
                            "all components. Counted per component row, O read median 6.8% "
                            "over 26 'lines' (615.6nm three times; 844.6nm 7.9/21.0/15.5% by "
                            "centring alone); per feature it is 28.9% over 13.",
        "unit_is_the_line": "the target row must lie within the reference set's OWN match "
                            "tolerance of the published wavelength. A second line of the "
                            "same species further out in the window is a BLEND PARTNER and "
                            "is counted as one — crediting it to the target is RYA-1195's "
                            "defect (12 claimed where 11 was right).",
        "why_depth_not_count": "counting neighbours makes 200 shallow lines look worse than "
                               "one saturated one. Summed central_depth is a share of "
                               "ABSORPTION, which is what a chi2 weights (RYA-1189).",
        "direction_of_the_bound": "catalogue-predicted depths only, so uncatalogued "
                                  "absorption would make the target's share SMALLER, never "
                                  "larger. The number is an UPPER BOUND on N I's share — the "
                                  "conservative direction for the claim.",
        "half_width_per_band": hws,
        "half_width_basis": "config.synth_bands - the SSOT the products themselves used. An "
                            "earlier pass of this script used a flat 0.5 A, which is not any "
                            "band's width: red-optical is 1.10 A and NIR 1.40 A.",
        "swept_A": list(SWEEP_A),
        "the_share_depends_on_the_width": "at +/-0.25 A a red-optical N I line carries the "
                                          "MAJORITY of its window, so this finding is not "
                                          "'N I is intrinsically weak' - it is 'N I is a "
                                          "minority absorber over the width the fit "
                                          "actually integrates'. Quoting the narrow number "
                                          "would describe a fit nobody ran; quoting only "
                                          "the wide one without this sentence would hide "
                                          "that a narrower window is a real lever.",
        "per_element": summary,
        "consequence": "where the target carries a few per cent, the fit is fitting the "
                       "blend and the abundance it returns for the target is whatever makes "
                       "the blend fit. This is why optical/red-optical N I reads +0.46 to "
                       "+0.52 dex above AGSS21 and why AGSS21's own nitrogen is an IR CN "
                       "measurement.",
        "rows": rows,
    }, indent=2) + "\n")
    print(f"\n  wrote n_blend_dominance.csv / _sweep.csv / .prov.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
