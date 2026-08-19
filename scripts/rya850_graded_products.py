#!/usr/bin/env python3
"""
RYA-850 — promote the graded (lab-gf) pool to the reported value
================================================================
Two products per cell, GRADED primary, UNGRADED kept and labelled.

WHAT WAS ACTUALLY MISSING. Nothing here is a new measurement. The lab-gf pools were built
by RYA-824 (VIS/IR) and RYA-836 (near-UV), `error_budget.gf_term(graded=...)` already
carries both terms (0.041 vs 0.170), and the near-UV lab-gf cell already reaches the
matrix. What never happened is the wiring: `derive_band_products.py:770` passes
`gf_graded=False` UNCONDITIONALLY, so every EW-route cell is charged the Kurucz 0.17
regardless of what its lines are, and RYA-824's two pools were never emitted as cells at
all. This script closes that, and reports the pair.

🔴 ONLY THREE OF EIGHTEEN CELLS CAN HAVE A GRADED TWIN TODAY, and that is a finding, not a
shortfall of this script. A lab-gf pool exists for Fe I in the near-UV (RYA-836, n=59),
Fe I VIS (RYA-824, n=9) and Fe I red-optical (RYA-824, n=20). There is NO primary-lab-gf
pool for Fe II in any band, nor for any ENGINE-A/ENGINE-B product — nobody has measured
one. Those cells fall back to UNGRADED and say so, rather than being relabelled graded on
the strength of a term they are not entitled to.

THE PUBLISHED BAR IS THE POOL'S OWN CITED LAB SIGMA (Ryan, 2026-08-17; RYA-850 amended).
`GRADED_GF_SYSTEMATIC_DEX` is a BOUND — the worst grade we accept — and a bound is the
right answer only while the actual sigmas are unknown. These pools publish per-line
uncertainties, which MEASURE the same quantity, and RYA-853 refereed them line-by-line
against the source papers (Ruffoni 142/142, Den Hartog 203/203 perfect on value AND cited
sigma) before they were allowed near a published bar. They are larger, not smaller:
0.0600 (VIS), 0.0524 (red-optical), 0.0522 (near-UV) against the generic 0.041, so the
bound was optimistic by 15-25% on the two EW-route bands. Both are printed side by side so
the widening is visible rather than quietly absorbed.

⚠️ NO COMBINED HEADLINE IS COMPUTED. The ticket asks for a "combined/headline Fe value".
RYA-712 forbids a cross-band combined product and `pipeline.band_products` deliberately has
no `combine()` — the per-band values are separate measurements whose spread IS the result.
So this emits the per-band graded table and names the candidate headline cell; collapsing
them into one number is a reporting decision for a human, not arithmetic to do quietly here.

⚠️ THE LINE SETS ARE PRE-RYA-847, and that is survivable rather than blocking. 847 is 6/7
done but a deliberate NUMERICAL NO-OP today — `SYNTH_CONSTRAINT = None`, its near-UV
control excludes 0 lines — because the threshold is set by a sweep that is still running.
So nothing here contradicts it. When the sweep ratifies a cut, the near-UV cell is the
exposed one (synthesis route) and must be regenerated; 847's own scope note puts the
EW-route products, `1D-LTE-LABGF in VIS/IR` among them BY NAME, outside its remit, so the
VIS and red-optical cells stand.

Usage:
    python3 scripts/rya850_graded_products.py
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.band_policy import resolve as resolve_band          # noqa: E402
from pipeline.error_budget import build as build_budget           # noqa: E402
from pipeline import harness_residual                            # noqa: E402  RYA-869
from pipeline.graded_reporting import (                           # noqa: E402
    GRADED_SUFFIX, annotate, element_table, format_value, is_graded,
    pair_products,
    stat_from_scatter, total_sigma)

MATRIX = ROOT / "data" / "results" / "rya783" / "fe_product_matrix.csv"
RYA824 = ROOT / "data" / "results" / "rya824" / "rya824_lab_gf_per_line.csv"
OUT = ROOT / "data" / "results" / "rya850"

#: The EW/profile-fit handler's MEASURED residual, as `derive_band_products` charges it.
#: Kept identical so a graded cell differs from its ungraded twin in the gf term ALONE.
#: RYA-869 — IMPORTED, not restated. This was the second of two hand-written copies of
#: one number (the other was `derive_band_products.py:77`); they agreed, which is what a
#: duplicated constant does right up until it does not (RYA-845).
PROFILE_FIT_RESIDUAL_DEX = harness_residual.HANDLER_RESIDUAL_DEX["ProfileFitHandler"]

#: The pool's OWN cited laboratory sigma, derived from the lines it actually contains
#: rather than hardcoded. RYA-853 refereed this table against the source papers and found
#: Ruffoni 142/142 and Den Hartog 203/203 PERFECT on both value and cited sigma, so these
#: numbers are audited data now, not a plausible-looking alternative.
LAB_CSV = ROOT / "data" / "reference" / "fe_gf_lab" / "fe1_lab_loggf.csv"
RYA836_PER_LINE = (ROOT / "data" / "results" / "rya836"
                   / "rya836_nearuv_lab_gf_per_line.csv")
CITED_MATCH_TOL_A = 0.05


#: A pool whose cited sigmas cover only part of it is not described by their RMS -- the
#: unmatched lines would silently inherit the matched ones' uncertainty. Below this the
#: cited term is REFUSED and the generic bound stands, which is the honest fallback.
#: IMPORTED, not restated (RYA-855): `pipeline.gf_rung` applies the same threshold inside
#: `derive_band_products`, and a number declared in two places is the RYA-845 defect.
from pipeline.gf_rung import CITED_COVERAGE_MIN                    # noqa: E402


def cited_sigma(waves) -> dict:
    """RMS of the per-line cited laboratory sigma over a pool.

    RMS rather than the median because these combine in quadrature into the budget, and
    the median throws away exactly the tail that a quadrature sum is sensitive to: the
    IR pool's median is 0.030 while its RMS is 0.0524. Reproduces RYA-824's published
    0.0524 / 0.0600 exactly, which is the check that this is the same quantity 824 meant.

    A wavelength window is not a unique line ID (RYA-853): where it catches more than one
    lab row the line is counted as UNMATCHED rather than resolved by `iloc[0]`, because
    picking the first match is how RYA-853 manufactured 12-dex "defects". Ambiguity is
    reported so it can never be mistaken for coverage.
    """
    lab = pd.read_csv(LAB_CSV)
    got, ambiguous = [], 0
    for w in waves:
        m = lab[np.abs(lab.wavelength_air_A - float(w)) <= CITED_MATCH_TOL_A]
        if len(m) == 1:
            got.append(float(m.iloc[0].e_loggf_dex))
        elif len(m) > 1:
            ambiguous += 1
    n_pool = len(list(waves))
    if not got:
        return {"n": 0, "n_pool": n_pool, "coverage": 0.0, "ambiguous": ambiguous,
                "cited_sigma": float("nan"), "usable": False}
    a = np.asarray(got, dtype=float)
    cov = len(a) / n_pool if n_pool else 0.0
    return {"n": len(a), "n_pool": n_pool, "coverage": cov, "ambiguous": ambiguous,
            "cited_sigma": float(np.sqrt((a ** 2).mean())),
            "usable": cov >= CITED_COVERAGE_MIN}


#: Named, because a budget term is never unsourced -- and RYA-853 refereed every one of
#: these against the source papers before this term was allowed to set a published bar.
CITED_SOURCE = ("Ruffoni+2014 MNRAS 441 3127; Den Hartog+2014 ApJS 215 23; "
                "Belmonte+2017 ApJ 848 125; Den Hartog+2019 ApJS 243 33 "
                "-- refereed line-by-line by RYA-853")


def cited_by_band_map() -> dict:
    """The cited sigma of each graded pool, keyed by the matrix's band names."""
    per_line_waves = {}
    d824 = pd.read_csv(RYA824)
    per_line_waves["red-optical"] = d824[d824.band == "IR"].wavelength_air_A.tolist()
    per_line_waves["VIS"] = d824[d824.band == "VIS"].wavelength_air_A.tolist()
    if RYA836_PER_LINE.exists():
        d836 = pd.read_csv(RYA836_PER_LINE)
        per_line_waves["near-UV"] = d836[
            np.isfinite(d836.a_labgf)].wavelength_air_A.tolist()
    return {band: cited_sigma(waves) for band, waves in per_line_waves.items()}


def apply_cited_gf(combined: pd.DataFrame, cited: dict) -> tuple[pd.DataFrame, list]:
    """Recharge every GRADED cell already in the matrix onto its pool's cited gf sigma.

    The cells built here take the cited term at construction; the near-UV graded cell
    arrives from the matrix carrying the generic bound, and is entitled to the same
    treatment. Doing it in one pass over the combined frame is what keeps the three
    graded cells on ONE gf convention -- if the near-UV kept the bound while the other
    two moved, the band-to-band spread would mix two definitions of the bar.

    The scatter is recovered EXACTLY as `stat * sqrt(N)` (the definition of the stat term)
    rather than re-read, so this cannot drift from the value the cell was built with.
    """
    out = combined.copy()
    recharged = []
    for i, r in out.iterrows():
        # RYA-906 — ask the `gf` AXIS, falling back to the suffix for rows that predate
        # the axis columns. A raw `.endswith` on a label is the RYA-869 shape: it is a
        # string test standing in for a property of the product, and it stops being true
        # the moment the naming grows a variant. `is_graded` owns that decision now, so
        # there is one place to change when the naming moves again.
        if not is_graded(r.get("treatment", ""), gf=r.get("gf")):
            continue
        c = cited.get(str(r["band"]))
        if not c or not c["usable"]:
            continue
        n = int(r["n_lines"])
        scatter = float(r["stat_dex"]) * np.sqrt(n)
        # 🔴 RYA-869 — the harness term follows the HANDLER THE CELL DECLARES. What stood
        # here was `near_uv = str(r["band"]) == "near-UV"`, an equality against a BAND
        # name used to pick a handler: right today only because the one near-UV graded
        # cell happens to be the synthesis route and the two RYA-824 cells happen to be
        # the profile fitter. It is the same shape as the deriver's
        # `treatment == "ENGINE-B"` — a handler property derived from a label — and it
        # would have been wrong the first time a VIS cell arrived from a flux fit, which
        # `ENGINE-B` already is.
        b = build_budget("Fe", BAND_PIVOT_A[str(r["band"])], n,
                         scatter_dex=scatter, gf_graded=True,
                         **harness_residual.for_handler(
                             str(r["handler"])).budget_kwargs(),
                         cited_gf_sigma_dex=c["cited_sigma"],
                         cited_gf_source=CITED_SOURCE)
        _, syst = b.total()
        before = float(r["syst_dex"])
        after = round(float(syst), 4)
        if abs(after - before) > 5e-5:
            recharged.append({"band": str(r["band"]), "engine": str(r["treatment"]),
                              "syst_generic": before, "syst_cited": after,
                              "total_generic": total_sigma(float(r["stat_dex"]), before),
                              "total_cited": total_sigma(float(r["stat_dex"]), after)})
        out.at[i, "syst_dex"] = after
        out.at[i, "dominant"] = b.dominant().name if b.dominant() else ""
    return out, recharged


#: 824 named its bands by instrument regime; the matrix names them by band policy.
BAND_FROM_824 = {"IR": "red-optical", "VIS": "VIS"}

#: One representative wavelength per band, only ever used to resolve the BAND POLICY
#: (which continuum/telluric terms apply) -- never as a measurement.
BAND_PIVOT_A = {"VIS": 5000.0, "red-optical": 8000.0, "near-UV": 3390.0}


def graded_cells_from_824() -> pd.DataFrame:
    """Emit RYA-824's two lab-gf pools as matrix cells, on the graded gf term.

    The pools were measured and published in RYA-824 and then never became cells. This
    re-aggregates the SAME per-line abundances — no re-inversion — and charges the budget
    the graded term, which is the entire difference from the ungraded twin.

    Built on the GENERIC bound deliberately: `apply_cited_gf` is the single place the
    cited sigma is applied, so all three graded cells travel the same code path and the
    published-vs-bound comparison covers every one of them rather than only the cell that
    happened to arrive from the matrix.
    """
    d = pd.read_csv(RYA824)
    rows = []
    for band824, grp in d.groupby("band"):
        band = BAND_FROM_824.get(str(band824))
        if band is None:
            continue
        ok = grp[np.isfinite(grp.a_labgf)]
        if len(ok) < 3:
            continue
        value = float(ok.a_labgf.median())
        scatter = float(ok.a_labgf.std(ddof=1))
        n = int(len(ok))

        pol = resolve_band(float(ok.wavelength_air_A.median()))
        # RYA-824's pools are profile-fit EW inversions re-aggregated on the laboratory
        # gf — no re-measurement, so the same handler and the same residual (RYA-869).
        hr = harness_residual.for_handler("ProfileFitHandler")
        b = build_budget("Fe", float(ok.wavelength_air_A.median()), n,
                         scatter_dex=scatter, gf_graded=True, **hr.budget_kwargs())
        stat, syst = b.total()
        rows.append({
            "element": "Fe", "ion": str(ok.ion.iloc[0]), "band": pol.name,
            "instrument": "kpno_solar_atlas", "treatment": "1D-LTE-LABGF",
            "handler": hr.handler,
            "A": round(value, 3), "n_lines": n, "n_excluded": int(len(grp) - len(ok)),
            "stat_dex": round(float(stat), 4), "syst_dex": round(float(syst), 4),
            "dominant": (b.dominant().name if b.dominant() else ""),
            "_src": "rya824 lab-gf pool (RYA-850 wiring)",
            "scatter_dex": round(scatter, 4),
        })
    return pd.DataFrame(rows)


#: The route token in a band-product file stem, e.g.
#: `.../FeI_3000_3780_kpno_solar_atlas_SYNTH_products.csv` -> `SYNTH`. Same stem grammar
#: as `scripts/rya855_rung_audit.py::_STEM`.
_SRC_ROUTE = re.compile(r"_(PROFILEFIT|SYNTH|LABGF)_products\.csv$")


def _route_from_src(src: str) -> str:
    """Which deriver route wrote this matrix row, read off its `_src` file name.

    🔴 THE ROUTE, NOT THE ENGINE. An `ENGINE-B` cell lives in a `*_PROFILEFIT_*` stem
    because that stem names where the LINE SET came from; the treatment is what tells
    `handler_of_banked_cell` that this particular cell re-fitted the flux.
    """
    m = _SRC_ROUTE.search(str(src))
    if m:
        return m.group(1)
    if "rya824 lab-gf pool" in str(src):
        return "PROFILEFIT"      # emitted by graded_cells_from_824, an EW re-aggregation
    raise ValueError(
        f"cannot tell which route produced matrix row {src!r}, so the handler that "
        f"earned its harness residual is unknown. Re-derive the matrix (the products "
        f"table carries a `handler` column since RYA-869) rather than guessing.")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    if not MATRIX.exists():
        raise SystemExit(f"Fe matrix absent: {MATRIX}")

    matrix = pd.read_csv(MATRIX)
    # RYA-869 — the matrix must carry the handler that produced each cell, because the
    # recharge below charges that handler's measured residual. `derive_band_products`
    # writes a `handler` column now and `rya783_fe_matrix_report` concatenates it
    # straight through, so a freshly derived matrix arrives with it. A matrix banked
    # BEFORE that column existed is recovered once, here, from the route recorded in
    # `_src` — never per call site, and never from the treatment name.
    if "handler" not in matrix.columns or matrix["handler"].isna().any():
        matrix["handler"] = [
            harness_residual.handler_of_banked_cell(
                route=_route_from_src(str(src)), treatment=str(t))
            for src, t in zip(matrix["_src"], matrix["treatment"])]
        print(f"[RYA-869] matrix predates the handler column — recovered from the route "
              f"in _src for {len(matrix)} rows: "
              f"{dict(matrix.handler.value_counts())}")

    # The pools' own cited laboratory sigma SETS the published bar (Ryan, 2026-08-17);
    # the generic 0.041 is a bound and is kept only as the documented comparison.
    cited = cited_by_band_map()
    print("=== the gf term each graded pool is entitled to ===")
    for band, c in sorted(cited.items()):
        mark = "USED" if c["usable"] else "REFUSED (coverage)"
        print(f"  {band:<12} n={c['n']:3d}/{c['n_pool']:<3d} "
              f"coverage {c['coverage']:.0%}  ambiguous {c['ambiguous']}  "
              f"cited {c['cited_sigma']:.4f} ({c['cited_sigma']/0.041:.1f}x generic)  {mark}")

    new = graded_cells_from_824()
    print(f"[wire] RYA-824 lab-gf pools -> {len(new)} new graded cells")
    for _, r in new.iterrows():
        print(f"       {r.band:<12} {r.ion:<3} n={r.n_lines:<3} "
              f"A={r.A:.3f}  stat={r.stat_dex:.4f} syst={r.syst_dex:.4f} "
              f"(scatter {r.scatter_dex:.4f})")

    combined = pd.concat([matrix, new.drop(columns=["scatter_dex"])], ignore_index=True)
    combined, recharged = apply_cited_gf(combined, cited)
    d = annotate(combined)

    # ── the pairing ───────────────────────────────────────────────────────────────
    pairs = pair_products(combined)
    n_graded = sum(1 for p in pairs if p.primary_is_graded)
    print(f"\n=== {len(pairs)} cells: {n_graded} GRADED, {len(pairs)-n_graded} "
          f"UNGRADED (no lab-gf pool measured) ===")
    print(f"{'band':<13}{'ion':<4}{'engine':<16}{'pool':<10}"
          f"{'A':>8}{'stat':>8}{'syst':>8}{'total':>8}   vs ungraded")
    for p in pairs:
        prim = p.primary
        if prim is None:
            continue
        pool = "GRADED" if p.primary_is_graded else "ungraded"
        beats = p.graded_beats_ungraded
        cmp_txt = ""
        if beats is not None:
            u = p.ungraded["total_dex"]
            g = p.graded["total_dex"]
            cmp_txt = (f"  {u:.4f} -> {g:.4f}  "
                       f"{'TIGHTER' if beats else '⚠️ WIDER — do not promote'}")
        print(f"{p.band:<13}{p.ion:<4}{p.base:<16}{pool:<10}"
              f"{prim['A']:>8.3f}{prim['stat_dex']:>8.4f}{prim['syst_dex']:>8.4f}"
              f"{prim['total_dex']:>8.4f}{cmp_txt}")

    # ── the main-table view ───────────────────────────────────────────────────────
    et = element_table(combined)
    print(f"\n=== main Solar table view (graded where it exists, +/- VISIBLE) ===")
    for _, r in et[et.primary_is_graded].iterrows():
        print(f"  A({r.element} {r.ion}; {r.band}, {r.engine}) = "
              f"{format_value(r.A, r.total_dex)}   "
              f"[stat {r.stat_dex:.3f}, sys {r.syst_dex:.3f}, n={int(r.n_lines)}]")

    # ── what the cited term cost, against the bound it replaced ──────────────────
    # Reported the way round it is now DECIDED: the cited sigma sets the bar, and the
    # generic 0.041 is shown as the comparison so the widening is visible rather than
    # quietly absorbed. Every graded cell appears -- a cell that did NOT move would mean
    # the recharge missed it, so the empty case is as informative as the populated one.
    print(f"\n=== the published bar is the CITED lab sigma; the generic 0.041 for scale ===")
    for r in recharged:
        print(f"  {r['band']:<12} generic 0.041 -> total {r['total_generic']:.4f}   |   "
              f"cited -> total {r['total_cited']:.4f}   "
              f"({r['total_cited']/r['total_generic']-1:+.1%})")
    if not recharged:
        print("  (none moved — the cited term did not reach any graded cell; investigate)")

    # ── do the graded bands agree? ────────────────────────────────────────────────
    # If they disagree by more than their bars, "the graded value" is not well defined and
    # picking one is a choice rather than a measurement. RYA-712 keeps them separate; this
    # says whether that separation is hiding a conflict.
    g = [p for p in pairs if p.primary_is_graded]
    print(f"\n=== do the {len(g)} graded cells agree? (pairwise, in sigma) ===")
    consistency = []
    for i in range(len(g)):
        for j in range(i + 1, len(g)):
            a, b_ = g[i].graded, g[j].graded
            d = float(a["A"]) - float(b_["A"])
            sig = float(np.hypot(a["total_dex"], b_["total_dex"]))
            n_sig = abs(d) / sig if sig > 0 else float("nan")
            consistency.append({"a": f"{g[i].band}", "b": f"{g[j].band}",
                                "delta": d, "combined_sigma": sig, "n_sigma": n_sig})
            print(f"  {g[i].band:<12} vs {g[j].band:<12} "
                  f"delta {d:+.3f}  combined sigma {sig:.3f}  => {n_sig:.1f} sigma"
                  f"{'   ⚠️ DISCREPANT' if n_sig > 2 else ''}")
    spread = (max(float(x.graded["A"]) for x in g)
              - min(float(x.graded["A"]) for x in g)) if g else float("nan")
    print(f"  band-to-band spread of the graded values: {spread:.3f} dex")
    print(f"  (the 1D reference is 7.516 — RYA-669/783; gold's 7.466 is on the 3D scale)")

    combined.to_csv(OUT / "rya850_two_product_matrix.csv", index=False)
    et.to_csv(OUT / "rya850_element_table.csv", index=False)
    summary = {
        "ticket": "RYA-850",
        "n_cells": len(pairs),
        "n_graded_primary": n_graded,
        "n_ungraded_fallback": len(pairs) - n_graded,
        "graded_cells": [
            {"band": p.band, "ion": p.ion, "engine": p.base,
             "A": p.graded["A"], "n_lines": int(p.graded["n_lines"]),
             "stat_dex": p.graded["stat_dex"], "syst_dex": p.graded["syst_dex"],
             "total_dex": p.graded["total_dex"],
             "ungraded_total_dex": (p.ungraded or {}).get("total_dex"),
             "graded_beats_ungraded": p.graded_beats_ungraded}
            for p in pairs if p.primary_is_graded],
        "gf_term_published": "cited",
        "cited_gf_sigma_by_band": cited,
        "cells_recharged_to_cited": recharged,
        "graded_consistency": consistency,
        "graded_band_spread_dex": spread,
        "caveats": [
            "line sets are PRE-RYA-847; the graded pool is lab-gf INTERSECT constrained "
            "once 847 lands, so these numbers are provisional",
            "no combined cross-band headline is computed — RYA-712 forbids it and "
            "band_products has no combine(); the headline is a SELECTED cell",
            "the generic graded term 0.041 understates RYA-824's own cited lab sigma "
            "(0.0524 IR / 0.0600 VIS)",
            "15 of 18 cells have NO lab-gf pool measured and fall back to ungraded",
        ],
    }
    (OUT / "rya850_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\n[out] {OUT}/rya850_two_product_matrix.csv"
          f"\n[out] {OUT}/rya850_element_table.csv"
          f"\n[out] {OUT}/rya850_summary.json")


if __name__ == "__main__":
    main()
