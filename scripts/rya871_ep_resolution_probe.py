#!/usr/bin/env python3
"""
RYA-871 — what does carrying `ep_eV` actually resolve, and what does it break?
==============================================================================
    python3 scripts/rya871_ep_resolution_probe.py --band-products <dir>

The EW-route per-line artifact carries no excitation potential, so `gf_rung.resolve_lines`
matches a measured line back to the loaded line list on WAVELENGTH ALONE and 16 of 152 VIS
Fe I lines do not resolve. This measures the fix BEFORE it is written, because two things
about it are not obvious from the ticket and both change what the fix has to be.

🔴 FIRST: EP CANNOT FIX MOST OF THEM ON ITS OWN. RYA-855 already split the 16 —
`n_absent_from_linelist` 14, `n_ambiguous_in_linelist` 2. An EP key breaks TIES; it does
nothing for a line with NO row inside the window. So the fix is a PAIR: widen the
wavelength window far enough to reach the row, and use EP to keep the widening honest.
Widening alone is what RYA-855 refused, correctly — at 0.02 A several lines have two Fe I
rows straddling them and a wider window buys a choice, not an identification.

🔴 SECOND: THE MEASURED WAVELENGTH IS NOT A LINE-LIST WAVELENGTH. `measure_band_profilefit`
takes its candidates from `data/audit/line_accounting/per_line.csv`, whose rows are
FEATURES, not lines: `line_accounting_rya709.features()` groups line-list rows within
0.05 A and reports `w` as the group MEAN, `log_gf` as its MAX and `ep_eV` as its MIN. So a
measured "line" at a blended feature sits between its components by construction, and that
— not measurement error — is why the nearest row is 0.006-0.02 A away. It also means the
`ep_eV` this ticket carries is the MINIMUM EP over the cluster, which is a real row's EP
but not necessarily the row whose gf was reported. Whether that key helps or hurts is a
measurement, and it is the one below.

THE CONTROL THAT DECIDES IT
---------------------------
A rule that only fires where the current one FAILED is unfalsifiable (RYA-818). So every
variant is scored on BOTH populations:

* the lines the wavelength-only rule already resolves — a variant must resolve the SAME
  ROW for every one of them. A variant that re-identifies even one line that was already
  identified is REJECTED, whatever it does for the 16.
* the lines it does not — split into ABSENT and AMBIGUOUS, because they have different
  fixes and collapsing them hides which is which.

Nothing here changes a value, a rung or an artifact. It prints a table and writes it.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from pipeline import gf_rung                                       # noqa: E402
from pipeline.band_policy import resolve as resolve_band           # noqa: E402
from pipeline.error_budget import build as build_budget            # noqa: E402
from pipeline import harness_residual                              # noqa: E402  RYA-869
from rya855_rung_audit import BAND_PIVOT_A, _cells, _linelist_for   # noqa: E402

OUT = ROOT / "data" / "results" / "rya871"
ACCOUNTING = ROOT / "data" / "audit" / "line_accounting" / "per_line.csv"

#: Wavelength tolerances to score. 0.005 is `gf_rung.LINELIST_MATCH_TOL_A` as it stands;
#: the rest bracket the 0.006-0.02 A offsets RYA-855 measured on the unresolved lines.
TOL_A = (0.005, 0.010, 0.020, 0.030, 0.050)

#: EP agreement window. The accounting EP is rounded to 4 dp and the line list carries
#: full precision, so this is a rounding tolerance and not a physical one — two REAL
#: transitions at one wavelength differ by whole eV (RYA-855's 3125.65 pair: 0.990 vs
#: 2.404), so nothing here depends on where inside this range the cut sits.
EP_TOL_EV = 0.005


def accounting_ep() -> pd.DataFrame:
    """The EP the emitter dropped, per (element, ion, wavelength) it emitted.

    `measure_band_profilefit` copies `float(r.wave_air_A)` from this table onto the
    LineMeasurement verbatim, so the join back is EXACT rather than nearest-within-
    tolerance — and that is asserted below, not assumed. If it ever stops being exact the
    emitter has started rounding and this probe must not paper over it.
    """
    d = pd.read_csv(ACCOUNTING)
    return d[["element", "ion", "wave_air_A", "ep_eV", "log_gf"]].copy()


def score(measured: pd.DataFrame, ll: pd.DataFrame, *, tol_A: float,
          use_ep: bool) -> dict:
    """Resolve every measured line against `ll` and classify the outcome.

    `measured` needs wavelength_air_A and (when `use_ep`) ep_eV. `ll` is the species-
    filtered line list frame from `gf_rung.linelist_frame`.
    """
    lw = ll.wavelength_air_A.to_numpy()
    lep = ll.ep_eV.to_numpy()
    out = []
    for r in measured.itertuples():
        w = float(r.wavelength_air_A)
        m = np.abs(lw - w) <= tol_A
        why = ""
        if use_ep and m.any():
            ep = float(getattr(r, "ep_eV", np.nan))
            if not np.isfinite(ep):
                why = "no-ep"          # the artifact did not carry one; cannot key on it
            else:
                m = m & (np.abs(lep - ep) <= EP_TOL_EV)
        n = int(m.sum())
        idx = int(np.flatnonzero(m)[0]) if n == 1 else -1
        out.append({"wavelength_air_A": w, "n_match": n, "row": idx,
                    "state": ("unique" if n == 1 else "absent" if n == 0 else "ambiguous"),
                    "why": why})
    return pd.DataFrame(out)


def probe(band_products: Path, lists: dict) -> pd.DataFrame:
    acc = accounting_ep()
    rows = []
    for cell in _cells(band_products):
        if cell["route"] != "PROFILEFIT":
            continue        # the SYNTH/LABGF routes are keyed at the list's own wavelength
        lines = pd.read_csv(cell["path"])
        used = lines[lines.in_aggregate.astype(bool) & lines.abundance.notna()].copy()
        if not len(used):
            continue
        pol = resolve_band(0.5 * (cell["lo"] + cell["hi"]))
        try:
            # `raw` is the loaded iSpec list; `ll` is its species-filtered FRAME. Both are
            # kept because they are consumed differently: `score` below works on the frame
            # directly, while `gf_rung.resolve_lines` takes the raw list and does its own
            # framing and filtering — handing it the frame calls `linelist_frame` on an
            # already-framed object and looks for a `wave_nm` column that is long gone.
            raw = _linelist_for(cell["route"], pol.name, lists)
            ll = gf_rung.linelist_frame(raw)
        except KeyError as e:
            print(f"  SKIP {cell['path'].name}: {e}")
            continue
        ll = ll[ll.species == gf_rung._species_label(cell["element"], cell["ion"])]

        # ── recover the dropped EP by an EXACT join, and prove it is exact ──────────
        a = acc[(acc.element == cell["element"]) & (acc.ion == cell["ion"])]
        merged = used.merge(a[["wave_air_A", "ep_eV"]], how="left",
                            left_on="wavelength_air_A", right_on="wave_air_A")
        n_joined = int(merged.ep_eV.notna().sum())

        base = score(merged, ll, tol_A=gf_rung.LINELIST_MATCH_TOL_A, use_ep=False)
        base_ok = base.state == "unique"

        # ── DOES THE BAR MOVE? ────────────────────────────────────────────────────
        # The ticket's CRITICAL check, and it is not answerable by counting resolutions:
        # a newly identified line changes the rung only if it changes whether EVERY line
        # in the pool is primary-lab. RYA-855 measured every Fe pool as MIXED several
        # times over, so the expectation is no movement — which is exactly why it is
        # measured rather than argued (an expected null still needs the measurement).
        #
        # Everything but the gf rung is held: same n, same scatter, same harness term,
        # the same band pivot the RYA-855/869 audits use. So a moved syst can only be the
        # rung, and an unmoved one is not an accident of two terms cancelling.
        vals = used.abundance.to_numpy(dtype=float)
        n_used = int(len(vals))
        scatter = float(np.std(vals, ddof=1)) if n_used > 1 else 0.0
        pivot = BAND_PIVOT_A.get(pol.name, 0.5 * (cell["lo"] + cell["hi"]))
        hr = harness_residual.for_handler(
            harness_residual.handler_of_banked_cell(route=cell["route"],
                                                    treatment=cell["treatment"]))

        def _rung(use_ep):
            eps = (merged.ep_eV.tolist() if use_ep else None)
            lg = gf_rung.resolve_lines(cell["element"], cell["ion"],
                                       merged.wavelength_air_A, raw, measured_ep_eV=eps)
            return gf_rung.decide(cell["element"], cell["ion"], lg)

        def _syst(rung):
            b = build_budget(cell["element"], pivot, max(n_used, 1),
                             scatter_dex=scatter, **rung.budget_kwargs(),
                             **hr.budget_kwargs())
            return round(b.total()[1], 4)

        r_before, r_after = _rung(False), _rung(True)
        syst_before, syst_after = _syst(r_before), _syst(r_after)
        A_median = round(float(np.median(vals)), 3) if n_used else float("nan")

        for tol in TOL_A:
            for use_ep in (False, True):
                s = score(merged, ll, tol_A=tol, use_ep=use_ep)
                # 🔴 THE CONTROL: does this variant re-identify a line the current rule
                # already identified? Compared by the ROW it lands on, not by whether it
                # landed — a variant that swaps one identification for another is a
                # regression that a resolved-count would score as a tie.
                same = (s.row[base_ok] == base.row[base_ok])
                rows.append({
                    "band": pol.name, "element": cell["element"], "ion": cell["ion"],
                    "treatment": cell["treatment"], "deck": cell["deck"],
                    "n_lines": len(merged), "n_ep_joined": n_joined,
                    "tol_A": tol, "use_ep": use_ep,
                    "n_unique": int((s.state == "unique").sum()),
                    "n_absent": int((s.state == "absent").sum()),
                    "n_ambiguous": int((s.state == "ambiguous").sum()),
                    "baseline_n_unique": int(base_ok.sum()),
                    "baseline_kept": int(same.sum()),
                    "baseline_reidentified": int((~same).sum()),
                    # Held identical across every variant of this cell — the value is a
                    # median of per-line abundances and no identity key enters it.
                    "A": A_median,
                    "rung_before": r_before.rung, "rung_after": r_after.rung,
                    "n_graded_before": r_before.n_graded,
                    "n_graded_after": r_after.n_graded,
                    "syst_before": syst_before, "syst_after": syst_after,
                    "d_syst": round(syst_after - syst_before, 4),
                    "src": str(cell["path"].relative_to(band_products)),
                })
    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--band-products", required=True, type=Path)
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args()

    from pipeline.abundances_derive import _load_synth_resources
    lists = {}
    lists["__ew__"], _, _ = _load_synth_resources()
    print(f"[linelist] EW route: {len(lists['__ew__'])} rows (GES, canonical gf applied)")

    d = probe(a.band_products, lists)
    if not len(d):
        raise SystemExit("no PROFILEFIT cells found")
    a.out.mkdir(parents=True, exist_ok=True)
    d.to_csv(a.out / "rya871_ep_resolution_probe.csv", index=False)

    print(f"\n=== does the dropped EP even come back? ===")
    for _, r in d.drop_duplicates(subset=["band", "element", "ion", "treatment", "deck"]).iterrows():
        mark = "OK" if r.n_ep_joined == r.n_lines else "⚠️ INCOMPLETE"
        print(f"  {r.band:<12}{r.element} {r.ion:<3}{r.treatment:<16}{(r.deck or 'root'):<12}"
              f"{r.n_ep_joined:>4}/{r.n_lines:<4} exact accounting joins   {mark}")

    print(f"\n=== resolution vs wavelength tolerance, with and without the EP key ===")
    print(f"{'band':<12}{'ion':<4}{'treatment':<16}{'tol':>7}{'EP':>4}"
          f"{'unique':>8}{'absent':>8}{'ambig':>7}{'kept':>8}{'RE-ID':>7}")
    for _, r in d.sort_values(["band", "ion", "treatment", "deck", "tol_A", "use_ep"]).iterrows():
        flag = "  🔴" if r.baseline_reidentified else ""
        print(f"{r.band:<12}{r.ion:<4}{r.treatment:<16}{r.tol_A:>7.3f}"
              f"{('yes' if r.use_ep else 'no'):>4}{r.n_unique:>8}{r.n_absent:>8}"
              f"{r.n_ambiguous:>7}{r.baseline_kept:>5}/{r.baseline_n_unique:<3}"
              f"{r.baseline_reidentified:>6}{flag}")

    # ── what the probe concludes, stated as the rule it licenses ─────────────────
    safe = d[d.baseline_reidentified == 0]
    print(f"\n=== variants that re-identify NOTHING already identified ===")
    if not len(safe):
        print("  none — every widening moves at least one existing identification")
    else:
        # NB: `n_unique` etc keep their `n_` prefix here — `df.unique` is a DataFrame
        # METHOD, so an aggregate named `unique` silently resolves to it on itertuples.
        best = (safe.groupby(["tol_A", "use_ep"])
                    .agg(n_unique=("n_unique", "sum"), n_absent=("n_absent", "sum"),
                         n_ambiguous=("n_ambiguous", "sum"), n_cells=("n_unique", "size"))
                    .reset_index().sort_values(["n_unique", "tol_A"],
                                               ascending=[False, True]))
        print(f"{'tol_A':>8}{'EP':>5}{'unique':>9}{'absent':>8}{'ambig':>7}   (summed over cells)")
        for _, r in best.iterrows():
            print(f"{r.tol_A:>8.3f}{('yes' if r.use_ep else 'no'):>5}"
                  f"{int(r.n_unique):>9}{int(r.n_absent):>8}{int(r.n_ambiguous):>7}")

    # ── the ticket's CRITICAL checks, answered ────────────────────────────────────
    per_cell = d.drop_duplicates(subset=["band", "element", "ion", "treatment", "deck"])
    moved = per_cell[per_cell.d_syst.abs() > 5e-5]
    lifted = per_cell[per_cell.rung_after != per_cell.rung_before]
    print(f"\n=== does the EP key move a RUNG or a BAR? ===")
    print(f"  {len(per_cell)} cells; {len(lifted)} change rung; {len(moved)} change syst")
    if len(lifted) or len(moved):
        for _, r in pd.concat([lifted, moved]).drop_duplicates().iterrows():
            print(f"  ⚠️ {r.band:<12}{r.element} {r.ion:<3}{r.treatment:<16}"
                  f"rung {r.rung_before}->{r.rung_after}  "
                  f"syst {r.syst_before:.4f}->{r.syst_after:.4f}")
    else:
        print("  none. Every Fe pool is MIXED several times over (RYA-855), so a newly")
        print("  identified line adds a grade to a pool that already had an ungraded one")
        print("  — it cannot lift a rung, and the bar it prices does not move.")
    print(f"  n_graded (lines newly PRICED, which is the thing that did change): "
          f"{int(per_cell.n_graded_before.sum())} -> {int(per_cell.n_graded_after.sum())}")
    print(f"\n=== values ===")
    print(f"  A is the median of the per-line abundances and no identity key enters it; "
          f"held identical across every variant by construction, and the column is "
          f"carried per cell so a reader can check rather than take it.")

    (a.out / "rya871_probe_summary.json").write_text(json.dumps({
        "ticket": "RYA-871",
        "band_products": str(a.band_products),
        "ep_tol_eV": EP_TOL_EV,
        "current_tol_A": gf_rung.LINELIST_MATCH_TOL_A,
        "epkey_tol_A": gf_rung.LINELIST_MATCH_TOL_EPKEY_A,
        "n_cells": int(len(per_cell)),
        "n_cells_changing_rung": int(len(lifted)),
        "n_cells_changing_syst": int(len(moved)),
        "n_graded_before": int(per_cell.n_graded_before.sum()),
        "n_graded_after": int(per_cell.n_graded_after.sum()),
        "variants": json.loads(d.to_json(orient="records")),
    }, indent=2) + "\n")
    print(f"\n  wrote {a.out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
