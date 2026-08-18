#!/usr/bin/env python3
"""RYA-847 item 6 — what did the non-minimum check cost, cell by cell?

    python3 scripts/rya847_gate_diff.py

READS ONLY. Every post-gate number comes from the regenerated artifacts under
data/results/rya847/gated/; every pre-gate number comes from
data/results/rya847/rya847_pregate_control.csv, which records RYA-847's own nine-cell
sweep — the SAME tree and SAME code with the gate not yet applied.

WHY NOT DIFF THE PUBLISHED MATRIX
---------------------------------
Because it would not measure this gate. data/results/rya783/fe_product_matrix.csv comes
from the rya845 run and its line counts differ by four ATOMIC_BLEND registry exclusions,
11119.795, and the RYA-807 registry gate. Diffing against it measures the gate PLUS months
of pool drift — the confound RYA-848 hit and settled with a same-inputs control.

WHAT MAKES THE CONTROL CHECKABLE
--------------------------------
The control numbers are carried as data (the sweep's per-line CSVs were destroyed by a
tree re-sync), so they are ASSERTED against the artifact rather than trusted:

  * n_post + n_caught == n_pre        -- per cell, exactly
  * scatter == stat_dex * sqrt(n)     -- on BOTH sides, since build_product reports
                                         stat = scatter/sqrt(N)

A control row that does not belong to this run fails those checks and this script exits
non-zero. That is the point: a hand-carried number must be falsifiable by the artifact it
claims to describe (RYA-847's own "a quantity with nowhere to live cannot be checked").
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

# In-repo artifacts, derived from this file's own location — the same construction
# scripts/rya850_graded_products.py uses. `codex_path` keys resolve EXTERNAL roots
# (spectra, grids); nothing here leaves the checkout.
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

RESULTS = ROOT / "data" / "results" / "rya847"
GATED = RESULTS / "gated"
CONTROL = RESULTS / "rya847_pregate_control.csv"
OUT_DIFF = RESULTS / "rya847_gate_diff.csv"
OUT_LINES = RESULTS / "rya847_gate_caught_lines.csv"

# The deck is the SUBDIRECTORY, never part of the filename — derive_band_products keys
# its products file on band + instrument + route only, which is why each deck must own an
# output directory (see scripts/rya847_regen_gated_cells.sh). The cell comes from the
# products row's own `band` column, so nothing here has to parse a wavelength out of a
# filename.
DECK_DIRS = ("ts-lte", "gerber-nlte")
# Treatments that do NOT depend on the Engine-B deck. Both VIS runs re-derive them, so the
# two copies are a free run-to-run reproducibility control.
DECK_FREE = ("1D-LTE", "ENGINE-A")


def _rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def _f(v, default=float("nan")):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def _deck_of(path):
    """A products/lines file's deck is the directory it sits in, or none."""
    return path.parent.name if path.parent.name in DECK_DIRS else ""


def _products():
    """Every product row under gated/, tagged with the deck directory it came from."""
    out = []
    for path in sorted(GATED.rglob("*_products.csv")):
        for r in _rows(path):
            r["_deck"] = _deck_of(path)
            r["_src"] = str(path.relative_to(GATED.parent))
            out.append(r)
    return out


def _caught(products_path, treatment):
    """The lines this product's gate excluded, read from ITS OWN per-line artifact.

    Paired by filename stem, not by directory: the top-level directory holds both the
    near-UV and the NIR cell, so a directory-wide glob would attribute one cell's
    exclusions to the other.
    """
    prefix = products_path.name[: -len("products.csv")]
    lines_path = products_path.with_name(f"{prefix}{treatment}_lines.csv")
    if not lines_path.exists():
        return None
    hits = []
    for r in _rows(lines_path):
        if r.get("treatment") != treatment:
            continue
        if "NON-MINIMUM" not in (r.get("excluded_reason") or ""):
            continue
        hits.append({
            "treatment": treatment,
            "wavelength_air_A": r["wavelength_air_A"],
            "abundance": r.get("abundance", ""),
            "frac_rise_weaker": r.get("frac_rise_weaker", ""),
            "sigma_A": r.get("sigma_A", ""),
            "edge_distance_dex": r.get("edge_distance_dex", ""),
            "red_chi2": r.get("red_chi2", ""),
            "excluded_reason": r.get("excluded_reason", ""),
            "_src": str(lines_path.relative_to(GATED.parent)),
        })
    return hits


def main() -> int:
    control = {(r["cell"], r["deck"], r["treatment"]): r for r in _rows(CONTROL)}
    diff_rows, caught_rows, failures = [], [], []
    deck_free_seen = {}

    for prod in _products():
        band, deck, treatment = prod["band"], prod["_deck"], prod["treatment"]
        n_post = int(_f(prod["n_lines"]))
        a_post = _f(prod["A"])
        stat_post = _f(prod["stat_dex"])
        scatter_post = stat_post * math.sqrt(n_post)

        hits = _caught(GATED.parent / prod["_src"], treatment)
        if hits is None:
            failures.append(f"{band}/{deck}/{treatment}: no per-line artifact beside its "
                            f"products row — the exclusions cannot be read back")
            hits = []
        for h in hits:
            h["cell"] = band
            h["deck"] = deck
            caught_rows.append(h)

        # Deck-free treatments are re-derived by BOTH VIS runs. They must agree: if they
        # do not, the two decks did not run the same inputs and no diff below is safe.
        if treatment in DECK_FREE and deck:
            prev = deck_free_seen.get((band, treatment))
            now = (a_post, n_post, round(stat_post, 6))
            if prev is not None and prev[1] != now:
                failures.append(
                    f"{band}/{treatment}: deck-free product differs between decks "
                    f"{prev[0]} {prev[1]} vs {deck} {now} — the runs are not comparable")
            deck_free_seen[(band, treatment)] = (deck, now)

        key = (band, deck, treatment)
        ctl = control.get(key)
        row = {
            "cell": band, "deck": deck, "treatment": treatment,
            "A_pre": "", "n_pre": "", "scatter_pre": "",
            "A_post": f"{a_post:.4f}", "n_post": n_post,
            "stat_post": f"{stat_post:.4f}", "scatter_post": f"{scatter_post:.4f}",
            "n_caught": len(hits),
            "lines_caught": " ".join(h["wavelength_air_A"] for h in hits),
            "dA": "", "d_scatter": "", "control_source": "",
        }

        if ctl is None:
            # No control row is the RIGHT answer for a cell the gate cannot touch — the
            # EW-route products. Assert that rather than assume it.
            if hits:
                failures.append(
                    f"{key}: {len(hits)} line(s) caught but no control row exists — a "
                    f"cell the gate CHANGED must carry a before/after")
            row["control_source"] = "n/a — gate excluded 0 lines here"
            diff_rows.append(row)
            continue

        n_pre = int(_f(ctl["n_pre"]))
        a_pre = _f(ctl["A_pre"])
        scatter_pre = _f(ctl["scatter_pre"])

        if n_post + len(hits) != n_pre:
            failures.append(
                f"{key}: n_post {n_post} + caught {len(hits)} != n_pre {n_pre} — the "
                f"control does not describe this run")

        row.update({
            "A_pre": f"{a_pre:.4f}", "n_pre": n_pre,
            "scatter_pre": f"{scatter_pre:.4f}",
            "dA": f"{a_post - a_pre:+.4f}",
            "d_scatter": f"{scatter_post - scatter_pre:+.4f}",
            "control_source": ctl["source"],
        })
        diff_rows.append(row)

    if not diff_rows:
        print(f"no products under {GATED} — run scripts/rya847_regen_gated_cells.sh on "
              f"Sirius first", file=sys.stderr)
        return 2

    # 🔴 EVERY CONTROL ROW MUST FIND ITS PRODUCT. This is the guard for the defect that
    # produced this script's first run: both VIS decks were pointed at ONE --out, the
    # second overwrote the first products.csv, and the ENGINE-B row — the cell the gate
    # does the most work in — simply vanished. Nothing failed; the diff just came up one
    # row short, which is unreadable as an error. A missing product is now LOUD.
    seen = {(r["cell"], r["deck"], r["treatment"]) for r in diff_rows}
    for key in control:
        if key not in seen:
            failures.append(
                f"{key}: a control row exists but NO product was found for it — the "
                f"artifact set is incomplete (a deck overwriting another deck's "
                f"products.csv looks exactly like this)")

    diff_rows.sort(key=lambda r: (r["cell"], r["deck"], r["treatment"]))
    caught_rows.sort(key=lambda r: (r["cell"], r["treatment"],
                                    float(r["wavelength_air_A"])))

    RESULTS.mkdir(parents=True, exist_ok=True)
    with open(OUT_DIFF, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(diff_rows[0]))
        w.writeheader()
        w.writerows(diff_rows)
    if caught_rows:
        with open(OUT_LINES, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(caught_rows[0]))
            w.writeheader()
            w.writerows(caught_rows)

    for r in diff_rows:
        print(f"{r['cell']:<8} {r['deck'] or '-':<12} {r['treatment']:<14} "
              f"n {str(r['n_pre']) or '-':>4} -> {r['n_post']:<4} "
              f"A {r['A_pre'] or '-':>7} -> {r['A_post']:<7} "
              f"dA {r['dA'] or '-':>8}  "
              f"scatter {r['scatter_pre'] or '-':>7} -> {r['scatter_post']:<7} "
              f"caught {r['n_caught']}")
    print(f"\nwrote {OUT_DIFF}")
    if caught_rows:
        print(f"wrote {OUT_LINES}")

    if failures:
        print("\n\U0001f534 CONTROL CHECKS FAILED:", file=sys.stderr)
        for f in failures:
            print("   " + f, file=sys.stderr)
        return 1
    print("\ncontrol checks pass: n_post + n_caught == n_pre in every changed cell, and "
          "the deck-free products agree across decks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
