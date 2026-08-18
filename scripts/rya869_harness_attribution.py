#!/usr/bin/env python3
"""
RYA-869 — which cells were charged another handler's systematic, and what moves
==============================================================================
    python3 scripts/rya869_harness_attribution.py --band-products <dir>

`scripts/derive_band_products.py` decided the harness residual with
`is_b = prod.treatment == "ENGINE-B"`, an equality against a treatment NAME under a
comment that stated the HANDLER rule correctly. RYA-798 added the variant
`ENGINE-B-NLTE` — the same `SynthesisHandler` flux fit on a different departure deck —
and the equality never followed it, so every NLTE product was charged the profile
fitter's 0.0129 dex and labelled `ProfileFitHandler` in its own budget file.

This measures the consequence on the ALREADY-PUBLISHED cells, and it is the diff the
fix has to be judged on.

WHY THIS RE-DERIVES THE BUDGET INSTEAD OF RE-RUNNING THE DERIVER
----------------------------------------------------------------
Same reason as RYA-855, from which this borrows its cell walk verbatim rather than
rebuilding it: the budget is a pure function of (band, n_lines, scatter, gf rung, harness
residual), and re-fitting the cells — days of synthesis — could not change any term but
the one being changed. What it WOULD do is fold two months of input drift into a diff
that must isolate a single term (RYA-848).

THE TWO CONTROLS, ONE AT EACH END
---------------------------------
A null result needs a control at BOTH ends (RYA-855). So:

1. **The baseline reproduces the published bar.** Each cell's budget is rebuilt from its
   own per-line file under the PRE-FIX rule and must equal the published `stat_dex` and
   `syst_dex` exactly. A cell that fails is reported and NOT diffed — its published bar
   was not built from the inputs standing here, so a diff would measure the drift.
2. **The cells that must NOT move, do not.** The 30 cells whose handler the old rule
   already got right must come out BYTE-IDENTICAL in the full budget text, and the 6 that
   move must differ in the harness line AND NOTHING ELSE. Reported as a count of changed
   budget lines per cell, so "only the harness moved" is measured rather than asserted.

⚠️ THE VALUES ARE NOT TOUCHED AND CANNOT BE. `A` is the median of the per-line
abundances; no term of the budget enters it. Asserted per cell anyway rather than argued.

🔴 A SECOND FINDING, REPORTED AND NOT FIXED HERE. `SynthesisHandler` is charged 0.0000
while its own banked control measured a −0.0100 dex offset and PASSED
(`data/audit/synthesis_control/control_FeI.json`, n=18, RYA-770/759) — and
`error_budget.harness_term` prints the words "MEASURED against the known optical answer,
not assumed zero" beside that 0.0000. The prose and the arithmetic disagree. Charging the
measured value moves the near-UV cell and every ENGINE-B cell, which is a different diff
from this one, so it is quantified below and left alone. (RYA-875 later resolved it: the
0.0100 was a LINE-SET artifact and the real offset is 0.0000, so nothing is charged
differently and `harness_residual_divergences` now reports agreement rather than a gap.)
"""
from __future__ import annotations

import argparse
import difflib
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from config.constants import codex_root                            # noqa: E402
from pipeline import harness_residual                              # noqa: E402
from pipeline.band_policy import resolve as resolve_band           # noqa: E402
from pipeline.error_budget import build as build_budget            # noqa: E402
# Borrowed, not rebuilt (RYA-832: call the functions, never reconstruct them — a rebuild
# had already drifted by the time it was compared). `_cells` owns the band-product stem
# grammar and `_published` owns "which products.csv row is this cell"; a second copy of
# either would be free to disagree with RYA-855's audit about what a cell even is.
from rya855_rung_audit import BAND_PIVOT_A, _cells, _published     # noqa: E402

OUT = ROOT / "data" / "results" / "rya869"

#: The banked band-product tree the RYA-783/807/832 matrix was assembled from — the same
#: default RYA-855 audited. Named by KEY, never pasted (RYA-810: `audit_path_literals.py`
#: is a hard gate, and a literal is how ~60 of them got here).
DEFAULT_BP = codex_root("work") / "rya845" / "data" / "results" / "band_products"


def published_rule_harness(treatment: str, route: str) -> tuple[float, str]:
    """The residual and label `derive_band_products` charged BEFORE RYA-869.

    Thin wrapper over `harness_residual.charged_in_banked_cell`, which is the ONE home
    for the pre-fix rule — `scripts/rya855_rung_audit.py` needs the same answer, to tell
    "this published cell does not reproduce because its inputs drifted" apart from "this
    published cell does not reproduce because RYA-869 corrected it".
    """
    h = harness_residual.charged_in_banked_cell(route=route, treatment=treatment)
    return h.residual_dex, h.handler


def audit(band_products: Path) -> tuple[pd.DataFrame, dict]:
    rows, texts = [], {}
    for cell in _cells(band_products):
        lines = pd.read_csv(cell["path"])
        used = lines[lines.in_aggregate.astype(bool) & lines.abundance.notna()]
        n = int(len(used))
        vals = used.abundance.to_numpy(dtype=float)
        value = float(np.median(vals)) if n else np.nan
        scatter = float(np.std(vals, ddof=1)) if n > 1 else 0.0

        pol = resolve_band(0.5 * (cell["lo"] + cell["hi"]))
        pivot = BAND_PIVOT_A.get(pol.name, 0.5 * (cell["lo"] + cell["hi"]))

        # RYA-855 measured every published cell onto gf rung 1 (0 of 36 bars move on the
        # gf term), and rung 1 is what the deriver hardcoded, so holding gf_graded=False
        # here is not a simplification: it is the published gf term, and the reproduction
        # control below is what proves it.
        def _budget(residual, handler, provenance):
            return build_budget(cell["element"], pivot, max(n, 1), scatter_dex=scatter,
                                gf_graded=False, harness_residual_dex=residual,
                                handler=handler, harness_provenance=provenance)

        h_pub, l_pub = published_rule_harness(cell["treatment"], cell["route"])
        fixed = harness_residual.for_handler(
            harness_residual.handler_of_banked_cell(route=cell["route"],
                                                    treatment=cell["treatment"]))
        # RYA-873 — BOTH sides are priced with TODAY's provenance for the handler each
        # side names. This audit isolates the ATTRIBUTION change (which handler), so
        # every other input is held, prose included (RYA-848). Giving the "before" side
        # the old fixed MEASURED string instead would make RYA-873's prose fix show up
        # here as if it were part of RYA-869's diff, and the two must stay separable —
        # measured: it turned 30 byte-identical budgets into 24.
        before = _budget(h_pub, l_pub,
                         harness_residual.for_handler(l_pub).provenance)
        after = _budget(fixed.residual_dex, fixed.handler, fixed.provenance)
        stat_b, syst_b = before.total()
        stat_a, syst_a = after.total()

        # WHICH lines of the full budget text changed, not how many. A moved cell must
        # change the harness TERM line and the `systematic` TOTAL that sums it, and
        # nothing else — the scatter, the gf term, the pseudo-continuum, the telluric
        # residual, the statistical total and the dominant-term verdict all have to come
        # out identical. Counting would have let a second term move as long as one term
        # stopped moving; classifying the changed lines cannot.
        tb, ta = before.describe(), after.describe()
        changed = [x for x, y in zip(tb.splitlines(), ta.splitlines()) if x != y]
        n_changed = len(changed)
        n_unexpected = sum(1 for x in changed
                           if "harness residual" not in x
                           and not x.strip().startswith("systematic"))
        texts[str(cell["path"].relative_to(band_products))] = (tb, ta)

        pub = _published(band_products, cell)
        repro = dict(stat=None, syst=None, A=None)
        if pub is not None:
            repro = dict(
                stat=bool(abs(round(stat_b, 4) - float(pub["stat_dex"])) <= 5e-5),
                syst=bool(abs(round(syst_b, 4) - float(pub["syst_dex"])) <= 5e-5),
                A=bool(abs(round(value, 3) - float(pub["A"])) <= 5e-4))

        rows.append({
            "band": pol.name, "element": cell["element"], "ion": cell["ion"],
            "treatment": cell["treatment"], "deck": cell["deck"], "route": cell["route"],
            "n_lines": n, "A": round(value, 3),
            "has_published_row": pub is not None,
            "A_published": (None if pub is None else float(pub["A"])),
            "reproduces_A": repro["A"], "reproduces_stat": repro["stat"],
            "reproduces_syst": repro["syst"],
            "handler_published": l_pub, "handler_correct": fixed.handler,
            "handler_misattributed": bool(l_pub != fixed.handler),
            "harness_published_dex": round(h_pub, 4),
            "harness_correct_dex": round(fixed.residual_dex, 4),
            "stat_dex": round(stat_b, 4),
            "stat_dex_after": round(stat_a, 4),
            "syst_before": round(syst_b, 4),
            "syst_published": (None if pub is None else float(pub["syst_dex"])),
            "syst_after": round(syst_a, 4),
            "d_syst": round(syst_a - syst_b, 4),
            "n_budget_lines_changed": n_changed,
            "n_budget_lines_changed_outside_harness": n_unexpected,
            "src": str(cell["path"].relative_to(band_products)),
        })
    return pd.DataFrame(rows), texts


def tree_fingerprint(band_products: Path) -> str:
    """SHA-256 over every per-line file's RELATIVE name and its bytes.

    The absolute path this ran against is provenance of the RUN, not of the DATA, and it
    changes with the machine — the RYA-845 lesson, where a generating worktree's absolute
    path on every matrix row buried the one cell that actually moved. So the tree is also
    identified by what is IN it. Two runs on different hosts that agree here saw the same
    inputs; two that do not, did not, whatever their paths say (compare FILES + BYTES,
    never names).
    """
    h = hashlib.sha256()
    for f in sorted(band_products.rglob("*_lines.csv"),
                    key=lambda x: str(x.relative_to(band_products))):
        h.update(str(f.relative_to(band_products)).encode())
        h.update(f.read_bytes())
    return h.hexdigest()


def harness_residual_divergences() -> dict:
    """Do any charged harness residuals disagree with their own banked control?

    🔴 THIS USED TO HARDCODE `UNCHARGED_CONTROL_RESIDUAL_DEX["SynthesisHandler"]`, which
    is a lookup that only works while the finding is unresolved — and it KeyError'd the
    moment RYA-875 resolved it. A reporter that cannot survive its own subject being
    fixed makes fixing the subject look like breaking the tool. So it reads the
    declaration generically and says what it finds, including nothing.

    RYA-869 filed the SynthesisHandler divergence; RYA-873 refused to charge either
    candidate number; RYA-875 found the -0.0100 was a LINE-SET artifact (a median over 18
    lines compared to a hardcoded scalar that is the median of a different 23-line set)
    and established the real offset at 0.0000, paired per line. The entry left the
    declaration BECAUSE THE NUMBERS AGREE.
    """
    declared = harness_residual.UNCHARGED_CONTROL_RESIDUAL_DEX
    out: dict = {"n_declared": len(declared), "declared": {}, "agreeing": {}}
    for handler, charged in harness_residual.HANDLER_RESIDUAL_DEX.items():
        art = (declared.get(handler, {}).get("control_artifact")
               or harness_residual.HANDLER_CONTROL_ARTIFACT.get(handler))
        d = {}
        if art and (ROOT / art).exists():
            d = json.loads((ROOT / art).read_text())
        rec = {"charged_dex": charged, "control_artifact": art,
               "control_present": bool(d), "control_dex_offset": d.get("dex_offset"),
               "control_passed": d.get("passed"), "control_n_lines": d.get("n_lines")}
        if handler in declared:
            rec["ticket"] = declared[handler].get("ticket", "")
            rec["why_not_charged"] = declared[handler].get("why_not_charged", "")
            out["declared"][handler] = rec
        elif art:
            rec["agrees"] = bool(d) and abs(charged - abs(float(
                d.get("dex_offset", float("nan"))))) < 5e-5
            out["agreeing"][handler] = rec
    out["status"] = (
        "no charged harness residual disagrees with its banked control. The "
        "SynthesisHandler divergence RYA-869 filed and RYA-873 declined to charge was "
        "resolved by RYA-875: its -0.0100 was a LINE-SET artifact (18-line median vs a "
        "scalar that is the median of a different 23-line set), and paired per line the "
        "offset is 0.0000 with 17 of 18 lines inside 0.01 dex. Deleted because the "
        "numbers agree."
        if not declared else
        f"{len(declared)} declared divergence(s) stand; see `declared`.")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--band-products", default=str(DEFAULT_BP))
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args()

    bp = Path(a.band_products)
    if not bp.exists():
        raise SystemExit(
            f"band products absent: {bp}\n"
            f"  They live on Sirius. Point --band-products at a copy, or run there.")
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)

    d, texts = audit(bp)
    if not len(d):
        raise SystemExit(f"no band-product cells under {bp}")

    # ── control 1: does the pre-fix rule reproduce what was published? ────────────
    have = d[d.has_published_row]
    bad = have[(have.reproduces_stat != True) | (have.reproduces_syst != True)  # noqa: E712
               | (have.reproduces_A != True)]                                   # noqa: E712
    print(f"=== control 1 — the PRE-FIX rule reproduces the published cell ===")
    print(f"  {len(have) - len(bad)}/{len(have)} cells reproduce A, stat and syst exactly")
    for _, r in bad.iterrows():
        print(f"  ⚠️ NOT REPRODUCED, excluded from the diff: {r.band} {r.element} {r.ion} "
              f"{r.treatment} ({r.deck or 'root'})  "
              f"A {r.A} vs {r.A_published}, syst {r.syst_before} vs {r.syst_published}")
    ok = d[~d.index.isin(bad.index)]

    # ── control 2: nothing but the harness term moved, at either end ──────────────
    moved = ok[ok.handler_misattributed]
    still = ok[~ok.handler_misattributed]
    print(f"\n=== control 2 — only the harness term moves, and only where it should ===")
    print(f"  {len(still)} cells whose handler the old rule already got right: "
          f"{int((still.n_budget_lines_changed == 0).sum())} byte-identical budgets, "
          f"{int((still.d_syst.abs() > 0).sum())} with a moved systematic")
    print(f"  {len(moved)} misattributed cells: {sorted(set(moved.n_budget_lines_changed))}"
          f" budget line(s) changed each — the harness TERM and the `systematic` TOTAL "
          f"that sums it; {int(moved.n_budget_lines_changed_outside_harness.sum())} "
          f"changed lines are anything else")
    assert (still.n_budget_lines_changed == 0).all(), \
        "a cell whose handler did not change produced a different budget"
    assert (moved.n_budget_lines_changed_outside_harness == 0).all(), \
        "a misattributed cell changed a budget line that is not the harness term or its total"
    assert (ok.stat_dex == ok.stat_dex_after).all(), "the harness term moved a STAT bar"
    assert (ok.A == ok.A_published.fillna(ok.A)).all(), "an ABUNDANCE moved"

    # ── the finding ───────────────────────────────────────────────────────────────
    print(f"\n=== {len(moved)} cells were charged another handler's systematic ===")
    if len(moved):
        print(f"{'band':<13}{'ion':<4}{'treatment':<16}{'deck':<14}{'n':>5}"
              f"{'charged':>10}{'label':>20}{'syst':>9}{'->':>4}{'syst':>9}")
        for _, r in moved.sort_values(["band", "ion", "deck"]).iterrows():
            print(f"{r.band:<13}{r.ion:<4}{r.treatment:<16}{(r.deck or 'root'):<14}"
                  f"{r.n_lines:>5}{r.harness_published_dex:>10.4f}"
                  f"{r.handler_published:>20}{r.syst_before:>9.4f}{'->':>4}"
                  f"{r.syst_after:>9.4f}")
        print(f"\n  every one is the SAME handler under two names: "
              f"{sorted(set(moved.handler_correct))} produced them, "
              f"{sorted(set(moved.handler_published))} was charged and printed.")
        print(f"  the bar was too LARGE in every case "
              f"(max |d_syst| {moved.d_syst.abs().max():.4f} dex), which is why nothing "
              f"was ever going to make it look wrong.")

    # ── item 4: the four cells' budgets, regenerated, diffed ──────────────────────
    # ONE file, not one per cell: the deriver itself writes all of a stem's treatments
    # into a single `*_budgets.txt`, and six near-identical artifacts would each need
    # their own GENERATORS.yaml entry saying the same sentence (RYA-686).
    diffs, after_texts = [], []
    for _, r in moved.sort_values(["band", "ion", "deck"]).iterrows():
        tb, ta = texts[r.src]
        after_texts.append(f"# {r.src}  ({r.deck or 'root'})\n{ta}")
        diffs.append("\n".join(difflib.unified_diff(
            tb.splitlines(), ta.splitlines(),
            fromfile=f"{r.src}  (as published)", tofile=f"{r.src}  (RYA-869)", lineterm="")))
    (out / "rya869_budgets_after.txt").write_text("\n\n".join(after_texts) + "\n")
    (out / "rya869_budget_diff.txt").write_text("\n\n".join(diffs) + "\n")

    # ── the published matrix cells this moves ─────────────────────────────────────
    matrix = ROOT / "data" / "results" / "rya783" / "fe_product_matrix.csv"
    mrows = []
    if matrix.exists():
        m = pd.read_csv(matrix)
        for _, r in m.iterrows():
            hit = moved[(moved.band == r.band) & (moved.ion == r.ion)
                        & (moved.treatment == r.treatment)
                        & (moved.src.apply(lambda s, x=str(r["_src"]): Path(x).name
                                           .replace("_products.csv", "") in s))]
            if len(hit):
                mrows.append({"band": r.band, "element": r.element, "ion": r.ion,
                              "treatment": r.treatment,
                              "syst_published": float(r.syst_dex),
                              "syst_rya869": float(hit.iloc[0].syst_after),
                              "A_unchanged": float(r.A)})
    print(f"\n=== published RYA-783/807/832 matrix cells that move ===")
    if not mrows:
        print("  none matched — the matrix at data/results/rya783 does not name these cells")
    for r in mrows:
        print(f"  {r['band']:<13}{r['element']} {r['ion']:<4}{r['treatment']:<16}"
              f"syst {r['syst_published']:.4f} -> {r['syst_rya869']:.4f}   "
              f"A {r['A_unchanged']:.3f} UNCHANGED")

    scd = harness_residual_divergences()
    print(f"\n=== do any charged harness residuals disagree with their control? ===")
    if scd["n_declared"]:
        for h, r in scd["declared"].items():
            print(f"  ⚠️ {h}: charged {r['charged_dex']:.4f}, control "
                  f"{r['control_dex_offset']} — {r.get('ticket', '')}")
    else:
        for h, r in scd["agreeing"].items():
            print(f"  {h}: charged {r['charged_dex']:.4f}, control "
                  f"{r['control_dex_offset']} over {r['control_n_lines']} lines — "
                  f"{'AGREE' if r.get('agrees') else '⚠️ DISAGREE'}")
    print(f"  {scd['status']}")

    d.to_csv(out / "rya869_harness_by_cell.csv", index=False)
    summary = {
        "ticket": "RYA-869",
        "band_products": str(bp),
        "band_products_fingerprint_sha256": tree_fingerprint(bp),
        "n_cells": int(len(d)),
        "n_reproducing_published": int(len(have) - len(bad)),
        "n_not_reproducing": int(len(bad)),
        "n_misattributed": int(len(moved)),
        "misattributed_cells": json.loads(
            moved[["band", "element", "ion", "treatment", "deck", "n_lines",
                   "handler_published", "handler_correct", "harness_published_dex",
                   "harness_correct_dex", "syst_before", "syst_after", "d_syst"]]
            .to_json(orient="records")),
        "published_matrix_cells_moved": mrows,
        "controls": {
            "pre_fix_rule_reproduces_published": int(len(have) - len(bad)),
            "unmoved_cells_byte_identical": int((still.n_budget_lines_changed == 0).sum()),
            "moved_cells_changing_only_the_harness_term_and_its_total":
                int((moved.n_budget_lines_changed_outside_harness == 0).sum()),
            "abundances_moved": 0,
            "stat_bars_moved": 0,
        },
        "harness_residual_divergences": scd,
    }
    (out / "rya869_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"\n  wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
