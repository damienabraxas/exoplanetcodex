"""RYA-1214 — `--lines-from-set`, and the two things that make it honest.

WHY IT EXISTS. `select_lines` ranks by theoretical central depth and refuses below 0.15.
Measured on the GES v6 list the VIS and red-optical bands synthesise, that floor excludes
BOTH of AGSS21's forbidden CNO indicators ([C I] 8727 at depth 0.030, [O I] 6300 at 0.040)
and ALL FIVE of its N I lines (max 0.040) — the production rule reaches 11 of the 23
atomic indicators, and `vis_N_*` / `vis_O_*` failed outright with "no line has theoretical
depth in [0.15, 0.9]". The light elements' indicators are weak BECAUSE they are
unsaturated, which is what makes them good diagnostics.

TWO INVARIANTS:

1. **A named set is not an EW comparison.** The selector is part of the product KEY
   (RYA-984), so an AGSS21-set run must not be tagged `_FROMEW`.
2. **A named set carries its own match window, or it is refused.** `_EW_MATCH_TOL_A` is
   0.005 A because it pairs two writers of the SAME 4-decimal list. AGSS21 prints lambda
   in nanometres to 2 dp, so at that window [C I] 8727.12 (list: 8727.139) is simply not
   found — the campaign's headline carbon line would vanish from its own product with no
   verdict to fall through to. Same shape as RYA-1211, one stage earlier.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
for p in (str(ROOT), str(ROOT / "scripts")):
    if p not in sys.path:
        sys.path.insert(0, p)

SETS = ROOT / "data" / "linelists" / "reference_sets"
PROV = SETS / "agss21_cno_lineset_rya1214.prov.json"


def _load_dbp():
    """Import `derive_band_products` WITHOUT iSpec, which it only needs at run time.

    The module imports `ispec` transitively through the pipeline, and this test is about
    two pure functions. Skipping where iSpec is absent would make the test invisible on
    every machine but Sirius — and these are exactly the checks a Mac edit can break.
    """
    spec = importlib.util.spec_from_file_location(
        "_dbp_rya1214", ROOT / "scripts" / "derive_band_products.py")
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    except Exception as exc:                      # pragma: no cover - env dependent
        pytest.skip(f"derive_band_products not importable here: {exc}")
    return mod


# ── 1. the selector says which set drove the run ─────────────────────────────
def test_a_named_set_is_not_tagged_as_an_ew_comparison():
    dbp = _load_dbp()
    a = argparse.Namespace(lines_from_set="AGSS21=/tmp/x.csv", lines_from_ew=None,
                           lines_tier="all", lines_deep_graded=False)
    assert dbp._selector_tag(a) == "_SET-AGSS21"
    # and the EW route is untouched
    b = argparse.Namespace(lines_from_set=None, lines_from_ew="/tmp/ew.csv",
                           lines_tier="all", lines_deep_graded=False)
    assert dbp._selector_tag(b) == "_FROMEW"


# ── 2. the window travels with the set ───────────────────────────────────────
def test_the_default_window_would_have_dropped_the_headline_carbon_line():
    """The failure this mechanism exists to prevent, asserted as arithmetic."""
    dbp = _load_dbp()
    census_A, list_A = 8727.12, 8727.139        # AGSS21's printing vs the GES v6 row
    sep = abs(list_A - census_A)
    assert sep > dbp._EW_MATCH_TOL_A, (
        "the default window would have found [C I] 8727 after all — re-derive this test "
        "from the current wavelengths rather than deleting it")
    w_sorted = np.array([list_A])
    keep, missing = dbp._match_into_list(np.array([census_A]), w_sorted, np.array([0]))
    assert keep == [] and len(missing) == 1, "the default window must MISS it"
    keep, missing = dbp._match_into_list(np.array([census_A]), w_sorted, np.array([0]),
                                         tol_A=0.05)
    assert list(keep) == [0] and not missing, "the set's own window must FIND it"


def test_widening_is_per_call_and_does_not_loosen_the_ew_comparison():
    dbp = _load_dbp()
    assert dbp._EW_MATCH_TOL_A == 0.005, (
        "the EW constant was measured for pairing two writers of one 4-decimal list "
        "(RYA-1036); a named set must pass its own window rather than widen this")


# ── 3. every shipped set declares a derived window ───────────────────────────
@pytest.mark.parametrize("path", sorted(SETS.glob("agss21_cno_*_rya1214.csv")))
def test_every_shipped_cno_set_declares_its_match_window_and_its_basis(path):
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert rows, f"{path.name} is empty"
    tol = {r["match_tol_A"] for r in rows}
    assert tol == {"0.05"}, (
        f"{path.name} declares {tol}; AGSS21 prints lambda in NANOMETRES to 2 dp = 0.1 A, "
        f"so the rounding half-width is 0.05 A — the same derivation "
        f"reference_lineset.SETS['asplund'] records for the Fe table of the same paper")
    for r in rows:
        assert "printed resolution" in r["match_tol_basis"], (
            "a window without its basis is a number somebody typed")


def test_the_adopted_and_grid_sets_are_never_merged():
    """The census says grid membership is NOT proof AGSS21 adopted a line. Merging the
    two would publish a set the paper does not claim."""
    import json
    prov = json.loads(PROV.read_text())
    statuses = {f["use_status"] for f in prov["files"]}
    assert len(prov["files"]) >= 2 and len(statuses) == 2, (
        f"expected the adopted and grid sets in separate files, got {prov['files']}")
    for f in prov["files"]:
        p = ROOT / f["file"]
        with p.open(newline="") as fh:
            got = {r["use_status"] for r in csv.DictReader(fh)}
        assert got == {f["use_status"]}, f"{p.name} mixes {got}"
