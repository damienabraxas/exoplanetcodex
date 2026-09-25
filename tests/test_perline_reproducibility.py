"""RYA-870 — the reproduce-or-fail guard, as a test (RYA-489 §6.1).

⚠️ WHY THIS VALIDATES AN ARTIFACT INSTEAD OF RE-RUNNING THE SYNTHESIS. The reproduction
needs iSpec, which lives in venv312; pytest lives in venv_ci, and the two are not
interchangeable on Sirius. A test that simply skipped when iSpec was absent would be green
in CI forever while proving nothing — the guard would have no teeth exactly where teeth are
required. So the expensive half runs as `scripts/rya870_reproduce_perline.py` under
venv312 and commits its result; this test refuses to pass unless that result exists, is
CURRENT for the product it describes, and reports no failures.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PRODUCT = ROOT / "data" / "products" / "solar" / "Fe_perline.csv"
REPORT = ROOT / "data" / "results" / "rya870" / "rya870_reproducibility.json"
QUANTISER_FLOOR_DEX = 0.01     # RYA-771


def _header(path: Path) -> dict:
    out = {}
    for line in path.read_text().splitlines():
        if not line.startswith("#"):
            break
        if ": " in line:
            k, v = line[1:].split(": ", 1)
            out[k.strip()] = v.strip()
    return out


@pytest.fixture(scope="module")
def report():
    if not REPORT.exists():
        pytest.fail(
            f"{REPORT.relative_to(ROOT)} is missing. The per-line product is not "
            f"replication-grade until the guard has RUN — on Sirius, venv312:\n"
            f"    python3 scripts/rya870_reproduce_perline.py --star solar --element Fe")
    return json.loads(REPORT.read_text())


def test_the_guard_describes_the_product_that_is_committed(report):
    """🔴 A stale pass is not a pass. If the product was regenerated after the guard ran,
    the result on disk describes a file that no longer exists, and treating it as evidence
    would be the RYA-848 banked-artifact confound in miniature."""
    assert PRODUCT.exists(), f"{PRODUCT.relative_to(ROOT)} missing"
    assert report["product_commit_sha"] == _header(PRODUCT)["commit_sha"], (
        "the reproducibility report was produced against a different revision of the "
        "product than the one committed — re-run the guard")


def _ew_route_rows_in_the_product() -> tuple[int, int]:
    """(testable, labelled-but-widthless) in-aggregate EW-route rows, off the product.

    ⚠️ TESTABLE IS NARROWER THAN LABELLED, and conflating them is a mistake I made
    writing this test: it asserted the guard should have tested 1270 rows that carry no
    `ew_mA` at all. The guard re-inverts an equivalent width, so a row with no width is
    not something it is declining to test -- there is nothing there to invert. The second
    number is reported separately because it is a real defect, just not this guard's:
    the ENGINE-A band products set `ew_inversion=True` while their own `ew_method` says
    the measurement was a synthesis flux-fit.
    """
    import csv
    import io
    body = "".join(l for l in PRODUCT.read_text().splitlines(keepends=True)
                   if not l.startswith("#"))
    testable = widthless = 0
    for r in csv.DictReader(io.StringIO(body)):
        if r["status"] != "in_aggregate" or r["method"] != "ew_integration":
            continue
        if (r["ew_mA"] or "").strip() in ("", "nan"):
            widthless += 1
        else:
            testable += 1
    return testable, widthless


def test_a_zero_coverage_verdict_is_true_of_the_product(report):
    """🔴 "NOTHING TO TEST" IS A CLAIM, AND IT GETS CHECKED (RYA-1229).

    The reproduction re-inverts an equivalent width, so it can only test rows whose
    abundance came from EW inversion. Since the product became a projection of the feed,
    no published Fe product with force-added per-line evidence uses the PROFILEFIT (EW)
    route — all ten are among the 21 unresolved — so the guard has zero coverage. That is
    a real loss and it is recorded as NO-COVERAGE rather than hidden: the guard's previous
    8/8 PASS was measured on rows from two SUPERSEDED artifacts (`1D-LTE (ts-lte)`) that
    back no published product.

    ⚠️ This test is what stops NO-COVERAGE becoming a silent skip. If the product does
    contain EW-route rows, a NO-COVERAGE verdict is a BUG in the guard, not a state of the
    repo — and the next test then demands they all reproduce.
    """
    testable, widthless = _ew_route_rows_in_the_product()
    if report.get("verdict") == "NO-COVERAGE":
        assert testable == 0, (
            f"the guard reported NO-COVERAGE but the product has {testable} in-aggregate "
            f"EW-route rows WITH an ew_mA that it could have tested — the guard is not "
            f"looking at them")
        assert report["n_tested"] == 0 and not report["results"]
        assert "NO EW-ROUTE ROWS TO TEST" in report["coverage_note"], (
            "a zero-coverage report must say WHY, in the report, not only in a log")
        # 🔴 AND THE OTHER HALF OF THE REASON MUST BE COUNTED, NOT PARAPHRASED. A report
        # that said only "no EW rows" would hide that 1270 rows claim the EW route and
        # carry no width — the source-side flag defect the note points at.
        assert report.get("n_ew_labelled_without_ew_mA") == widthless, (
            f"the report says {report.get('n_ew_labelled_without_ew_mA')} EW-labelled "
            f"rows have no ew_mA; the product has {widthless}")
    else:
        assert testable > 0, (
            f"the guard claims to have tested {report['n_tested']} rows but the product "
            f"has no EW-route row with a width — the report describes another product")


def test_every_sampled_row_reproduces_its_own_number(report):
    """The deliverable's teeth: a row that cannot reproduce its A(X) from its own published
    constants is a FAILURE of the product, never a warning."""
    if report.get("verdict") == "NO-COVERAGE":
        # asserted for real in test_a_zero_coverage_verdict_is_true_of_the_product
        pytest.skip("the guard has no EW-route rows to test; the verdict is checked "
                    "against the product in test_a_zero_coverage_verdict_is_true_of_"
                    "the_product")
    assert report["n_tested"] > 0, "a guard that tested nothing has not run"
    failures = [r for r in report["results"] if r.get("outcome") == "FAIL"]
    assert not failures, (
        "rows failed to reproduce their own published abundance within "
        f"{QUANTISER_FLOOR_DEX} dex: "
        + "; ".join(f"{r['wavelength_air_A']} delta={r.get('delta_dex')}"
                    for r in failures))


def test_uncovered_rows_are_reported_not_counted_as_passing(report):
    """An untested row must never look like a tested one. synthesis_fit rows cannot be
    reproduced from the row alone (they need the observed spectrum), so they are counted
    separately and named in the report."""
    assert "n_synthesis_rows_NOT_COVERED" in report
    assert report["n_tested"] == report["n_passed"] + report["n_failed"]
    assert report["n_tested"] <= report["n_ew_route_in_aggregate"]
