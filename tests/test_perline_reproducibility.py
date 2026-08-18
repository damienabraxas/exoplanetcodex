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


def test_every_sampled_row_reproduces_its_own_number(report):
    """The deliverable's teeth: a row that cannot reproduce its A(X) from its own published
    constants is a FAILURE of the product, never a warning."""
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
