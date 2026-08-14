"""RYA-809 — the per-line RCA verdicts, pinned.

RYA-808 gave twelve lines a LABEL; this gave them a CAUSE, or documented why none could be
established. These pin both halves, because the half that is easy to erode later is the
second one: an undiagnosed line quietly becoming an excluded one is precisely the RYA-161
failure this ticket was written to avoid.
"""
from __future__ import annotations

import csv

import pytest

from pipeline import problem_children as pc

#: diagnosed, cause established -> removed from the aggregate
CONVICTED = {4769.812, 4975.412, 5783.907, 6604.585, 7941.832}
#: every test run, no cause established -> stays in the aggregate, flagged
AMBIGUOUS = {5609.961, 6185.693, 8090.325, 8592.951, 4932.084, 7366.370, 8615.311}
#: clean on position, isolation AND depth — the line is exonerated, the inversion is not
EXONERATED = {8615.311, 6185.693}


@pytest.fixture(scope="module")
def rows():
    with pc.REGISTRY_CSV.open() as fh:
        return {r["lambda_or_scope"].split()[0]: r for r in csv.DictReader(fh)}


def test_all_twelve_carry_a_verdict(rows):
    for lam in CONVICTED | AMBIGUOUS:
        r = rows.get(f"{lam:.3f}")
        assert r is not None, f"{lam} must be in the registry"
        assert "RYA-809 RCA:" in r["notes"], f"{lam} must carry its RCA evidence"


@pytest.mark.parametrize("lam", sorted(CONVICTED))
def test_convicted_lines_are_excluded_with_a_named_cause(lam, rows):
    r = rows[f"{lam:.3f}"]
    assert r["required_treatment"] == "exclude"
    assert r["status"] == "active", "a cause was established, so it is no longer owed"
    assert r["problem_class"], "an exclusion must name its class"
    assert pc.aggregate_action(r) == "exclude"


@pytest.mark.parametrize("lam", sorted(AMBIGUOUS))
def test_ambiguous_lines_are_kept_not_excluded(lam, rows):
    """RYA-161: no exclusion without an established cause. These ran every test and none
    returned a positive signal, so they stay in the aggregate."""
    r = rows[f"{lam:.3f}"]
    assert r["required_treatment"] == "investigate"
    assert r["status"] == "owed"
    assert pc.aggregate_action(r) == "flag"


@pytest.mark.parametrize("lam", sorted(EXONERATED))
def test_the_exonerated_lines_implicate_the_inversion_not_the_line(lam, rows):
    """These match position to <1 mA, are isolated, and absorb LESS than their catalogued
    central_depth supports — yet the inversion returns A ~ 9.1. Excluding a clean line
    would be backwards; the note must say the method is implicated."""
    n = rows[f"{lam:.3f}"]["notes"]
    assert "EXONERATED" in n and "INVERSION" in n
    assert rows[f"{lam:.3f}"]["required_treatment"] == "investigate"


def test_registry_stays_self_consistent(rows):
    bad = [k for k, r in rows.items()
           if r["required_treatment"] == "exclude" and r["status"] == "owed"]
    assert not bad, f"exclude+owed reintroduced (RYA-808 invariant): {bad}"


def test_blend_verdicts_are_abundance_weighted():
    """The trap this RCA nearly fell into: `log_gf - EP*theta` is a PER-ATOM strength, so
    comparing Fe I (A=7.46) against Ce II (A~1.6) without the abundance term makes every
    rare-earth neighbour look dominant. A first pass scored neighbours at +3.5 to +6 dex
    and would have excluded ten lines on it."""
    import scripts.rya809_fe_rca as rca
    assert rca._abund("Fe I") == pytest.approx(7.46)
    assert rca._abund("C2") is None, "a molecule's density is not its constituent's abundance"
    assert rca._abund("Ce II") is None, "an untabulated abundance must be inconclusive"
    src = __import__("inspect").getsource(rca.verdict)
    assert "abundance weighting" in src or "a_nb" in src
