"""RYA-808 — registry hygiene, ratified.

RYA-807 had to pick a discriminator BECAUSE the registry held rows asserting two
incompatible things at once: `required_treatment=exclude` (the prescribed treatment is
removal) together with `status=owed` (we have not established why). RYA-808 removes the
contradiction from the data rather than leaving the code to arbitrate it, and pins the
semantic so it cannot drift back.
"""
from __future__ import annotations

import csv

import pytest

from pipeline import problem_children as pc


@pytest.fixture(scope="module")
def registry_rows():
    with pc.REGISTRY_CSV.open() as fh:
        return list(csv.DictReader(fh))


# ── §5 self-consistency ──────────────────────────────────────────────────────
def test_no_row_prescribes_exclusion_for_an_undiagnosed_cause(registry_rows):
    """THE INVARIANT. `exclude` says what to do; `owed` says we do not know why. A
    treatment cannot be prescribed from an undiagnosed cause (RYA-161), so the pair is a
    contradiction and `investigate` is the honest value."""
    bad = [r["lambda_or_scope"] for r in registry_rows
           if r["required_treatment"].strip() == "exclude"
           and r["status"].strip() == "owed"]
    assert not bad, f"exclude+owed is contradictory; use investigate: {bad}"


def test_investigate_is_the_treatment_for_owed_rows(registry_rows):
    """Ryan ratified 2026-08-13: `investigate`, not blank. Blank would read as 'no
    treatment needed', which is the opposite of what an owed row means."""
    owed = [r for r in registry_rows if r["status"].strip() == "owed"]
    assert owed, "the registry should still carry owed rows"
    for r in owed:
        t = r["required_treatment"].strip()
        assert t and t != "exclude", (r["lambda_or_scope"], t)


# ── §3 the semantic, pinned ──────────────────────────────────────────────────
def test_investigate_never_excludes():
    """The whole point of the relabel: it must change what the registry SAYS without
    changing what the harness DOES."""
    assert pc.aggregate_action(
        dict(required_treatment="investigate", status="owed")) == "flag"
    assert pc.aggregate_action(
        dict(required_treatment="investigate", status="active")) == "flag"


def test_only_exclude_plus_active_removes_a_line():
    """Exhaustive over the treatments the registry actually uses."""
    treatments = {"exclude", "investigate", "synthesis", "HFS_sum", "upper_limit",
                  "NLTE_grid", "none", "cited_substitution", "per_region_source",
                  "astrophysical_gf_differential"}
    for t in treatments:
        for status in ("active", "owed", "resolved", "predicted"):
            act = pc.aggregate_action(dict(required_treatment=t, status=status))
            expected = "exclude" if (t == "exclude" and status == "active") else (
                "none" if status in ("resolved", "predicted") else "flag")
            assert act == expected, (t, status, act, expected)


# ── §2 the four that nothing was tracking ────────────────────────────────────
@pytest.mark.parametrize("lam", [7941.832, 7366.370, 8090.325, 8592.951])
def test_the_untracked_outliers_are_registered_and_kept(lam):
    """RYA-807 surfaced these above A = 8.0 with no registry row anywhere. Registered so
    they are TRACKED — and flagged, not excluded, because the cause is unknown."""
    d = pc.disposition_for_line("Fe", "I", lam)
    assert d is not None, f"{lam} must be tracked"
    assert d["problem_class"] == "ABUNDANCE_OUTLIER"
    assert d["status"] == "owed"
    assert d["required_treatment"] == "investigate"
    assert pc.aggregate_action(d) == "flag", "an undiagnosed outlier must NOT be excluded"


# ── §4 the ceiling clause, declined on purpose ───────────────────────────────
def test_the_aggregate_is_never_gated_on_the_abundance_VALUE():
    """RYA-807's acceptance asked that 'no line above the pool ceiling remains in any
    aggregate'. RYA-808 drops that: a physical-plausibility ceiling is a DIFFERENT
    mechanism from a curated registry, and implementing it here would exclude lines whose
    cause is not established — the exact thing `investigate` exists to prevent.

    There is no such gate in the code today (verified across pipeline/ and scripts/), and
    this test exists so one cannot be added inside the registry path without a ticket that
    argues for it. If a plausibility gate is ever wanted, it is its own mechanism with its
    own name.
    """
    import inspect
    src = inspect.getsource(pc.aggregate_action)
    for token in ("abundance", "8.0", "ceiling", "plausib"):
        assert token not in src, (
            f"aggregate_action must decide on PROVENANCE, not on the measured value "
            f"(found {token!r})")
