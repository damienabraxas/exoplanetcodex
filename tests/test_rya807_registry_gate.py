"""RYA-807 — the curated registry gates the band-product aggregate.

RYA-463 built the registry, RYA-774 landed the harness, and nothing connected them for
months. These pin the connection AND the rule that decides exclude-vs-flag, because that
rule is the part a later edit is most likely to get wrong in a way nothing else catches.
"""
from __future__ import annotations

import pytest

from pipeline import problem_children as pc
from pipeline.band_products import LineMeasurement


@pytest.fixture(scope="module")
def table():
    return pc.line_dispositions()


def test_registry_explodes_multi_wavelength_rows(table):
    """One registry row can name several lines: '7052.715 / 7120.021 (IR)'. Both must
    become matchable, or half a quarantine silently does nothing."""
    w = set(round(x, 3) for x in table.wavelength_air_A)
    assert {7052.715, 7120.021} <= w
    # ...and a non-numeric scope ("Procyon gf outliers") must NOT produce a line row:
    # every exploded row carries a real wavelength, never a parsed fragment of prose.
    assert table.wavelength_air_A.between(1000.0, 30000.0).all()


def test_exclude_requires_active_status():
    """THE RULE. `required_treatment` says what to do once the cause is known; `status`
    says whether it is. The registry holds rows that are BOTH exclude AND owed, and
    removing those would be tuning on an undiagnosed cause (RYA-161)."""
    diagnosed = dict(required_treatment="exclude", status="active")
    undiagnosed = dict(required_treatment="exclude", status="owed")
    assert pc.aggregate_action(diagnosed) == "exclude"
    assert pc.aggregate_action(undiagnosed) == "flag", "owed must never be excluded"


def test_non_exclude_treatments_flag_but_never_exclude():
    """'synthesis' / 'HFS_sum' / 'upper_limit' say MEASURE IT DIFFERENTLY, not drop it."""
    for treat in ("synthesis", "HFS_sum", "upper_limit", "NLTE_grid", "none"):
        assert pc.aggregate_action(dict(required_treatment=treat,
                                        status="active")) == "flag"


def test_resolved_and_predicted_are_no_action():
    assert pc.aggregate_action(dict(required_treatment="exclude",
                                    status="resolved")) == "none"
    assert pc.aggregate_action(dict(required_treatment="exclude",
                                    status="predicted")) == "none"
    assert pc.aggregate_action(None) == "none"


def test_the_line_that_started_this_is_excluded(table):
    """Fe I 8024.543 — ATOMIC_BLEND / exclude / critical / active, 'MISIDENTIFIED, the
    absorber is another species' — entered the merged RYA-783 product at A = 9.678."""
    d = pc.disposition_for_line("Fe", "I", 8024.543, table=table)
    assert d is not None, "8024.543 must be findable by the band path"
    assert pc.aggregate_action(d) == "exclude"
    assert "782" in str(d["governing_tickets"])


def test_rya764_outliers_are_kept_not_excluded(table):
    """The eight ABUNDANCE_OUTLIER rows are `exclude`+`owed`. They stay IN the aggregate
    with a flag: their cause is not established, and RYA-807 §2 is explicit that removing
    them would be tuning."""
    for w in (8615.311, 6604.585, 5783.907, 4769.812):
        d = pc.disposition_for_line("Fe", "I", w, table=table)
        assert d is not None and d["problem_class"] == "ABUNDANCE_OUTLIER"
        assert pc.aggregate_action(d) == "flag"


def test_rya780_scale_evidence_is_retained(table):
    """7190.122 / 7330.137 are RYA-780's own signal. Recorded so the harness can SEE the
    quarantine, with required_treatment='none' so they are NOT removed — deleting them
    would delete the evidence for the gf-scale argument."""
    for w in (7190.122, 7330.137):
        d = pc.disposition_for_line("Fe", "I", w, table=table)
        assert d is not None, f"{w} must be machine-readable, not prose-only"
        assert pc.aggregate_action(d) == "flag"


def test_ion_is_respected(table):
    """Fe I and Fe II at neighbouring wavelengths are different transitions."""
    assert pc.disposition_for_line("Fe", "II", 8024.543, table=table) is None


def test_tolerance_is_tight(table):
    """Both sides are CATALOGUED wavelengths — the largest real offset measured is
    0.0004 A — so a loose tolerance would start capturing neighbours."""
    assert pc.LINE_MATCH_TOL_A <= 0.05
    assert pc.disposition_for_line("Fe", "I", 8024.543 + 0.5, table=table) is None


def test_measurement_carries_the_verdict():
    """RYA-711 quarantine-not-cull / RYA-429 reported-never-dropped: an excluded line must
    still appear in the per-line output, saying why."""
    lm = LineMeasurement(element="Fe", ion="I", wavelength_air_A=8024.543,
                         instrument="kpno_solar_atlas", ew_mA=49.3, ew_method="profile")
    for f in ("problem_class", "problem_status", "problem_tickets", "problem_action"):
        assert hasattr(lm, f), f"{f} must reach the per-line CSV via asdict()"
        assert getattr(lm, f) == "", "defaults must leave existing artifacts unchanged"


def test_driver_consults_the_registry_for_any_element():
    """Harness-level, not Fe-special: the read is unconditional on element."""
    src = (pc.ROOT / "scripts" / "derive_band_products.py").read_text()
    assert "problem_children" in src
    assert "aggregate_action" in src or "_pc.aggregate_action" in src
    # The lookup MUST pass a.element -- it has to match the right species. What would make
    # this Fe-special is a CONDITIONAL: gating the read itself on which element is running.
    assert "_pc.line_dispositions()" in src
    # Scope the check to the registry block itself: an unrelated element conditional
    # elsewhere in the driver (the Engine-B deck lookup) is not this ticket's concern.
    import re as _re
    head = src.split("_pc_table = _pc.line_dispositions()")[0]
    block = head[-800:]
    assert not _re.search(r"if\s+a\.element\s*==", block), \
        "the registry read must not be gated on which element is running"
