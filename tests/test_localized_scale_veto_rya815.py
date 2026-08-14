"""
tests/test_localized_scale_veto_rya815.py — RYA-815
===================================================
A reference self-contradiction must stop ONE element, not the channel.

RYA-681 was right to make the contradiction LOUD — it converted a silent
0.05 dex/cycle ratchet (RYA-669 measured A(Fe I) drifting to 7.416 with all nine
gates green) into a refusal. What was wrong was the BLAST RADIUS: one stale cell
vetoed all 28 species, and since the RYA-654 tracker generates from the artifact
phase_c can then no longer produce, every element's status froze behind it.

These tests pin both halves: the refusal still exists, and it is now local.
"""
from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import pytest

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.solar_scale_provenance import (  # noqa: E402
    SCALE_OFF_SCALE, SCALE_SELF_CONTRADICTION, SCALE_UNDECLARED,
    ScaleProvenanceError, classify_gold_scale, resolve_gold_scale,
    try_apply_reported_scale_correction)

CONSISTENT = {"method_scale": "3D-NLTE (Fe I)"}      # 3D value, 3D label
CONTRADICTORY = {"method_scale": "1D-NLTE (Fe I)"}   # 3D value, 1D label — gold v3


def _phase_c():
    sys.argv = ["phase_c"]
    spec = importlib.util.spec_from_file_location(
        "pc_rya815", ROOT / "scripts" / "phase_c_verdict_rya371.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# ── RYA-681's refusal is PRESERVED, not softened ─────────────────────────────
def test_resolve_gold_scale_still_raises_on_contradiction():
    with pytest.raises(ScaleProvenanceError, match="CONTRADICTS ITSELF"):
        resolve_gold_scale("Fe", CONTRADICTORY, 7.466)


def test_resolve_still_raises_on_off_scale_and_undeclared():
    with pytest.raises(ScaleProvenanceError, match="NEITHER recognised"):
        resolve_gold_scale("Fe", CONSISTENT, 7.416)      # the doubled-correction value
    with pytest.raises(ScaleProvenanceError, match="declares NO abundance scale"):
        resolve_gold_scale("Fe", {}, 7.466)


def test_the_two_surfaces_share_one_body_and_agree():
    """A guard that disagrees with itself is worse than no guard."""
    for row, val in ((CONTRADICTORY, 7.466), (CONSISTENT, 7.416), ({}, 7.466)):
        state, contra = classify_gold_scale("Fe", row, val)
        assert state is None and contra is not None
        with pytest.raises(ScaleProvenanceError) as exc:
            resolve_gold_scale("Fe", row, val)
        assert contra.message == str(exc.value)          # verbatim, not paraphrased


# ── the classifier names the offending cell ──────────────────────────────────
def test_contradiction_kinds_are_named_and_distinct():
    assert classify_gold_scale("Fe", CONTRADICTORY, 7.466)[1].kind == SCALE_SELF_CONTRADICTION
    assert classify_gold_scale("Fe", CONSISTENT, 7.416)[1].kind == SCALE_OFF_SCALE
    assert classify_gold_scale("Fe", {}, 7.466)[1].kind == SCALE_UNDECLARED


def test_contradiction_names_the_offending_cell():
    contra = classify_gold_scale("Fe", CONTRADICTORY, 7.466)[1]
    assert contra.element == "Fe"
    assert "Fe." in contra.cell            # e.g. Fe.method_scale
    assert contra.declared == "1D-NLTE" and contra.from_value == "3D-NLTE"


def test_consistent_row_classifies_cleanly():
    state, contra = classify_gold_scale("Fe", CONSISTENT, 7.466)
    assert state == "3D-NLTE" and contra is None


# ── the non-raising apply refuses to publish a number it cannot vouch for ────
def test_try_apply_returns_value_uncorrected_and_flags_the_contradiction():
    a, rec, contra = try_apply_reported_scale_correction("Fe", 7.466, CONTRADICTORY)
    assert contra is not None
    assert rec["applied"] is False
    assert a == 7.466                       # NOT silently corrected
    assert rec["contradiction_cell"] == contra.cell


def test_try_apply_matches_the_raising_form_when_consistent():
    from pipeline.solar_scale_provenance import apply_reported_scale_correction
    a1, r1 = apply_reported_scale_correction("Fe", 7.516, {"method_scale": "1D-NLTE (Fe I)"})
    a2, r2, contra = try_apply_reported_scale_correction(
        "Fe", 7.516, {"method_scale": "1D-NLTE (Fe I)"})
    assert contra is None
    assert (a1, r1["applied"], r1["correction_dex"]) == (a2, r2["applied"], r2["correction_dex"])


# ── THE ACCEPTANCE CASE: one element withheld, the rest proceed ──────────────
def test_injected_contradiction_withholds_one_element_and_the_rest_proceed():
    m = _phase_c()
    ab, ew, pa, _resolved = m._load("CURRENT")

    baseline = m.build_verdicts(ab, ew, pa)
    n_total = len(baseline)
    assert n_total >= 26
    assert not [r for r in baseline if r["verdict"] == "INDETERMINATE"], (
        "gold v4 is self-consistent (RYA-811), so the baseline must have none")

    # inject gold v3's exact shape: the 3D value under a 1D label
    poisoned = ab.copy()
    fe = poisoned["element"] == "Fe"
    assert fe.any()
    poisoned.loc[fe, "method_scale"] = "1D-NLTE (Fe I)"

    rows = m.build_verdicts(poisoned, ew, pa)
    assert len(rows) == n_total, "no element may be DROPPED from the counts"

    indet = [r for r in rows if r["verdict"] == "INDETERMINATE"]
    assert len(indet) == 1, f"expected exactly one withheld element, got {len(indet)}"
    fe_row = indet[0]
    assert fe_row["element"] == "Fe"

    # LOUD + NAMED
    assert "reference self-contradiction" in fe_row["owed"]
    assert "verdict withheld" in fe_row["owed"]
    assert "Fe." in fe_row["owed"]                       # the offending cell is named
    assert "CONTRADICTS ITSELF" in fe_row["owed"]        # the RYA-681 message survives
    assert fe_row["A_measured"] is None                  # no number we cannot vouch for

    # and every OTHER element is untouched, verdict for verdict
    others_before = {r["element"]: r["verdict"] for r in baseline if r["element"] != "Fe"}
    others_after = {r["element"]: r["verdict"] for r in rows if r["element"] != "Fe"}
    assert others_after == others_before, "one bad row must not disturb any other element"


def test_a_self_consistent_reference_produces_no_indeterminate():
    m = _phase_c()
    ab, ew, pa, _ = m._load("CURRENT")
    rows = m.build_verdicts(ab, ew, pa)
    assert not [r for r in rows if r["verdict"] == "INDETERMINATE"]
    fe = [r for r in rows if r["element"] == "Fe"][0]
    assert fe["A_measured"] == 7.466        # RYA-811's value, uncorrected twice


# ── INDETERMINATE is a DISTINCT, COUNTED, NAMED state downstream ─────────────
def test_indeterminate_is_declared_not_left_to_fall_through():
    """
    Falling through both `in` checks is only ACCIDENTALLY correct. Declaring the
    state makes its neutrality deliberate, and stops the next person extending
    MEASURED_VERDICTS from silently turning a withheld verdict into a data claim.
    """
    from pipeline.ledger_consistency_guard import (
        MEASURED_VERDICTS, NO_DATA_VERDICTS, WITHHELD_VERDICTS)
    assert WITHHELD_VERDICTS == {"INDETERMINATE"}
    assert not (WITHHELD_VERDICTS & MEASURED_VERDICTS)
    assert not (WITHHELD_VERDICTS & NO_DATA_VERDICTS)


def test_indeterminate_is_not_absorbed_into_pass_or_owed():
    """It must be its own bucket in the counts, never folded into a neighbour."""
    m = _phase_c()
    ab, ew, pa, _ = m._load("CURRENT")
    poisoned = ab.copy()
    poisoned.loc[poisoned["element"] == "Fe", "method_scale"] = "1D-NLTE (Fe I)"
    rows = m.build_verdicts(poisoned, ew, pa)

    import collections
    counts = collections.Counter(r["verdict"] for r in rows)
    base = collections.Counter(r["verdict"] for r in m.build_verdicts(ab, ew, pa))
    assert counts["INDETERMINATE"] == 1
    assert counts["PASS"] == base["PASS"] - 1          # Fe LEFT pass, not silently kept
    assert sum(counts.values()) == sum(base.values())  # nothing dropped from the totals
    for other in ("CURATION-OWED", "NLTE-OWED", "DATA-GAP"):
        assert counts[other] == base[other], f"{other} absorbed the withheld element"
