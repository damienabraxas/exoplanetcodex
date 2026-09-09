"""RYA-1206 — the CRIRES+ H ENGINE-B-NLTE product and why it is n=12, not n=25.

EMIT + DIAGNOSE. These pin the three things that would otherwise rot:

  * an NLTE product must contain ONLY genuinely-NLTE lines. Partial coverage is not a
    reason to refuse the engine (Ryan, RYA-1206, superseding the RYA-1203 refusal) -- the
    trap is fitting all 25 and letting 13 fall back to departure = 1 under an NLTE label;
  * the count is per LINE, not per neighbourhood. `pool_label_coverage` says 13 because it
    credits Fe I 16179.583 with a label 0.026 A away belonging to a line 2.5 dex weaker;
  * the shortfall is a DECK limit, not a fixable labelling defect -- so nobody re-opens it
    as a bug. Both levels exist; the transition does not.
"""
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "data/results/rya1206/h_nlte_coverage.json"


@pytest.fixture(scope="module")
def doc():
    if not DOC.exists():
        pytest.skip("RYA-1206 coverage artifact absent")
    return json.loads(DOC.read_text())


def test_the_h_pool_splits_twelve_genuinely_nlte_and_thirteen_lte_fallback(doc):
    lines = doc["lines"]
    assert len(lines) == 25, "the graded H Ruffoni pool is 25 lines"
    good = [l for l in lines if l["genuinely_NLTE"]]
    assert len(good) == 12, "12 genuinely-NLTE lines -- not the 13 a neighbourhood count gives"
    assert doc["headline"]["genuinely_NLTE_lines"] == 12
    assert doc["headline"]["LTE_fallback_lines"] == 13


def test_every_genuinely_nlte_line_is_the_e7D_multiplet(doc):
    """The 12 are not a scatter of lines that happened to match -- they are one multiplet,
    which is what makes 'the deck carries this transition family and not that one' a
    physical statement rather than a matching artefact."""
    good = [l for l in doc["lines"] if l["genuinely_NLTE"]]
    assert all(l["label_low"].startswith("e7D") for l in good), [l["label_low"] for l in good]
    assert all(l["label_up"].startswith("n7D") for l in good), [l["label_up"] for l in good]
    assert all(l["nlte_flag"] == "T" for l in good)


def test_the_dropped_lines_carry_no_label_at_all(doc):
    """They must be genuinely unlabelled in the list the synthesis reads -- not merely
    'unmatched by our search'."""
    bad = [l for l in doc["lines"] if not l["genuinely_NLTE"]]
    assert all(l["nlte_flag"] == "F" for l in bad)
    assert all(l["label_low"] == "none" and l["label_up"] == "none" for l in bad)


def test_the_split_is_by_excitation_and_the_bands_do_not_touch(doc):
    """The discriminator is EP, and the gap between the two populations is what makes the
    boundary a statement about the deck rather than a threshold someone chose."""
    good = [l["EP_eV"] for l in doc["lines"] if l["genuinely_NLTE"]]
    bad = [l["EP_eV"] for l in doc["lines"] if not l["genuinely_NLTE"]]
    assert max(good) < min(bad), "the EP bands must be disjoint"
    assert min(bad) - max(good) > 0.5, "and separated by a real gap, not a rounding edge"


def test_the_limit_is_not_a_wavelength_cut(doc):
    """🔴 THE POINT OF THE WHOLE DIAGNOSTIC. If the covered and dropped sets were separated
    in wavelength, this would be the same kind of limit ENGINE-A has and the two engines
    would be co-limited. They INTERLEAVE, so it is not."""
    lines = sorted(doc["lines"], key=lambda l: l["wave_A"])
    flags = [l["genuinely_NLTE"] for l in lines]
    flips = sum(1 for a, b in zip(flags, flags[1:]) if a != b)
    assert flips >= 4, (
        "covered and dropped lines must interleave in wavelength; a clean split would mean "
        f"a wavelength boundary after all (flips={flips})")


def test_no_dropped_line_is_a_labelling_defect(doc):
    """The classification Ryan asked for: reach vs labelling vs level-absent. If any line
    were IN the atom but unlabelled, that is a BUG to file -- and this says there are none."""
    bad = [l for l in doc["lines"] if not l["genuinely_NLTE"]]
    assert all("TRANSITION-ABSENT" in l["classification"] for l in bad)
    assert doc["headline"]["fixable_by_code"] is False


def test_the_engine_a_overlap_is_reported_as_degenerate_not_as_agreement(doc):
    """🔴 A 100% overlap with a UNIVERSAL set is not evidence of a shared cause. ENGINE-A
    serves nothing in H, so every subset overlaps it perfectly; reading that as 'both grids
    stop in the same place' would be the saturation trap. The artifact must say so."""
    x = doc["cross_check_vs_engine_a_rya1165"]
    assert "DEGENERATE" in x["verdict"]
    assert x["informative_comparison"]["lines Gerber covers that ENGINE-A cannot"] == 12
    assert x["informative_comparison"]["lines ENGINE-A covers that Gerber cannot"] == 0
    assert "COMPLEMENTARY" in x["answer"]


# ── the selection primitive itself ───────────────────────────────────────────────────

def _pair():
    """The real 16179 pair, in the schema `_select_species_rows` filters on.

    `turbospectrum_species` is REQUIRED and collapses the ionisation stages (every Fe row
    is 26.0), so `element` is what separates Fe I from Fe II -- a fixture missing either
    field silently selects nothing and the test passes for the wrong reason.
    """
    return np.array(
        [(16179.583, 26.0, "Fe 1", "F", "none", "none"),
         (16179.609, 26.0, "Fe 1", "T", "w5G4*", "f5G4")],
        dtype=[("wave_A", "f8"), ("turbospectrum_species", "f8"), ("element", "U8"),
               ("nlte", "U2"), ("nlte_label_low", "U8"), ("nlte_label_up", "U8")])


def test_the_fixture_actually_selects_both_rows():
    """Control: if the fixture schema drifts, `_select_species_rows` returns nothing and
    every assertion below passes vacuously."""
    from pipeline import gerber_nlte as G
    assert len(G._select_species_rows(_pair(), 26, 1)) == 2


def test_pool_label_mask_reads_the_line_not_its_neighbour():
    """The defect that made 12 look like 13, pinned on a synthetic pair: a target with no
    label and a labelled neighbour inside the tolerance. The mask must say False."""
    from pipeline import gerber_nlte as G
    ll = _pair()
    m = G.pool_label_mask(ll, 26, [16179.583], ion=1)
    assert m.tolist() == [False], (
        "the target is nlte=F; the label 0.026 A away belongs to another transition")
    m2 = G.pool_label_mask(ll, 26, [16179.609], ion=1)
    assert m2.tolist() == [True], "and the neighbour itself IS labelled"


def test_pool_label_coverage_still_over_counts_that_pair():
    """The control that keeps the finding honest: the defect is in the COVERAGE function,
    and this ticket did not change it (that is Ryan's to file). If this ever starts
    agreeing with the mask, the over-count was fixed and the note above can be retired."""
    from pipeline import gerber_nlte as G
    ll = _pair()
    lab, tot = G.pool_label_coverage(ll, 26, [16179.583], ion=1)
    assert (lab, tot) == (1, 1), "pool_label_coverage credits the neighbour's label"
