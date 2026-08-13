"""
tests/test_grade_subject_naming_rya711.py
=========================================
RYA-711 items 1 and 2 — a grade must name its subject, and the gf cut must be derived.

WHAT THESE PIN, AND WHY IT IS WORTH PINNING
-------------------------------------------
The failure this prevents is not a crash. It is a READER — human or downstream table —
seeing a `B` and not knowing whether it means "log gf good to <=10 %" (NIST, about the
atomic data) or "composite line_score 0.60-0.80" (ours, about our measurement). Those are
unrelated claims wearing one glyph, and the repo already contains a table that puts both
in adjacent columns.

So:
  * our grade VALUES must carry the MQ- prefix, not merely our column name -- a value
    outlives the header it was emitted under;
  * NIST's letters must stay unprefixed, because they are someone else's published scale;
  * the two vocabularies must not overlap at all;
  * and the >25 % gf cut must be DERIVABLE from the numbers it claims to follow, not
    asserted -- the whole point of item 2 is that `# >25 % -- cull` read as arbitrary.
"""
import math
import re
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from config.constants import LINE_GRADE_THRESHOLDS  # noqa: E402
from pipeline.curate_nonfe_pools import NIST_GRADE_HIGH, NIST_GRADE_CULL  # noqa: E402

RYA561_GATE_DEX = 0.10          # the ratification tolerance the cut corresponds to


# ── item 1: the grade names its subject ──────────────────────────────────────

def test_our_grade_values_are_mq_prefixed():
    """The prefix lives in the VALUE. A bare 'B' copied into another table is ambiguous
    forever; 'MQ-B' explains itself wherever it lands."""
    assert set(LINE_GRADE_THRESHOLDS) == {"MQ-A", "MQ-B", "MQ-C", "MQ-D"}


def test_our_grades_never_collide_with_a_nist_grade():
    """The two vocabularies must be disjoint -- that is the whole of item 1."""
    nist = NIST_GRADE_HIGH | NIST_GRADE_CULL | {"C+", "C"}
    assert not (set(LINE_GRADE_THRESHOLDS) & nist)


def test_nist_letters_are_left_unprefixed():
    """We must NOT 'fix' the collision by renaming someone else's published scale."""
    for g in NIST_GRADE_HIGH | NIST_GRADE_CULL:
        assert not g.startswith("MQ"), f"{g} is a NIST grade and must not be re-labelled"


def test_thresholds_are_ordered_and_span_the_score_range():
    t = LINE_GRADE_THRESHOLDS
    assert t["MQ-A"] > t["MQ-B"] > t["MQ-C"] > t["MQ-D"] == 0.0
    assert t["MQ-A"] <= 1.0


def test_the_grader_emits_prefixed_values():
    from pipeline.abundances_derive import _compute_line_scores
    n = 6
    df = pd.DataFrame({
        "element": ["Fe"] * n, "ion": ["I"] * n,
        "wavelength_air_A": [5000.0 + i for i in range(n)],
        "ew_mA": [50.0] * n, "ew_err_mA": [1.0] * n,
        "a_1dlte": [7.5, 7.5, 7.5, 7.5, 7.5, 7.5],
        "fit_chi2_red": [1.0] * n, "vald_proximity_flag": [1.0] * n,
        "delta_nlte": [0.0] * n,
    })
    ew = pd.DataFrame({"element": ["Fe"] * n, "ion": ["I"] * n,
                       "wavelength_air_A": df["wavelength_air_A"],
                       "ew_mA": df["ew_mA"], "ew_err_mA": df["ew_err_mA"],
                       "fit_chi2_red": df["fit_chi2_red"]})
    out = _compute_line_scores(df, ew)
    assert "mq_grade" in out.columns, "the grade column must name its subject"
    assert "line_grade" not in out.columns, "the ambiguous name must be gone"
    assert out["mq_grade"].str.match(r"^MQ-[ABCD]$").all()


def test_no_bare_line_grade_survives_in_code():
    """A rename that leaves half the codebase on the old name is worse than none: two
    names for one thing is how the ambiguity comes back."""
    hits = []
    for d in ("pipeline", "scripts"):
        for p in (ROOT / d).rglob("*.py"):
            txt = p.read_text(encoding="utf-8", errors="ignore")
            for m in re.finditer(r"""['"]line_grade['"]|['"]element_grade['"]""", txt):
                hits.append(f"{p.relative_to(ROOT)}:{txt[:m.start()].count(chr(10)) + 1}")
    assert not hits, "bare grade names still present: " + ", ".join(hits)


# ── item 2: the gf cut is derived, not asserted ──────────────────────────────

def _dex(pct: float) -> float:
    return math.log10(1.0 + pct / 100.0)


@pytest.mark.parametrize("pct,expected", [
    (0.3, 0.0013), (1, 0.0043), (2, 0.0086), (3, 0.0128), (7, 0.0294), (10, 0.0414),
    (18, 0.0719), (25, 0.0969), (40, 0.1461), (50, 0.1761), (75, 0.2430),
])
def test_the_documented_percent_to_dex_table_re_derives(pct, expected):
    """Every number written into SCIENCE_STANDARDS and the code comment, recomputed.
    A documented derivation that does not reproduce is worse than an undocumented one."""
    assert round(_dex(pct), 4) == expected


def test_the_cut_is_the_last_rung_inside_the_gate():
    """The actual claim: C (25 %) still fits inside +-0.10 dex and D+ (40 %) does not.
    That, not roundness, is why the cut sits where it does."""
    assert _dex(25) < RYA561_GATE_DEX < _dex(40)
    assert round(_dex(25) / RYA561_GATE_DEX, 2) == 0.97


def test_every_culled_grade_exceeds_the_gate_and_every_high_grade_does_not():
    pct = {"AAA": 0.3, "AA": 1, "A+": 2, "A": 3, "B+": 7, "B": 10,
           "C+": 18, "C": 25, "D+": 40, "D": 50, "E": 75}
    for g in NIST_GRADE_HIGH:
        assert _dex(pct[g]) < RYA561_GATE_DEX, f"{g} is HIGH but exceeds the gate"
    for g in NIST_GRADE_CULL - {"F"}:            # F is a local legacy code, not a NIST tier
        assert _dex(pct[g]) > RYA561_GATE_DEX, f"{g} is culled but fits inside the gate"


def test_c_and_c_plus_fall_through_deliberately():
    """The documented middle tier: at the gate, not past it, so reference-tiered rather
    than excluded. Pinned because it is behaviour defined by ABSENCE from both sets --
    exactly the kind that gets 'tidied up' by someone adding C to the cull list."""
    assert "C" not in NIST_GRADE_HIGH and "C" not in NIST_GRADE_CULL
    assert "C+" not in NIST_GRADE_HIGH and "C+" not in NIST_GRADE_CULL
    assert _dex(18) < RYA561_GATE_DEX and _dex(25) < RYA561_GATE_DEX


def test_the_derivation_is_written_where_a_reader_will_look():
    doc = (ROOT / "docs" / "SCIENCE_STANDARDS.md").read_text()
    assert "log10(1 + ε)" in doc or "log10(1 + eps)" in doc
    assert "0.0969" in doc and "RYA-561" in doc
    assert "MQ-" in doc and "name its subject" in doc.lower()

    src = (ROOT / "pipeline" / "curate_nonfe_pools.py").read_text()
    assert "log10(1 + eps)" in src, "the code site must carry the derivation too"
    assert "0.0969" in src
    # the thing item 2 exists to remove
    assert "# >25 % — cull\n" not in src
