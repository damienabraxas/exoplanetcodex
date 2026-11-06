"""RYA-1233: the published grade rule (Fe's axis) applies to every element but Fe."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline.cno_grade import grade_for  # noqa: E402


def test_a_published_set_is_reference_grade_for_any_element():
    assert grade_for({"element": "Si", "selector": "SET-SI_AGSS21"}) == "Reference Grade"
    assert grade_for({"element": "Si", "selector": "SET-SI_ELGUETA2026"}) == "Reference Grade"
    assert grade_for({"element": "C", "selector": "SET-AGSS21"}) == "Reference Grade"
    assert grade_for({"element": "Al", "selector": None}) == "Codex Grade"


def test_fe_stamps_its_own():
    assert grade_for({"element": "Fe", "selector": "SET-AGSS21"}) is None


def test_the_deep_pool_is_deep_grade():
    assert grade_for({"element": "Si", "selector": "DEEPGRADED"}) == "Deep Grade"
    assert grade_for({"element": "Si", "selector": "GRADED"}) == "Codex Grade"


def test_the_laboratory_pool_is_reference_grade():
    # RYA-1311: Ca published 44 lab-gf products as Codex Grade
    assert grade_for({"element": "Ca", "selector": "REFERENCE", "tier": "REFERENCE"}) == "Reference Grade"


def test_publish_keeps_the_reference_tier():
    # RYA-1311: publish sent the Reference pool as tier ALL while its budget was scoped
    # REFERENCE, so every Reference product quarantined (scope mismatch)
    src = (Path(__file__).resolve().parents[1] / "pipeline" / "element_publish.py").read_text()
    assert "(\"GRADED\", \"DEEPGRADED\", \"REFERENCE\") else \"ALL\"" in src
