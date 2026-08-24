"""An element run writes ONE page and cannot reach another's (RYA-1031).

🔴 WHY THIS EXISTS. On 2026-08-23 regenerating the single combined `live_status.json`
from an Fe-only checkout would have written a page with 3 elements and 10 sections over
one with 26 and 19 -- silently dropping 23 elements, the telluric block, the graded
block and the reporting contract. It did not happen because the output was diffed
against the served copy first and Ryan said stop, which is luck plus vigilance, not a
property of the system.

Splitting the page per element makes it a property. An Fe run derives its path from the
`--element` argument and nothing else, so there is no code path from "write Fe" to a
filename containing "al". These tests pin that, and pin the index staying DERIVED --
because the moment the index is authored by a run, the shared mutable file is back and
so is the race.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _tracker():
    spec = importlib.util.spec_from_file_location(
        "rya935_live_status", ROOT / "scripts" / "rya935_live_status.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _status(elements=("Fe", "Al")):
    """A minimal status doc shaped like the real one."""
    return {
        "generated": "2026-08-24T00:00:00+00:00",
        "generator": "scripts/rya935_live_status.py",
        "refresh_seconds": 5,
        "system": "solar",
        "elements": ["FeI", "AlI", "CI"],
        "bands": ["VIS"],
        "instruments": [{"holding": "solar_kpno_molecfit_corrected",
                         "instrument": "kpno_solar_atlas",
                         "telluric_applied": "applied"},
                        {"holding": "solar_kpno", "instrument": "kpno_solar_atlas",
                         "telluric_applied": "not-applied"}],
        "telluric": [], "telluric_summary": {}, "reference": {},
        "reporting_contract": {}, "model_matrix": {"cells": []},
        "products": [{"element": e, "band": "VIS", "instrument": "kpno_solar_atlas",
                      "holding": "solar_kpno_molecfit_corrected", "A": 7.4}
                     for e in elements],
        "graded": [{"element": e} for e in elements],
        "reachability": {}, "run_contexts": [], "variant_products": 0,
        "unattributed_products": 0, "pre_correction_products": 0,
    }


def test_an_element_run_writes_only_its_own_page(tmp_path):
    """🔴 THE INVARIANT: writing Fe must leave Al's bytes untouched."""
    t = _tracker()
    t.write_element_page(_status(), "Al", tmp_path)
    al = tmp_path / "elements" / "al.json"
    before = al.read_bytes()

    t.write_element_page(_status(), "Fe", tmp_path)

    assert al.read_bytes() == before, (
        "an Fe run modified al.json. The whole reason the page is split per element is "
        "that this must be impossible, not merely unlikely."
    )
    assert (tmp_path / "elements" / "fe.json").exists()


def test_an_element_page_carries_only_that_element(tmp_path):
    t = _tracker()
    dest = t.write_element_page(_status(), "Fe", tmp_path)
    page = json.loads(dest.read_text())
    assert page["element"] == "Fe"
    assert {r["element"] for r in page["products"]} == {"Fe"}, \
        "another element's products leaked onto the Fe page"
    assert {r["element"] for r in page["graded"]} == {"Fe"}


@pytest.mark.parametrize("bad", ["../al", "Fe/../Al", "fe.json", "", "FE", "Fe I", "*"])
def test_a_path_traversing_element_is_refused(bad):
    """The output path is built from this argument, so it must be an element or nothing."""
    t = _tracker()
    with pytest.raises(SystemExit):
        t._element_slug(bad)


def test_the_index_holds_no_abundance_values(tmp_path):
    """A headline number on an overview page is how a stale value gets quoted."""
    t = _tracker()
    t.write_element_page(_status(), "Fe", tmp_path)
    index = json.loads(t.write_index(_status(), tmp_path).read_text())
    assert "products" not in index and "graded" not in index, \
        "the index must carry status, not values -- the numbers live on element pages"
    assert any("A" in str(k) and isinstance(v, float)
               for k, v in index.items()) is False


def test_the_index_is_derived_from_whatever_pages_exist(tmp_path):
    """Derived, never authored: adding a page changes the index with no index write."""
    t = _tracker()
    t.write_element_page(_status(), "Fe", tmp_path)
    one = json.loads(t.write_index(_status(), tmp_path).read_text())
    assert [d["element"] for d in one["dashboard"]] == ["Fe"]

    t.write_element_page(_status(), "Al", tmp_path)
    two = json.loads(t.write_index(_status(), tmp_path).read_text())
    assert sorted(d["element"] for d in two["dashboard"]) == ["Al", "Fe"], \
        "the index did not pick up a page written by a separate run"


def test_an_element_with_no_page_is_a_DECLARED_gap(tmp_path):
    """RYA-833: absence must be visible as absence, never as a blank reading as zero."""
    t = _tracker()
    t.write_element_page(_status(), "Fe", tmp_path)
    index = json.loads(t.write_index(_status(), tmp_path).read_text())
    # reported in the ROSTER's vocabulary (species "CI"), not the page's ("C") --
    # see the note in write_index; the comparison strips the ion, the report does not.
    assert "CI" in index["elements_without_a_page"], (
        "an element in the roster with no page must be listed as missing, so the reader "
        "can tell 'not measured yet' from 'measured and empty'."
    )
