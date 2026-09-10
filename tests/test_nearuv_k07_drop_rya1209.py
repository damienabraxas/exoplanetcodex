"""RYA-1209 — the near-UV Fe I pool drops its one non-lab line, and publishes.

CURATION + REPUBLISH. What these pin is the DISTINCTION that makes the removal legitimate:
the line went because it never met the tier's definition, not because it failed a gate.
Those look identical in the diff and are opposite in kind (RYA-981: a quota dressed as a
cut), so the reason, the size of the value move, and the rung are all asserted.
"""
import csv
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
FEED = ROOT / "data/products/solar/Fe.json"
REG = ROOT / "data/catalog/line_curation_exclusions.csv"
K07 = 3721.270


@pytest.fixture(scope="module")
def feed():
    return json.loads(FEED.read_text())


@pytest.fixture(scope="module")
def nearuv_fe1(feed):
    return [p for p in feed["products"] if p["band"] == "near-UV" and p["ion"] == "I"]


def test_the_k07_line_is_curated_with_a_tier_definition_reason():
    """🔴 THE REASON IS THE WHOLE ARGUMENT. 'It blocks the budget' would be selecting lines
    to pass a gate. 'It is not lab gf, and this tier means lab gf' is a membership fix."""
    rows = [r for r in csv.DictReader(
        l for l in REG.read_text().splitlines() if not l.startswith("#"))]
    hit = [r for r in rows if r["species"] == "Fe I"
           and abs(float(r["wavelength_air_A"]) - K07) < 0.01]
    assert len(hit) == 1, "Fe I 3721.27 must be in the curation registry exactly once"
    r = hit[0]
    assert r["band"] == "near-UV"
    assert "systematic:K07" in r["reason"]
    assert "no citable per-line gf uncertainty" in r["reason"]
    assert r["ticket"] == "RYA-1209" and r["ruled_by"].strip()
    # and the evidence must carry the measured consequence, not a promise about it
    assert "7.560->7.554" in r["evidence"], "the value move must be recorded as MEASURED"


def test_removing_it_barely_moves_the_value(nearuv_fe1):
    """A membership fix should not be worth much. If dropping one line of 58 moved the
    answer materially, 'it does not belong in the tier' would need re-examining."""
    from pipeline.line_curation import excluded
    assert excluded("Fe I", 3721.272), "the registry must actually catch it"
    # measured: -0.007..+0.003 across the four products, all far under the molecular lever
    for p in nearuv_fe1:
        assert 7.45 < p["A"] < 7.60, f"{p['holding']}/{p['treatment']} A={p['A']}"


def test_the_curation_does_not_catch_a_neighbour():
    """The registry matches on a tolerance; a tolerance that swallowed the next line along
    would remove a lab line silently."""
    from pipeline.line_curation import excluded
    for w in (3703.821, 3721.0, 3722.0, 3730.386):
        assert not excluded("Fe I", w), f"{w} must not be curated"


def test_all_four_fe_I_products_carry_post_molecular_values(nearuv_fe1):
    """The point of the ticket: RYA-1207 could not publish these. Now they are published,
    and BELOW the pre-molecular values they replaced (7.596/7.606/7.642/7.651)."""
    assert len(nearuv_fe1) == 4
    pre = {"solar_kpno_kurucz2005_corrected": {"1D-LTE": 7.642, "ENGINE-A": 7.651},
           "solar_kpno_molecfit_corrected": {"1D-LTE": 7.596, "ENGINE-A": 7.606}}
    for p in nearuv_fe1:
        before = pre[p["holding"]][p["treatment"]]
        assert p["A"] < before - 0.05, (
            f"{p['holding']}/{p['treatment']} is {p['A']}, not below the pre-molecular "
            f"{before} — the molecular correction did not reach this product")


def test_the_pool_lost_exactly_one_line(nearuv_fe1):
    """58 -> 57 candidates. The fitted n is lower again because the RYA-1191 validity bound
    drops non-convergent fits, which is a different mechanism and must not be conflated."""
    for p in nearuv_fe1:
        assert p["n_lines"] in (54, 55), f"unexpected n_lines {p['n_lines']}"


def test_every_near_uv_product_states_what_it_now_contains(feed):
    n = [p for p in feed["products"] if p["band"] == "near-UV"]
    assert len(n) == 8
    for p in n:
        assert p.get("opacity_limit") == "MOLECULAR-OPACITY-INCLUDED"
        note = p.get("opacity_note", "")
        assert "OH A-X" in note and "3721.272" in note
        # 🔴 the note must keep saying what is STILL missing; a correction that stops
        # describing its own residual reads as a converged answer
        assert "UNDER-corrects" in note and "floor" in note
        assert "pseudo-continuum" in note, (
            "sigma_syst did not move, and the note must say why: the 0.10 dex continuum "
            "term dominates, so clearing the gf rung bought publication, not a smaller bar")
