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
    """The four DEEPGRADED near-UV Fe I products RYA-1209 unblocked.

    🔴 SCOPED TO THE TIER, NOT WIDENED. RYA-1213 added Reference Grade rows in this band
    and they are a DIFFERENT POOL: Reference applies no depth gate, so it also takes
    3026.056 A -- the one lab line below the 0.05 depth floor -- which Deep excludes by
    construction. Its n is 56, not 54/55, and asserting this ticket's pool size on it
    would be asserting the wrong claim about the right band. The properties that ARE
    about the band rather than the pool are asserted over every near-UV Fe I product; see
    `nearuv_fe1_all`.
    """
    return [p for p in feed["products"] if p["band"] == "near-UV" and p["ion"] == "I"
            and p["tier"] == "DEEPGRADED"]


@pytest.fixture
def nearuv_fe1_all(feed):
    """EVERY near-UV Fe I product, whatever its pool — for band-level claims."""
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


def test_EVERY_near_uv_fe_I_product_is_below_its_pre_molecular_value(nearuv_fe1_all):
    """RYA-1213 — the molecular correction is a property of the BAND, so it must reach
    every pool measured in it, not only the four this ticket published. The Reference
    rows are a different pool and a different n; what they share with the Deep rows is
    the opacity, and that is what this asserts."""
    pre = {"solar_kpno_kurucz2005_corrected": {"1D-LTE": 7.642, "ENGINE-A": 7.651},
           "solar_kpno_molecfit_corrected": {"1D-LTE": 7.596, "ENGINE-A": 7.606}}
    for p in nearuv_fe1_all:
        before = pre[p["holding"]][p["treatment"]]
        assert p["A"] < before - 0.05, (
            f"{p['holding']}/{p['treatment']} tier={p['tier']} is {p['A']}, not below "
            f"the pre-molecular {before} — molecular opacity did not reach this pool")


def test_the_pool_lost_exactly_one_line(nearuv_fe1):
    """58 -> 57 candidates. The fitted n is lower again because the RYA-1191 validity bound
    drops non-convergent fits, which is a different mechanism and must not be conflated."""
    for p in nearuv_fe1:
        assert p["n_lines"] in (54, 55), f"unexpected n_lines {p['n_lines']}"


def test_every_near_uv_product_states_what_it_now_contains(feed):
    n = [p for p in feed["products"] if p["band"] == "near-UV"]
    # 🔴 10, NOT 8. RYA-1208 added the Fe II `synth-1D-LTE-gerber` leg on both near-UV KP
    # holdings. They are Fe II, so the K07 drop this ticket is about does not touch them
    # (K07 was an Fe I line) -- but they ARE near-UV products and the per-product
    # assertions below apply to them like any other, which is why they are swept in
    # rather than filtered out.
    # 🔴 RYA-1213 — THE PIN MOVES TO THE DEPTH-SPLIT TIERS AND THE LOOP COVERS ALL.
    # The REFERENCE tier adds six near-UV Fe II products, and they must carry the opacity
    # metadata exactly like every other near-UV row -- which is what the loop below
    # asserts, over all of them. Re-pinning the total each time a Reference cell lands
    # would turn a vanish-detector into a number that gets bumped without being read, so
    # the literal now pins the ten Codex/Deep rows it was written about.
    assert len([q for q in n if q["tier"] != "REFERENCE"]) == 10
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
