"""RYA-1088 — a measured zero must be sayable, and an undeclared delta_p must not be one.

🔴 THE CONFUSION THIS ENDS. For the Sun the Type B term is ~0 by construction, so the
⟨3D⟩ product carried no parameter term at all — and an ABSENT term reads as an omission.
A readiness checklist reading "sigma_params never measured" on the solar product
manufactured a false alarm that way, while the term had existed since RYA-158.

RYA-907 ratified one half of this: UNMEASURED IS NOT ZERO. This is the mirror —
A MEASURED ZERO IS NOT AN ABSENCE. Both states have to be sayable, and these tests pin
that the product can say each.
"""
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PRODUCT = ROOT / "data" / "products" / "solar" / "Fe.json"
MEAN3D = ("synth-mean3D-NLTE-gerber-stagger", "synth-mean3D-LTE-gerber-stagger")

from config.constants import STAR_PARAMS                      # noqa: E402
from pipeline.uncertainty_stack import params_and_deltas       # noqa: E402


#: 🔴 UNBLOCKED 2026-08-28 (RYA-1089) -- AND THE BLOCK WAS REAL WHILE IT LASTED.
#: RYA-1088 requires sigma_params to be recorded IN data/products/solar/Fe.json. RYA-1080's
#: `test_no_published_value_was_edited_to_reconcile` USED TO assert that every published
#: field of every product was byte-for-byte what `origin/main` published. Its DOCSTRING
#: scoped that to RYA-1080; its IMPLEMENTATION was global, so no ticket could add a field
#: to a product at all. RYA-1088 refused to widen another ticket's CRITICAL guard to pass
#: its own change, reverted the stamp, and left these three tests as STRICT xfail -- the
#: requirement kept on the record rather than deleted, ready to flip.
#:
#: They flipped. RYA-1100 rewrote that guard along exactly the line RYA-1088 asked for --
#: its option (a), refuse CHANGED values but allow ADDED keys -- and its
#: `published_value_edits` now says so in its own docstring, naming RYA-1088's
#: `sigma_params` as the case that forced the fix. VERIFIED, not assumed: with the stamp
#: applied, `tests/test_feed_repo_reconciliation_rya1080.py` is 15/15 green and these three
#: XPASS(strict). The xfail markers are therefore RETIRED, not the tests.
#:
#: ⚠️ The lesson worth keeping: the xfail REASON went stale on `main` for several merges.
#: RYA-1100 removed the blocker; nobody flipped the tests that named it, so the repo kept
#: asserting "blocked, Ryan's call" about a conflict that no longer existed. A test parked
#: against another ticket's guard has to be re-checked when that guard moves.
@pytest.fixture(scope="module")
def mean3d_products():
    d = json.loads(PRODUCT.read_text())
    ps = [p for p in d["products"] if p.get("treatment") in MEAN3D]
    assert ps, "no ⟨3D⟩ product in the solar Fe record"
    return ps


def test_sigma_params_is_PRESENT_not_absent(mean3d_products):
    """The whole point. `None` and `missing` are different from a measured value."""
    for p in mean3d_products:
        assert "sigma_params" in p, f"{p['treatment']} has NO parameter term"
        assert p["sigma_params"] is not None
        assert p.get("sigma_params_reason"), "a zero without a reason is still an absence"


def test_sigma_reported_is_attached_and_correctly_combined(mean3d_products):
    """sigma_reported = sqrt(sigma_stat^2 + sigma_params^2) — the RYA-282 convention."""
    import math
    for p in mean3d_products:
        exp = math.sqrt(p["sigma_stat"] ** 2 + p["sigma_params"] ** 2)
        assert p["sigma_reported"] == pytest.approx(exp, abs=5e-4)


def test_the_per_parameter_terms_are_recorded_so_the_zero_is_LEGIBLE(mean3d_products):
    """0.012 alone does not tell a reader that logg and [Fe/H] contributed exactly nothing
    BY DEFINITION rather than by accident."""
    for p in mean3d_products:
        t = p["sigma_params_terms"]
        assert t["logg"] == 0.0 and t["FeH"] == 0.0
        assert t["vmic"] > 0 and t["Teff"] > 0


def test_solar_logg_and_FeH_deltas_are_ZERO_by_definition():
    """🔴 CRITICAL per the ticket: a non-zero solar sigma_B for logg or [Fe/H] is wrong.
    The Sun DEFINES the [Fe/H] zero point and its logg is fixed by IAU mass and radius."""
    _, d = params_and_deltas("solar")
    assert d["logg"] == 0.0
    assert d["feh"] == 0.0


def test_the_solar_e_feh_in_stars_yaml_is_NOT_used_as_a_delta():
    """⚠️ THE TRAP. `stars.yaml` carries `e_feh: 0.03` for the Sun — that is the
    uncertainty on the solar iron SCALE, not on the Sun's [Fe/H] relative to itself.
    Reading it would produce a non-zero solar sigma_B([Fe/H])."""
    assert float(STAR_PARAMS["solar"]["e_feh"]) == 0.03      # it is there
    _, d = params_and_deltas("solar")
    assert d["feh"] == 0.0                                   # and it is NOT used


def test_the_solar_delta_xi_is_MEASURED_not_a_literal_and_not_borrowed():
    """🔴 RYA-1089 THEN RYA-311. Two retirements, and the second went further than the first.

    The original 0.05 came from the RYA-158 spec line "vturb = 0.9 +/- 0.05" -- no citation,
    and doubly stale, since that line's central value is 0.9 against the adopted 1.0
    (RYA-196). RYA-1089 replaced it with Jofré+2014 Table 2 sigma_vmic = 0.18: honest
    provenance, but the wrong QUANTITY -- an inter-node dispersion across seven GBS
    pipelines, not the formal error of ours.

    RYA-311 measured it: the Delta-chi2 = 1 error on our own Fe I reduced-EW FITEXY slope,
    propagated through the solve. Neither borrowed number survives for the Sun; Jofré's
    per-star sigma_vmic stays the interim delta_p for stars this solve has not been run on.
    """
    assert float(STAR_PARAMS["solar"]["e_xi"]) == pytest.approx(0.0588)
    _, d = params_and_deltas("solar")
    assert d["vturb_kms"] == pytest.approx(0.0588)
    src = (ROOT / "pipeline" / "uncertainty_stack.py").read_text()
    assert "'vturb_kms': 0.05" not in src, "the bare literal must not come back"
    y = (ROOT / "config" / "stars.yaml").read_text()
    assert "RYA-311" in y and "FITEXY" in y, "the measured value must say where it came from"


def test_the_measured_delta_xi_RELAXES_the_bar_and_the_larger_candidates_are_ON_THE_RECORD():
    """⚠️ RYA-161 evidence, and this one cuts the uncomfortable way -- which is why it is
    a test and not a paragraph.

    RYA-1089's sourced 0.18 moved sigma_reported TOWARDS the 0.05 gate (0.0214 -> 0.0467),
    and that was easy to defend: nobody picks a delta_p that makes their own bar worse.
    RYA-311's measured 0.0588 moves it back AWAY (0.0467 -> 0.0227), which is exactly the
    direction a tuned number would move. Two things make it not that, and both are checked
    here rather than asserted in a comment:

      1. The METHOD was fixed by the ticket before the number existed -- the formal FITEXY
         fit error, no alternative offered.
      2. The SAME measurement produced two LARGER candidates -- the 0.2695 km/s ceiling
         sensitivity and the 0.3107 chi2-scaled error, either of which puts the solar bar
         OUTSIDE its gate. They are recorded beside e_xi instead of being discarded, so the
         smaller number is not the only one on the record and the choice between them stays
         visible and reversible.
    """
    _, d = params_and_deltas("solar")
    assert d["vturb_kms"] < 0.18, "RYA-311's measured error is smaller than Jofré's"
    rec = STAR_PARAMS["solar"]
    assert float(rec["e_xi_selection"]) > float(rec["e_xi"]), \
        "the ceiling sensitivity must be on the record beside the formal error"
    assert float(rec["e_xi_chi2_scaled"]) > float(rec["e_xi"]), \
        "the chi2-scaled alternative must be on the record beside the formal error"


def test_the_measured_solar_xi_is_recorded_even_though_the_pin_stands():
    """RYA-311 solved xi = 0.7088 and did NOT adopt it -- moving to it raises solar
    A(Fe I) by +0.0595 dex, which is a RYA-196/register reconciliation, not a silent edit.
    A measured value that disagrees with the pin must still be WRITTEN DOWN; dropping it
    because it was not adopted would leave the pin looking unchallenged."""
    rec = STAR_PARAMS["solar"]
    assert float(rec["xi"]) == pytest.approx(1.0), "the pin is unchanged"
    assert float(rec["xi_measured_rya311"]) == pytest.approx(0.7088)
    assert "measured, RYA-311" in str(rec["xi_provenance"])


@pytest.mark.parametrize("star,e_xi", [("procyon", 0.11),
                                       ("alpha_cen_a", 0.07),
                                       ("alpha_cen_b", 0.31),
                                       ("solar", 0.0588)])
def test_per_star_delta_p_are_declared_and_sourced(star, e_xi):
    """Every delta_p from a cited source, never typed from memory (the ticket's firewall).
    The three off-solar xi uncertainties are Jofré+2014 Table 2 σvmic, verified by
    extracting the paper's own text rather than a search summary. The SOLAR one is no
    longer Jofré's: RYA-311 measured it on our own data, and the two provenances are
    deliberately different -- Jofré's stays the interim value for stars the solve has not
    yet been run on."""
    rec = STAR_PARAMS[star]
    for f in ("e_teff", "e_logg", "e_feh", "e_xi"):
        assert rec.get(f) is not None, f"{star} is missing {f}"
    assert float(rec["e_xi"]) == pytest.approx(e_xi)
    _, d = params_and_deltas(star)
    assert d["vturb_kms"] == pytest.approx(e_xi)


def test_the_citation_travels_with_the_value():
    """A sourced number whose source is not written next to it is not sourced."""
    y = (ROOT / "config" / "stars.yaml").read_text()
    assert y.count("Jofré+2014 (A&A 564, A133) Table 2") >= 3
    assert "standard deviation of this mean" in y


@pytest.mark.parametrize("star", ["tau_ceti", "eps_eri", "55cnc_a"])
def test_a_star_with_NO_declared_delta_p_still_REFUSES(star):
    """The fence moved from 'is this the Sun' to 'is this sourced'. An undeclared
    uncertainty must raise, never become a silent zero (RYA-907)."""
    with pytest.raises(NotImplementedError) as e:
        params_and_deltas(star)
    assert "CITED" in str(e.value) or "cited" in str(e.value)


def test_declared_stars_no_longer_raise():
    """The other half: the fence must actually be open for the run-order targets."""
    for star in ("procyon", "alpha_cen_a", "alpha_cen_b"):
        p, d = params_and_deltas(star)
        assert all(v is not None for v in d.values())
        assert d["teff_K"] > 0, f"{star} Teff delta should be a real uncertainty"


# ── RYA-1089: the gate-margin artifact must not drift away from the record ──────────
def test_the_gate_margin_artifact_tracks_the_RECORDED_delta_xi():
    """A table of consequences is only useful while it describes the value actually carried.

    🔴 The artifact names a `recorded_e_xi` and marks that row. If someone changes
    `solar.e_xi` in stars.yaml and does not re-run the generator, the artifact keeps
    asserting a gate verdict for a delta_xi the budget no longer uses -- a status surface
    that stopped reading the thing it reports on (RYA-1064). This fails loudly instead.
    """
    import json
    from pathlib import Path
    p = Path(__file__).resolve().parent.parent / "data/results/rya1089/rya1089_delta_xi_gate_margin.json"
    if not p.exists():
        pytest.skip("gate-margin artifact not generated in this tree")
    art = json.loads(p.read_text())
    assert float(art["recorded_e_xi"]) == pytest.approx(float(STAR_PARAMS["solar"]["e_xi"])), (
        "rya1089_delta_xi_gate_margin.json was generated against a different solar e_xi "
        "than stars.yaml now carries -- re-run scripts/rya1089_delta_xi_gate_margin.py")
    assert art["adopted_xi_kms"] == pytest.approx(float(STAR_PARAMS["solar"]["xi"])), (
        "the artifact's adopted xi no longer matches the pin it was computed against")
    # the honest finding itself: the candidates must still straddle the gate, or the
    # conclusion reported to Ryan has silently changed.
    verdicts = {m["verdict"] for c in art["candidates"] for m in c["members"].values()}
    assert verdicts == {"INSIDE", "OVER"}, (
        f"the gate verdict is no longer straddled ({verdicts}) -- the RYA-1089 finding "
        f"has changed and the ticket's conclusion needs restating")
