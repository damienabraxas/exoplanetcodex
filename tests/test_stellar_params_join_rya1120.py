"""RYA-1120 — the join between the two uncertainty budgets, and what it refuses.

🔴 WHAT THIS PROTECTS. `uncertainty_stack` derived the stellar-parameter systematics and
`error_budget` published the product bar, and nothing connected them: every live VIS Fe
product carried a `sigma_syst` with no Teff, no log g and no microturbulence in it. These
tests pin the join AND the three refusals that make it honest — no borrowed element-level
number, no partial sum, no applicability read off an engine name.
"""
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

import sys  # noqa: E402

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline import error_budget as EB          # noqa: E402
from pipeline import stellar_params as SP        # noqa: E402


@pytest.fixture(scope="module")
def solar_deltas():
    from pipeline.uncertainty_stack import params_and_deltas
    return params_and_deltas("solar")[1]


def test_delta_p_comes_from_uncertainty_stack_and_is_never_retyped(solar_deltas):
    """The join READS the other budget's numbers; it does not carry its own copy."""
    t = SP.for_product("1D-LTE")
    assert t.deltas == solar_deltas
    # the Sun's three pinned parameters contribute exactly zero and need no run
    assert solar_deltas["logg"] == 0.0 and solar_deltas["feh"] == 0.0
    assert t.required() == ["teff_K", "vturb_kms"]
    src = (ROOT / "pipeline" / "stellar_params.py").read_text()
    assert "0.2912" not in src, (
        "the solar delta_xi has been re-typed into the join — it must come from "
        "uncertainty_stack, or the two budgets drift apart again (RYA-845)")


def test_it_reproduces_RYA_1089s_stamped_Fe_I_number(solar_deltas):
    """An independent cross-check that the join computes the SAME quantity RYA-1089
    stamped: |dA/dvmic| = 0.24 dex/(km/s) x delta_xi -> sigma_B_vmic = 0.0699."""
    t = SP.for_product("1D-LTE", responses={"vturb_kms": 0.24,
                                            "teff_K": 0.0665 / 100.0})
    assert t.sigma_dex() == pytest.approx(0.0699, abs=5e-5)
    c = t.contributions()
    assert c["vturb_kms"]["dex"] == pytest.approx(0.0699, abs=5e-5)
    # and the Sun's Teff is so well pinned that its term is three orders down
    assert c["teff_K"]["dex"] < 1e-3
    assert c["logg"]["dex"] == 0.0 and c["feh"]["dex"] == 0.0


def test_a_pool_with_no_measured_response_gets_NOTHING_not_a_borrowed_number():
    """🔴 The central refusal. RYA-1089's -0.24 was measured on ONE 62-line pool, and
    RYA-1093 showed xi-sensitivity is a property of the LINE SET. A product that has not
    measured its own response must publish a FLOOR, not somebody else's derivative."""
    t = SP.for_product("1D-LTE")           # no responses supplied
    assert t.missing() == ["teff_K", "vturb_kms"]
    assert t.sigma_dex() is None
    assert "NOT MEASURED" in t.source() and "LINE SET" in t.source()


def test_a_partial_measurement_is_not_a_total():
    """Measuring one of two required parameters must not publish their partial sum —
    it would be smaller than the truth and indistinguishable from a complete budget."""
    t = SP.for_product("1D-LTE", responses={"vturb_kms": 0.24})
    assert t.missing() == ["teff_K"]
    assert t.sigma_dex() is None


def test_applicability_is_a_declared_route_not_a_substring():
    """🔴 The RYA-1092 pattern. `"3D" in treatment` swept the <3D> MEAN in with full 3D."""
    assert SP.FULL_3D_TREATMENTS == frozenset({"ENGINE-A-3DNLTE"})
    assert SP.for_product("ENGINE-A-3DNLTE").applicable is False
    for mean3d in ("synth-mean3D-NLTE-gerber-stagger",
                   "synth-mean3D-LTE-gerber-stagger"):
        assert SP.for_product(mean3d).applicable is True, (
            "the <3D> mean was exempted again — RYA-1099 measured it at xi=0 as "
            "+0.137 dex WORSE, because a mean atmosphere averages the velocity "
            "structure OUT and still runs on an inherited xi")


def test_the_full_3D_exemption_is_recorded_with_its_contrary_measurement():
    """An exemption stated bare is indistinguishable from one nobody checked."""
    s = SP.for_product("ENGINE-A-3DNLTE").source()
    assert "NOT APPLICABLE" in s
    assert "+0.0985" in s, (
        "the contrary measurement has been dropped — vmic IS an input axis of the "
        "Amarsi MLP and its correction responds; the exemption rests on the published "
        "value being a SUM whose halves move oppositely, and must say so")


def test_an_exempt_product_may_not_also_carry_responses():
    with pytest.raises(ValueError, match="cannot also carry measured responses"):
        SP.StellarParamSystematic(star_id="solar", deltas={"vturb_kms": 0.3},
                                  responses={"vturb_kms": 0.2}, applicable=False)


def test_a_response_for_an_undeclared_parameter_is_a_STOP():
    """RYA-282 §2: an undeclared delta_p is a stop-and-report, never a zero."""
    with pytest.raises(SP.UnknownParameter, match="no declared delta_p"):
        SP.StellarParamSystematic(star_id="solar", deltas={"vturb_kms": 0.3},
                                  responses={"nonesuch": 1.0})


# ── and the term actually reaches the published budget ───────────────────────
def test_budget_kwargs_flow_into_error_budget_and_reach_sigma_syst():
    t = SP.for_product("1D-LTE", responses={"vturb_kms": 0.24,
                                            "teff_K": 0.0665 / 100.0})
    b = EB.build("Fe", 5500.0, 67, scatter_dex=0.18, gf_graded=True,
                 harness_residual_dex=0.01, handler="SynthesisHandler",
                 harness_provenance="test", **t.budget_kwargs())
    term = next(x for x in b.terms if x.name == "stellar parameters")
    assert term.measured and term.dex == pytest.approx(0.0699, abs=5e-5)
    assert not term.averages_down
    # it is IN the published systematic, not merely attached
    without = EB.build("Fe", 5500.0, 67, scatter_dex=0.18, gf_graded=True,
                       harness_residual_dex=0.01, handler="SynthesisHandler",
                       harness_provenance="test")
    assert b.systematic() > without.systematic()


def test_an_unmeasured_join_makes_the_bar_an_honest_floor_not_a_silent_absence():
    b = EB.build("Fe", 5500.0, 67, scatter_dex=0.18, gf_graded=True,
                 harness_residual_dex=0.01, handler="SynthesisHandler",
                 harness_provenance="test",
                 **SP.for_product("1D-LTE").budget_kwargs())
    names = [t.name for t in b.unmeasured_terms()]
    assert "stellar parameters" in names, (
        "the term vanished instead of being declared unmeasured — that is the F1 defect")
    assert not b.inapplicable_terms()


def test_full_3D_reports_inapplicable_and_its_budget_is_COMPLETE():
    b = EB.build("Fe", 5500.0, 50, scatter_dex=0.18, gf_graded=True,
                 harness_residual_dex=0.01, handler="SynthesisHandler",
                 harness_provenance="test",
                 **SP.for_product("ENGINE-A-3DNLTE").budget_kwargs())
    assert [t.name for t in b.inapplicable_terms()] == ["stellar parameters"]
    assert not b.unmeasured_terms(), (
        "an inapplicable term was counted as a gap — a full-3D bar would read as a "
        "floor forever")
    assert "n/a" in b.describe()
