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


#: 🔴 BLOCKED BY RYA-1080's GUARD, AND NOT WORKED AROUND.
#: RYA-1088 requires sigma_params to be recorded IN data/products/solar/Fe.json.
#: RYA-1080's `test_no_published_value_was_edited_to_reconcile` asserts every published
#: field of every product is byte-for-byte what origin/main published -- "copied_to is the
#: only key this ticket may write". Its DOCSTRING scopes that to RYA-1080; its
#: IMPLEMENTATION is global, so no ticket can ever add a field to a product.
#:
#: Stamping the fields makes RYA-1080's guard fail. Measured: clean origin/main fails 1
#: test (the pre-existing rya1080 generator gate); with the stamp applied the branch fails
#: 2. That is a regression I introduced.
#:
#: Widening another ticket's CRITICAL guard to make my own change pass is the exact
#: antipattern this project ratifies against, so the stamp is REVERTED and these tests are
#: xfail rather than deleted -- they flip to passing the moment the conflict is resolved,
#: and deleting them would erase the requirement.
#:
#: Resolution is Ryan's: either (a) make RYA-1080's guard precise -- refuse CHANGED values,
#: allow ADDED keys, which preserves its stated intent ("reconcile the artifacts, not the
#: numbers") -- or (b) record sigma_params somewhere other than the product JSON, which
#: contradicts RYA-1088's own critical condition.
_BLOCKED = ("RYA-1080's published-field guard forbids adding any key to "
            "data/products/solar/Fe.json; see the note above -- Ryan's call")


@pytest.fixture(scope="module")
def mean3d_products():
    d = json.loads(PRODUCT.read_text())
    ps = [p for p in d["products"] if p.get("treatment") in MEAN3D]
    assert ps, "no ⟨3D⟩ product in the solar Fe record"
    return ps


@pytest.mark.xfail(reason=_BLOCKED, strict=True)
def test_sigma_params_is_PRESENT_not_absent(mean3d_products):
    """The whole point. `None` and `missing` are different from a measured value."""
    for p in mean3d_products:
        assert "sigma_params" in p, f"{p['treatment']} has NO parameter term"
        assert p["sigma_params"] is not None
        assert p.get("sigma_params_reason"), "a zero without a reason is still an absence"


@pytest.mark.xfail(reason=_BLOCKED, strict=True)
def test_sigma_reported_is_attached_and_correctly_combined(mean3d_products):
    """sigma_reported = sqrt(sigma_stat^2 + sigma_params^2) — the RYA-282 convention."""
    import math
    for p in mean3d_products:
        exp = math.sqrt(p["sigma_stat"] ** 2 + p["sigma_params"] ** 2)
        assert p["sigma_reported"] == pytest.approx(exp, abs=5e-4)


@pytest.mark.xfail(reason=_BLOCKED, strict=True)
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


def test_the_solar_delta_xi_is_SOURCED_not_a_literal():
    """🔴 RYA-1089. The solar delta_xi was 0.05 from the RYA-158 spec line
    "vturb = 0.9 +/- 0.05" -- no citation, and doubly stale, since that line's central
    value is 0.9 against the adopted 1.0 (RYA-196). It drove the ENTIRE solar sigma_params.
    It is now Jofré+2014 Table 2 sigma_vmic = 0.18, read from the record beside its
    citation, and an ABSENT value refuses rather than becoming a number again."""
    assert float(STAR_PARAMS["solar"]["e_xi"]) == pytest.approx(0.18)
    _, d = params_and_deltas("solar")
    assert d["vturb_kms"] == pytest.approx(0.18)
    src = (ROOT / "pipeline" / "uncertainty_stack.py").read_text()
    assert "'vturb_kms': 0.05" not in src, "the bare literal must not come back"


def test_the_sourced_delta_xi_made_the_bar_WORSE_not_better():
    """⚠️ RYA-161 evidence, as a test. The sourced value is 3.6x the literal it replaced
    and moves sigma_reported TOWARDS the 0.05 gate (0.0214 -> 0.0467). A delta_p chosen to
    keep a bar inside a gate would have moved it the other way."""
    _, d = params_and_deltas("solar")
    assert d["vturb_kms"] > 0.05, "the sourced value is larger than the literal it retired"


@pytest.mark.parametrize("star,e_xi", [("procyon", 0.11),
                                       ("alpha_cen_a", 0.07),
                                       ("alpha_cen_b", 0.31),
                                       ("solar", 0.18)])
def test_per_star_delta_p_are_declared_and_sourced(star, e_xi):
    """Every delta_p from a cited source, never typed from memory (the ticket's firewall).
    The xi uncertainties are Jofré+2014 Table 2 σvmic, verified by extracting the paper's
    own text rather than a search summary."""
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
