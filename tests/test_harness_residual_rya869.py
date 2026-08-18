"""
tests/test_harness_residual_rya869.py — RYA-869
===============================================
`derive_band_products.py` decided the harness residual with
`is_b = prod.treatment == "ENGINE-B"`, under a comment that stated the rule correctly:
*"The harness residual is a property of the HANDLER that produced the number."* The
comment was right and the code tested a LABEL, so when RYA-798 added the variant
`ENGINE-B-NLTE` — the same `SynthesisHandler` flux fit on a different departure deck —
four published Fe matrix bars were charged the profile fitter's measured residual and
LABELLED `ProfileFitHandler` in their own budget files.

WHAT IS PINNED HERE, AND WHAT IS DELIBERATELY NOT
-------------------------------------------------
Not a number. `assert syst == 0.1700` would have been green through the entire life of
this defect for the ENGINE-B cell beside the broken one, and RYA-845 was pinned in place
by exactly that kind of assertion. What is pinned is the RELATIONSHIPS the defect
violated, each stated so that the PRE-FIX code is red on it:

1. two products of the SAME handler get the SAME harness term whatever their treatment
   names are — the assertion `treatment == "ENGINE-B"` cannot satisfy;
2. the SAME treatment name under two handlers gets two DIFFERENT harness terms — which
   no widening of a treatment set can satisfy, because `1D-LTE` really is a flux fit in
   the near-UV and a profile fit in VIS;
3. the residual and the LABEL are one decision and are handed over together, so a budget
   file cannot name one handler in prose and charge another in arithmetic;
4. nothing may be charged by default: an undeclared handler raises;
5. the banked-artifact recovery table is TOTAL over `band_products.TREATMENTS`, so a
   seventh treatment cannot inherit a handler by falling off the end of a set;
6. every charged residual either equals its handler's own banked control or is declared
   as a known divergence naming the artifact — so "assumed zero" cannot hide inside a
   term whose printed prose says "MEASURED".
"""
from __future__ import annotations

import ast
import re
import sys
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from pipeline import harness_residual as hr                        # noqa: E402
from pipeline.band_products import (                               # noqa: E402
    TREATMENTS, LineMeasurement, Product, build_product)
from pipeline.error_budget import build as build_budget            # noqa: E402


def _lines(treatment: str, n: int = 4) -> list[LineMeasurement]:
    return [LineMeasurement(element="Fe", ion="I", wavelength_air_A=5000.0 + i,
                            instrument="kpno_solar_atlas", ew_mA=float("nan"),
                            ew_method="test", abundance=7.50 + 0.01 * i,
                            treatment=treatment, ew_inversion=False)
            for i in range(n)]


# ── 1. one handler, two treatment names, one harness term ────────────────────────────

def test_engine_b_and_engine_b_nlte_are_charged_identically():
    """THE REGRESSION. `ENGINE-B-NLTE` is `ENGINE-B`'s handler on a different deck.

    RED on the pre-fix code by construction: `prod.treatment == "ENGINE-B"` answers
    False for the NLTE product and True for the LTE one, so the two sides of this
    assertion were 0.0129/ProfileFitHandler and 0.0000/SynthesisHandler — the published
    defect, in one line.
    """
    lte = build_product("Fe", "I", "kpno_solar_atlas", "VIS", "ENGINE-B",
                        _lines("ENGINE-B"), handler="SynthesisHandler")
    nlte = build_product("Fe", "I", "kpno_solar_atlas", "VIS", "ENGINE-B-NLTE",
                         _lines("ENGINE-B-NLTE"), handler="SynthesisHandler")
    assert (hr.for_product(lte).budget_kwargs()
            == hr.for_product(nlte).budget_kwargs()), (
        "two products of the same handler are charged different harness terms — the "
        "residual is following the treatment name again (RYA-869)")


def test_a_new_flux_fit_variant_follows_without_being_listed():
    """A treatment name this module has never heard of still gets its handler's term.

    The ticket's own instruction: do NOT fix it by adding `"ENGINE-B-NLTE"` to the
    equality, because the next variant breaks it again. So the decider is asked with a
    handler and a treatment it does not know, and must still answer.
    """
    made_up = Product(element="Fe", ion="I", instrument="kpno_solar_atlas", band="VIS",
                      treatment="ENGINE-B-SOME-FUTURE-DECK", value=7.5, sigma=0.1,
                      n_lines=4, n_excluded=0, handler="SynthesisHandler")
    assert hr.for_product(made_up).residual_dex == \
        hr.HANDLER_RESIDUAL_DEX["SynthesisHandler"]


# ── 2. one treatment name, two handlers, two harness terms ───────────────────────────

def test_same_treatment_under_two_handlers_is_charged_differently():
    """`1D-LTE` is a flux fit in the near-UV and a profile fit in VIS — really, today.

    This is the assertion that no set of treatment names can satisfy, which is why the
    ticket refused a membership test as the fix. `derive_band_products` builds both of
    these products under the treatment name `1D-LTE`.
    """
    vis = build_product("Fe", "I", "kpno_solar_atlas", "VIS", "1D-LTE",
                        _lines("1D-LTE"), handler="ProfileFitHandler")
    nuv = build_product("Fe", "I", "kpno_solar_atlas", "near-UV", "1D-LTE",
                        _lines("1D-LTE"), handler="SynthesisHandler")
    assert hr.for_product(vis).residual_dex != hr.for_product(nuv).residual_dex
    assert hr.for_product(vis).handler == "ProfileFitHandler"
    assert hr.for_product(nuv).handler == "SynthesisHandler"


# ── 3. the value and the label are ONE decision ──────────────────────────────────────

def test_budget_prose_and_arithmetic_name_the_same_handler():
    """The published artifact said `ProfileFitHandler MEASURED ...` beside 0.0129 on a
    cell the SynthesisHandler produced. The label is part of the decision, so it travels
    with the number rather than being passed separately."""
    for handler in hr.HANDLER_RESIDUAL_DEX:
        kw = hr.for_handler(handler).budget_kwargs()
        assert set(kw) == {"harness_residual_dex", "handler"}, (
            "budget_kwargs must carry the residual AND its label together; splitting "
            "them is how one call site drifts from the other (RYA-845/855)")
        text = build_budget("Fe", 5000.0, 40, scatter_dex=0.10, gf_graded=False,
                            **kw).describe()
        harness_line = [x for x in text.splitlines() if "harness residual" in x]
        assert len(harness_line) == 1
        assert handler in harness_line[0]
        assert f"{kw['harness_residual_dex']:.4f}" in harness_line[0]
        for other in hr.HANDLER_RESIDUAL_DEX:
            if other != handler:
                assert other not in harness_line[0]


# ── 4. nothing is charged by default ─────────────────────────────────────────────────

def test_a_product_with_no_handler_is_refused():
    p = Product(element="Fe", ion="I", instrument="kpno_solar_atlas", band="VIS",
                treatment="ENGINE-B", value=7.5, sigma=0.1, n_lines=4, n_excluded=0)
    with pytest.raises(ValueError, match="declares no handler"):
        hr.for_product(p)


def test_an_unknown_handler_raises_rather_than_returning_zero():
    """A handler with no control has no measured systematic (`MeasurementHandler.
    systematic_dex`'s own rule). Defaulting to 0.0 here would understate every frontier
    bar it touched while printing the word MEASURED next to it."""
    with pytest.raises(KeyError, match="no measured harness residual"):
        hr.for_handler("SomeHandlerNobodyControlled")


def test_build_product_requires_a_handler():
    with pytest.raises(TypeError):
        build_product("Fe", "I", "kpno_solar_atlas", "VIS", "ENGINE-B",
                      _lines("ENGINE-B"))          # noqa: no handler=
    with pytest.raises(KeyError):
        build_product("Fe", "I", "kpno_solar_atlas", "VIS", "ENGINE-B",
                      _lines("ENGINE-B"), handler="NotAHandler")


def test_build_product_records_the_handler_on_the_product():
    p = build_product("Fe", "I", "kpno_solar_atlas", "VIS", "ENGINE-B",
                      _lines("ENGINE-B"), handler="SynthesisHandler")
    assert p.handler == "SynthesisHandler"
    assert "handler" in p.to_frame().columns or p.handler        # per-line frame unchanged


# ── 5. the banked-artifact recovery is TOTAL, and is real handlers only ──────────────

def test_banked_recovery_covers_every_treatment():
    """A seventh treatment must not inherit a handler by falling off the end of a set.

    `_BANKED_PROFILEFIT_HANDLER` is the only inference left in the repo and it is for
    files written before `Product.handler` existed; making it total means adding a
    treatment forces someone to say which handler produces it.
    """
    missing = [t for t in TREATMENTS if t not in hr._BANKED_PROFILEFIT_HANDLER]
    assert not missing, (
        f"treatments with no recorded handler: {missing}. Add them to "
        f"pipeline.harness_residual._BANKED_PROFILEFIT_HANDLER — inheriting one by "
        f"default is the RYA-869 defect.")
    for t in TREATMENTS:
        assert hr.handler_of_banked_cell(route="PROFILEFIT", treatment=t) \
            in hr.HANDLER_RESIDUAL_DEX


def test_banked_recovery_names_registered_handlers_only():
    """Control against the primary: every handler NAME the recovery can return is a real
    registered `MeasurementHandler`, so a typo cannot invent a handler with a residual."""
    from pipeline.measure import HANDLERS
    real = {h.name for h in HANDLERS.values()}
    assert real, "no measurement handlers registered — the control cannot discriminate"
    assert set(hr.HANDLER_RESIDUAL_DEX) <= real, (
        f"residuals are declared for handlers that do not exist: "
        f"{set(hr.HANDLER_RESIDUAL_DEX) - real}")
    assert set(hr._BANKED_PROFILEFIT_HANDLER.values()) <= real


def test_banked_recovery_is_by_route_not_by_engine():
    """An ENGINE-B cell is written into a `*_PROFILEFIT_*` stem — that stem names where
    the LINE SET came from, not who measured it. A route-only rule would call it a
    profile fit, which is the published defect by another road."""
    assert hr.handler_of_banked_cell(route="PROFILEFIT",
                                     treatment="ENGINE-B") == "SynthesisHandler"
    assert hr.handler_of_banked_cell(route="PROFILEFIT",
                                     treatment="1D-LTE") == "ProfileFitHandler"
    assert hr.handler_of_banked_cell(route="SYNTH",
                                     treatment="1D-LTE") == "SynthesisHandler"
    with pytest.raises(KeyError):
        hr.handler_of_banked_cell(route="PROFILEFIT", treatment="ENGINE-Z")


# ── 6. a charged residual is never quietly different from its own control ────────────

def test_every_charged_residual_matches_its_control_or_is_a_declared_divergence():
    """`harness_term()` prints "MEASURED against the known optical answer, not assumed
    zero" beside whatever number it is handed. So a charged value that disagrees with the
    handler's banked control has to be DECLARED, with the artifact named — otherwise the
    prose is an assertion nobody can check.

    🔴 This currently passes BECAUSE of the declaration, not because the numbers agree:
    `SynthesisHandler` is charged 0.0 and its control measured 0.0100 (RYA-770/759).
    """
    for handler, charged in hr.HANDLER_RESIDUAL_DEX.items():
        note = hr.UNCHARGED_CONTROL_RESIDUAL_DEX.get(handler)
        if note is None:
            continue
        p = ROOT / note["control_artifact"]
        assert p.exists(), (
            f"{handler}'s divergence names {note['control_artifact']}, which is absent; "
            f"a declared divergence with no artifact is an unfalsifiable claim")
        d = json.loads(p.read_text())
        assert d["handler"] == handler
        assert abs(abs(float(d["dex_offset"])) - float(note["control_dex"])) < 5e-5, (
            f"{handler}'s declared control residual has drifted from the artifact: "
            f"declared {note['control_dex']}, banked {abs(float(d['dex_offset']))}")
        assert abs(charged - float(note["control_dex"])) > 1e-9, (
            f"{handler} is declared as a DIVERGENCE but the charged value now equals "
            f"its control. Charge it and delete the declaration.")
        assert note["why_not_charged"].strip(), "a divergence must state why"


# ── 7. no caller decides this for itself any more ────────────────────────────────────

_HARNESS_CALLERS = [
    ROOT / "scripts" / "derive_band_products.py",
    ROOT / "scripts" / "rya850_graded_products.py",
    ROOT / "scripts" / "rya836_nearuv_lab_gf_subpool.py",
]


@pytest.mark.parametrize("path", _HARNESS_CALLERS, ids=lambda p: p.name)
def test_no_production_route_states_a_harness_term_at_a_budget_call(path):
    """The pair `harness_residual_dex=` / `handler=` may only reach `error_budget.build`
    through `harness_residual`. RYA-845's double-count and this ticket's misattribution
    are the same shape: a rule written at its call sites drifts between them, and each
    copy is internally consistent while the pair is wrong.

    🔴 IT IS THE BUDGET CALL THAT IS FORBIDDEN, NOT THE WORD. `build_product(...,
    handler="SynthesisHandler")` is the DECLARATION — the route that ran the handler is
    the only code that knows which one it was, and saying so there is the fix. Checked by
    walking the AST rather than grepping, because a text rule could not tell the two
    apart and would have had to permit both.
    """
    tree = ast.parse(path.read_text())
    offenders = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        fn = node.func
        name = fn.attr if isinstance(fn, ast.Attribute) else getattr(fn, "id", "")
        if name not in {"build_budget", "build"}:
            continue
        for kw in node.keywords:
            if kw.arg in {"harness_residual_dex", "handler"}:
                offenders.append((node.lineno, kw.arg))
    assert not offenders, (
        f"{path.name} states the harness term itself at a budget call "
        f"{offenders}; hand it `**harness_residual.for_product(...).budget_kwargs()` "
        f"so the value and its label cannot be separated")


def test_the_guard_is_red_on_the_shape_it_forbids():
    """VERIFY THE TEST DISCRIMINATES (RYA-805/845). A source guard that cannot fail is a
    guard nobody can trust, and this one passes trivially on a file with no budget call
    at all."""
    bad = ast.parse('build_budget("Fe", 5000.0, 4, harness_residual_dex=0.0129,\n'
                    '             handler="ProfileFitHandler")')
    found = [kw.arg for n in ast.walk(bad) if isinstance(n, ast.Call)
             for kw in n.keywords if kw.arg in {"harness_residual_dex", "handler"}]
    assert sorted(found) == ["handler", "harness_residual_dex"]
    good = ast.parse('build_budget("Fe", 5000.0, 4, '
                     '**harness_residual.for_product(p).budget_kwargs())')
    assert not [kw.arg for n in ast.walk(good) if isinstance(n, ast.Call)
                for kw in n.keywords if kw.arg in {"harness_residual_dex", "handler"}]


def test_the_deriver_no_longer_tests_a_treatment_name_for_the_handler():
    """The defect, spelled out, must not come back in any spelling.

    Both the equality and the widened set are rejected: the ticket's instruction was that
    a bigger set rebuilds the same defect one variant later.
    """
    src = (ROOT / "scripts" / "derive_band_products.py").read_text()
    code = "\n".join(l for l in src.splitlines() if not l.lstrip().startswith("#"))
    assert not re.search(r'treatment\s*==\s*"ENGINE-B"', code)
    assert not re.search(r'treatment\s+in\s*[\({\[][^)\]}]*ENGINE-B', code)
    assert "harness_residual.for_product(" in code, (
        "the deriver must ASK the product which handler made it")
