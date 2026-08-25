"""RYA-1044: the Engine-B leg on the SYNTHESIS route, and the blast radius it must not have.

THE DEFECT, AND WHY IT IS A "PAPER-DONE"
----------------------------------------
`derive_band_products.py` routes `--force-synthesis` straight into `synthesis_route()`,
and that function had **no Engine-B code at all**. So the leg that RYA-784 ("WIRE Engine-B
into derive_band_products") and RYA-798 ("wire the Gerber decks into the production
flux-fit path") both closed as *Done* was unreachable from every synthesis-only band. The
wiring landed in the MAIN route; the bands that can only synthesise take this one.

Measured across every committed `*_products.csv` before this change: **1D-LTE ×12,
ENGINE-A ×7, and not one Engine-B product** — for any element, any band, any deck. Two
tickets, both Done, zero products.

WHAT THESE TESTS PIN
--------------------
The ticket makes one guard mandatory: existing `--force-synthesis` outputs must be
UNCHANGED except for the intended Engine-B addition. Most of these tests are that guard,
because the risk here is not that the new leg is wrong — it is that turning it on quietly
moves a band that was already right.
"""
import ast
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DRIVER = (ROOT / "scripts" / "derive_band_products.py").read_text()
HARNESS = (ROOT / "scripts" / "rya759_nearuv_fe_product.py").read_text()


def _func(src: str, name: str) -> ast.FunctionDef:
    for node in ast.walk(ast.parse(src)):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return node
    raise AssertionError(f"{name} not found")


def _src(src: str, name: str) -> str:
    return ast.get_source_segment(src, _func(src, name)) or ""


# ── the defect is closed ────────────────────────────────────────────────────

def test_synthesis_route_now_has_an_engine_b_leg():
    """The whole ticket: 0 hits for `engine_b_deck`/`ENGINE-B` in this function was the
    defect. It is reachable now."""
    s = _src(DRIVER, "synthesis_route")
    assert "engine_b_deck" in s
    assert "eb_lines = _fit_lines(" in s


def test_the_synthesis_leg_has_NO_measured_EW_precondition():
    """🔴 THE FALSE COUPLING THIS TICKET REMOVES. The main route refuses without a
    profile-fit EW file. A flux fit does not consume an EW — that file is a precondition
    of the EW ROUTE — and requiring one to emit a synthesis product is what made the leg
    unreachable for exactly the bands where synthesis is the only permitted method."""
    s = _src(DRIVER, "synthesis_route")
    assert "no measured EWs" not in s
    assert "measure_band_profilefit" not in s


# ── the blast-radius guard: the mandatory half ──────────────────────────────

def test_the_leg_runs_by_default_and_only_skip_engine_b_stops_it():
    """A leg that produces a product should produce it (Ryan, 2026-08-25). An opt-in flag
    would leave RYA-784/798 half-closed in a NEW way — reachable but never reached — which
    is the same paper-done shape this ticket exists to end.

    ⚠️ It shipped opt-in first, justified as a CI guard, and THAT JUSTIFICATION WAS WRONG:
    `check_result_generators` scores TRACKED artifacts only, so a new output file breaks
    nothing until it is committed. The real cost is RUNTIME — every `--force-synthesis`
    band now fits its lines twice — and `--skip-engine-b` is what that is for."""
    s = _src(DRIVER, "synthesis_route")
    assert "if not a.skip_engine_b:" in s
    assert "synth_engine_b" not in DRIVER, \
        "the opt-in flag is gone -- a vestigial no-op flag is worse than none"


def test_skip_engine_b_still_wins():
    """The one escape hatch. A band that must not pay the second fit still has one, and it
    is the SAME flag the EW route honours — one spelling, both routes."""
    s = _src(DRIVER, "synthesis_route")
    assert "a.skip_engine_b" in s


def test_the_1D_LTE_leg_is_the_same_call_it_always_was():
    """🔴 THE ONE THING THAT MUST NOT MOVE. RYA-759's published near-UV value came off
    this leg. Its body now lives in `_fit_lines` so the Engine-B leg can share it, but the
    1D-LTE call passes NO fit kwargs — so it reaches `fit_one` with exactly the arguments
    it reached before this ticket."""
    s = _src(DRIVER, "synthesis_route")
    assert 'lines = _fit_lines("1D-LTE")' in s, \
        "the 1D-LTE leg must take no fit kwargs at all"


def test_the_engine_b_leg_writes_its_OWN_files_like_ENGINE_A():
    """🔴 THE BLAST-RADIUS DECISION, AND IT FOLLOWS THIS FUNCTION'S OWN PRECEDENT.

    `synthesis_route` already runs THREE legs — 1D-LTE, ENGINE-A and now Engine-B — and
    ENGINE-A writes `{stem}_ENGINE-A_products.csv`, not a second row in the shared file.
    Appending to `{stem}_products.csv` would change the ROW COUNT of a file every existing
    band writes with exactly one row, which is a shift in an existing output. Separate
    files are purely additive: a band that does not opt in is byte-identical, and one that
    does gains files rather than altering them."""
    s = _src(DRIVER, "synthesis_route")
    # ⚠️ matched on the PATH, not on `to_csv(<path>)` -- the lines.csv write wraps across
    # two source lines, so a call-shaped substring misses it and fails for the wrong reason.
    for artifact in ("_products.csv", "_lines.csv", "_budgets.txt", "_provenance.txt"):
        assert 'f"{stem}_{eb_treatment}%s"' % artifact in s, artifact
    # the shared products.csv is written exactly once, by the 1D-LTE leg, unconditionally
    assert s.count('to_csv(out / f"{stem}_products.csv"') == 1


def test_the_shared_products_csv_is_untouched_by_the_leg():
    """The 1D-LTE write must sit OUTSIDE the leg's branch, so opting in cannot move it."""
    s = _src(DRIVER, "synthesis_route")
    head, _, _tail = s.partition("if eb_lines is not None:")
    assert 'to_csv(out / f"{stem}_products.csv"' in head, \
        "the 1D-LTE products.csv must be written before the leg's branch, unconditionally"


# ── the fit harness is SHARED, never copied ─────────────────────────────────

def test_both_legs_go_through_ONE_fit_loop():
    """🔴 RYA-701: one Ba->Al copy of a fitting routine produced thirteen defects. Here
    there is a sharper reason — the two legs are DIFFERENCED against each other, so drift
    between them lands directly in the reported NLTE effect rather than in one product's
    value. A copy could not be wrong in a way that cancels."""
    s = _src(DRIVER, "synthesis_route")
    assert s.count("def _fit_lines(") == 1
    assert s.count("res = fit_one(") == 1, \
        "a second fit_one call site means the loop was copied rather than shared"


def test_fit_one_pass_throughs_default_to_None_and_are_keyword_only():
    """`fit_one` is RYA-759's published harness. The new parameters must be inert when
    unset, or the near-UV value moves by way of an argument nobody passed — and
    keyword-only so no positional caller can acquire one by accident."""
    fn = _func(HARNESS, "fit_one")
    kwonly = {a.arg for a in fn.args.kwonlyargs}
    assert {"nlte_deck", "nlte_deck_key", "atmosphere_layers_file",
            "atmosphere"} <= kwonly
    for d in fn.args.kw_defaults:
        assert isinstance(d, ast.Constant) and d.value is None
    # ...and unset, they must not appear in the call at all
    s = _src(HARNESS, "fit_one")
    assert "_extra = {k: v for k, v in (" in s
    assert "**_extra)" in s


def test_the_atmosphere_override_falls_back_to_ctx():
    """With `atmosphere=None` the fitter must read `ctx['atmosphere']` exactly as before."""
    s = _src(HARNESS, "fit_one")
    assert "ctx['atmosphere'] if atmosphere is None else atmosphere" in s


# ── what the leg emits ──────────────────────────────────────────────────────

def test_the_treatment_tokens_come_from_the_axis_registry():
    """Never retyped. RYA-798 emitted a treatment `TREATMENTS` had never heard of and the
    product died at `build_product` AFTER the synthesis had run."""
    s = _src(DRIVER, "synthesis_route")
    assert "taxes.MEAN3D_NLTE_STAGGER" in s and "taxes.MEAN3D_LTE_STAGGER" in s
    assert "ENGINE-B-MEAN3D" not in s


def test_the_3D_leg_runs_both_gates_before_fitting():
    """Depths must pair index-for-index AND the two τ scales must be the same scale.
    iSpec overwrites the departure τ with the atmosphere's, so a disagreement is applied
    silently rather than raised."""
    s = _src(DRIVER, "synthesis_route")
    assert "assert_depth_match(_dep, _layers)" in s
    assert "assert_tau_consistent(_dep, _model3d)" in s


def test_the_NLTE_label_check_runs_only_where_departures_are_applied():
    """⚠️ iSpec fails SOFT: an unlabelled element lands in `nlte_ignored` and synthesises
    in LTE WITHOUT raising (RYA-764). But demanding NLTE labels of the ⟨3D⟩-LTE comparand
    — a run that deliberately applies none — would refuse it for lacking what it does not
    use, and that comparand is mandatory (RYA-1040)."""
    s = _src(DRIVER, "synthesis_route")
    block = s.split("if _nlte:\n            # iSpec fails SOFT")[1]
    assert "assert_linelist_supports_nlte" in block


def test_the_3D_leg_refuses_an_unregistered_deck_by_name():
    """A missing deck must say which key it looked for and where the aux verdict lives —
    Fe and Mn require the plain aux, the other 15 may use either (RYA-1035)."""
    s = _src(DRIVER, "synthesis_route")
    assert "no <3D> deck registered for" in s
    assert "mean3d_aux_defect_sweep.csv" in s


def test_the_leg_emits_a_separate_product_never_a_correction():
    """RYA-712: engines are separate data products and are never combined. The leg calls
    `build_product` for its own row; it does not adjust the 1D-LTE value."""
    s = _src(DRIVER, "synthesis_route")
    assert "eb_product = build_product(" in s
    # TWO products, built independently -- not one product adjusted into another.
    # THREE legs on this route -- 1D-LTE, ENGINE-A, Engine-B -- each building its own
    # product independently rather than adjusting another's.
    assert s.count("build_product(") == 3
    # ⚠️ Checked on ASSIGNMENT, not on the characters. The leg PRINTS
    # `eb_product.value - product.value`, which is the differential RYA-1040 exists to
    # report; a substring ban on "product.value -" forbids the diagnostic rather than the
    # correction. What must never happen is a WRITE back into either product.
    fn = _func(DRIVER, "synthesis_route")
    for node in ast.walk(fn):
        if isinstance(node, ast.Assign):
            for tgt in node.targets:
                if isinstance(tgt, ast.Attribute) and tgt.attr in ("value", "sigma"):
                    raise AssertionError(
                        f"a product field is assigned at line {node.lineno} -- engines "
                        f"are separate products and are never combined (RYA-712)")


@pytest.mark.parametrize("deck,expected", [
    ("gerber-mean3d", "synth-mean3D-NLTE-gerber-stagger"),
    ("gerber-mean3d-lte", "synth-mean3D-LTE-gerber-stagger"),
])
def test_the_deck_choice_selects_the_matching_treatment(deck, expected):
    """The CLI choice and the emitted treatment must agree, or a product is labelled as
    something it is not."""
    from pipeline import treatment_axes as T
    nlte = deck in ("gerber-nlte", "gerber-mean3d")
    token = (T.MEAN3D_NLTE_STAGGER if nlte else T.MEAN3D_LTE_STAGGER).token
    assert token == expected
