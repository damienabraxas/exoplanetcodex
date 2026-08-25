"""RYA-1040: the ⟨3D⟩ route through the driver, the fit and the shared generator.

Three layers had to learn the same two facts, and each had its own version of the same
defect:

  * `derive_band_products` resolved the deck from the ELEMENT (`"Fe" if a.element == "Fe"
    else a.element` — an identity expression), so it picked `Fe` and could never reach
    `Fe@mean3D`;
  * `_fit_synth_flux` did the same one layer down, so fixing only the driver would have
    left the fit running 1D departures under a ⟨3D⟩ label — a well-formed product wrong
    in exactly the axis it claims to be about;
  * `_synth_flux_at_abund` had no way to hand iSpec an atmosphere as a FILE, which is the
    RYA-798 shape one axis over: the one generator both routes call had no parameter for
    the thing that needed passing.

These tests do not synthesise. They pin the CONTRACT each layer now honours, and above
all that the defaults are untouched — RYA-770 stabilised the LTE path at −0.026 dex
against a banked answer and it must not move.
"""
import ast
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

# ⚠️ READ AS SOURCE, NOT IMPORTED. `pipeline.abundances_derive` imports iSpec at module
# level, which exists only on Sirius -- so importing it here would make this file SKIP on
# the Mac and run only in CI. That is precisely how RYA-1040's own predecessor shipped two
# assertions nobody had ever executed (PR #384, caught by CI on the 7.46 -> 7.50 change).
# Parsing the source instead means these run EVERYWHERE, which is what a contract test
# about signatures and call shapes actually needs.
_SRC = {name: (ROOT / rel).read_text()
        for name, rel in (("derive", "pipeline/abundances_derive.py"),
                          ("handler", "pipeline/measure/synthesis.py"),
                          ("driver", "scripts/derive_band_products.py"))}


def _func(module_key: str, name: str) -> ast.FunctionDef:
    tree = ast.parse(_SRC[module_key])
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return node
    raise AssertionError(f"{name} not found in {module_key}")


def _defaults(fn: ast.FunctionDef) -> dict:
    """arg name -> default node, for keyword-or-positional args with defaults."""
    args = fn.args.args + fn.args.kwonlyargs
    defaults = ([None] * (len(fn.args.args) - len(fn.args.defaults))
                + list(fn.args.defaults) + list(fn.args.kw_defaults))
    return {a.arg: d for a, d in zip(args, defaults)}


def _src_of(module_key: str, name: str) -> str:
    return ast.get_source_segment(_SRC[module_key], _func(module_key, name)) or ""


# ── the shared generator: a new parameter that is inert by default ──────────

def test_the_generator_accepts_an_atmosphere_FILE():
    """A mul23 ⟨3D⟩ model has five columns and iSpec's MARCS writer needs Depth, Pₑ, P_g
    and P_rad — so for a ⟨3D⟩ atmosphere there is nothing for `write_atmosphere` to write.
    The file is the only route in."""
    d = _defaults(_func("derive", "_synth_flux_at_abund"))
    assert "atmosphere_layers_file" in d
    assert isinstance(d["atmosphere_layers_file"], ast.Constant)
    assert d["atmosphere_layers_file"].value is None


def test_the_LTE_call_is_unchanged_when_no_file_is_given():
    """🔴 THE INVARIANT RYA-770 BOUGHT. The parameter is splatted from a dict that is
    EMPTY when unused, so with `atmosphere_layers_file=None` both calls are exactly the
    calls they were before it existed — not "the same plus a default-valued keyword".

    `test_rya798_gerber_nlte_wiring` pins the NLTE half of this; this pins the half that
    RYA-1040 could have broken."""
    src = _src_of("derive", "_synth_flux_at_abund")
    assert "_atm_file = ({} if atmosphere_layers_file is None" in src
    # both branches splat the same dict — neither hard-codes the keyword
    assert src.count("**_atm_file") == 2


def test_the_fit_takes_a_deck_KEY_and_defaults_it_to_the_element():
    """🔴 THE DECK IS NAMED, NOT INFERRED. Defaulting to the element keeps every existing
    caller byte-identical while making `<El>@mean3D` REACHABLE, which it was not."""
    d = _defaults(_func("derive", "_fit_synth_flux"))
    for k in ("nlte_deck_key", "atmosphere_layers_file"):
        assert isinstance(d[k], ast.Constant) and d[k].value is None, k
    src = _src_of("derive", "_fit_synth_flux")
    assert "_deck_key = nlte_deck_key or element" in src
    # and every deck lookup goes through it — none may still key on `element`
    assert "for_node(element," not in src
    assert "has_abundance_axis(element)" not in src
    assert src.count("_deck_key") >= 4


def test_the_handler_forwards_both_without_requiring_them():
    """`.get` rather than `[...]`: a context that has never heard of the ⟨3D⟩ route must
    keep working, because every 1D product still flows through this handler."""
    src = _SRC["handler"]
    assert 'nlte_deck_key=context.get("nlte_deck_key")' in src
    assert 'atmosphere_layers_file=context.get("atmosphere_layers_file")' in src


# ── the driver: the pair, and the deck it selects ───────────────────────────

def test_the_cli_offers_both_members_of_the_pair():
    """The comparand is not optional. A ⟨3D⟩-NLTE number whose only comparand is 1D-LTE
    reports the 1D→mean-3D ATMOSPHERE shift as non-LTE physics (RYA-542)."""
    src = _SRC["driver"]
    assert '"gerber-mean3d", "gerber-mean3d-lte"' in src
    # the LTE member takes the <3D> branch too -- same atmosphere, departures withheld
    assert 'mean3d = a.engine_b_deck in ("gerber-mean3d", "gerber-mean3d-lte")' in src
    assert 'nlte = a.engine_b_deck in ("gerber-nlte", "gerber-mean3d")' in src


def test_the_driver_builds_the_deck_key_and_never_infers_it():
    """The identity expression `"Fe" if a.element == "Fe" else a.element` is gone.

    ⚠️ CHECKED ON THE CODE, NOT ON THE TEXT. That expression still appears in the file --
    in a COMMENT, documenting what it used to be, which is exactly where it should still
    appear. A substring test would either fail on the comment or force deleting the
    history to satisfy the test. So the AST is walked instead: no `for_node` call in the
    driver may take a conditional as its deck argument."""
    src = _SRC["driver"]
    assert 'deck_key = f"{a.element}@mean3D"' in src
    calls = [n for n in ast.walk(ast.parse(src))
             if isinstance(n, ast.Call)
             and isinstance(n.func, ast.Attribute) and n.func.attr == "for_node"]
    assert calls, "the driver must still call for_node somewhere"
    for c in calls:
        assert not isinstance(c.args[0], ast.IfExp), \
            "the deck argument is a conditional -- it is being inferred, not selected"


def test_the_driver_refuses_an_unregistered_3D_deck_by_name():
    """A missing deck must say WHICH key it looked for and where the aux verdict lives —
    Fe and Mn require the plain aux, the other 15 may use either (RYA-1035)."""
    src = _SRC["driver"]
    assert "no <3D> deck registered for" in src
    assert "mean3d_aux_defect_sweep.csv" in src


def test_both_gates_run_on_the_3D_leg():
    """Depths must pair index-for-index AND the two τ scales must be the same scale —
    iSpec overwrites the departure τ with the atmosphere's, so a disagreement is applied
    silently rather than raised."""
    src = _SRC["driver"]
    assert "assert_depth_match(dep, layers)" in src
    assert "assert_tau_consistent(dep, model3d)" in src


def test_the_LTE_member_withholds_departures_and_skips_the_NLTE_LABEL_CHECK():
    """⚠️ `assert_linelist_supports_nlte` belongs ONLY to the NLTE member. Running it on
    the comparand would demand NLTE labels of a run that deliberately applies none, and
    refuse the product for lacking what it does not use."""
    src = _SRC["driver"]
    # ⚠️ THE SECOND `if mean3d:` -- the first selects the treatment token, this one is the
    # atmosphere/deck block. Splitting on the first occurrence silently tests the wrong
    # block and passes for the wrong reason.
    block = src.split("if mean3d:")[2].split("elif nlte:")[0]
    assert "dep = None" in block
    assert "if nlte:" in block and "assert_linelist_supports_nlte" in block
    # the withholding happens in the else, i.e. only for the LTE member
    assert block.index("if nlte:") < block.index("dep = None")


def test_the_treatment_tokens_come_from_the_axis_registry():
    """Never retyped — RYA-798 emitted a treatment `TREATMENTS` had never heard of and the
    product died at `build_product` after the synthesis had already run."""
    src = _SRC["driver"]
    assert "taxes.MEAN3D_NLTE_STAGGER" in src and "taxes.MEAN3D_LTE_STAGGER" in src
    assert "ENGINE-B-MEAN3D" not in src


# ── the atmosphere array is inert everywhere but the three fields iSpec reads ──

def test_the_3D_layer_array_survives_the_three_reads_iSpec_makes():
    """The array is NOT the atmosphere — babsma reads the real structure from the mul23
    file. iSpec touches it in exactly three places, and this reproduces all three against
    the real committed model."""
    from pipeline import mean3d_atmosphere as M
    layers, model = M.load(ROOT / "data/atmospheres/stagger_avg3d_rya442/sun_avg3d_stagger.mod")
    assert len(layers[0]) != 11                      # :114 -> MARCS-FILE=.false.
    assert np.isfinite(layers[:, 7]).all()           # :159 -> the departure tau
    assert np.isfinite(layers[0][-1])                # :85  -> radius / plane-parallel
    # everything else is NaN, so a fourth read fails loudly instead of using a plausible 0
    assert np.isnan(layers[:, [0, 1, 2, 3, 4, 5, 6, 8]]).all()
    assert len(layers) == model["ndep"] == 101
