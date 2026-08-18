"""
tests/test_synth_ew_rya878.py — RYA-878
=======================================
The SynthesisHandler control's ANGLE 1 compares a SYNTHETIC EW to a MEASURED one. It was
never a measured angle, for a structural reason: the handler computed a synthetic EW
inline and the production synth-v2 path computed none at all, so there was nothing on the
engine side to compare against and "the engine shares the offset" could only be inferred.

Fixing that means TWO callers of one quantity, which is the moment a second definition
appears and the two sides quietly stop meaning the same thing (RYA-845, and again in
855/869/873). So the definition lives in `pipeline.synth_ew` and these tests hold it there.

They pin structure and invariants, not the offset — the offset is a property of today's
line set and is measured in `data/results/rya878/`.
"""
from __future__ import annotations

import ast
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import synth_ew                                      # noqa: E402


# ── one definition ───────────────────────────────────────────────────────────────────

def _builds_the_difference(path: Path) -> bool:
    """Does this file CONSTRUCT the difference-of-syntheses, as code?

    Checked on the AST: `scripts/rya873_control_rca.py` quotes the expression in its
    docstring to explain the mechanism, and a text search cannot tell a citation from a
    copy — the same trap that bit the RYA-873 ownership guard.
    """
    tree = ast.parse(path.read_text())
    for n in ast.walk(tree):
        if (isinstance(n, ast.BinOp) and isinstance(n.op, ast.Sub)
                and isinstance(n.left, ast.Name) and isinstance(n.right, ast.Name)
                and {n.left.id, n.right.id} == {"without", "with_el"}):
            return True
    return False


def test_only_synth_ew_builds_the_difference_of_syntheses():
    """The EW is `trapz(clip(without - with_el))`. Exactly one module may construct it;
    a second copy is how the handler and the deriver would come to disagree about what
    'the synthetic EW' means while each stayed internally consistent."""
    hits = [f.relative_to(ROOT)
            for f in sorted((ROOT / "pipeline").rglob("*.py"))
            + sorted((ROOT / "scripts").rglob("*.py"))
            if _builds_the_difference(f)]
    assert hits == [Path("pipeline/synth_ew.py")], hits


@pytest.mark.parametrize("mod", ["pipeline/measure/synthesis.py",
                                 "pipeline/abundances_derive.py"])
def test_both_sides_call_the_shared_definition(mod):
    """The handler and the production path must both go through it — otherwise ANGLE 1
    is back to comparing two things measured differently, which is the defect."""
    src = (ROOT / mod).read_text()
    assert "synth_ew.synthetic_ew_mA(" in src, f"{mod} does not call the shared EW"


def test_the_production_path_banks_the_column():
    """A quantity with nowhere to live is a quantity nobody can check (RYA-843). The
    engine side of ANGLE 1 exists only if it reaches the artifact."""
    src = (ROOT / "pipeline" / "abundances_derive.py").read_text()
    assert "'ew_synth_mA'" in src or '"ew_synth_mA"' in src


def test_the_engine_computes_no_ew_for_a_line_that_produced_no_abundance():
    """An EW evaluated at a failed or edge-pinned fit is the EW of a number the fit did
    not determine — it would look like a measurement and be one of the search's failures.
    """
    src = (ROOT / "pipeline" / "abundances_derive.py").read_text()
    i = src.index("ew_synth_mA = float('nan')")
    guard = src[i:i + 260]
    assert "if status == 'ok':" in guard, guard[:160]


# ── invariants of the quantity itself ────────────────────────────────────────────────

def _flat(w_nm, trial_A, **kw):
    """A stand-in generator: a Gaussian whose depth scales with 10**A, continuum 1."""
    w = np.asarray(w_nm, float) * 10.0
    depth = min(0.9, 10.0 ** (trial_A - 7.5) * 0.5)
    return 1.0 - depth * np.exp(-0.5 * ((w - w.mean()) / 0.05) ** 2)


def test_depleting_the_element_leaves_no_absorption_to_integrate():
    """The 'without' spectrum is the baseline the difference is taken against; if the
    depletion did not remove the element the EW would be a difference of two similar
    lines and systematically too small."""
    assert synth_ew.EW_DEPLETION_DEX >= 3.0
    deep = _flat(np.linspace(499.0, 501.0, 50), 7.5)
    gone = _flat(np.linspace(499.0, 501.0, 50), 7.5 - synth_ew.EW_DEPLETION_DEX)
    assert (1.0 - gone).max() < 0.01 * (1.0 - deep).max()


def test_ew_is_positive_and_grows_with_abundance():
    kw = dict(centre_A=5000.0, fit_half_width_A=0.2, synth_kwargs={})
    lo = synth_ew.synthetic_ew_mA(_flat, abundance=7.0, **kw)
    hi = synth_ew.synthetic_ew_mA(_flat, abundance=7.5, **kw)
    assert 0 < lo < hi, (lo, hi)


def test_a_failing_generator_yields_nan_not_zero():
    """Zero is a measurement; NaN is 'no measurement'. A synthesis that raised must not
    report a line as having no absorption."""
    def boom(*a, **k):
        raise RuntimeError("iSpec fell over")
    v = synth_ew.synthetic_ew_mA(boom, centre_A=5000.0, abundance=7.5,
                                 fit_half_width_A=0.2, synth_kwargs={})
    assert np.isnan(v)


# ── the integration range ────────────────────────────────────────────────────────────

def test_the_range_is_wider_than_the_fit_window_and_floored():
    """Integrating over the fit window truncates the damping wings and drove the EW
    ratio to 0.250 by construction — the reason the range is deliberately wider."""
    assert synth_ew.ew_half_width_A(0.05) == synth_ew.EW_HALF_WIDTH_FLOOR_A
    assert synth_ew.ew_half_width_A(1.0) == synth_ew.EW_HALF_WIDTH_FACTOR * 1.0
    assert synth_ew.EW_HALF_WIDTH_FACTOR > 1.0


def test_production_never_overrides_the_range():
    """`integration_half_width_A` exists for RYA-878's definitional sweep. If a
    production caller started passing it, the banked EW would stop meaning what the
    definition says and the sweep's own baseline would move under it."""
    for mod in ("pipeline/measure/synthesis.py", "pipeline/abundances_derive.py"):
        assert "integration_half_width_A" not in (ROOT / mod).read_text(), mod
    sweep = (ROOT / "scripts" / "rya878_ew_definition_sweep.py").read_text()
    assert "integration_half_width_A=hw" in sweep, "the sweep must be the one that sweeps"


def test_the_override_actually_changes_the_range():
    """VERIFY THE KNOB DISCRIMINATES — a sweep whose parameter did nothing would report a
    flat ratio and be read as 'the range is not the cause'."""
    kw = dict(centre_A=5000.0, abundance=7.5, fit_half_width_A=1.0, synth_kwargs={})
    wide = synth_ew.synthetic_ew_mA(_flat, **kw)
    narrow = synth_ew.synthetic_ew_mA(_flat, integration_half_width_A=0.1, **kw)
    assert narrow < wide, (narrow, wide)


def test_the_signature_defaults_to_the_derived_range():
    import inspect
    p = inspect.signature(synth_ew.synthetic_ew_mA).parameters
    assert p["integration_half_width_A"].default is None
