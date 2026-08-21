"""RYA-767 Seam 2: a run is an object, and its preconditions are DATA.

Everything pinned here is something an operator had to know by heart during
RYA-933/934, where getting it wrong produced a plausible wrong answer rather
than an error. The value of this layer is not that it runs anything -- it runs
nothing -- but that the knowledge stops living in a person.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))


@pytest.fixture(scope="module", autouse=True)
def _kp(tmp_path_factory):
    kp = tmp_path_factory.mktemp("kp")
    (kp / "lm0296").touch()
    os.environ.setdefault("CODEX_KP_ATLAS", str(kp))


def _d(**kw):
    from pipeline.run_descriptor import RunDescriptor
    base = dict(element="Fe", ion="I", instrument="harps", holding="solar_harps",
                lo_A=3800.0, hi_A=6910.0)
    base.update(kw)
    return RunDescriptor(**base)


def _resolve(descriptor, **kw):
    from pipeline.run_descriptor import resolve
    opts = dict(interpreter="/x/venv312/bin/python", ispec_dir="/x/ispec")
    opts.update(kw)
    return resolve(descriptor, **opts)


def test_the_key_carries_the_holding_not_just_the_instrument():
    """Two products of one instrument differing by telluric state must not collide."""
    a = _d(holding="solar_harps").key
    b = _d(holding="solar_harps_molecfit_corrected").key
    assert a != b
    assert "solar_harps_molecfit_corrected" in b


def test_band_and_method_come_from_policy_not_from_the_caller():
    """near-UV forbids profile-fit; that is not a preference anyone gets to express."""
    assert _d(lo_A=3800.0, hi_A=6910.0).band == "VIS"
    assert _resolve(_d()).method == "profile-fit"
    kp = _d(instrument="kpno_solar_atlas", holding="solar_kpno", lo_A=3000.0, hi_A=3800.0)
    assert kp.band == "near-UV"
    assert _resolve(kp).method == "synthesis"


def test_a_holding_that_does_not_cover_the_band_is_blocked_before_anything_runs():
    kp = _d(instrument="kpno_solar_atlas", holding="solar_kpno_kurucz2005_corrected",
            lo_A=10000.0, hi_A=24000.0)
    r = _resolve(kp)
    assert not r.runnable
    assert "does not cover" in r.blocked_reason


def test_not_holding_a_window_is_reported_as_such_not_as_a_failure():
    """'We do not hold this' and 'this failed' are different answers (RYA-796/833)."""
    kp = _d(instrument="kpno_solar_atlas", holding="solar_kpno_kurucz2005_corrected",
            lo_A=10000.0, hi_A=24000.0)
    assert "NOT a failure" in _resolve(kp).blocked_reason


def test_an_unwired_holding_names_the_rya904_shape():
    r = _resolve(_d(holding="solar_harps_imaginary"))
    assert not r.runnable
    assert "not wired" in r.blocked_reason
    assert "like having no data" in r.blocked_reason


def test_the_numpy_ceiling_is_a_precondition_not_a_discovery():
    """RYA-682: above the ceiling iSpec writes a zero-row artifact and EXITS 0."""
    r = _resolve(_d(), interpreter=None)
    check = next(p for p in r.preconditions if p.name == "numpy_below_ceiling")
    assert not check.satisfied
    # The message must name the FAILURE MODE, not merely the version bound:
    # the danger is a silent empty product, not a crash.
    assert "silent empty product" in check.detail.lower()
    assert not r.runnable


def test_the_engine_environment_is_a_precondition():
    r = _resolve(_d(), ispec_dir=None)
    check = next(p for p in r.preconditions if p.name == "ispec_dir_set")
    assert not check.satisfied
    assert not r.runnable


def test_the_ew_step_is_ordered_first_and_says_why():
    r = _resolve(_d())
    names = [s["name"] for s in r.steps]
    assert names == ["measure_ew", "derive_products"]
    assert "does not measure" in r.steps[0]["why_ordered_first"]


def test_synthesis_bands_have_no_separate_ew_step():
    kp = _d(instrument="kpno_solar_atlas", holding="solar_kpno", lo_A=3000.0, hi_A=3800.0)
    assert [s["name"] for s in _resolve(kp).steps] == ["derive_products"]


def test_every_step_carries_a_postcondition_that_is_not_an_exit_code():
    r = _resolve(_d())
    for step in r.steps:
        assert step["postcondition"]
        assert "exit" not in step["postcondition"].lower()
    assert "zero-row product is a FAILURE" in r.steps[-1]["postcondition"]


def test_harps_declares_its_extent_rather_than_leaving_it_discovered():
    """Without a declared span, `covers()` answered True for the NIR."""
    import measure_band_ew as M
    for holding in ("solar_harps", "solar_harps_molecfit_corrected"):
        spec = next(h for h in M._INSTRUMENT_HOLDINGS["harps"]
                    if h.holding_id == holding)
        assert spec.span_A is not None, f"{holding} must declare its extent"
        assert not spec.covers(11000.0, 2.0), f"{holding} must not claim the NIR"
        assert spec.covers(5000.0, 2.0)


def test_the_resolver_decides_nothing_about_science():
    """It resolves availability. It must not import, or become, a science engine.

    Checked on the IMPORTS rather than on words in prose -- the docstring
    legitimately mentions `ispec/abundances.py` when explaining the numpy
    ceiling, and a substring match on "abundance" would fail on an accurate
    comment while missing an actual dependency.
    """
    import ast as _ast
    tree = _ast.parse((ROOT / "pipeline" / "run_descriptor.py").read_text())
    imported = set()
    for node in _ast.walk(tree):
        if isinstance(node, _ast.Import):
            imported.update(a.name for a in node.names)
        elif isinstance(node, _ast.ImportFrom) and node.module:
            imported.add(node.module)
    for engine in ("ispec", "pipeline.abundances_derive", "pipeline.cno_synthesis",
                   "pipeline.band_products", "pipeline.lines_fit"):
        assert engine not in imported, (
            f"run_descriptor imports {engine} -- orchestration must resolve what is "
            f"possible, never compute the science itself")
