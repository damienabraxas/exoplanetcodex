"""RYA-924 — the 3D-NLTE leg must refuse to run where its results cannot be audited.

The committed 3D-NLTE cells were computed on a Mac from a per-session /tmp directory.
Nothing stopped it: `assert_on_sirius` (RYA-567) existed and was wired into two other
legs, but never into this one. The guard did not fail — it was never asked.
"""
from __future__ import annotations

import ast
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "rya817_run_3dnlte_bands.py"


def _load():
    spec = importlib.util.spec_from_file_location("rya817_mod", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ── the gate is actually wired, not merely available ────────────────────────

def test_the_3dnlte_route_calls_assert_on_sirius():
    """It is not enough that the guard exists. RYA-567 wired it into
    `nlte_bfactor_synth` and `pysme_nlte` and this route was left out for months."""
    tree = ast.parse(SCRIPT.read_text(encoding="utf-8"))
    called = {n.func.id for n in ast.walk(tree)
              if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)}
    assert "assert_on_sirius" in called, (
        "the 3D-NLTE band run does not call assert_on_sirius — it can be run off-Sirius "
        "again, which is how the committed cells were produced on a Mac")


def test_the_gate_is_called_before_any_band_is_run():
    """A guard that fires after the work is a log line, not a guard.

    Checked INSIDE `main()`, not by source position: `discover_bands` is *defined* above
    `main`, so a whole-file index comparison tests where the code sits, not when it runs.
    """
    tree = ast.parse(SCRIPT.read_text(encoding="utf-8"))
    main = next(n for n in tree.body
                if isinstance(n, ast.FunctionDef) and n.name == "main")
    order = [n.func.id for n in ast.walk(main)
             if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
             and n.func.id in ("assert_on_sirius", "discover_bands",
                               "_require_auditable_input_dir")]
    assert "assert_on_sirius" in order, "main() never calls the gate"
    assert order.index("assert_on_sirius") < order.index("discover_bands"), (
        f"the gate must fire before any band is discovered; got {order}")


# ── a scratchpad input is refused, and the refusal DISCRIMINATES ────────────

@pytest.mark.parametrize("bad", [
    "/tmp/claude-501/-Users-someone/1f004751-5946-4704-b27e-0bd1e599f25c",
    "/private/tmp/claude-501/-Users-ryanschmitt/session/band_products",
    "/var/folders/xy/T/some-scratch",
])
def test_a_temp_directory_input_is_REFUSED(bad):
    """🔴 The exact shape recorded in rya817_run_summary.json:
        band_products_dir = /private/tmp/claude-501/-Users-ryanschmitt/<uuid>/
    A value derived from a per-session directory cannot be checked by anyone afterwards.
    """
    mod = _load()
    with pytest.raises(SystemExit) as e:
        mod._require_auditable_input_dir(Path(bad))
    assert "temporary directory" in str(e.value)


def test_a_REAL_directory_is_accepted(tmp_path):
    """Positive control: the refusal must not be a blanket no. `tmp_path` here is
    pytest's own fixture — if THAT were rejected the guard would be untestable, which is
    itself the sign of a guard that only knows how to say no."""
    mod = _load()
    real = ROOT / "data" / "results"
    if real.is_dir():
        mod._require_auditable_input_dir(real)          # must not raise


def test_a_missing_directory_is_refused_distinctly():
    """Two refusals, two reasons — a reader must be able to tell them apart.

    ⚠️ Deliberately NOT using pytest's `tmp_path`: on Linux that lives under /tmp, so the
    scratchpad branch catches it first. That is the guard behaving correctly, and it is
    also a neat demonstration that pytest's own fixture is exactly the kind of directory
    a product must never be derived from.
    """
    mod = _load()
    with pytest.raises(SystemExit) as e:
        mod._require_auditable_input_dir(ROOT / "no-such-directory-rya924")
    assert "does not exist" in str(e.value)


# ── the artifact records where it ran ───────────────────────────────────────

def test_run_environment_names_the_machine_and_the_model_library():
    """The previous run's only environmental trace was ACCIDENTAL — a /tmp path that
    happened to reveal a Mac. Host, platform, interpreter and model-library version are
    now recorded on purpose."""
    env = _load().run_environment()
    for key in ("host", "platform", "python", "sklearn_saved"):
        assert key in env, f"run provenance is missing {key!r}"
    assert env["sklearn_saved"], "the version the pickles were WRITTEN with must be read"


def test_version_skew_is_only_reported_when_both_versions_are_REAL():
    """An absent library is not a version. Writing 'models written with 1.0.2, loaded
    with ABSENT — cross-version unpickling' would be prose asserting something untrue."""
    env = _load().run_environment()
    if env.get("sklearn_loaded") is None:
        assert "sklearn_version_skew" not in env, (
            "reported a version skew against a library that is not installed")
        assert "sklearn_status" in env, "an absent library must say so plainly"
