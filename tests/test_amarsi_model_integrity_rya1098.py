"""
RYA-1098 — the Amarsi MLP load must be validated, not merely uneventful.

🔴 THE CLASS OF DEFECT. The vendored pickles were written by scikit-learn 1.0.2 and are
loaded under 1.9.0. Cross-version unpickling of sklearn estimators is unsupported and does
NOT reliably raise: it can deserialize cleanly and then predict subtly WRONG. These models
emit a physics correction applied to Fe abundances, so a silent mis-load is an abundance
error with no stack trace. `nlte_corrections._load_models` was wrapping the load in
`warnings.simplefilter('ignore')` -- silencing the one signal that existed.

Every check here re-MEASURES; none of them trusts a docstring, and none asks the loaded
object to vouch for itself. That last point is not pedantry: unpickling under a newer
sklearn REWRITES `_sklearn_version` on the estimator, so the object reports 1.9.0 for a
file written by 1.0.2 -- a confidently wrong answer to the only question that matters.
"""
from __future__ import annotations

import ast
import warnings
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sklearn = pytest.importorskip("sklearn", reason="the Amarsi engine needs scikit-learn")

from pipeline import amarsi_model_integrity as ami          # noqa: E402


@pytest.fixture(scope="module")
def report():
    if not (ami.VENDOR_DIR / ami.MODEL_FILES[0]).exists():
        pytest.skip("vendor/1L-3NErrors not present")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return ami.validate()


# ── the version fact, read from the right place ──────────────────────────────────

def test_the_written_version_is_read_from_the_BYTES_not_the_object():
    """The loaded estimator reports the RUNTIME version, so asking it is worse than not
    asking. Pinned by measurement: the file says 1.0.2 and the object says otherwise."""
    if not (ami.VENDOR_DIR / ami.MODEL_FILES[0]).exists():
        pytest.skip("vendor absent")
    for f in ami.MODEL_FILES:
        assert ami.written_version(ami.VENDOR_DIR / f) == ami.WRITTEN_SKLEARN

    import pickle
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        with (ami.VENDOR_DIR / ami.MODEL_FILES[1]).open("rb") as fh:
            scaler, mlp = pickle.load(fh)
    reported = mlp.__getstate__().get("_sklearn_version")
    assert reported == sklearn.__version__, (
        "the premise of `written_version` is that the unpickled object reports the "
        "RUNTIME version. If that stops being true the function's reason to exist has "
        "changed and its docstring is now wrong.")
    if sklearn.__version__ != ami.WRITTEN_SKLEARN:
        assert reported != ami.WRITTEN_SKLEARN, "the object must NOT know it was 1.0.2"


# ── the two validations, re-measured ─────────────────────────────────────────────

def test_the_authors_published_pair_reproduces(report):
    """CHECK 1 — ground truth that never passed through our environment. sklearn 1.0.2
    predates Python 3.12 and is not installable here, so capturing our own 1.0.2 pairs is
    impossible; the authors' worked example (README §2.4) is the substitute, and it is a
    better one because it is independent of us entirely."""
    assert report.reference_ok, (
        f"predicted {report.reference_predicted:+.6f} vs the authors' published "
        f"{report.reference_published:+.3f} -- beyond {ami.REFERENCE_TOL}. The load is "
        f"NOT reproducing the 1.0.2 environment and every correction is suspect.")


def test_an_independent_forward_pass_reproduces_predict(report):
    """CHECK 2 — the learned parameters are plain numpy arrays and pickle round-trips
    those exactly whatever the sklearn version, so the parameters cannot have been
    corrupted; only their USE could be. This tests the use."""
    assert report.forward_ok, (
        f"a numpy reimplementation of the documented forward pass disagrees with "
        f"predict() by {report.forward_max_abs_diff:.3e} -- sklearn {report.runtime_sklearn} "
        f"is not applying these arrays the way the architecture says it should.")


def test_the_reimplementation_is_actually_independent():
    """A second opinion that calls the thing it is checking is not a second opinion."""
    src = (ROOT / "pipeline" / "amarsi_model_integrity.py").read_text()
    fn = next(n for n in ast.walk(ast.parse(src))
              if isinstance(n, ast.FunctionDef) and n.name == "forward")
    body = ast.get_source_segment(src, fn)
    assert ".predict(" not in body, "forward() must not call predict()"
    assert ".transform(" not in body, "forward() must not call the scaler's transform()"


def test_the_forward_pass_REFUSES_an_architecture_it_does_not_implement():
    """POSITIVE CONTROL. `forward` hardcodes relu hidden layers and an identity output. If
    a future model used something else, silently applying relu anyway would produce a
    confident wrong second opinion -- worse than no check at all."""
    class _Fake:
        activation = "tanh"
        out_activation_ = "identity"
        coefs_ = [np.eye(2)]
        intercepts_ = [np.zeros(2)]

    class _Scaler:
        mean_ = np.zeros(2)
        scale_ = np.ones(2)

    with pytest.raises(ami.ModelVersionSkew, match="activation"):
        ami.forward(_Scaler(), _Fake(), np.zeros((1, 2)))


# ── the recurrence guard ─────────────────────────────────────────────────────────

def test_an_unvalidated_runtime_is_REFUSED_not_ridden():
    """🔴 THE GUARD. A future sklearn bump must fail loud rather than silently reintroduce
    a skewed load. Simulated by asking the loader for a runtime that is not in the
    validated set -- the same code path a real bump would take."""
    if not (ami.VENDOR_DIR / ami.MODEL_FILES[0]).exists():
        pytest.skip("vendor absent")
    real = ami.VALIDATED_RUNTIMES
    try:
        ami.VALIDATED_RUNTIMES = frozenset({"0.0.0"})
        with pytest.raises(ami.ModelVersionSkew, match="VALIDATED_RUNTIMES"):
            ami.load_models(require_validated=True)
    finally:
        ami.VALIDATED_RUNTIMES = real


def test_the_current_runtime_IS_in_the_validated_set(report):
    """The other half: the guard must be open for the runtime we actually use, or the
    engine is simply off and the guard is an outage rather than a check."""
    assert report.validated_runtime, (
        f"sklearn {report.runtime_sklearn} is not in VALIDATED_RUNTIMES. Run "
        f"scripts/rya1098_validate_amarsi_models.py and, if both checks pass, add it "
        f"WITH the numbers -- never as a convenience.")


def test_the_production_loader_no_longer_suppresses_the_version_warning():
    """The original defect, pinned. `_load_models` must route through the checked loader
    and must not blanket-ignore warnings around a pickle load.

    ⚠️ SCANNED AS CODE, NOT AS TEXT. A first version of this test matched substrings in
    the function's source and went red on the DOCSTRING, which quotes the very defect it
    describes ("the body here was `warnings.simplefilter('ignore')` wrapped around three
    `pickle.load` calls"). Prose about a defect is not the defect; the AST separates them
    and a grep cannot.
    """
    src = (ROOT / "pipeline" / "nlte_corrections.py").read_text()
    fn = next(n for n in ast.walk(ast.parse(src))
              if isinstance(n, ast.FunctionDef) and n.name == "_load_models")
    called = set()
    for node in ast.walk(fn):
        if isinstance(node, ast.Call):
            f = node.func
            if isinstance(f, ast.Attribute):
                base = getattr(f.value, "id", "")
                called.add(f"{base}.{f.attr}" if base else f.attr)
            elif isinstance(f, ast.Name):
                called.add(f.id)
    assert "simplefilter" not in {c.split(".")[-1] for c in called}, (
        f"a blanket warning filter around the model load is exactly the defect RYA-1098 "
        f"fixes -- InconsistentVersionWarning was the only evidence there was. Calls: "
        f"{sorted(called)}")
    assert "pickle.load" not in called, (
        f"the load belongs in one place, behind the check. Calls: {sorted(called)}")
    imported = {n.module for n in ast.walk(fn) if isinstance(n, ast.ImportFrom)}
    assert any("amarsi_model_integrity" in str(m) for m in imported), (
        f"the production load must go through the version check; imports: {imported}")


def test_every_model_array_is_finite():
    """The scary `divide by zero / overflow / invalid value encountered in matmul` seen
    around these models is a BLAS artifact, not a corrupted estimator -- and the way to
    know that rather than hope it is to check the arrays."""
    if not (ami.VENDOR_DIR / ami.MODEL_FILES[0]).exists():
        pytest.skip("vendor absent")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        models = ami.load_models(require_validated=False)
    for key, (scaler, mlp) in models.items():
        assert np.isfinite(scaler.mean_).all() and np.isfinite(scaler.scale_).all(), key
        assert all(np.isfinite(c).all() for c in mlp.coefs_), key
        assert all(np.isfinite(b).all() for b in mlp.intercepts_), key
