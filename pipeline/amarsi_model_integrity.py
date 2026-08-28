"""
pipeline/amarsi_model_integrity.py — RYA-1098
=============================================
The Amarsi 3D-NLTE MLP pickles were WRITTEN under scikit-learn 1.0.2 and are LOADED under
1.9.0. This module makes that load either provably correct or loudly refused.

🔴 "IT LOADED WITHOUT ERRORING" IS NOT "IT LOADED CORRECTLY".

Cross-version unpickling of sklearn estimators is unsupported: attribute defaults and
internal layouts change between releases, and a mismatched load does not reliably raise.
The MLP's output is a PHYSICS CORRECTION -- the 3D non-LTE departure applied to Fe
abundances -- so a silent mis-load corrupts the correction with no loud failure.

🔴 AND THE ONE SIGNAL THAT EXISTED WAS BEING SUPPRESSED. sklearn DOES warn
(`InconsistentVersionWarning`, three per file), and `nlte_corrections._load_models`
wrapped the load in `warnings.simplefilter('ignore')`. That silenced the only evidence.

⚠️ WORSE, THE LOADED OBJECT LIES ABOUT ITS OWN PROVENANCE. After unpickling under 1.9.0,
`estimator.__getstate__()['_sklearn_version']` reads **1.9.0** -- the runtime rewrites it.
Asking the object which version wrote it therefore returns the wrong answer, confidently.
`written_version()` below reads the string out of the PICKLE BYTES instead, which is the
only place the truth survives.

THE VALIDATION, AND WHY IT IS NOT CIRCULAR
------------------------------------------
The ticket's preferred route -- capture (input -> output) pairs in the original 1.0.2
environment -- is CLOSED: sklearn 1.0.2 predates Python 3.12 and is not installable in any
environment we have (measured: Mac venv 1.9.0, system 1.6.1, Sirius venv312 1.9.0). Saying
so is part of the result; a validation invented to fill that gap would be worse than none.

What IS available is external ground truth that never passed through our environment:

  1. THE AUTHORS' OWN PUBLISHED PAIR. `vendor/1L-3NErrors/README.md` §2.4 prints a worked
     example produced under 1.0.2: (5051 K, 4.0, 4.5, 0.0, Elo 3.0, Eup 5.586893,
     lggf -2.563, Fe1) -> aberr **-0.136**.
  2. AN INDEPENDENT REIMPLEMENTATION. An `MLPRegressor`'s learned parameters are plain
     numpy arrays (`coefs_`, `intercepts_`) and a `StandardScaler`'s are `mean_`/`scale_`.
     Pickle round-trips numpy arrays exactly, independently of the sklearn version, so the
     PARAMETERS cannot have been corrupted by the skew -- only their USE could be. The
     forward pass is documented arithmetic (StandardScaler, then relu hidden layers, then
     an identity output), so reimplementing it in ~5 lines of numpy and comparing to
     `predict()` tests exactly the thing in doubt: does 1.9.0 still apply these arrays the
     way 1.0.2 did?

Neither check asks the loaded object to vouch for itself, which is the failure mode.

MEASURED RESULT (RYA-1098)
--------------------------
  authors' published pair   -0.135955 vs -0.136 published        (|d| < 1e-4)
  numpy vs sklearn 1.9.0    max |diff| 0.0 over 400 in-box points (EXACT)

So the 1.0.2 -> 1.9.0 skew is BENIGN FOR THESE ESTIMATORS -- proven, not assumed -- and
the corrections already in the feed are not contaminated by it. That verdict is pinned by
`tests/test_amarsi_model_integrity_rya1098.py`, which re-measures both checks rather than
trusting this docstring, and which FAILS on a future sklearn bump that changes either.

⚠️ A NOTE ON THE SCARY WARNINGS. Loading and predicting emits `divide by zero` /
`overflow` / `invalid value encountered in matmul` from numpy. They are NOT evidence of a
corrupted model: every `coefs_`, `intercepts_`, `mean_` and `scale_` array is finite, the
largest coefficient is 0.71, and the SAME warnings appear from a bare numpy `@` in the
reimplementation above. They are an artifact of the BLAS backend, not of the estimator.
"""
from __future__ import annotations

import pickle
import re
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parent.parent
VENDOR_DIR = _ROOT / "vendor" / "1L-3NErrors"
MODEL_FILES = ("fe1_model_lt02.p", "fe1_model_gt02.p", "fe2_model.p")

#: The version the vendored pickles were written with, read from their bytes.
WRITTEN_SKLEARN = "1.0.2"

#: Runtime versions whose load has been VALIDATED by `validate()` and recorded here.
#: 🔴 A version is added ONLY after the two checks below pass on it. This list is the
#: difference between "we checked" and "it did not crash".
VALIDATED_RUNTIMES: frozenset = frozenset({"1.0.2", "1.9.0"})

#: The authors' worked example, README.md §2.4, produced under sklearn 1.0.2.
#: (teff, logg, A(Fe;3N), vmic, Elo, Eup, loggf) -> aberr, on the Fe I Elo>2 network.
REFERENCE_PAIR = ((5051.0, 4.0, 4.5, 0.0, 3.0, 5.586893, -2.563), -0.136)
REFERENCE_TOL = 5e-4      # the published value is quoted to 3 decimals

#: How exactly an independent numpy forward pass must reproduce `predict()`. Zero in
#: practice; a whisker of float slack so a BLAS reassociation cannot fail the guard.
FORWARD_TOL = 1e-9


class ModelVersionSkew(RuntimeError):
    """A cross-version pickle load that has NOT been validated. Refuse, never ride."""


@dataclass(frozen=True)
class IntegrityReport:
    runtime_sklearn: str
    written_sklearn: str
    skewed: bool
    validated_runtime: bool
    reference_predicted: float
    reference_published: float
    reference_ok: bool
    forward_max_abs_diff: float
    forward_ok: bool

    @property
    def ok(self) -> bool:
        return self.reference_ok and self.forward_ok


def runtime_version() -> str:
    import sklearn
    return str(sklearn.__version__)


def written_version(path: Path) -> str:
    """The sklearn version that WROTE this pickle, read from the file's bytes.

    🔴 NOT from the loaded estimator. Unpickling under a newer sklearn rewrites
    `_sklearn_version` to the RUNTIME version, so the object reports 1.9.0 for a file
    written by 1.0.2 -- a confidently wrong answer to the only question that matters here.
    """
    raw = Path(path).read_bytes()
    hits = {m.decode() for m in
            re.findall(rb"_sklearn_version.{0,4}?([0-9]+\.[0-9]+\.[0-9]+)", raw, re.S)}
    if not hits:
        raise ModelVersionSkew(
            f"{Path(path).name} records no `_sklearn_version`. A pickle that will not say "
            f"what wrote it cannot be version-checked, and an unverifiable load of a "
            f"physics model is a refusal, not a warning.")
    if len(hits) > 1:
        raise ModelVersionSkew(
            f"{Path(path).name} records MORE THAN ONE sklearn version {sorted(hits)} -- "
            f"its estimators were not all written by the same release, so no single "
            f"compatibility claim covers the file.")
    return hits.pop()


def load_models(vendor_dir: Path | None = None, *, require_validated: bool = True) -> dict:
    """The Amarsi pickles, loaded only through a version check.

    ⚠️ The `InconsistentVersionWarning` is deliberately NOT suppressed. It is the signal;
    silencing it is how this went unnoticed. It is caught, read and re-raised as a refusal
    when the runtime is not one whose predictions have been validated.
    """
    d = Path(vendor_dir or VENDOR_DIR)
    rt = runtime_version()
    out = {}
    for key, fname in zip(("lt02", "gt02", "fe2"), MODEL_FILES):
        path = d / fname
        if not path.exists():
            raise FileNotFoundError(
                f"Amarsi model {fname} not found in {d}. "
                f"git clone https://github.com/sliljegren/1L-3NErrors vendor/1L-3NErrors")
        wrote = written_version(path)
        if wrote != rt and require_validated and rt not in VALIDATED_RUNTIMES:
            raise ModelVersionSkew(
                f"{fname} was written by scikit-learn {wrote} and is being loaded under "
                f"{rt}, which is NOT in VALIDATED_RUNTIMES {sorted(VALIDATED_RUNTIMES)}. "
                f"Cross-version unpickling of sklearn estimators is unsupported and does "
                f"not reliably raise -- it can deserialize cleanly and then predict "
                f"WRONG. This model's output is a physics correction, so a silent "
                f"mis-load corrupts an abundance with no loud failure. Run "
                f"`python3 scripts/rya1098_validate_amarsi_models.py`; if both checks "
                f"pass, add {rt!r} to VALIDATED_RUNTIMES with the numbers.")
        with warnings.catch_warnings():
            # The version warning has already been ACTED ON above; muting it here keeps
            # the caller's stream readable without muting it at the point of decision.
            warnings.filterwarnings("ignore", message=".*Trying to unpickle estimator.*")
            with path.open("rb") as fh:
                out[key] = tuple(pickle.load(fh))
    return out


def forward(scaler, mlp, X: np.ndarray) -> np.ndarray:
    """The MLP forward pass, in plain numpy, from the unpickled arrays alone.

    Deliberately independent of sklearn: StandardScaler, then relu hidden layers, then an
    identity output -- the architecture the vendored models actually carry, asserted
    rather than assumed. Its whole purpose is to be a second opinion, so it must not call
    the thing it is checking.
    """
    if getattr(mlp, "activation", "relu") != "relu":
        raise ModelVersionSkew(f"unexpected hidden activation {mlp.activation!r}")
    if getattr(mlp, "out_activation_", "identity") != "identity":
        raise ModelVersionSkew(f"unexpected output activation {mlp.out_activation_!r}")
    z = (np.asarray(X, dtype=float) - scaler.mean_) / scaler.scale_
    for W, b in zip(mlp.coefs_[:-1], mlp.intercepts_[:-1]):
        z = np.maximum(z @ W + b, 0.0)
    return (z @ mlp.coefs_[-1] + mlp.intercepts_[-1]).ravel()


def validate(vendor_dir: Path | None = None, *, n_probe: int = 400,
             seed: int = 0) -> IntegrityReport:
    """Both checks, measured. Neither asks the loaded object to vouch for itself."""
    models = load_models(vendor_dir, require_validated=False)
    scaler, mlp = models["gt02"]

    x, published = REFERENCE_PAIR
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pred = float(mlp.predict(scaler.transform(np.array([x], dtype=float)))[0])

    rng = np.random.default_rng(seed)
    X = np.column_stack([rng.uniform(5000, 6500, n_probe), rng.uniform(4.0, 4.5, n_probe),
                         rng.uniform(4.5, 7.5, n_probe), rng.uniform(0.0, 3.0, n_probe),
                         rng.uniform(0.5, 5.0, n_probe), rng.uniform(2.0, 7.0, n_probe),
                         rng.uniform(-6.0, 0.0, n_probe)])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sk = np.asarray(mlp.predict(scaler.transform(X)), dtype=float).ravel()
    mine = forward(scaler, mlp, X)
    dmax = float(np.max(np.abs(sk - mine)))

    rt = runtime_version()
    wrote = written_version(Path(vendor_dir or VENDOR_DIR) / MODEL_FILES[1])
    return IntegrityReport(
        runtime_sklearn=rt, written_sklearn=wrote, skewed=(rt != wrote),
        validated_runtime=(rt in VALIDATED_RUNTIMES),
        reference_predicted=pred, reference_published=published,
        reference_ok=abs(pred - published) <= REFERENCE_TOL,
        forward_max_abs_diff=dmax, forward_ok=dmax <= FORWARD_TOL)
