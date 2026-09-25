"""
RYA-1220 -- the PySME invariant the native NLTE results depend on, and a way to check it.

🔴 THE FIX FOR THE ORIGINAL DEFECT LIVES OUTSIDE THIS REPOSITORY. SMElib stores departure
coefficients GLOBALLY in the C library. If an LTE synthesis runs after an NLTE one without
resetting them, it silently reuses the previous NLTE coefficients, and the "LTE" leg is not
LTE. Every NLTE correction is A(NLTE) - A(LTE), so a contaminated LTE leg drives the
correction toward zero -- which is exactly what a near-zero correction looks like.

`git grep ResetDepartureCoefficients` over this repo finds only PROSE inside two JSON
artifacts. The actual guarantee is a property of the PySME install on the synthesis host:

    pysme.nlte.Nlte.update_coefficients must call dll.ResetDepartureCoefficients()
    BEFORE its "no NLTE to do" early return.

Placement is the whole point. A reset that sits after the early return never runs on the
LTE path, which is the path that needs it.

⚠️ AND THERE IS MORE THAN ONE COPY ON SIRIUS. /srv/codex/engines/pysme_src and the
installed /mnt/codex-data/venv_pysme tree carry DIFFERENT code. Recording which one runs
is not bookkeeping: the results depend on it.

This module does not patch anything. It states the invariant and checks it, so a
reinstall that quietly reverts the behaviour fails a test instead of shifting the science.
"""
from __future__ import annotations

import inspect
import pathlib

#: Fingerprint of the install that produced the RYA-1220 native N I corrections.
RECORDED_ENVIRONMENT = {
    "host": "sirius",
    "interpreter": "/mnt/codex-data/venv_pysme/bin/python",
    "module_path": "/mnt/codex-data/venv_pysme/lib/python3.12/site-packages/pysme/nlte.py",
    "pysme_version": "1.0.2",
    "nlte_py_sha256_prefix": "31e43a179bf8f8cf",
    "note": ("A SECOND, DIFFERENT copy of pysme exists at /srv/codex/engines/pysme_src; "
             "it is not the one that runs. Check module_path before trusting a result."),
}

RESET_CALL = "ResetDepartureCoefficients"
EARLY_RETURN_MARKER = "No NLTE to do"


class NlteEnvironmentError(RuntimeError):
    """The install cannot be trusted to produce a genuine LTE leg."""


def _carrier():
    """The class in pysme.nlte that owns update_coefficients.

    Discovered, not hard-coded: the class is spelled `NLTE` in the installed 1.0.2 and a
    guessed `Nlte` imports cleanly as a NAME ERROR that a broad except would swallow into
    "pysme is not installed" -- the check would then pass by being silently skipped on the
    one host where it matters.
    """
    import pysme.nlte as module
    for obj in vars(module).values():
        if inspect.isclass(obj) and hasattr(obj, "update_coefficients"):
            return module, obj
    raise NlteEnvironmentError(
        "pysme.nlte exposes no class with update_coefficients; the invariant cannot be "
        "checked and must not be assumed")


def update_coefficients_source() -> str | None:
    """Source of pysme's update_coefficients, or None when pysme is not IMPORTABLE.

    Only a genuine ImportError returns None. Anything else raises, so a renamed symbol
    cannot masquerade as an absent install.
    """
    try:
        import pysme.nlte  # noqa: F401
    except ImportError:
        return None
    _, carrier = _carrier()
    return inspect.getsource(carrier.update_coefficients)


def reset_precedes_early_return(source: str) -> bool:
    """Is the departure-coefficient reset reached on the LTE path?"""
    if RESET_CALL not in source:
        return False
    if EARLY_RETURN_MARKER not in source:
        # No early return at all: an unconditional reset is still correct.
        return True
    return source.index(RESET_CALL) < source.index(EARLY_RETURN_MARKER)


def verify() -> dict:
    """Report whether the running install satisfies the invariant."""
    source = update_coefficients_source()
    if source is None:
        return {"state": "UNAVAILABLE", "ok": None,
                "why": "pysme is not importable here; run this on the synthesis host"}
    ok = reset_precedes_early_return(source)
    module = pathlib.Path(inspect.getfile(_carrier()[0]))
    return {"state": "CHECKED", "ok": ok, "module_path": str(module),
            "why": ("reset precedes the no-NLTE early return, so the LTE leg is genuinely "
                    "LTE" if ok else
                    "THE RESET DOES NOT RUN ON THE LTE PATH. LTE syntheses will reuse the "
                    "previous NLTE departure coefficients and every NLTE correction "
                    "derived here is suspect.")}


def require() -> None:
    got = verify()
    if got["state"] == "CHECKED" and not got["ok"]:
        raise NlteEnvironmentError(got["why"])
