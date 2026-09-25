"""RYA-1220 -- the departure-coefficient reset must stay on the LTE path.

The repair that produced the native N I corrections is a property of the PySME install,
not of this repository, so nothing in CI could see it. These tests make the invariant
checkable: they verify the logic here always, and the live install where pysme exists.
"""
from __future__ import annotations

import pytest

from pipeline import nlte_environment as env


def test_reset_before_the_early_return_is_accepted():
    src = """
        dll.ResetDepartureCoefficients()
        if np.all(self.grids == ""):
            # No NLTE to do
            return sme
    """
    assert env.reset_precedes_early_return(src) is True


def test_reset_after_the_early_return_is_rejected():
    """This is the original defect: the LTE path returns before resetting."""
    src = """
        if np.all(self.grids == ""):
            # No NLTE to do
            return sme
        dll.ResetDepartureCoefficients()
    """
    assert env.reset_precedes_early_return(src) is False


def test_a_missing_reset_is_rejected():
    assert env.reset_precedes_early_return("if x:\n    # No NLTE to do\n    return sme") is False


def test_an_unconditional_reset_with_no_early_return_is_accepted():
    assert env.reset_precedes_early_return("dll.ResetDepartureCoefficients()") is True


def test_the_recorded_environment_names_the_ambiguity():
    rec = env.RECORDED_ENVIRONMENT
    assert rec["module_path"].startswith("/mnt/codex-data/venv_pysme")
    assert "pysme_src" in rec["note"], "the second, different copy must stay recorded"


@pytest.mark.skipif(env.update_coefficients_source() is None,
                    reason="pysme not installed here; this check belongs on the synthesis host")
def test_the_live_install_satisfies_the_invariant():
    got = env.verify()
    assert got["ok"] is True, got["why"]
