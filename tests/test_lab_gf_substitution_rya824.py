"""
tests/test_lab_gf_substitution_rya824.py — RYA-824
==================================================
The lab-gf sub-pool rests on one operation: replacing a single physical line's total
oscillator strength in the synth line list and changing nothing else. Three ways that
goes wrong quietly, one test each:

  * it moves an HFS component instead of the physical line, changing the BRANCHING as
    well as the total — a different physical claim from the one the lab paper makes;
  * it matches more than one physical line and silently picks the first, so the
    substituted line is not the line the census named;
  * it perturbs rows outside the target cluster, which would make the re-inverted
    abundance a function of the whole line list rather than of the one gf that changed.

The synthesis itself is not exercised here (it needs iSpec and Turbospectrum, which live
on Sirius). What is exercised is the input to it, which is where the reasoning lives.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from pipeline.gf_grades import lab_lines  # noqa: E402

pytest.importorskip("pandas")


def _fake_linelist(rows):
    """A minimal synth-array stand-in with the fields the substitution reads."""
    dt = np.dtype([("element", "U8"), ("ion", "i4"), ("molecule", "U8"),
                   ("wave_A", "f8"), ("lower_state_eV", "f8"), ("loggf", "f8")])
    return np.array(rows, dtype=dt)


def _sub(*args, **kw):
    from rya824_lab_gf_subpool import substitute_gf
    return substitute_gf(*args, **kw)


def test_the_lab_table_is_present_and_positive():
    """RYA-799's positivity guard, restated: this ticket re-reads the same table."""
    lab = lab_lines()
    assert len(lab) > 0
    assert (lab.e_loggf_dex > 0).all()
    assert lab.e_loggf_dex.max() <= 0.5


def test_substitution_moves_the_total_and_preserves_hfs_branching():
    """Two components of one physical line must both shift by the SAME amount, so their
    ratio is untouched while the total lands on the lab value."""
    arr = _fake_linelist([
        ("Fe", 1, "", 6000.000, 3.0, -1.30),
        ("Fe", 1, "", 6000.004, 3.0, -1.70),   # same physical line, HFS component
    ])
    before = arr["loggf"].copy()
    ratio_before = before[0] - before[1]
    out = _sub(arr, 6000.002, 3.0, -0.50)

    total = float(np.log10(np.sum(10.0 ** out["loggf"])))
    assert total == pytest.approx(-0.50, abs=1e-9)
    assert (out["loggf"][0] - out["loggf"][1]) == pytest.approx(ratio_before, abs=1e-12)
    # the same additive shift on every component is what preserves the branching
    shifts = out["loggf"] - before
    assert shifts[0] == pytest.approx(shifts[1], abs=1e-12)


def test_substitution_is_a_copy_and_leaves_the_input_untouched():
    """The caller inverts the ORIGINAL list as its control immediately afterwards. If the
    substitution mutated in place, the control would silently become a second lab-gf
    run and the two numbers would agree for the wrong reason."""
    arr = _fake_linelist([("Fe", 1, "", 6000.000, 3.0, -1.30)])
    before = arr["loggf"].copy()
    out = _sub(arr, 6000.000, 3.0, -0.50)
    assert out["loggf"][0] != before[0]
    assert arr["loggf"][0] == pytest.approx(before[0])


def test_substitution_refuses_when_no_physical_line_matches():
    arr = _fake_linelist([("Fe", 1, "", 6000.000, 3.0, -1.30)])
    with pytest.raises(ValueError, match="matched 0 physical lines"):
        _sub(arr, 7000.000, 3.0, -0.50)


def test_substitution_refuses_an_ambiguous_match_rather_than_picking_one(monkeypatch):
    """Two DIFFERENT physical lines inside the tolerance must be refused, not resolved by
    taking the first — that is how a blend inherits a neighbour's gf (RYA-785: a
    same-species neighbour does not cancel).

    The ambiguity is injected rather than built out of line data, and deliberately so: I
    tried to construct it with two rows at different excitation potentials and got ZERO
    matches, because a fixture close enough in wavelength to be ambiguous is also close
    enough for `cluster_physical_lines` to merge into ONE physical line, and one far
    enough apart to stay separate falls outside the 0.02 eV window. That is a real
    property of the clustering and it means this branch is close to unreachable through
    the data — which is a reason to pin the refusal logic directly, not a reason to skip
    it. The guard is cheap and the failure it prevents is silent.
    """
    import rya824_lab_gf_subpool as mod
    import pipeline.gf_resolver as res

    arr = _fake_linelist([
        ("Fe", 1, "", 6000.000, 3.00, -1.30),
        ("Fe", 1, "", 6000.004, 3.00, -2.10),
    ])
    monkeypatch.setattr(res, "cluster_physical_lines",
                        lambda keys, wls, eps: [[0], [1]])
    with pytest.raises(ValueError, match="matched 2 physical lines"):
        _sub(arr, 6000.002, 3.00, -0.50)


def test_substitution_does_not_touch_rows_outside_the_target():
    """A neighbouring species in the same synthesis window must come out bit-identical,
    or the re-inverted abundance is a function of more than the one gf that changed."""
    arr = _fake_linelist([
        ("Fe", 1, "", 6000.000, 3.0, -1.30),
        ("Ni", 1, "", 6000.050, 4.1, -2.40),
        ("Fe", 1, "", 6002.000, 2.2, -0.90),
    ])
    out = _sub(arr, 6000.000, 3.0, -0.50)
    assert out["loggf"][1] == pytest.approx(arr["loggf"][1])
    assert out["loggf"][2] == pytest.approx(arr["loggf"][2])


def test_a_zero_change_substitution_is_a_no_op():
    """Substituting a line's own current total must return it unchanged. If this drifts,
    every reported dA carries an offset that has nothing to do with the lab value."""
    arr = _fake_linelist([
        ("Fe", 1, "", 6000.000, 3.0, -1.30),
        ("Fe", 1, "", 6000.004, 3.0, -1.70),
    ])
    cur_total = float(np.log10(np.sum(10.0 ** arr["loggf"])))
    out = _sub(arr, 6000.002, 3.0, cur_total)
    assert np.allclose(out["loggf"], arr["loggf"], atol=1e-12)
