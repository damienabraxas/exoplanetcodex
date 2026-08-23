"""An empty census is a configuration error, not a result (RYA-1020).

`rya824_lab_gf_subpool.census()` used to end in `return pd.DataFrame(rows)`. When no
pool matched, `rows` was empty, and `pd.DataFrame([])` carries **no columns** -- so a
wrong `--pool-dir` did not fail there. It failed roughly a hundred lines later in the
census *printer*, as:

    AttributeError: 'DataFrame' object has no attribute 'band'

which reads as a schema problem in the pool CSV and sends the reader to the data. The
actual fault was upstream and mundane: the function matched nothing, and said so only in
`[skip]` lines that scroll past.

This pins the fix at the fault: refuse, name the directory, and say how many pools were
sought and how many were missing. The distinction the message must preserve is
**"not present" vs "present but zero matches"** -- the first is a path mistake, the
second is a real (and much more interesting) scientific absence, and collapsing them
into one error would send the next person down the wrong road.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _mod():
    spec = importlib.util.spec_from_file_location(
        "rya824_lab_gf_subpool", ROOT / "scripts" / "rya824_lab_gf_subpool.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    return m


def test_empty_census_raises_instead_of_returning_a_columnless_frame(tmp_path):
    """🔴 The regression: an empty pool dir must fail HERE, not as an AttributeError later."""
    m = _mod()
    with pytest.raises(SystemExit) as e:
        m.census(tmp_path)                      # empty dir — nothing can match
    msg = str(e.value)

    # names the directory, so the reader knows which path was wrong
    assert str(tmp_path) in msg, f"error does not name the pool dir: {msg}"

    # and distinguishes a path mistake from a real scientific absence
    assert "not present" in msg and "0 matches" in msg, (
        "the message must separate 'pool file missing' (a path mistake) from 'pool "
        f"present but nothing matched' (a real absence): {msg}")


def test_the_old_failure_mode_cannot_return(tmp_path):
    """Guard the shape, not the string: whatever census() does on an empty dir, it must
    never hand back an object whose `.band` access is what finally raises."""
    m = _mod()
    try:
        out = m.census(tmp_path)
    except SystemExit:
        return                                   # refused up front — correct
    assert hasattr(out, "band"), (
        "census() returned a frame with no `band` column. That is the RYA-1020 bug: the "
        "failure surfaces later, in the printer, as an AttributeError about the pool "
        "schema rather than about the empty census."
    )
