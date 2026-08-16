"""
tests/test_gf_column_honesty_rya825.py — RYA-825
================================================
A `log_gf` column must report the gf the inversion USED, not the one the linelist
happened to carry at intake. The two diverged silently for a long time and the cost was
not abstract: RYA-799 read the intake column, concluded the Fe pool was "100 % Kurucz
K14", and built a 48-line SCALE-MISMATCH finding on it. RYA-824 re-inverted and found 18
of 29 lab-covered lines were already on the canonical value; re-grading on the corrected
column turns 37 of those 48 into pure metadata artifacts and takes the graded count from
2 to 25.

So this file pins three things:

  * the live accounting table is honest — every row that claims `gf_canonical` really
    does carry the resolver's value (the drift guard, asserted on the artifact rather
    than on a fixture, because a fixture cannot drift);
  * the intake value is preserved, so nothing was destroyed to achieve that;
  * the guard can FAIL — a column nudged off the resolver's answer is caught.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.gf_resolver import (  # noqa: E402
    GfResolutionError, annotate_used_gf, assert_gf_column_is_honest, resolve)
from pipeline.species import species_key  # noqa: E402

ACCOUNTING = ROOT / "data" / "audit" / "line_accounting" / "per_line.csv"
KW = dict(wl_col="wave_air_A", ep_col="ep_eV", gf_col="log_gf")


@pytest.fixture(scope="module")
def acc():
    return pd.read_csv(ACCOUNTING)


# ── the live artifact ─────────────────────────────────────────────────────────

def test_the_accounting_table_carries_the_provenance_columns(acc):
    for col in ("log_gf", "log_gf_intake", "gf_source_intake", "gf_canonical"):
        assert col in acc.columns, f"{col} missing — the column cannot be audited"


def test_the_reported_gf_is_the_resolver_value_wherever_it_claims_to_be(acc):
    """THE DRIFT GUARD. Asserted on the live table, not a fixture."""
    stats = assert_gf_column_is_honest(acc, **KW)
    assert stats["n_checked"] > 0, "no row claims gf_canonical — the guard checked nothing"
    assert stats["n_bad"] == 0


def test_the_intake_value_was_preserved_not_destroyed(acc):
    """RYA-825 forbids dropping the intake gf. It must be finite everywhere the primary
    column is, or the correction traded one blind spot for another."""
    both = acc[np.isfinite(acc.log_gf)]
    assert np.isfinite(both.log_gf_intake).all()
    assert (acc.gf_source_intake.astype(str).str.strip() != "").all()


def test_unresolved_rows_keep_their_intake_value_unchanged(acc):
    """Outside canonical's range the resolver does not run, and the honest report there
    is the intake value untouched — not a gap, and not a guess."""
    un = acc[~acc.gf_canonical.astype(bool)]
    assert len(un) > 0, "expected some rows outside canonical coverage"
    assert np.allclose(un.log_gf.values, un.log_gf_intake.values, equal_nan=True)


def test_the_correction_actually_changed_something(acc):
    """A guard that passes because nothing was corrected proves nothing. The whole
    premise is that the two columns disagree on a substantial minority of rows."""
    d = np.abs(acc.log_gf - acc.log_gf_intake)
    assert int((d > 1e-6).sum()) > 100


# ── the guard discriminates ───────────────────────────────────────────────────

def test_the_drift_guard_catches_a_nudged_column(acc):
    """The failure mode is a column that drifts off the resolver by a small amount. If
    the guard cannot see that, it is decoration."""
    bad = acc.copy()
    j = bad.index[bad.gf_canonical.astype(bool)][0]
    bad.loc[j, "log_gf"] = float(bad.loc[j, "log_gf"]) + 0.05
    with pytest.raises(GfResolutionError, match="drifted"):
        assert_gf_column_is_honest(bad, **KW)


def test_the_guard_refuses_a_table_with_no_canonical_flag(acc):
    """Without the flag there is no way to know which rows are claiming to be
    resolver-sourced, so the guard must refuse rather than pass vacuously."""
    with pytest.raises(GfResolutionError, match="cannot be checked"):
        assert_gf_column_is_honest(acc.drop(columns=["gf_canonical"]), **KW)


def test_the_guard_catches_COVERAGE_GROWTH_not_just_drift(acc):
    """The half the first version of this guard was missing.

    Validating only the rows that CLAIM canonical is blind to canonical GAINING lines:
    rows correctly flagged False go stale and nothing notices. That is not hypothetical.
    RYA-822 extended canonical_gf blueward the same day RYA-825 landed, and 3,413
    accounting rows became resolvable while the freshly-corrected table still reported
    them as outside coverage — a one-sided guard let a just-fixed table rot within hours.

    Simulated by flipping a resolvable row's flag to False, which is exactly the state
    coverage growth produces.
    """
    stale = acc.copy()
    j = stale.index[stale.gf_canonical.astype(bool)][0]
    stale.loc[j, "gf_canonical"] = False
    with pytest.raises(GfResolutionError, match="gained coverage"):
        assert_gf_column_is_honest(stale, **KW)


def test_the_live_table_is_current_with_todays_canonical(acc):
    """Not just self-consistent — CURRENT. This is the assertion that would have failed
    the moment RYA-822 landed, and it is the point of the two-sided guard."""
    stats = assert_gf_column_is_honest(acc, **KW)
    assert stats["n_newly_resolvable"] == 0


# ── the annotator ─────────────────────────────────────────────────────────────

def test_annotate_used_gf_round_trips_a_known_line(acc):
    """Take a row the table says is canonical and re-derive it from scratch."""
    row = acc[acc.gf_canonical.astype(bool)].iloc[0]
    key = species_key(str(row.element), str(row.ion))
    assert float(row.log_gf) == pytest.approx(
        resolve(key, float(row.wave_air_A), float(row.ep_eV)), abs=1e-9)


def test_annotate_used_gf_is_idempotent(acc):
    """Running the correction twice must not move anything — otherwise regenerating the
    table would keep nudging the column."""
    once, _ = annotate_used_gf(acc[["element", "ion", "wave_air_A", "log_gf", "ep_eV"]]
                               .head(400), **KW)
    twice, _ = annotate_used_gf(once, **KW)
    assert np.allclose(once.log_gf.values, twice.log_gf.values, equal_nan=True)
