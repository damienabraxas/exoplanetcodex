"""
tests/test_round_dex_determinism_rya1084.py — RYA-1084
======================================================
A published uncertainty moved 0.0217 -> 0.0218 between two runs that agreed on every other
field. The ticket's prior was a Python `round()` vs numpy `np.round()` half-tie. It is not:
they agree at every point. The real mechanism is that 0.02175 as a float is
0.0217499999999999985 — one ULP BELOW the decimal tie — so `round()` returns 0.0217 while
the very next representable float returns 0.0218, and `stat_basis` prints five decimals and
says "0.02175" for both.

These tests pin the mechanism, so nobody re-derives the wrong hypothesis from the artifact.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.error_budget import round_dex, DEX_PLACES  # noqa: E402

TIE = 0.02175
TIE_NEXT = float(np.nextafter(TIE, 1))


def test_the_float_literal_sits_below_the_decimal_tie():
    """The whole mechanism in one assertion. 0.02175 is not 0.02175."""
    assert TIE < 0.02175000000000000  or f"{TIE:.20f}".startswith("0.02174999")
    assert f"{TIE:.20f}".startswith("0.021749999")


def test_python_round_and_numpy_round_agree_AT_THIS_VALUE():
    """🔴 The ticket's leading prior, refuted precisely rather than sweepingly.

    round() and np.round() DO disagree in general — 0.12345 gives 0.1235 vs 0.1234, and
    5e-05 gives 0.0001 vs 0.0 — which is exactly why the standing rule (RYA-1033) exists.
    But at 0.02175 and its ULP neighbour they AGREE, so a round()/np.round() disagreement
    cannot be what moved this bar. Getting that distinction wrong would have sent the next
    person after the wrong defect."""
    for x in (TIE, TIE_NEXT, 0.0217499, 0.0217501):
        assert round(x, 4) == float(np.round(x, 4)), x
    # and the general disagreement is real, which is why the rule had to be pinned at all
    assert round(0.12345, 4) != float(np.round(0.12345, 4))


def test_one_ulp_flips_the_published_digit_under_plain_round():
    """The defect, reproduced: two floats one ULP apart round to different published bars."""
    assert round(TIE, 4) == 0.0217
    assert round(TIE_NEXT, 4) == 0.0218


def test_stat_basis_precision_cannot_distinguish_them():
    """Why it looked like a tie: the basis string is 5 dp and shows 0.02175 for both, so the
    artifact carries no evidence of which float produced it."""
    assert f"{TIE:.5f}" == f"{TIE_NEXT:.5f}" == "0.02175"


def test_round_dex_is_stable_across_that_ulp_pair():
    """THE FIX. Half-up on the shortest round-tripping decimal gives the same answer for
    both members of the pair that `round()` splits."""
    assert round_dex(TIE) == round_dex(TIE_NEXT) == 0.0218


def test_round_dex_still_separates_genuinely_different_values():
    """It must not paper over a real difference — only the tie neighbourhood is affected."""
    assert round_dex(0.0217499) == 0.0217
    assert round_dex(0.0217501) == 0.0218


def test_round_dex_changes_nothing_but_ties():
    """A rounding change that moved ordinary values would be a silent re-publication of
    every bar in the repo. It differs from round() on ties only."""
    rng = np.random.default_rng(0)
    v = rng.uniform(0, 2, 200_000)
    assert not [x for x in v if abs(round_dex(float(x)) - round(float(x), 4)) > 1e-12]


def test_the_rule_is_single_sourced_at_the_write_sites():
    """🔴 The bug was possible because `round(x, 4)` was written at five separate call
    sites: a rule declared five times is a rule nobody owns (the RYA-845 lesson)."""
    src = (ROOT / "scripts" / "derive_band_products.py").read_text()
    assert "round_dex(" in src
    for bad in ("round(stat, 4)", "round(syst, 4)", "round(eb_stat, 4)",
                "round(stat_a, 4)"):
        assert bad not in src, f"a raw {bad} came back — route it through round_dex"
    assert DEX_PLACES == 4


def test_the_ten_blocked_products_were_re_ingested_not_papered_over():
    """The resolution must be visible in the feed, with a reason attached to each row."""
    import json
    doc = json.loads((ROOT / "data/products/solar/Fe.json").read_text())
    re_ing = [p for p in doc["products"]
              if p["provenance"].get("reingested_by") == "RYA-1084"]
    assert len(re_ing) == 10
    assert all("deterministic" in p["provenance"]["reingest_reason"] for p in re_ing)


def test_no_value_was_special_cased():
    """🔴 CRITICAL condition: pin one rounding rule, never special-case 0.02175."""
    import ast
    for f in ("pipeline/error_budget.py", "scripts/derive_band_products.py",
              "scripts/rya1084_reingest_blocked.py"):
        tree = ast.parse((ROOT / f).read_text())
        # Constants only: prose may name the value, executable code may not branch on it.
        consts = [n.value for n in ast.walk(tree)
                  if isinstance(n, ast.Constant) and isinstance(n.value, (int, float))]
        assert not [c for c in consts if abs(float(c) - 0.02175) < 1e-9], (
            f"{f} carries 0.02175 as a literal — pin one rounding rule, never a value")
