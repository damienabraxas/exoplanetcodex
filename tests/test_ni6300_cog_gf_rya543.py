"""
tests/test_ni6300_cog_gf_rya543.py
==================================
RYA-543 — the Ni I 6300.34 gf used by the [O I] 6300 EW/COG subtraction
(`NI6300_COG` → `lines_fit._predict_ni6300_ew`) must be single-sourced through
gf_resolver (canonical = Johansson+2003 log gf −2.11), NOT a hardcoded constant.

The old `NI6300_COG['log_gf'] = −2.841` (stale VALD3) was a duplicate that
diverged 0.73 dex from the RYA-365-adjudicated canonical −2.11, so the EW-path
[O I] Ni-subtraction used the wrong gf. This test pins:
  1. NI6300_COG no longer hardcodes a log gf (SSOT — resolved at use).
  2. the resolver returns −2.11 for that line.
  3. `_predict_ni6300_ew` resolves via gf_resolver (no hardcoded gf in the path).
  4. the stewardship guard now FAILS loudly if a divergent hardcoded gf is
     reintroduced into NI6300_COG (the guard-gap closer).
"""
import inspect

import pytest

import config.constants as const
import scripts.check_stewardship as sc
from pipeline import lines_fit
from pipeline.gf_resolver import resolve as resolve_gf, _index

NI_KEY = (28, 1)
NI_WL, NI_EP = 6300.342, 4.266
JOHANSSON_GF = -2.11
STALE_VALD_GF = -2.841


# ── 1. SSOT: no hardcoded gf in the constants dict ───────────────────────────

def test_ni6300_cog_has_no_hardcoded_loggf():
    assert 'log_gf' not in const.NI6300_COG, (
        "NI6300_COG must not hardcode a log gf — it is single-sourced to "
        "canonical_gf.csv and resolved via gf_resolver (RYA-543).")
    # the metadata the COG still needs is present and mirrors canonical
    assert abs(const.NI6300_COG['excitation_potential_eV'] - NI_EP) < 1e-6
    assert abs(const.NI6300_COG['wavelength_air_A'] - NI_WL) <= 0.02


def test_resolver_returns_johansson_value():
    _index.cache_clear()
    assert abs(resolve_gf(NI_KEY, NI_WL, NI_EP) - JOHANSSON_GF) < 1e-6


# ── 2. the EW/COG path resolves via gf_resolver, no hardcoded gf ──────────────

def test_predict_ni6300_ew_uses_resolver_not_hardcoded_gf():
    src = inspect.getsource(lines_fit._predict_ni6300_ew)
    assert '_gr.resolve(' in src or 'gf_resolver.resolve(' in src, (
        "_predict_ni6300_ew must resolve the Ni 6300.34 gf via gf_resolver.")
    assert "NI6300_COG['log_gf']" not in src, (
        "the hardcoded NI6300_COG['log_gf'] reference must be gone from the COG path.")


# ── 3. stewardship guard closes the gap (the RYA-543 regression) ──────────────

def test_guard_passes_in_target_state():
    # with no hardcoded gf, the const-gf invariant is a clean no-op
    assert sc.check_constants_gf_duplicates() == []


def test_guard_fails_on_reintroduced_stale_gf(monkeypatch):
    # inject the exact stale VALD3 duplicate the ticket describes
    monkeypatch.setitem(const.NI6300_COG, 'log_gf', STALE_VALD_GF)
    viols = sc.check_constants_gf_duplicates()
    assert len(viols) == 1
    v = viols[0]
    assert v.invariant == 'const_gf'
    assert not v.tracked, "a divergent hardcoded gf must be UNTRACKED → FAIL the build"
    assert 'NI6300_COG' in v.locus
    assert f"{STALE_VALD_GF:+.3f}" in v.value and f"{JOHANSSON_GF:+.3f}" in v.value


def test_guard_fails_on_any_divergent_gf(monkeypatch):
    # a small divergence above the 0.05-dex threshold still fails
    monkeypatch.setitem(const.NI6300_COG, 'log_gf', JOHANSSON_GF - 0.2)
    assert len(sc.check_constants_gf_duplicates()) == 1


def test_guard_allows_canonical_value(monkeypatch):
    # a hardcoded copy that MATCHES canonical is not a divergence (no false positive)
    monkeypatch.setitem(const.NI6300_COG, 'log_gf', JOHANSSON_GF)
    assert sc.check_constants_gf_duplicates() == []
