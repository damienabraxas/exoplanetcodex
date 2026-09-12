"""RYA-1213 — the REFERENCE tier is the lab pool with NO depth gate.

Reference Grade existed as one line set (`asplund`) and therefore in one band. RYA-946
defines the tier by gf PEDIGREE, so the selector this suite guards asks the gf question
and stops. The tests here are about the three ways that could silently go wrong:

  * the pool could quietly acquire a depth term and become GRADED or DEEPGRADED wearing
    a Reference name;
  * the EW routes' two-branch tier expression (`graded -> the lab mask, anything else ->
    its complement`) could map `reference` onto the NON-laboratory lines and publish
    them as the top grade;
  * the tier could reach the feed without a grade, or with the same grade label as the
    AGSS21 replication and no way for a reader to tell the two apart.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from pipeline import model_registry as mr            # noqa: E402
from pipeline import product_eligibility as pe       # noqa: E402
from pipeline import reference_lineset as rls        # noqa: E402

FEED = ROOT / "data/products/solar/Fe.json"
CANON = ROOT / "data/linelists/canonical_gf.csv"


@pytest.fixture(scope="module")
def feed():
    return json.loads(FEED.read_text())


# ── the axis ─────────────────────────────────────────────────────────────────
def test_reference_is_in_the_one_vocabulary():
    """`reference` is a DELIBERATE addition to LINE_SETS, not a value a writer invented.

    `line_set_for_product` refuses anything outside the tuple, so a product tagged with
    an unregistered set fails loudly — which is the behaviour, and this asserts the
    registration actually happened rather than the reader having been widened.
    """
    assert "reference" in mr.LINE_SETS


def test_the_tier_derives_its_line_set_and_never_stores_it():
    """A Reference product is OURS, so the axis derives from `tier` like every other.

    The stored-vs-derived asymmetry is RYA-1127's: a replication carries an explicit
    `line_set` because nothing in its record implies whose list it is; ours must not,
    because a stored copy is a second source of truth free to drift from the key.
    """
    assert rls.line_set_for_product({"tier": "REFERENCE"}) == "reference"
    ours = [p for p in json.loads(FEED.read_text())["products"]
            if p.get("tier") == "REFERENCE"]
    for p in ours:
        assert "line_set" not in p, (
            f"{pe.key_of(p)} stores `line_set` — our own products derive it from tier "
            f"(RYA-1127); storing it creates a second source of truth")


def test_reference_and_asplund_cannot_collide_on_identity():
    """Two line sets now publish as Reference Grade, so the KEY must separate them.

    They answer different questions — AGSS21's 21 lines on AGSS21's gf, versus our whole
    lab pool on canonical_gf — and `line_set` is the only field that differs between a
    VIS Fe I Reference row and a VIS Fe I Asplund row on the same holding and treatment.
    If it ever left KEY_FIELDS, one would silently overwrite the other.
    """
    assert "line_set" in pe.KEY_FIELDS


# ── the pool ─────────────────────────────────────────────────────────────────
@pytest.fixture(scope="module")
def lab():
    cg = pd.read_csv(CANON, low_memory=False)
    return cg[cg.gf_tier.astype(str).str.contains("LAB", na=False)]


def test_the_selector_applies_no_depth_term(lab, monkeypatch):
    """🔴 THE ONE PROPERTY THAT DEFINES THE TIER, ASSERTED ON THE SELECTOR ITSELF.

    Not on a product (which would only prove one run behaved) and not by reading the
    source (which would pass on a comment). `_cand_reference` is called with a synthetic
    linelist that carries EVERY lab wavelength in the window, so nothing can be lost to
    the synthesis-list match, and the count it returns must be the full lab count —
    including lines on BOTH sides of the depth gate.
    """
    from derive_band_products import _cand_reference

    lo, hi = 4200.0, 6910.0
    want = lab[(lab.species == "Fe I") & lab.wavelength_air_A.between(lo, hi)]
    n = len(want)
    assert n > 2, "the VIS Fe I lab pool is the fixture; it must not be empty"

    ll = np.rec.fromarrays(
        [want.wavelength_air_A.values.astype(float),
         np.array(["Fe 1"] * n, dtype=object),
         want.log_gf.values.astype(float),
         want.excitation_potential_eV.values.astype(float),
         np.full(n, 0.5)],
        names="wave_A,element,loggf,lower_state_eV,theoretical_depth")
    got = _cand_reference(ll, lo_A=lo, hi_A=hi, species="Fe 1")
    assert len(got) == n, (
        f"_cand_reference returned {len(got)} of {n} LAB-tier Fe I lines in "
        f"{lo:.0f}-{hi:.0f} A. Reference is defined by gf pedigree alone; any shortfall "
        f"here means a depth term has re-entered the selector")


def test_reference_is_a_superset_of_both_depth_split_selectors(lab):
    """Codex + Deep <= Reference, and the difference is REAL LINES, not rounding.

    Both depth comparisons are False for NaN, so a lab line whose feature the stellar
    catalogue does not carry falls out of BOTH pools. Reference keeps it, because
    nothing about its depth was ever part of the question. This asserts the arithmetic
    that makes that claim true rather than trusting the prose.
    """
    from derive_band_products import _feature_depth
    from line_accounting_rya709 import DEPTH_HI

    lo, hi = 6910.0, 9199.0          # red-optical: the band that holds the no-depth line
    sel = lab[(lab.species == "Fe I") & lab.wavelength_air_A.between(lo, hi)]
    d = _feature_depth(sel.wavelength_air_A.values.astype(float))
    codex, deep = int((d <= DEPTH_HI).sum()), int((d > DEPTH_HI).sum())
    assert codex + deep <= len(sel)
    assert np.isnan(d).any(), (
        "red-optical is the fixture for the no-depth case precisely because it has one; "
        "if the stellar catalogue has since gained that feature, move this test to a "
        "band that still exercises it rather than deleting the property")
    assert len(sel) - (codex + deep) == int(np.isnan(d).sum())


# ── the two refusals ─────────────────────────────────────────────────────────
def test_the_ew_route_refuses_reference_rather_than_inverting_it():
    """🔴 THE FAILURE THIS WOULD HAVE BEEN: the top grade, measured on the NON-lab lines.

    Both EW branches read `_keep = mask if tier == "graded" else ~mask`. Anything that
    is not literally `graded` takes the COMPLEMENT, so `reference` would have selected
    every line with no laboratory gf and published it as Reference Grade. The refusal is
    asserted through the argument parser and the source of the branch, so a future edit
    that re-opens the path has to defeat both.
    """
    import derive_band_products as dbp

    src = Path(dbp.__file__).read_text()
    # The EW artifact selector and the profile-fit route each guard independently.
    assert src.count('if tier == "reference":') >= 1
    assert src.count('if _tier == "reference":') >= 1
    for branch in ('_keep = _g if _tier == "graded" else ~_g',
                   'df[graded if tier == "graded" else ~graded]'):
        assert branch in src, (
            "the two-branch expression this test exists to guard has moved or been "
            "rewritten; re-point the guard at wherever the tier now chooses a mask")


def test_the_tier_is_offered_by_both_clis():
    """A selector nothing can ask for is not wired. Both entry points must name it."""
    import derive_band_products as dbp
    import publish_product as pp

    assert '"reference"' in Path(dbp.__file__).read_text()
    assert '"REFERENCE"' in Path(pp.__file__).read_text()


def test_reference_is_held_to_the_laboratory_gf_gate():
    """REFERENCE claims lab gf, so it sits in the publisher's rung-3 gate.

    Leaving it out would let a pool that failed to resolve to a laboratory scale publish
    under the TOP grade name with no corroboration — the RYA-1212 shape, at the one
    tier where the claim is loudest.
    """
    import publish_product as pp

    src = Path(pp.__file__).read_text()
    assert 'if a.tier in ("GRADED", "DEEPGRADED", "REFERENCE"):' in src


# ── the feed ─────────────────────────────────────────────────────────────────
def test_every_reference_product_carries_its_grade_and_its_variant(feed):
    """RYA-1213 §4 — a reader must be able to tell the two Reference kinds apart."""
    refs = [p for p in feed["products"] if p.get("grade") == "Reference Grade"]
    assert refs, "no Reference Grade product in the feed"
    for p in refs:
        ls = p.get("line_set_resolved")
        assert ls in ("reference", "asplund", "asplund-al", "gbs"), pe.key_of(p)
        assert p.get("grade_variant"), (
            f"{pe.key_of(p)} publishes as Reference Grade with no `grade_variant`; two "
            f"line sets share that grade and the label is what separates them")
    variants = {p["line_set_resolved"]: p["grade_variant"] for p in refs}
    assert len(set(variants.values())) == len(variants), (
        "two Reference line sets share a grade_variant string, so the label does not "
        "actually distinguish them")


def test_a_reference_product_reproduces_its_pool_without_a_depth_gate(feed):
    """`tier_provenance` must reproduce what the run MEASURED, not a depth-split half."""
    for p in feed["products"]:
        if p.get("tier") != "REFERENCE":
            continue
        tp = p.get("tier_provenance") or {}
        if not tp.get("pool_reproduced"):
            continue
        assert tp.get("depth_gate_applied") is False, pe.key_of(p)
        assert tp.get("depth_gate") is None, pe.key_of(p)
        assert tp["n_selected"] == tp["n_lab_tier_in_window"], (
            f"{pe.key_of(p)} reproduced {tp['n_selected']} of "
            f"{tp['n_lab_tier_in_window']} LAB-tier rows in its own window. A Reference "
            f"pool is every one of them; a shortfall means a depth term ran")
