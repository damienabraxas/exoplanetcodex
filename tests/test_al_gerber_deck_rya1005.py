"""tests/test_al_gerber_deck_rya1005.py — RYA-1005

Al's TS-native Gerber deck, and the abundance axis that Fe's deck does not have.

The thing worth guarding here is not that Al is registered. It is that `gerber_nlte` was
built on a property of the Fe grid — "the departures do not depend on A(X)" — which the
module measured, stated, and correctly relied on, and which is FALSE for Al. Adding Al
without the axis would have returned well-formed departure files at one pinned abundance
for every trial in a chi2 loop, which is invisible in the output.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline import gerber_nlte as g          # noqa: E402

_GT = Path(g.GT)
needs_grids = pytest.mark.skipif(
    not (_GT / g.DECKS["Al"]["aux"]).exists(),
    reason="Gerber grids are Sirius-only (external drive)")


def test_al_is_registered_with_its_staged_assets():
    d = g.DECKS["Al"]
    assert d["Z"] == 13
    assert d["atom"] == "atom.al_qmh"
    assert "MARCS" in d["grid"]


def test_the_provenance_record_exists_and_refuses_a_scalar_deck_abundance():
    """`deck_abundance()` is the Fe concept. Al has 31 abundances, so the record carries
    `a_sun: null` DELIBERATELY and says why — a number there would be a fiction that
    every caller would then interpolate at."""
    import json
    p = ROOT / "data/nlte_grids/gerber_ts/Al_gerber2023.prov.json"
    assert p.exists(), "an unregistered deck has not passed the RYA-534/785 gate"
    rec = json.loads(p.read_text())
    assert rec["deck_abundance"]["a_sun"] is None
    assert rec["abundance_axis"]["n"] == 31


@needs_grids
def test_al_has_an_abundance_axis_and_fe_does_not():
    """🔴 THE STRUCTURAL DIFFERENCE. Fe: 15,229 aux rows, ONE A(X) at [Fe/H]=0.
    Al: 454,466 rows, THIRTY-ONE. Everything else in this file follows from it."""
    al, fe = g.abundance_axis("Al"), g.abundance_axis("Fe")
    assert len(fe) == 1 and not g.has_abundance_axis("Fe")
    assert len(al) == 31 and g.has_abundance_axis("Al")
    assert al[0] == pytest.approx(4.43) and al[-1] == pytest.approx(7.43)
    # the solar value must be ON the grid, not extrapolated to
    assert al[0] < 6.43 < al[-1]


@needs_grids
def test_al_refuses_without_an_abundance():
    """POSITIVE CONTROL on the guard. Silently pinning the chi2 loop to one arbitrary
    abundance is the failure this exists to prevent, and it would look fine."""
    with pytest.raises(g.GerberDeckError, match="abundance is required"):
        g.for_node("Al", 5750.0, 4.50, 0.0)


@needs_grids
def test_al_refuses_off_the_axis():
    """Extrapolating departures off the grid is not a correction, it is an invention."""
    with pytest.raises(g.GerberDeckError, match="outside this deck's abundance axis"):
        g.for_node("Al", 5750.0, 4.50, 0.0, abundance=9.9)


@needs_grids
def test_fe_still_works_with_no_abundance_passed():
    """INERTNESS. Every pre-RYA-1005 caller passes no abundance; Fe must be unchanged —
    same deck abundance and same shape.

    🔴 CORRECTED RYA-1035 — 7.46 → 7.50, and the parenthetical was backwards. It read
    "7.46, NOT the 7.50 in its aux table"; the aux table was right and the record was
    wrong. `abu_ref` is read from stdin (interpol_modeles_nlte.f:206), written verbatim
    into the departure file (:761) and printed back by bsyn (:988), and `gerber_nlte` fed
    that stdin from `deck_abundance()` itself — so the 7.46 was our own input echoing
    round a closed loop. The grid says 7.50: `atom.fe607a` line 2 (`7.50  55.85`), both
    aux tables' A(X) = 7.50 + [Fe/H], and Turbospectrum's own `metal = abund(15) - 7.50`.

    The INERTNESS this test is really about is unaffected: what changed is a label, not a
    departure. `deck_abundance` now raises if the record and the aux ever disagree again."""
    d = g.for_node("Fe", 5750.0, 4.50, 0.0)
    assert d["deck_abundance"] == pytest.approx(7.50)
    assert d["ndep"] == 56 and d["nk"] == 607


@needs_grids
@pytest.mark.skipif(not __import__("os").environ.get("RYA1005_SLOW"), reason="3 interpolator runs; set RYA1005_SLOW=1")
def test_al_departures_actually_differ_between_abundances():
    """The measurement the whole design rests on, kept runnable.

    If this ever returns identical blocks, Al's axis is nominal after all and the Fe
    approximation should be restored — so the assertion is the claim, not a formality.
    """
    import hashlib
    shas = set()
    for A in (6.20, 6.43, 6.70):
        d = g.for_node("Al", 5750.0, 4.50, 0.0, abundance=A)
        shas.add(hashlib.sha256(d["departures"].tobytes()).hexdigest())
    assert len(shas) == 3, "Al's departures did NOT vary with abundance"
