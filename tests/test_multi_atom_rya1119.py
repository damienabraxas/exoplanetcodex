"""RYA-1119 — the MULTI -> Lightweaver atom converter.

🔴 WHAT THESE TESTS ARE FOR, AND WHAT THEY ARE NOT.

They pin the FORMAT and the UNIT CONVENTIONS, because those are objectively checkable and
because every one of them was got wrong on the first attempt, and three of the four did
NOT raise:

  * continuum grids are DESCENDING in MULTI, ascending in Lightweaver  -> raises
  * wavelengths are ANGSTROM in MULTI, NANOMETRE in Lightweaver        -> raises, far away
  * cross-sections are CM^2 in MULTI, M^2 in Lightweaver               -> SILENT, 1e4 error
  * collisional rates are CM^3 in MULTI, M^3 in Lightweaver            -> SILENT, 1e6 error

The cm^3 one is the instructive failure: it made every collisional rate a million times too
strong, which thermalised the atom completely and returned b = 1.0000 at every level and
every depth. It converged in five iterations and reported no error. A converged LTE answer
wearing an NLTE label — the RYA-1118 striding lesson in a different costume, caught by
looking at the departure STRUCTURE rather than at whether the solver finished.

⚠️ THEY DO NOT VALIDATE THE PHYSICS. No test here shows the restricted atom's departures are
correct. That gate is stated on RYA-1137 and is a separate job: reproduce fe607a's own
departures on the 40 AGSS21 lines through the Gerber route, and land A(Fe) on 7.46. Passing
this file means the conversion is faithful, not that the atom is right.
"""
from pathlib import Path

import numpy as np
import pytest

from pipeline.multi_atom import (CM1_TO_EV, read_multi_atom, levels_touched, restrict)

ATOM = Path.home() / "Documents/Exoplanet Codex/grids/nlte/gerber_ts/atom.fe607a"
pytestmark = pytest.mark.skipif(not ATOM.exists(),
                                reason="atom.fe607a is a vendor deck, not in the repo")


@pytest.fixture(scope="module")
def atom():
    return read_multi_atom(ATOM)


def test_the_header_counts_are_what_was_parsed(atom):
    """A MULTI file declares its own sizes. Parsing fewer is a silent truncation."""
    assert len(atom.levels) == 607
    assert len(atom.lines) == 12635
    assert len(atom.continua) == 606
    assert len(atom.collisions) == 29983


def test_continua_are_walked_not_indexed(atom):
    """🔴 Continua are VARIABLE length — header then NLAMB pairs.

    Indexing them at a fixed stride reads (wavelength, alpha) pairs as continuum headers
    and every later block is misaligned, which shows up as garbage collision indices
    rather than as an exception.
    """
    for c in atom.continua:
        assert len(c.wavelength) == c.nlamb
        assert len(c.alpha) == c.nlamb
    # and the block after them parsed cleanly, which is what proves the walk landed right
    kinds = {c.kind for c in atom.collisions}
    assert kinds == {"CE", "CH", "CI", "CH0"}, kinds


def test_stage_is_one_based_in_the_file(atom):
    """MULTI counts 1 = neutral. Lightweaver counts 0 = neutral (checked separately)."""
    stages = {l.stage for l in atom.levels}
    assert stages == {1, 2, 3}
    assert atom.level(1).stage == 1 and atom.level(1).E_cm == 0.0


def test_fe607a_covers_the_agss21_line_set(atom):
    """The premise of the restriction: RYA-1137 measured 40/40 at 1 meV."""
    import csv
    p = Path(__file__).resolve().parents[1] / "data/reference/asplund2021_fe/asplund2021_fe_lines.csv"
    if not p.exists():
        pytest.skip("AGSS21 reference table not present in this worktree")
    ag = [r for r in csv.DictReader(p.open()) if r["ion"] == "I"]
    touched, missed = levels_touched(atom, [(float(r["elo_eV"]), float(r["eup_eV"]))
                                            for r in ag])
    assert missed == [], f"{len(missed)} AGSS21 lines not found in fe607a"
    assert len(touched) == 58, len(touched)


def test_restriction_keeps_the_ion_ground_or_there_is_no_ionisation(atom):
    """🔴 Without the ion ground there is no continuum and no ionisation balance, and the
    atom becomes a closed neutral system that converges happily to the wrong answer."""
    sub, _, rep = restrict(atom, [1, 2, 3])
    assert any(l.stage > 1 for l in sub.levels)
    assert rep.n_levels_after == 4                      # three asked for, plus the ion


def test_the_restriction_reports_what_it_dropped(atom):
    """Restriction is a physics change; it has to leave a record, not be implicit."""
    sub, _, rep = restrict(atom, [1, 2, 3])
    assert rep.n_levels_before == 607 and rep.n_levels_after == 4
    assert len(rep.dropped_levels) == 603
    text = rep.describe()
    assert "physics change" in text and "NOT validated" in text


# ── the unit conventions, which is where every silent error lived ──────────────

@pytest.fixture(scope="module")
def converted(atom):
    lw = pytest.importorskip("lightweaver")
    import csv
    p = Path(__file__).resolve().parents[1] / "data/reference/asplund2021_fe/asplund2021_fe_lines.csv"
    if not p.exists():
        pytest.skip("AGSS21 reference table not present in this worktree")
    from pipeline.multi_atom import to_lightweaver
    ag = [r for r in csv.DictReader(p.open()) if r["ion"] == "I"]
    touched, _ = levels_touched(atom, [(float(r["elo_eV"]), float(r["eup_eV"])) for r in ag])
    sub, remap, _ = restrict(atom, touched)
    return sub, to_lightweaver(sub, remap)


def test_stage_is_zero_based_after_conversion(converted):
    """🔴 Off by one here does not crash — it produces a convergent atom with the wrong
    ionisation balance."""
    _, lw = converted
    stages = sorted({l.stage for l in lw.levels})
    assert stages == [0, 1]
    assert sum(1 for l in lw.levels if l.stage == 0) == 58


def test_continuum_wavelengths_are_nanometres_and_ascending(converted):
    """Angstrom in put every point above `lambdaEdge`; the grid filtered to EMPTY and
    surfaced as `IndexError: index -1` deep inside `compute_wavelength_grid`."""
    sub, lw = converted
    src = {(c.j, c.i): c for c in sub.continua}
    kept = sorted(l.index for l in sub.levels)
    for c in lw.continua:
        w = np.asarray(c.wavelengthGrid)
        assert np.all(np.diff(w) > 0), "wavelength grid must be strictly ascending"
        # Referenced to the SOURCE FILE, not to a magic bound. A first version asserted
        # `< 1000 nm` and failed on a legitimate 1011.8 nm edge from a high-lying level
        # sitting 1.23 eV below the ionisation limit — the threshold was wrong, not the
        # conversion. Lightweaver also appends the computed edge to the grid, so compare
        # against the tabulated points only.
        s = src[(kept[c.j], kept[c.i])]
        tabulated = np.sort(s.wavelength) * 0.1
        assert np.allclose(w[:tabulated.size], tabulated, rtol=1e-12), \
            "tabulated grid is not the source wavelengths in nm"


def test_continuum_cross_sections_are_converted_to_m2(converted):
    """🔴 SILENT 1e4. Checked against the source file rather than a magic number."""
    sub, lw = converted
    src = {(c.j, c.i): c for c in sub.continua}
    kept = sorted(l.index for l in sub.levels)
    for c in lw.continua:
        s = src[(kept[c.j], kept[c.i])]
        # 🔴 atol=0 IS LOAD-BEARING. `np.isclose` defaults to atol=1e-8, and these values
        # are ~1e-20, so the default makes ANY two of them "close" — the first version of
        # this test passed with the conversion reverted, i.e. it asserted nothing. Caught
        # by mutation-testing the constant, not by reading the test.
        assert np.isclose(max(c.alphaGrid), s.alpha.max() * 1e-4, rtol=1e-12, atol=0.0)


def test_collisional_rates_are_converted_to_m3(converted):
    """🔴 SILENT 1e6, and the one that produced a fully thermalised b = 1.0000 atom."""
    sub, lw = converted
    src = {(c.kind, c.j, c.i): c for c in sub.collisions}
    kept = sorted(l.index for l in sub.levels)
    n = 0
    for c in lw.collisions:
        key = type(c).__name__
        kinds = {"CE": "CE", "CI": "CI", "CH": "CH", "ChargeExchangeNeutralH": "CH0"}
        s = src.get((kinds[key], kept[c.j], kept[c.i]))
        if s is None:
            continue
        # atol=0 for the same reason as the cross-sections above: these rates reach
        # ~1e-14 and the default atol=1e-8 would swamp them. This one happened to catch
        # its mutation anyway (raw values are ~1e-7, above the default atol) — happening
        # to work is not the same as being right.
        assert np.allclose(np.asarray(c.rates), s.rates * 1e-6, rtol=1e-12, atol=0.0)
        n += 1
    assert n > 100, f"only {n} rate sets cross-checked"


def test_CH0_is_charge_exchange_not_plain_neutral_H_collisions(converted):
    """CH0 is downward-only charge exchange. Mapping it to CH invents an upward rate."""
    from lightweaver.collisional_rates import CH, ChargeExchangeNeutralH
    sub, lw = converted
    n_ch0 = sum(1 for c in sub.collisions if c.kind == "CH0")
    if n_ch0:
        assert sum(1 for c in lw.collisions
                   if isinstance(c, ChargeExchangeNeutralH)) == n_ch0
