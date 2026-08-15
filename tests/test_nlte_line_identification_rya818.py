"""
tests/test_nlte_line_identification_rya818.py — RYA-818
=======================================================
The arithmetic that decides whether an atom is term-resolved, and the matching
that follows from the answer.

These run WITHOUT Sirius: the atom-shaped fixtures below are synthetic. The tests
that need the real `atom.cr374` are marked and skip when the grid volume is not
mounted, so CI on the runner exercises them and a laptop run still passes.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.model_atom import (  # noqa: E402
    FINE_STRUCTURE, TERM_RESOLVED, ModelAtomError, atom_resolution,
    classify_label, ion_stage_histogram, level_j, read_gerber_atom,
    term_statistical_weight)
from pipeline.nlte_line_identification import (  # noqa: E402
    FLAG_ENERGY, FLAG_LABEL, FLAG_UNMATCHED, LevelResolver,
    LineIdentificationError, air_to_vacuum, identification_provenance,
    identify_lines, ionisation_offset_eV, level_energies_cm, reach_report,
    read_species_lines, render_identification_fields)

# ── fixtures ─────────────────────────────────────────────────────────────────
TERM_ATOM = """Cr
    5.74   52.000
    4     9     3     0
*  EC          G       LABEL                   ION
       0.000      7.0  '    Level    1 =        a7S '   1
    8090.208     25.0  '    Level    2 =        a5D '   1
   24079.555     33.0  '    Level    3 =        a3H '   1
   54580.617      6.0  '    Level    4 =        a6S '   2
"""

FINE_ATOM = """Fe
    7.50   56.000
    3     5     2     0
*  EC          G       LABEL                   ION
       0.000      9.0  '    Level    1 =       a5D4 '   1
     415.933      7.0  '    Level    2 =       a5D3 '   1
   26874.550      7.0  '    Level    3 =      z3D3* '   1
"""

LINELIST = (
    "'Cr I    LTE'\n"
    "  4200.101  3.079 -1.035   -7.800    7.0  2.14E+07 's' 'p'   0.0    1.0"
    " 'Cr I LS:3d5.(4G).4s a5D LS:3d4.(a3F).4s.4p.(3P*) a3H'\n"
    "  4201.000  1.000 -2.000   -7.800    5.0  1.00E+07 's' 'p'   0.0    1.0"
    " 'Cr I LS:3d5.(4G).4s a7S LS:3d4.(a3F).4s.4p.(3P*) zzZ'\n"
    "'Fe I    NLTE'\n"
    "  5000.000  1.000 -1.000   -7.000    9.0  1.00E+07 's' 'p'   0.0    1.0"
    " 'Fe I LS:x a5D4 LS:y z3D3*'\n"
)


@pytest.fixture
def term_atom(tmp_path) -> Path:
    p = tmp_path / "atom.cr4"
    p.write_text(TERM_ATOM)
    return p


@pytest.fixture
def linelist(tmp_path) -> Path:
    p = tmp_path / "ll.list"
    p.write_text(LINELIST)
    return p


# ── the identity that decides term vs fine structure ─────────────────────────
def test_term_weight_is_the_exact_2s1_2l1_identity():
    """sum_J (2J+1) over an LS term telescopes to (2S+1)(2L+1)."""
    assert term_statistical_weight("a5D") == 25      # 1+3+5+7+9,  J=0..4
    assert term_statistical_weight("a3H") == 33      # 9+11+13,    J=4..6
    assert term_statistical_weight("z7F*") == 49     # 1+3+..+13,  J=0..6
    assert term_statistical_weight("a6D") == 30      # 2+4+6+8+10, J=1/2..9/2
    assert term_statistical_weight("b4P") == 12      # 2+4+6
    assert term_statistical_weight("a7S") == 7       # single level


def test_term_weight_declines_labels_the_question_does_not_apply_to():
    for label in ("z3D3*", "3s2.3p2P*", "", "none", "a2[11/2]"):
        assert term_statistical_weight(label) is None


def test_level_j_comes_from_the_label_never_from_g():
    """Deriving J as (g-1)/2 is precisely what fabricates J=12 for a super-level."""
    assert level_j("z3D3*") == 3.0
    assert level_j("a2D1") == 1.0
    assert level_j("a5D") is None            # a bare term has no single J


def test_classify_label_requires_g_to_confirm_the_label():
    assert classify_label("a5D", 25.0) == "term"
    assert classify_label("a5D", 9.0) == "unparsed"     # bare term, g disagrees
    assert classify_label("z3D3*", 7.0) == "j-resolved"
    assert classify_label("z3D3*", 25.0) == "unparsed"  # J label, g disagrees


# ── reading the atom ─────────────────────────────────────────────────────────
def test_element_is_read_from_the_file_not_assumed(term_atom):
    """The Fe hardcode regression: a Cr atom must not report itself as Fe.

    `rya763_level_mapping.read_gerber_atom` sets `species=f"Fe {ion}"` for every
    atom it reads, and RYA-776 imports it for the whole Engine-B deck.
    """
    lv = read_gerber_atom(term_atom)
    assert set(lv["element"]) == {"Cr"}


def test_g_is_preserved_and_not_converted_to_j(term_atom):
    lv = read_gerber_atom(term_atom)
    assert float(lv.loc[lv["term"] == "a5D", "g"].iloc[0]) == 25.0
    assert "J" not in lv.columns, (
        "a J column on a term-resolved atom would be (g-1)/2 = 12 for a5D, a "
        "coordinate no real line carries")


def test_a_short_read_refuses_rather_than_returning_a_partial_atom(tmp_path):
    p = tmp_path / "atom.trunc"
    p.write_text("\n".join(TERM_ATOM.splitlines()[:-1]) + "\n")   # 3 of 4 levels
    with pytest.raises(ModelAtomError, match="declares 4 levels, parsed 3"):
        read_gerber_atom(p)


def test_ion_stages_are_reported(term_atom):
    assert ion_stage_histogram(read_gerber_atom(term_atom)) == {1: 3, 2: 1}


def test_resolution_verdicts(term_atom, tmp_path):
    assert atom_resolution(read_gerber_atom(term_atom)).verdict == TERM_RESOLVED
    p = tmp_path / "atom.fe3"
    p.write_text(FINE_ATOM)
    assert atom_resolution(read_gerber_atom(p)).verdict == FINE_STRUCTURE


# ── energies ─────────────────────────────────────────────────────────────────
def test_air_to_vacuum_is_a_small_positive_shift():
    for w in (4000.0, 6000.0, 9000.0):
        assert 1.0002 < air_to_vacuum(w) / w < 1.0004


def test_upper_level_sits_above_the_lower_by_the_photon_energy():
    lo, up = level_energies_cm(5000.0, 1.0, 0.0)
    assert up > lo
    assert abs((up - lo) - 1.0e8 / air_to_vacuum(5000.0)) < 1e-6


def test_ion_offset_puts_stage_two_on_the_neutral_ground_scale():
    assert ionisation_offset_eV("Cr", 1) == 0.0
    assert abs(ionisation_offset_eV("Cr", 2) - 6.76651) < 1e-9


def test_ion_offset_refuses_rather_than_defaulting_to_zero():
    with pytest.raises(LineIdentificationError, match="no first ionisation potential"):
        ionisation_offset_eV("Xx", 2)
    with pytest.raises(LineIdentificationError, match="not supported"):
        ionisation_offset_eV("Cr", 3)


# ── resolving levels ─────────────────────────────────────────────────────────
def test_resolver_is_built_per_stage_and_cannot_see_the_other(term_atom):
    lv = read_gerber_atom(term_atom)
    r1 = LevelResolver(levels=lv, ion=1)
    assert r1.by_label("a6S")[0] is None          # a6S is a stage-2 level
    assert r1.by_label("a5D")[0]["index"] == 2
    r2 = LevelResolver(levels=lv, ion=2)
    assert r2.by_label("a5D")[0] is None
    assert r2.by_label("a6S")[0]["index"] == 4


def test_resolver_refuses_an_empty_stage(term_atom):
    with pytest.raises(LineIdentificationError, match="no levels for ion stage"):
        LevelResolver(levels=read_gerber_atom(term_atom), ion=3)


def test_ambiguity_refuses_instead_of_taking_the_first(tmp_path):
    dup = TERM_ATOM.replace("Level    3 =        a3H", "Level    3 =        a5D")
    p = tmp_path / "atom.dup"
    p.write_text(dup)
    r = LevelResolver(levels=read_gerber_atom(p), ion=1)
    lvl, why = r.by_label("a5D")
    assert lvl is None and "ambiguous" in why


def test_energy_route_is_OFF_by_default(term_atom):
    """🔴 The finding, pinned. Controlled against the label route where the truth
    is known, the energy route names the WRONG level 18% of the time for Cr I and
    35% for Cr II — this atom's energies are term averages, and neighbouring terms
    sit closer together than the fine structure they average over. It must not run
    unless a caller asks for it and owns the error rate."""
    r = LevelResolver(levels=read_gerber_atom(term_atom), ion=1, energy_tol_cm=50.0)
    lvl, flag, reason = r.resolve("nosuchterm", 8090.0)       # would match by energy
    assert lvl is None and flag == FLAG_UNMATCHED
    assert reason.endswith("energy-route-disabled")


def test_energy_fallback_fires_only_after_the_label_fails_when_enabled(term_atom):
    lv = read_gerber_atom(term_atom)
    r = LevelResolver(levels=lv, ion=1, energy_tol_cm=50.0)
    lvl, flag, reason = r.resolve("a5D", 0.0, allow_energy=True)
    assert flag == FLAG_LABEL and lvl["index"] == 2      # label wins outright
    lvl, flag, reason = r.resolve("nosuchterm", 8090.0, allow_energy=True)
    assert flag == FLAG_ENERGY and lvl["index"] == 2
    assert reason.startswith("label-absent->")
    lvl, flag, _ = r.resolve("nosuchterm", 999999.0, allow_energy=True)
    assert lvl is None and flag == FLAG_UNMATCHED


def test_the_energy_control_can_actually_fail(term_atom):
    """The control must be able to REPORT disagreement, or it is decoration.

    Two terms 284 cm-1 apart — the real a3P/a3H spacing — with a computed energy
    nearer the wrong one is exactly the production failure, so the control has to
    catch it here.
    """
    lv = read_gerber_atom(term_atom)
    # a3H sits at 24079.555; ask the energy route about an energy 12 cm-1 away
    # while the LABEL says a7S (level 1, at 0.0). They must disagree.
    r = LevelResolver(levels=lv, ion=1, energy_tol_cm=50.0)
    truth, _ = r.by_label("a7S")
    guess, _ = r.by_energy(24091.8)
    assert truth is not None and guess is not None
    assert int(truth["index"]) != int(guess["index"])


# ── end to end on the fixture ────────────────────────────────────────────────
def test_reads_only_the_requested_species(linelist):
    cr = read_species_lines(linelist, "Cr", 1)
    assert len(cr) == 2
    assert list(cr["term_low"]) == ["a5D", "a7S"]
    assert list(cr["term_up"]) == ["a3H", "zzZ"]
    assert read_species_lines(linelist, "Cr", 2).empty      # no Cr II block


def test_identify_marks_the_unmatched_line_as_lte(term_atom, linelist):
    lv = read_gerber_atom(term_atom)
    ident = identify_lines(read_species_lines(linelist, "Cr", 1), lv, "Cr", 1)
    assert list(ident["nlte"]) == [True, False]
    good, bad = ident.iloc[0], ident.iloc[1]
    assert (good["level_low"], good["level_up"]) == (2, 3)
    assert (good["flag_low"], good["flag_up"]) == (FLAG_LABEL, FLAG_LABEL)
    assert bad["level_up"] == 0 and bad["label_up"] == "none"
    assert bad["flag_up"] == FLAG_UNMATCHED


def test_reach_report_counts_lte_lines_as_lte_not_as_coverage(term_atom, linelist):
    lv = read_gerber_atom(term_atom)
    rep = reach_report(identify_lines(read_species_lines(linelist, "Cr", 1), lv, "Cr", 1))
    assert rep["n_lines"] == 2
    assert rep["n_nlte"] == 1
    assert rep["reach_pct"] == 50.0
    assert rep["n_lte_despite_nlte_block"] == 1
    assert rep["upper_by_flag"][FLAG_UNMATCHED] == 1


def test_emitted_fields_carry_all_six_in_order(term_atom, linelist):
    lv = read_gerber_atom(term_atom)
    ident = identify_lines(read_species_lines(linelist, "Cr", 1), lv, "Cr", 1)
    assert render_identification_fields(ident.iloc[0]) == "  2 3  'a5D' 'a3H'  'c' 'c'"
    assert render_identification_fields(ident.iloc[1]).endswith("'none'  'c' 'x'")


# ── 🔴 the trap ──────────────────────────────────────────────────────────────
def test_provenance_records_whether_the_energy_route_ran(term_atom):
    res = atom_resolution(read_gerber_atom(term_atom))
    off = identification_provenance("Cr", "atom.cr4", res, 50.0)
    assert off["energy_fallback_used"] is False
    assert off["energy_tolerance_cm-1"] is None      # no tolerance to report
    assert "18%" in off["energy_fallback_note"]
    on = identification_provenance("Cr", "atom.cr4", res, 50.0, energy_fallback=True)
    assert on["energy_fallback_used"] is True
    assert on["energy_tolerance_cm-1"] == 50.0
    assert "ENABLED" in on["matching"]


def test_multiplicity_matching_is_disabled_and_says_why(term_atom):
    """Requiring the line's 2J+1 to equal a super-level's g rejects ~95% of
    CORRECT matches. It must stay off, and the reason must travel with the
    product so it is not switched on as a tightened cross-check."""
    res = atom_resolution(read_gerber_atom(term_atom))
    prov = identification_provenance("Cr", "atom.cr4", res, 50.0)
    assert prov["multiplicity_matching_used"] is False
    assert prov["wavelength_matching_used"] is False
    assert "one fine-structure component" in prov["multiplicity_note"]


def test_term_resolved_product_declares_the_approximation(term_atom, tmp_path):
    res = atom_resolution(read_gerber_atom(term_atom))
    assert "TERM-SHARED" in identification_provenance("Cr", "a", res, 50.0)["approximation"]
    p = tmp_path / "atom.fe3"
    p.write_text(FINE_ATOM)
    fine = atom_resolution(read_gerber_atom(p))
    assert "approximation" not in identification_provenance("Fe", "a", fine, 50.0)


# ── against the real atom, when it is mounted ────────────────────────────────
def _real_atom() -> Path | None:
    try:
        from config.constants import codex_path
        p = codex_path('grids.gerber_ts') / "atom.cr374"
        return p if p.exists() else None
    except Exception:
        return None


needs_atom = pytest.mark.skipif(_real_atom() is None,
                                reason="atom.cr374 not mounted (grid volume offline)")


@needs_atom
def test_the_real_cr_atom_is_term_resolved_and_that_is_why_j_keying_fails():
    lv = read_gerber_atom(_real_atom())
    res = atom_resolution(lv)
    assert res.is_term_resolved, res.describe()
    assert res.n_term > 250
    # the exact level that breaks (J, energy) keying
    a5d = lv[(lv["term"] == "a5D") & (lv["ion"] == 1)].iloc[0]
    assert a5d["g"] == 25.0
    assert (a5d["g"] - 1) / 2 == 12.0, "the fabricated J no Cr line can ever carry"


@needs_atom
def test_the_real_cr_stage_split_is_what_the_driver_asserts():
    assert ion_stage_histogram(read_gerber_atom(_real_atom())) == {1: 148, 2: 225, 3: 1}
