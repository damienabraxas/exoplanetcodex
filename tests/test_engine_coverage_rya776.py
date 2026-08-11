"""
RYA-776 -- the engine x wavelength coverage reference.

What is actually at risk here, and therefore what these test:

  * THE THREE (FOUR) STATES STAYING APART. The whole ticket exists because SERVED,
    REACHABLE-NOT-EXTRACTED and UNCOVERED had been collapsing into one blurry "no
    coverage". The state machine is tested directly, including the case that forced a
    fourth state: reach that is not locally decidable must NOT read as absence.
  * THE REDUCTION NEVER MANUFACTURING AN ABSENCE. An element can hold two grids per
    engine. Reducing them to one answer must never let UNCOVERED win over a grid whose
    reach is unknown -- that would invent the exact false verdict this table prevents.
  * DETERMINISM (RYA-768). The row order must be a TOTAL order over the full key,
    grid_id included, or the artifact never byte-diffs clean.
  * THE READER REFUSING BAD INPUT. An unrecognised state or a missing file must raise,
    never degrade to "uncovered" -- an unreadable reference and a real gap must not be
    indistinguishable (the RYA-708 rule, applied to the engine half).

These are hermetic: they build tables in tmp_path rather than depending on the generated
artifact, which is produced on Sirius and is not present in a Mac checkout.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline import coverage  # noqa: E402
from pipeline.coverage import (  # noqa: E402
    REACH_UNKNOWN, REACHABLE_NOT_EXTRACTED, SERVED, UNCOVERED, CoverageError,
    EngineCoverage, engine_reach, engine_summary, load_engine_coverage,
)

HEADER = ("element,ion,engine,grid_id,band,band_lo_A,band_hi_A,state,n_lines_served,"
          "n_lines_reachable,n_lines_catalogued,level_asset,grid_asset,note\n")


def _row(element="Fe", ion="I", engine="A", grid_id="G", band="VIS",
         lo=3800.0, hi=6910.0, state=SERVED, served=5, reach=5, cat=10,
         level="label_Fe.txt", grid="G.csv", note=""):
    return (f"{element},{ion},{engine},{grid_id},{band},{lo},{hi},{state},"
            f"{served},{reach},{cat},{level},{grid},{note}\n")


def _table(tmp_path: Path, *rows: str) -> Path:
    p = tmp_path / "engine_coverage.csv"
    p.write_text("# generated\n" + HEADER + "".join(rows), encoding="utf-8")
    return p


# ── the state machine ────────────────────────────────────────────────────────

def test_classify_keeps_the_four_states_apart():
    """The generator's state machine, exercised on the four distinct situations."""
    gen = pytest.importorskip("scripts.generate_engine_coverage_rya776",
                              reason="generator imports pipeline.nlte_corrections")
    # served: the extract resolves lines here
    assert gen.classify(3, 0, 10, "label_Fe.txt", "")[0] == SERVED
    # reachable but not extracted: the levels cover it, no extract does
    assert gen.classify(0, 7, 10, "label_Fe.txt", "")[0] == REACHABLE_NOT_EXTRACTED
    # uncovered: catalogued lines exist, the table answers for some endpoints (so the
    # key works) and carries neither level of any line
    assert gen.classify(0, 0, 10, "label_Fe.txt", "", 4)[0] == UNCOVERED
    # NOT decidable -- no local level asset. This is Fe/MPIA, and calling it UNCOVERED
    # would be the false "no coverage" the ticket exists to end.
    assert gen.classify(0, 0, 10, "", "web-service supplier", 0)[0] == REACH_UNKNOWN
    # NOT decidable -- no catalogued line in band. A zero here measures the LINELIST's
    # span, not the grid's (the 9199.9 A wall is GES; atom.fe607a reaches 20000 A).
    state, note = gen.classify(0, 0, 0, "atom.fe607a", "", 0)
    assert state == REACH_UNKNOWN
    assert "LINELIST" in note


def test_uncovered_requires_a_decidable_reach():
    """UNCOVERED must never be reachable without a level asset AND a live denominator."""
    gen = pytest.importorskip("scripts.generate_engine_coverage_rya776")
    for level, cat in (("", 10), ("", 0), ("atom.fe607a", 0)):
        assert gen.classify(0, 0, cat, level, "why", 0)[0] != UNCOVERED


def test_a_key_that_addresses_nothing_is_not_an_absence():
    """The Fe II super-level defect, caught on the first real run.

    The Gerber atoms pack higher ionisation stages into super-levels whose statistical
    weight is the sum of the merged levels', so J = (G-1)/2 is not a J (Fe II reads
    14.5, 27.5) and not one endpoint resolves. Writing UNCOVERED there would have
    claimed we cannot model Fe II in the optical -- the ionisation arbiter, measured in
    production. Zero partial matches means the key failed, not that the physics is gone.
    """
    gen = pytest.importorskip("scripts.generate_engine_coverage_rya776")
    state, note = gen.classify(0, 0, 8870, "atom.fe607a", "", 0)
    assert state == REACH_UNKNOWN
    assert "super-level" in note
    # ...but a table that answers for SOME endpoints and lacks the others IS reporting a
    # real absence, and must still say so. This is RYA-763's Fe I IR result.
    state, note = gen.classify(0, 0, 5871, "atom.fe607a", "", 4189)
    assert state == UNCOVERED
    assert "resolve one endpoint" in note


def test_higher_stages_are_rebased_to_their_own_ground_state():
    """The Mn II false absence: a Gerber atom counts EVERY stage's energies from the
    NEUTRAL ground state, GES counts each ion's from its own. Compared raw, not one of
    3386 Mn II endpoints resolved and the generator wrote UNCOVERED over 1693 catalogued
    optical lines; rebased, 171 resolve both."""
    import pandas as pd
    gen = pytest.importorskip("scripts.generate_engine_coverage_rya776")
    atom = pd.DataFrame({"ion": [1, 1, 2, 2], "J": [0.5, 1.5, 0.5, 1.5],
                         "energy_eV": [0.0, 2.1, 7.434, 9.534]})
    neutral, off0 = gen._ion_filter(atom, "I")
    assert off0 == 0.0 and list(neutral.energy_eV) == [0.0, 2.1]
    ionised, off = gen._ion_filter(atom, "II")
    assert off == pytest.approx(7.434)
    assert list(ionised.energy_eV) == pytest.approx([0.0, 2.1])


def test_an_implausible_offset_is_refused_not_applied():
    """The rebase's one failure mode is a stage whose GROUND STATE is missing -- then the
    minimum is an excited level and shifting by it invents a coordinate. Bounds-checked
    against the range real ionisation potentials occupy, and refused outside it."""
    import pandas as pd
    gen = pytest.importorskip("scripts.generate_engine_coverage_rya776")
    atom = pd.DataFrame({"ion": [2, 2], "J": [0.5, 1.5], "energy_eV": [0.4, 1.4]})
    sub, off = gen._ion_filter(atom, "II")
    assert off == 0.0                       # 0.4 eV is no ionisation potential
    assert list(sub.energy_eV) == [0.4, 1.4]


def test_ion_code_does_not_fold_every_stage_onto_two():
    """`1 if ion == 'I' else 2` mapped Fe III onto Fe 2, so the Fe III rows carried
    Fe II's catalogued count. Caught on the first real run."""
    gen = pytest.importorskip("scripts.generate_engine_coverage_rya776")
    assert gen.ges_ion_code("I") == 1
    assert gen.ges_ion_code("II") == 2
    assert gen.ges_ion_code("III") == 3
    assert gen.ges_ion_code("III") != gen.ges_ion_code("II")


# ── the reader helper ────────────────────────────────────────────────────────

def test_engine_reach_answers_the_ticket_question(tmp_path):
    """'Fe Engine-A at 8000 A' -- the lookup that used to be a re-derivation."""
    p = _table(
        tmp_path,
        _row(band="VIS", lo=3800.0, hi=6910.0, state=SERVED, served=81),
        _row(band="red-optical", lo=6910.0, hi=10000.0, state=REACH_UNKNOWN,
             served=0, reach=0, cat=462, level="",
             note="web-service supplier; not locally decidable"),
    )
    tab = load_engine_coverage(p)
    assert engine_reach("Fe", "I", "A", 5000.0, tab).state == SERVED
    ans = engine_reach("Fe", "I", "A", 8000.0, tab)
    assert ans.state == REACH_UNKNOWN
    assert ans.band == "red-optical"
    # The load-bearing assertion: this is NOT a data gap, and must not read as one.
    assert not ans.is_data_gap
    assert "NOT a claim of absence" in ans.why()


def test_band_edges_belong_to_exactly_one_band(tmp_path):
    """Bands are half-open, so a boundary wavelength lands in one band, never two."""
    p = _table(tmp_path,
               _row(band="VIS", lo=3800.0, hi=6910.0),
               _row(band="red-optical", lo=6910.0, hi=10000.0, state=UNCOVERED,
                    served=0, reach=0))
    tab = load_engine_coverage(p)
    assert len(engine_reach("Fe", "I", "A", 6910.0, tab).rows) == 1
    assert engine_reach("Fe", "I", "A", 6910.0, tab).band == "red-optical"


def test_missing_row_is_not_a_coverage_verdict(tmp_path):
    """A species/engine pair that was never generated is an UNBUILT reference."""
    tab = load_engine_coverage(_table(tmp_path, _row()))
    ans = engine_reach("V", "I", "A", 5000.0, tab)
    assert ans.state == REACH_UNKNOWN
    assert not ans.is_data_gap
    assert "NO ROW" in ans.why()


def test_reduction_never_invents_an_absence(tmp_path):
    """Two grids for one engine: UNCOVERED must not outrank an undecided one.

    Mg and Si each hold both a Bergemann/MPIA and an Amarsi/PySME extract, and they do
    not have the same reach. Asserting 'genuinely absent' requires EVERY grid to have
    been decidable, so one REACH-UNKNOWN is enough to withhold the verdict.
    """
    p = _table(
        tmp_path,
        _row(element="Mg", grid_id="Mg_Amarsi2020_PySME", state=UNCOVERED,
             served=0, reach=0, cat=12),
        _row(element="Mg", grid_id="Mg_Bergemann_MPIA", state=REACH_UNKNOWN,
             served=0, reach=0, cat=12, level=""),
    )
    ans = engine_reach("Mg", "I", "A", 5000.0, load_engine_coverage(p))
    assert ans.state == REACH_UNKNOWN
    assert not ans.is_data_gap
    assert len(ans.rows) == 2          # the reduction never hides the per-grid detail


def test_served_outranks_everything(tmp_path):
    p = _table(tmp_path,
               _row(grid_id="A1", state=UNCOVERED, served=0, reach=0),
               _row(grid_id="A2", state=SERVED, served=4))
    ans = engine_reach("Fe", "I", "A", 5000.0, load_engine_coverage(p))
    assert ans.state == SERVED and ans.is_data_gap is False


def test_a_real_gap_is_still_reported_as_one(tmp_path):
    """The refusal to over-claim absence must not become a refusal to report it."""
    p = _table(tmp_path, _row(element="V", state=UNCOVERED, served=0, reach=0, cat=31))
    ans = engine_reach("V", "I", "A", 5000.0, load_engine_coverage(p))
    assert ans.state == UNCOVERED and ans.is_data_gap
    assert "real modelling gap" in ans.why()


# ── refusing bad input ───────────────────────────────────────────────────────

def test_unknown_state_raises_rather_than_degrading(tmp_path):
    p = tmp_path / "engine_coverage.csv"
    p.write_text(HEADER + _row(state="PROBABLY-FINE"), encoding="utf-8")
    with pytest.raises(CoverageError, match="not one of"):
        load_engine_coverage(p)


def test_unknown_engine_raises(tmp_path):
    p = tmp_path / "engine_coverage.csv"
    p.write_text(HEADER + _row(engine="C"), encoding="utf-8")
    with pytest.raises(CoverageError):
        load_engine_coverage(p)
    with pytest.raises(CoverageError):
        engine_reach("Fe", "I", "C", 5000.0, [])


def test_missing_reference_raises_and_says_how_to_build_it(tmp_path):
    with pytest.raises(CoverageError, match="GENERATED"):
        load_engine_coverage(tmp_path / "nope.csv")


def test_empty_reference_is_not_no_coverage(tmp_path):
    p = tmp_path / "engine_coverage.csv"
    p.write_text(HEADER, encoding="utf-8")
    with pytest.raises(CoverageError, match="not 'no coverage'"):
        load_engine_coverage(p)


# ── determinism (RYA-768) ────────────────────────────────────────────────────

def test_sort_key_is_a_total_order_over_the_full_key():
    """grid_id is IN the key. RYA-768's defect was a key that omitted a field and so
    left rows ordered by whatever the filesystem happened to return."""
    gen = pytest.importorskip("scripts.generate_engine_coverage_rya776")
    rows = [
        dict(element="Mg", ion="I", engine="A", grid_id="Mg_Bergemann_MPIA",
             band="VIS", band_lo_A=3800.0),
        dict(element="Mg", ion="I", engine="A", grid_id="Mg_Amarsi2020_PySME",
             band="VIS", band_lo_A=3800.0),
    ]
    key = lambda r: (r["element"], r["ion"], r["engine"], r["grid_id"], r["band_lo_A"])  # noqa: E731
    assert len({key(r) for r in rows}) == len(rows), "key must separate these two rows"
    assert [r["grid_id"] for r in sorted(rows, key=key)] == \
        ["Mg_Amarsi2020_PySME", "Mg_Bergemann_MPIA"]
    assert sorted(rows, key=key) == sorted(list(reversed(rows)), key=key)


def test_render_is_byte_stable_across_calls():
    gen = pytest.importorskip("scripts.generate_engine_coverage_rya776")
    rows = [dict(element="Fe", ion="I", engine="A", grid_id="G", band="VIS",
                 band_lo_A=3800.0, band_hi_A=6910.0, state=SERVED, n_lines_served=1,
                 n_lines_reachable=1, n_lines_catalogued=2, level_asset="l",
                 grid_asset="g", note="")]
    assert gen.render(rows) == gen.render(rows)


def test_ion_normalisation_joins_the_two_conventions():
    """Al_Amarsi2020_PySME.csv writes ion `1` where every other extract writes `I`.
    Unnormalised, Al would split into two species that never join -- a silently
    half-empty answer, which is worse than a loud one."""
    gen = pytest.importorskip("scripts.generate_engine_coverage_rya776")
    assert gen._norm_ion("1") == gen._norm_ion("I") == "I"
    assert gen._norm_ion(2) == gen._norm_ion("II") == "II"
    assert gen._norm_ion("1.0") == "I"
    with pytest.raises(SystemExit):
        gen._norm_ion("neutral")


# ── the tracker join ─────────────────────────────────────────────────────────

def test_engine_summary_marks_reach_only_bands(tmp_path):
    """The compact tracker cell: served bands plain, reach-only bands with `?`."""
    p = _table(
        tmp_path,
        _row(engine="A", band="VIS", lo=3800.0, hi=6910.0, state=SERVED, served=81),
        _row(engine="A", band="red-optical", lo=6910.0, hi=10000.0,
             state=REACH_UNKNOWN, served=0, reach=0, level=""),
        _row(engine="B", grid_id="atom.fe607a", band="VIS", lo=3800.0, hi=6910.0,
             state=SERVED, served=79),
        _row(engine="B", grid_id="atom.fe607a", band="red-optical", lo=6910.0,
             hi=10000.0, state=REACHABLE_NOT_EXTRACTED, served=0, reach=1620),
    )
    s = engine_summary("Fe", "I", load_engine_coverage(p))
    assert s == "A:VIS · B:VIS,red-optical?"


def test_tracker_join_degrades_visibly_not_silently():
    """A tracker regenerated where the reference was never built must SAY so.

    The reach table is generated on Sirius; a Mac regeneration must not fail, but it must
    also not emit a blank that reads as "this engine reaches nothing".
    """
    gen = pytest.importorskip("scripts.generate_element_status_tracker_rya654")
    cell = gen._engine_reach_cell(None, "Fe", "I")
    assert "not generated" in cell and "Sirius" in cell


def test_engine_reach_is_in_the_generated_columns():
    """It is derived from a generated sibling, so hand-editing it must be a build break
    like every other generated column -- not an analyst field."""
    gen = pytest.importorskip("scripts.generate_element_status_tracker_rya654")
    assert "engine_reach" in gen.COLUMNS
    assert "engine_reach" in gen.GENERATED_COLUMNS


def test_instrument_half_is_untouched():
    """RYA-776 COMPLEMENTS instrument coverage; it must not have replaced it."""
    for name in ("coverage_at", "load_registry", "instruments_for", "verify"):
        assert hasattr(coverage, name)
    assert coverage.ENGINE_COVERAGE.name == "engine_coverage.csv"
    assert coverage.CATALOG.name == "instrument_catalog.csv"
