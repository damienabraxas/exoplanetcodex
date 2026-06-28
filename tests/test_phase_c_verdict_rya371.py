"""
tests/test_phase_c_verdict_rya371.py
====================================
RYA-371 Phase C — guard the solar 27-element verdict assembler.
VALIDATE-DON'T-TUNE: the verdict is a CLASSIFICATION of the measured baseline
against Asplund 2021; the assembler must never mutate a measured value toward the
anchor, must source C/N/O from the Phase A synthesis (not the uncurated EW C),
must keep Fe the validated PASS leg, and must flag P as the DATA-GAP.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

import phase_c_verdict_rya371 as P  # noqa: E402


def _rows():
    ab, ew, phase_a = P._load()
    return {r['element']: r for r in P.build_verdicts(ab, ew, phase_a)}, phase_a


def test_fe_is_the_validated_pass_leg():
    rows, _ = _rows()
    fe = rows['Fe']
    assert fe['verdict'] == 'PASS'
    assert fe['n_lines'] >= 60                      # the deep Fe I pool
    # +0.05 vs Asplund is the documented scale offset, well inside the band
    assert abs(fe['delta_vs_asplund']) <= P.TOL_PASS


def test_oxygen_passes_cross_arm_not_tuned():
    rows, _ = _rows()
    o = rows['O']
    assert o['verdict'] == 'PASS'
    # measured O reconciles near (NOT exactly equal to) the Asplund anchor — proof
    # it is a measurement carried with its real +0.05 offset, not snapped to 8.69
    assert o['A_measured'] != 8.69
    assert abs(o['delta_vs_asplund']) <= P.TOL_PASS


def test_cno_sourced_from_phase_a_not_ew_artifact():
    rows, phase_a = _rows()
    # The EW path yields C = 10.26 (uncurated artifact). The verdict must report the
    # Phase A synthesis value (~8.49), never the EW garbage.
    assert rows['C']['A_measured'] == phase_a['cross_arm']['C']['primary_mean']
    assert rows['C']['A_measured'] < 9.0
    assert rows['C']['channel'].startswith('synthesis')


def test_nitrogen_is_nlte_owed_with_named_grid():
    rows, _ = _rows()
    n = rows['N']
    assert n['verdict'] == 'NLTE-OWED'
    assert 'N I' in n['owed'] and 'RYA-369' in n['owed']   # the owed grid is named


def test_phosphorus_is_the_data_gap():
    rows, _ = _rows()
    p = rows['P']
    assert p['verdict'] == 'DATA-GAP'
    assert 'HST' in p['owed'] or 'RYA-119' in p['owed']


def test_measured_metals_are_curation_owed_not_silently_passed():
    rows, _ = _rows()
    # RYA-456: the RYA-395/398 graded curation is now wired into the run. The
    # gf-limited metals that carry an independent-gf line (Ti/Cr/Ni/Si) now produce
    # A(X) AT THE DOCUMENTED FLOOR and must stay CURATION-OWED — never PASS.
    for el in ('Ti', 'Cr', 'Ni', 'Si'):
        assert rows[el]['verdict'] == 'CURATION-OWED', el
        assert rows[el]['A_measured'] is not None          # produced at the floor, off-anchor
    # Mg's pool has no independent-gf line → the graded firewall culls it to nothing;
    # still CURATION-OWED (gf-data-limited), but with no produced value.
    assert rows['Mg']['verdict'] == 'CURATION-OWED'
    # The Cr canary: the graded pool lands Cr ~+0.40 vs Asplund (RYA-398), NOT at the
    # anchor. A Cr PASS / Cr≈0 would mean tuning crept into the wiring.
    cr = rows['Cr']
    assert cr['verdict'] != 'PASS'
    assert cr['delta_vs_asplund'] is not None and cr['delta_vs_asplund'] > 0.25, cr


def test_counts_cover_all_metals_once():
    rows, _ = _rows()
    assert len(rows) == 26                                   # Asplund metals (no H/He)
    counts = {}
    for r in rows.values():
        counts[r['verdict']] = counts.get(r['verdict'], 0) + 1
    assert sum(counts.values()) == 26
    assert counts['PASS'] == 3                               # O, C, Fe
    assert set(v for v in counts) <= {'PASS', 'NLTE-OWED', 'CURATION-OWED', 'DATA-GAP'}


def test_no_value_mutated_toward_anchor():
    # Cardinal rule: delta = A_measured - asplund must be the raw arithmetic, never
    # shrunk. Recompute independently and assert identity.
    rows, _ = _rows()
    for r in rows.values():
        if r['A_measured'] is not None:
            assert r['delta_vs_asplund'] == round(r['A_measured'] - r['asplund2021'], 3), r['element']


def test_differential_backbone_registers_optical_instruments():
    bb = P.differential_backbone()
    insts = bb['instruments']
    # the four optical/IR legs + the RYA-459/460 reference legs (KITT_PEAK measured,
    # UV_COMPOSITE cited)
    assert {'HARPS', 'ESPRESSO', 'UVES', 'CRIRES+'} <= set(insts)
    assert {'KITT_PEAK', 'UV_COMPOSITE'} <= set(insts)      # RYA-460 registered the KP leg
    assert 'PARKED' in insts['CRIRES+']['status']           # Phase B telluric gate
    assert insts['ESPRESSO']['R'] == 140000
