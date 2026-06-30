"""
tests/test_scope_cno_redo_rya486.py
===================================
RYA-486 is a SCOPE decision (docs/design/scope_cno_redo_rya486.md): the CNO redo
after RYA-485 is full-CNO solar + surgical-CNO Procyon, and the 27-element Procyon
run (RYA-404) is NOT a CNO side-effect. That scope rests on four facts about the
code. These tests pin those facts so the scope is flagged if the code drifts out
from under it. CI-safe: no external grid/spectra needed.
"""
from pathlib import Path

from config.constants import get_acceptance_profile, STAR_SPECTRAL_TYPE
from pipeline import nlte_cno

ROOT = Path(__file__).resolve().parents[1]


def test_decision_A_cno_is_coupled_so_o_alone_is_impossible():
    # The synthesis engine must declare the molecular coupling that makes a solar
    # O-only redo impossible (O redo => CNO redo).
    doc = (ROOT / 'pipeline' / 'cno_synthesis.py').read_text().lower()
    assert 'molecular equilibrium' in doc
    assert 'co' in doc and 'cn' in doc and 'ch' in doc


def test_decision_A_ir_unlock_oi_844_926_in_amarsi_grid():
    # The IR unlock is real: O I 844 and 926 nm are grid-routable multiplets.
    spans = nlte_cno._OI_MULTIPLET_SPANS
    assert '844nm' in spans and '926nm' in spans
    lo844, hi844 = spans['844nm']
    lo926, hi926 = spans['926nm']
    assert lo844 <= 8446.4 <= hi844            # O I 8446 routes to 844nm
    assert lo926 <= 9263.0 <= hi926            # O I 9263 routes to 926nm


def test_decision_C_metals_ew_path_is_independent_of_cno_synthesis():
    # RYA-404 metals come from the EW path, which must not import the CNO engine —
    # so re-deriving CNO cannot perturb a banked metal.
    src = (ROOT / 'pipeline' / 'abundances_derive.py').read_text()
    assert 'cno_synthesis' not in src
    assert 'run_cno' not in src


def test_decision_C_404_is_its_own_fstar_gate_not_a_cno_consequence():
    # The F-star (Procyon / RYA-404) acceptance profile is its own NLTE-unavailable
    # gate, unrelated to the CNO denominators.
    assert STAR_SPECTRAL_TYPE['procyon'] == 'F'
    f = get_acceptance_profile('F')
    assert f['nlte_available'] is False        # Fe I NLTE grid runs out > 6500 K
    assert 'fe1_scatter_max' in f              # F-star scatter floor (RYA-281)
