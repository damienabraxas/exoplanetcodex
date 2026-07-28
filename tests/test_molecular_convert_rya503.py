"""
tests/test_molecular_convert_rya503.py
======================================
RYA-503 — the molecular converter + the acquired OH/NH/CH mid-IR lists.

Offline throughout: the converter physics is exercised on a tiny synthetic ExoMol
pair (no network); the CO round-trip acceptance is validated live against the ExoMol
source at build time (reported on the ticket) and its recipe is pinned in the manifest.
The vendored mid-IR lists + manifest + the [molecular] guard are checked from the repo.
"""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from scripts import molecular_linelist_convert as mc      # noqa: E402
from pipeline import molecular_lists as ml                 # noqa: E402
import scripts.check_stewardship as sc                     # noqa: E402

_ACQUIRED = {'OH_MYTHOS_midIR', 'NH_kNigHt_midIR', 'CH_MoLLIST_midIR'}


# ── converter physics (offline; synthetic ExoMol pair) ────────────────────────
def _write_synthetic(tmp, states_rows, trans_rows):
    (tmp / 's.states').write_text("\n".join(states_rows) + "\n")
    (tmp / 't.trans').write_text("\n".join(trans_rows) + "\n")


def test_gf_lambda_chi_formulas(tmp_path):
    # two ground-state states: v0 J8 (E=138.39) and v11 J9 (E=22286.66) → the vendored
    # CO v11-0 J9-8 line: λ 4515.03 (vac), gu 19, A 8.35e-10, loggf -16.314, χ 0.01716
    _write_synthetic(
        tmp_path,
        ["           1   138.390000      17       8   0 e X",
         "           2 22286.663200      19       9  11 e X"],
        ["           2           1 8.3500E-10"])
    edef = mc.ExomolDef('', '', ['v', 'e/f', 'ElecState'], False, False, False)
    states = mc.parse_exomol_states(tmp_path / 's.states', edef)
    lines = mc.convert_exomol(states, tmp_path / 't.trans', tag='X')
    assert len(lines) == 1
    r = lines[0]
    assert abs(r.lam_A - 4515.03) < 0.02
    assert r.g_u == 19.0 and abs(r.chi_low_eV - 0.01716) < 1e-4
    assert abs(r.loggf - (-16.314)) < 0.01
    assert r.Ju == 9 and r.Jl == 8 and r.desc == 'v11-0_J9-8_X'


def test_upper_lower_assigned_by_energy(tmp_path):
    # trans lists (lower_id, upper_id) reversed → converter must still pick higher-E as upper
    _write_synthetic(
        tmp_path,
        ["           1   138.390000      17       8   0 e X",
         "           2 22286.663200      19       9  11 e X"],
        ["           1           2 8.3500E-10"])         # reversed order
    edef = mc.ExomolDef('', '', ['v', 'e/f', 'ElecState'], False, False, False)
    lines = mc.convert_exomol(mc.parse_exomol_states(tmp_path / 's.states', edef),
                              tmp_path / 't.trans', tag='X')
    assert lines[0].g_u == 19.0 and lines[0].loggf < 0        # upper = the v11 state


def test_ground_only_filter(tmp_path):
    # an A-state (electronic) upper must be dropped by ground_only; X-X kept
    _write_synthetic(
        tmp_path,
        ["           1   138.390000      17       8   0 e X",
         "           2 22286.663200      19       9  11 e X",
         "           3 32500.000000      19       9   0 e A"],
        ["           2           1 8.3500E-10",           # X-X kept
         "           3           1 1.0000E+06"])          # A-X dropped
    edef = mc.ExomolDef('', '', ['v', 'e/f', 'ElecState'], False, False, False)
    states = mc.parse_exomol_states(tmp_path / 's.states', edef)
    alll = mc.convert_exomol(states, tmp_path / 't.trans', tag='X')
    gnd = mc.convert_exomol(states, tmp_path / 't.trans', tag='X', ground_only=True)
    assert len(alll) == 2 and len(gnd) == 1


def test_nu_min_filter(tmp_path):
    _write_synthetic(
        tmp_path,
        ["           1     0.000000       1       0   0 e X",
         "           2     5.000000       3       1   0 e X"],    # ν=5 cm⁻¹ (far-IR)
        ["           2           1 1.0E-05"])
    edef = mc.ExomolDef('', '', ['v', 'e/f', 'ElecState'], False, False, False)
    states = mc.parse_exomol_states(tmp_path / 's.states', edef)
    with pytest.raises(mc.MolecularConvertError):     # everything below ν_min → 0 lines → loud
        mc.convert_exomol(states, tmp_path / 't.trans', tag='X', nu_min=1000.0)


def test_format_parse_roundtrip(tmp_path):
    _write_synthetic(
        tmp_path,
        ["           1   138.390000      17       8   0 e X",
         "           2 22286.663200      19       9  11 e X"],
        ["           2           1 8.3500E-10"])
    edef = mc.ExomolDef('', '', ['v', 'e/f', 'ElecState'], False, False, False)
    lines = mc.convert_exomol(mc.parse_exomol_states(tmp_path / 's.states', edef),
                              tmp_path / 't.trans', tag='X')
    out = tmp_path / 'o.bsyn'
    out.write_text(mc.format_bsyn(lines, '0608.012016', 'test'))
    back = mc.parse_bsyn(out)
    assert 'v11-0_J9-8_X' in back
    r = back['v11-0_J9-8_X']
    assert abs(r.lam_A - lines[0].lam_A) < 1e-3 and abs(r.loggf - lines[0].loggf) < 1e-3


def test_selftest_passes():
    assert mc.selftest() is True


def test_constants_are_cited_values():
    assert abs(mc.CM_PER_EV - 8065.543937) < 1e-6      # CODATA
    assert abs(mc.GF_CONST - 1.49919e-16) < 1e-22


# ── the acquired mid-IR lists (from the repo; offline) ────────────────────────
@pytest.fixture(scope='module')
def manifest():
    return ml.load_manifest()


def test_co_conversion_recipe_pinned(manifest):
    # Phase-1 gap-closer: the CO entry now carries the reproducible recipe + round-trip
    conv = manifest['molecules']['CO'].get('conversion')
    assert conv and 'PASS' in conv['roundtrip']
    assert conv['version'] == '20170101' and 'ν ≥ 1000' in conv['recipe']
    assert conv['doi'].startswith('10.1088')


@pytest.mark.parametrize('key', sorted(_ACQUIRED))
def test_acquired_entry_complete(manifest, key):
    e = manifest['molecules'][key]
    assert e['origin'] == 'acquired (RYA-503)' and e['in_ispec'] is False
    assert e['coverage_gate'] == 'midir'
    assert e['doi'] and e['version'] and e['source'] and e['distribution']
    wc = e['wavelength_coverage']
    assert wc['regime'] == 'mid-IR-rovibrational' and wc['min_A'] and wc['max_A']
    mw = e['midir_window']
    assert mw['count'] > 0 and 'PRESENT' in mw['verdict']


@pytest.mark.parametrize('key', sorted(_ACQUIRED))
def test_vendored_file_matches_recorded_counts(manifest, key):
    e = manifest['molecules'][key]
    sub = ml.VENDORED_DIR / e['vendored_subdir']
    live = sum(ml.count_bsyn_lines(sub / f) for f in e['files'])
    assert live == e['line_count']
    # recompute the mid-IR window count from the vendored bytes → matches the manifest
    clo, chi = e['midir_window']['range_cm-1']
    n = 0
    for f in e['files']:
        for w in ml.bsyn_wavelengths(sub / f):
            nu = 1.0e8 / w
            if clo <= nu <= chi:
                n += 1
    assert n == e['midir_window']['count']


def test_midir_present_where_electronic_was_empty(manifest):
    # the point of the acquisition: non-zero mid-IR rows, vs the 0 the held .bsyn carry
    for mol, key in (('OH', 'OH_MYTHOS_midIR'), ('NH', 'NH_kNigHt_midIR'), ('CH', 'CH_MoLLIST_midIR')):
        held = manifest['molecules'][mol].get('midir_window', {}).get('count')
        assert held == 0                                   # RYA-360 measured absence
        assert manifest['molecules'][key]['midir_window']['count'] > 0


# ── the [molecular] guard covers the acquired lists ───────────────────────────
def test_molecular_guard_passes_with_acquisitions():
    assert sc.check_molecular_lists() == []


def test_negative_missing_midir_list_fails_loud(monkeypatch, tmp_path):
    fake = {'molecules': {'OH_MYTHOS_midIR': {
        'vendored_subdir': 'OH', 'files': ['16O-1H__MYTHOS_rovib.bsyn'], 'line_count': 72762,
        'source': 's', 'distribution': 'd', 'origin': 'acquired (RYA-503)', 'in_ispec': False,
        'coverage_gate': 'midir', 'wavelength_coverage': {'min_A': 3133, 'max_A': 99939, 'regime': 'mid-IR-rovibrational'},
        'midir_window': {'label': 'OH 1-0', 'range_cm-1': [2600, 3600], 'count': 5779},
    }}}
    monkeypatch.setattr(ml, 'load_manifest', lambda: fake)
    monkeypatch.setattr(ml, 'VENDORED_DIR', tmp_path)                # file absent
    monkeypatch.setattr(ml, 'ISPEC_MOLECULES_DIR', tmp_path / 'no_ispec')
    viol = sc.check_molecular_lists()
    assert any(v.invariant == 'molecular' and not v.tracked and 'present' in v.quantity
               for v in viol)


def test_negative_empty_midir_window_fails_loud(monkeypatch, tmp_path):
    (tmp_path / 'OH').mkdir()
    (tmp_path / 'OH' / 'x.bsyn').write_text("'0108.000016' 1 1\n'lbl'\n 5000.0 0.0 -1.0 0 3 1e0 'X' 'X' 1 0 'X' 'X' 0 0 'l'\n")
    fake = {'molecules': {'OH_MYTHOS_midIR': {
        'vendored_subdir': 'OH', 'files': ['x.bsyn'], 'line_count': 1,
        'source': 's', 'distribution': 'd', 'origin': 'acquired (RYA-503)', 'in_ispec': False,
        'coverage_gate': 'midir', 'wavelength_coverage': {'min_A': 5000, 'max_A': 5000, 'regime': 'mid-IR-rovibrational'},
        'midir_window': {'label': 'OH 1-0', 'range_cm-1': [2600, 3600], 'count': 0},   # empty!
    }}}
    monkeypatch.setattr(ml, 'load_manifest', lambda: fake)
    monkeypatch.setattr(ml, 'VENDORED_DIR', tmp_path)
    monkeypatch.setattr(ml, 'ISPEC_MOLECULES_DIR', tmp_path / 'no_ispec')
    viol = sc.check_molecular_lists()
    assert any(v.invariant == 'molecular' and not v.tracked
               and 'mid-IR window' in v.quantity for v in viol)
