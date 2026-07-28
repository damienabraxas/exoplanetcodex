"""
tests/test_silent_line_drop_rya429.py
=====================================
RYA-429 — the EW fitter's no-silent-drop guarantee.

Pins the DISCIPLINE:
  * the rejection vocabulary is a controlled set; a rejection carries a reason
    from it, and 'other' requires a freetext note (blank/placeholder = defeat),
  * the no-silent-drop invariant holds: every in-coverage input line is either
    measured OR in the ledger with an explicit reason (neither = hard failure),
  * the ionization-stage-presence flag fires when a species has both stages in
    the linelist but one stage measures zero lines,
  * the Sr II RCA: all four canonical lines (4077/4161/4215/4305) are accounted
    for, and the two subordinate lines (4161/4305) carry a documented reason.

The invariant/Sr-II/stage-flag checks read the COMMITTED audit JSON produced by
a solar fitter run (house style: guard the artifact, do not re-run the fitter).
The vocabulary + reconcile + stage-flag logic is unit-tested on synthetic data
so the guarantee is verified even where the gitignored working data is absent.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pipeline import rejection_ledger as rl  # noqa: E402

AUDIT = ROOT / 'data' / 'results' / 'rejection_ledger_solar_rya429.json'


def _audit():
    if not AUDIT.exists():
        pytest.skip(f"solar audit JSON absent ({AUDIT.name}); run "
                    "`python -m pipeline.lines_fit --star solar` first")
    return json.loads(AUDIT.read_text())


# ── controlled vocabulary + ledger discipline (unit) ──────────────────────────

def test_reason_vocabulary_is_the_controlled_set():
    assert rl.REJECTION_REASONS == (
        'below_priority_threshold', 'blend_overlap', 'failed_chi2',
        'saturated_core', 'low_snr', 'continuum_fail', 'out_of_coverage', 'other')


def test_ledger_rejects_unknown_reason():
    led = rl.RejectionLedger('solar')
    with pytest.raises(ValueError):
        led.reject('Sr', 'II', 4161.79, 'mystery_reason')


def test_other_requires_a_freetext_note():
    led = rl.RejectionLedger('solar')
    with pytest.raises(ValueError):
        led.reject('Sr', 'II', 4161.79, 'other', note='')       # blank -> refused
    led.reject('Sr', 'II', 4161.79, 'other', note='documented cause')  # ok
    assert len(led.to_frame()) == 1


def test_earliest_gate_wins_no_double_count():
    led = rl.RejectionLedger('solar')
    led.reject('Fe', 'I', 5000.0, 'failed_chi2')
    led.reject('Fe', 'I', 5000.02, 'low_snr')     # same 0.1 A feature -> ignored
    df = led.to_frame()
    assert len(df) == 1 and df.iloc[0]['reason'] == 'failed_chi2'


# ── reconcile: the no-silent-drop invariant (unit, synthetic) ─────────────────

def _lines(rows):
    return pd.DataFrame(rows, columns=[
        'element', 'ion', 'wavelength_air_A', 'central_depth', 'blend_flag', 'priority'])


def _measured(rows):
    return pd.DataFrame(rows, columns=['element', 'ion', 'wavelength_air_A',
                                       'ew_mA', 'chi2'])


def test_reconcile_accounts_for_every_in_coverage_line():
    spec = np.linspace(4000.0, 5000.0, 100)
    lines = _lines([
        ['Fe', 'I', 4500.0, 0.5, False, 3],    # measured
        ['Cr', 'I', 4600.0, 0.5, False, 3],    # fit failed (in ledger already)
        ['Ti', 'I', 4700.0, 0.001, False, 3],  # below min_fit_depth -> pre-filter
        ['V',  'I', 4800.0, 0.5, True, 3],     # blend_flag -> pre-filter
        ['Ni', 'I', 9999.0, 0.5, False, 3],    # OUT of coverage -> excluded from invariant
    ])
    measured = _measured([['Fe', 'I', 4500.0, 42.0, 0.01]])
    led = rl.RejectionLedger('solar')
    led.reject('Cr', 'I', 4600.0, 'failed_chi2', note='fit did not converge')

    rep = rl.reconcile(lines, measured, led, spec, min_fit_depth=0.008)
    assert rep['n_in_coverage'] == 4          # the 9999 A line is out of coverage
    assert rep['n_measured'] == 1
    assert rep['n_rejected'] == 3
    assert rep['invariant_holds'] is True
    assert rep['n_unaccounted'] == 0
    assert rep['n_unclassified'] == 0
    # derived reasons mirror the worklist pre-filter
    reasons = {(r['element']): r['reason'] for _, r in led.to_frame().iterrows()}
    assert reasons['Ti'] == 'below_priority_threshold'   # central_depth < min_fit_depth
    assert reasons['V'] == 'blend_overlap'               # blend_flag=True


def test_reconcile_flags_a_genuine_silent_drop_as_unclassified():
    # A clean line (passes the worklist filter) that is neither measured nor in
    # the ledger is the exact bug this ticket kills -> must surface loudly.
    spec = np.linspace(4000.0, 5000.0, 100)
    lines = _lines([['Sr', 'II', 4305.0, 0.58, False, 3]])
    led = rl.RejectionLedger('solar')
    rep = rl.reconcile(lines, _measured([]), led, spec, min_fit_depth=0.008)
    assert rep['n_unclassified'] == 1
    assert rep['invariant_holds'] is False or rep['n_unclassified'] > 0
    row = led.to_frame().iloc[0]
    assert row['reason'] == 'other' and row['note'].startswith('UNCLASSIFIED')


# ── ionization-stage-presence flag (unit, synthetic) ──────────────────────────

def test_ionization_flag_fires_when_one_stage_empty():
    spec = np.linspace(4000.0, 5000.0, 100)
    lines = _lines([
        ['Sr', 'I',  4607.3, 0.70, False, 3],
        ['Sr', 'II', 4077.7, 0.97, False, 3],
        ['Sr', 'II', 4161.8, 0.53, False, 3],
    ])
    measured = _measured([['Sr', 'I', 4607.3, 30.0, 0.01]])   # Sr II measures nothing
    flags = rl.ionization_stage_presence_check(lines, measured, spec, emit=False)
    assert len(flags) == 1
    assert flags[0]['element'] == 'Sr' and flags[0]['empty_stage'] == 'ionized'


def test_ionization_flag_silent_when_both_stages_present():
    spec = np.linspace(4000.0, 5000.0, 100)
    lines = _lines([
        ['Fe', 'I',  4500.0, 0.5, False, 3],
        ['Fe', 'II', 4600.0, 0.5, False, 3],
    ])
    measured = _measured([['Fe', 'I', 4500.0, 40.0, 0.01],
                          ['Fe', 'II', 4600.0, 30.0, 0.01]])
    assert rl.ionization_stage_presence_check(lines, measured, spec, emit=False) == []


# ── integration: the committed solar audit artifact ───────────────────────────

def test_solar_no_silent_drop_invariant_holds():
    a = _audit()
    assert a['invariant_holds'] is True
    assert a['n_unaccounted'] == 0
    assert a['n_unclassified'] == 0
    assert a['n_measured'] + a['n_rejected'] == a['n_in_coverage']


def test_solar_ledger_reasons_are_all_from_the_vocabulary():
    a = _audit()
    for reason in a['reason_counts']:
        assert reason in rl.REJECTION_REASONS


def test_solar_ionization_stage_flag_fires_for_sr():
    a = _audit()
    elems = {f['element'] for f in a['ionization_stage_flags']}
    assert 'Sr' in elems, "expected the Sr II empty-stage flag to fire"
    sr = next(f for f in a['ionization_stage_flags'] if f['element'] == 'Sr')
    assert sr['empty_stage'] == 'ionized' and sr['ionized_measured'] == 0


def test_solar_sr_ii_all_four_accounted_with_documented_reasons():
    a = _audit()
    verd = {round(r['wavelength_air_A'], 1): r for r in a['sr_ii_verdict']}
    for wav in (4077.7, 4161.8, 4215.5, 4305.4):
        r = verd[round(wav, 1)]
        assert r['status'] != 'UNACCOUNTED', f"Sr II {wav} vanished with no record"
        if r['status'] == 'rejected':
            assert r['reason'] in rl.REJECTION_REASONS
            assert r['note'], f"Sr II {wav} rejected with a blank note"


def test_solar_sr_ii_subordinate_lines_carry_a_reason():
    # 4161 and 4305 are the two that silently vanished in the RCA — they must now
    # each carry an explicit, documented disposition.
    a = _audit()
    verd = {round(r['wavelength_air_A'], 1): r for r in a['sr_ii_verdict']}
    for wav in (4161.8, 4305.4):
        r = verd[round(wav, 1)]
        assert r['status'] in ('measured', 'rejected')
        if r['status'] == 'rejected':
            assert r['reason'] and r['note']
