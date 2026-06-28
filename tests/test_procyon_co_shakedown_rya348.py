"""
tests/test_procyon_co_shakedown_rya348.py
=========================================
RYA-348 — guard the Procyon C/O shakedown harness's reporting logic (NOT the synthesis
run, which is a heavy Turbospectrum job). Pins the things that make the deliverable
honest: the differential denominator is OUR measured Sun, the gap table records the
multi-instrument wiring reality (UVES O I 777 primary-O NOT wired for Procyon), and the
clean-fit threshold the synth-quality read uses.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import importlib.util  # noqa: E402

spec = importlib.util.spec_from_file_location(
    'procyon_co_shakedown', ROOT / 'scripts' / 'procyon_co_shakedown.py')
sh = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sh)


def test_differential_denominator_is_our_sun():
    # RYA-348 §4: differential anchor is OUR measured Sun (Phase C), NOT Asplund.
    assert sh.SOLAR_OURS['C'] == 8.491
    assert sh.SOLAR_OURS['O'] == 8.735
    assert 'Phase C' in sh.SOLAR_PROV


def test_clean_fit_threshold():
    assert sh.CHI2_CLEAN == 10.0      # brief Step 3: χ²ᵣ < 10 = fits cleanly


def test_gap_table_records_unwired_primary_O():
    w = sh.gap_table()
    # The UVES O I 777 primary-O arm exists as a config but is NOT wired for Procyon.
    oi777 = w[w['diagnostic'].str.contains('7771')].iloc[0]
    assert 'NOT WIRED' in oi777['wiring']
    assert oi777['role'] == 'O PRIMARY'
    # HARPS arms are the ones actually run here.
    harps = w[w['instrument'] == 'HARPS VIS']
    assert (harps['wiring'] == 'WIRED (run here)').all()
    assert len(harps) >= 4
    # UV and IR arms are recon-confirmed not-wired (build / telluric-gated).
    assert (w['wiring'].str.contains('NOT WIRED')).sum() >= 4


def test_registry_preflight_runs():
    # The RYA-463 pre-flight must produce a heads-up string for Procyon's params.
    txt = sh.preflight_registry(6554, 0.01)
    assert isinstance(txt, str) and 'Procyon' in txt
