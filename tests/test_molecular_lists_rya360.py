"""
tests/test_molecular_lists_rya360.py
====================================
RYA-360 — the C/N/O molecular line lists are vendored + guarded. Verifies the
vendored secure record, the provenance manifest (incl. the RYA-499 wavelength_coverage
regime field), the reproduced RYA-236 HARPS-window baselines, the measured mid-IR
absence (OH/NH/CH) vs CO's mid-IR reach, the [molecular] stewardship invariant passing
today, and — the point — that a missing/emptied list fails the guard loudly.
"""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pipeline import molecular_lists as ml           # noqa: E402
import scripts.check_stewardship as sc               # noqa: E402

_REQUIRED = {'CH', 'CN', 'C2', 'OH', 'NH', 'CO'}
_HEADLINE = {'CH': 583, 'CN': 3534, 'C2': 1019}       # RYA-236 baselines to reproduce


@pytest.fixture(scope='module')
def manifest():
    return ml.load_manifest()


# ── vendored secure record present + complete ─────────────────────────────────
def _held(manifest):
    # the RYA-360 held iSpec-bundle lists (RYA-503 acquisitions have their own tests)
    return {m: e for m, e in manifest['molecules'].items()
            if not str(e.get('origin', '')).startswith('acquired')}


def test_all_required_molecules_vendored(manifest):
    assert _REQUIRED <= set(_held(manifest))
    for mol, e in manifest['molecules'].items():        # every entry (held + acquired) present
        sub = ml.VENDORED_DIR / e['vendored_subdir']
        assert e['files'], mol
        for f in e['files']:
            p = sub / f
            assert p.exists() and p.stat().st_size > 0, f"{mol}/{f} missing/empty"


def test_manifest_provenance_and_coverage_complete(manifest):
    for mol, e in _held(manifest).items():
        assert e['source'] and e['distribution'], mol
        wc = e['wavelength_coverage']
        assert wc['min_A'] and wc['max_A'] and wc['regime'], mol
        # electronic vs mid-IR regime is machine-recorded (RYA-499)
        if mol == 'CO':
            assert 'mid-IR' in wc['regime']
        else:
            assert wc['regime'] == 'electronic-optical'


def test_harps_headline_counts_reproduce_rya236(manifest):
    for mol, expected in _HEADLINE.items():
        assert manifest['molecules'][mol]['harps_window']['count'] == expected, mol


def test_line_counts_match_live_files(manifest):
    # the recorded baseline equals a fresh count of the vendored bytes (guard basis)
    for mol, e in manifest['molecules'].items():
        sub = ml.VENDORED_DIR / e['vendored_subdir']
        live = sum(ml.count_bsyn_lines(sub / f) for f in e['files'])
        assert live == e['line_count'], f"{mol}: {live} != {e['line_count']}"


# ── the RYA-499 mid-IR verdict (measured, not inferred) ───────────────────────
def test_ohnhch_midir_absent_by_measurement(manifest):
    for mol in ('OH', 'NH', 'CH'):
        w = manifest['molecules'][mol]['midir_window']
        assert w['count'] == 0, f"{mol} unexpectedly has mid-IR rows"
        assert 'ABSENT' in w['verdict']                      # RYA-503 unblocked


def test_co_reaches_midir(manifest):
    reach = manifest['molecules']['CO']['midir_reach']
    assert reach['CO 1-0 fundamental ~4.6 µm (46000 Å)'] is True
    assert reach['CO 2-0 overtone ~2.3 µm (23000 Å)'] is True


# ── the stewardship [molecular] invariant ─────────────────────────────────────
def test_molecular_invariant_passes_today():
    viol = sc.check_molecular_lists()
    assert viol == [], f"expected no molecular violations, got {[v.locus for v in viol]}"


def test_negative_missing_list_fails_loud(monkeypatch, tmp_path):
    """Remove/empty a required list → UNTRACKED violation (exit-1 class). Simulated by
    pointing the guard at a manifest whose vendored file is absent."""
    fake = {'molecules': {'CH': {
        'vendored_subdir': 'CH', 'files': ['12CH_400-450.bsyn'], 'line_count': 3238,
        'source': 'src', 'distribution': 'dist', 'harps_range_count': 100,
        'harps_window': {'name': 'CH A-X G-band', 'range_A': [4290, 4315], 'count': 583},
        'wavelength_coverage': {'min_A': 4200.0, 'max_A': 9200.0, 'regime': 'electronic-optical'},
    }}}
    monkeypatch.setattr(ml, 'load_manifest', lambda: fake)
    monkeypatch.setattr(ml, 'VENDORED_DIR', tmp_path)                 # CH/12CH... absent
    monkeypatch.setattr(ml, 'ISPEC_MOLECULES_DIR', tmp_path / 'no_ispec')
    viol = sc.check_molecular_lists()
    untracked = [v for v in viol if v.invariant == 'molecular' and not v.tracked]
    assert untracked, "a missing molecular list must fail loudly (untracked)"
    assert any('present' in v.quantity for v in untracked)


def test_negative_emptied_list_fails_loud(monkeypatch, tmp_path):
    """A present-but-truncated list (count < baseline) is also a loud failure."""
    (tmp_path / 'CH').mkdir()
    (tmp_path / 'CH' / '12CH_400-450.bsyn').write_text(
        "'         0106.000012'    1      3238\n'CH PGopher'\n")  # header only, 0 data rows
    fake = {'molecules': {'CH': {
        'vendored_subdir': 'CH', 'files': ['12CH_400-450.bsyn'], 'line_count': 3238,
        'source': 'src', 'distribution': 'dist', 'harps_range_count': 100,
        'wavelength_coverage': {'min_A': 4200.0, 'max_A': 9200.0, 'regime': 'electronic-optical'},
    }}}
    monkeypatch.setattr(ml, 'load_manifest', lambda: fake)
    monkeypatch.setattr(ml, 'VENDORED_DIR', tmp_path)
    monkeypatch.setattr(ml, 'ISPEC_MOLECULES_DIR', tmp_path / 'no_ispec')
    viol = sc.check_molecular_lists()
    assert any(v.invariant == 'molecular' and not v.tracked
               and 'non-empty' in v.quantity for v in viol)


# ── Finding B: cno_synthesis.py comment corrected ─────────────────────────────
def test_cno_synthesis_comment_is_regime_qualified():
    src = (ROOT / 'pipeline' / 'cno_synthesis.py').read_text()
    # the old unqualified "held" list is gone
    assert 'CH Masseron, CN Brooke+Sneden,\nC2, OH, NH, CO' not in src
    assert 'ELECTRONIC bands only' in src
    assert 'RYA-503' in src and 'CO_IR_Li2015' in src
