"""
tests/test_ace_co_3d_rya442.py
==============================
RYA-442 — guard the 1D->3D determination. Fast, deterministic checks of the
literature determination, the verdict logic, and the model-acquisition finding.
The concrete <3D> synthesis probe (~45 s, needs iSpec + the <3D> model) is the CLI
smoke test `python scripts/ace_co_3d_probe.py --validate`, not CI.
"""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))
import ace_co_3d_probe as p3  # noqa: E402


def test_full3d_literature_all_within_band():
    s1 = p3.step1_literature()
    # Every cited FULL-3D molecular-C abundance lands within 0.10 dex of 8.46.
    assert s1['all_full3d_within_band'] is True
    lo, hi = s1['full_3d_min_max']
    assert 8.35 <= lo <= 8.55 and 8.45 <= hi <= 8.55


def test_implied_correction_is_negative_and_closes_gap():
    s1 = p3.step1_literature()
    lo, hi = s1['implied_1D_to_3D_correction_band']
    assert hi < 0                        # correction is NEGATIVE (toward reference)
    # applying the band to the 1D value lands within the reopen band
    for corr in (lo, hi):
        assert abs((p3.A_C_1D_MU1_441 + corr) - p3.ace.A_C_REF) <= p3.NEAR_REF + 1e-9


def test_verdict_is_reopen():
    rep = p3.run(do_probe=False)         # literature-only, fast
    assert rep['reopen'] is True
    assert rep['verdict'].startswith('REOPEN')
    assert "wrong sign" in rep['441_record_correction']


def test_constants_single_source():
    assert p3.ace.A_C_REF == 8.46
    assert p3.A_C_1D_MU1_441 == 8.646    # RYA-441 result
    assert p3.NEAR_REF == 0.10


def test_3d_model_acquired_with_provenance():
    assert p3.MOD3D.exists(), "the <3D> STAGGER solar model should be on disk"
    prov = p3.MOD3D.parent / 'PROVENANCE.md'
    assert prov.exists()
    txt = prov.read_text()
    assert 'STAGGER' in txt and 'Magic' in txt and 'INGESTION CAVEAT' in txt
    head = p3.MOD3D.read_text().splitlines()
    assert head[0].strip().startswith('5777')        # Teff
    assert any('TAU5000 SCALE' in ln for ln in head[:12])
