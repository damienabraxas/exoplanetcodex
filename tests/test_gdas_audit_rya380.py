"""
tests/test_gdas_audit_rya380.py
===============================
RYA-380 — the codex-data-audit GDAS coverage check. Wavelength gate (red-optical/IR
needs GDAS, blue/UV does not), and the flag logic: a gated dataset with a MISSING night
is flagged; a gated dataset with all nights available is OK; an un-enumerated gated
dataset is flagged for back-fill.
"""
import datetime as dt

import pytest

from pipeline.telluric import gdas_audit as ga
from pipeline.telluric import gdas_fetch as gf

_HAS_ESO = (gf._eso_gdas_dir() is not None
            and (gf._eso_gdas_dir() / "gdas_profiles_C-70.4-24.6.tar.gz").exists())


# ── wavelength gate ──────────────────────────────────────────────────────────
def test_gate_red_optical_and_ir_need_gdas():
    assert ga.needs_gdas({'max_wave_A': 7699.0})    # K I 7699 red-optical
    assert ga.needs_gdas({'max_wave_A': 23000.0})   # K-band IR
    assert ga.needs_gdas({'telluric_gated': True})  # explicit


def test_gate_blue_uv_not_gated():
    assert not ga.needs_gdas({'max_wave_A': 3800.0})
    assert not ga.needs_gdas({'max_wave_A': 5000.0})  # mid-optical, telluric-clean
    assert not ga.needs_gdas({})                       # no info → not gated


# ── flag logic ───────────────────────────────────────────────────────────────
def test_blue_dataset_not_flagged():
    rep = ga.audit_gdas_coverage([{'name': 'UV', 'site': 'paranal',
                                   'nights': [], 'max_wave_A': 3000.0}])
    assert rep['ok'] and rep['n_flagged'] == 0
    assert 'NOT-GATED' in rep['datasets'][0]['verdict']


def test_gated_without_nights_is_flagged():
    rep = ga.audit_gdas_coverage([{'name': 'IR-no-nights', 'site': 'paranal',
                                   'nights': [], 'max_wave_A': 23000.0}])
    assert not rep['ok'] and rep['n_flagged'] == 1
    assert 'FLAG' in rep['datasets'][0]['verdict']


def test_missing_gdas_night_is_flagged(tmp_path):
    # gated dataset, an unknown site → every night MISSING/ERROR → flagged
    rep = ga.audit_gdas_coverage(
        [{'name': 'IR-unknown-site', 'site': 'atlantis',
          'nights': [dt.date(2022, 11, 23)], 'max_wave_A': 23000.0}],
        cache_dir=tmp_path)
    assert not rep['ok'] and rep['n_flagged'] == 1
    assert rep['datasets'][0]['nights'][0]['status'] in ('MISSING', 'ERROR')


@pytest.mark.skipif(not _HAS_ESO, reason="ESO GDAS tarball not installed")
def test_paranal_gated_nights_ok(tmp_path):
    rep = ga.audit_gdas_coverage(
        [{'name': 'Vesta IR', 'site': 'paranal',
          'nights': [dt.date(2022, 11, d) for d in (21, 23, 25)], 'max_wave_A': 23000.0}],
        cache_dir=tmp_path)
    assert rep['ok'] and rep['n_flagged'] == 0
    assert all(n['status'] in ('cached', 'fetchable') for n in rep['datasets'][0]['nights'])
