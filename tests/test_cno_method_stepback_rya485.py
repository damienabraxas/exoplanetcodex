"""
tests/test_cno_method_stepback_rya485.py
========================================
RYA-485 — guard the three CNO-method corrections.
  Issue 2 (solar O modeling verification, the correctness-critical one): the live nlte_cno
    path gives the Sun 3D-NLTE on O I 777 (not silent-1D); [O I] 6300's 3D term is ~0 (so
    forbidden-LTE is correct); and the differential regime-mismatch (Procyon-1D vs Sun-3D)
    is the named finding. Asserts against the live grid, plus the committed report JSON.
  Issue 3 (3 independent solar references): the fetch manifests are well-formed + cited
    (Baker Zenodo upstream-md5; Reiners CDS + Wallace NSO computed-md5).
"""
import json
import sys
import warnings
from pathlib import Path

warnings.filterwarnings('ignore')
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pipeline import nlte_cno as N                                # noqa: E402
from config.constants import get_star_params                     # noqa: E402

REPORT = ROOT / 'data' / 'results' / 'solar_o_modeling_verification_rya485.json'
MANI = ROOT / 'data' / 'sirius_manifest'


def _p(star):
    r = get_star_params(star)
    return float(r['teff']), float(r['logg']), float(r['feh_ref']), float(r.get('xi', 1.0))


# ── Issue 2: solar O is 3D/synthesis, not silent-LTE ──────────────────────────

def test_sun_gets_3d_nlte_on_oi777_not_silent_1d():
    te, lg, fe, xi = _p('solar')
    assert N.select_leg(te) == '3D'                      # Sun is under the 6500 K 3D ceiling
    d3 = N.cno_nlte_delta('OI', '777nm', te, lg, fe, xi, 8.9, leg='3D')
    assert d3 == d3 and abs(d3) > 0.1                    # 3D leg fires (non-NaN, ~-0.17)


def test_procyon_is_1d_leg_3d_unavailable_above_ceiling():
    te, lg, fe, xi = _p('procyon')
    assert te > N.TEFF_3D_CEILING and N.select_leg(te) == '1D'
    import math
    assert math.isnan(N.cno_nlte_delta('OI', '777nm', te, lg, fe, xi, 9.36, leg='3D'))


def test_oi6300_forbidden_lte_is_correct_3d_term_negligible():
    te, lg, fe, xi = _p('solar')
    d = N.cno_nlte_delta('OI', '630.0nm', te, lg, fe, xi, 8.73, leg='3D')
    assert abs(d) < 0.01                                 # ~+0.001 -> forbidden-LTE is right


def test_report_records_no_silent_lte_and_the_regime_finding():
    d = json.loads(REPORT.read_text())
    assert d['no_silent_lte'] is True
    assert d['sun_gets_3d_on_777'] is True
    assert d['forbidden_lte_is_correct'] is True
    # the named finding: differential regime-mismatch (~0.02 dex), not silently absorbed
    rm = d['regime_mismatch']
    assert abs(rm['regime_term_3d_minus_1d']) > 0.01
    assert rm['procyon_oh_vs_sun3d'] != rm['procyon_oh_vs_sun1d']
    # solar O spread does NOT mirror Procyon's 0.18 (Procyon continuum is UVES-specific)
    assert d['per_arm_solar_O']['mirrors_procyon'] is False


# ── Issue 3: three independent solar 777 references, cited ────────────────────

def test_three_solar_reference_manifests_present_and_cited():
    expect = {
        'solar_atlas_iag_baker2020.json': ('zenodo', 'upstream'),
        'solar_atlas_iag_reiners2016.json': ('cds', 'computed'),
        'solar_atlas_wallace2011_kpno.json': ('nso', 'computed'),
    }
    for fn, (src, md5kind) in expect.items():
        m = json.loads((MANI / fn).read_text())
        assert m['ticket'].startswith('RYA-485')
        assert m['items'] and m['n_items'] == len(m['items'])
        assert m['source_cite'] and m['purpose']
        if md5kind == 'upstream':                        # Zenodo: every item carries an md5
            assert all(len(it['md5']) == 32 for it in m['items'])
        else:                                            # CDS/NSO: no upstream md5 -> computed at fetch
            assert all(it.get('md5') is None for it in m['items'])
            assert 'computed-at-fetch' in m['checksum_provenance']
