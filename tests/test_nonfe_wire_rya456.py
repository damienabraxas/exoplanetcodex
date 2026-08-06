"""
tests/test_nonfe_wire_rya456.py
===============================
RYA-456 — wire curate_nonfe_pools (RYA-395/398) into the default solar run.

These pin the WIRING CONTRACT, not the science numbers:
  * the curated kept pool flows into production A(X) by REUSING curate_nonfe_pools
    (the cull is not reimplemented; the abundance-blindness firewall stays live),
  * the run's non-Fe line set + A(X) reproduce the standalone curation
    (line-set identity + A(X) ≤0.05 dex vs the committed graded diagnostics),
  * the Cr canary lands at the RYA-398 floor (~+0.40), NOT at PASS / Asplund,
  * the Phase C classifier MAPS the curation's blind verdict (VALIDATED→PASS,
    RESIDUAL/LOW_CONFIDENCE→CURATION-OWED) — it never re-derives a threshold.

Heavy iSpec/MOOG and run-output tests skip cleanly when the engine or the
(gitignored) solar_abundances.csv is absent, matching the suite convention.
"""
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

import phase_c_verdict_rya371 as P  # noqa: E402
from pipeline import data_namespace as ns  # noqa: E402

# RYA-469: the solar baseline is now the FROZEN gold reference (committed), not the
# gitignored working file — so these assertions run in CI instead of skipping.
# Pinned to gold v1 ON PURPOSE (RYA-313). These three tests assert that the RYA-456
# wiring reproduced the RYA-398 GRADED CULL DIAGNOSTICS
# (curation_diagnostics_graded_rya398.csv, below) — a historical reproduction claim, so
# it must name the reference it reproduces. Following `current_version()` silently
# repointed them at gold v2 when RYA-522 re-froze, where the comparison is meaningless:
# v2 drops the `curation_verdict` column and deliberately HOLDS Ca/Ti/Ni/Na/Al blank at
# the ratified `owed` tier (RYA-522/596), so a v1-era value check reads nan-vs-6.324.
# Tests that should track the CURRENT pointer live in test_data_namespacing_rya469.py.
AB = ns.reference_path('v1')
GRADED_DIAG = ROOT / 'data' / 'curation' / 'nonfe_pools' / 'curation_diagnostics_graded_rya398.csv'
NONFE_DIR = ROOT / 'data' / 'curation' / 'nonfe_pools'


def _moog_available():
    try:
        import pipeline.abundances_derive as ad  # noqa: F401
        sys.path.insert(0, str(ad.ISPEC_DIR))
        import ispec
        return bool(ispec.is_moog_support_enabled())
    except Exception:
        return False


# ── Wiring locus / reuse (no engine needed) ──────────────────────────────────

def test_nonfe_wire_excludes_fe_cno_li_p():
    import pipeline.abundances_derive as ad
    assert {'Fe', 'C', 'N', 'O', 'Li', 'P'} <= ad._CURATED_NONFE_EXCLUDE
    for el in ('Cr', 'Si', 'Ti', 'Ni', 'Ca', 'Na', 'Al', 'Mn', 'Sr', 'Ba', 'Eu'):
        assert el not in ad._CURATED_NONFE_EXCLUDE, el


def test_nonfe_wire_reuses_curation_module_not_reimplemented():
    """The wiring helper must call into curate_nonfe_pools (cull stays there)."""
    import inspect
    import pipeline.abundances_derive as ad
    src = inspect.getsource(ad._curate_nonfe_pool) + inspect.getsource(ad._curated_nonfe_rows)
    assert 'curate_nonfe_pools' in src
    assert 'apply_cull' in src and 'grade_restrict=True' in src   # RYA-398 firewall
    assert 'compute_lte_abundances' in src                        # curation's own A(X)
    # no parallel cull thresholds reinvented here
    assert 'REW_LINEAR_CEILING' not in src and 'BLEND_FRAC_MAX' not in src


def test_nonfe_wire_keeps_blindness_firewall_live():
    """The abundance-blindness firewall must still raise on the wired path."""
    import pipeline.curate_nonfe_pools as cur
    with pytest.raises(ValueError):
        cur.assert_no_astrophysical_gf(pd.DataFrame({'log_gf_astro': [0.0]}), 'wired')


# ── Phase C verdict mapping (pure classification) ────────────────────────────

def _classify(cverdict, el='Cr', a_meas=6.02, asp=5.62, delta=0.40, n=7):
    return P._classify(el, asp, a_meas, delta, 0.3, n, None, False, '',
                        True, True, None, cverdict)[0]


def test_phase_c_maps_curation_verdict():
    assert _classify('VALIDATED', el='Ca', a_meas=6.31, asp=6.30, delta=0.01, n=8) == 'PASS'
    assert _classify('RESIDUAL') == 'CURATION-OWED'
    assert _classify('LOW_CONFIDENCE', el='Ni', a_meas=6.95, asp=6.20, delta=0.75, n=2) == 'CURATION-OWED'


def test_phase_c_validate_dont_tune_no_pass_without_validated():
    """A produced metal at a large offset must NOT be PASS unless the blind cull
    said VALIDATED — the verdict never fits the value to the anchor."""
    assert _classify('RESIDUAL', a_meas=6.02, asp=5.62, delta=0.40) != 'PASS'
    assert _classify('LOW_CONFIDENCE', a_meas=4.96, asp=2.83, delta=2.13) != 'PASS'


# ── Faithfulness vs the committed standalone curation (needs the run output) ──

_CURATED_ION_I = ['Si', 'Ca', 'Ni', 'Ti', 'Cr']


@pytest.mark.skipif(not AB.exists(), reason="needs a solar run (solar_abundances.csv is gitignored)")
def test_curation_faithful_abundance_reproduces_diagnostics():
    ab = pd.read_csv(AB, comment='#')
    diag = pd.read_csv(GRADED_DIAG).set_index('element')
    for el in _CURATED_ION_I:
        row = ab[(ab['element'] == el) & (ab['ion'] == 'I')]
        assert len(row) == 1, el
        a_run = float(row.iloc[0]['A_X'])
        a_cur = float(diag.loc[el, 'A_lte_curated'])
        assert abs(a_run - a_cur) <= 0.05, (el, a_run, a_cur)        # A(X) consistency


@pytest.mark.skipif(not AB.exists(), reason="needs a solar run (gitignored output)")
def test_cr_canary_lands_at_floor_not_pass():
    ab = pd.read_csv(AB, comment='#')
    cr = ab[(ab['element'] == 'Cr') & (ab['ion'] == 'I')]
    assert len(cr) == 1
    cr = cr.iloc[0]
    assert str(cr['curation_verdict']) == 'RESIDUAL'                 # not VALIDATED
    a_nlte = float(cr['A_X_nlte'])
    assert abs(a_nlte - 6.02) < 0.05, a_nlte                          # ~+0.40 vs Asplund 5.62
    assert (a_nlte - 5.62) > 0.25, a_nlte                             # NOT snapped to the anchor


@pytest.mark.skipif(not AB.exists(), reason="needs a solar run (gitignored output)")
def test_fe_legs_present_and_validated():
    ab = pd.read_csv(AB, comment='#')
    fe1 = ab[(ab['element'] == 'Fe') & (ab['ion'] == 'I')].iloc[0]
    fe2 = ab[(ab['element'] == 'Fe') & (ab['ion'] == 'II')]
    assert int(fe1['n_lines']) >= 60                                 # deep Fe I pool untouched
    assert abs(float(fe1['A_X_nlte']) - 7.516) < 0.01                # validated leg (RYA-407/446)
    assert len(fe2) == 1
    # Fe carries no curation verdict (its own pool, RYA-279/347 — not RYA-456-wired)
    assert not isinstance(fe1.get('curation_verdict', None), str) or fe1.get('curation_verdict') != fe1.get('curation_verdict') or pd.isna(fe1.get('curation_verdict'))


# ── Line-set identity vs the committed graded cull (needs MOOG) ───────────────

@pytest.mark.skipif(not _moog_available(), reason="needs iSpec MOOG engine")
def test_curation_faithful_lineset_identity():
    import pipeline.abundances_derive as ad
    pool = ad._curate_nonfe_pool(_CURATED_ION_I)
    for el in _CURATED_ION_I:
        cull = pd.read_csv(NONFE_DIR / f'{el}_cull_graded_rya398.csv', comment='#')
        committed = sorted(round(float(w), 3) for w in
                           cull[(cull['kept']) & (cull['ion'] == 'I') & (cull['A_lte'].notna())]
                           ['wavelength_air_A'])
        mine = sorted(round(float(w), 3) for w in
                      pool[(pool['element'] == el) & (pool['kept']) & (pool['ion'] == 'I')
                           & (pool['A_lte'].notna())]['wavelength_air_A'])
        assert committed == mine, (el, committed, mine)             # SAME kept lines
