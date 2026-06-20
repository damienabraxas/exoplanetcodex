"""
tests/test_gf_store_consistency_rya368.py
=========================================
RYA-368 — all-stores gf consistency: store-#2 (linelist_solar.csv) is reconciled
to the single canonical source at the consumer, and the CI invariant now covers
every gf store (orphan = hard fail; raw divergence = tracked-visible; the RYA-367
trigger lines must resolve to canonical).
"""
import numpy as np
import pandas as pd
import pytest

from pipeline import gf_resolver as gr
from pipeline.audit import gf_store_consistency as gs


# ── resolve_df_gf (the safe per-line consumer route) ─────────────────────────

class TestResolveDfGf:
    def test_resolves_trigger_lines_to_canonical(self):
        df = pd.DataFrame([
            {'element': 'O', 'ion': 'I', 'wavelength_air_A': 6300.304,
             'excitation_potential_eV': 0.0, 'log_gf': -9.776},   # raw VALD3
            {'element': 'Ni', 'ion': 'I', 'wavelength_air_A': 6300.337,
             'excitation_potential_eV': 4.266, 'log_gf': -2.841},  # raw VALD3 comp
        ])
        out, stats = gr.resolve_df_gf(df, keep_unresolved=True)
        assert stats['n_resolved'] == 2 and stats['n_changed'] == 2
        oi = out[out['element'] == 'O'].iloc[0]
        ni = out[out['element'] == 'Ni'].iloc[0]
        assert abs(oi['log_gf'] - (-9.717)) < 0.005      # → canonical (S&Z 2000)
        assert abs(ni['log_gf'] - (-2.11)) < 0.01        # → canonical (Johansson 2003)
        assert bool(oi['gf_canonical']) and bool(ni['gf_canonical'])

    def test_keep_unresolved_keeps_raw(self):
        df = pd.DataFrame([{'element': 'Xx', 'ion': 'I',
                            'wavelength_air_A': 1234.567,
                            'excitation_potential_eV': 0.0, 'log_gf': -1.0}])
        out, stats = gr.resolve_df_gf(df, keep_unresolved=True)
        assert stats['n_unresolved'] == 1
        assert out.iloc[0]['log_gf'] == -1.0            # kept raw
        assert not bool(out.iloc[0]['gf_canonical'])

    def test_keep_unresolved_false_raises(self):
        df = pd.DataFrame([{'element': 'Xx', 'ion': 'I',
                            'wavelength_air_A': 1234.567,
                            'excitation_potential_eV': 0.0, 'log_gf': -1.0}])
        with pytest.raises((gr.GfResolutionError, ValueError)):
            gr.resolve_df_gf(df, keep_unresolved=False)


# ── store registry + per-store report ────────────────────────────────────────

class TestStoreReports:
    @pytest.fixture(scope='class')
    def reports(self):
        return {s.label: gs.store_report(s) for s in gs.STORES}

    def test_all_three_stores_registered(self):
        labels = {s.label for s in gs.STORES}
        assert {'atomic_lines.tsv', 'iSpec line regions',
                'linelist_solar.csv'} <= labels

    def test_no_orphans_anywhere(self, reports):
        # every store line must have a canonical home (the hard guarantee)
        for label, rep in reports.items():
            assert rep['n_orphan'] == 0, f"{label} has {rep['n_orphan']} orphan(s)"

    def test_store2_has_raw_divergence_reported(self, reports):
        # the landmine is visible, not silent
        s2 = reports['linelist_solar.csv']
        assert s2['n_divergent'] > 0
        assert s2['n_overlap'] > 0

    def test_store2_not_load_bearing(self):
        s2 = next(s for s in gs.STORES if s.label == 'linelist_solar.csv')
        assert s2.is_load_bearing_gf is False


# ── target resolution + the audit run() gate ─────────────────────────────────

class TestTargetsAndRun:
    def test_targets_resolve_to_canonical(self):
        for t in gs.check_targets():
            assert t['ok'], f"{t['label']} resolved to {t['got']}, expected {t['expect']}"

    def test_run_passes(self):
        out = gs.run(verbose=False)
        assert out['pass'] is True


# ── CI invariant in check_stewardship ────────────────────────────────────────

class TestStewardshipInvariant:
    def test_all_stores_invariant_tracked_only(self):
        import scripts.check_stewardship as cs
        vs = cs.check_all_stores_resolve()
        untracked = [v for v in vs if not v.tracked]
        assert untracked == [], f"unexpected untracked store violations: {untracked}"
        # raw-divergence summaries are present and tracked
        assert any(v.quantity.startswith('log gf (raw') for v in vs)

    def test_target_regression_guard_fires_untracked(self, monkeypatch):
        # simulate the [O I]/Ni landmine returning → must produce an UNTRACKED fail
        import scripts.check_stewardship as cs
        from pipeline.audit import gf_store_consistency as _gs
        monkeypatch.setattr(_gs, 'check_targets',
                            lambda: [{'label': '[O I] 6300.304', 'expect': -9.717,
                                      'got': -9.776, 'ok': False}])
        vs = cs.check_all_stores_resolve()
        untracked = [v for v in vs if not v.tracked]
        assert any('6300' in v.locus for v in untracked)
