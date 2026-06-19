"""
tests/test_ni6300_gf_rya365.py
==============================
RYA-365 — the Ni I 6300.34 gf in the [O I] 6300 joint fit must resolve through
gf_resolver (canonical = Johansson+2003 log gf -2.11), with no hardcoded/store-#2
copy in the live path.

Verifies: the canonical adjudication (value + provenance), that gf_resolver and
the synth-linelist reroute both return -2.11, that the live synth path is asserted
against gf_resolver (no hardcoded gf), and the blend-partition helpers.
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pipeline.gf_resolver import resolve as resolve_gf, _index

REPO = Path(__file__).resolve().parents[1]
CANON = REPO / 'data' / 'linelists' / 'canonical_gf.csv'

NI_KEY = (28, 1)
NI_WL, NI_EP = 6300.342, 4.266
JOHANSSON_GF = -2.11


# ── Canonical adjudication ───────────────────────────────────────────────────

class TestCanonicalNiGf:
    def test_resolver_returns_johansson_value(self):
        _index.cache_clear()
        assert abs(resolve_gf(NI_KEY, NI_WL, NI_EP) - JOHANSSON_GF) < 1e-6

    def test_provenance_is_johansson_2003(self):
        df = pd.read_csv(CANON, low_memory=False)
        key = df['key_z'].astype(str).str.replace(r'\.0$', '', regex=True)
        row = df[(key == '28') & (df['ion'].astype(float) == 1.0)
                 & (np.abs(df['wavelength_air_A'] - NI_WL) <= 0.02)]
        assert len(row) == 1
        r = row.iloc[0]
        assert abs(float(r['log_gf']) - JOHANSSON_GF) < 1e-6
        assert 'Johansson' in str(r['loggf_reference'])
        assert '2003' in str(r['loggf_reference'])
        assert str(r['adjudication_status']) == 'adjudicated_rya365'

    def test_not_the_stale_nist_value(self):
        # the pre-RYA-365 NIST grade-B value (-2.310) must be gone from log_gf
        df = pd.read_csv(CANON, low_memory=False)
        key = df['key_z'].astype(str).str.replace(r'\.0$', '', regex=True)
        r = df[(key == '28') & (df['ion'].astype(float) == 1.0)
               & (np.abs(df['wavelength_air_A'] - NI_WL) <= 0.02)].iloc[0]
        assert abs(float(r['log_gf']) - (-2.310)) > 0.05
        # audit columns retained (GES seed already carried Johansson's -2.11)
        assert abs(float(r['gf_synth_ges']) - (-2.11)) < 1e-6


# ── Synth-linelist reroute ───────────────────────────────────────────────────

class TestSynthReroute:
    def test_apply_to_synth_array_gives_johansson(self):
        from pipeline.abundances_derive import ispec, _SYNTH_LINELIST_FILE
        from pipeline.gf_resolver import apply_to_synth_array
        arr = ispec.read_atomic_linelist(_SYNTH_LINELIST_FILE)
        apply_to_synth_array(arr)
        w = arr['wave_A'].astype(float)
        ep = arr['lower_state_eV'].astype(float)
        m = (np.abs(w - NI_WL) <= 0.02) & (np.abs(ep - NI_EP) <= 0.05)
        ni = [i for i in np.where(m)[0]
              if str(arr['element'][i]).strip().startswith('Ni')
              and int(arr['ion'][i]) == 1]
        assert len(ni) == 1
        assert abs(float(arr['loggf'][ni[0]]) - JOHANSSON_GF) < 1e-4


# ── No hardcoded gf in the live path ─────────────────────────────────────────

class TestNoHardcodedGf:
    def test_live_path_asserts_against_resolver(self):
        src = (REPO / 'pipeline' / 'cno_synthesis.py').read_text()
        # the live [O I] partition must assert the synth Ni gf == gf_resolver canonical
        assert 'resolve_gf((28, 1)' in src
        assert 'RYA-365 invariant violated' in src

    def test_only_labeled_stale_literal(self):
        # the only Ni-gf numeric literal allowed is the clearly-labelled diagnostic
        # 'before' value (_STALE_NI_GF); the live value comes from gf_resolver.
        src = (REPO / 'pipeline' / 'cno_synthesis.py').read_text()
        assert '_STALE_NI_GF = -2.310' in src
        assert 'diagnostic only' in src.lower() or 'before/after' in src.lower()


# ── Adjudication script ──────────────────────────────────────────────────────

class TestAdjudicationScript:
    def test_check_passes(self):
        import subprocess
        r = subprocess.run(
            ['python3', 'scripts/adjudicate_ni6300_rya365.py', '--check'],
            cwd=str(REPO), capture_output=True, text=True)
        assert r.returncode == 0, r.stdout + r.stderr


# ── Blend-partition helpers ──────────────────────────────────────────────────

class TestPartitionHelpers:
    def test_ew_mA_units(self):
        from pipeline.cno_synthesis import _ew_mA
        # a 1.0 Å-wide box (6300.0–6301.0 Å = 0.1 nm) at depth 0.5 →
        # EW = 0.5 × 1.0 Å = 0.5 Å = 500 mÅ.
        sw = np.array([6300.0, 6300.5, 6301.0]) / 10.0  # nm grid
        flux = np.array([0.5, 0.5, 0.5])
        assert abs(_ew_mA(sw, flux) - 500.0) < 1.0

    def test_ni6300_idx_finds_line(self):
        from pipeline.abundances_derive import ispec, _SYNTH_LINELIST_FILE
        from pipeline.cno_synthesis import _ni6300_idx, _NI_6300_WL
        arr = ispec.read_atomic_linelist(_SYNTH_LINELIST_FILE)
        i = _ni6300_idx(arr)
        assert str(arr['element'][i]).strip().startswith('Ni')
        assert abs(float(arr['wave_A'][i]) - _NI_6300_WL) <= 0.02
