"""
tests/test_oi_6300_audit_rya367.py
==================================
RYA-367 — solar [O I] 6300 atomic-data audit (Step 0) + the Step-1 [O I] gf
provenance pin.

Verifies the two CRITICAL gate conditions hold on the live branch:
  1. Ni I 6300.34 resolves through gf_resolver to canonical -2.11 and the
     load-bearing synth carries that value (not the store-#2 VALD3 pair).
  2. [O I] 6300.304 gf is a cited primary (NIST grade A / Storey & Zeippen 2000),
     value -9.717, NOT generic VALD3.
"""
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pipeline.audit import oi_6300_region as audit
from pipeline.gf_resolver import resolve

REPO = Path(__file__).resolve().parents[1]
CANON = REPO / 'data' / 'linelists' / 'canonical_gf.csv'


# ── Step 0 gate ──────────────────────────────────────────────────────────────

class TestStep0Gate:
    @pytest.fixture(scope='class')
    def findings(self):
        return audit.run(verbose=False)

    def test_step0_passes(self, findings):
        assert findings['step0_pass'] is True

    def test_ni_resolves_to_canonical(self, findings):
        assert abs(findings['ni_resolver_gf'] - (-2.11)) < 0.01
        assert abs(findings['ni_synth_gf'] - (-2.11)) < 0.01

    def test_oi_is_cited_primary_not_generic(self, findings):
        assert abs(findings['oi_gf'] - (-9.717)) < 0.005
        assert findings['oi_matches_sz2000']
        assert findings['oi_ref'].upper() not in ('VALD3', 'VALD', '', 'NAN')

    def test_store2_divergence_detected_but_not_fatal(self, findings):
        # store-#2 still carries the raw VALD3 pair (2 comps, ~-2.70) — reported,
        # not a STOP (it is not the synthesis path).
        assert findings['ni_store2_ncomp'] == 2
        assert findings['store2_diverges'] is True


def test_resolver_direct():
    assert abs(resolve((28, 1), 6300.337, 4.266) - (-2.11)) < 0.01
    assert abs(resolve((8, 1), 6300.304, 0.0) - (-9.717)) < 0.005


# ── Step 1 provenance pin ────────────────────────────────────────────────────

class TestOiProvenancePin:
    def test_canonical_cites_storey_zeippen(self):
        df = pd.read_csv(CANON, low_memory=False)
        row = df[(df['key_z'].astype(str).str.replace(r'\.0$', '', regex=True) == '8')
                 & (np.abs(df['wavelength_air_A'] - 6300.304) <= 0.01)]
        assert len(row) == 1
        r = row.iloc[0]
        assert abs(float(r['log_gf']) - (-9.717)) < 1e-6      # value unchanged
        assert 'Storey & Zeippen 2000' in str(r['loggf_reference'])
        assert str(r['adjudication_status']) == 'adjudicated_rya367'

    def test_provenance_check_script_passes(self):
        r = subprocess.run(
            ['python3', 'scripts/adjudicate_oi6300_provenance_rya367.py', '--check'],
            cwd=str(REPO), capture_output=True, text=True)
        assert r.returncode == 0, r.stdout + r.stderr


# ── Ni structural consistency (Step 2) ───────────────────────────────────────

def test_ni_single_component_canonical():
    df = pd.read_csv(CANON, low_memory=False)
    row = df[(df['key_z'].astype(str).str.replace(r'\.0$', '', regex=True) == '28')
             & (df['ion'].astype(float) == 1.0)
             & (np.abs(df['wavelength_air_A'] - 6300.342) <= 0.02)]
    assert len(row) == 1                                   # single combined line
    assert int(row.iloc[0]['hfs_n_components']) == 1
    assert abs(float(row.iloc[0]['log_gf']) - (-2.11)) < 1e-6
