"""
tests/test_ir_triage_rya376.py
==============================
RYA-376 (re-baselined 2026-06-30) — guard the 27-element IR triage + the Group-A
depth-confirmation artifacts produced by scripts/audit_vald_inventory_rya376.py.

Fast checks only: imports the static triage table and reads the small generated CSVs
(ir_triage.csv, groupA_depth_confirm.csv). Does NOT load the large LFS line lists.
"""
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'data' / 'audit' / 'vald_inventory'
sys.path.insert(0, str(ROOT / 'scripts'))

import audit_vald_inventory_rya376 as audit  # noqa: E402

# The 27 target elements the audit must cover (matches elements_master.json).
TARGETS = {'C', 'N', 'O', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'K', 'Ca', 'Sc', 'Ti',
           'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Sr', 'Y', 'Zr', 'Ba',
           'Eu', 'Li'}


def test_triage_covers_all_27_elements_with_valid_tags():
    assert set(audit.IR_TRIAGE) == TARGETS, "every target element must carry an IR tag"
    valid = {'IR-PRIMARY', 'IR-CROSSCHECK', 'IR-NONE'}
    for el, (group, tag, _why) in audit.IR_TRIAGE.items():
        assert group in {'A', 'B', 'C'}, el
        assert tag in valid, el


def test_group_membership_matches_ticket_triage():
    by_group = {g: {el for el, (gg, _t, _w) in audit.IR_TRIAGE.items() if gg == g}
                for g in ('A', 'B', 'C')}
    # Group A = flagship + IR-best channels (extraction priority)
    assert by_group['A'] == {'C', 'N', 'O', 'S', 'P', 'K', 'Na', 'Mg', 'Al'}
    # Group C = neutron-capture heavies + Li → UV/near-blue leg
    assert by_group['C'] == {'Ba', 'Y', 'Zr', 'Sr', 'Eu', 'Li'}
    # Fe is the deliberate Group-B exception (IR sigma ruler); Zn is the unlisted 27th
    assert audit.IR_TRIAGE['Fe'][:2] == ('B', 'IR-CROSSCHECK')
    assert audit.IR_TRIAGE['Zn'][:2] == ('B', 'IR-CROSSCHECK')


@pytest.mark.skipif(not (OUT / 'ir_triage.csv').exists(), reason="artifact not generated")
def test_ir_triage_artifact_confirms_groupAB_have_nir():
    df = pd.read_csv(OUT / 'ir_triage.csv')
    assert len(df) == 27
    ab = df[df['group'].isin(['A', 'B'])]
    # Every IR-PRIMARY / IR-CROSSCHECK element must actually hold NIR data.
    assert ab['nir_held'].all(), "Group A/B element tagged IR but no NIR held → defect"


@pytest.mark.skipif(not (OUT / 'groupA_depth_confirm.csv').exists(), reason="artifact not generated")
def test_groupA_named_diagnostics_all_covered_with_depth():
    df = pd.read_csv(OUT / 'groupA_depth_confirm.csv')
    assert df['covered'].all(), "a named Group-A IR diagnostic is missing from the held list"
    # The S I 1.045 um channel must be materially deeper than the weak optical S I 6757.
    s_triplet = df[df['diagnostic'].str.startswith('S I 1.045')]
    assert s_triplet['max_central_depth'].max() > 0.10
