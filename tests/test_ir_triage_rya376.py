"""
tests/test_ir_triage_rya376.py
==============================
RYA-376 (re-baselined 2026-06-30, corrected per the 16:26 review) — guard the
CANONICAL-27 IR triage + the Group-A depth-confirmation artifacts produced by
scripts/audit_vald_inventory_rya376.py.

The canonical 27 (RYA-109) is 26 distinct elements with **Fe I and Fe II counted as
SEPARATE species** = 27, and **no Zn**. These tests lock that in so the earlier silent
miscount (collapse Fe I+Fe II → one "Fe", backfill Zn to 27) cannot recur.

Fast checks only: imports the static triage table and reads the small generated CSVs.
Does NOT load the large LFS line lists.
"""
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'data' / 'audit' / 'vald_inventory'
sys.path.insert(0, str(ROOT / 'scripts'))

import audit_vald_inventory_rya376 as audit  # noqa: E402

# The canonical 27 SPECIES (RYA-109) — Fe split, no Zn.
CANONICAL_27 = {
    'C', 'N', 'O', 'S', 'P', 'K', 'Na', 'Mg', 'Al',                       # A (9)
    'Si', 'Ca', 'Ti', 'Cr', 'Mn', 'V', 'Co', 'Ni', 'Sc', 'Cu', 'Fe I', 'Fe II',  # B (12)
    'Ba', 'Y', 'Zr', 'Sr', 'Eu', 'Li',                                   # C (6)
}


def test_triage_is_the_canonical_27_species_fe_split_no_zn():
    assert set(audit.IR_TRIAGE) == CANONICAL_27
    assert len(audit.IR_TRIAGE) == 27
    # Fe I and Fe II are SEPARATE species; bare "Fe" is not a key.
    assert 'Fe I' in audit.IR_TRIAGE and 'Fe II' in audit.IR_TRIAGE
    assert 'Fe' not in audit.IR_TRIAGE
    # Zn is non-canonical and must NOT be in the triage.
    assert 'Zn' not in audit.IR_TRIAGE


def test_tags_valid_and_group_membership():
    valid = {'IR-PRIMARY', 'IR-CROSSCHECK', 'IR-NONE'}
    for sp, (group, tag, _why) in audit.IR_TRIAGE.items():
        assert group in {'A', 'B', 'C'} and tag in valid, sp
    by_group = {g: {sp for sp, (gg, _t, _w) in audit.IR_TRIAGE.items() if gg == g}
                for g in ('A', 'B', 'C')}
    assert by_group['A'] == {'C', 'N', 'O', 'S', 'P', 'K', 'Na', 'Mg', 'Al'}
    assert by_group['B'] == {'Si', 'Ca', 'Ti', 'Cr', 'Mn', 'V', 'Co', 'Ni', 'Sc', 'Cu',
                             'Fe I', 'Fe II'}
    assert by_group['C'] == {'Ba', 'Y', 'Zr', 'Sr', 'Eu', 'Li'}
    # Fe I = sigma ruler; Fe II = sparse-in-IR cross-check (tag from data, both cross-check).
    assert audit.IR_TRIAGE['Fe I'][:2] == ('B', 'IR-CROSSCHECK')
    assert audit.IR_TRIAGE['Fe II'][:2] == ('B', 'IR-CROSSCHECK')


def test_zn_is_non_canonical_only():
    assert 'Zn' not in audit.CANONICAL_SYMBOLS  # Zn is NOT a canonical element symbol


@pytest.mark.skipif(not (OUT / 'ir_triage.csv').exists(), reason="artifact not generated")
def test_ir_triage_artifact_27_species_fe_split_no_zn():
    df = pd.read_csv(OUT / 'ir_triage.csv')
    assert len(df) == 27
    species = set(df['species'])
    assert species == CANONICAL_27
    assert {'Fe I', 'Fe II'} <= species and 'Zn' not in species
    # Every IR-PRIMARY / IR-CROSSCHECK species must actually hold NIR data.
    ab = df[df['group'].isin(['A', 'B'])]
    assert ab['nir_held'].all(), "Group A/B species tagged IR but no NIR held → defect"
    # Fe I is the sigma ruler (most NIR lines); Fe II is sparse by comparison.
    fe1 = df.loc[df['species'] == 'Fe I', 'solar_nir_lines'].iloc[0]
    fe2 = df.loc[df['species'] == 'Fe II', 'solar_nir_lines'].iloc[0]
    assert fe1 > 1000 and fe2 < fe1 / 10, "Fe I should dwarf Fe II in IR line count"


@pytest.mark.skipif(not (OUT / 'non_canonical_holdings.csv').exists(), reason="artifact not generated")
def test_zn_reported_as_non_canonical_holding():
    df = pd.read_csv(OUT / 'non_canonical_holdings.csv')
    assert 'Zn' in set(df['element'])
    assert not df['canonical'].any(), "non-canonical report must contain only non-canonical species"


@pytest.mark.skipif(not (OUT / 'groupA_depth_confirm.csv').exists(), reason="artifact not generated")
def test_groupA_named_diagnostics_all_covered_with_depth():
    df = pd.read_csv(OUT / 'groupA_depth_confirm.csv')
    assert df['covered'].all(), "a named Group-A IR diagnostic is missing from the held list"
    # The S I 1.045 um channel must be materially deeper than the weak optical S I 6757.
    s_triplet = df[df['diagnostic'].str.startswith('S I 1.045')]
    assert s_triplet['max_central_depth'].max() > 0.10
