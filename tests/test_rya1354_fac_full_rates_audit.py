"""Complete full-valence rates, mapped-state provenance and zero corrections."""
from pathlib import Path
import sys
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import rya1354_fac_full_rates_audit as audit


def test_actual_eight_full_space_rate_cases_and_paired_comparison():
    r=audit.summarize(ROOT/'data/audit/rya1354')
    assert len(r['cases'])==8 and r['accepted_lifetime'] is None and r['scientific_admission']=='HOLD'
    for case in r['cases'].values():
        v=case['audit'];assert v['channels']==13
        assert v['printed_rows_verified']+len(v['writer_suppressed'])==13
        for cfg in case['configuration_state'].values():
            assert cfg['groups']==31 and cfg['configurations']==189
            assert cfg['mode']==-1 and cfg['max_abs_delta']==cfg['max_abs_shift']==0
    assert set(r['paired_SD_full_comparison_coherent100'])=={'ground','equal'}
    assert all(set(v)=={'B','C'} for v in r['paired_SD_full_comparison_coherent100'].values())


def test_unqualified_full_mapping_hash_rejected(monkeypatch):
    monkeypatch.setattr(audit,'AUDIT_SHA','0'*64)
    with pytest.raises(ValueError,match='mapping not qualified'):audit.summarize(ROOT/'data/audit/rya1354')
