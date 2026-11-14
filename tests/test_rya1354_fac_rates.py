"""Explicit gauge/operator inventory from actual mapped FAC levels."""
from pathlib import Path
import sys
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_native_parts import load_parts
from rya1354_fac_recoupling_audit import summarize
from rya1354_fac_oxygen_levels import inspect_inventory
from rya1354_fac_rates import rate_input


def test_native_mapping_drives_thirteen_explicit_relativistic_queries():
    base=ROOT/'data/audit/rya1354';report=summarize(base)
    files,_=load_parts(base/'fac_coherent_complete_v18')
    p='runs/fac_coherent_controls_v18/coherent100/'
    levels=inspect_inventory(files[p+'o.lev'].decode())
    for gauge in (1,2):
        text,requests=rate_input(.03,100,gauge,report['cases']['coherent100']['recoupling'],list(reversed(levels)))
        assert len(requests)==13 and all(r['upper']==6 for r in requests)
        assert f'SetTransitionOptions({gauge},0)' in text and 'SetTransitionCut(0,0)' in text
        assert {r['operator'] for r in requests}=={'E1','M2','E3','M1','E2','M3'}
        assert text.count("TransitionTable('rates.tr.b'")==6
        assert text.count("PrintTable('rates.tr.b','rates.tr',1)")==1


def test_extra_lower_state_cannot_be_silently_omitted():
    base=ROOT/'data/audit/rya1354';report=summarize(base)
    files,_=load_parts(base/'fac_coherent_complete_v18')
    levels=inspect_inventory(files['runs/fac_coherent_controls_v18/coherent100/o.lev'].decode())
    levels[7]['energy_eV']=levels[6]['energy_eV']-.01
    with pytest.raises(ValueError,match='lower-state inventory differs'):
        rate_input(.03,100,2,report['cases']['coherent100']['recoupling'],levels)


def test_archived_raw_inventory_retains_unresolved_channels():
    from rya1354_fac_rate_inventory import summarize_inventory
    report=summarize_inventory(ROOT/'data/audit/rya1354')
    assert len(report['cases'])==4 and report['accepted_lifetime'] is None
    for case in report['cases'].values():
        assert case['queries']==13 and case['finite_native_rates']==11
        assert {(r['lower'],r['multipole'],r['status']) for r in case['unresolved']}=={(0,-3,-1),(3,-3,-1)}
        assert case['lifetime_ns'] is None


def test_raw_missing_wrong_gauge_and_zero_filled_refusal_rejected():
    from rya1354_fac_rate_inventory import inspect_raw
    files,_=load_parts(ROOT/'data/audit/rya1354/fac_rates_complete_v21')
    text=files['runs/fac_model_rates_v21/coherent100_B/raw.tsv'].decode()
    import json
    receipt=json.loads(files['runs/fac_model_rates_v21/receipt.json'])
    requests=receipt['cases']['coherent100_B']['requests']
    for bad in ('\n'.join(text.splitlines()[:-1]),text.replace('nan nan nan','0 0 0')):
        with pytest.raises(ValueError):inspect_raw(bad,requests,2)
    with pytest.raises(ValueError,match='gauge'):inspect_raw(text,requests,1)
