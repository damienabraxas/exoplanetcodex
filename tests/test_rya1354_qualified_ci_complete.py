"""Replay native matched DC/Breit results and refuse model/input substitution."""
import json
from pathlib import Path
import sys
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_qualified_ci_audit import summarize,RUN


def test_complete_matched_native_ci_and_hamiltonian_sensitivity():
    base=ROOT/'data/audit/rya1354'
    result=summarize(base,base/'qualified_ci_complete_v17')
    assert result==json.loads((base/'source_curation_v2/qualified_ci_v17.json').read_text())
    assert all(len(c['channels'])==13 for c in result['cases'].values())
    assert result['cases']['DC']['model']['B']['model_one_photon_lifetime_ns']==pytest.approx(2.28448374775)
    assert result['cases']['Breit']['model']['C']['model_one_photon_lifetime_ns']==pytest.approx(1.71254016485)
    assert not result['adoption_allowed'] and result['accepted_lifetime'] is result['uncertainty'] is None


@pytest.mark.parametrize('change',['hamiltonian','source_wave','transition_mixing'])
def test_ci_replay_refuses_substituted_physics_and_inputs(monkeypatch,change):
    import rya1354_qualified_ci_audit as audit
    original=audit.load_parts
    def changed(path):
        files,manifest=original(path)
        if path.name=='qualified_ci_complete_v17':
            if change=='hamiltonian':
                files[RUN+'Breit/ground/rci.stdin']=files[RUN+'DC/ground/rci.stdin']
            elif change=='source_wave':
                files[RUN+'DC/triplet/state.w']=files[RUN+'DC/ground/state.w']
            else:
                files[RUN+'DC/to_ground/upper.cm']=files[RUN+'DC/ground/state.cm']
        return files,manifest
    monkeypatch.setattr(audit,'load_parts',changed)
    with pytest.raises(ValueError,match='Hamiltonian/root/mixing differs|matched source input differs|transition state differs'):
        audit.summarize(ROOT/'data/audit/rya1354',ROOT/'data/audit/rya1354/qualified_ci_complete_v17')
