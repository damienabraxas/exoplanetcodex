"""Convergence-only variant must be explicit and must preserve native updates."""
from pathlib import Path
import json
import sys
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_strict_scf import ENERGY_STOP,strict_source
from rya1354_two_method_summary import capsule_files
from rya1354_error_resolution import summarize,qualified_scf
from rya1354_triplet_correlation import require_triplet_quintet_order


def test_strict_variant_changes_only_the_declared_stop():
    files=capsule_files(ROOT/'data/audit/rya1354/strict_controls_n34_v1.tar.gz')
    original=files['toolchains/grasp-strict-scf-v1/scf.original.f90'].decode()
    changed=files['toolchains/grasp-strict-scf-v1/scf.f90'].decode()
    assert strict_source(original)==changed
    assert sum(a!=b for a,b in zip(original.splitlines(),changed.splitlines()))==1
    assert len(original.splitlines())==len(changed.splitlines())
    for name in ('rwfn.out','rmix.out'):
        assert files['runs/strict_control_original_v1/'+name]==files['runs/strict_control_variant_v1/'+name]
    with pytest.raises(ValueError):
        strict_source('unexpected revision')
    with pytest.raises(ValueError):
        strict_source(original+original)


def test_native_replay_resolves_numerics_without_claiming_accuracy(tmp_path):
    out=tmp_path/'resolution.json'
    summarize(out)
    assert out.read_bytes()==(ROOT/'data/audit/rya1354/source_curation_v2/error_resolution_v3.json').read_bytes()
    d=json.loads(out.read_text())
    for r in d['numerical_continuations']:
        assert r['native_SCF_diagnostics']['all_reported_functions_below_ACC']
        assert len(r['channels'])==16
    assert len(d['lifetime']['channels'])==13
    assert d['lifetime']['triplet_SCF_diagnostics']['all_reported_functions_below_ACC']
    assert d['lifetime_n4']['triplet_SCF_diagnostics']['all_reported_functions_below_ACC']
    assert len(d['lifetime_n4']['channels'])==13
    assert len(d['numerical_continuations'])==3
    assert d['finer_grid']['native_N']==439
    assert d['finer_grid']['diagnostics']['native_function_ACC']==.03125**6
    assert set(d['ground_n5_hamiltonian'])=={'DC','Breit_static_limit'}
    for case in d['ground_n5_hamiltonian'].values():
        assert len(case['channels'])==16
        assert {r['operator'] for r in case['channels']}=={'M1','E2','M3','E4'}
    assert d['calibrated_sigma'] is None and d['covariance'] is None
    assert not d['adoption_allowed']


def test_unqualified_or_reversed_model_cannot_enter_lifetime():
    f=capsule_files(ROOT/'data/audit/rya1354/strict_controls_n34_v1.tar.gz')
    root='runs/strict_n3_v1'
    broken=dict(f)
    broken[root+'/rmcdhf.stderr']+=b'\nMaximum iterations in SCF Exceeded.\n'
    with pytest.raises(ValueError):
        qualified_scf(broken,root)
    triplet=f['runs/triplet_correlation_v1/to_5S/upper.lsj.lbl'].decode()
    quintet=f['runs/triplet_correlation_v1/to_5S/lower.lsj.lbl'].decode()
    with pytest.raises(ValueError,match='energy ordering'):
        require_triplet_quintet_order(triplet,quintet)
    with pytest.raises(ValueError,match='identity'):
        require_triplet_quintet_order(triplet.replace('3s_3S','3s_5S'),quintet)


def test_custom_grid_does_not_silently_recompute_native_accy():
    f=capsule_files(ROOT/'data/audit/rya1354/strict_extensions_v3.tar.gz')
    with pytest.raises(ValueError,match='native ACCY prompt'):
        qualified_scf(f,'runs/strict_n3_halfgrid_v1')


def test_split_native_capsule_requires_all_receipted_parts():
    with pytest.raises(ValueError,match='missing or changed'):
        capsule_files(ROOT/'data/audit/rya1354/ground_n5_hamiltonian_v1_parts.DC.tar.gz')
