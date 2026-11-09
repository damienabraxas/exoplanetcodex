"""Actual compiled complete-coupling controls, independently replayed from arrays."""
from pathlib import Path
import numpy as np
from tests.test_bride_transition_operator_rya1281 import discrepancy

ROOT=Path(__file__).resolve().parents[1]/'data/results/rya1281'
ROUND=4096*np.finfo(float).eps


def test_full_operator_preserves_source_rates_and_differs_from_pairwise_map():
    with np.load(ROOT/'solver_acceleration_v1/rya1281-mali-operator-audit-v3/operator_audit.npz')as pair,np.load(ROOT/'full_coupling_operator_v1/rya1281-mali-operator-audit-v4/operator_audit.npz')as full:
        n=pair['populations'];np.testing.assert_array_equal(n,full['populations'])
        for i in range(3):
            key=str(i);r=pair['physical_'+key]
            np.testing.assert_array_equal(r,full['physical_'+key])
            assert discrepancy(full['gamma_'+key],r,n)<ROUND
        assert not np.array_equal(full['gamma_1'],pair['gamma_1'])
        assert not np.array_equal(full['gamma_1'],full['physical_1'])


def test_full_operator_repairs_nonlte_overlap_without_changing_physical_rates():
    with np.load(ROOT/'solver_acceleration_v1/rya1281-continuum-overlap-v3/operators.npz')as native,np.load(ROOT/'full_coupling_operator_v1/rya1281-continuum-overlap-v6/operators.npz')as full:
        for i in range(9):
            key=str(i);n=native['populations_'+key];r=native['physical_'+key]
            np.testing.assert_array_equal(n,full['populations_'+key])
            np.testing.assert_array_equal(r,full['physical_'+key])
            assert discrepancy(full['operator_'+key],r,n)<ROUND
        assert discrepancy(native['operator_7'],native['physical_7'],native['populations_7'])>1e-6


def test_compiled_linear_solve_agrees_with_population_scaled_coordinate_system():
    with np.load(ROOT/'full_coupling_operator_v1/rya1281-population-linear-audit-v1/linear_solutions.npz')as d:
        assert np.all(d['compiled']>0)and np.all(d['numpy_relative']>0)
        assert np.max(abs(d['compiled']/d['numpy_relative']-1))<ROUND


def test_cap_response_changes_only_top_preconditioner_and_preserves_physics():
    base=ROOT/'full_coupling_operator_v1/rya1281-mali-operator-audit-v4/operator_audit.npz'
    cap=ROOT/'solar_full_cap_v1/rya1281-mali-operator-audit-v5/operator_audit.npz'
    with np.load(base)as b,np.load(cap)as c:
        np.testing.assert_array_equal(b['populations'],c['populations'])
        for i in range(3):
            key=str(i)
            np.testing.assert_array_equal(b['physical_'+key],c['physical_'+key])
            np.testing.assert_array_equal(b['gamma_'+key][:,:,1:],c['gamma_'+key][:,:,1:])
            assert discrepancy(c['gamma_'+key],c['physical_'+key],c['populations'])<ROUND
        assert not np.array_equal(b['gamma_1'],c['gamma_1'])
