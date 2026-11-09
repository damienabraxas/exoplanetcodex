"""Independent numerical checks of actual compiled operator-control artifacts."""
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]/'data/results/rya1281/solver_acceleration_v1'
ROUND=4096*np.finfo(float).eps


def discrepancy(operator,physical,n):
    physical_flow=np.einsum('ijd,jd->id',physical,n)
    operator_flow=np.einsum('ijd,jd->id',operator,n)
    gross=np.einsum('ijd,jd->id',abs(physical),n)
    return np.max(abs(operator_flow-physical_flow)/gross)


def test_compiled_nonlte_overlap_discriminates_single_continua_and_operator_fix():
    with np.load(ROOT/'rya1281-continuum-overlap-v3/operators.npz')as native,np.load(ROOT/'rya1281-continuum-overlap-v5/operators.npz')as fixed:
        for i in range(9):
            key=str(i);n=native['populations_'+key];r=native['physical_'+key]
            np.testing.assert_array_equal(n,fixed['populations_'+key])
            np.testing.assert_array_equal(r,fixed['physical_'+key])
            assert discrepancy(fixed['operator_'+key],r,n)<ROUND
            if i==7:assert discrepancy(native['operator_'+key],r,n)>1e-6
            else:assert discrepancy(native['operator_'+key],r,n)<ROUND


def test_full_source_fix_preserves_physical_rates_and_population_state_exactly():
    old=ROOT.parent/'solar_solver_controls_v1/rya1281-mali-operator-audit-v1/operator_audit.npz'
    with np.load(old)as native,np.load(ROOT/'rya1281-mali-operator-audit-v3/operator_audit.npz')as fixed:
        n=native['populations'];np.testing.assert_array_equal(n,fixed['populations'])
        for i in range(3):
            key=str(i);r=native['physical_'+key]
            np.testing.assert_array_equal(r,fixed['physical_'+key])
            assert discrepancy(fixed['gamma_'+key],r,n)<ROUND
        # A physical-matrix fallback would also preserve the action, but would
        # not exercise the intended approximate MALI operator.
        assert not np.array_equal(fixed['gamma_1'],fixed['physical_1'])
        assert discrepancy(native['gamma_1'],native['physical_1'],n)>1e-3
