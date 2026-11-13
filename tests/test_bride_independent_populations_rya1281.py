import numpy as np
import pytest
from scripts.audit_bride_independent_oxygen_populations import compare_departures


def test_reference_mapping_preserves_total_without_hiding_raw_departure_discrepancy():
    n=np.array([[2.,3.],[8.,7.]]);s=np.array([[1.,2.],[9.,8.]]);b=np.array([[2.1,1.6],[.9,.9]]);before=[a.copy()for a in [n,s,b]];r=compare_departures(n,s,b,[.01,1.])
    np.testing.assert_allclose(r['reference_mapped_populations'].sum(axis=0),n.sum(axis=0));np.testing.assert_array_equal(r['delta_b'],n/s-b)
    assert np.any(r['reference_mapping_scale']!=1)
    for a,old in zip([n,s,b],before):np.testing.assert_array_equal(a,old)
    with pytest.raises(ValueError):compare_departures(n,s,b.T[:1],[.01,1.])
    with pytest.raises(ValueError):compare_departures(n,s,b*0,[.01,1.])
