import numpy as np
import pytest
from scripts.audit_bride_cn_interface import cn_partition_record,project_partition,isotope_weights


def test_partition_selects_neutral_species_and_rejects_ambiguous_records():
    row="63 '607.00' 'C N ' -7.724 2 0D0 1D0 2 6 1 7 1"
    assert cn_partition_record(row)['coefficients']==[0.,1.]
    np.testing.assert_allclose(project_partition([1.,10.,100.],[0.,1.]),[1.,10.,100.],rtol=5e-7)
    for text in [row+'\n'+row,row.replace('2 6 1 7 1','2 6 2 7 1'),row.replace('C N ','C N -')]:
        with pytest.raises(ValueError):cn_partition_record(text)
    with pytest.raises(ValueError):project_partition([0.],[1.])


def test_isotope_products_use_serialized_fractions_once():
    source='6 12 12 .9893\n6 13 13 .0107\n7 14 14 .99632\n7 15 15 .00368'
    _,weights=isotope_weights(source)
    assert weights[0]['serialized_product']==.989*.996
    assert weights[2]['serialized_product']==.989*.004
    assert sum(x['serialized_product']for x in weights)==pytest.approx(1.)
    assert weights[2]['relative_change']>0.08
    with pytest.raises(ValueError):isotope_weights(source+'\n6 12 12 .9893')
