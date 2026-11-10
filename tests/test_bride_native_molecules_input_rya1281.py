import numpy as np
import pytest
from scripts.run_bride_native_molecules import lte_export,native_status


def test_lte_export_clears_blank_and_populated_nlte_labels_without_changing_source_or_gf():
    data=np.array([(0.1,'',''),(-2.,'lower','upper')],dtype=[('loggf',float),('nlte_label_low','U20'),('nlte_label_up','U20')])
    old=data.copy();out=lte_export(data)
    np.testing.assert_array_equal(data,old);np.testing.assert_array_equal(out['loggf'],old['loggf'])
    assert (out['nlte_label_low']=='none').all()and(out['nlte_label_up']=='none').all()


def test_native_input_failure_refuses_even_with_zero_exit_status():
    native_status(0,'normal native output')
    for code,text in [(1,''),(0,'input.f; bad keyword: NLTE:'),(0,'INPUT: bad answer:')]:
        with pytest.raises(RuntimeError):native_status(code,text)
