from types import SimpleNamespace
import numpy as np
import pytest
from pipeline.native_strength_precision import strength_factors


def test_consumed_strength_precision_is_source_determined_and_identity_bound(tmp_path):
    p=tmp_path/'atom.txt';p.write_text('i j lam_air\n0 0 1 8446 0 -0.463 1 1 0 9 5 0 0\n')
    factors,rows=strength_factors([SimpleNamespace(i=0,j=1,f=.1148)],[3.,1.],p)
    np.testing.assert_allclose(factors,[.9998546189411687],rtol=1e-15)
    assert rows[0]['consumed_loggf']==-.463
    with pytest.raises(ValueError):strength_factors([SimpleNamespace(i=0,j=2,f=.1148)],[3.,1.,1.],p)
