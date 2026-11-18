import numpy as np
import pytest
from pipeline.rate_semantics import gross_rate_relative_difference

def test_effective_native_se_weights_cannot_be_compared_as_gross_rates():
    for reference,replay in [('native-preconditioned-effective','physical-gross'),('physical-gross','unknown'),('native-preconditioned-effective','native-preconditioned-effective')]:
        with pytest.raises(ValueError,match='unpreconditioned'):
            gross_rate_relative_difference([1.],[1000.],reference_semantics=reference,replay_semantics=replay)
    assert gross_rate_relative_difference([0.,1.],[0.,1.25],reference_semantics='physical-gross',replay_semantics='physical-gross')==.25
    with pytest.raises(ValueError,match='support'):
        gross_rate_relative_difference([0.,1.],[1.,1.],reference_semantics='physical-gross',replay_semantics='physical-gross')
