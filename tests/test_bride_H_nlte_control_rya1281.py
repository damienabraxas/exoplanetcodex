import numpy as np
import pytest
from pipeline.bride_H_nlte_control import scale_initial_fractions


def test_same_positive_fraction_field_maps_to_different_totals_without_mutation():
    pop=np.array([[1.,4.],[3.,2.]]);before=pop.copy();a=scale_initial_fractions(pop,[20.,30.]);b=scale_initial_fractions(pop,[22.,27.])
    np.testing.assert_array_equal(pop,before);np.testing.assert_allclose(a.sum(axis=0),[20.,30.]);np.testing.assert_allclose(b.sum(axis=0),[22.,27.]);np.testing.assert_allclose(a/a.sum(axis=0),b/b.sum(axis=0))


def test_invalid_fraction_or_target_refused():
    for pop,total in [([[0.,1.],[1.,1.]],[1.,1.]),([[1.,1.],[1.,np.nan]],[1.,1.]),([[1.,1.],[1.,1.]],[1.]),([[1.,1.],[1.,1.]],[0.,1.])]:
        with pytest.raises(ValueError):scale_initial_fractions(pop,total)
