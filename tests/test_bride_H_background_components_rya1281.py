import numpy as np
import pytest
from scripts.check_bride_H_background_components import pair_tables,interpolate_pair


def test_source_table_axes_and_bilinear_interpolation():
    text='double lambdaFFPlus[2]={0,10};double tempFFPlus[2]={100,200};double kappaFFPlus[4]={/* ignored9 */0,10,20,30};'
    tables=pair_tables(text);np.testing.assert_allclose(interpolate_pair([5],[150],tables),[[15]])
    with pytest.raises(ValueError):pair_tables(text+text)
    with pytest.raises(ValueError):interpolate_pair([10],[150],tables)


def test_nonfinite_interpolation_coordinates_refused():
    tables=(np.array([0.,10.]),np.array([100.,200.]),np.ones((2,2)))
    for wave,temp in [([np.nan],[150]),([5],[np.inf])]:
        with pytest.raises(ValueError):interpolate_pair(wave,temp,tables)
