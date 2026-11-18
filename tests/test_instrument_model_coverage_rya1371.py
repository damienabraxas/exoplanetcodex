import numpy as np
import pytest
from pipeline.instrument_model_coverage import catalog_interval_vac_nm,interval_intersections


def test_catalog_wavelength_bounds_are_not_reinterpreted_as_raw_wavenumber():
    np.testing.assert_array_equal(catalog_interval_vac_nm(404.7,1065.,'vacuum_wavenumber'),[404.7,1065.])
    air=catalog_interval_vac_nm(296.,1300.,'air');assert air[0]>296. and air[1]>1300.
    with pytest.raises(ValueError):catalog_interval_vac_nm(90.,320.,'air')
    with pytest.raises(ValueError):catalog_interval_vac_nm(90.,320.,'unspecified')


def test_setup_gaps_are_preserved_as_independent_intersections():
    assert interval_intersections([100.,500.],[[90.,110.],[300.,310.],[600.,700.]])==[[100.,110.],[300.,310.]]
