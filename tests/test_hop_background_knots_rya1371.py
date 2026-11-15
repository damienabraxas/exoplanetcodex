import numpy as np
import pytest
from pipeline.hop_background_knots import read_background_knots


def test_source_wavelength_knots_keep_close_edges_and_convert_once(tmp_path):
    p=tmp_path/'absdat';p.write_text('FeI source\n ILOGL=0 KVADL=0 NLATB=4\n 900.0 2511.0\n 2513.0 2624.0\n ILOGT=1\n')
    table=read_background_knots(p)[0]
    np.testing.assert_array_equal(table['wave_vac_nm'],np.array([900.,2511.,2513.,2624.])*.1)
    p.write_text('bad\n NLATB=3\n 900 900 1000\n')
    with pytest.raises(ValueError):read_background_knots(p)


@pytest.mark.parametrize('count',['-3','unknown'])
def test_unsupported_knot_count_refuses_instead_of_omitting_component(tmp_path,count):
    p=tmp_path/'absdat';p.write_text(f'FeI source\n NLATB={count}\n 900 1000 1100\n')
    with pytest.raises(ValueError):read_background_knots(p)
