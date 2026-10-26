"""Mixed TEMP blocks, exact coefficient preservation and fail-closed conversion."""
from dataclasses import replace
import numpy as np
import pytest
from pipeline.multi_atom import read_multi_atom, to_lightweaver

SOURCE = """O
8.70 15.999
3 1 0 0
0.0 5.0 'lower' 1
12000.0 7.0 'upper' 1
110000.0 4.0 'ion' 2
2 1 0.4 9 5.0 0.0 0 3.89e7 453.234 0.0 VOIGT
GENCOL
TEMP
3 2000 5000 10000
CE
2 1 1e-9 2e-9 3e-9
TEMP
4 1000 2000 3000 4000
CH
2 1 1e-10 2e-10 3e-10 4e-10
END
"""

def parse(tmp_path, text=SOURCE):
    p=tmp_path/'atom.synthetic';p.write_text(text)
    return read_multi_atom(p)


def test_each_collision_retains_its_active_temperature_grid_and_all_rates(tmp_path):
    a=parse(tmp_path)
    np.testing.assert_array_equal(a.collisions[0].temperature,[2000,5000,10000])
    np.testing.assert_array_equal(a.collisions[1].temperature,[1000,2000,3000,4000])
    np.testing.assert_array_equal(a.collisions[1].rates,[1e-10,2e-10,3e-10,4e-10])


@pytest.mark.parametrize('old,new', [('CH\n','MISSING_PHYSICS\n'),
    ('1e-10 2e-10 3e-10 4e-10','1e-10 2e-10 3e-10'),
    ('1e-10 2e-10 3e-10 4e-10','1e-10 2e-10 3e-10 4e-10 5e-10'),
    ('1000 2000 3000 4000','1000 3000 2000 4000'),
    ('1e-10 2e-10 3e-10 4e-10','1e-10 -2e-10 3e-10 4e-10'),
    ('END\n',''),('END\n','END\nCE\n2 1 1 2 3 4\n')])
def test_malformed_collision_data_is_never_skipped_or_truncated(tmp_path,old,new):
    with pytest.raises(ValueError):parse(tmp_path,SOURCE.replace(old,new))


def test_native_conversion_preserves_abo_radiative_width_and_per_rate_grid(tmp_path):
    pytest.importorskip('lightweaver')
    a=parse(tmp_path);native=to_lightweaver(a)
    line=native.lines[0]
    assert type(line.broadening.elastic[0]).__name__=='VdwBarklem'
    np.testing.assert_allclose(line.broadening.elastic[0].vals,[453.,.234],rtol=1e-12,atol=0)
    assert line.broadening.natural[0].gamma==3.89e7
    for src,dst in zip(a.collisions,native.collisions):
        np.testing.assert_array_equal(dst.temperature,src.temperature)
        np.testing.assert_allclose(dst.rates,src.rates*1e-6,rtol=1e-15,atol=0)


@pytest.mark.parametrize('field,value',[('gStark',1.2e-8),('gVdw',-7.),('profile','GAUSS')])
def test_unvalidated_broadening_or_profile_cannot_be_silently_replaced(tmp_path,field,value):
    pytest.importorskip('lightweaver')
    a=parse(tmp_path);a.lines[0]=replace(a.lines[0],**{field:value})
    with pytest.raises(ValueError):to_lightweaver(a)


def test_unmapped_rate_cannot_disappear_during_conversion(tmp_path):
    pytest.importorskip('lightweaver')
    a=parse(tmp_path);a.collisions[0]=replace(a.collisions[0],kind='UNMAPPED')
    with pytest.raises(ValueError,match='Unsupported'):to_lightweaver(a)


def test_fixed_radiative_rates_are_not_skipped(tmp_path):
    with pytest.raises(ValueError,match='fixed radiative'):
        parse(tmp_path,SOURCE.replace('3 1 0 0','3 1 0 1'))


def test_explicit_source_core_spacing_is_preserved(tmp_path):
    pytest.importorskip('lightweaver')
    a=parse(tmp_path,SOURCE.replace('0.4 9 5.0 0.0','0.4 251 100.0 3.0'))
    native=to_lightweaver(a)
    assert native.lines[0].quadrature.qCore == 3.
    assert native.lines[0].quadrature.qWing == 100.
    assert native.lines[0].quadrature.Nlambda == 251


def test_multi_explicit_continuum_keeps_table_points_beyond_energy_edge(tmp_path):
    """Source MULTI FREQC does not clip the absolute table to energy-derived FRQ0."""
    pytest.importorskip('lightweaver')
    from pipeline.multi_continuum import MultiExplicitContinuum
    from lightweaver.atomic_model import AtomicModel,AtomicLevel
    from lightweaver.atomic_table import Element
    c=MultiExplicitContinuum(j=1,i=0,wavelengthGrid=[100.,105.,110.],
                             alphaGrid=[1e-22,2e-22,3e-22])
    a=AtomicModel(element=Element(8),levels=[AtomicLevel(E=0.,g=1.,label='fixture',stage=0),
        AtomicLevel(E=120000.,g=1.,label='fixtureIon',stage=1)],lines=[],continua=[c],collisions=[])
    assert c.energy_edge_nm < 100.
    assert c.lambda0 == c.energy_edge_nm  # energy normalization is retained
    assert c.lambdaEdge == 110.
    np.testing.assert_array_equal(c.wavelength(),[100.,105.,110.])
    np.testing.assert_allclose(c.alpha(np.array([100.,105.,110.])),[1e-22,2e-22,3e-22],rtol=1e-15,atol=0)
    # Interpolation support is the source table, not silently extended cross-sections.
    np.testing.assert_array_equal(c.alpha(np.array([99.,111.])),[0.,0.])


@pytest.mark.parametrize('qmax',[-1.,10.])
def test_explicit_continuum_format_is_preserved_and_analytic_form_refused(tmp_path,qmax):
    text=SOURCE.replace('3 1 0 0','3 1 1 0').replace('GENCOL',
        f'3 1 1e-18 3 {qmax}\n1100 1e-18\n1050 2e-18\n1000 3e-18\nGENCOL')
    if qmax>=0:
        with pytest.raises(ValueError,match='analytic MULTI'):parse(tmp_path,text)
        return
    a=parse(tmp_path,text)
    np.testing.assert_array_equal(a.continua[0].wavelength,[1100.,1050.,1000.])
