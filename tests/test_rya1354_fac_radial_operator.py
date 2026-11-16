from pathlib import Path
import sys
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from rya1354_fac_radial_operator import ACTIVE,pairs,input_text,raw_rows
from rya1354_fac_radial_preflight import input_text as radial


def fixture(gauge=2,request=.3):
    rows=[]
    for (n,k),(n2,k2) in pairs():
        e=request if request else .2
        rows.append(f'{n} {k} {n2} {k2} {gauge} 0 {int(request==0)} 1 {request} {e} -1 -1.2 {e*.01:.17g} {e*.01:.17g} 2 2')
    return '\n'.join(rows)


def test_complete_ordered_valence_e1_pairs_and_frozen_radial_prefix():
    assert len(ACTIVE)==10 and len(pairs())==len(set(pairs()))==42
    assert ((3,-1),(2,1)) in pairs() and ((3,-1),(2,-2)) in pairs()
    assert all((b,a) in pairs() for a,b in pairs())
    for scale in (10,100):
        for gauge in (1,2):
            for energy in (0,.3):
                text=input_text(scale,gauge,energy)
                assert text.startswith(radial(scale)) and text.count('RawRadialE1(')==42
                assert 'Structure(' not in text and 'OptimizeRadial(' not in text
                assert f'SetTransitionOptions({gauge},0)' in text
    with pytest.raises(ValueError,match='Undeclared'):input_text(10,3,.3)
    explicit=input_text(10,2,.3,1e-8)
    assert explicit.startswith(radial(10))
    assert "SetOption('radial:awmin',1e-08)" in explicit
    assert explicit.replace("SetOption('radial:awmin',1e-08)\n",'')==input_text(10,2,.3)
    with pytest.raises(ValueError,match='floor'):input_text(10,2,.3,0)


def test_radial_rows_require_complete_finite_separate_gauge_context():
    c={'FINE_STRUCTURE_CONST':.01}
    assert len(raw_rows(fixture(),2,.3,c))==42
    with pytest.raises(ValueError,match='Incomplete'):raw_rows('\n'.join(fixture().splitlines()[:-1]),2,.3,c)
    with pytest.raises(ValueError,match='pair/gauge/context'):raw_rows(fixture()+'\n'+fixture().splitlines()[0],2,.3,c)
    with pytest.raises(ValueError,match='pair/gauge/context'):raw_rows(fixture(1),2,.3,c)
    bad=fixture().splitlines();f=bad[0].split();f[13]='.02';bad[0]=' '.join(f)
    with pytest.raises(ValueError,match='photon grid'):raw_rows('\n'.join(bad),2,.3,c)


def test_orbital_gap_is_native_energy_difference_not_an_arbitrary_frequency():
    c={'FINE_STRUCTURE_CONST':.01}
    # Use exact binary subtraction, as the native wrapper does.
    text=fixture(request=0).replace('0.2 -1 -1.2 0.002 0.002',f'{abs(-1.-(-1.2))} -1 -1.2 {abs(-1.-(-1.2))*.01:.17g} {abs(-1.-(-1.2))*.01:.17g}')
    assert len(raw_rows(text,2,0,c))==42
    with pytest.raises(ValueError,match='frequency/orbital'):raw_rows(fixture(request=0),2,0,c)


def test_rate_sensitivity_separates_photon_power_and_amplitude():
    from rya1354_fac_radial_operator_audit import factorization
    for m in (-1,2,-3):
        exponent=2*abs(m)+1
        first=dict(multipole=m,energy_au=1.,amplitude=2.,A_s1=5.)
        second=dict(multipole=m,energy_au=2.,amplitude=6.,A_s1=5.*2**exponent*9)
        result=factorization(first,second)
        assert result['energy_factor']==2**exponent and result['amplitude_squared_factor']==9
        assert result['floating_factorization_residual']==0
    first.update(amplitude=0.,A_s1=0.)
    assert factorization(first,second)['actual_rate_ratio'] is None


def test_actual_diagnostic_preserves_failure_and_all_exact_contexts():
    import json
    from rya1354_fac_radial_operator_audit import summarize
    result=summarize(ROOT/'data/audit/rya1354')
    assert result==json.loads((ROOT/'data/audit/rya1354/source_curation_v2/fac_radial_operator_v43.json').read_text())
    assert len(result['probes'])==24 and sum(v['pairs'] for v in result['probes'].values())==1008
    assert all(v['physical_orbitals_unchanged'] for v in result['probes'].values())
    assert result['retained_failure']['status']=='FAILED_RETAINED'
    assert result['retained_failure']['qualified_cases']==2
    assert result['retained_failure']['previous_successful_controls_exact']
    assert result['scientific_admission']=='HOLD' and not result['adoption_allowed']
    assert result['accepted_lifetime'] is None and result['uncertainty'] is None
    for modes in result['gauges'].values():
        assert all(len(rows)==42 for rows in modes.values())
        assert all(row['orbital_gap_over_photon']==1 for row in modes['orbital_gap'])
