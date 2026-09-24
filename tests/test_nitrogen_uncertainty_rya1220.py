"""N-specific coupling/admission checks against the shared contract, not a second recipe."""
import math
import pytest
from pipeline import uncertainty_contract as uc
from pipeline.nitrogen_uncertainty import molecular_coupling, attach_budget


def coupling(**kw):
    return molecular_coupling({'CN': 8.}, {'CN': 8.1}, {'CN': 7.9},
        {'CN': 7.98}, {'CN': 8.02}, carbon_step=.1, oxygen_step=.1,
        source='synthetic CN response control', **kw)


def test_missing_c_o_covariance_keeps_measured_responses_but_holds_total():
    term = coupling()
    assert term['state'] == 'HOLD' and term['sigma_dex'] is None
    assert term['evidence']['responses']['C']['jacobian'] == pytest.approx(-1)
    p = attach_budget(dict(element='N', star='solar', tier='ALL'), ['CN'], [term],
                      covariance=[], covariance_source='no resolved numerical components',
                      assumptions='C/O covariance unknown')
    assert p['uncertainty']['schema'] == uc.SCHEMA
    assert p['sigma_reported'] is None
    assert set(uc.COMPONENTS) == {c['name'] for c in p['uncertainty']['components']}
    assert uc.publication_problems(p)


def test_c_o_cross_term_propagates_through_canonical_arithmetic():
    term = coupling(abundance_covariance=[[.01,.003],[.003,.0025]],
                    covariance_source='synthetic correlated C/O',
                    response_assessment='Synthetic linear response, analytically exact')
    assert term['sigma_dex'] == pytest.approx(math.sqrt(.01+.04*.0025-2*.2*.003))
    independent = coupling(abundance_covariance=[[.01,0],[0,.0025]],
                    covariance_source='synthetic independent C/O',
                    response_assessment='Synthetic linear response')
    assert term['sigma_dex'] < independent['sigma_dex']


def test_changed_indicator_pool_is_not_silently_restricted():
    with pytest.raises(uc.UncertaintyError, match='pools differ'):
        molecular_coupling({'CN':8}, {'other':8.1}, {'CN':7.9}, {'CN':8}, {'CN':8},
                           carbon_step=.1, oxygen_step=.1, source='synthetic')


def test_invalid_c_o_covariance_is_rejected():
    with pytest.raises(uc.UncertaintyError, match='positive semidefinite'):
        coupling(abundance_covariance=[[.01,.02],[.02,.01]],
                 covariance_source='invalid control', response_assessment='linear control')
