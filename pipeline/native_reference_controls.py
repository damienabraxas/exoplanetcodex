"""Explicit one-factor native reference controls; no undeclared substitution."""
import copy


def validate_observer_plan(plan,qualified):
    if qualified['controls']['nq_fac']!=2.:raise ValueError('Require qualified NQFAC2 base')
    factor=plan['controls']['nq_fac']
    if isinstance(factor,bool)or factor not in [2.,4.]:raise ValueError('Only predeclared NQFAC2/4 observer controls allowed')
    candidate=copy.deepcopy(plan);candidate['controls']['nq_fac']=2.;candidate['scope']=qualified['scope']
    if candidate!=qualified:raise ValueError('Native observer refinement changed more than NQFAC and scope')
    return factor


def validate_shared_radiation_plan(plan,qualified_frequency4):
    if qualified_frequency4['controls']['nq_fac']!=4. or 'decouple_continuum'in qualified_frequency4['controls']:raise ValueError('Require default-decoupled qualified NQFAC4 base')
    candidate=copy.deepcopy(plan)
    if candidate['controls'].pop('decouple_continuum',None)is not False:raise ValueError('Require explicit shared-radiation control')
    candidate['scope']=qualified_frequency4['scope']
    if candidate!=qualified_frequency4:raise ValueError('Shared-radiation control changed another physical/numerical axis')
