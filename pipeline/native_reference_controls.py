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


def validate_shared_refinement_plan(plan,qualified_shared,*,target_factor=8.):
    """Only predeclared consecutive shared-field NQFAC4->8 or8->16 refinements."""
    base_factor=qualified_shared['controls']['nq_fac']
    if (base_factor,target_factor)not in [(4.,8.),(8.,16.)]or qualified_shared['controls'].get('decouple_continuum')is not False:raise ValueError('Require the preceding qualified shared-field frequency case')
    candidate=copy.deepcopy(plan)
    if isinstance(candidate['controls']['nq_fac'],bool)or candidate['controls']['nq_fac']!=target_factor:raise ValueError('Require the explicit predeclared frequency factor')
    candidate['controls']['nq_fac']=base_factor;candidate['scope']=qualified_shared['scope']
    if candidate!=qualified_shared:raise ValueError('Shared frequency refinement changed another physical/numerical axis')
