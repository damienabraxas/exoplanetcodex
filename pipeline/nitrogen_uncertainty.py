"""Nitrogen-specific evidence adapters for the canonical RYA-587 contract.

No base component vocabulary, total, validation or differential implementation
lives here. Those belong exclusively to pipeline.uncertainty_contract.
"""
from __future__ import annotations

import math
from pipeline import uncertainty_contract as uc


def molecular_coupling(nominal, carbon_minus, carbon_plus, oxygen_minus,
                       oxygen_plus, *, carbon_step, oxygen_step, source,
                       abundance_covariance=None, covariance_source=None,
                       response_assessment=None):
    """Propagate measured CN/NH C/O responses when their covariance is established.

    Steps are numerical probes, not adopted C/O errors. Keep a HOLD with all
    measured responses when the abundance covariance or assessment is missing.
    The common pool is an indicator set, not a count of model transitions/pixels.
    """
    responses = {}
    for element, minus, plus, step in (
        ('C', carbon_minus, carbon_plus, carbon_step),
        ('O', oxygen_minus, oxygen_plus, oxygen_step),
    ):
        responses[element] = uc.paired_response(
            nominal, minus, plus, delta=step, source=source,
            parameter_source=f'Numerical A({element}) probe only; not adopted sigma')
    evidence = {'responses': responses, 'abundance_order': ['C', 'O'],
                'abundance_covariance': abundance_covariance,
                'covariance_source': covariance_source,
                'response_assessment': response_assessment,
                'pool_sha256': uc.pool_digest(nominal)}
    if abundance_covariance is None or not covariance_source or not response_assessment:
        evidence['hold_reason'] = 'C/O covariance or nonlinear/numerical response assessment unresolved'
        return uc.component('molecular_coupling', None, state='HOLD',
                            source=source, evidence=evidence)
    j = [responses[k]['jacobian'] for k in ('C', 'O')]
    sigma = math.sqrt(uc.covariance_variance(j, abundance_covariance))
    return uc.component('molecular_coupling', sigma, state='MEASURED',
                        source=source, evidence=evidence)


def attach_budget(product, indicator_ids, components, *, covariance,
                  covariance_source, assumptions):
    """Attach the shared contract to one nitrogen product; omitted terms stay HOLD."""
    if product.get('element') != 'N':
        raise ValueError('Nitrogen evidence adapter received another element')
    p = {**product, 'uncertainty_indicator_ids': list(indicator_ids)}
    scope = uc.product_scope(p, star=p['star'], indicator_ids=indicator_ids)
    p['uncertainty'] = uc.assemble(scope, components, covariance=covariance,
                                 covariance_source=covariance_source,
                                 assumptions=assumptions)
    p['sigma_reported'] = p['uncertainty']['sigma_reported']
    return p
