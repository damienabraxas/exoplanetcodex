"""Emit the measured/HOLD CN candidate under the merged RYA-587 contract.

This is a review artifact. It never writes data/products/solar and is expected
to fail the publication gate while transition-data, C/O covariance, blends and
other material components remain unresolved.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline import uncertainty_contract as uc


ROOT = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    ids = ['CN_AX_IR:12C14N:A-X(0-0):windows_sha256:c4460615564ca9f01a14bc0bf80da543b453b56c9b347a43205ec2130a91babb']
    product = dict(element='N', ion='I', band='NIR', instrument='iag_fts_solar_atlas',
                   holding='solar_iag', tier='ALL', selector='MOL-CN_AX_IR', route='SYNTH',
                   treatment='1D-LTE', star='solar', A=7.716)
    scope = uc.product_scope(product, star='solar', indicator_ids=ids)
    run = 'data/output/rya1220/cn_iag_stellar_v2/fits.json'
    measurement = uc.fit_measurement(.016, source=run,
        likelihood='pipeline.fit_constraint.measure_constraint; profile curvature rescaled to red_chi2=1',
        correlation_treatment='UNRESOLVED: 2211 correlated pixels; no independent residual covariance yet',
        n_pixels=2211)
    terms = [measurement]
    for name, sigma, state, evidence in [
        ('stellar.teff', 0., 'MEASURED', dict(pool_sha256=scope['pool_sha256'],
            parameter='teff_K', parameter_source='pipeline.uncertainty_stack.params_and_deltas(solar): e_teff=1 K',
            response_assessment='Both ±1 K refits return A(N)=7.716 at 0.001 dex optimizer precision',
            signed_response_dex=0., delta_minus_dex=0., delta_plus_dex=0., delta_parameter=1.)),
        ('stellar.logg', 0., 'DEFINED', dict(pool_sha256=scope['pool_sha256'],
            parameter='logg', parameter_source='solar log g exact by definition in uncertainty_stack',
            response_assessment='N/A: declared delta_p=0', delta_parameter=0.)),
        ('stellar.xi', .0955, 'MEASURED', dict(pool_sha256=scope['pool_sha256'],
            parameter='vturb_kms', parameter_source='pipeline.uncertainty_stack.params_and_deltas(solar): adopted xi allowance 0.2912 km/s',
            response_assessment='Asymmetric ±xi refits retained; central linear response used provisionally',
            signed_response_dex=.0955, delta_minus_dex=0., delta_plus_dex=.191,
            delta_parameter=.2912)),
        ('stellar.metallicity', 0., 'DEFINED', dict(pool_sha256=scope['pool_sha256'],
            parameter='feh', parameter_source='solar [Fe/H] exact by definition in uncertainty_stack',
            response_assessment='N/A: declared delta_p=0', delta_parameter=0.)),
    ]:
        terms.append(uc.component(name, sigma, state=state,
                                  source='RYA-1220 CN response run; shared RYA-587 contract',
                                  evidence=evidence))
    holds = {
        'transition_data': 'Exact molecular quantum identities and line-specific uncertainty/covariance not yet complete',
        'continuum': 'Source-normalized holding; local continuum perturbation not measured',
        'profile_ew': 'Fit residual/profile covariance unresolved',
        'pseudo_continuum': 'Molecular pseudo-continuum contribution unresolved',
        'telluric': 'Holding is declared clean, but residual uncertainty is not quantified',
        'holding_instrument': 'Repeatability/conditioning response not measured',
        'nlte': 'CN route is LTE-by-design; model-form uncertainty unresolved',
        'model_atmosphere': '3D molecular correction not established for this pool',
        'hfs_isotopes': 'Isotope sensitivity not measured',
        'blends': 'Species-specific blend decomposition unresolved',
        'molecular_coupling': 'C/O covariance absent despite measured dA/dC and dA/dO responses',
    }
    for name, reason in holds.items():
        terms.append(uc.component(name, None, state='HOLD',
            source='RYA-1220 CN response run; unresolved evidence',
            evidence={'pool_sha256': scope['pool_sha256'], 'hold_reason': reason}))
    numeric = [c for c in terms if c['sigma_dex'] is not None]
    cov = [[0.0]*len(numeric) for _ in numeric]
    for i, c in enumerate(numeric):
        cov[i][i] = c['sigma_dex']**2
    budget = uc.assemble(scope, terms, covariance=cov,
        covariance_source='Diagonal only for measured terms; unresolved shared covariance is not assumed zero',
        assumptions='Measured terms are provisional; no covariance cancellation assumed for unresolved components')
    candidate = {**product, 'uncertainty_indicator_ids': ids,
                 'uncertainty': budget, 'sigma_reported': budget['sigma_reported']}
    output = {'schema': uc.SCHEMA, 'candidate': candidate,
              'publication_problems': uc.publication_problems(candidate),
              'status': 'HOLD_REVIEW_ONLY', 'live_feed_modified': False}
    (args.out/'candidate.json').write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps({'status': output['status'], 'holds': budget['holds'],
                      'publication_problems': output['publication_problems']}, indent=2))


if __name__ == '__main__':
    main()
