"""Review-only migration of existing N evidence to the shared RYA-587 contract.

Never writes the live feed or calls synthesis. Preserves original measurements,
matches accepted atomic lines through the shared wavelength+EP matcher, and
retains missing covariance/fit evidence as HOLD. Outputs require a new directory.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import pandas as pd
from pipeline import uncertainty_contract as uc
from pipeline.line_match import match, require_resolved
from pipeline.gf_grades import nist_sigma_dex
from pipeline.nitrogen_uncertainty import attach_budget
from pipeline.uncertainty_stack import params_and_deltas


def file_evidence(path):
    return {'path': str(path.relative_to(ROOT)),
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    feed = ROOT/'data/products/solar/N.json'
    original = json.loads(feed.read_text())
    canonical_path = ROOT/'data/linelists/canonical_gf.csv'
    canonical = pd.read_csv(canonical_path, low_memory=False)
    canonical = canonical[(canonical.key_z == 7) & (canonical.ion == 1)].reset_index(drop=True)
    inv_path = ROOT/'data/reference/nitrogen_method_rya1220/indicator_inventory.json'
    inventory = json.loads(inv_path.read_text())
    _, deltas = params_and_deltas('solar')
    candidates, rows = [], []
    for index, product in enumerate(original['quarantine']):
        p = {**product, 'star': 'solar'}
        artifact = ROOT/product['provenance']['copied_to']
        evidence = {'product_artifact': file_evidence(artifact),
                    'historical_sigma_stat': p.get('sigma_stat'),
                    'historical_sigma_syst': p.get('sigma_syst'),
                    'historical_stat_basis': p.get('stat_basis')}
        is_molecular = p['selector'].startswith('MOL-')
        molecular_definition = None
        if is_molecular:
            key = p['selector'].removeprefix('MOL-')
            definitions = [r for r in inventory['indicators']
                           if r['key'] == key and r['holding'] == p['holding']]
            if len(definitions) != 1:
                raise ValueError(f'No unique molecular indicator definition: {key}/{p["holding"]}')
            molecular_definition = definitions[0]
            # Stable band + exact window identity, not a claim of resolved rotational lines.
            payload = {'key': key, 'kind': molecular_definition['kind'],
                       'windows_air_A': molecular_definition['windows_air_A']}
            digest = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
            ids = [f'molecular-indicator:{key}:windows_sha256:{digest}']
            evidence['indicator_definition'] = molecular_definition
            evidence['definition_source'] = file_evidence(inv_path)
        else:
            lines_path = artifact.with_name(artifact.name.replace('_products.csv', '_lines.csv'))
            if not lines_path.exists() and p['treatment'] == '1D-LTE':
                lines_path = artifact.with_name(artifact.name.replace('_products.csv', '_1D-LTE_lines.csv'))
            lines = pd.read_csv(lines_path)
            accepted = lines[lines.in_aggregate.astype(str).str.lower().eq('true')]
            joined = match(accepted.wavelength_air_A, canonical.wavelength_air_A,
                           want_ep=accepted.ep_eV, src_ep=canonical.excitation_potential_eV,
                           require_ep=True)
            indices = require_resolved(joined, what=str(lines_path), species='N I')
            matched = canonical.iloc[indices]
            ids = matched.physical_id.astype(str).tolist()
            if len(ids) != int(p['n_lines']):
                raise ValueError('Historical product count differs from its accepted physical pool')
            evidence['line_artifact'] = file_evidence(lines_path)
            evidence['canonical_artifact'] = file_evidence(canonical_path)
            evidence['line_data'] = [
                {'physical_id': str(r.physical_id), 'loggf': float(r.log_gf),
                 'grade': str(r.nist_grade), 'sigma_dex': nist_sigma_dex(r.nist_grade),
                 'source': str(r.loggf_reference)} for r in matched.itertuples()]
            evidence['accepted_abundances'] = accepted.abundance.tolist()

        reasons = {
            'measurement': 'Existing statistic retained; correlated pixels/lines and aggregation likelihood not established under RYA-587',
            'transition_data': ('Molecular transition identities/uncertainties and common-source covariance unresolved'
                                if is_molecular else 'Per-line canonical grades recovered; aggregate response weights and common-source covariance unresolved'),
            'continuum': 'Source normalization is not a measured continuum error',
            'profile_ew': 'Separation from fit likelihood and correlated residuals unresolved',
            'pseudo_continuum': 'Blended spectrum; pseudo-continuum effect not quantified',
            'telluric': 'Verified correction or generic legacy allowance does not quantify this product residual',
            'holding_instrument': 'Matched-pool repeatability and conditioning response unresolved',
            'nlte': 'Treatment recorded; formation-model uncertainty/applicability unresolved',
            'model_atmosphere': 'Model discrepancy is not a sourced uncertainty distribution',
            'hfs_isotopes': 'N isotope/component applicability and sensitivity not yet established',
            'blends': 'Species-specific blend/continuum separation unresolved',
            'molecular_coupling': 'C/O covariance and exact-indicator responses unresolved; CN also blends N I',
        }
        terms = []
        for name in uc.COMPONENTS:
            parameter = {'stellar.teff': 'teff_K', 'stellar.logg': 'logg',
                         'stellar.xi': 'vturb_kms', 'stellar.metallicity': 'feh'}.get(name)
            source = 'RYA-1220 existing-evidence audit: ' + str(artifact.relative_to(ROOT))
            ev = {'pool_sha256': uc.pool_digest(ids)}
            if parameter and deltas[parameter] == 0:
                terms.append(uc.component(name, 0., state='DEFINED',
                    source='pipeline.uncertainty_stack.params_and_deltas(solar): declared exact solar parameter',
                    evidence={**ev, 'parameter': parameter, 'adopted_sigma': 0.}))
                continue
            ev['hold_reason'] = reasons.get(name, 'Exact-pool symmetric adopted-sigma responses not yet supplied')
            if parameter:
                ev.update(parameter=parameter, adopted_sigma=deltas[parameter],
                          parameter_source='pipeline.uncertainty_stack.params_and_deltas(solar)')
            if name in {'measurement', 'transition_data'}:
                ev['existing_evidence'] = evidence
            terms.append(uc.component(name, None, state='HOLD', source=source, evidence=ev))
        numeric = [c for c in terms if c['sigma_dex'] is not None]
        p = attach_budget(p, ids, terms, covariance=[[0.]*len(numeric) for _ in numeric],
            covariance_source='Only definitionally exact solar parameters are numerical in this migration',
            assumptions='No independence assumed for unresolved components; no legacy subtotal promoted')
        problems = uc.publication_problems(p)
        if not problems:
            raise RuntimeError('Incomplete historical nitrogen product unexpectedly eligible')
        candidates.append({'historical_index': index, 'product': p,
                           'publication_problems': problems, 'evidence': evidence})
        rows.append({'index': index, 'product_key': p['uncertainty']['scope']['product_key'],
                     'A': p['A'], 'n_indicators': len(ids), 'n_model_lines': p['n_lines'],
                     'state': p['uncertainty']['state'], 'sigma_reported': p['sigma_reported'],
                     'holds': ';'.join(p['uncertainty']['holds'])})
    (args.out/'candidates.json').write_text(json.dumps({'schema': uc.SCHEMA,
        'source_feed': file_evidence(feed), 'candidates': candidates}, indent=2, allow_nan=False)+'\n')
    pd.DataFrame(rows).to_csv(args.out/'summary.csv', index=False)
    print(f'{len(rows)} historical candidates mapped; all HOLD, original feed unchanged')


if __name__ == '__main__':
    main()
