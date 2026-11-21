"""Build an exposed O I source-grid matrix; no adoption or laboratory claims.

Wavelength proximity yields identity candidates, never canonical physical joins.
Retain every candidate (including distinct M1/E2 operators) and every missing
prediction. Only explicitly checked native quintet identities receive residuals.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import re

ROOT = Path(__file__).resolve().parents[1]


def rows(path, delimiter=','):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def numeric(value):
    # ASD '+' marks a perturbed level; preserve the original string alongside.
    return float(re.sub(r'[^0-9.]', '', value))


def triplet_reference(matches):
    """Refuse wavelength-only matches to different configurations/terms/J."""
    checked = [r for r in matches if r['Lower level'].split('|')[0].strip() == '2s2.2p3.(4S*).3s'
               and r['Lower level'].split('|')[1].strip() == '5S*'
               and r['Lower level'].split('|')[2].strip() == '2'
               and r['Upper level'].split('|')[0].strip() == '2s2.2p3.(4S*).3p'
               and r['Upper level'].split('|')[1].strip() == '5P'
               and r['Upper level'].split('|')[2].strip() in ('1','2','3')]
    if len(checked) != 1:
        raise ValueError('Triplet identity ambiguous or absent')
    return checked[0]


def build(root=ROOT, destination=None):
    out = root/'data/audit/rya1354/nonblind_v1'
    census_path = root/'data/audit/rya1136_cno_intake/atomic_source_census.csv'
    canonical_path = root/'data/linelists/canonical_gf.csv'
    nist_path = out/'nist_OI_3000_10000.csv'
    candidate_path = root/'data/audit/rya1354/exploratory_candidates_v3.json'
    product_path = root/'data/products/solar/O.json'
    native_path = out/'native_summary.json'
    inputs = [census_path, canonical_path, nist_path, candidate_path,
              product_path, root/'data/config/elements_master.json', native_path]
    canonical = {r['line_id']: r for r in rows(canonical_path)}
    nist = rows(nist_path)
    predictions = json.loads(candidate_path.read_text())['candidates']
    native = json.loads(native_path.read_text())
    product_scopes = []
    for p in json.loads(product_path.read_text())['products']:
        scope = {k:p.get(k) for k in ('grade','selector','tier','band','instrument','holding','route','treatment')}
        scope['recorded_indicator_ids'] = sorted({identity
            for c in p.get('uncertainty',{}).get('components',[])
            for identity in c.get('evidence',{}).get('indicator_ids',[])})
        product_scopes.append(scope)
    origins = {'T5803': dict(doi='10.1088/0953-4075/24/18/010',origin='THEORY_CIV3_WITH_EMPIRICAL_ENERGY_ADJUSTMENTS',
                            inspected='NIST citation and publisher abstract; full line tables not inspected'),
               'T4539': dict(doi='10.1103/PhysRevA.28.3169',origin='THEORY_MCHF_BREIT_PAULI',
                            inspected='NIST citation/title only'),
               'T5081': dict(doi='10.1088/0953-4075/21/9/007',origin='UNVERIFIED_PRIMARY_METHOD',
                            inspected='NIST citation/title only; not accepted as laboratory measurement')}
    for code in origins:
        inputs.append(out/f'nist_{code}.html')
    matrix = []
    for grid in rows(census_path):
        if grid['species'] != 'O I':
            continue
        wave = float(grid['wavelength_vac_A'])
        # Source grid uses rounded component wavelengths. All nearby records are
        # retained. Dense multiplets remain ambiguous pending source-level review.
        matches = [r for r in nist if r['Ritz'] and abs(numeric(r['Ritz'])-wave) <= .03]
        canon = canonical.get(grid['canonical_line_id'], {})
        entry = dict(source_grid=grid, actual_canonical_tier=canon.get('gf_tier'),
                     canonical_physical_id=canon.get('physical_id'),
                     product_pool_membership_via_existing_canonical_tuple=[s for s in product_scopes
                         if canon.get('physical_id') and canon['physical_id'] in s['recorded_indicator_ids']],
                     actual_canonical_nist_grade=canon.get('nist_grade'),
                     actual_canonical_gf_reference=canon.get('loggf_reference'),
                     reference_use='SOURCE_ANALYSIS_GRID_NOT_PROOF_OF_FINAL_ADOPTION',
                     grade_scope='REFERENCE_PRODUCT_DESIGNATION_IS_NOT_LAB_GF',
                     identity_status='PROVISIONAL_PROXIMITY_CANDIDATES_NOT_CANONICAL_JOIN',
                     nist_candidates=matches,
                     nist_gf_provenance={code:info for code,info in origins.items()
                                         if any(code in m['TP'] for m in matches)},
                     gf_origin='COMPILED_ANALYSIS_GRID_UPSTREAM_PRIMARY_PROVENANCE_UNRESOLVED',
                     uncertainty=None, covariance=None,
                     dependence_group=grid['line_label'],
                     dependence_note='Shared multiplet/upper-level and common source scale must be curated before covariance scoring',
                     isotope_hfs='UNSPECIFIED_NOT_RESOLVED_BY_WAVELENGTH',
                     exposure='EXPOSED_NONBLIND', prediction_status='NOT_CALCULATED',
                     experimental_gf_accuracy_score=None)
        if 7770 < float(grid['wavelength_air_A']) < 7780:
            # Independent checks of configuration, parent term, term and J;
            # no closest-wavelength assignment of calculated energies.
            reference = triplet_reference(matches)
            upper_j = int(reference['Upper level'].split('|')[2])
            prediction, = [p for p in predictions if p['two_J_upper'] == 2*upper_j]
            ei, ek = map(float, reference['Ei           Ek'].split('  -'))
            entry.update(identity_status='CONFIGURATION_TERM_J_PARITY_CHECKED_ISOTOPE_UNSPECIFIED',
                         prediction_status='SINGLE_CONFIGURATION_EXPLORATORY',
                         prediction_id=prediction['candidate_id'],
                         energy_comparison=dict(reference_level_difference_cm1=ek-ei,
                                                native_prediction_cm1=prediction['transition_energy_cm1_native'],
                                                residual_cm1=prediction['transition_energy_cm1_native']-(ek-ei),
                                                uncertainty=None, purpose='LEVEL_SEPARATION_NOT_EXPERIMENTAL_GF_TEST'),
                         exposed_grid_loggf_residual_dex=math.log10(prediction['gf_B'])-float(grid['published_loggf']),
                         comparison_kind='THEORY_VERSUS_COMPILED_GRID_NOT_LAB_ACCURACY')
        elif len(matches)==1:
            reference=matches[0]
            lower_conf, lower_term, lower_j=map(str.strip,reference['Lower level'].split('|'))
            upper_conf, upper_term, upper_j=map(str.strip,reference['Upper level'].split('|'))
            family=None
            # Full source-qualified parent/configuration/term checks precede J selection.
            if lower_conf=='2s2.2p3.(4S*).3p' and lower_term=='5P' and upper_term=='5D*':
                family={'2s2.2p3.(4S*).3d':'IR_3d_E1',
                        '2s2.2p3.(4S*).4d':'optical_4d_E1'}.get(upper_conf)
            if (lower_conf,lower_term,upper_conf,upper_term)==('2s2.2p3.(4S*).3s','3S*','2s2.2p3.(4S*).3p','3P'):
                family='triplet_8446_E1'
            if family:
                prediction, = [p for p in native['coverage'][family]
                                if p['lower_J']==int(lower_j) and p['upper_J']==int(upper_j)]
                ei,ek=map(float,reference['Ei           Ek'].split('  -'))
                entry.update(identity_status='CONFIGURATION_TERM_J_PARITY_CHECKED_ISOTOPE_UNSPECIFIED',
                             prediction_status='REFERENCE_CONFIGURATION_EXPLORATORY',native_prediction=prediction,
                             native_table=f'runs/coverage_v8/{family}/lower.upper.t',
                             energy_comparison=dict(reference_level_difference_cm1=ek-ei,
                                                    native_prediction_cm1=prediction['native_energy_cm1'],
                                                    residual_cm1=prediction['native_energy_cm1']-(ek-ei),
                                                    uncertainty=None,purpose='LEVEL_SEPARATION_NOT_EXPERIMENTAL_GF_TEST'),
                             exposed_grid_loggf_residual_dex=math.log10(prediction['gf_B'])-float(grid['published_loggf']),
                             comparison_kind='THEORY_VERSUS_COMPILED_GRID_NOT_LAB_ACCURACY')
        elif matches and {m['Type'] for m in matches}=={'M1','E2'}:
            entry.update(prediction_status='NATIVE_OPERATORS_AVAILABLE_REVIEW_NOT_SCORED',
                         native_table=native['forbidden']['table'],
                         operator_warning=native['forbidden']['warning'])
        if entry['prediction_status']=='NOT_CALCULATED' and matches:
            entry['prediction_status']='HOLD_AMBIGUOUS_GRID_COMPONENT_IDENTITY'
            entry['related_native_table']='runs/coverage_v8/optical_4d_E1/lower.upper.t'
        matrix.append(entry)
    if len(matrix) != 26:
        raise ValueError('O reference grid count changed; review scope explicitly')
    targets = json.loads((root/'data/config/elements_master.json').read_text())['elements']
    all_canon = list(canonical.values())
    prospect = []
    for target in targets:
        symbol = target['symbol']
        species = 'Fe II' if symbol == 'Fe II' else symbol+' I'
        held = [r for r in all_canon if r['species'] == species]
        prospect.append(dict(target=symbol, priority=target['priority'], canonical_species=species,
                             canonical_lines=len(held), held_lab_tier_rows=sum(r['gf_tier']=='LAB' for r in held),
                             exposure='PROJECT_HOLDINGS_EXPOSED_NOT_BLIND',
                             prospective_status='REQUIRES_DISTINCT_CURATOR_ROLE_AND_GENUINELY_WITHHELD_DATA',
                             partition='upper-level + multiplet + lifetime/branching normalization + common calibration scale'))
    result = dict(schema='rya1354.nonblind-reference-matrix.v1', scientific_admission='HOLD',
                  inputs={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
                  scope='ALL_26_SOURCE_GRID_COMPONENTS_NOT_ALL_PROVEN_ADOPTED_REFERENCE_LINES',
                  nist_wavelength_frame='VACUUM_ANGSTROM_SHOW_AV_3',
                  matching_tolerance_vac_A=.03, matching_status='DISCOVERY_ONLY_EXCEPT_EXPLICIT_TRIPLET_LABEL_CHECK',
                  rows=matrix, canonical_target_screen=prospect,
                  actual_product_scopes=product_scopes,
                  product_scope_warning='Actual Reference Grade products use SET-AGSS21; this does not establish which of all 26 grid components enter each fitted pool. No grade propagated from wavelength proximity.',
                  prospective_shortlist=['Na','Mg'],
                  shortlist_basis='Canonical priority-1 light atoms with simpler valence structures than oxygen; selection is methodological, not a blindness claim. Curator must verify new independent measurements and coverage before registration.',
                  scoring_summary=dict(grid_rows=26, native_prediction_rows=sum(r['prediction_status'] not in ('NOT_CALCULATED','HOLD_AMBIGUOUS_GRID_COMPONENT_IDENTITY') for r in matrix),
                                       scored_nonblind_compiled_gf_rows=sum('exposed_grid_loggf_residual_dex' in r for r in matrix),
                                       unresolved_grid_identity_rows=sum(r['prediction_status']=='HOLD_AMBIGUOUS_GRID_COMPONENT_IDENTITY' for r in matrix),
                                       experimental_gf_accuracy_rows=0, calibrated_chi_square=None,
                                       missing_predictions_retained=True))
    (destination or out/'reference_matrix.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result['scoring_summary']))


if __name__ == '__main__':
    build()
