"""Mutation tests for scientific scope/provenance errors, not just JSON shape."""
import copy,json,unittest
from jsonschema.exceptions import ValidationError
from pathlib import Path
from pipeline.frontier_atlas import ROOT,build,validate,safe_path,derive
class FrontierAtlasTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.data=build()
    def mutate(self,fn):
        data=copy.deepcopy(self.data);fn(data)
        with self.assertRaises((ValueError,ValidationError)):validate(data)
    def test_generated_snapshot_matches_fresh_source_inventory(self):
        self.assertEqual(self.data,json.loads((ROOT/'data/products/frontier/coverage.v1.json').read_text()))
    def test_feii_2019_lte_cannot_be_relabelled_nlte(self):
        self.mutate(lambda d:next(e for e in d['entries'] if e['id']=='fe-ii-amarsi2019-table4').update(treatment='FULL3D_NLTE'))
    def test_zero_feii_transitions_cannot_be_promoted(self):
        self.mutate(lambda d:next(e for e in d['entries'] if e['id']=='fe-ii-gerber-atom').update(coverage_status='unvalidated'))
    def test_no_scientific_adoption_from_inventory(self):
        self.mutate(lambda d:d['entries'][0].update(coverage_status='validated'))
        self.mutate(lambda d:d['entries'][0].update(execution_approved=True))
    def test_no_claim_of_literature_absence_from_catalog_gap(self):
        self.mutate(lambda d:d['entries'][0].update(coverage_status='not_published'))
    def test_model_atom_is_not_an_atmosphere(self):
        self.mutate(lambda d:d['entries'][0]['components']['atmosphere_ids'].append('fe-i-gerber-atom'))
    def test_molecular_list_does_not_imply_nlte(self):
        self.mutate(lambda d:next(e for e in d['entries'] if e['component_type']=='line_list').update(treatment='FULL3D_NLTE'))
    def test_hash_drift_and_dangling_sources_are_rejected(self):
        self.mutate(lambda d:d['sources']['cno-lines'].update(sha256='0'*64))
        self.mutate(lambda d:d['entries'][0]['source_ids'].append('invented-source'))
    def test_future_audit_and_impossible_domain_rejected(self):
        self.mutate(lambda d:d['entries'][0]['literature_audit'].update(date='2099-01-01'))
        self.mutate(lambda d:d['entries'][0]['domain']['axes'][0].update(min=1e9))
    def test_no_repository_escape_or_unsafe_evidence(self):
        with self.assertRaises(ValueError):safe_path(ROOT,'../../outside.json')
        self.mutate(lambda d:d['sources']['cno-lines'].update(evidence_url='javascript:alert(1)'))
    def test_unknown_species_does_not_enter_atomic_inventory(self):
        self.mutate(lambda d:d['entries'][0]['species'].update(element='Invented'))
    def test_species_and_scope_remain_distinct(self):
        entries={e['id']:e for e in self.data['entries']}
        self.assertFalse(entries['fe-ii-gerber-atom']['coverage']['atom_reach']['nlte_capable'])
        for key in ['o-i-amarsi2019-3n','o-i-amarsi2019-1n']:
            e=entries[key];self.assertTrue(any(label.startswith('777') for label in e['coverage']['source_labels']));self.assertTrue(e['domain']['axes'])
        for molecule in ['h2o','tio','feh','cah','mgh','sio']:
            e=entries[molecule+'-discovery'];self.assertEqual(e['coverage_status'],'not_yet_searched');self.assertIsNone(e['literature_audit']['date'])
if __name__=='__main__':unittest.main()
