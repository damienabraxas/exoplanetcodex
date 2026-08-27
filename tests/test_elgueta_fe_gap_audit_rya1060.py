import csv,json
from pathlib import Path
R=Path(__file__).resolve().parents[1];P=R/'data/audit/rya1060_elgueta_fe_gap'
def rows(n):
 with (P/n).open() as f:return list(csv.DictReader(f))
def test_every_fe_transition_partitioned_once():
 r=rows('elgueta_fe_transition_audit.csv');assert r and {x['species'] for x in r}=={'FeI','FeII'}
 assert all(x['gap_class'] and x['gap_class']!='OTHER' for x in r)
def test_known_elgueta_counts_recomputed():
 r=rows('elgueta_fe_transition_audit.csv')
 assert sum(x['band']=='Y' and x['species']=='FeI' and x['gd_robust']=='Y' for x in r)==0
 assert sum(x['band']=='J' and x['species']=='FeI' and x['gd_robust']=='Y' for x in r)==21
 assert sum(x['band']=='H' and x['species']=='FeI' and x['gd_robust']=='Y' for x in r)==53
 assert sum(x['band']=='Y' and x['species']=='FeI' and x['elgueta_status']=='unassessed' and float(x['wavelength_A'])<10280 for x in r)==71
def test_reconciliation_closes_and_is_read_only():
 assert all(x['unexplained_delta']=='0' for x in rows('band_species_reconciliation.csv'))
 src=(R/'scripts/rya1060_elgueta_fe_gap_audit.py').read_text();assert 'write_text' not in src.split("canonical_gf.csv")[0]
def test_five_recent_branch_lines_are_explicit():
 r=rows('elgueta_fe_transition_audit.csv');w={round(float(x['wavelength_A']),3) for x in r if x['in_recent_rya1054_commit']=='True'}
 assert w=={9913.180,9944.207,10142.844,10216.313,10435.355}
