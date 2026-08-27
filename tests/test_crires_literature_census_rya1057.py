import csv,json
from pathlib import Path
P=Path(__file__).resolve().parents[1]/'data/audit/rya1057_crires_literature'
def read(n):
 with (P/n).open() as f:return list(csv.DictReader(f))
def test_lanes_and_firewall():
 r=read('publication_program_census.csv'); lanes={x['lane'] for x in r}
 assert 'STELLAR ABUNDANCE — PRIMARY' in lanes and 'PLANET/ATMOSPHERE' in lanes and 'LEGACY CRIRES — HISTORICAL' in lanes
 assert all(x['crires_plus_confirmed'].startswith('no') for x in r if x['lane'].startswith('LEGACY'))
def test_elgueta_and_direct_overlaps():
 r=read('publication_program_census.csv'); e=next(x for x in r if x['reference'].startswith('Elgueta'))
 assert all(x in e['codex_overlap'] for x in ('Sun','Procyon','eps Eri'))
 assert 'already vendored' in e['machine_readable_lines']
def test_planet_rows_never_enter_element_index():
 idx=read('per_element_reverse_index.csv'); assert idx and all('planet' not in x['source'].lower() for x in idx)
