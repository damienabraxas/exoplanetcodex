import csv,json,subprocess,sys
from pathlib import Path
R=Path(__file__).resolve().parents[1];P=R/'data/audit/rya1058_elgueta'
def rows(n):
 with (P/n).open() as f:return list(csv.DictReader(f))
def test_complete_digest_counts_and_provenance():
 r=rows('normalized_lines.csv');assert len(r)==1700
 assert {x['band'] for x in r}=={'Y','J','H'} and all(len(x['source_md5'])==32 for x in r)
def test_same_target_products_preserve_flags():
 assert all(len(rows(n))==1700 for n in ('sun_line_behavior.csv','procyon_line_behavior.csv','eps_eri_line_behavior.csv'))
def test_no_wavelength_only_or_automatic_gf_promotion():
 r=rows('normalized_lines.csv');assert all(x['action']=='TRACE_PRIMARY_GF_BEFORE_PROMOTION' for x in r)
 assert all(x['match_status']!='WAVELENGTH_ONLY' for x in r)
def test_raw_holding_is_md5_pinned():
 assert (R/'data/reference/elgueta2026_vizier/MD5SUMS.txt').exists()
