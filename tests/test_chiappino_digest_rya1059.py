from __future__ import annotations
import json, subprocess, sys
from pathlib import Path
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'data/audit/rya1059_chiappino'; RAW=ROOT/'data/reference/chiappino2026'

def test_reproducible(tmp_path):
 subprocess.run([sys.executable,str(ROOT/'scripts/rya1059_digest_chiappino.py'),'--out',str(tmp_path)],check=True)
 for p in OUT.iterdir(): assert (tmp_path/p.name).read_bytes()==p.read_bytes()

def test_source_package_is_checksum_guarded():
 m=pd.read_csv(RAW/'SHA256SUMS.csv'); assert len(m)==4; assert (m.sha256.str.len()==64).all()

def test_atomic_and_molecular_systems_stay_separate():
 d=pd.read_csv(OUT/'normalized_lines.csv'); assert {'atomic','molecular'}==set(d.atomic_or_molecular)
 assert {'C12O','C13O','CN','OH'} <= set(d.species)
 assert not d[d.atomic_or_molecular=='molecular'].canonical_line_id.notna().any()

def test_no_empirical_line_is_promoted():
 d=pd.read_csv(OUT/'normalized_lines.csv'); assert not d.promotion_allowed.astype(bool).any()
 assert d[d.atomic_or_molecular=='atomic'].traced_primary_source.eq('UNTRACED_BEYOND_VALD3').all()

def test_match_is_physical_not_wavelength_only():
 d=pd.read_csv(OUT/'normalized_lines.csv'); assert {'PHYSICAL_KEY_UNIQUE','NO_PHYSICAL_KEY_MATCH','MOLECULAR_SCHEMA_SEPARATE'} <= set(d.canonical_match_status)
 assert d[d.canonical_match_status=='PHYSICAL_KEY_UNIQUE'].canonical_line_id.notna().all()

def test_summary_and_required_routing():
 s=json.loads((OUT/'summary.json').read_text()); assert s['mutation_policy']=='DISCOVERY_REFERENCE_ONLY'; assert s['lines']>200
 assert (OUT/'fe_jhk_delta.csv').exists(); assert (OUT/'al_completeness_delta.csv').exists(); assert (OUT/'cno_diagnostic_source_map.csv').exists(); assert (OUT/'empirical_lines_failing_provenance.csv').exists()
