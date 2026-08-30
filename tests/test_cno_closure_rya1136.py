import csv, json, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "data/audit/rya1136_cno_intake"

def test_closure_build_is_reproducible_and_never_wavelength_only():
    subprocess.run([sys.executable,"scripts/build_cno_intake_rya1136.py"],cwd=ROOT,check=True)
    rows=list(csv.DictReader((AUDIT/"molecular_physical_crossmatch.csv").open()))
    assert len(rows)==408
    assert all(r["identity_basis"] in {"", "wavelength+lower_energy+loggf"} for r in rows)
    assert not any(r["identity_basis"]=="wavelength" for r in rows)
    assert sum(r["join_status"]=="PHYSICAL_TUPLE_MATCH" for r in rows)==80

def test_atomic_and_band_ledgers_are_explicit():
    atomic=list(csv.DictReader((AUDIT/"atomic_source_census.csv").open()))
    assert {r["element"] for r in atomic}=={"C","N","O"}
    coverage=list(csv.DictReader((AUDIT/"combined_coverage_matrix.csv").open()))
    assert len(coverage)==12
    assert {r["band"] for r in coverage}=={"FUV","NUV","VIS","RED_OPTICAL","NIR","IR"}

def test_final_verdict_is_complete_but_not_falsely_ready():
    verdict=json.loads((AUDIT/"intake_verdict.json").read_text())
    assert verdict["intake_census_complete"] is True
    assert verdict["frozen_ready_for_measurement"] is False
    assert verdict["verdict"]=="BLOCKED_MOLECULAR_DATA"
    assert "No abundance derived" in verdict["safety"]
