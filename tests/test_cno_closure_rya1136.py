import csv, json, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "data/audit/rya1136_cno_intake"

def test_closure_build_is_reproducible_and_never_wavelength_only():
    subprocess.run([sys.executable,"scripts/build_cno_intake_rya1136.py"],cwd=ROOT,check=True)
    rows=list(csv.DictReader((AUDIT/"molecular_physical_crossmatch.csv").open()))
    assert len(rows)==408
    assert all(r["identity_basis"] in {"", "wavelength+lower_energy+loggf",
                                             "wavenumber+band+lower_energy+gf"} for r in rows)
    assert not any(r["identity_basis"]=="wavelength" for r in rows)
    assert sum(r["join_status"]=="PHYSICAL_TUPLE_MATCH" for r in rows)==80
    # 🔴 364, NOT 390. RYA-1144 found that 26 of 32 sum-matches were argmin-resolved and
    # counted as coverage; >1 viable subset now yields AMBIGUOUS_SUM_MATCH, which is
    # deliberately NOT an accepted join. Its own words: "accepted coverage falls 390 ->
    # 364. That drop is the point: those 26 features were fitted, not identified."
    # The fix landed in e4b52d56 and this assertion never followed it, so the suite has
    # been red on main ever since -- asserting a number the builder deliberately stopped
    # producing. Verified: the committed artifact and a fresh rebuild both give 364.
    assert sum(r["join_status"] in {"PHYSICAL_TUPLE_MATCH", "PRIMARY_TUPLE_MATCH",
                                    "PRIMARY_UNRESOLVED_SUM_MATCH"} for r in rows)==364

def test_atomic_and_band_ledgers_are_explicit():
    atomic=list(csv.DictReader((AUDIT/"atomic_source_census.csv").open()))
    # 🔴 RYA-1183. The SOURCE elements are still exactly C/N/O. Ni appears only as the
    # [O I] 6300.30 blend component — Ni I 6300.34 is the contaminant that makes the
    # diagnostic unusable alone, and a census listing the [O I] line without it describes
    # a feature that does not exist in isolation. It is marked
    # BLEND_COMPONENT_NOT_A_SOURCE_LINE so it can never be counted as an AGSS21 line.
    source_rows = [r for r in atomic if r["use_status"] != "BLEND_COMPONENT_NOT_A_SOURCE_LINE"]
    assert {r["element"] for r in source_rows} == {"C", "N", "O"}
    blend = [r for r in atomic if r["use_status"] == "BLEND_COMPONENT_NOT_A_SOURCE_LINE"]
    assert {r["element"] for r in blend} == {"Ni"}
    assert len(blend) == 1 and blend[0]["blend_role"] == "CONTAMINANT_COMPONENT"
    coverage=list(csv.DictReader((AUDIT/"combined_coverage_matrix.csv").open()))
    assert len(coverage)==12
    assert {r["band"] for r in coverage}=={"FUV","NUV","VIS","RED_OPTICAL","NIR","IR"}
    # Same RYA-1144 correction, read through the coverage matrix. RYA-1150 additionally
    # made the builder own summary.json and assert loudly if the two ever disagree on the
    # molecular total again -- "Both read 364 accepted / 44 review."
    assert sum(int(r["matched"]) for r in coverage if r["domain"]=="molecular")==364

def test_final_verdict_is_complete_but_not_falsely_ready():
    verdict=json.loads((AUDIT/"intake_verdict.json").read_text())
    assert verdict["intake_census_complete"] is True
    assert verdict["frozen_ready_for_measurement"] is False
    assert verdict["verdict"]=="INTAKE_COMPLETE_REVIEW_REQUIRED"
    assert "No abundance derived" in verdict["safety"]
