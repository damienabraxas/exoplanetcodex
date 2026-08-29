"""RYA-1109 - the AGSS21 Fe reference line set is what the paper published.

⚠️ WHAT RUNS WHERE. Every test below reads the COMMITTED artifact and runs anywhere,
including CI. The one test that re-extracts from the PDF is skipped when the PDF is absent,
because the source lives outside the repo (`Reference documents/Apslund 2021.pdf`) and is
not redistributable here. That skip is declared rather than silent: a skipped test is not a
passing test, so the artifact-side controls below are written to stand on their own and to
fail if the CSV is edited by hand.
"""
from __future__ import annotations

import csv
import json
import math
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

CSV_PATH = ROOT / "data" / "reference" / "asplund2021_fe" / "asplund2021_fe_lines.csv"
PROV_PATH = ROOT / "data" / "reference" / "asplund2021_fe" / "asplund2021_fe_lines.prov.json"
PDF = Path("/Users/ryanschmitt/Documents/Exoplanet Codex/Reference documents/"
           "Apslund 2021.pdf")

HC_EV_A = 12398.419843320026
PUBLISHED = {"I": 7.46, "II": 7.47}


@pytest.fixture(scope="module")
def rows():
    return list(csv.DictReader(CSV_PATH.open()))


@pytest.fixture(scope="module")
def prov():
    return json.loads(PROV_PATH.read_text())


def _by_ion(rows, ion):
    return [r for r in rows if r["ion"] == ion]


def test_counts_are_the_published_table(rows):
    assert len(rows) == 53
    assert len(_by_ion(rows, "I")) == 40
    assert len(_by_ion(rows, "II")) == 13


def test_fe2_mean_and_the_recommended_average_reproduce_the_paper(rows):
    """The acceptance control. AGSS21 Sect. 4: 7.46 and 7.47 from Fe I and Fe II
    'based on unweighted means', recommended value the average of the two."""
    m = {ion: sum(float(r["a_nlte_3d"]) for r in _by_ion(rows, ion))
         / len(_by_ion(rows, ion)) for ion in ("I", "II")}
    assert m["II"] == pytest.approx(PUBLISHED["II"], abs=0.005)
    assert (m["I"] + m["II"]) / 2 == pytest.approx(7.46, abs=0.005)


def test_the_fe1_discrepancy_is_recorded_not_tuned_away(rows, prov):
    """🔴 Fe I's per-ion mean is 7.4549 against a published 7.46. That is REPORTED.

    It is asserted here so it cannot quietly disappear: if someone later 'fixes' it by
    dropping a line or nudging a value, this fails. The transcription is validated
    independently (energies, row count, the paper's own <3D>/MARCS claims), so the delta is
    a finding about the paper's quoted per-ion value, not a licence to edit the table."""
    fe1 = _by_ion(rows, "I")
    mean = sum(float(r["a_nlte_3d"]) for r in fe1) / len(fe1)
    assert mean == pytest.approx(7.4549, abs=5e-4)
    assert prov["controls"]["Fe I"]["mean_reproduces_published"] is False


def test_energies_are_self_consistent(rows):
    """Eup - Elo == hc/lambda_vac. Validates three columns at once in a quantity the table
    does not tabulate, so a column shift cannot satisfy it."""
    for r in rows:
        lam = float(r["wavelength_air_A"])
        s2 = (1e4 / lam) ** 2
        n = 1 + 8342.13e-8 + 2406030e-8 / (130 - s2) + 15997e-8 / (38.9 - s2)
        dE = HC_EV_A / (lam * n)
        assert abs((float(r["eup_eV"]) - float(r["elo_eV"])) - dE) < 5e-4, r


def test_column_mapping_matches_the_papers_own_claims(rows):
    """AGSS21: with <3D> 'ionisation balance is not fulfilled' while 3D non-LTE shows
    'excellent agreement'; the 1D theoretical (MARCS) mean is 'unrealistically low'."""
    def mean(ion, key):
        s = _by_ion(rows, ion)
        return sum(float(r[key]) for r in s) / len(s)
    gap_3d = abs(mean("I", "a_nlte_3d") - mean("II", "a_nlte_3d"))
    gap_m3d = abs(mean("I", "a_nlte_mean3d") - mean("II", "a_nlte_mean3d"))
    assert gap_m3d > gap_3d
    for ion in ("I", "II"):
        assert mean(ion, "a_nlte_marcs") < mean(ion, "a_nlte_3d")


def test_rew_is_in_the_weak_line_regime(rows):
    """Guards a real bug that was in this extractor: REW computed as mA/A instead of A/A
    lands every line near -2 instead of -5, which would silently break any weak-line cut
    downstream (Amarsi's own is REW < -4.9)."""
    for r in rows:
        rew = float(r["rew"])
        assert -7.0 < rew < -4.0, r
        assert rew == pytest.approx(
            math.log10((float(r["ew_mA"]) / 1000.0) / float(r["wavelength_air_A"])),
            abs=1e-5)


def test_gf_provenance_is_collective_and_per_line_is_empty(rows):
    """RYA-161 firewall. AGSS21 Table A.2 publishes NO per-line gf source, so ours is
    empty and the collective attribution carries the paper's own sentence. Filling it in
    would present OUR cross-match as Asplund's attribution."""
    for r in rows:
        assert r["gf_source_per_line"] == ""
        assert "Den Hartog" in r["gf_source_collective"]
        assert "Belmonte" in r["gf_source_collective"]
        assert r["source_key"] == "asplund2021"
        assert r["line_set"] == "asplund_agss21"


def test_it_is_not_the_ap2002_list(prov):
    """The set we already held is a DIFFERENT paper's. Replicating AGSS21 on AP2002's
    lines reproduces Amarsi's solar row, not Asplund's result."""
    assert "amarsi2022_solar_control_lines.csv" in prov["not_the_same_as"]
    ov = prov["ap2002_overlap"]
    assert ov["Fe I"]["agss21"] == 40 and ov["Fe I"]["ap2002"] == 50
    assert ov["Fe I"]["shared_lambda_ep"] < 40      # they are not the same list


def test_elo_reaches_below_our_graded_floor(rows):
    """The whole reason this list matters: our graded pool floors at Elo 2.85 eV because
    the lab sources measured high-excitation Fe I. AGSS21's set reaches far below it."""
    lo = min(float(r["elo_eV"]) for r in _by_ion(rows, "I"))
    assert lo < 0.1
    assert sum(1 for r in _by_ion(rows, "I") if float(r["elo_eV"]) < 2.85) >= 8


@pytest.mark.skipif(not PDF.exists(),
                    reason="AGSS21 PDF lives outside the repo and is not redistributable; "
                           "the artifact-side controls above cover CI")
def test_the_extractor_reproduces_the_committed_artifact():
    """The control: re-extract from the PDF and require the committed CSV to be what comes
    out. Without it, a hand-edit to the CSV would pass every test above."""
    r = subprocess.run([sys.executable, "scripts/rya1109_extract_asplund_fe_lines.py",
                        "--check-only"], cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "all controls passed" in r.stdout
    before = CSV_PATH.read_text()
    w = subprocess.run([sys.executable, "scripts/rya1109_extract_asplund_fe_lines.py"],
                       cwd=ROOT, capture_output=True, text=True)
    assert w.returncode == 0, w.stdout + w.stderr
    assert CSV_PATH.read_text() == before, (
        "re-extracting the PDF does not reproduce the committed CSV byte-for-byte")
