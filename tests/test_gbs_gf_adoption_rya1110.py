"""RYA-1110 - the adopted GBS gf are exactly the 12 Ryan approved, and nothing else moved.

Ryan ratified adopting 12 of the 21 GBS lines Jofre+2014 never published (2026-08-30), on
RYA-1111's held-out validation. These tests pin the ADOPTION ITSELF -- what was taken, what
was refused, whose value it is, and that taking it did not disturb anything Jofre DID
publish.
"""
from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402
from pipeline import line_match, reference_lineset as rls  # noqa: E402


def _src_rows():
    return list(csv.DictReader(GBS.open()))


def _match_into_source(lams, eps):
    """Resolve reference wavelengths into the GBS artifact on the canonical lambda+EP key.

    ⚠️ NOT a rounded-wavelength dict. CI's RYA-1037 guard caught exactly that in the first
    draft of this work, and it was right to: a rounded wavelength is not an identity
    (RYA-1033), so a test built on one can pass while joining the wrong rows.
    """
    src = _src_rows()
    idx = line_match.match(
        np.asarray(lams, float),
        np.array([float(r["wavelength_air_A"]) for r in src]),
        want_ep=np.asarray(eps, float),
        src_ep=np.array([float(r["excitation_potential_eV"]) for r in src]),
        require_ep=True, tol_A=rls.SETS["gbs"].match_tol_A).index
    return src, np.asarray(idx)

GBS = ROOT / "data" / "linelists" / "reference_sets" / "gbs_solar_fe_rya1110.csv"
SIDECAR = ROOT / "data" / "linelists" / "reference_sets" / "gbs_gf_adoption_rya1110.csv"
REPORT = ROOT / "data" / "results" / "rya1111" / "gbs_missing_gf_options.json"


@pytest.fixture(scope="module")
def gbs():
    return rls.load("gbs")


@pytest.fixture(scope="module")
def report():
    return json.loads(REPORT.read_text())


def test_exactly_the_twelve_are_adopted(gbs, report):
    """12 in, 9 refused, 150 measurable. Not 13, not 21."""
    assert int(gbs.gf_adopted.sum()) == 12
    assert int(gbs.gf_missing.sum()) == 9
    assert len(rls.measurable(gbs)) == 150
    assert len(gbs) == 159


def test_the_adopted_set_is_the_report_s_ADOPTABLE_set_exactly(gbs, report):
    """🔴 The selection is the audited one. Not a re-derivation, not a superset: the THIN
    case (BKK, exact but n=2) and all 8 RISKY must stay out."""
    adopted = gbs[gbs.gf_adopted]
    _, idx = _match_into_source(adopted.wavelength_air_A, adopted.elo_eV)
    assert (idx >= 0).all(), "an adopted line does not resolve into the source set"
    want = sorted(round(r["wavelength_air_A"], 2)
                  for r in report["per_line"] if r["verdict"] == "ADOPTABLE")
    got = sorted(round(float(x), 2) for x in adopted.wavelength_air_A)
    assert got == want
    refused = {round(r["wavelength_air_A"], 2)
               for r in report["per_line"] if r["verdict"] in ("THIN", "RISKY")}
    assert set(got) & refused == set(), "a THIN or RISKY line was adopted"


def test_the_thin_case_is_not_adopted(gbs, report):
    """"2 of 2 exact" is not validation (RYA-282's small-n refusal). BKK stays out."""
    thin = [r for r in report["per_line"] if r["verdict"] == "THIN"]
    assert len(thin) == 1
    _, idx = _match_into_source([thin[0]["wavelength_air_A"]], [thin[0]["elo_eV"]])
    assert idx[0] >= 0
    m = line_match.match(
        np.array([thin[0]["wavelength_air_A"]]),
        gbs.wavelength_air_A.to_numpy(float),
        want_ep=np.array([thin[0]["elo_eV"]]),
        src_ep=gbs.elo_eV.to_numpy(float),
        require_ep=True, tol_A=rls.SETS["gbs"].match_tol_A).index
    row = gbs.iloc[int(np.asarray(m)[0])]
    assert not bool(row.gf_adopted)
    assert bool(row.gf_missing)


def test_log_gf_gbs_is_still_empty_for_all_21(gbs):
    """🔴 THE ADOPTED VALUE IS NOT WRITTEN INTO `log_gf_gbs`. That column means "what Jofre
    published", and these 21 lines have no such value. Writing Heiter's number there would
    make the file assert something false about a publication and make the adopted twelve
    indistinguishable from the 138."""
    raw = list(csv.DictReader(GBS.open()))
    blank = [r for r in raw if not r["log_gf_gbs"].strip()]
    assert len(blank) == 21, "the source artifact must still show 21 unpublished gf"


def test_the_rya1110_artifact_was_not_rewritten():
    """The adoption is a SIDECAR. `gbs_solar_fe_rya1110.csv` must be byte-identical to what
    is committed -- the safest edit to a data file is none (RYA-1084)."""
    r = subprocess.run(["git", "status", "--porcelain",
                        "data/linelists/reference_sets/gbs_solar_fe_rya1110.csv"],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.stdout.strip() == "", "the GBS reference artifact was modified"


def test_every_adopted_value_is_heiters_own(gbs):
    """Not a rounded echo: re-read Heiter's column from the source and require a match."""
    adopted = gbs[gbs.gf_adopted]
    src, idx = _match_into_source(adopted.wavelength_air_A, adopted.elo_eV)
    for k, j in enumerate(idx):
        assert j >= 0
        assert abs(float(src[j]["heiter2021_log_gf"])
                   - float(adopted.iloc[k].loggf)) < 5e-4


def test_each_adopted_line_carries_its_provenance_and_says_it_is_adopted(gbs):
    """A reader must be able to tell an adopted gf from a published one, from the row."""
    for r in gbs[gbs.gf_adopted].itertuples():
        assert "ADOPTED" in r.gf_source
        assert "Heiter" in r.gf_source
        assert "RYA-1110" in r.gf_source
        assert r.gf_adopted_source.strip()
    for r in gbs[~gbs.gf_adopted & ~gbs.gf_missing].itertuples():
        assert "ADOPTED" not in r.gf_source, "a published line was labelled adopted"


def test_the_138_published_lines_are_untouched(gbs):
    """RYA-161: the adoption adds values where there were none. It changes nothing Jofre
    published."""
    pub = gbs[~gbs.gf_adopted & ~gbs.gf_missing]
    src, idx = _match_into_source(pub.wavelength_air_A, pub.elo_eV)
    n = 0
    for k, j in enumerate(idx):
        assert j >= 0
        assert src[j]["log_gf_gbs"].strip(), "a line with no published gf is being measured"
        assert abs(float(src[j]["log_gf_gbs"]) - float(pub.iloc[k].loggf)) < 1e-9
        n += 1
    assert n == 138


def test_the_sidecar_matches_the_report():
    """--check: the committed sidecar is what the audited report says to adopt."""
    r = subprocess.run([sys.executable, "scripts/rya1110_adopt_validated_gf.py", "--check"],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "MATCHES" in r.stdout


def test_the_sidecar_never_overrides_a_published_gf():
    """Structural: every row in the sidecar names a line whose `log_gf_gbs` is empty."""
    side = list(csv.DictReader(SIDECAR.open()))
    src, idx = _match_into_source([r["wavelength_air_A"] for r in side],
                                  [r["elo_eV"] for r in side])
    for j in idx:
        assert j >= 0
        assert not src[j]["log_gf_gbs"].strip()
