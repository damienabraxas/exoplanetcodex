"""RYA-1047 — canonical_gf extended from 12975.9 A into the H band.

RYA-1046 ingested 28 primary-lab Fe I log gf at 14680-21386 A and they graded NOTHING,
because canonical_gf stopped at 12975.9 A. This is the extension that makes them reachable.

🔴 The source is a file marked `hfsoff_quarantine`. That is admissible for Fe and for Fe
ONLY, and the reason is physics rather than convenience: 56Fe is 91.75% of iron and has
nuclear spin I = 0, so Fe I has no hyperfine structure to switch off. These tests pin both
halves -- that Fe got in, and that the hfs-sensitive species did not.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.gf_grades import canonical_species, grade_line, lab_lines  # noqa: E402

AUDIT = ROOT / "data" / "audit" / "rya1047_hband_gf"
# 🔴 TWO DIFFERENT QUESTIONS, TWO DIFFERENT NUMBERS — AND I CONFLATED THEM.
#   (a) "inside the H arm's measured extent"     -> 14557.3-17500.0 A, 27 of 28
#   (b) "in Elgueta's CERTIFIED H Fe line list"  -> Fe starts 15009.4 A, 25 of 28
# The 2 lines that separate them (14679.844, 14826.412 A) are graded via canonical +
# Ruffoni but have NO Elgueta entry within 1 A, so they are not certified lines.
#
# I first wrote (b)'s count against a bound of my own, 15007-17494 A, which is not in
# the code -- `normalize_vesta_ir.ARMS` sets H's window to None deliberately, because H
# is Elgueta's most telluric-contaminated arm and they skip regions INSIDE it (RYA-787).
# Then I "corrected" it to (a) and silently changed which question the test asked. Both
# tests below now say which one they mean.
H_ARM = (14557.3, 17500.0)          # atomich.dat's measured span


@pytest.fixture(scope="module")
def summary() -> dict:
    return json.loads((AUDIT / "hband_gf_summary.json").read_text())


@pytest.fixture(scope="module")
def fe1() -> pd.DataFrame:
    return canonical_species("Fe I")


def test_canonical_now_covers_the_h_arm(fe1):
    """The whole point: the bound that made 28 lab measurements unreachable is gone."""
    assert fe1.wavelength_air_A.max() > H_ARM[1]
    in_h = fe1[(fe1.wavelength_air_A >= H_ARM[0]) & (fe1.wavelength_air_A <= H_ARM[1])]
    assert len(in_h) > 500, f"only {len(in_h)} Fe I lines in the H arm"


def test_the_lab_lines_actually_grade_now(fe1):
    """A lab table that reaches a band its line list does not is worth nothing. 27 of the
    28 Ruffoni lines must come back as primary-lab graded, carrying their MEASURED sigma
    rather than the 0.17 dex Kurucz blanket (RYA-850)."""
    ruff = lab_lines("Fe I")
    ruff = ruff[ruff.source == "Ruffoni2013"]
    graded = 0
    for _, row in ruff.iterrows():
        m = fe1[(fe1.wavelength_air_A - row.wavelength_air_A).abs() < 0.05]
        if not len(m):
            continue
        v = grade_line(float(m.wavelength_air_A.iloc[0]),
                       float(m.excitation_potential_eV.iloc[0]),
                       float(m.log_gf.iloc[0]), species="Fe I")
        if v.is_graded:
            graded += 1
            assert v.gf_sigma_dex <= 0.13, "a lab line must carry its own measured sigma"
    assert graded == 27


def test_the_ep_check_refused_the_one_line_wavelength_alone_would_have_taken(fe1):
    """🔴 THE SAFETY NET, ON A REAL CASE. Lab line 21386.166 A has a VALD neighbour just
    0.005 A away -- wavelength alone matches it instantly -- but their excitation
    potentials differ by 0.554 eV, so they are DIFFERENT TRANSITIONS. RYA-780 manufactured
    a 2.85 dex discrepancy exactly this way. It must stay unadjudicated."""
    adj = pd.read_csv(AUDIT / "hband_fe1_gf_adjudication.csv")
    near = adj[(adj.wavelength_air_A - 21386.166).abs() < 0.02]
    assert len(near) == 1, "expected the close VALD neighbour to still be there"
    assert not bool(near.adjudicated.iloc[0]), (
        "21386.16 A was adjudicated — the wavelength+EP matcher has been loosened to "
        "wavelength-alone behaviour and is assigning lab gf to the wrong transition")
    assert abs(float(near.excitation_potential_eV.iloc[0]) - 5.6207) > 0.5


def test_matching_beats_a_randomised_null(summary):
    """A match rate is only interpretable against the rate the same matcher achieves on
    randomised wavelengths. 96.4% vs 0.0% is a result; 96.4% on its own is a number."""
    assert summary["match_rate"] > 0.9
    assert summary["randomised_null_rate"] < 0.05
    assert summary["match_rate"] - summary["randomised_null_rate"] > 0.5
    assert summary["ambiguous_left_unadjudicated"] == 0


def test_the_hfs_quarantine_was_checked_not_waved_through(summary):
    """The provenance argument, as data. Fe must be measured untouched by hfs, and every
    species that IS discordant must be one hfs explains -- otherwise the discordance has
    some other cause and the source is not admissible."""
    p = summary["provenance_check"]
    assert p["fe_discordant"] == 0, "the entire hfs argument rests on this"
    assert p["fe_compared"] > 200
    assert set(p["discordant_species"]) <= {"K", "Na", "Al", "Co", "Mn", "V", "Sc", "Cu"}


def test_only_fe_was_admitted(fe1):
    """Fe is admissible from an hfs-off pull; K/Na/Al/Co/Mn are NOT, and must not have
    ridden along. Above the old red edge there should be Fe and nothing else from this
    extension."""
    canon = pd.read_csv(ROOT / "data" / "linelists" / "canonical_gf.csv",
                        comment="#", low_memory=False)
    new = canon[canon.seed_source.astype(str).str.contains("RYA-1047", na=False)]
    assert len(new) > 2000
    assert set(new.species.str.split().str[0]) == {"Fe"}, (
        f"non-Fe species entered from an hfs-off source: "
        f"{sorted(set(new.species.str.split().str[0]) - {'Fe'})}")


def test_adjudicated_rows_carry_their_citation(fe1):
    """RYA-850's measured sigma is only usable if the DOI and sigma travelled with it."""
    canon = pd.read_csv(ROOT / "data" / "linelists" / "canonical_gf.csv",
                        comment="#", low_memory=False)
    new = canon[canon.adjudication_status.astype(str) == "lab_rya1047"]
    assert len(new) == 27
    assert new.gf_sigma_dex.notna().all()
    assert (new.gf_source_doi.astype(str) == "10.1088/0004-637X/779/1/17").all()
    assert new.lab_source_tag.astype(str).eq("Ruffoni2013").all()
