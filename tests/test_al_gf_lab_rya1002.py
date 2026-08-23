"""tests/test_al_gf_lab_rya1002.py — RYA-1002

Al's primary-laboratory gf table, and the generalisation of `gf_grades` that lets a
non-Fe species reach `GF-LAB` at all.

The two things worth guarding here are not the numbers — they are:
  1. that generalising the module did NOT re-grade Fe (a silent regression on the
     element every other product depends on), and
  2. that the sigma/lambda column trap which defeated RYA-835 cannot come back.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline import gf_grades as gg                      # noqa: E402
from pipeline.wavelength_util import vac_to_air           # noqa: E402

AL_CSV = ROOT / "data" / "reference" / "al_gf_lab" / "al1_lab_loggf.csv"
CM1_TO_EV = 1.0 / 8065.543937


@pytest.fixture(scope="module")
def al() -> pd.DataFrame:
    if not AL_CSV.exists():
        pytest.skip(f"Al lab table absent at {AL_CSV}")
    return pd.read_csv(AL_CSV)


def test_table_is_the_twelve_lines_the_paper_reports(al):
    """Burheim states 12 lines over 670-4200 nm. A different count means the parse
    drifted, and the count is the only end-to-end check on that."""
    assert len(al) == 12
    assert set(al.source) == {"Burheim2023"}
    assert al.lam_vac_A.min() > 6600 and al.lam_vac_A.max() < 42000


def test_every_line_carries_a_measured_sigma(al):
    """The entire reason to prefer this source over a compilation is the per-line bar.
    A blank one would make the table look graded and not be."""
    assert al.e_loggf_dex.notna().all()
    assert (al.e_loggf_dex > 0).all()
    # the paper's stated 2-11% accuracy band, in dex
    assert al.e_loggf_dex.max() <= np.log10(1.11) + 1e-9
    assert al.e_loggf_dex.min() >= np.log10(1.015) - 1e-9


def test_the_column_trap_cannot_return(al):
    """⚠️ THE RYA-835 DEFECT, pinned.

    Burheim prints `sigma [cm^-1]` beside `lambda_vac [A]` and over this range the two
    are the same order of magnitude. RYA-835 read the union of both columns as a
    wavelength list and concluded the paper covered none of our lines, when it covers
    eight. `sigma * lambda == 1e8` is the identity that makes the confusion impossible,
    so it is asserted rather than trusted.
    """
    assert np.allclose(al.sigma_cm1 * al.lam_vac_A, 1e8, rtol=1e-4)
    assert not (set(al.sigma_cm1.round(3)) & set(al.lam_vac_A.round(3)))


def test_air_wavelengths_come_from_the_ssot_converter(al):
    """A second vac<->air formula is the SSOT defect this project forbids."""
    assert np.allclose(al.wavelength_air_A, vac_to_air(al.lam_vac_A.values), atol=1e-6)


def test_levels_are_independent_and_reproduce_the_papers_sigma(al):
    """Levels come from NIST ASD, not from Burheim's own sigma — so `Ek - Ei == sigma`
    is a real cross-source test that the level columns and the wavelength column
    describe the SAME transition, and not an identity."""
    d = al.dropna(subset=["elo_eV", "eup_eV"])
    assert len(d) == 12, "every line should have tied to a NIST level pair"
    implied = (d.eup_eV - d.elo_eV) / CM1_TO_EV
    assert np.allclose(implied, d.sigma_cm1, atol=0.5)


def test_lines_sharing_a_lower_level_share_its_energy(al):
    """Cross-row consistency on the level assignment. 6696/6698 and 13123/13150 are all
    4s 2S1/2 lines; 12749 is 3d 2D5/2 and 12757 is 3d 2D3/2 — the pair must NOT be
    swapped, which a wavelength-only match would not catch."""
    by_label = al.groupby(al.lower_level.str.replace(r"\s+", " ", regex=True)).elo_eV
    for label, g in by_label:
        assert g.nunique() == 1, f"{label} resolves to several energies: {sorted(set(g))}"


# ── the generalisation ────────────────────────────────────────────────────────────────

def test_fe_is_still_the_default_and_is_unchanged():
    """The generalisation is ADDITIVE. Every pre-RYA-1002 caller passes no species and
    must keep the Fe I behaviour exactly."""
    assert gg.DEFAULT_SPECIES == "Fe I"
    assert gg.LAB_CSV == gg.LAB_TABLES["Fe I"]
    fe = gg.lab_lines()
    assert fe is gg.lab_lines("Fe I")
    r = fe.iloc[0]
    v_default = gg.grade_line(float(r.wavelength_air_A), float(r.elo_eV), float(r.loggf))
    v_explicit = gg.grade_line(float(r.wavelength_air_A), float(r.elo_eV), float(r.loggf),
                               species="Fe I")
    assert (v_default.gf_grade, v_default.gf_sigma_dex) == \
           (v_explicit.gf_grade, v_explicit.gf_sigma_dex)
    assert v_default.gf_grade == gg.GRADE_LAB


def test_every_lab_source_is_citable():
    """`grade_line` names the paper in its verdict. A source with no CITATIONS entry
    would raise there, mid-grade, on whichever line happened to hit it first."""
    for species in gg.LAB_TABLES:
        for src in set(gg.lab_lines(species).source):
            assert src in gg.CITATIONS, f"{species}: {src} has no citation"


def test_unregistered_species_raises_rather_than_grading_nothing():
    """POSITIVE CONTROL on the absence. An empty table would send every line to the
    systematic and look identical to 'no lab measurement exists' (RYA-833) — so the
    two must not look alike."""
    with pytest.raises(KeyError, match="no primary-lab gf table registered"):
        gg.lab_lines("Ba II")


def test_al_reaches_GF_LAB_on_the_value_the_linelist_carries(al):
    """The point of the whole ticket: an Al line can now be graded."""
    r = al[al.wavelength_air_A.round(2) == 13123.42].iloc[0]
    v = gg.grade_line(float(r.wavelength_air_A), float(r.elo_eV), float(r.loggf),
                      species="Al I")
    assert v.gf_grade == gg.GRADE_LAB
    assert v.is_graded
    assert v.gf_sigma_dex == pytest.approx(float(r.e_loggf_dex))
    assert "Burheim" in v.gf_grade_source


def test_al_scale_mismatch_fires_on_the_kurucz_vis_values(al):
    """RYA-1001 found `canonical_gf` carries Kurucz-1995 for the VIS pair, 0.113 dex from
    the experiment. A lab value existing is NOT the same as the pool having used it — the
    line must be bounded by the systematic and FLAGGED, never handed Burheim's sigma."""
    r = al[al.wavelength_air_A.round(2) == 6696.02].iloc[0]
    v = gg.grade_line(float(r.wavelength_air_A), float(r.elo_eV), -1.347, species="Al I")
    assert v.gf_grade == gg.GRADE_MISMATCH
    assert not v.is_graded
    assert v.gf_sigma_dex == gg.K07_SYSTEMATIC_DEX
    assert v.gf_delta_dex == pytest.approx(-0.113, abs=1e-3)
