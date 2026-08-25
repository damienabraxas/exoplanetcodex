"""RYA-1046 — Ruffoni et al. 2013 H-band Fe I laboratory log gf.

WHY THIS TABLE EXISTS. Our Fe I lab pool stopped at 11316.06 A (RYA-834 checked that and
stated it). CRIRES+ is the science instrument and its J/H arms sit at 14000-17500 A, so
every CRIRES Fe I line was rung 1 -- the 0.17 dex Kurucz blanket -- for want of a lab
measurement, not because none exists. This is the only laboratory Fe I gf we hold above
11316 A.

🔴 THE COLUMN CHOICE IS THE WHOLE TICKET. Ruffoni's Table 6 carries THREE log gf columns:
"BF & Effective tau", "Ladenburg", and "Recommended". The first and third refine the
radiative lifetimes by fitting the SOLAR spectrum, and their quoted uncertainty folds in
12% from Asplund's solar Fe abundance. Adopting those to measure a solar Fe abundance
would be circular -- RYA-161's firewall, the same test that keeps `melendez1999` out.
We take LADENBURG, the lab emission/absorption network, and these tests pin that.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.gf_grades import CITATIONS, lab_lines  # noqa: E402

TABLE = ROOT / "data" / "reference" / "fe_gf_lab" / "fe1_lab_loggf_ruffoni2013.csv"


@pytest.fixture(scope="module")
def ruff() -> pd.DataFrame:
    return pd.read_csv(TABLE)


def test_the_table_is_registered_and_cited(ruff):
    """A lab table whose source is missing from CITATIONS grades lines uncited."""
    assert TABLE.exists()
    assert set(ruff.source.unique()) == {"Ruffoni2013"}
    cite, doi = CITATIONS["Ruffoni2013"]
    assert "Ruffoni" in cite and "2013" in cite
    assert doi == "10.1088/0004-637X/779/1/17"
    # The citation must SAY which column, because the other two are disqualified.
    assert "Ladenburg" in cite


def test_it_reaches_where_nothing_else_does(ruff):
    """The point of the ingest: lab reach above 11316.06 A, where there was none."""
    lab = lab_lines("Fe I")
    older = lab[lab.source != "Ruffoni2013"]
    assert older.wavelength_air_A.max() == pytest.approx(11316.06, abs=0.1)
    assert ruff.wavelength_air_A.min() > older.wavelength_air_A.max()
    # ... and it covers the CRIRES+ H arm (Elgueta 15007-17494 A).
    in_h = ruff[(ruff.wavelength_air_A >= 15007.0) & (ruff.wavelength_air_A <= 17494.0)]
    assert len(in_h) == 25


def test_the_two_fe1_tables_never_claim_the_same_line(ruff):
    """Two laboratory values for one line is an adjudication (RYA-780), not a tie-break.
    `lab_lines` refuses on overlap; this asserts the tables actually are disjoint so that
    refusal is never reached in normal operation."""
    lab = lab_lines("Fe I")
    assert not lab.duplicated("wavelength_air_A").any()


def test_every_line_carries_a_measured_sigma(ruff):
    """RYA-850: a MEASURED per-line sigma supersedes the consistency bound. A fallback
    sigma dressed as a measurement is the failure that rule exists to stop."""
    assert (ruff.e_loggf_dex > 0).all()
    assert ruff.e_loggf_dex.max() <= 0.5
    # Ruffoni quotes 0.05-0.13 dex; a constant column would mean we invented it.
    assert ruff.e_loggf_dex.nunique() > 1


def test_wavelengths_are_air_not_vacuum(ruff):
    """RYA-944: a ticket said VACUUM and the atlas was AIR. In the H band the two differ
    by ~4.4 A -- far more than any matching tolerance -- so getting this wrong would
    match every line to the wrong transition or to nothing."""
    # Air wavelength recomputed from the level energies, via vacuum + the air index.
    # RYA-264/501: the refractive-index formula lives in exactly ONE place. Re-deriving
    # it here would be a second copy to drift -- and this test would then be validating
    # the ingest against a formula the pipeline does not use.
    from pipeline.wavelength_util import vac_to_air
    sigma = ruff.eup_cm1 - ruff.elo_cm1              # cm^-1
    lam_vac = 1.0e8 / sigma                          # Angstrom, vacuum
    lam_air = vac_to_air(lam_vac.to_numpy())
    assert (abs(lam_air - ruff.wavelength_air_A) < 0.05).all()
    # And the control: had we recorded vacuum, this would be the ~4 A error we avoided.
    assert (lam_vac - ruff.wavelength_air_A).min() > 3.0
