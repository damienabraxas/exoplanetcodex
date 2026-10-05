"""RYA-1233: IR Si I gf from Pehlivan Rhodin et al. 2024 (CDS J/A+A/682/A184)."""
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import rya1233_adopt_pr2024_si as pr  # noqa: E402


def _canon(w):
    rows = [r for r in csv.DictReader(open(ROOT / "data/linelists/canonical_gf.csv"))
            if r["species"] == "Si I" and abs(float(r["wavelength_air_A"]) - w) < 0.01]
    assert len(rows) == 1
    return rows[0]


def test_table_b3_is_air_and_b4_is_vacuum_as_measured():
    """The ReadMe says vacuum for both. B3's 1210.354 nm is Bergemann's AIR 12103.54 A."""
    b3 = pr.read_b3()
    assert abs(b3[b3.lam_nm.round(3) == 1210.354].air_A.iloc[0] - 12103.54) < 0.001
    b4 = pr.read_b4()
    inside = (b4.lam_nm * 10 > 2000) & (b4.lam_nm * 10 < 20000)
    assert (b4[inside].air_A < b4[inside].lam_nm * 10).all()      # vacuum -> air there
    assert (b4[~inside].air_A == b4[~inside].lam_nm * 10).all()   # air outside (ReadMe)


def test_the_calculated_gf_floor_is_measured_from_the_papers_own_table():
    f = pr.calc_floor_dex()
    assert 0.04 < f < 0.05                          # rms(exp - calc) over Table B3 = 0.046


def test_experimental_gf_is_lab_graded_and_calculated_is_priced_not_lab():
    exp = _canon(12103.534)
    assert exp["gf_tier"] == "LAB" and exp["lab_source_tag"] == "PR2024_exp"
    calc = _canon(16434.9265)
    assert calc["gf_tier"] == "PR2024-CALC"
    assert float(calc["gf_sigma_dex"]) >= pr.calc_floor_dex()
    assert abs(float(calc["log_gf"]) - (-1.597)) < 0.001   # was NIST-C+ -1.063: read A=7.06


def test_optical_si_stays_on_asplunds_garz_gf():
    assert _canon(5645.613)["lab_source_tag"] == "AGSS21_Si_Amarsi2017"
