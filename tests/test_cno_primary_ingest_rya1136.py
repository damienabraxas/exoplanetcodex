"""RYA-1136 primary molecular-source normalization checks."""
import importlib.util
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "primary_cno", ROOT / "scripts/ingest_cno_molecular_primary_rya1136.py")
MOD = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MOD
SPEC.loader.exec_module(MOD)


def nearest(rows, species, wavelength_A):
    return min((row for row in rows if row.species == species),
               key=lambda row: abs(row.wavelength_vac_A - wavelength_A))


def test_ch_source_gf_is_already_weighted():
    # Amarsi Table 2 row 40: 425.420 nm vacuum, Elow=0.523 eV, loggf=-1.506.
    row = nearest(MOD.parse_ch(), "CH", 4254.20)
    assert abs(row.wavelength_vac_A - 4254.20) < 0.02
    assert abs(row.lower_energy_eV - 0.523) < 0.002
    assert abs(row.loggf - (-1.506)) < 0.002


def test_cn_source_f_requires_lower_statistical_weight():
    # Amarsi Table 2 CN row at 1087.516 nm; Brooke tabulates f, not gf.
    row = nearest(MOD.parse_cn(), "CN", 10875.16)
    assert row.system == "A-X"
    assert abs(row.lower_energy_eV - 0.056) < 0.002
    assert abs(row.loggf - (-2.601)) < 0.002
    assert math.isclose(row.gf, (2 * row.j_lower + 1) * 7.827300e-5,
                        rel_tol=1e-8)


def test_primary_crossmatch_is_not_wavelength_only():
    text = (ROOT / "scripts/ingest_cno_molecular_primary_rya1136.py").read_text()
    assert "abs(x.lower_energy_eV - energy)" in text
    assert "abs(x.loggf - loggf)" in text
