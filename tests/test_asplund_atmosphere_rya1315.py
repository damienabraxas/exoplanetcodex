"""RYA-1315: Asplund+2021's atmosphere terms come from Scott+2015's per-line 3D/<3D>/HM74 table."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "scripts"), str(ROOT)]

import rya1230_cno_budget as b  # noqa: E402
import rya1230_cno_budget_legs as legs  # noqa: E402


def test_na_pool_terms_are_half_the_scott_differences():
    got = b.scott_atmosphere("Na", "I", [4751.822, 5148.838, 6154.225, 6160.747])
    assert got["n_scott_lines"] == 4 and got["basis"] == "pool lines"
    # Scott+2015a Table 2: <3D> - HM = -0.042,-0.042,-0.041,-0.041; 3D - <3D> = -0.022,-0.022,-0.020,-0.016
    assert abs(got["temperature_structure_dex"] - 0.02075) < 1e-9
    assert abs(got["inhomogeneity_dex"] - 0.01) < 1e-9


def test_pool_without_scott_lines_takes_the_species_value_and_unknown_species_is_none():
    assert b.scott_atmosphere("Ti", "II", [9999.0])["basis"].startswith("species mean")
    assert b.scott_atmosphere("C", "I", [5380.3]) is None


def test_no_rya587_leg_is_run_by_default():
    assert legs.DEFAULT_LEGS == ("nominal",)
    assert legs.FE_EXTRA_LEGS == ("cscale",)
