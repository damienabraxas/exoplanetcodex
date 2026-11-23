"""RYA-731 -- an ion whose only lines are an adopted Theoretical Codex set (Ni II, Clear+2026
gf) gets set cells and NO lab-pool cell; ions with graded lines are unchanged."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def test_tc_only_ion_gets_set_cells_only():
    from pipeline import run_matrix as rm
    assert "II" not in rm.graded_ions("Ni")
    assert "II" in rm.theoretical_codex_ions("Ni")
    ii = [d for d in rm.expand("solar", "Ni") if d.ion == "II"]
    assert ii, "Ni II has an adopted Theoretical Codex set and must get cells"
    assert {d.pool for d in ii} == {rm.SET_POOL_PREFIX + "NI_CLEAR2026"}


def test_graded_ion_keeps_its_lab_pools():
    from pipeline import run_matrix as rm
    pools = {d.pool for d in rm.expand("solar", "Ni") if d.ion == "I"}
    assert pools & set(rm.POOLS), "Ni I must keep its lab-pool cells"
