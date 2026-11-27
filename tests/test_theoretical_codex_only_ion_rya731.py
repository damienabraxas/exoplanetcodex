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


def test_clear_set_applies_source_gf_to_the_matching_transition(tmp_path):
    import numpy as np
    from scripts.derive_band_products import _cand_from_ew_artifact

    source = tmp_path / "clear.csv"
    source.write_text("wavelength_air_A,ep_eV,match_tol_A,loggf_published\n"
                      "3200,2.9,0.005,-2.5\n")
    linelist = np.array(
        [("Ni 2", 3200.0, 4.0, -1.0, 0.2), ("Ni 2", 3200.002, 2.9, -2.0, 0.1),
         ("Fe 1", 3200.002, 2.9, -3.0, 0.3)],
        dtype=[("element", "U4"), ("wave_A", "f8"), ("lower_state_eV", "f8"),
               ("loggf", "f8"), ("theoretical_depth", "f8")],
    )

    candidates = _cand_from_ew_artifact(
        linelist, source, lo_A=3000, hi_A=3780, species="Ni 2", tier="all",
        source_label="lines-from-set:NI_CLEAR2026",
    )

    assert candidates.loggf.tolist() == [-2.5]
    assert linelist["loggf"].tolist() == [-1.0, -2.5, -3.0]
