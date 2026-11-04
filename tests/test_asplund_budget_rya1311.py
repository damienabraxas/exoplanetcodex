"""RYA-1311 (Ryan 2026-11-03): the published uncertainty is Asplund+2021's, read off the cell's
own nominal products -- no RYA-587 perturbation leg is run."""
import importlib.util
import inspect
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))


def _budget():
    spec = importlib.util.spec_from_file_location("rya1230_cno_budget",
                                                  ROOT / "scripts" / "rya1230_cno_budget.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _lines(path, deltas, accepted=True):
    pd.DataFrame({"wavelength_air_A": list(deltas), "nlte_delta_dex": list(deltas.values()),
                  "in_aggregate": accepted}).to_csv(path, index=False)


def test_the_3d_effect_comes_from_the_units_engine_a_pair_for_any_treatment(tmp_path):
    b = _budget()
    base = "CaI_4200_6908_harps_solar_harps_molecfit_corrected_SYNTH_SET-CA_AGSS21"
    _lines(tmp_path / f"{base}_ENGINE-A-3DNLTE_lines.csv", {6162.17: 0.10, 6439.07: 0.06})
    _lines(tmp_path / f"{base}_ENGINE-A_lines.csv", {6162.17: 0.02, 6439.07: 0.00})
    for treatment in ("1D-LTE", "ENGINE-B-NLTE", "synth-mean3D-NLTE-gerber-stagger", "ENGINE-A-3DNLTE"):
        got = b._unit_3d_effect(tmp_path, f"{base}_{treatment}_lines.csv", np.mean)
        assert np.isclose(got, 0.07), treatment


def test_no_engine_a_pair_means_no_3d_term(tmp_path):
    assert _budget()._unit_3d_effect(tmp_path, "CaII_1_2_x_y_SYNTH_REFERENCE_1D-LTE_lines.csv",
                                     np.mean) is None


def test_publish_defaults_to_the_asplund_budget_and_publishes_every_grade_but_new_deep():
    from pipeline import element_publish as ep
    assert inspect.signature(ep.publish).parameters["budget"].default == "asplund"
    src = inspect.getsource(ep.units_for)
    assert 'pools=["set", "reference", "codex", "deep"]' in src
