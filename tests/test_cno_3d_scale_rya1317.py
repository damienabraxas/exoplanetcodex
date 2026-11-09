"""RYA-1317: CNO terms come from the indicator's own Amarsi group, on Asplund's 3D scale."""
import json
from pathlib import Path

import pytest

from pipeline import molecular_model_form as mm
from pipeline.uncertainty_contract import asplund_sigma

ROOT = Path(__file__).resolve().parents[1]


def doc(**terms):
    return {"components": [{"name": n, "state": s, "sigma_dex": v} for n, (s, v) in terms.items()]}


def test_g_band_uses_ch_a_x_not_the_species_row():
    g = mm.group_terms("CH_Gband")
    assert g["table3_row"] == "CH (A-X)"
    assert g["shift_3d_dex"] == pytest.approx(8.459 - 8.393)
    assert float(mm._row("CH").model_form_3D_minus_MARCS) == pytest.approx(0.035)


@pytest.mark.parametrize("indicator", [k for k, v in mm.INDICATOR_GROUP.items() if v[1]])
def test_table3_shift_matches_the_per_line_cds_group(indicator):
    """The transcribed Table 3 row and Amarsi's own per-line table are the same group."""
    species, system, dv = mm.INDICATOR_GROUP[indicator][1]
    d = mm._lines()
    g = d[(d.species == species) & (d.system == system)]
    if dv is not None:
        g = g[g.delta_nu.abs() == dv]
    cds = g.abundance_3d.mean() - g.abundance_marcs.mean()
    assert mm.group_terms(indicator)["shift_3d_dex"] == pytest.approx(cds, abs=0.003)


def test_a_molecular_group_charges_amarsis_atmosphere_term():
    d = doc(measurement=("MEASURED", 0.0), model_atmosphere=("MEASURED", 0.001),
            nlte=("N/A", None), molecular_coupling=("MEASURED", 0.0))
    g = mm.group_terms("CH_Gband")
    want = max(abs(g["shift_3d_dex"]) / 2, g["atmos_sem_dex"])
    assert asplund_sigma(d, "CH", "CH_Gband") == pytest.approx(want)
    assert asplund_sigma(d, "CH") == pytest.approx(0.035 / 2)      # no group: species row


@pytest.mark.parametrize("element", ["C", "N", "O"])
def test_1d_cno_products_carry_their_3d_shift_and_keep_their_a(element):
    products = json.loads((ROOT / f"data/products/solar/{element}.json").read_text())["products"]
    for p in products:
        g = mm.group_terms(mm.indicator_key(p)) if p.get("treatment") == "1D-LTE" else None
        if g is None:
            assert "asplund_3d_shift" not in p
        else:
            assert p["asplund_3d_shift"] == g["shift_3d_dex"]
            assert p["asplund_3d_shift_basis"].startswith(mm.TABLE3_SOURCE)
