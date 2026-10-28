"""RYA-1251: the published sigma follows Asplund, Amarsi & Grevesse 2021 Sect. 2.1."""
import json
import math
from pathlib import Path

import pytest

from pipeline.uncertainty_contract import ASPLUND_NLTE_FLOOR, asplund_sigma

ROOT = Path(__file__).resolve().parents[1]


def doc(**terms):
    return {"components": [{"name": n, "state": s, "sigma_dex": v} for n, (s, v) in terms.items()]}


def test_only_statistics_3d_and_nlte_enter():
    d = doc(measurement=("MEASURED", 0.02), model_atmosphere=("MEASURED", 0.01),
            nlte=("MEASURED", 0.05), transition_data=("MEASURED", 0.2),
            continuum=("MEASURED", 0.3), molecular_coupling=("N/A", None))
    assert asplund_sigma(d) == pytest.approx(math.sqrt(0.02**2 + 0.01**2 + 0.05**2))


def test_nlte_floor_applies_to_every_atomic_product():
    d = doc(measurement=("MEASURED", 0.0), model_atmosphere=("MEASURED", 0.0), nlte=("N/A", None))
    assert asplund_sigma(d) == pytest.approx(ASPLUND_NLTE_FLOOR)


def test_molecules_carry_no_nlte_term():
    d = doc(measurement=("MEASURED", 0.02), model_atmosphere=("MEASURED", 0.0),
            nlte=("N/A", None), molecular_coupling=("MEASURED", 0.08))
    assert asplund_sigma(d) == pytest.approx(0.02)


def test_hold_or_missing_measurement_gives_none():
    assert asplund_sigma(doc(measurement=("HOLD", None), nlte=("MEASURED", 0.03))) is None
    assert asplund_sigma(doc(measurement=("MEASURED", 0.01), nlte=("HOLD", None))) is None


@pytest.mark.parametrize("element", ["Fe", "C", "N", "O", "Si", "Al", "Na"])
def test_every_budgeted_product_is_stamped_with_its_own_recompute(element):
    products = json.loads((ROOT / f"data/products/solar/{element}.json").read_text())["products"]
    budgeted = [p for p in products if "uncertainty" in p]
    assert budgeted
    for p in budgeted:
        assert p["sigma_asplund"] == pytest.approx(asplund_sigma(p["uncertainty"]), abs=1e-12)
        assert p["sigma_asplund_basis"].startswith("Asplund, Amarsi & Grevesse 2021")
