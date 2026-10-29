"""RYA-1251: the published sigma follows Asplund, Amarsi & Grevesse 2021 Sect. 2.1."""
import json
import math
from pathlib import Path

import pytest

from pipeline.uncertainty_contract import ASPLUND_NLTE_FLOOR, asplund_sigma, molecule_of

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


def test_molecules_take_the_compositional_term_not_nlte():
    d = doc(measurement=("MEASURED", 0.02), model_atmosphere=("MEASURED", 0.0),
            nlte=("MEASURED", 0.2), molecular_coupling=("MEASURED", 0.08))
    assert asplund_sigma(d) == pytest.approx(math.hypot(0.02, 0.08))


def test_hold_or_missing_measurement_gives_none():
    assert asplund_sigma(doc(measurement=("HOLD", None), nlte=("MEASURED", 0.03))) is None
    assert asplund_sigma(doc(measurement=("MEASURED", 0.01), nlte=("HOLD", None))) is None


@pytest.mark.parametrize("element", ["Fe", "C", "N", "O", "Si", "Al", "Na", "Mg"])
def test_every_budgeted_product_is_stamped_with_its_own_recompute(element):
    products = json.loads((ROOT / f"data/products/solar/{element}.json").read_text())["products"]
    budgeted = [p for p in products if "uncertainty" in p]
    assert budgeted
    for p in budgeted:
        assert p["sigma_asplund"] == pytest.approx(asplund_sigma(p["uncertainty"], molecule_of(p)), abs=1e-12)
        assert p["sigma_asplund"] == pytest.approx(math.hypot(p["sigma_asplund_stat"],
                                                              p["sigma_asplund_syst"]), abs=1e-12)
        assert p["sigma_asplund_basis"].startswith("Asplund, Amarsi & Grevesse 2021")


def test_a_molecules_3d_term_is_at_least_half_amarsis_shift():
    from pipeline import molecular_model_form
    shift = abs(float(molecular_model_form._row("CH").model_form_3D_minus_MARCS))
    d = doc(measurement=("MEASURED", 0.0), model_atmosphere=("MEASURED", 0.001),
            nlte=("N/A", None), molecular_coupling=("MEASURED", 0.0))
    assert asplund_sigma(d, "CH") == pytest.approx(shift / 2)
    assert asplund_sigma(d) == pytest.approx(0.001)


def test_legacy_3dnlte_fe_anchor_gets_its_own_asplund_sigma():
    """The Fe top number (HARPS 3D NLTE, Asplund lines) has no budget, but its own files give the sigma."""
    import csv
    from pathlib import Path
    from pipeline.uncertainty_contract import legacy_3dnlte_asplund_terms
    doc = json.loads((ROOT / "data/products/solar/Fe.json").read_text())
    rows = [p for p in doc["products"] if "uncertainty" not in p
            and p["treatment"] == "ENGINE-A-3DNLTE" and p["selector"] == "ASPLUND_AGSS21"
            and p["instrument"] == "harps"]
    assert len(rows) == 1
    p = rows[0]
    per_line = ROOT / p["provenance"]["copied_to"].replace("_products.csv", "_per_line.csv")
    ab = [float(r["aberr"]) for r in csv.DictReader(per_line.open()) if r["in_domain"] == "True"]
    stat, syst = legacy_3dnlte_asplund_terms(p)
    assert stat == p["sigma_stat"]
    assert syst == pytest.approx(math.hypot(abs(sum(ab) / len(ab)) / 2, 0.03))
    assert p["sigma_asplund"] == pytest.approx(math.hypot(stat, syst))
    assert 0.03 < p["sigma_asplund"] < p["sigma_reported"]


def test_legacy_stamp_route_refuses_a_wrong_stamp_and_any_other_change():
    from pipeline.uncertainty_contract import UncertaintyError, assert_publication_feed
    doc = json.loads((ROOT / "data/products/solar/Fe.json").read_text())
    p = next(q for q in doc["products"] if "uncertainty" not in q and "sigma_asplund" in q)
    previous = {"products": [{k: v for k, v in p.items() if not k.startswith("sigma_asplund")}]}
    assert_publication_feed({"products": [p]}, previous=previous)
    with pytest.raises(UncertaintyError):
        assert_publication_feed({"products": [dict(p, sigma_asplund=0.01)]}, previous=previous)
    with pytest.raises(UncertaintyError):
        assert_publication_feed({"products": [dict(p, A=p["A"] + 0.001)]}, previous=previous)
