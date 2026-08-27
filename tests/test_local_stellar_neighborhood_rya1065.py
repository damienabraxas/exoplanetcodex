import json
import math
from pathlib import Path

import pytest

from scripts.generate_local_stellar_neighborhood import SCHEMA_VERSION, generate, validate

ROOT = Path(__file__).resolve().parents[1]
SITE = Path("/Users/ryanschmitt/codex/rya1065-site")


def test_product_is_data_driven_finite_and_geometrically_truthful():
    product = generate(SITE)
    assert product["schema_version"] == SCHEMA_VERSION
    assert len(product["targets"]) == 19
    assert {t["name"] for t in product["targets"]} == {
        row["system_name"] for row in __import__("pipeline.system_catalog", fromlist=["load_catalog"]).load_catalog()
    }
    for target in product["targets"]:
        xyz = target["cartesian_pc"].values()
        assert math.isclose(math.sqrt(sum(v * v for v in xyz)), target["distance_pc"], abs_tol=2e-6)


def test_links_exist_or_are_explicitly_unpublished():
    product = generate(SITE)
    published = {t["id"] for t in product["targets"] if t["publication_status"] == "published"}
    assert published == {
        "solar", "alpha_cen_a", "55cnc_a", "hd209458", "hd189733", "tau_ceti",
        "eps_eri", "gliese581", "hd89307", "proxima", "61-vir", "51-peg",
    }
    assert all(t["url"] is None for t in product["targets"] if t["id"] not in published)


def test_validator_rejects_fabricated_unpublished_route():
    product = generate(SITE)
    procyon = next(t for t in product["targets"] if t["id"] == "procyon")
    procyon["url"] = "/systems/procyon/"
    with pytest.raises(ValueError, match="fabricated URL"):
        validate(product, SITE)
