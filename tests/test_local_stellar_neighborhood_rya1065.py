import json
import math
from pathlib import Path

import pytest

from scripts.generate_local_stellar_neighborhood import SCHEMA_VERSION, generate, validate
from pipeline.system_catalog import load_catalog

ROOT = Path(__file__).resolve().parents[1]


def registry_site(tmp_path: Path) -> Path:
    for row in load_catalog():
        route = tmp_path / "systems" / row["website_slug"] / "index.html"
        route.parent.mkdir(parents=True, exist_ok=True)
        route.write_text(f'<title>{row["system_name"]}</title>', encoding="utf-8")
    return tmp_path


def test_product_is_data_driven_finite_and_geometrically_truthful(tmp_path):
    product = generate(registry_site(tmp_path))
    assert product["schema_version"] == SCHEMA_VERSION
    assert len(product["targets"]) == 19
    assert {t["name"] for t in product["targets"]} == {
        row["system_name"] for row in __import__("pipeline.system_catalog", fromlist=["load_catalog"]).load_catalog()
    }
    for target in product["targets"]:
        xyz = target["cartesian_pc"].values()
        assert math.isclose(math.sqrt(sum(v * v for v in xyz)), target["distance_pc"], abs_tol=2e-6)


def test_every_registry_target_has_a_resolving_canonical_page(tmp_path):
    site = registry_site(tmp_path)
    product = generate(site)
    assert all(t["publication_status"] == "published" for t in product["targets"])
    for target in product["targets"]:
        assert target["url"] == f'/systems/{target["slug"]}/'
        assert (site / target["url"].lstrip("/") / "index.html").is_file()


def test_validator_rejects_a_broken_published_route(tmp_path):
    site = registry_site(tmp_path)
    product = generate(site)
    procyon = next(t for t in product["targets"] if t["id"] == "procyon")
    procyon["url"] = "/systems/not-procyon/"
    with pytest.raises(ValueError, match="published route does not resolve"):
        validate(product, site)
