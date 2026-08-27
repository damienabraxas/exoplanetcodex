from scripts.generate_release_roadmap import generate


def test_release_progress_is_derived_from_current_publishable_products():
    product = generate()
    assert product["canonical_product_count"] == 27
    releases = {release["id"]: release for release in product["releases"]}
    assert releases["foundation"]["progress_fraction"] == 1.0
    solar = releases["alpha"]["targets"][0]
    assert solar["complete_product_ids"] == ["Fe", "Fe-II"]
    assert solar["complete_products"] == 2
    assert solar["total_products"] == 27
    assert [target["complete_products"] for target in releases["beta"]["targets"]] == [0, 0]


def test_public_progress_declares_its_gate_and_excludes_noncurrent_channels():
    product = generate()
    assert "current `products` array" in product["sources"]["publication_gate"]
    for release in product["releases"]:
        for target in release.get("targets", []):
            assert target["gate"] == "current_publishable_product"
            assert target["complete_products"] == len(target["complete_product_ids"])
