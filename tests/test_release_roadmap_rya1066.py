from scripts.generate_release_roadmap import generate


def test_release_progress_is_derived_from_current_publishable_products():
    product = generate()
    assert product["canonical_product_count"] == 27
    releases = {release["id"]: release for release in product["releases"]}
    assert releases["foundation"]["progress_fraction"] == 1.0
    solar = releases["alpha"]["targets"][0]
    # RYA-1270: pinned to the store, not to a count -- the count of 2 outlived C/N/O/Si/Al.
    from scripts.generate_release_roadmap import PRODUCTS
    on_disk = {p.stem for p in (PRODUCTS / "solar").glob("*.json")}
    assert {"Fe", "Fe-II"} <= set(solar["complete_product_ids"])
    assert {i for i in solar["complete_product_ids"] if i != "Fe-II"} <= on_disk
    assert solar["complete_products"] == len(solar["complete_product_ids"])
    assert solar["total_products"] == 27
    note = releases["alpha"]["notes"][0]
    assert note.startswith("Published for the Sun: Fe (I and II)")
    assert all(i in note for i in solar["complete_product_ids"] if i != "Fe-II")
    assert [target["complete_products"] for target in releases["beta"]["targets"]] == [0, 0]


def test_public_progress_declares_its_gate_and_excludes_noncurrent_channels():
    product = generate()
    assert "current `products` array" in product["sources"]["publication_gate"]
    for release in product["releases"]:
        for target in release.get("targets", []):
            assert target["gate"] == "current_publishable_product"
            assert target["complete_products"] == len(target["complete_product_ids"])
