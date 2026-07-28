from pathlib import Path

from data.catalog.instruments import (
    load_catalog, render_markdown_table, validate_all, validate_catalog,
    validate_holdings, validate_modes,
)


def test_instrument_catalog_contract():
    assert validate_catalog() == []
    assert validate_modes() == []
    assert validate_holdings() == []
    assert validate_all() == []


def test_required_multiwavelength_sources_are_registered():
    ids = {row["instrument_id"] for row in load_catalog()}
    assert {
        "harps", "harps_n", "espresso", "uves", "feros", "crires_plus",
        "hst_stis", "hst_cos", "hires", "nirspec", "kpf", "spirou",
        "espadons", "narval", "phoenix", "ghost", "chiron",
        "flames_giraffe", "kpno_solar_atlas", "calspec_solar",
    } <= ids


def test_capability_is_not_per_system_holdings():
    rows = load_catalog()
    assert all("file_count" not in row for row in rows)
    assert all("target_count" not in row for row in rows)


def test_readme_instrument_table_is_generated_from_catalog():
    readme = Path(__file__).resolve().parents[1] / "README.md"
    text = readme.read_text()
    start = "<!-- instrument-table:start -->"
    end = "<!-- instrument-table:end -->"
    embedded = text.split(start, 1)[1].split(end, 1)[0].strip()
    assert embedded == render_markdown_table()
