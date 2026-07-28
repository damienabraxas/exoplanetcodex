"""Keep the documented installation and quick-start entry points executable."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_documented_entry_points_exist():
    for relative in (
        "scripts/validate_installation.py",
        "scripts/quickstart_example.py",
        "docs/setup/environment.md",
        "docs/models/assets.md",
        "docs/data/spectra.md",
        "docs/data/instruments.md",
        "docs/reproduction/workflows.md",
        "docs/reproducibility.md",
        "docs/troubleshooting.md",
        "docs/references.md",
        "data/catalog/system_catalog.csv",
        "data/catalog/instrument_catalog.csv",
    ):
        assert (ROOT / relative).is_file(), relative
