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
        "docs/reproduction/workflows.md",
        "docs/reproducibility.md",
        "docs/troubleshooting.md",
    ):
        assert (ROOT / relative).is_file(), relative
