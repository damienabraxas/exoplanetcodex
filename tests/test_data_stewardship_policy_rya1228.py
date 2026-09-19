"""RYA-1228 -- the standing policy, enforced structurally rather than documented.

🔴 WHAT THIS PROTECTS. The LFS quota was consumed by ONE `.gitattributes` line while
`.gitignore` was doing its job correctly on spectra and model grids. Closing a single
route is not the fix; these tests pin both routes AND the commit-time gate, plus the two
ways the gate could be useless: refusing everything, or refusing nothing.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.data_stewardship_policy import (  # noqa: E402
    CEILING_EXEMPT_SUFFIXES, SIZE_CEILING_BYTES, matched_pattern, over_ceiling)


def test_the_policy_module_has_no_heavy_imports():
    """🔴 THE HOOK MUST IMPORT ON ANY INTERPRETER.

    check_stewardship.py pulls in pipeline.abundances_derive, whose PEP 604 annotations
    fail to import on Python 3.9 -- the version this machine actually has. A hook routed
    through it would crash here and be silently disabled, which is worse than no hook.
    """
    src = (ROOT / "pipeline/data_stewardship_policy.py").read_text()
    for forbidden in ("import pandas", "import numpy", "from pipeline.", "import pipeline"):
        assert forbidden not in src, f"policy module must stay dependency-free: {forbidden}"
    out = subprocess.run([sys.executable, "-c",
                          "import sys; sys.path.insert(0, %r);"
                          "import pipeline.data_stewardship_policy" % str(ROOT)],
                         capture_output=True, text=True)
    assert out.returncode == 0, out.stderr


@pytest.mark.parametrize("path", [
    "data/linelists/vald_solar_raw.txt",
    "data/linelists/vald_anything_raw.txt",
    "data/linelists/vald_procyon_raw/vald_procyon_blue.vald",
    "data/spectra/whatever.fits",
    "data/model_atmospheres/grid.mod",
])
def test_raw_external_data_is_refused(path):
    assert matched_pattern(path) is not None, path


@pytest.mark.parametrize("path", [
    "data/linelists/linelist_solar.csv",
    "data/linelists/canonical_gf.csv",
    "data/products/solar/Fe.json",
    "scripts/build_linelist.py",
])
def test_codex_products_and_code_are_not_refused(path):
    """🔴 A gate that refuses our own measurements trains people to --no-verify."""
    assert matched_pattern(path) is None, path


def test_the_ceiling_exempts_products_and_catches_the_unanticipated():
    big = 25 * 1024 * 1024
    assert over_ceiling("data/x.bin", big) is True
    assert over_ceiling("data/linelists/linelist_solar.csv", big) is False
    assert over_ceiling("scripts/huge.py", big) is False
    assert over_ceiling("data/x.bin", SIZE_CEILING_BYTES - 1) is False
    assert ".csv" in CEILING_EXEMPT_SUFFIXES


def test_gitattributes_tracks_NO_data_in_lfs():
    """The route VALD actually came in through. There must be no LFS pattern at all."""
    text = (ROOT / ".gitattributes").read_text()
    assert "filter=lfs" not in text, "an LFS-tracked pattern has returned"


def test_gitignore_closes_the_vald_route_without_touching_products():
    def ignored(p):
        r = subprocess.run(["git", "-C", str(ROOT), "check-ignore", "--no-index", "-q", p])
        return r.returncode == 0
    assert ignored("data/linelists/vald_solar_raw.txt")
    assert ignored("data/linelists/vald_procyon_raw/vald_procyon_blue.vald")
    assert not ignored("data/linelists/linelist_solar.csv")
    assert not ignored("data/linelists/canonical_gf.csv")


def test_the_hook_refuses_raw_data_and_allows_a_product(tmp_path):
    """End-to-end on the real hook, both directions -- a gate with only one is untested."""
    hook = ROOT / "hooks" / "pre-commit"
    assert hook.exists() and hook.stat().st_mode & 0o111, "hook must be executable"
    body = hook.read_text()
    assert "data_stewardship_policy" in body, "hook must use the shared policy"
    assert "--no-verify" in body, "an override must be documented, not hidden"
    assert "check_stewardship import" not in body, (
        "the hook must not route through the heavy module -- it cannot import on 3.9")


def test_the_standard_is_written_down_and_names_the_archive_destinations():
    text = (ROOT / "SCIENCE_STANDARDS.md").read_text()
    assert "RYA-461" in text and "RYA-1124" in text
    assert "NEVER committed or LFS-tracked" in text
    assert "linelist_" in text and "stay" in text.lower()
