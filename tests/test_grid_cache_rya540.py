"""
tests/test_grid_cache_rya540.py
===============================
RYA-540 anti-drift guard for the persistent Gerber-TS grid cache. CI-safe: it reads only
the committed provenance JSONs + the grid_cache module logic — no grids, no network, no
Sirius. Guards that every Gerber grid stays md5-pinned and uniquely named so a cached
.bin can always be verified on load (the defect RYA-540 exists to prevent).
"""
import importlib.util
import json
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[1]
_PROV_DIR = _REPO / "data" / "nlte_grids" / "gerber_ts"
_PROVS = sorted(_PROV_DIR.glob("*_gerber2023.prov.json"))

# Load scripts/grid_cache.py by path (scripts/ is not a package).
_spec = importlib.util.spec_from_file_location(
    "grid_cache", _REPO / "scripts" / "grid_cache.py")
grid_cache = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(grid_cache)

_ELEMENTS = [p.name.split("_")[0] for p in _PROVS]


def test_gerber_provs_present():
    """The Family-A Gerber set is fully provenanced (>=11 elements)."""
    assert len(_PROVS) >= 11, f"expected >=11 gerber provs, found {len(_PROVS)}: {_ELEMENTS}"


@pytest.mark.parametrize("prov", _PROVS, ids=_ELEMENTS)
def test_prov_pins_zip_md5_and_bytes(prov):
    """Every grid pins a 32-hex zip md5, a .bin.zip name, and an integer byte size
    (the RYA-540 capacity-guard backfill) — no unpinned re-fetch is possible."""
    d = json.loads(prov.read_text())["files"]
    z = d.get("grid_1d_bin_zip") or d.get("grid_1d_bin")
    assert z is not None, f"{prov.name}: no grid_1d_bin[_zip] entry"
    assert z["name"].endswith(".bin.zip"), f"{prov.name}: name is not a .bin.zip"
    assert len(z["md5"]) == 32 and all(c in "0123456789abcdef" for c in z["md5"]), \
        f"{prov.name}: zip md5 not 32-hex"
    assert isinstance(z.get("bytes"), int) and z["bytes"] > 0, \
        f"{prov.name}: missing/invalid exact zip bytes (RYA-540 backfill)"


@pytest.mark.parametrize("el", _ELEMENTS)
def test_grid_cache_reads_every_prov(el):
    """grid_cache._prov resolves every element (schema-tolerant: grid_1d_bin OR
    grid_1d_bin_zip) and derives a .bin name matching the pinned zip."""
    pv = grid_cache._prov(el)
    assert pv["bin_name"].endswith(".bin")
    assert pv["zip_name"] == pv["bin_name"] + ".zip"
    assert len(pv["zip_md5"]) == 32


def test_bin_names_and_md5s_unique():
    """No two grids share a .bin name or a zip md5 — a collision would let the cache
    serve the wrong grid or mask a drift."""
    provs = [grid_cache._prov(el) for el in _ELEMENTS]
    names = [p["bin_name"] for p in provs]
    md5s = [p["zip_md5"] for p in provs]
    assert len(set(names)) == len(names), f"duplicate .bin name: {names}"
    assert len(set(md5s)) == len(md5s), f"duplicate zip md5: {md5s}"
