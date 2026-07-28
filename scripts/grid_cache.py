#!/usr/bin/env python3
"""
scripts/grid_cache.py — RYA-540: persistent md5-pinned NLTE grid cache.

Governing rule (Ryan, standing): Sirius HOLDS the working set locally. Re-downloading an
already-acquired grid is a DEFECT. So: download a Gerber TS departure grid ONCE, keep the
unzipped .bin on local disk (canonical store, RYA-461), and on every subsequent need VERIFY
by md5 and reuse — NEVER re-download a grid whose md5 matches. There is NO free-after-gate
(that RYA-534 policy is replaced); eviction is an explicit admin action only.

Layout (RYA-540 STEP-0): grids live on the M.2 at /srv/codex/grids, exposed at the deck's
hardcoded path /mnt/codex-data/grids via a compat symlink. This module uses the canonical
gerber_ts dir and is agnostic to which physical drive backs it.

Integrity + capacity:
  * On acquire: download the .bin.zip -> verify against the ZIP md5 pinned in the element's
    provenance JSON (data/nlte_grids/gerber_ts/<El>_gerber2023.prov.json) -> unzip -> record
    the unzipped .bin md5 in the cache index -> keep the .bin, delete the .zip.
  * On load: md5 the cached .bin against the recorded index md5. MISMATCH RAISES (corrupt /
    truncated = the RYA-534 concurrent-unzip class) and re-fetches THAT grid only.
  * Capacity: before a download, check free space vs the expected grid size; if it will not
    fit, LOUD-FAIL with the shortfall and point at RYA-477 (external SSD / Orion). A capacity
    problem is a hardware decision, never a silent re-download loop.

No silent fallbacks anywhere: a miss downloads (once); a mismatch raises; a full disk raises.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

# Canonical gerber_ts store (deck path; symlink -> M.2 /srv/codex/grids/nlte/gerber_ts).
GERBER_DIR = Path(os.environ.get(
    "GERBER_TS_DIR", "/mnt/codex-data/grids/nlte/gerber_ts"))
CACHE_INDEX = GERBER_DIR / "_cache_index.json"
# Provenance JSONs live in the repo (committed). Resolve from repo roots, like the deck.
_PROV_CANDIDATES = [
    os.environ.get("TSGERBER_REPO_ROOT", ""),
    str(Path(__file__).resolve().parents[1]),
    "/mnt/codex-data/codex/rya519", "/mnt/codex-data/codex/repo",
]


def _prov_path(element: str) -> Path:
    # Primary: co-located with the grids (self-contained cache on Sirius). Then repo roots.
    local = GERBER_DIR / f"{element}_gerber2023.prov.json"
    if local.exists():
        return local
    for root in _PROV_CANDIDATES:
        if not root:
            continue
        p = Path(root) / "data" / "nlte_grids" / "gerber_ts" / f"{element}_gerber2023.prov.json"
        if p.exists():
            return p
    raise SystemExit(f"no provenance JSON for {element} (looked in {GERBER_DIR} + {_PROV_CANDIDATES}); "
                     f"REFUSING to fetch a grid without pinned url+md5.")


def _md5(path: Path, buf=1 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(buf), b""):
            h.update(chunk)
    return h.hexdigest()


def _load_index() -> dict:
    if CACHE_INDEX.exists():
        return json.loads(CACHE_INDEX.read_text())
    return {}


def _save_index(idx: dict):
    CACHE_INDEX.write_text(json.dumps(idx, indent=2, sort_keys=True) + "\n")


def _free_bytes(path: Path) -> int:
    return shutil.disk_usage(path).free


def _prov(element: str) -> dict:
    d = json.loads(_prov_path(element).read_text())
    # provenance schema pins the .bin.ZIP (name ends .bin.zip; md5/bytes are the zip's; the
    # unzipped grid is name without .zip). The key is inconsistent across tickets: RYA-533
    # (Na) used "grid_1d_bin"; RYA-534 (the other 10) used "grid_1d_bin_zip". Accept either.
    files = d["files"]
    z = files.get("grid_1d_bin_zip") or files.get("grid_1d_bin")
    if z is None:
        raise SystemExit(f"{element}: provenance has no grid_1d_bin[_zip] entry — cannot fetch.")
    return {"zip_name": z["name"], "zip_url": z["url"], "zip_md5": z["md5"],
            "zip_bytes": z.get("bytes"),
            "bin_name": z["name"][:-4] if z["name"].endswith(".zip") else z["name"] + ".bin"}


def status(element: str | None = None):
    idx = _load_index()
    print(f"cache dir: {GERBER_DIR}  (backing device: {_device(GERBER_DIR)})")
    print(f"free here: {_free_bytes(GERBER_DIR)/1e9:.1f} GB")
    els = [element] if element else sorted(idx)
    if not els:
        print("cache index empty (no grids acquired yet).")
    for el in els:
        rec = idx.get(el)
        binf = GERBER_DIR / rec["bin_name"] if rec else None
        present = binf.exists() if binf else False
        print(f"  {el}: index={'yes' if rec else 'no'} bin_present={present}"
              + (f" bin_md5={rec['bin_md5'][:12]} bytes={rec.get('bin_bytes','?')}" if rec else ""))


def _device(path: Path) -> str:
    try:
        return subprocess.run(["df", "--output=source", str(path)], capture_output=True,
                              text=True).stdout.splitlines()[-1].strip()
    except Exception:
        return "?"


def ensure_grid(element: str, verify: bool = True) -> Path:
    """Return the local path to element's unzipped Gerber departure .bin, downloading ONLY on
    a cache miss/mismatch. Never deletes. RAISES on md5 mismatch or insufficient capacity."""
    GERBER_DIR.mkdir(parents=True, exist_ok=True)
    prov = _prov(element)
    binf = GERBER_DIR / prov["bin_name"]
    idx = _load_index()
    rec = idx.get(element)

    # CACHE HIT: bin present + (fast) recorded, verify md5 if asked.
    if binf.exists() and rec and rec.get("bin_md5"):
        if verify:
            actual = _md5(binf)
            if actual != rec["bin_md5"]:
                print(f"[grid_cache] {element}: md5 MISMATCH ({actual[:12]} != "
                      f"{rec['bin_md5'][:12]}) — corrupt/truncated; re-fetching THIS grid only.")
                binf.unlink()
            else:
                print(f"[grid_cache] {element}: cache HIT (md5 {actual[:12]}), skipping download.")
                return binf
        else:
            print(f"[grid_cache] {element}: cache hit (bin present; md5 not re-checked).")
            return binf

    # MISS -> download once, with a capacity guard.
    zipf = GERBER_DIR / prov["zip_name"]
    # Pre-download capacity guard: we know the exact zip size (prov bytes). Need room for the
    # zip AND the unzipped .bin simultaneously (.bin ~= zip; Gerber grids are already ~un-
    # compressible binaries). Reserve zip + ~1.5x for the bin as a floor. LOUD-FAIL if short —
    # a capacity shortfall is a hardware decision (RYA-477), never a re-download loop.
    free = _free_bytes(GERBER_DIR)
    zbytes = prov.get("zip_bytes") or (1 << 30)
    need = int(zbytes * 2.5)
    if free < need:
        raise SystemExit(
            f"[grid_cache] CAPACITY: fetching {element} needs ~{need/1e9:.0f} GB (zip {zbytes/1e9:.0f} "
            f"GB + unzip headroom) but only {free/1e9:.1f} GB free at {GERBER_DIR}. ADD STORAGE "
            f"(external SSD / Orion, RYA-477); do NOT thrash re-downloads. Working set exceeds "
            f"disk = a hardware decision, not a re-download.")
    print(f"[grid_cache] {element}: cache MISS -> downloading {prov['zip_name']} "
          f"(Sirius-only, RYA-526) from {prov['zip_url']}")
    r = subprocess.run(["curl", "-sL", "-C", "-", "-o", str(zipf), prov["zip_url"]])
    if r.returncode != 0 or not zipf.exists():
        raise SystemExit(f"[grid_cache] {element}: download failed (curl rc={r.returncode}).")
    zmd5 = _md5(zipf)
    if zmd5 != prov["zip_md5"]:
        raise SystemExit(f"[grid_cache] {element}: ZIP md5 mismatch ({zmd5} != {prov['zip_md5']}) "
                         f"— truncated/corrupt download; NOT unzipping. Re-run to resume/re-fetch.")
    # capacity for unzip: need room for the .bin (~3-4x zip) on top of the zip
    zsize = zipf.stat().st_size
    if _free_bytes(GERBER_DIR) < 4 * zsize:
        raise SystemExit(
            f"[grid_cache] CAPACITY: unzip of {element} needs ~{4*zsize/1e9:.0f} GB but only "
            f"{_free_bytes(GERBER_DIR)/1e9:.1f} GB free — ADD STORAGE (RYA-477). Zip kept for resume.")
    print(f"[grid_cache] {element}: zip md5 OK; unzipping...")
    if subprocess.run(["unzip", "-o", "-q", str(zipf), "-d", str(GERBER_DIR)]).returncode != 0:
        raise SystemExit(f"[grid_cache] {element}: unzip failed (disk? corrupt?). Zip kept.")
    if not binf.exists():
        # find the produced .bin
        cands = sorted(GERBER_DIR.glob(f"NLTEgrid*{element.upper()}*.bin")) + \
                sorted(GERBER_DIR.glob(f"NLTEgrid*{element}*.bin"))
        if not cands:
            raise SystemExit(f"[grid_cache] {element}: no .bin after unzip — RAISE (no default).")
        binf = cands[0]
    bmd5 = _md5(binf)
    idx[element] = {"bin_name": binf.name, "bin_md5": bmd5, "bin_bytes": binf.stat().st_size,
                    "zip_md5": prov["zip_md5"], "zip_name": prov["zip_name"]}
    _save_index(idx)
    zipf.unlink()   # keep the .bin (ready to use); drop the re-downloadable .zip
    print(f"[grid_cache] {element}: acquired + cached (.bin md5 {bmd5[:12]}, "
          f"{binf.stat().st_size/1e9:.1f} GB). Zip removed; .bin retained (no free-after-gate).")
    return binf


def evict(element: str):
    """Explicit admin disk-reclaim ONLY (never automated / never mid-survey)."""
    idx = _load_index()
    rec = idx.get(element)
    if not rec:
        print(f"[grid_cache] {element}: not in cache."); return
    binf = GERBER_DIR / rec["bin_name"]
    if binf.exists():
        binf.unlink()
    idx.pop(element, None)
    _save_index(idx)
    print(f"[grid_cache] {element}: EVICTED (admin). Re-acquire via ensure_grid (md5-pinned).")


def _cli():
    ap = argparse.ArgumentParser(description="RYA-540 persistent md5-pinned NLTE grid cache")
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--ensure", metavar="EL", help="acquire (cache-or-download) a grid")
    ap.add_argument("--verify", metavar="EL", help="md5-verify a cached grid, re-fetch on mismatch")
    ap.add_argument("--evict", metavar="EL", help="ADMIN disk-reclaim (never automated)")
    a = ap.parse_args()
    if a.evict:
        evict(a.evict)
    elif a.ensure:
        print(ensure_grid(a.ensure, verify=True))
    elif a.verify:
        print(ensure_grid(a.verify, verify=True))
    else:
        status(a.status if isinstance(a.status, str) else None)


if __name__ == "__main__":
    _cli()
