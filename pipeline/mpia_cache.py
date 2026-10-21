"""On-disk cache + retry in front of the live MPIA NLTE query (RYA-1232).

WHY. ENGINE-A for Fe is a live per-line query to nlte.mpia.de. The orchestrator issues it from
every cell and every budget leg in parallel: on 2026-10-20 70 of 73 requests came back 504 /
truncated and ENGINE-A was "NOT PRODUCED" on every Asplund-set product. MPIA's answer for a
line at a grid node does not change between requests, so each (species, node, wavelength) is
fetched ONCE and kept; a failed request is retried with backoff instead of dropping the leg.

`install()` replaces `scripts.build_nlte_grids_mpia._submit_batch` (the function
derive_band_products' engine_a_delta imports at call time) with a cached version. Values are
returned exactly as MPIA served them; only lines MPIA actually served are stored, so an
unserved line is asked again next time. One JSON file per answer: concurrent workers never
race on a shared file. Dir: $CODEX_MPIA_CACHE_DIR, default data/cache/mpia (gitignored).
"""
from __future__ import annotations

import hashlib
import json
import os
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RETRIES = 6
_INSTALLED = False


def cache_dir() -> Path:
    d = Path(os.environ.get("CODEX_MPIA_CACHE_DIR") or ROOT / "data" / "cache" / "mpia")
    d.mkdir(parents=True, exist_ok=True)
    return d


def _path(code: str, node: dict, wave: float) -> Path:
    tag = hashlib.sha1(json.dumps(node, sort_keys=True).encode()).hexdigest()[:12]
    return cache_dir() / f"{code}_{tag}_{float(wave):.4f}.json"


def install() -> bool:
    global _INSTALLED
    if _INSTALLED:
        return True
    import sys
    sys.path.insert(0, str(ROOT))
    import scripts.build_nlte_grids_mpia as M
    orig = M._submit_batch

    def cached(lines_angstrom, species_code, grid_params):
        out: dict = {g["name"]: {} for g in grid_params}
        missing = []
        for w in lines_angstrom:
            hit = True
            for g in grid_params:
                f = _path(species_code, g, w)
                if f.exists():
                    rec = json.loads(f.read_text())
                    out[g["name"]].update({float(k): v for k, v in rec.items()})
                else:
                    hit = False
            if not hit:
                missing.append(w)
        if missing:
            last = None
            for k in range(RETRIES):
                try:
                    r = orig(missing, species_code, grid_params)
                    break
                except Exception as e:                       # noqa: BLE001
                    last = e
                    time.sleep(min(60, 5 * 2 ** k))
            else:
                raise last
            for g in grid_params:
                served = r.get(g["name"], {}) or {}
                out[g["name"]].update(served)
                for w in missing:
                    if not served:
                        continue
                    kk = min(served, key=lambda x: abs(x - w))
                    if abs(kk - w) < 0.6 and served[kk] is not None:
                        fd, tp = tempfile.mkstemp(dir=cache_dir(), suffix=".tmp")
                        with os.fdopen(fd, "w") as fh:
                            json.dump({str(kk): served[kk]}, fh)
                        os.replace(tp, _path(species_code, g, w))
        return out

    M._submit_batch = cached
    _INSTALLED = True
    return True
