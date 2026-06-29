#!/usr/bin/env python3
"""
scripts/stage_amarsi_mn_grid_rya476.py
======================================
RYA-476 — source + stage the Amarsi GALAH Mn NLTE departure grid
(`nlte_Mn_scatt_pysme.grd`) so RYA-473's HFS-resolved synthesis can use the
TRIPLET-EXACT δ on 6013/16/21 instead of the high-EP-reference MPIA +0.108.

The grid is a multi-GB binary — gitignored, never committed (the RYA-410 disk-safe
pattern). This stages it from the cited Zenodo deposit, md5-verified:
  Amarsi, Lind, Osorio et al. 2020, A&A 642, A62 (GALAH non-LTE departure grids)
  Zenodo record 3982506 (v3, 2020-06-11), DOI 10.5281/zenodo.3982506
  nlte_Mn_scatt_pysme.tar.gz  md5 ba01596031b054386a79a992141ebf72 (RYA-410 verified)

Disk-safe: download → md5-verify the tar → extract the .grd → REMOVE the tar (the
.grd is KEPT on disk; re-fetchable from the verified tar md5 if freed later). No
binary is committed; the committed products are the prov JSON + the per-line δ that
RYA-473 emits once the grid is live.

    python -m scripts.stage_amarsi_mn_grid_rya476            # download + extract + verify
    python -m scripts.stage_amarsi_mn_grid_rya476 --check    # report presence only, no fetch
"""
from __future__ import annotations

import argparse
import hashlib
import sys
import tarfile
import urllib.request
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
GRID_DIR = _REPO / 'data' / 'nlte_grids' / 'amarsi_galah'
TAR_NAME = 'nlte_Mn_scatt_pysme.tar.gz'
GRD_NAME = 'nlte_Mn_scatt_pysme.grd'
TAR_MD5 = 'ba01596031b054386a79a992141ebf72'   # RYA-410-verified Zenodo tar
URL = 'https://zenodo.org/records/3982506/files/nlte_Mn_scatt_pysme.tar.gz?download=1'


def _md5(path: Path, chunk=1 << 20) -> str:
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for blk in iter(lambda: fh.read(chunk), b''):
            h.update(blk)
    return h.hexdigest()


def check() -> bool:
    grd = GRID_DIR / GRD_NAME
    if grd.exists():
        print(f"  PRESENT: {grd.relative_to(_REPO)} ({grd.stat().st_size/1e9:.2f} GB)")
        return True
    print(f"  ABSENT:  {grd.relative_to(_REPO)} — fetch from Zenodo 3982506 (tar md5 {TAR_MD5})")
    return False


def stage():
    GRID_DIR.mkdir(parents=True, exist_ok=True)
    grd = GRID_DIR / GRD_NAME
    if grd.exists():
        print(f"  already staged: {grd.relative_to(_REPO)} ({grd.stat().st_size/1e9:.2f} GB)")
        return grd
    tar = GRID_DIR / TAR_NAME
    print(f"  downloading {TAR_NAME} from Zenodo 3982506 ...")
    urllib.request.urlretrieve(URL, tar)
    print(f"  downloaded {tar.stat().st_size/1e9:.2f} GB; verifying md5 ...")
    got = _md5(tar)
    if got != TAR_MD5:
        tar.unlink(missing_ok=True)
        raise SystemExit(f"FATAL: md5 mismatch (got {got}, want {TAR_MD5}) — refusing, tar removed")
    print(f"  md5 PASS ({got})")
    print(f"  extracting {GRD_NAME} ...")
    with tarfile.open(tar) as t:
        members = [m for m in t.getmembers() if m.name.endswith(GRD_NAME) or m.name.endswith('.grd')]
        if not members:
            raise SystemExit(f"FATAL: no .grd member in {TAR_NAME}")
        for m in members:
            m.name = Path(m.name).name           # flatten into GRID_DIR
            t.extract(m, GRID_DIR)
    if not grd.exists():
        raise SystemExit(f"FATAL: extraction did not produce {grd}")
    print(f"  extracted {grd.relative_to(_REPO)} ({grd.stat().st_size/1e9:.2f} GB)")
    tar.unlink()
    print(f"  removed {TAR_NAME} (re-fetchable from verified md5); .grd KEPT on disk")
    return grd


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true', help='report presence only, no fetch')
    args = ap.parse_args(argv)
    print(f"\n  RYA-476 — Amarsi GALAH Mn NLTE grid staging\n  {'='*52}")
    if args.check:
        sys.exit(0 if check() else 1)
    stage()


if __name__ == '__main__':
    main()
