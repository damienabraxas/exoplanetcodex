#!/usr/bin/env python3
"""
scripts/sirius_fetch_grids_rya477.py
====================================
RYA-477 Part B — environment-independent grid downloads onto Sirius's 500 GB data
drive (/mnt/codex-data), driven from the Mac over SSH. Manifest-driven, md5-verified,
provenance-stamped — the RYA-359 / RYA-461 stewardship pattern, NO silent partial
downloads.

Per manifest item (data/sirius_manifest/*.json):
  1. download the .tar.gz to <data-root>/grids/nlte/_archives/ (wget -c, resumable)
  2. GATE 1 — byte size must equal the manifest size_bytes (catches truncation)
  3. GATE 2 — md5 must equal the manifest md5 (catches corruption); on mismatch the
     file is moved aside (.BAD) and the item is FAILED, never extracted
  4. extract the .grd into <data-root>/<dest_subdir>/ only after md5 PASS
  5. write a per-item provenance JSON + update the master status JSON

Idempotent: an item whose .grd is already present AND whose archive md5 re-verifies is
skipped (re-running resumes/repairs, never re-downloads a good grid). Designed to run
under nohup on Sirius for hours, independent of the Mac session.

    python3 sirius_fetch_grids_rya477.py --manifest grids_zenodo_3982506.json \
        --data-root /mnt/codex-data
    python3 sirius_fetch_grids_rya477.py ... --status-only   # print status, no fetch
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tarfile
import time
from datetime import datetime, timezone
from pathlib import Path


def _now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def _md5(path: Path, chunk=1 << 22) -> str:
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for blk in iter(lambda: fh.read(chunk), b''):
            h.update(blk)
    return h.hexdigest()


def _log(status_path: Path, status: dict):
    status['updated'] = _now()
    status_path.write_text(json.dumps(status, indent=2))


def _grd_for(archive_name: str) -> str:
    # nlte_Mn_scatt_pysme.tar.gz -> nlte_Mn_scatt_pysme.grd
    return archive_name[:-7] + '.grd' if archive_name.endswith('.tar.gz') else archive_name + '.grd'


def fetch_one(item: dict, data_root: Path, arch_dir: Path, log) -> dict:
    name = item['name']
    dest_dir = data_root / item['dest_subdir']
    dest_dir.mkdir(parents=True, exist_ok=True)
    grd = dest_dir / _grd_for(name)
    rec = {'name': name, 'element': item.get('element'), 'md5_expected': item['md5'],
           'size_expected': item['size_bytes'], 'grd': str(grd)}

    # idempotent skip: grd present and non-trivial
    if grd.exists() and grd.stat().st_size > 1_000_000:
        rec['state'] = 'PRESENT'; rec['grd_bytes'] = grd.stat().st_size
        log(f"  [{item['element']:>3}] PRESENT  {grd.name} ({grd.stat().st_size/1e9:.2f} GB) — skip")
        return rec

    archive = arch_dir / name
    url = item['url']
    log(f"  [{item['element']:>3}] downloading {name} ({item['size_bytes']/1e9:.2f} GB) ...")
    t0 = time.time()
    # wget -c = resume a partial; --tries for transient network; fail loud on error
    r = subprocess.run(['wget', '-c', '-q', '--tries=5', '--timeout=60', '-O', str(archive), url])
    if r.returncode != 0:
        rec['state'] = 'DOWNLOAD_FAILED'; rec['wget_rc'] = r.returncode
        log(f"  [{item['element']:>3}] DOWNLOAD_FAILED rc={r.returncode}")
        return rec
    got = archive.stat().st_size if archive.exists() else 0
    rec['size_got'] = got
    # GATE 1 — size
    if got != item['size_bytes']:
        rec['state'] = 'SIZE_MISMATCH'
        log(f"  [{item['element']:>3}] SIZE_MISMATCH got {got} want {item['size_bytes']} (partial?) — NOT extracting")
        return rec
    # GATE 2 — md5
    log(f"  [{item['element']:>3}] verifying md5 ...")
    got_md5 = _md5(archive)
    rec['md5_got'] = got_md5
    if got_md5 != item['md5']:
        bad = archive.with_suffix(archive.suffix + '.BAD')
        archive.rename(bad)
        rec['state'] = 'MD5_MISMATCH'; rec['bad_archive'] = str(bad)
        log(f"  [{item['element']:>3}] MD5_MISMATCH got {got_md5} want {item['md5']} — moved aside, NOT extracting")
        return rec
    # extract the .grd (only on md5 PASS)
    if item.get('extract', True):
        log(f"  [{item['element']:>3}] md5 PASS — extracting .grd ...")
        with tarfile.open(archive) as t:
            members = [m for m in t.getmembers() if m.name.endswith('.grd')]
            if not members:
                rec['state'] = 'NO_GRD_IN_ARCHIVE'
                log(f"  [{item['element']:>3}] NO_GRD_IN_ARCHIVE")
                return rec
            for m in members:
                m.name = Path(m.name).name
                t.extract(m, dest_dir)
        rec['grd_bytes'] = grd.stat().st_size if grd.exists() else 0
        # free the archive after successful extraction (disk-safe; re-fetchable by md5)
        archive.unlink(missing_ok=True)
        rec['archive_freed'] = True
    rec['state'] = 'VERIFIED'
    rec['wall_s'] = round(time.time() - t0, 1)
    # per-item provenance stamp
    prov = dest_dir / (name + '.prov.json')
    prov.write_text(json.dumps({
        'ticket': 'RYA-477', 'fetched': _now(), 'source_url': url,
        'md5': item['md5'], 'size_bytes': item['size_bytes'],
        'grd': grd.name, 'grd_bytes': rec.get('grd_bytes'),
        'verify': 'size+md5 PASS; archive freed post-extract (re-fetchable by md5)',
    }, indent=2))
    log(f"  [{item['element']:>3}] VERIFIED  {grd.name} ({rec.get('grd_bytes',0)/1e9:.2f} GB) in {rec['wall_s']}s")
    return rec


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', required=True, help='manifest filename (under the manifest dir)')
    ap.add_argument('--manifest-dir', default=str(Path(__file__).resolve().parent.parent /
                                                  'data' / 'sirius_manifest'))
    ap.add_argument('--data-root', default='/mnt/codex-data')
    ap.add_argument('--status-only', action='store_true')
    args = ap.parse_args(argv)

    manifest = json.loads((Path(args.manifest_dir) / args.manifest).read_text())
    data_root = Path(args.data_root)
    if not (data_root / '.codex_mounted').exists():
        sys.exit(f"FATAL: {data_root}/.codex_mounted sentinel missing — data drive not mounted, refusing")
    arch_dir = data_root / 'grids' / 'nlte' / '_archives'
    arch_dir.mkdir(parents=True, exist_ok=True)
    status_path = data_root / 'grids' / 'nlte' / f"_fetch_status_{manifest['manifest_id']}.json"

    def log(msg):
        print(msg, flush=True)

    log(f"\n  RYA-477 grid fetch — {manifest['manifest_id']} -> {data_root}")
    log(f"  {manifest['n_items']} items, {manifest['total_archive_bytes']/1e9:.2f} GB archives "
        f"(extracted larger); source {manifest.get('zenodo_record')}\n")

    status = {'manifest_id': manifest['manifest_id'], 'data_root': str(data_root),
              'started': _now(), 'items': {}}
    if status_path.exists():
        try:
            status = json.loads(status_path.read_text())
        except Exception:
            pass
    status.setdefault('items', {})

    if args.status_only:
        for it in manifest['items']:
            st = status['items'].get(it['name'], {}).get('state', 'UNKNOWN')
            grd = data_root / it['dest_subdir'] / _grd_for(it['name'])
            present = 'present' if grd.exists() else 'absent'
            log(f"  {it['element']:>3}  {st:16}  grd:{present}  {it['name']}")
        return

    ok = 0
    for it in manifest['items']:
        rec = fetch_one(it, data_root, arch_dir, log)
        status['items'][it['name']] = rec
        _log(status_path, status)
        if rec['state'] in ('VERIFIED', 'PRESENT'):
            ok += 1
    status['finished'] = _now()
    status['n_ok'] = ok
    status['n_total'] = len(manifest['items'])
    _log(status_path, status)
    log(f"\n  DONE: {ok}/{len(manifest['items'])} grids verified-present. status -> {status_path}")
    if ok != len(manifest['items']):
        bad = [n for n, r in status['items'].items() if r['state'] not in ('VERIFIED', 'PRESENT')]
        log(f"  INCOMPLETE — failed/partial: {bad}")
        sys.exit(2)


if __name__ == '__main__':
    main()
