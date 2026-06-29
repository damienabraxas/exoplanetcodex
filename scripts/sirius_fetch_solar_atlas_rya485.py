#!/usr/bin/env python3
"""
scripts/sirius_fetch_solar_atlas_rya485.py
==========================================
RYA-485 Issue 3 — fetch an independent (non-Kitt-Peak) solar reference atlas onto Sirius's
data drive, md5-verified, no silent partials. The standing rule (RYA-477): reference atlases
live on Sirius; the resolver reads from there. The RYA-481 discipline: a new solar reference
is DECLARED with provenance, never silently substituted for Kitt Peak.

Direct-file fetcher (no tar extraction): per manifest item, download to
<data-root>/<dest_subdir>/, GATE on byte size then md5 against the cited Zenodo checksum;
a mismatch moves the file aside (.BAD) and FAILs the item — never accepted silently. Writes
a per-record provenance manifest + a status JSON. Idempotent (skip md5-verified files).

    python3 sirius_fetch_solar_atlas_rya485.py --manifest solar_atlas_iag_baker2020.json \
        --data-root /mnt/codex-data
    python3 sirius_fetch_solar_atlas_rya485.py ... --status-only
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
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


def _gz_ok(path: Path) -> bool:
    if not str(path).endswith('.gz'):
        return True
    return subprocess.run(['gzip', '-t', str(path)]).returncode == 0


def fetch_one(item, dest_dir: Path, log):
    name = item['name']
    out = dest_dir / name
    exp_md5 = item.get('md5') or None              # None => no upstream checksum published
    rec = {'name': name, 'md5_expected': exp_md5, 'size_expected': item.get('size_bytes')}
    if out.exists() and (exp_md5 is None or _md5(out) == exp_md5) and _gz_ok(out) and out.stat().st_size > 1000:
        rec['state'] = 'PRESENT'; rec['md5_got'] = _md5(out)
        log(f"  PRESENT  {name} ({'md5 ok' if exp_md5 else 'recorded md5'}) — skip")
        return rec
    sz = (f"{item['size_bytes']/1e6:.1f} MB" if item.get('size_bytes') else 'size unknown')
    log(f"  downloading {name} ({sz}) ...")
    t0 = time.time()
    r = subprocess.run(['wget', '-c', '-q', '--tries=5', '--timeout=60', '-O', str(out), item['url']])
    if r.returncode != 0:
        rec['state'] = 'DOWNLOAD_FAILED'; rec['wget_rc'] = r.returncode
        log(f"  DOWNLOAD_FAILED {name} rc={r.returncode}")
        return rec
    got = out.stat().st_size if out.exists() else 0
    rec['size_got'] = got
    if got < 1000:
        rec['state'] = 'TOO_SMALL'
        log(f"  TOO_SMALL {name} ({got} B) — likely an error page, not the atlas")
        return rec
    # size gate only when an expected size is known
    if item.get('size_bytes') and got != item['size_bytes']:
        rec['state'] = 'SIZE_MISMATCH'
        log(f"  SIZE_MISMATCH {name} got {got} want {item['size_bytes']} (partial?)")
        return rec
    # .gz integrity (catches truncation/HTML when no md5 is published)
    if not _gz_ok(out):
        out.rename(out.with_suffix(out.suffix + '.BAD'))
        rec['state'] = 'GZ_CORRUPT'
        log(f"  GZ_CORRUPT {name} — gzip -t failed, moved aside")
        return rec
    got_md5 = _md5(out); rec['md5_got'] = got_md5
    if exp_md5 is not None and got_md5 != exp_md5:
        out.rename(out.with_suffix(out.suffix + '.BAD'))
        rec['state'] = 'MD5_MISMATCH'
        log(f"  MD5_MISMATCH {name} got {got_md5} want {exp_md5} — moved aside")
        return rec
    rec['state'] = 'VERIFIED' if exp_md5 else 'FETCHED_MD5_RECORDED'
    rec['wall_s'] = round(time.time() - t0, 1)
    log(f"  {rec['state']} {name} ({got/1e6:.1f} MB, "
        f"{'md5 ok' if exp_md5 else 'md5 '+got_md5[:8]+'… recorded, no upstream checksum'}) in {rec['wall_s']}s")
    return rec


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--manifest-dir', default=str(Path(__file__).resolve().parent.parent /
                                                  'data' / 'sirius_manifest'))
    ap.add_argument('--data-root', default='/mnt/codex-data')
    ap.add_argument('--status-only', action='store_true')
    args = ap.parse_args(argv)

    man = json.loads((Path(args.manifest_dir) / args.manifest).read_text())
    data_root = Path(args.data_root)
    if not (data_root / '.codex_mounted').exists():
        sys.exit(f"FATAL: {data_root}/.codex_mounted missing — data drive not mounted, refusing")
    dest = data_root / man['dest_subdir']
    dest.mkdir(parents=True, exist_ok=True)
    status_path = dest / f"_fetch_status_{man['manifest_id']}.json"

    def log(m):
        print(m, flush=True)

    log(f"\n  RYA-485 solar-atlas fetch — {man['manifest_id']} -> {dest}")
    log(f"  {man['n_items']} files, {man['total_bytes']/1e6:.0f} MB; {man.get('zenodo_record')}\n")

    if args.status_only:
        for it in man['items']:
            p = dest / it['name']
            log(f"  {'present' if p.exists() else 'absent ':8} {it['name']}")
        return

    OK_STATES = ('VERIFIED', 'PRESENT', 'FETCHED_MD5_RECORDED')
    status = {'manifest_id': man['manifest_id'], 'started': _now(), 'items': {}}
    ok = 0
    for it in man['items']:
        rec = fetch_one(it, dest, log)
        status['items'][it['name']] = rec
        status['updated'] = _now()
        status_path.write_text(json.dumps(status, indent=2))
        if rec['state'] in OK_STATES:
            ok += 1
    status['finished'] = _now(); status['n_ok'] = ok; status['n_total'] = len(man['items'])
    status_path.write_text(json.dumps(status, indent=2))
    # provenance manifest alongside the data (RYA-481 declare-don't-substitute). Record the
    # OBSERVED md5 (= upstream md5 where published, else computed-at-fetch for CDS/NSO).
    (dest / 'PROVENANCE.json').write_text(json.dumps(
        {'ticket': man['ticket'], 'source_cite': man['source_cite'],
         'source_url': man.get('zenodo_record') or man.get('source_url'),
         'zenodo_doi': man.get('zenodo_doi'), 'purpose': man['purpose'], 'fetched': _now(),
         'checksum_provenance': man.get('checksum_provenance', 'upstream md5'),
         'files': [{'name': i['name'], 'md5': status['items'].get(i['name'], {}).get('md5_got'),
                    'md5_source': 'upstream' if i.get('md5') else 'computed-at-fetch',
                    'size': status['items'].get(i['name'], {}).get('size_got')} for i in man['items']],
         'note': 'INDEPENDENT solar reference (non-Kitt-Peak). Declare explicitly in the '
                 'resolver before any 777 swap-test — never silently substitute for Kitt Peak.'},
        indent=2))
    log(f"\n  DONE: {ok}/{len(man['items'])} files verified. status -> {status_path}")
    if ok != len(man['items']):
        sys.exit(2)


if __name__ == '__main__':
    main()
