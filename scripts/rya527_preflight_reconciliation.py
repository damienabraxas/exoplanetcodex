#!/usr/bin/env python3
"""
scripts/rya527_preflight_reconciliation.py
==========================================
RYA-527 MANDATORY PREFLIGHT (the anti-repetition gate, 2026-07-16).

Before ANY verdict/two-engine run: reconcile the COMMITTED Engine-B provisioning
truth (data/nlte_grids/gerber_ts/*.prov.json, md5-pinned) against the physical
grids on Sirius. LOUD-FAIL, NEVER a silent scope reduction: if the live view
disagrees with committed state, RAISE and STOP naming element + expected md5 +
path. A missing/mismatched grid is a re-fetch from the pinned provenance (or an
RYA-547 storage-eviction to fix) — never "report a lower provisioned count", never
fall back to an overlay.

Verification chain (the grids are stored UNZIPPED on Sirius; the pinned md5 is the
.bin.zip download hash, the zip is deleted post-extract):
  1. model_atom + grid_1d_aux  : recompute md5 on Sirius, MUST equal the prov pin
     (these are stored as-pinned — a direct pin-level md5 check).
  2. grid_1d_bin(.zip)         : prov zip_md5 MUST equal the Sirius RYA-540 cache
     index zip_md5 (the download-time hash the extract was verified against), AND
     the unzipped bin MUST be physically present (symlinks resolved) with a size
     matching the cache index bin_bytes.

Green preflight => the engines/grids are real and provisioned; proceed to the run.
"""
import glob
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PROV_DIR = ROOT / 'data' / 'nlte_grids' / 'gerber_ts'
SIRIUS = 'sirius'
SIRIUS_GRID_DIR = '/mnt/codex-data/grids/nlte/gerber_ts'


def _load_pins():
    pins = {}
    for f in sorted(glob.glob(str(PROV_DIR / '*.prov.json'))):
        d = json.loads(Path(f).read_text())
        el = d.get('element') or os.path.basename(f).split('_')[0]
        files = d.get('files', {})
        atom = files.get('model_atom', {})
        aux = files.get('grid_1d_aux', {})
        binf = files.get('grid_1d_bin') or files.get('grid_1d_bin_zip', {})
        pins[el] = {'atom': (atom.get('name'), atom.get('md5')),
                    'aux': (aux.get('name'), aux.get('md5')),
                    'bin_zip_name': binf.get('name'), 'bin_zip_md5': binf.get('md5'),
                    'bin_zip_bytes': binf.get('bytes')}
    return pins


def _sirius(cmd):
    r = subprocess.run(['ssh', '-o', 'BatchMode=yes', SIRIUS, cmd],
                       capture_output=True, text=True)
    return r.stdout, r.returncode


def _sirius_state():
    """One round trip: atom/aux md5, bin sizes (symlinks resolved), cache index."""
    cmd = (f"cd {SIRIUS_GRID_DIR} && "
           "md5sum atom.* auxData_* 2>/dev/null && "
           "echo '---BINS---' && "
           "for b in NLTEgrid*_*.bin; do "
           "sz=$(stat -Lc '%s' \"$b\" 2>/dev/null || echo MISSING); "
           "echo \"$b $sz\"; done && "
           "echo '---INDEX---' && cat _cache_index.json")
    out, rc = _sirius(cmd)
    if rc != 0:
        raise SystemExit(f"PREFLIGHT LOUD-FAIL: cannot reach Sirius or grid dir ({rc}):\n{out}")
    md5s, bins, index_txt = {}, {}, ''
    section = 'md5'
    for ln in out.splitlines():
        if ln.strip() == '---BINS---':
            section = 'bins'; continue
        if ln.strip() == '---INDEX---':
            section = 'index'; continue
        if section == 'md5' and ln.strip():
            h, _, name = ln.partition('  ')
            md5s[name.strip()] = h.strip()
        elif section == 'bins' and ln.strip():
            name, _, sz = ln.rpartition(' ')
            bins[name.strip()] = sz.strip()
        elif section == 'index':
            index_txt += ln + '\n'
    index = json.loads(index_txt) if index_txt.strip() else {}
    return md5s, bins, index


def main():
    pins = _load_pins()
    print(f"[preflight] {len(pins)} committed Engine-B grids (gerber_ts prov pins): "
          f"{', '.join(sorted(pins))}")
    md5s, bins, index = _sirius_state()

    failures, rows = [], []
    for el in sorted(pins):
        p = pins[el]
        el_fail = []
        # 1. atom + aux: direct md5 must equal the pin
        for kind in ('atom', 'aux'):
            name, want = p[kind]
            got = md5s.get(name)
            if got is None:
                el_fail.append(f"{kind} {name} MISSING on Sirius ({SIRIUS_GRID_DIR}/{name})")
            elif got != want:
                el_fail.append(f"{kind} {name} md5 {got} != pinned {want}")
        # 2. bin: prov zip_md5 == cache-index zip_md5, and bin present + size match
        idx = index.get(el, {})
        if not idx:
            el_fail.append(f"no RYA-540 cache-index entry (expected zip_md5 {p['bin_zip_md5']})")
        else:
            if idx.get('zip_md5') != p['bin_zip_md5']:
                el_fail.append(f"bin zip_md5 {idx.get('zip_md5')} != pinned {p['bin_zip_md5']} "
                               f"({p['bin_zip_name']})")
            bin_name = idx.get('bin_name')
            got_sz = bins.get(bin_name)
            if got_sz in (None, 'MISSING'):
                el_fail.append(f"bin {bin_name} MISSING/unresolved on Sirius "
                               f"({SIRIUS_GRID_DIR}/{bin_name})")
            elif str(got_sz) != str(idx.get('bin_bytes')):
                el_fail.append(f"bin {bin_name} size {got_sz} != index {idx.get('bin_bytes')}")
        status = 'OK' if not el_fail else 'FAIL'
        rows.append((el, status, p['atom'][0], idx.get('bin_name', '?')))
        for m in el_fail:
            failures.append(f"{el}: {m}")

    print("\n  element  status  atom            bin")
    for el, st, atom, binn in rows:
        print(f"  {el:>4s}    {st:<5s}  {atom:<14s}  {binn}")

    if failures:
        raise SystemExit(
            "\nPREFLIGHT LOUD-FAIL — committed provisioning truth disagrees with Sirius. "
            "Do NOT proceed, do NOT reduce the count, do NOT overlay. Re-fetch from the "
            "pinned provenance (or fix RYA-547 storage eviction), then re-run:\n  - "
            + "\n  - ".join(failures))
    print(f"\n  PREFLIGHT GREEN: all {len(pins)} Engine-B grids present + md5-reconciled on "
          f"Sirius. Ti ships atom.ti503b (Engine-B xfail RYA-548; production Ti NLTE is "
          f"Engine-A Mallinson-2024 — not a blocker). Safe to run the two-engine verdict.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
