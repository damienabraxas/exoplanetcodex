#!/usr/bin/env python3
"""
stage_target_data.py — transfer audited KEEPER spectra for ONE target from the
local data tree to the Sirius data drive, with sha256 verification on both ends.

RYA-420. The reusable replacement for one-off rsync / download-use-delete. It:
  * refuses to run unless the remote data drive is mounted (runs the RYA-419
    mount guard on the remote first);
  * transfers ONLY the keepers named in an audit-derived manifest (relative paths,
    one per line, '#' comments ok) — never a whole folder, never a list from memory;
  * sha256-hashes every present keeper locally, rsyncs by --files-from, then runs
    `sha256sum -c` on the REMOTE — staged bytes must equal source bytes;
  * NEVER deletes any local source (reclaim is a separate, gated step);
  * loudly reports keepers that are listed but absent on disk (e.g. a URL that was
    never downloaded) and never counts them as staged.

All paths are arguments — no hardcoded source/dest/manifest.

Usage:
  python scripts/stage_target_data.py --target procyon \
    --source-root "/abs/.../Procyon" --keeper-manifest procyon_keepers.present \
    --dest-host sirius [--dest-root /mnt/codex-data] [--dry-run]
"""
import argparse
import hashlib
import subprocess
import sys
from pathlib import Path
# Standalone-script bootstrap (RYA-313): repo root on sys.path BEFORE importing
# config/pipeline, so this runs from any cwd. Derived from __file__, never cwd.
import os as _os_boot, sys as _sys_boot
_sys_boot.path.insert(0, _os_boot.path.dirname(_os_boot.path.dirname(
    _os_boot.path.abspath(__file__))))
from config.constants import codex_root  # RYA-810 path register


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run(cmd, **kw):
    print("+ " + " ".join(map(str, cmd)), flush=True)
    return subprocess.run(cmd, check=True, **kw)


def check_remote_mount(dest_host: str, guard: str, dest_root: str):
    print(f"[precondition] verifying data drive mounted on {dest_host} ...", flush=True)
    with open(guard, "rb") as g:
        r = subprocess.run(["ssh", dest_host, f"CODEX_DATA_ROOT={dest_root} bash -s"],
                           stdin=g, capture_output=True, text=True)
    print(r.stdout.strip())
    if r.returncode != 0:
        sys.exit(f"CRITICAL: remote mount guard failed:\n{r.stderr.strip()}\n"
                 f"Aborting — drive not mounted on {dest_host} (see RYA-419).")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True)                      # procyon | tau_boo
    ap.add_argument("--source-root", required=True, type=Path)      # local data root for this target
    ap.add_argument("--keeper-manifest", required=True, type=Path)  # rel paths, one per line, # ok
    ap.add_argument("--dest-host", required=True)                   # e.g. sirius
    ap.add_argument("--dest-root", default=str(codex_root("data")))
    ap.add_argument("--dest-subdir", default=None,                  # default spectra/<target>
                    help="dest path under dest-root (default: spectra/<target>)")
    ap.add_argument("--mount-guard", default="scripts/check_data_mount.sh")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    check_remote_mount(a.dest_host, a.mount_guard, a.dest_root)

    keepers = [ln.strip() for ln in a.keeper_manifest.read_text().splitlines()
               if ln.strip() and not ln.lstrip().startswith("#")]
    present = [k for k in keepers if (a.source_root / k).is_file()]
    missing = [k for k in keepers if k not in present]
    if missing:
        print(f"WARNING: {len(missing)} keeper(s) listed but NOT on disk "
              f"(e.g. URL-not-downloaded):")
        for m in missing:
            print(f"   MISSING: {m}")
        print("   -> NOT transferred. Resolve before calling that arm complete.")
    if not present:
        sys.exit("CRITICAL: no keeper files present on disk — nothing to stage.")

    src_sums = {k: sha256_file(a.source_root / k) for k in present}
    sums_path = Path(f"{a.target}_src.sha256")
    sums_path.write_text("".join(f"{src_sums[k]}  {k}\n" for k in present))
    files_from = Path(f"{a.target}_keepers.present")
    files_from.write_text("".join(k + "\n" for k in present))
    print(f"[checksum] {len(present)} files hashed -> {sums_path}")

    subdir = a.dest_subdir or f"spectra/{a.target}"
    dest = f"{a.dest_host}:{a.dest_root}/{subdir}/"
    # ensure the dest dir exists (skeleton has it, but be robust for new targets)
    if not a.dry_run:
        run(["ssh", a.dest_host, f"mkdir -p '{a.dest_root}/{subdir}'"])
    rsync = ["rsync", "-av", "--checksum", f"--files-from={files_from}",
             f"{a.source_root}/", dest]
    if a.dry_run:
        rsync.insert(1, "--dry-run")
        run(rsync)
        print("[dry-run] no verification; re-run without --dry-run to transfer.")
        return
    run(rsync)

    # verify on the remote by sha256 (paths in the manifest are relative; -c matches them)
    run(["scp", str(sums_path), f"{a.dest_host}:{a.dest_root}/{subdir}/.staged.sha256"])
    r = subprocess.run(["ssh", a.dest_host,
                        f"cd '{a.dest_root}/{subdir}' && sha256sum -c .staged.sha256"],
                       capture_output=True, text=True)
    print(r.stdout)
    print(r.stderr)
    if r.returncode != 0:
        sys.exit("CRITICAL: remote checksum verification FAILED — staged data != source. "
                 "Do NOT reclaim local copies.")
    print(f"[done] {a.target}: {len(present)} staged + verified, {len(missing)} missing. "
          f"Local source left INTACT.")


if __name__ == "__main__":
    main()
