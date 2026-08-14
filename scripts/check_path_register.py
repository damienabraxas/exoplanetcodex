"""
scripts/check_path_register.py
==============================
RYA-810 — verify every entry in `config/path_register.yaml` actually resolves.

    python3 scripts/check_path_register.py [--strict]

Uses the REAL resolver (config.constants.codex_root/codex_path) rather than
reimplementing the expansion. An earlier throwaway checker did reimplement it, and
promptly reported two GOOD entries as BROKEN because it did not know about the
`{repo_parent}` token -- a checker that disagrees with the thing it checks is worse
than no checker.

Entries are skipped, not failed, when their root is not applicable to this machine:
`local_data` is workstation-side and absent on Sirius; the Sirius roots are absent on
the Mac. `--strict` fails on any skip, for use where everything is expected present.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from config.constants import _path_register, codex_root, codex_path  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--strict", action="store_true")
    a = ap.parse_args()
    reg = _path_register()

    print("roots:")
    live = {}
    for name in sorted(reg["roots"]):
        r = codex_root(name)
        live[name] = r.is_dir()
        rc = reg["roots"][name]
        extra = ""
        if rc.get("removable"):
            m = rc.get("mount")
            extra = f"  [removable: mounted={os.path.ismount(m) if m else '?'}]"
        print(f"  {'present' if live[name] else 'ABSENT ':8s} {name:12s} {r}{extra}")

    print("\nentries:")
    ok = bad = skip = 0
    for key in sorted(reg["entries"]):
        e = reg["entries"][key]
        p = codex_path(key)
        if not live[e["root"]]:
            print(f"  skip    {key:32s} (root {e['root']!r} not on this machine)")
            skip += 1
            continue
        good = p.is_dir() if e["kind"] == "dir" else p.is_file()
        print(f"  {'OK    ' if good else 'BROKEN'} {key:32s} -> {p}")
        ok, bad = ok + good, bad + (not good)

    print(f"\n{ok} ok, {bad} broken, {skip} skipped (root absent here)")
    if bad or (a.strict and skip):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
