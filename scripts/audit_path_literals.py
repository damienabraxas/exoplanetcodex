"""
scripts/audit_path_literals.py
=====================================
RYA-800 — report hardcoded absolute path literals that should be register keys.

    python3 scripts/audit_path_literals.py [--check]

WHY. The repo carried 113 absolute path literals against 10 resolver calls, and the
RYA-800 grid migration had to audit every one by hand to learn whether it still
resolved. Literals are also a LEAK: they bake mount points, drive names and
usernames into a history that outlives the machine.

`--check` exits non-zero if the literal count RISES above the recorded baseline, so
the number can only go down. It deliberately does NOT demand zero: migrating 113
call sites is incremental work, and a gate that can never pass gets disabled.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCAN = ("pipeline", "scripts", "tests", "config")
# machine-specific prefixes that must not be literals in committed code
PAT = re.compile(r"""['"](/mnt/|/srv/|/home/|/Users/)[^'"]*['"]""")
# the register itself, the migration tooling that must name real devices, and the
# resolver's own docstrings are legitimately allowed to contain them
EXEMPT = {
    "config/path_register.yaml",
    "scripts/audit_path_literals.py",
    "scripts/migrate_grids_to_ext_rya800.py",
    "scripts/verify_grid_migration_rya800.py",
    "scripts/fetch_gerber_grid.py",
}
BASELINE = 119   # RYA-810 batch 1: 135 -> 119 (16 retired, incl. ALL 15 that
# carried a username). Ratchet DOWN only.
# NB two earlier numbers appear in the RYA-810 ticket's history and are both wrong
# for THIS branch: an ad-hoc grep said 113 (it matched only /mnt/codex-data and
# /srv/codex, missing the /home/ and /Users/ literals that actually leak a
# username), and this scanner said 130 on the older RYA-800 branch point. main has
# since advanced and gained literals. Always re-measure on the branch you are on --
# the count grows on its own, which is exactly why the ratchet exists.
# /srv/codex. This scanner also catches /home/ and /Users/ literals, which are
# precisely the ones that leak a username into git. 130 is the honest number.


def scan():
    hits = []
    for d in SCAN:
        for p in sorted((ROOT / d).rglob("*.py")):
            rel = str(p.relative_to(ROOT))
            if rel in EXEMPT:
                continue
            for i, line in enumerate(p.read_text(errors="replace").splitlines(), 1):
                st = line.strip()
                if st.startswith("#"):
                    continue
                if PAT.search(line):
                    hits.append((rel, i, st[:100]))
    return hits


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    hits = scan()
    by_file: dict[str, int] = {}
    for rel, _i, _s in hits:
        by_file[rel] = by_file.get(rel, 0) + 1
    print(f"hardcoded absolute path literals: {len(hits)}  (baseline {BASELINE})\n")
    for rel, n in sorted(by_file.items(), key=lambda x: -x[1]):
        print(f"  {n:4d}  {rel}")
    print("\nReplace with config.constants.codex_path('<key>') / require_codex_path(...)")
    print("Keys live in config/path_register.yaml")
    if a.check and len(hits) > BASELINE:
        sys.exit(f"\nFAIL: literals rose to {len(hits)} (baseline {BASELINE}). "
                 f"New code must use the register.")
