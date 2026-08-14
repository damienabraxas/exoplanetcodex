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
    # THE root definitions themselves. A register has to bottom out in a real path
    # somewhere, and this is that somewhere -- config/path_register.yaml holds the
    # per-root defaults and constants.py holds SIRIUS_DATA_ROOT for the older RYA-567
    # resolver. Exempting them is not a loophole: it is the single place a machine
    # truth is allowed to be written down.
    # NOTE a real duplication remains -- SIRIUS_DATA_ROOT and the register's `data`
    # root default the same value independently. Folding one into the other needs the
    # register block to move above SIRIUS_DATA_ROOT in constants.py; deferred rather
    # than done blind, because the old resolver still has live call sites.
    "config/constants.py",
    "config/path_register.yaml",
    "scripts/audit_path_literals.py",
    "scripts/migrate_grids_to_ext_rya800.py",
    "scripts/verify_grid_migration_rya800.py",
    "scripts/fetch_gerber_grid.py",
}
BASELINE = 0     # RYA-810 COMPLETE. Every convertible literal is gone, so this is now
# a HARD GATE, not a ratchet: any NEW absolute path literal fails CI. The ticket
# originally said "do not demand zero" -- correct at 130 remaining, wrong now. A gate
# that CAN pass should be enforced. If a genuine machine-truth path is ever needed,
# add it to EXEMPT with a reason rather than raising this number.


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
