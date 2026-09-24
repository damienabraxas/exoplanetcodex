#!/usr/bin/env python3
"""
Register Version-pointer guard (RYA-1184 smoke 3, enforcing RYA-690).

CODEX_STATE_REGISTER.md carries a `**Version: vNNN**` header line whose ONLY job
is to name the newest Changelog row. RYA-690 is the defect this prevents: the file
once carried TWO Version header lines and a changelog row with no matching header,
so version numbers stopped being unique and every branch inherited the ambiguity.

Checks, all loud-fail:
  1. exactly ONE Version header line exists
  2. every Changelog row version is unique
  3. the header names the NUMERICALLY HIGHEST changelog row

Rows are not required to appear in sorted order -- the file is append-oriented and
v152 currently sits above v151 and v147. Order is not the invariant; the POINTER is.
"""

from __future__ import annotations

import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGISTER = ROOT / "CODEX_STATE_REGISTER.md"

HEADER_RE = re.compile(r"^\*\*Version:\s*v(?P<n>\d+)\*\*", re.M)
ROW_RE = re.compile(r"^-\s+\*\*v(?P<n>\d+)\*\*", re.M)


def main() -> int:
    if not REGISTER.exists():
        print(f"FAIL: {REGISTER} is missing", file=sys.stderr)
        return 2
    text = REGISTER.read_text(encoding="utf-8")

    headers = [int(m.group("n")) for m in HEADER_RE.finditer(text)]
    rows = [int(m.group("n")) for m in ROW_RE.finditer(text)]

    print("Register Version-pointer guard (RYA-1184 / RYA-690)")
    print(f"  Version header lines : {len(headers)} {headers}")
    print(f"  Changelog rows       : {len(rows)}")

    failed = False

    if len(headers) != 1:
        print(f"FAIL: expected exactly ONE Version header line, found {len(headers)} "
              f"-- this is the RYA-690 defect", file=sys.stderr)
        failed = True

    dupes = [v for v, c in Counter(rows).items() if c > 1]
    if dupes:
        print(f"FAIL: duplicate Changelog versions: {sorted(dupes)} "
              f"-- version numbers must be unique (RYA-690)", file=sys.stderr)
        failed = True

    if not rows:
        print("FAIL: no Changelog rows found -- cannot validate the pointer", file=sys.stderr)
        return 2

    newest = max(rows)
    if headers and headers[0] != newest:
        print(f"FAIL: Version header names v{headers[0]} but the newest Changelog row "
              f"is v{newest} -- the header is a POINTER, not a record (RYA-690)",
              file=sys.stderr)
        failed = True

    if failed:
        return 1
    print(f"OK: Version header names v{newest}, the newest of {len(rows)} unique Changelog rows.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
