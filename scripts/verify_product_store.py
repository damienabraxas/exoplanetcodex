#!/usr/bin/env python3
"""Do the store's provenance claims still hold? — RYA-1039.

Every published row carries the sha256 of the artifact it came from. That is a claim,
and an unchecked claim rots: RYA-1006 silently overwrote the rya984_graded_163 anchor
pools on disk, and the provenance files were BYTE-IDENTICAL either way, so a conditioned
product was undetectable. The published number stayed the same while the thing it came
from changed underneath it.

This re-hashes each source and reports one of four states per row:

  OK       the artifact is present and its bytes match what was published
  DRIFTED  the artifact is present and its bytes DIFFER -- something re-ran and
           overwrote it after publication. The stored VALUE may be stale.
  ABSENT   the artifact is gone. Not an error by itself: band_products is gitignored
           and a Sirius-produced row is expected to be absent on the Mac. It means the
           row can no longer be re-verified HERE, which is worth knowing.
  FOREIGN  the row was produced on another host, so its origin path is not ours to check

DRIFTED is the one that matters. It is the only state that says a published number and
its evidence have come apart, and it is silent in every other view.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import socket
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STORE = ROOT / "data" / "products"


def _sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify(doc: dict, *, here: str) -> list[tuple]:
    out = []
    for section in ("products", "quarantine"):
        for row in doc.get(section, []):
            pr = row.get("provenance") or {}
            host, path, want = pr.get("host"), pr.get("path"), pr.get("sha256")
            # A copied artifact is checked at the COPY, which is the only thing we hold.
            local = pr.get("copied_to") or path
            key = "|".join(str(row.get(f) or "") for f in
                           ("band", "holding", "tier", "route", "treatment"))
            # 🔴 THE ARTIFACT'S EXISTENCE IS THE FACT; the host label is metadata.
            # Keying on `host != hostname` marked 44 locally-produced rows FOREIGN,
            # because they were published under the nickname `mac` while the machine
            # calls itself something else. A verifier that cannot verify what is sitting
            # in front of it is worse than none -- it reports "not mine" for rows it
            # holds. So: if the file is here, CHECK IT, whatever the label says.
            p = Path(local) if local else None
            if p is None or not p.exists():
                state = "FOREIGN" if (host and host != here) else "ABSENT"
                out.append((state, section, key, row.get("A"), local or host)); continue
            got = _sha256(p)
            out.append((("OK" if got == want else "DRIFTED"), section, key,
                        row.get("A"), local))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--star", default="solar")
    ap.add_argument("--element", default=None)
    ap.add_argument("--quiet", action="store_true", help="only show non-OK rows")
    a = ap.parse_args()

    here = socket.gethostname().split(".")[0]
    files = ([STORE / a.star / f"{a.element}.json"] if a.element
             else sorted((STORE / a.star).glob("*.json")))
    rc = 0
    for f in files:
        if not f.exists():
            print(f"no such element product: {f}"); return 2
        doc = json.loads(f.read_text())
        rows = verify(doc, here=here)
        tally = Counter(r[0] for r in rows)
        print(f"{f.name}  v{doc.get('version')}  " +
              "  ".join(f"{k} {v}" for k, v in sorted(tally.items())))
        for state, section, key, A, where in rows:
            if a.quiet and state == "OK":
                continue
            mark = "🔴" if state == "DRIFTED" else "  "
            print(f"  {mark} {state:<8} {section:<10} A={str(A):<8} {key}")
            if state == "DRIFTED":
                print(f"       artifact changed since publication: {where}")
        if tally.get("DRIFTED"):
            rc = 1
    if rc:
        print("\n🔴 DRIFTED rows exist: a published value and its evidence have come "
              "apart. Re-derive and republish with --reason, or confirm the artifact "
              "was rewritten with identical science.")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
