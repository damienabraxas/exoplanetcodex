#!/usr/bin/env python3
"""RYA-1172 — group the CNO species by gf authority so a C grade is never read as an O grade.

Every string in `pipeline.cno_gf_pedigree` is quoted from Wiese & Fuhr, NASA LAW 2006
(`Reference documents/20060052476.pdf`). This prints the grouping and writes it as a ledger.
"""
from __future__ import annotations

import collections
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.cno_gf_pedigree import (  # noqa: E402
    PEDIGREE_BY_SPECIES, SOURCE_CITATION, SOURCE_DOCUMENT, UNESTABLISHED,
    is_cno, pedigree_for)

CANONICAL = ROOT / "data" / "linelists" / "canonical_gf.csv"
OUT = ROOT / "data" / "audit" / "rya1129_atomic_intake" / "cno_gf_pedigree_rya1172.csv"


def main() -> int:
    with CANONICAL.open(newline="") as fh:
        rows = [r for r in csv.DictReader(fh) if is_cno(r["species"])]
    counts = collections.Counter(r["species"].strip() for r in rows)

    by_key = collections.defaultdict(list)
    for sp in sorted(counts, key=lambda s: (s.split()[0], len(s))):
        by_key[pedigree_for(sp).key].append(sp)

    print(f"=== RYA-1172 — CNO gf authority is not uniform ===")
    print(f"  source: {SOURCE_DOCUMENT}")
    print(f"          {SOURCE_CITATION}\n")
    ledger = []
    for key in ("NIST_MCHF_PARTIAL_UPDATE_2006", "WFD1996_MONOGRAPH7_OPACITY_PROJECT",
                UNESTABLISHED.key):
        species = by_key.get(key, [])
        if not species:
            continue
        p = pedigree_for(species[0])
        n = sum(counts[s] for s in species)
        print(f"  {key}")
        print(f"    species  : {', '.join(species)}   ({n} rows)")
        print(f"    vintage  : {p.vintage or '(none stated)'}")
        print(f"    method   : {(p.method or '(none stated)')[:88]}")
        # ⚠️ Only say "theory" where a source says so. For the unestablished stages
        # is_laboratory=False is the SAFE DEFAULT (never tier LAB, RYA-1005) and not a
        # claim about the method -- calling it theory would be the same over-claim this
        # ticket exists to remove, one group down.
        note = ("<- critically-evaluated THEORY, not laboratory" if p.method
                else "<- safe default; the method is NOT established, not 'theory'")
        print(f"    lab gf?  : {p.is_laboratory}   {note}")
        print(f"    evidence : {p.evidence[:96]}...\n")
        for s in species:
            ledger.append({"species": s, "rows": counts[s], "gf_authority": key,
                           "compilation": p.compilation, "method": p.method,
                           "vintage": p.vintage or "", "is_laboratory": p.is_laboratory,
                           "evidence": p.evidence, "source_document": SOURCE_DOCUMENT})

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(ledger[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(ledger)
    print(f"  wrote {OUT.relative_to(ROOT)}  ({len(ledger)} species)")

    mixed = {k: v for k, v in by_key.items()}
    print(f"\n  🔴 {len(mixed)} DISTINCT PEDIGREES across CNO. A 'NIST grade' on a C row and")
    print(f"     on an O row are ten years and a different method apart.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
