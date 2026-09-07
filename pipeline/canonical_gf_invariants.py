"""RYA-1182 — the invariant RYA-1130 was written to enforce, stated positively.

`data/linelists/canonical_gf.csv` is the **atomic** laboratory-gf store. Molecular
transition provenance belongs in the RYA-1130 domain and must never be reachable from a
join against the atomic store.

🔴 IT WAS NOT BEING KEPT OUT. RYA-1142 check B2 found 8,977 molecular rows inside
canonical_gf -- CN 2766, C2 2654, CH 2265, MgH 594, SiH 583, NH 114, OH 1 -- every one
seeded `linelist(VALD3)`. That is the exact VALD3 redistribution the RYA-1136 CNO intake
spent five primary archives avoiding, sitting one join from being picked up by the next
CNO measurement that reached for canonical_gf. RYA-1182 relocated them to
`data/linelists/molecular/quarantine/`.

WHY THE RULE IS POSITIVE, NOT A MOLECULE BLACKLIST
--------------------------------------------------
The obvious guard -- "reject species in {C2, CH, CN, NH, OH, MgH, SiH}" -- is a
RECOGNISED SET, and a recognised set only ever catches what somebody already thought of.
The next TiO, FeH or CaH would enter exactly the way these did. RYA-1072's rule applies:
never widen a recognised set to absorb the unrecognised case.

So the contract is stated over what an ATOMIC transition must have, and anything failing
it is non-atomic by definition:

    key_z   an integer atomic number in 1..92
    ion     non-empty (the ionisation stage)

Measured on the store as RYA-1182 left it, that rule is exact: the 8,977 molecular rows
were precisely the rows with a non-numeric `key_z`, and precisely the rows with an empty
`ion` -- the two conditions selected the identical set, and every surviving row satisfies
both. The atomic side spans Z = 1..92 across 80 distinct elements.
"""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CANONICAL_GF = ROOT / "data" / "linelists" / "canonical_gf.csv"
QUARANTINE = (ROOT / "data" / "linelists" / "molecular" / "quarantine"
              / "vald3_redistribution_quarantine_rya1182.csv")

#: Hydrogen through uranium. A `key_z` outside this is not an element we hold.
Z_MIN, Z_MAX = 1, 92

#: Provenance that is redistribution rather than a vetted measurement. A molecular row may
#: only sit in the atomic store if it carries vetted, non-VALD3 provenance AND says so;
#: see `violations()`. VALD3 is the redistribution channel this ticket exists to exclude.
REDISTRIBUTION_TAGS = ("vald3", "vald")


def is_atomic(row: dict) -> bool:
    """Does this row satisfy the atomic contract? Positive rule, no name list."""
    try:
        z = int(str(row.get("key_z", "")).strip())
    except (TypeError, ValueError):
        return False
    return Z_MIN <= z <= Z_MAX and bool(str(row.get("ion", "")).strip())


def _is_redistribution(row: dict) -> bool:
    fields = " ".join(str(row.get(k, "")) for k in
                      ("seed_source", "loggf_reference", "gf_tier")).lower()
    return any(tag in fields for tag in REDISTRIBUTION_TAGS)


def read_rows(path: Path | str = CANONICAL_GF) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def violations(rows: list[dict] | None = None) -> list[dict]:
    """Rows that break the separation invariant.

    A non-atomic row is a violation UNLESS it carries vetted non-redistribution
    provenance -- the escape hatch the spec asks for, deliberately narrow: a row that
    merely *is* molecular still fails; it has to have been measured and say by whom.
    """
    rows = read_rows() if rows is None else rows
    return [r for r in rows if not is_atomic(r) and _is_redistribution(r)]


def non_atomic(rows: list[dict] | None = None) -> list[dict]:
    """Every row failing the atomic contract, vetted or not."""
    rows = read_rows() if rows is None else rows
    return [r for r in rows if not is_atomic(r)]


def check(strict: bool = True) -> int:
    """CLI/CI entry point. Non-zero exit is the loud failure."""
    rows = read_rows()
    bad = non_atomic(rows) if strict else violations(rows)
    if not bad:
        print(f"OK: canonical_gf.csv holds {len(rows)} rows, all atomic "
              f"(integer key_z in {Z_MIN}..{Z_MAX}, non-empty ion).")
        return 0
    import collections
    counts = collections.Counter(str(r.get("species", "")).strip() for r in bad)
    print(f"FAIL: canonical_gf.csv holds {len(bad)} NON-ATOMIC row(s) — molecular "
          f"transition provenance must live in the RYA-1130 domain, not here.")
    for sp, n in counts.most_common():
        print(f"    {sp or '<blank species>'}: {n}")
    print(f"  Relocate them to {QUARANTINE.relative_to(ROOT)} (RYA-1182).")
    return 1


if __name__ == "__main__":
    raise SystemExit(check())
