#!/usr/bin/env python3
"""RYA-1170 — every DOI in an intake artifact must agree with the SSOT, or say it isn't in it.

RYA-355's single-source rule, enforced. `data/refs/bibliography.csv` is the one
authoritative reference list (RYA-854). An intake that COPIES a DOI out of it has made a
second copy that can drift -- and one did: the RYA-1136 CNO intake carried
10.1051/0004-6361/201936179 for `Amarsi2019_Table1`, which Crossref resolves to Curran &
Moss, "Quasi-stellar object redshift estimates", A&A 629. A different paper, cited as the
source of the entire C I / O I atomic census, while the SSOT held the right DOI all along.

TWO CHECKS, and the second is the RYA-1170 amendment:

  1. DOI AGREEMENT.  A row naming an `ssot_key` must carry exactly that key's DOI.
  2. HELD-DOCUMENT LINK.  If the SSOT declares a `local_file` for the key, the row must
     point at that same path and carry a non-empty sha256. An empty checksum on a key the
     SSOT records as HELD is a broken link, not a missing source -- which is what left the
     AGSS21 reconciliation with no acquired referent.

⚠️ WHAT THIS DELIBERATELY DOES NOT DO. It does not verify the checksum's VALUE unless the
reference library is present. `local_file` is relative to the Codex root, one level above
the repo, and that library is a working directory rather than a repo artifact -- it
differs per machine and is legitimately absent on CI. `generate_sources_page.audit_library`
already says so and treats absence as non-fatal; failing here would only teach people to
disable the check.

Usage:
    python3 scripts/check_intake_ssot_rya1170.py           # report
    python3 scripts/check_intake_ssot_rya1170.py --check   # CI: exit 1 on divergence
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SSOT = ROOT / "data" / "refs" / "bibliography.csv"
LIBRARY_ROOT = Path(
    os.environ.get("CODEX_LIBRARY_ROOT", Path.home() / "Documents" / "Exoplanet Codex"))

#: The intake bibliographies in scope. RYA-1129/1136/1160 per the ticket, plus the two
#: siblings built on the same pattern -- the defect is the pattern, not the ticket.
INTAKE_BIBLIOGRAPHIES = (
    "data/audit/rya1136_cno_intake/source_bibliography.csv",
    "data/audit/rya1132_al_intake/source_bibliography.csv",
)

#: Directories swept for any loose DOI string, so a copy in a file nobody thought of is
#: still seen. A DOI here that the SSOT does not hold is REPORTED, never failed: most of
#: these sources are simply not in bibliography.csv yet.
SWEEP_DIRS = (
    "data/audit/rya1129_atomic_intake",
    "data/audit/rya1136_cno_intake",
    "data/audit/rya1160_cno_nist_gf",
)

DOI_RE = re.compile(r"10\.\d{4,9}/[^\s,\"'\]\}]+")


def ssot_rows() -> dict[str, dict]:
    with SSOT.open(newline="") as fh:
        return {r["key"]: r for r in csv.DictReader(fh)}


def divergences() -> list[str]:
    """Hard failures: a row that names an SSOT key and disagrees with it."""
    ssot = ssot_rows()
    bad: list[str] = []
    for rel in INTAKE_BIBLIOGRAPHIES:
        path = ROOT / rel
        if not path.exists():
            continue
        with path.open(newline="") as fh:
            rows = list(csv.DictReader(fh))
        for row in rows:
            key = (row.get("ssot_key") or "").strip()
            if not key:
                continue
            if key not in ssot:
                bad.append(f"{rel}: {row['source_id']} names ssot_key={key!r}, "
                           f"which bibliography.csv does not hold")
                continue
            want = ssot[key]["doi"].strip()
            got = (row.get("doi") or "").strip()
            if want != got:
                bad.append(f"{rel}: {row['source_id']} DOI {got!r} != SSOT[{key}] {want!r}")
            local = ssot[key]["local_file"].strip()
            if local:
                if (row.get("asset") or "").strip() != local:
                    bad.append(f"{rel}: {row['source_id']} does not point at the held "
                               f"document {local!r} (asset={row.get('asset')!r})")
                if not (row.get("sha256") or "").strip():
                    bad.append(f"{rel}: {row['source_id']} has an EMPTY sha256 while "
                               f"bibliography.csv records {key} as HELD — a broken link, "
                               f"not a missing source")
    return bad


def unlinked() -> list[tuple[str, str, str]]:
    """(file, source_id, doi) for rows carrying a DOI with no SSOT key. Reported only."""
    out = []
    for rel in INTAKE_BIBLIOGRAPHIES:
        path = ROOT / rel
        if not path.exists():
            continue
        with path.open(newline="") as fh:
            for row in csv.DictReader(fh):
                if not (row.get("ssot_key") or "").strip() and (row.get("doi") or "").strip():
                    out.append((rel, row["source_id"], row["doi"].strip()))
    return out


def loose_dois() -> dict[str, set[str]]:
    """Every DOI string in the swept dirs -> the files carrying it."""
    found: dict[str, set[str]] = {}
    for d in SWEEP_DIRS:
        for p in (ROOT / d).rglob("*"):
            if not p.is_file():
                continue
            try:
                text = p.read_text()
            except (UnicodeDecodeError, OSError):
                continue
            for doi in set(DOI_RE.findall(text)):
                found.setdefault(doi.rstrip(".").lower(), set()).add(
                    str(p.relative_to(ROOT)))
    return found


def checksum_status(key: str, recorded: str) -> str:
    rel = ssot_rows()[key]["local_file"].strip()
    if not rel:
        return "no local_file declared"
    path = LIBRARY_ROOT / rel
    if not path.is_file():
        return f"UNVERIFIABLE-HERE (reference library absent at {LIBRARY_ROOT})"
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    return "VERIFIED" if actual == recorded else f"MISMATCH (on disk {actual})"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true", help="exit 1 on any divergence")
    a = ap.parse_args()

    ssot = ssot_rows()
    bad = divergences()
    print(f"=== RYA-1170 intake/SSOT agreement — {len(bad)} divergence(s) ===")
    for b in bad:
        print(f"  🔴 {b}")
    if not bad:
        print("  every linked row matches data/refs/bibliography.csv")

    link = [r for r in ssot.values() if r["local_file"].strip()]
    print(f"\n  SSOT: {len(ssot)} keys, {len(link)} declaring a held document")
    for rel in INTAKE_BIBLIOGRAPHIES:
        path = ROOT / rel
        if not path.exists():
            continue
        with path.open(newline="") as fh:
            for row in csv.DictReader(fh):
                key = (row.get("ssot_key") or "").strip()
                if key and ssot.get(key, {}).get("local_file", "").strip():
                    print(f"    {row['source_id']}: checksum "
                          f"{checksum_status(key, (row.get('sha256') or '').strip())}")

    un = unlinked()
    print(f"\n  {len(un)} row(s) carry a DOI that bibliography.csv does not hold "
          f"(reported, not failed):")
    for rel, sid, doi in un:
        print(f"    {sid:20s} {doi}")

    loose = loose_dois()
    known = {r["doi"].strip().lower() for r in ssot.values() if r["doi"].strip()}
    in_ssot = set(loose) & known
    print(f"\n  loose-DOI sweep over {len(SWEEP_DIRS)} intake dir(s): "
          f"{len(loose)} distinct DOI(s), {len(in_ssot)} of them in the SSOT")

    return 1 if (a.check and bad) else 0


if __name__ == "__main__":
    sys.exit(main())
