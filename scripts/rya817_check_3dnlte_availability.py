#!/usr/bin/env python3
"""Verify the 3D-NLTE availability matrix — RYA-817 Deliverable B.

    python3 scripts/rya817_check_3dnlte_availability.py [--offline]

`data/curation/threednlte_availability.csv` is a HAND-AUTHORED literature survey: it
records, for each canonical species, whether a 3D non-LTE model or grid exists, whether
it covers stars other than the Sun, and what we actually hold. A survey is only as good
as its citations, so this checks two things a reader cannot check by eye:

  COVERAGE   the matrix covers the canonical set EXACTLY — no element invented, none
             dropped. The canonical set is 28 SPECIES, not 28 elements: RYA-109 counts
             Fe I and Fe II separately (Fe I constrains Teff, Fe II constrains log g)
             and RYA-757 adds Zn, so 27 ELEMENTS carry 28 species. Getting that wrong
             is the collision `pipeline/element_freeze.py` warns about, and a survey
             that quietly lists 28 elements would be wrong in a way nobody notices.

  DOI        every DOI resolves at Crossref AND the title Crossref returns matches the
             reference string in the row. RYA-758 found 4 of 22 DOIs wrong in an
             earlier literature pass; a DOI that resolves to the wrong paper is worse
             than no DOI, because it looks checked.

`--offline` runs COVERAGE only and says so, rather than passing silently without the
half of the check that needs the network.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

MATRIX = ROOT / "data" / "curation" / "threednlte_availability.csv"
TRACKER = ROOT / "data" / "audit" / "element_status_tracker.csv"

#: RYA-757. Zn is canonical by intake but has no tracker row yet (never measured here),
#: so it cannot come from the tracker and is named explicitly.
EXTRA_CANONICAL = {"Zn"}
EXPECTED_SPECIES = 28

CROSSREF = "https://api.crossref.org/works/{doi}"


def canonical_elements() -> tuple[set[str], int]:
    """The canonical element set and species count, read from the generated tracker."""
    df = pd.read_csv(TRACKER, comment="#")
    elements = set(df["element"].astype(str)) | EXTRA_CANONICAL
    # species = one row per (element, ion) in the tracker, plus Zn I
    species = len(df.drop_duplicates(["element", "ion"])) + len(EXTRA_CANONICAL)
    return elements, species


def _norm(s: str) -> set[str]:
    return {w for w in re.findall(r"[a-z]+", s.lower()) if len(w) > 3}


def check_doi(doi: str, reference: str) -> tuple[bool, str]:
    """Does this DOI resolve, and does the paper it resolves to match the reference?"""
    try:
        req = urllib.request.Request(
            CROSSREF.format(doi=doi),
            headers={"User-Agent": "exoplanetcodex/RYA-817 (mailto:ryan.damien.schmitt@gmail.com)"})
        with urllib.request.urlopen(req, timeout=30) as r:   # noqa: S310 (fixed https host)
            msg = json.load(r)["message"]
    except urllib.error.HTTPError as e:
        return False, f"Crossref HTTP {e.code} — DOI does not resolve"
    except Exception as e:                                   # noqa: BLE001
        return False, f"Crossref lookup failed: {e}"

    title = " ".join(msg.get("title") or [])
    surnames = {a.get("family", "").lower() for a in msg.get("author", []) or []}
    # A journal reference cites the ISSUE year; Crossref also carries the online-first
    # year, and for a January issue those differ by one (Wang et al. 2021 = MNRAS 500,
    # issue of Jan 2021, published online Nov 2020). Accept ANY of the dates Crossref
    # reports rather than picking one and calling a correct DOI wrong.
    years = set()
    for key in ("published-print", "published-online", "issued", "created"):
        parts = msg.get(key, {}).get("date-parts") or []
        if parts and parts[0]:
            years.add(parts[0][0])

    # WHAT MAKES THIS DISCRIMINATE. The obvious check — first-author surname plus year —
    # does NOT: hand it the Amarsi+2022 DOI against the Wang et al. reference and it
    # passes, because "Amarsi" is a co-author on Wang and the years are one apart. So
    # the check is anchored on VOLUME and PAGE, which are what a journal reference
    # actually pins down and which two different papers essentially never share.
    #
    # The year is deliberately NOT required: Crossref frequently carries only the
    # online-first date while the astronomical citation uses the ISSUE year (Wang et al.
    # is online 2020, MNRAS 500 issue 2 is January 2021, ADS says 2021MNRAS.500.2159W).
    ref_l = reference.lower()
    authors = msg.get("author") or [{}]
    first = authors[0].get("family", "").lower()
    volume = str(msg.get("volume") or "").strip()
    page = str(msg.get("page") or "").split("-")[0].strip()
    problems = []
    if first and first not in ref_l:
        problems.append(f"first author {first!r} not in the reference string")
    if volume and not re.search(rf"\b{re.escape(volume)}\b", reference):
        problems.append(f"volume {volume} not in the reference string")
    if page and not re.search(rf"\b{re.escape(page)}\b", reference, re.IGNORECASE):
        problems.append(f"page/article {page} not in the reference string")
    if not volume and not page:
        problems.append("Crossref returned neither volume nor page — cannot pin the "
                        "reference; verify by hand")
    if problems:
        return False, f"resolves to {title!r} {sorted(years)} — " + "; ".join(problems)
    return True, f"{first} {volume}:{page}"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--offline", action="store_true",
                    help="skip the Crossref DOI verification and say so")
    args = ap.parse_args(argv)

    m = pd.read_csv(MATRIX)
    canon, n_species = canonical_elements()
    failures = []

    print("COVERAGE")
    listed = set(m["element"].astype(str))
    missing, extra = sorted(canon - listed), sorted(listed - canon)
    ok = not missing and not extra
    print(f"  elements: matrix {len(listed)}  canonical {len(canon)}  "
          f"{'PASS' if ok else 'FAIL'}")
    if missing:
        failures.append(f"canonical elements absent from the matrix: {missing}")
    if extra:
        failures.append(f"matrix rows that are not canonical: {extra}")

    # the species count, stated so the 27-elements/28-species distinction is explicit
    matrix_species = sum(len(str(s).split("+")) for s in m["canonical_species"])
    sp_ok = matrix_species == n_species == EXPECTED_SPECIES
    print(f"  species : matrix {matrix_species}  tracker+Zn {n_species}  "
          f"expected {EXPECTED_SPECIES}  {'PASS' if sp_ok else 'FAIL'}")
    if not sp_ok:
        failures.append(f"species count {matrix_species} != {EXPECTED_SPECIES} "
                        f"(Fe carries TWO species — RYA-109)")

    print("\nSTATUS TALLY")
    for col in ("solar_3d_nlte", "offsolar_3d_nlte"):
        print(f"  {col}:")
        for k, v in m[col].value_counts().items():
            print(f"    {k:<32s} {v:2d}   "
                  f"{', '.join(sorted(m.loc[m[col] == k, 'element']))}")

    if args.offline:
        print("\nDOI VERIFICATION: SKIPPED (--offline). The matrix is NOT fully checked.")
    else:
        print("\nDOI VERIFICATION (Crossref)")
        seen = {}
        for col_doi, col_ref in (("solar_doi", "solar_reference"),
                                 ("offsolar_doi", "offsolar_reference")):
            for _, row in m.iterrows():
                doi = str(row[col_doi]).strip()
                if not doi or doi == "nan":
                    continue
                key = (doi, str(row[col_ref])[:60])
                if key in seen:
                    continue
                good, detail = check_doi(doi, str(row[col_ref]))
                seen[key] = good
                if not good:
                    failures.append(f"{row['element']} {col_doi} {doi}: {detail}")
                    print(f"  FAIL {row['element']:<3s} {doi:<34s} {detail}")
                time.sleep(0.2)          # be polite to Crossref
        print(f"  {sum(seen.values())}/{len(seen)} distinct DOIs verified")

    if failures:
        print("\nFAILURES")
        for f in failures:
            print(f"  ! {f}")
        return 1
    print("\nALL CHECKS PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
