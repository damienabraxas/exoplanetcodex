#!/usr/bin/env python3
"""Extract Lodders, Bergemann & Palme (2025) Table 6 into a citable CSV (RYA-935).

The tracker needed a second literature comparator. The paper is held
(`Reference documents/AbundancesLodders.pdf`, bibliography key `lodders2025`) but
its NUMBERS were not: nothing in the repository carried them. Typing them from
memory would have broken the rule the whole project runs on -- a value cites its
source -- so they are extracted here, reproducibly, and the extractor is
committed beside the table.

Table 6 is "Recommended Atomic Solar System Abundances: Present and Proto-Solar".
Its final four numeric columns per row are
    log(E/H)+12 present, sigma, log(E/H)+12 proto-solar, sigma
and the PRESENT column is the one a solar photospheric measurement compares
against. Proto-solar is carried too, and labelled, because the two differ by
~0.09 dex and quoting the wrong one would look like a real discrepancy.

The parse is checked against values fixed by definition or long-settled
(H = 12.000 exactly) and fails loudly rather than emitting a table that silently
drifted when the PDF text layer changes.
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

#: Anchors the parse. H is 12.000 BY DEFINITION; the rest are values this project
#: already holds independently, so a mis-parsed column cannot pass unnoticed.
SANITY = {"H": 12.000, "O": 8.76, "Mg": 7.57, "Si": 7.55, "Fe": 7.49}
SANITY_TOL = 0.02

ROW = re.compile(r"^\s*(\d{1,3})\s+([A-Z][a-z]?)\s+(.*)$")
NUM = re.compile(r"-?\d+\.?\d*(?:[Ee][+-]?\d+)?")


def extract(pdf: Path, pages: range) -> list[dict]:
    from pypdf import PdfReader
    reader = PdfReader(str(pdf))
    out: list[dict] = []
    for index in pages:
        if index >= len(reader.pages):
            continue
        for line in (reader.pages[index].extract_text() or "").splitlines():
            m = ROW.match(line)
            if not m:
                continue
            nums = NUM.findall(m.group(3))
            if len(nums) < 4:
                continue
            try:
                present, sig_p, proto, sig_pr = (float(x) for x in nums[-4:])
            except ValueError:
                continue
            out.append({"Z": int(m.group(1)), "element": m.group(2),
                        "A_present": present, "sigma_present": sig_p,
                        "A_protosolar": proto, "sigma_protosolar": sig_pr})
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("pdf", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--pages", default="74-76", help="0-based page range holding Table 6")
    args = ap.parse_args()
    lo, hi = (int(x) for x in args.pages.split("-"))
    rows = extract(args.pdf, range(lo, hi + 1))
    if len(rows) < 60:
        raise SystemExit(f"only {len(rows)} element rows parsed; Table 6 has ~83. "
                         f"The PDF text layer likely changed -- fix the parse rather "
                         f"than shipping a partial table.")
    found = {r["element"]: r["A_present"] for r in rows}
    for element, expect in SANITY.items():
        got = found.get(element)
        if got is None or abs(got - expect) > SANITY_TOL:
            raise SystemExit(f"sanity check failed: {element} parsed as {got}, expected "
                             f"~{expect}. Refusing to emit a table whose columns may have "
                             f"shifted.")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        fh.write("# Lodders, K., Bergemann, M. & Palme, H. 2025, Space Science Reviews 221\n")
        fh.write("# Table 6 -- Recommended Atomic Solar System Abundances\n")
        fh.write("# A_present is log(E/H)+12 for the PRESENT Sun; compare photospheric\n")
        fh.write("# measurements against THAT, not against A_protosolar (~0.09 dex higher).\n")
        fh.write("# Extracted by scripts/rya935_extract_lodders.py -- do not hand-edit.\n")
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: r["Z"]))
    print(f"{len(rows)} elements -> {args.out}")
    for e in ("Al", "Fe"):
        r = next((x for x in rows if x["element"] == e), None)
        if r:
            print(f"  {e}: present {r['A_present']:.2f} +/- {r['sigma_present']:.2f}  "
                  f"(proto-solar {r['A_protosolar']:.2f})")


if __name__ == "__main__":
    main()
