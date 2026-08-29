#!/usr/bin/env python3
"""RYA-1110 — acquire the Jofré et al. (2014) GBS metallicity line data.

    python3 scripts/rya1110_fetch_jofre2014.py              # fetch + transcribe + verify
    python3 scripts/rya1110_fetch_jofre2014.py --verify-only

"Gaia FGK benchmark stars: Metallicity" (Jofré+, 2014), A&A 564 A133,
2014A&A...564A.133J, VizieR J/A+A/564/A133. Bibliography key `jofre2014` (RYA-854).

THE HOLDING HAS TWO HALVES AND NEITHER ONE IS SUFFICIENT
--------------------------------------------------------
The published data for this paper is split across two carriers that do not overlap:

  * VizieR carries the SELECTION and the MEASUREMENTS — which lines were used for which
    star (table6.dat), and the per-method equivalent widths and abundances keyed by
    (star, species, wavelength, excitation potential) (ew.dat / abund.dat).
    🔴 IT CARRIES NO log gf AT ALL. Nothing in the VizieR tree gives an oscillator
    strength, so a line list built from VizieR alone has no gf column to fill.

  * The PAPER carries the ATOMIC DATA — Tables 4 and 5 list, for the "golden" lines only,
    wavelength / lower-level energy / log gf / van der Waals / a gf reference code.
    🔴 IT CARRIES NO PER-STAR SELECTION. The table columns are per stellar GROUP, so the
    paper alone cannot say which lines the SUN used.

So the ticket's "wavelength, Elo, log gf, gf provenance" needs BOTH, joined. This script
acquires both halves and pins them; `rya1110_build_gbs_fe_lineset.py` does the join.

WHY THE PDF IS TRANSCRIBED HERE RATHER THAN READ AT BUILD TIME
--------------------------------------------------------------
The PDF lives OUTSIDE the repo (the RYA-854 reference library) and is not something CI
can open. Tables 4 and 5 are transcribed once, to TSV, committed, and md5-pinned against
the PDF they came from — so the build is reproducible on a machine that has no PDF, and
the transcription can still be re-derived and diffed against its source on this one.

FIREWALL (RYA-161): every value here is transcribed AS PUBLISHED. No gf substitution, no
value edits, no filling of a blank from another source. What the paper does not publish
stays empty and is reported as empty.

🔴 WHAT THE PUBLISHED RECORD DOES NOT CONTAIN, RECORDED SO NOBODY RE-SEARCHES FOR IT
------------------------------------------------------------------------------------
1. THE gf REFERENCE CODES CANNOT BE DECODED FROM THE COPY WE HOLD. Tables 4/5 cite gf
   sources by integer code (102, 114, 129, 156, 158, 166, 167, 186, 187, 190) and the
   decoder is the tables' own "References:" footnote. In the arXiv copy in our library
   (1309.1099v2) that footnote failed to typeset and reads literally

       References: 102: ????????. 114: ??. 129: ??. 156: ?. 167: ??. 186: ??. 187: ??.

   The codes are transcribed anyway — they are real published values and they group the
   lines correctly — but they are NOT resolvable to a bibliography key from this source.
   arXiv has only v1 and v2 (checked); www.aanda.org returns HTTP 403 to a scripted GET.
   Resolving them needs the publisher's PDF, fetched by a human. Flagged, not guessed.
2. TABLES 4/5 DO NOT COVER EVERY GOLDEN LINE. table6.dat flags 193 lines golden; the
   paper tables list 183. The 10 without a row (Fe I 4939.7 5083.3 5166.3 5506.8;
   Fe II 5256.9 5316.6 5316.8 6113.3 6149.3 6369.5) have NO published gf here. Verified
   absent by searching the whole extracted PDF text for each wavelength, not merely by
   the table parse failing to find them.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import re
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

CATALOG = "J/A+A/564/A133"
BASE = "https://cdsarc.cds.unistra.fr/ftp/J/A+A/564/A133"
BIBCODE = "2014A&A...564A.133J"
BIB_KEY = "jofre2014"          # data/refs/bibliography.csv
DEST = ROOT / "data" / "reference" / "jofre2014_gbs"
VIZIER = DEST / "vizier"

#: The reference-library PDF Tables 4/5 are transcribed from. OUTSIDE the repo on purpose
#: (RYA-854 owns the library); pinned by md5 so the transcription can be re-derived and
#: proved to have come from this exact file.
PDF = (Path.home() / "Documents" / "Exoplanet Codex" / "Reference documents"
       / "1309.1099v2.pdf")
#: Absence of any one of these means the fetch FAILED, not that the catalog lacks them.
REQUIRED = ("ReadMe", "table6.dat", "stars.dat", "ew.dat", "abund.dat")

#: Row shape of Tables 4 and 5: lambda, Elow, log gf, Waals, reference code. Everything
#: after the code is per-GROUP scatter, which this ticket does not use and does not
#: transcribe (it is not per-line atomic data).
_ROW = re.compile(r'^\s*(\d{4}\.\d{2})\s+(\d+\.\d{4})\s+(-?\d+\.\d{3})\s+(\d+\.\d{3})\s+(\d+)\b')


def md5(p: Path) -> str:
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def fetch_vizier() -> None:
    VIZIER.mkdir(parents=True, exist_ok=True)
    # ew.dat and abund.dat are served gzipped; the rest plain. Store BOTH uncompressed so
    # the committed holding is greppable and byte-diffable.
    for name in ("ReadMe", "table6.dat", "stars.dat"):
        _get(f"{BASE}/{name}", VIZIER / name)
    for name in ("ew.dat", "abund.dat"):
        gz = VIZIER / f"{name}.gz"
        _get(f"{BASE}/{name}.gz", gz)
        raw = gzip.decompress(gz.read_bytes())
        (VIZIER / name).write_bytes(raw)
        gz.unlink()
        print(f"  {(VIZIER / name).relative_to(ROOT)}  {len(raw)} B (ungzipped)")


def _get(url: str, dest: Path) -> None:
    with urllib.request.urlopen(url, timeout=120) as r:
        body = r.read()
    if not body:
        raise RuntimeError(f"{url} returned 0 bytes — a partial holding is worse than none")
    dest.write_bytes(body)
    print(f"  {dest.relative_to(ROOT)}  {len(body)} B")


def _pdf_text() -> list[str]:
    if not PDF.exists():
        raise FileNotFoundError(
            f"{PDF} not found. Tables 4/5 are transcribed from the RYA-854 reference "
            f"library copy of arXiv:1309.1099v2; --verify-only works without it.")
    from pypdf import PdfReader
    r = PdfReader(str(PDF))
    return "\n".join(p.extract_text() or "" for p in r.pages).splitlines()


def transcribe_paper_tables() -> None:
    """Tables 4 (Fe I golden) and 5 (Fe II golden) -> TSV, as published."""
    lines = _pdf_text()
    i4 = next(i for i, l in enumerate(lines) if l.startswith("Table 4. List of"))
    i5 = next(i for i, l in enumerate(lines) if l.startswith("Table 5. List of"))

    def grab(lo, hi):
        out = []
        for l in lines[lo:hi]:
            m = _ROW.match(l)
            if m:
                out.append(tuple(m.groups()))
        return out

    for species, rows, path in (
            ("Fe I", grab(i4, i5), DEST / "paper_table4_fe1_golden.tsv"),
            ("Fe II", grab(i5, len(lines)), DEST / "paper_table5_fe2_golden.tsv")):
        lam = [r[0] for r in rows]
        if len(set(lam)) != len(lam):
            raise RuntimeError(f"{species}: duplicate wavelength in the parse — "
                               f"the table shape changed, do not trust the transcription")
        with path.open("w", newline="\n") as fh:
            fh.write(f"# Jofré et al. 2014, A&A 564 A133 ({BIBCODE}) — "
                     f"{'Table 4' if species == 'Fe I' else 'Table 5'}, "
                     f"golden {species} lines. TRANSCRIBED AS PUBLISHED (RYA-161).\n")
            fh.write(f"# Source PDF: {PDF.name}  md5 {md5(PDF)}\n")
            fh.write("# loggf_ref_code is the paper's integer gf reference code. The "
                     "decoder footnote did NOT typeset in this copy — see the module "
                     "docstring. NOT resolvable to a bibliography key from here.\n")
            fh.write("wavelength_air_A\telow_eV\tlog_gf\tvdw_abo\tloggf_ref_code\n")
            for r in rows:
                fh.write("\t".join(r) + "\n")
        print(f"  {path.relative_to(ROOT)}  {len(rows)} {species} rows")


def write_md5sums() -> None:
    out = []
    for p in sorted(DEST.rglob("*")):
        if p.is_file() and p.name != "MD5SUMS.txt":
            out.append(f"{md5(p)}  {p.relative_to(DEST)}")
    with (DEST / "MD5SUMS.txt").open("w", newline="\n") as fh:
        fh.write("\n".join(out) + "\n")
    print(f"  MD5SUMS.txt  {len(out)} files")


def verify() -> int:
    fails = []
    for name in REQUIRED:
        p = VIZIER / name
        if not p.exists() or p.stat().st_size == 0:
            fails.append(f"missing or empty: {p.relative_to(ROOT)}")
    for p in (DEST / "paper_table4_fe1_golden.tsv", DEST / "paper_table5_fe2_golden.tsv"):
        if not p.exists():
            fails.append(f"missing: {p.relative_to(ROOT)}")
    sums = DEST / "MD5SUMS.txt"
    if not sums.exists():
        fails.append("missing: MD5SUMS.txt")
    else:
        for line in sums.read_text().splitlines():
            want, _, rel = line.partition("  ")
            p = DEST / rel
            if not p.exists():
                fails.append(f"pinned file gone: {rel}")
            elif md5(p) != want:
                fails.append(f"md5 CHANGED: {rel}")
    for f in fails:
        print(f"  FAIL {f}")
    if not fails:
        print(f"  OK — holding intact under {DEST.relative_to(ROOT)}")
    return 1 if fails else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--verify-only", action="store_true")
    a = ap.parse_args()
    if not a.verify_only:
        print(f"VizieR {CATALOG}:")
        fetch_vizier()
        print("Paper Tables 4/5:")
        transcribe_paper_tables()
        write_md5sums()
    print("Verify:")
    return verify()


if __name__ == "__main__":
    raise SystemExit(main())
