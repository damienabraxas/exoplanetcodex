#!/usr/bin/env python3
"""RYA-1110 — stage the gf-provenance DECODERS: Heiter+2021 (GES) and Jofré's footnote.

    python3 scripts/rya1110_stage_heiter2021.py                 # stage + verify
    python3 scripts/rya1110_stage_heiter2021.py --verify-only
    python3 scripts/rya1110_stage_heiter2021.py --source-dir ~/Downloads

Closes the gap the first RYA-1110 pass had to leave open: Jofré Tables 4/5 cite gf by
INTEGER CODE, and the arXiv copy we held (1309.1099v2) failed to typeset the footnote that
decodes them — it printed `References: 102: ????????. 114: ??. …`.

TWO DECODERS, AND THEY ARE NOT INTERCHANGEABLE
----------------------------------------------
1. **Heiter et al. 2021, A&A 645 A106** (VizieR J/A+A/645/A106) — the GES line list itself.
   `geslines.dat` carries a PER-LINE `r_loggf`, and `refs.dat` maps every code to an author
   and a bibcode. This is the finer instrument: it says which source THIS line's gf came
   from.
   🔴 BUT IT DESCRIBES GES v6. Jofré used GES **v3** (2014, in prep.). Where the two
   versions carry a different log gf, Heiter's `r_loggf` is the provenance of a DIFFERENT
   NUMBER, and attaching it to Jofré's value would be a fabricated pedigree — the
   `gf_grades` SCALE-MISMATCH defect. The builder therefore uses this route ONLY where
   Heiter's value equals the published GBS value.
2. **Jofré et al. 2014's own Table 4/5 footnote**, transcribed here from the PUBLISHED A&A
   PDF (`aa22440-13.pdf`), which does typeset it. Coarser — each integer code names a LIST
   of sources rather than one — but it is attached to the right number, so it is the
   fallback wherever route 1 is not valid.

WHAT IS STAGED
--------------
    heiter2021_ges/geslines_Fe_4700_6900.tsv   Fe I/II/III rows over the GBS window, the
                                               columns the join uses. Not the whole 80,612-
                                               row table: the full file is ~29 MB and the
                                               other 65,000 rows are other elements.
    heiter2021_ges/refs.tsv                    all 327 reference codes, unfiltered
    jofre2014_gbs/paper_table45_refcodes.tsv   the published footnote, transcribed verbatim

The window is 4700-6900 A, wider than the GBS span (4787.83-6820.37) on purpose: the
displaced-null control in the builder shifts lines by up to 0.5 A and must see the real
line DENSITY, not an edge.

FIREWALL (RYA-161): this stages DECODERS. No gf value is read from here into any product.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import unicodedata
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

CATALOG = "J/A+A/645/A106"
BIBCODE = "2021A&A...645A.106H"
BIB_KEY = "heiter2021"                       # data/refs/bibliography.csv
DEST = ROOT / "data" / "reference" / "heiter2021_ges"
JOFRE_DEST = ROOT / "data" / "reference" / "jofre2014_gbs"

#: The VizieR FITS renderings Ryan downloaded. OUTSIDE the repo (they are 29 MB and 84 kB
#: of FITS); md5-pinned in SOURCES.json so the staged TSVs can be proved to come from them.
SRC_LINES = "J_A+A_645_A106_geslines.dat.gz.fits.gz"
SRC_REFS = "J_A+A_645_A106_refs.dat.fits.gz"

#: The PUBLISHED Jofré PDF, in the RYA-854 reference library. The arXiv copy will NOT do:
#: this is the whole reason the codes were undecodable in the first pass.
JOFRE_PDF = (Path.home() / "Documents" / "Exoplanet Codex" / "Reference documents"
             / "aa22440-13.pdf")

WINDOW_A = (4700.0, 6900.0)
_COLS = ("Element", "Ion", "Isotope", "lambda", "loggf", "e_loggf", "r_loggf",
         "gfflag", "synflag", "Elow")

#: Jofré et al. 2014, Tables 4 and 5, "References." footnote — TRANSCRIBED AS PUBLISHED
#: from aa22440-13.pdf (RYA-161: verbatim, no substitution, no expansion of an "et al.").
#:
#: 🔴 190 IS ABSENT FROM THE PUBLISHED FOOTNOTE. Table 4's body uses it (Fe I 4985.55) and
#: the footnote does not define it. That is a gap in the paper, not in this transcription —
#: checked against the published PDF, both the table row and the footnote line.
JOFRE_FOOTNOTE = {
    102: ("Table 4", "Bard et al. (1991); Bard & Kock (1994); Blackwell et al. "
                     "(1979a,b, 1982a,b, 1995); O'Brian et al. (1991)"),
    114: ("Table 4", "Bridges & Kornblith (1974); Fuhr & Wiese (2006)"),
    129: ("Table 4", "Garz & Kock (1969); Fuhr et al. (1988)"),
    156: ("Table 4", "May et al. (1974)"),
    167: ("Table 4", "Richter & Wulff (1970); Fuhr et al. (1988)"),
    186: ("Table 4", "Wolnik et al. (1970); Fuhr et al. (1988)"),
    187: ("Table 4", "Wolnik et al. (1971); Fuhr et al. (1988)"),
    158: ("Table 5", "Meléndez & Barbuy (2009)"),
    166: ("Table 5", "Raassen & Uylings (1998)"),
}

#: The first-author surname of each source in a footnote entry, for the cross-check against
#: Heiter's per-line codes. Derived by hand from the strings above and asserted against them
#: (`_check_surnames`) so the two cannot drift apart.
JOFRE_SURNAMES = {
    102: ("bard", "blackwell", "obrian"),
    114: ("bridges", "fuhr"),
    129: ("garz", "fuhr"),
    156: ("may",),
    167: ("richter", "fuhr"),
    186: ("wolnik", "fuhr"),
    187: ("wolnik", "fuhr"),
    158: ("melendez",),
    166: ("raassen",),
}

#: The footnote line as it appears, so a reader can see the transcription is complete.
FOOTNOTE_VERBATIM = (
    "Table 4 — References. 102: Bard et al. (1991); Bard & Kock (1994); Blackwell et al. "
    "(1979a,b, 1982a,b, 1995); O'Brian et al. (1991). 114: Bridges & Kornblith (1974); "
    "Fuhr & Wiese (2006). 129: Garz & Kock (1969); Fuhr et al. (1988). 156: May et al. "
    "(1974). 167: Richter & Wulff (1970); Fuhr et al. (1988). 186: Wolnik et al. (1970); "
    "Fuhr et al. (1988). 187: Wolnik et al. (1971); Fuhr et al. (1988).\n"
    "Table 5 — References. 158: Meléndez & Barbuy (2009). 166: Raassen & Uylings (1998).")


class StageError(RuntimeError):
    """A decoder could not be staged in a form the builder can trust."""


def md5(p: Path) -> str:
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _fold(s: str) -> str:
    """ASCII-fold for surname comparison ONLY. The transcription keeps its accents —
    "Meléndez" is how the paper spells it and how this repo must quote it; the ASCII form
    exists so a surname summary can be checked against it, never to replace it."""
    return "".join(c for c in unicodedata.normalize("NFD", s)
                   if unicodedata.category(c) != "Mn").lower().replace("'", "")


def _check_surnames() -> None:
    """JOFRE_SURNAMES must actually appear in the verbatim strings it summarises."""
    for code, names in JOFRE_SURNAMES.items():
        text = _fold(JOFRE_FOOTNOTE[code][1])
        missing = [n for n in names if n not in text]
        if missing:
            raise StageError(
                f"code {code}: surname(s) {missing} are not in the transcribed footnote "
                f"{JOFRE_FOOTNOTE[code][1]!r}. The summary and the verbatim text have "
                f"drifted; fix the summary, never the transcription.")


def stage_geslines(src: Path) -> int:
    from astropy.io import fits
    with fits.open(src) as h:
        d = h[1].data
    el = np.char.strip(d["Element"].astype(str))
    keep = (el == "Fe") & (d["lambda"] >= WINDOW_A[0]) & (d["lambda"] <= WINDOW_A[1])
    rows = d[keep]
    if not len(rows):
        raise StageError(f"no Fe rows in {WINDOW_A} A — wrong file or wrong column names")
    out = DEST / "geslines_Fe_4700_6900.tsv"
    with out.open("w", newline="\n") as fh:
        fh.write(f"# Heiter et al. 2021, A&A 645 A106 ({BIBCODE}) — VizieR {CATALOG}, "
                 f"geslines.dat.\n")
        fh.write(f"# SUBSET: Element == 'Fe', {WINDOW_A[0]:.0f} <= lambda <= "
                 f"{WINDOW_A[1]:.0f} A; {len(_COLS)} of 24 columns. Rows and values are "
                 f"AS PUBLISHED (RYA-161).\n")
        fh.write(f"# Source: {src.name}  md5 {md5(src)}\n")
        fh.write("\t".join(_COLS) + "\n")
        for r in rows:
            fh.write("\t".join(
                (r[c].strip() if isinstance(r[c], str) else
                 f"{r[c]:.4f}" if c == "lambda" else
                 f"{r[c]:.3f}" if c in ("loggf", "e_loggf", "Elow") else
                 str(int(r[c]))) for c in _COLS) + "\n")
    print(f"  {out.relative_to(ROOT)}  {len(rows)} Fe rows")
    return len(rows)


def stage_refs(src: Path) -> int:
    from astropy.io import fits
    with fits.open(src) as h:
        d = h[1].data
    out = DEST / "refs.tsv"
    with out.open("w", newline="\n") as fh:
        fh.write(f"# Heiter et al. 2021, A&A 645 A106 ({BIBCODE}) — VizieR {CATALOG}, "
                 f"refs.dat. THE DECODER: reference code -> author + bibcode.\n")
        fh.write(f"# Source: {src.name}  md5 {md5(src)}\n")
        fh.write("Ref\tAut\tBibCode\tCom\n")
        for r in d:
            fh.write("\t".join(str(r[c]).strip().replace("\t", " ")
                               for c in ("Ref", "Aut", "BibCode", "Com")) + "\n")
    print(f"  {out.relative_to(ROOT)}  {len(d)} reference codes")
    return len(d)


def stage_jofre_footnote() -> int:
    _check_surnames()
    out = JOFRE_DEST / "paper_table45_refcodes.tsv"
    with out.open("w", newline="\n") as fh:
        fh.write("# Jofré et al. 2014, A&A 564 A133 — the Tables 4/5 'References.' "
                 "footnote, TRANSCRIBED AS PUBLISHED (RYA-161).\n")
        fh.write(f"# Source: {JOFRE_PDF.name} (the PUBLISHED A&A PDF). The arXiv copy "
                 f"1309.1099v2 does NOT typeset this footnote — it prints "
                 f"'References: 102: ????????. 114: ??. ...' — which is why the first "
                 f"RYA-1110 pass could not decode the codes.\n")
        fh.write("# 190 IS ABSENT FROM THE PUBLISHED FOOTNOTE although Table 4's body "
                 "uses it (Fe I 4985.55). A gap in the paper, not in this transcription.\n")
        fh.write("code\ttable\tsources_as_published\tfirst_author_surnames\n")
        for code in sorted(JOFRE_FOOTNOTE):
            tbl, txt = JOFRE_FOOTNOTE[code]
            fh.write(f"{code}\t{tbl}\t{txt}\t{','.join(JOFRE_SURNAMES[code])}\n")
    print(f"  {out.relative_to(ROOT)}  {len(JOFRE_FOOTNOTE)} codes")
    # This file lands in the Jofré holding, whose manifest that script owns. Re-pin through
    # ITS writer rather than growing a second one here — an unpinned file in a pinned
    # directory reads as "verified" without being verified.
    sys.path.insert(0, str(ROOT / "scripts"))
    import rya1110_fetch_jofre2014 as J
    J.write_md5sums()
    return len(JOFRE_FOOTNOTE)


def write_sources(src_dir: Path) -> None:
    doc = {
        "ticket": "RYA-1110",
        "heiter2021": {
            "bibcode": BIBCODE, "bibliography_key": BIB_KEY, "vizier_catalog": CATALOG,
            "doi": "10.1051/0004-6361/201936291",
            "source_files": {n: {"path": str(src_dir / n), "md5": md5(src_dir / n)}
                             for n in (SRC_LINES, SRC_REFS)},
            "note": "VizieR FITS renderings, downloaded by Ryan 2026-08-29. Held outside "
                    "the repo; the staged TSVs are the committed derivative.",
        },
        "jofre2014_published_pdf": {
            "path": str(JOFRE_PDF),
            "md5": md5(JOFRE_PDF) if JOFRE_PDF.exists() else None,
            "note": "The PUBLISHED A&A PDF. Carries the Table 4/5 reference footnote that "
                    "arXiv:1309.1099v2 fails to typeset.",
        },
        "ges_refs_bib": {
            "path": str(Path.home() / "Documents" / "Exoplanet Codex"
                        / "Reference documents" / "ges_refs.bib"),
            "note": "Heiter+2021's BibTeX companion. NOT staged: refs.tsv already carries "
                    "the code -> author + bibcode mapping this ticket needs, and the .bib "
                    "would be a second declaration of the same fact.",
        },
    }
    p = DEST / "SOURCES.json"
    p.write_text(json.dumps(doc, indent=2, ensure_ascii=False) + "\n")
    print(f"  {p.relative_to(ROOT)}")


def write_md5sums() -> None:
    for base in (DEST,):
        out = []
        for p in sorted(base.rglob("*")):
            if p.is_file() and p.name != "MD5SUMS.txt":
                out.append(f"{md5(p)}  {p.relative_to(base)}")
        with (base / "MD5SUMS.txt").open("w", newline="\n") as fh:
            fh.write("\n".join(out) + "\n")
        print(f"  {(base / 'MD5SUMS.txt').relative_to(ROOT)}  {len(out)} files")


def verify() -> int:
    fails = []
    for p in (DEST / "geslines_Fe_4700_6900.tsv", DEST / "refs.tsv",
              DEST / "SOURCES.json", JOFRE_DEST / "paper_table45_refcodes.tsv"):
        if not p.exists() or p.stat().st_size == 0:
            fails.append(f"missing or empty: {p.relative_to(ROOT)}")
    # EVERY file in a pinned holding must be pinned. A file present but absent from
    # MD5SUMS is the worst state: the directory verifies and the file is unchecked.
    for base in (DEST, JOFRE_DEST):
        m = base / "MD5SUMS.txt"
        if not m.exists():
            continue
        pinned = {ln.partition("  ")[2] for ln in m.read_text().splitlines() if ln.strip()}
        have = {str(p.relative_to(base)) for p in base.rglob("*")
                if p.is_file() and p.name != "MD5SUMS.txt"}
        for rel in sorted(have - pinned):
            fails.append(f"present but NOT PINNED in {base.name}/MD5SUMS.txt: {rel}")
    sums = DEST / "MD5SUMS.txt"
    if not sums.exists():
        fails.append("missing: heiter2021_ges/MD5SUMS.txt")
    else:
        for line in sums.read_text().splitlines():
            want, _, rel = line.partition("  ")
            p = DEST / rel
            if not p.exists():
                fails.append(f"pinned file gone: {rel}")
            elif md5(p) != want:
                fails.append(f"md5 CHANGED: {rel}")
    try:
        _check_surnames()
    except StageError as e:
        fails.append(str(e))
    for f in fails:
        print(f"  FAIL {f}")
    if not fails:
        print(f"  OK — decoders intact under {DEST.relative_to(ROOT)} "
              f"and {JOFRE_DEST.relative_to(ROOT)}")
    return 1 if fails else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--source-dir", type=Path, default=Path.home() / "Downloads",
                    help="where the VizieR FITS renderings live (default ~/Downloads)")
    ap.add_argument("--verify-only", action="store_true")
    a = ap.parse_args()
    if not a.verify_only:
        src = a.source_dir.expanduser()
        for n in (SRC_LINES, SRC_REFS):
            if not (src / n).exists():
                raise SystemExit(
                    f"{src / n} not found. Both VizieR renderings of {CATALOG} are needed:\n"
                    f"  {SRC_LINES}  (geslines.dat — the per-line atomic data)\n"
                    f"  {SRC_REFS}   (refs.dat — the code -> source decoder)\n"
                    f"Pass --source-dir if they are elsewhere.")
        DEST.mkdir(parents=True, exist_ok=True)
        print(f"Heiter+2021 {CATALOG}:")
        stage_geslines(src / SRC_LINES)
        stage_refs(src / SRC_REFS)
        write_sources(src)
        print("Jofré+2014 Table 4/5 footnote (from the PUBLISHED PDF):")
        stage_jofre_footnote()
        write_md5sums()
    print("Verify:")
    return verify()


if __name__ == "__main__":
    raise SystemExit(main())
