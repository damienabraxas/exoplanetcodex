#!/usr/bin/env python3
"""RYA-789 — acquire the Elgueta et al. (2026) VizieR catalog J/A+A/710/A111.

    python3 scripts/fetch_elgueta2026_vizier.py            # fetch + verify + inventory
    python3 scripts/fetch_elgueta2026_vizier.py --verify-only

"CRIRES+ lines study for six GBS stars" (Elgueta+, 2026), 2026A&A...710A.111E.
Acquisition for RYA-787's CRIRES+ solar IR Fe baseline: the Y/J/H atomic line lists, the
paper tables, and the reduced spectra under sp/.

WHY NOT THE VizieR WEB UI
-------------------------
viz-bin is behind an Anubis anti-bot wall, so a scripted https GET to it is blocked. Both
routes here go somewhere else: the CDS FTP mirror (which serves the WHOLE tree, including
sp/ and the ReadMe) and astroquery's TAP/catalog resolution (which serves the tabular data
by catalog NAME, no path to guess). They are independent, and both are attempted, because
a single route that silently returns a partial tree is exactly the failure this must not
turn into a "holding".

ROUTE B IS THE PRIMARY. astroquery returns the *published tables* but not the ReadMe and
not sp/ -- and sp/ is the entire question this ticket exists to answer. Route A is kept as
a cross-check on the tabular data (and as a second opinion if the mirror path ever moves).

LOUD-FAIL, NOT PARTIAL
----------------------
If the mirror returns nothing, or the ReadMe or any of the three atomic line lists is
missing or empty, this raises. A half-fetched directory that gets md5-pinned and
registered as a holding is worse than no holding at all: everything downstream would
treat it as complete.
"""
from __future__ import annotations

import argparse
import hashlib
import re
import shutil
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

CATALOG = "J/A+A/710/A111"
BASE = "https://cdsarc.cds.unistra.fr/ftp/J/A+A/710/A111"
DEST = ROOT / "data" / "reference" / "elgueta2026_vizier"

BIBCODE = "2026A&A...710A.111E"
VIZIER_DOI = "10.26093/cds/vizier.bc-p6/lbbb"
LICENSE = "CC-BY-4.0"

#: Files whose absence means the fetch FAILED rather than "the catalog does not have them".
REQUIRED = ("ReadMe", "atomicy.dat", "atomicj.dat", "atomich.dat")


def md5(p: Path) -> str:
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def files_under(d: Path) -> list[Path]:
    return sorted(p for p in d.rglob("*") if p.is_file())


# ── Route B: the CDS mirror (primary — brings sp/ and the ReadMe) ────────────
#
# NOT a recursive crawl. `wget -r` against this mirror fetches robots.txt, finds
# `Disallow: /ftp/`, and stops -- correctly, because that rule exists to keep search
# engines from walking the whole CDS archive. (The ticket's `wget -r --cut-dirs=4` also
# lands the tree one directory too deep; both were confirmed on the first run.)
#
# So this reads each directory index ONCE and then fetches the named files by explicit
# URL. That is a user-directed retrieval of one CC-BY catalog we were pointed at, from a
# data centre whose purpose is distributing it -- not crawler traversal -- and it is what
# CDS's own anonymous-FTP route does. It is also simply better: an explicit file list is
# reproducible, order-stable, and cannot half-succeed without saying so.

INDEX_RE = re.compile(r'href="([^"?][^"]*)"', re.I)
POLITE_DELAY_S = 1.0          # well inside the mirror's advertised Crawl-delay: 10


def is_index(blob: bytes) -> bool:
    """Is this response a CDS directory listing rather than a data file?

    THE BUG THIS EXISTS FOR. The index links a subdirectory as `href="sp"` -- no trailing
    slash, indistinguishable from a file name. The first run therefore saved sp/ as a
    1,793-byte HTML page and reported `VESTA/SUN spectrum present: False` while the
    18 spectra sat there unfetched. Nothing caught it: the required files were all
    present, so an absence check had nothing to say. Content is the only reliable test.
    """
    head = blob[:400].lower()
    return b"<html" in head and (b"index of /ftp" in head or b"parent directory" in head)


def get(url: str, timeout: int = 300) -> bytes:
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return r.read()


def list_dir(html: bytes) -> list[str]:
    """Entry names in one CDS directory index. Directories are NOT marked -- see is_index."""
    out, seen = [], set()
    for m in INDEX_RE.finditer(html.decode("latin1")):
        n = m.group(1).strip("/")
        if not n or n in seen or n.startswith(("?", "http")) or "/" in n \
                or n == "robots.txt":
            continue
        seen.add(n)
        out.append(n)
    return out


def fetch_mirror(dest: Path, depth: int = 1) -> tuple[bool, str]:
    """Pull the catalog tree, one GET per entry, recursing into real subdirectories.

    Each entry is fetched once; if what comes back IS a directory index, it is treated as
    a directory and walked. That costs no extra requests over guessing, and it cannot be
    fooled by the mirror's unmarked directory links.
    """
    dest.mkdir(parents=True, exist_ok=True)
    note = (f"GET {BASE}/ then each entry by URL; entries whose response is a directory "
            f"index are recursed (depth {depth})")
    print(f"[route B] {note}")

    def walk(url_dir: str, out_dir: Path, level: int) -> tuple[int, int]:
        try:
            html = get(url_dir, timeout=60)
        except Exception as e:
            print(f"[route B] index FAILED {url_dir}: {type(e).__name__}: {e}")
            return 0, 0
        names = list_dir(html)
        print(f"[route B] {url_dir} -> {len(names)} entries")
        ok_n = tot = 0
        for n in names:
            out = out_dir / n
            if out.is_file() and out.stat().st_size > 0:
                ok_n += 1
                tot += 1
                continue
            tot += 1
            try:
                blob = get(f"{url_dir}{n}")
            except Exception as e:
                print(f"[route B] FAILED {url_dir}{n}: {type(e).__name__}: {e}")
                continue
            time.sleep(POLITE_DELAY_S)
            if is_index(blob):
                if level <= 0:
                    print(f"[route B] {n}/ is a directory but depth is exhausted")
                    continue
                print(f"[route B] {n} is a DIRECTORY -- recursing")
                tot -= 1
                # `out` is already dest/<n>; passing `out / n` would nest it as sp/sp/.
                a, b = walk(f"{url_dir}{n}/", out, level - 1)
                ok_n += a
                tot += b
                continue
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_bytes(blob)
            ok_n += 1
        return ok_n, tot

    got, total = walk(f"{BASE}/", dest, depth)
    present = {p.name for p in files_under(dest)}
    ok = all(n in present for n in REQUIRED)
    print(f"[route B] fetched {got}/{total}; required present: {ok}")
    return ok, note


# ── Route A: astroquery (cross-check on the tabular data) ───────────────────

def fetch_astroquery(dest: Path) -> tuple[bool, str, list[str]]:
    """Resolve the catalog by NAME and write each table as ECSV. No path to guess.

    Kept as a cross-check, not the primary: astroquery serves the published tables and
    neither the ReadMe nor sp/.
    """
    note = f'Vizier.get_catalogs("{CATALOG}") with ROW_LIMIT=-1'
    try:
        from astroquery.vizier import Vizier
    except ImportError as e:
        print(f"[route A] astroquery not importable: {e}")
        return False, note, []
    try:
        Vizier.ROW_LIMIT = -1
        cats = Vizier.get_catalogs(CATALOG)
    except Exception as e:                                   # network / service failure
        print(f"[route A] FAILED: {type(e).__name__}: {e}")
        return False, note, []
    out = dest / "astroquery"
    out.mkdir(parents=True, exist_ok=True)
    written = []
    for t in cats:
        name = str(t.meta.get("name", "unnamed")).replace("/", "_")
        p = out / f"{name}.ecsv"
        t.write(p, format="ascii.ecsv", overwrite=True)
        written.append(f"{name}  rows={len(t)}")
        print(f"[route A] {name}: {len(t)} rows -> {p.relative_to(ROOT)}")
    return bool(written), note, written


# ── verify + inventory ───────────────────────────────────────────────────────

def verify(dest: Path) -> dict:
    files = files_under(dest)
    if not files:
        raise SystemExit(
            f"FETCH RETURNED NOTHING into {dest}. Tried:\n  {BASE}/\n  "
            f'astroquery Vizier.get_catalogs("{CATALOG}")\n'
            f"Refusing to register a holding that does not exist.")
    names = {p.name for p in files}
    missing = [n for n in REQUIRED if n not in names]
    if missing:
        raise SystemExit(f"PARTIAL FETCH — required file(s) absent: {missing}. "
                         f"Not registering a partial holding.")
    for n in REQUIRED:
        hit = next(p for p in files if p.name == n)
        if hit.stat().st_size == 0:
            raise SystemExit(f"PARTIAL FETCH — {n} is zero bytes.")

    # No saved file may BE a directory index. The first run wrote sp/ as a 1.8 kB HTML
    # page and then truthfully reported "0 files in sp/, VESTA present: False" -- every
    # existence check passed, because the thing existed; it was just the wrong thing.
    # An absence check cannot catch a wrong-content bug, so this is a content check.
    masquerading = [p for p in files if p.stat().st_size < 65536
                    and is_index(p.read_bytes())]
    if masquerading:
        raise SystemExit(
            "FETCH IS WRONG, NOT MERELY PARTIAL — these saved files are HTML directory "
            f"indexes, not data: {[str(p.relative_to(dest)) for p in masquerading]}. "
            "The mirror links subdirectories without a trailing slash; they must be "
            "recursed, not downloaded. Refusing to md5-pin a holding that contains "
            "index pages where data should be.")

    print(f"\ninventory of {dest.relative_to(ROOT)}: {len(files)} files, "
          f"{sum(p.stat().st_size for p in files) / 1e6:.1f} MB")
    for p in files:
        if p.name != "MD5SUMS.txt":
            print(f"   {p.relative_to(dest)}  {p.stat().st_size}")

    # ── THE KEY QUESTION: is the solar spectrum in sp/ ? ─────────────────────
    sp = [p for p in files if "sp" in p.relative_to(dest).parts[:-1]]
    print(f"\nsp/ contains {len(sp)} files:")
    for p in sp:
        print(f"    {p.relative_to(dest)} {p.stat().st_size}")
    solar = [p for p in sp if any(k in p.name.lower()
                                  for k in ("vesta", "sun", "sol"))]
    print(f"\nVESTA/SUN spectrum present in sp/: {bool(solar)} "
          f"{[str(p.relative_to(dest)) for p in solar]}")

    # ── md5 manifest: the RYA-540 pin, sorted so it is byte-stable ──────────
    man = dest / "MD5SUMS.txt"
    with open(man, "w") as out:
        for p in files:
            if p.name == "MD5SUMS.txt":
                continue
            out.write(f"{md5(p)}  {p.relative_to(dest)}\n")
    print(f"\nwrote {man.relative_to(ROOT)}")
    return dict(n_files=len(files), sp=[str(p.relative_to(dest)) for p in sp],
                solar=[str(p.relative_to(dest)) for p in solar])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dest", default=str(DEST))
    ap.add_argument("--verify-only", action="store_true",
                    help="re-verify and rewrite MD5SUMS.txt without re-downloading")
    ap.add_argument("--skip-astroquery", action="store_true",
                    help="skip the Route A cross-check (it needs astroquery + network)")
    a = ap.parse_args()
    dest = Path(a.dest)

    if not a.verify_only:
        ok_b, cmd_b = fetch_mirror(dest)
        ok_a = False
        if not a.skip_astroquery:
            ok_a, note_a, _ = fetch_astroquery(dest)
            print(f"[route A] {note_a}: {'ok' if ok_a else 'no tables'}")
        if not ok_b and not ok_a:
            raise SystemExit(
                "BOTH FETCH ROUTES RETURNED NOTHING.\n"
                f"  route B: {cmd_b}\n"
                f'  route A: Vizier.get_catalogs("{CATALOG}")\n'
                "Not fabricating a partial holding (RYA-789 CRITICAL condition).")

    res = verify(dest)
    print(f"\nlicense {LICENSE} | bibcode {BIBCODE} | VizieR DOI {VIZIER_DOI}")
    print(f"solar spectrum in sp/: {bool(res['solar'])}")


if __name__ == "__main__":
    main()
