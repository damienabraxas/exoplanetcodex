#!/usr/bin/env python3
"""
scripts/generate_sources_page.py
================================
RYA-854 -- generate the public /sources page from `data/refs/bibliography.csv`.

    python3 scripts/generate_sources_page.py            # write the page
    python3 scripts/generate_sources_page.py --check     # regenerate + diff, exit 1 on drift
    python3 scripts/generate_sources_page.py --site-root ../exoplanetcodex-site

WHY GENERATED, NEVER HAND-AUTHORED (the RYA-775 pattern). A hand-written reference
page rots: it drifts from the line lists, from the tickets, and from the papers it
claims to cite, and nothing anywhere reports that it has. So there is exactly ONE
source of truth -- `data/refs/bibliography.csv` -- and this script is the only thing
allowed to write the page. `--check` makes that enforceable: it regenerates into
memory and diffs against the committed page, so a hand-edit to the HTML is a build
break rather than a silent divergence.

NO CITATION DATA LIVES IN THIS FILE. Authors, years, titles, venues, DOIs, URLs,
roles, tickets and notes all come from the CSV. What lives here is *schema and
presentation*: the category vocabulary with its display labels and ordering, and the
markup. If a row's category is not in that vocabulary the script FAILS -- a new
category has to be a deliberate edit here, not a typo that silently vanishes from the
page.

THE LINK DISCIPLINE (RYA-854 spec). A row whose `url` is blank renders WITHOUT a link
and carries a visible "no verified link" marker plus its `license_note`. That is
deliberate: several sources in the pool are pre-DOI literature with no registered DOI,
and inventing a plausible-looking DOI would be a fabricated provenance record. An
absent link that says why it is absent is honest; a guessed one is not.

THE HASH IN THE HEADER pins the INPUT (the CSV), not the commit. Stamping
`git rev-parse HEAD` into an artifact makes the artifact change on every unrelated
commit and tells a reader nothing about what the page was built from.
"""
from __future__ import annotations

import argparse
import csv
import difflib
import hashlib
import html
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BIB_CSV = ROOT / "data" / "refs" / "bibliography.csv"
PAGE = ROOT / "docs" / "site" / "sources" / "index.html"

#: repo-relative default for a sibling checkout of the website repo
SITE_ROOT_DEFAULT = ROOT.parent / "exoplanetcodex-site"
SITE_PAGE_REL = Path("sources") / "index.html"

REQUIRED_COLS = ("key", "authors", "year", "title", "venue", "doi", "url", "category",
                 "role_in_codex", "tickets", "local_file", "license_note", "verified")

#: The category vocabulary: (csv value, section heading, section blurb).
#: Order here IS the order on the page -- the solar scale first, because every Codex
#: number is differential against it; the project's own documents last.
CATEGORIES: tuple[tuple[str, str, str], ...] = (
    ("solar_reference", "Solar reference scales &amp; atlases",
     "The abundance scales every Codex measurement is referenced against, and the "
     "solar spectra used as the observational anchor."),
    ("atomic_data_gf", "Atomic data &mdash; oscillator strengths",
     "Where each log gf comes from. The distinction that matters is PRIMARY LABORATORY "
     "(branching fractions from a Fourier-transform spectrometer times a measured "
     "radiative lifetime, with no reference to the solar spectrum) versus "
     "semi-empirical or solar-fitted values, which cannot referee a solar abundance "
     "without restating their own input."),
    ("atomic_data_broadening", "Atomic data &mdash; line broadening",
     "Collisional broadening theory behind the van der Waals damping the line lists "
     "carry."),
    ("nlte", "Non-LTE departures",
     "The model atoms and departure-coefficient grids that move an abundance off its "
     "LTE value &mdash; including the papers that define where a grid does NOT reach."),
    ("nlte_synthesis", "Non-LTE synthesis codes",
     "Radiative-transfer codes that solve statistical equilibrium during synthesis."),
    ("model_atmosphere", "Model atmospheres",
     "The 1D and 3D atmospheric structures the synthesis runs on top of."),
    ("benchmark_params", "Benchmark stellar parameters",
     "Independently determined temperatures, gravities, metallicities and "
     "microturbulences &mdash; parameters the Codex adopts rather than fits."),
    ("instrument_line_selection", "Instruments &amp; line selection",
     "External line-selection catalogues and instrument-specific abundance work the "
     "Codex cross-checks its own pools against."),
    ("cno", "Carbon, nitrogen &amp; oxygen",
     "The C/N/O literature behind the C/O ratio &mdash; the flagship derived product "
     "and the hardest measurement in the project."),
    ("motivation", "Why this project exists",
     "The results that make a host-star abundance worth measuring at all."),
    ("method_tool", "Methods, codes &amp; tools",
     "The software and methodological references the pipeline is built on."),
    ("project_doc", "Exoplanet Codex documents",
     "The project&rsquo;s own methodology, glossary and architecture documents."),
)

#: `verified` vocabulary -> (badge text, explanation). A row's verification state is
#: part of its provenance, so it is shown rather than hidden.
VERIFIED_LABELS: dict[str, tuple[str, str]] = {
    "extracted": ("PDF", "Metadata read directly off the paper held in the local "
                         "reference library."),
    "crossref": ("DOI", "DOI record fetched from Crossref; title, authors and venue "
                        "matched against this row."),
    "publisher": ("SITE", "Publisher or database landing page resolved and matched."),
    "unconfirmed": ("UNCONFIRMED", "No DOI or stable URL could be confirmed &mdash; "
                                   "nothing is guessed. See the note."),
}

#: A row whose role_in_codex opens with this string is flagged as a reference standard,
#: consistent with the abundance pages (RYA-851). Driven by the CSV, not a name list.
GOLD_MARKER = "GOLD reference standard"

NAV = (
    ("/systems/", "Systems"),
    ("/glossary/", "Glossary"),
    ("/method/", "Method"),
    ("/sources/", "Sources"),
    ("/mission-log/", "Mission Log"),
    ("/origin/", "Origin"),
    ("/about/", "About"),
    ("/roadmap.html", "Roadmap"),
)


# ─────────────────────────────────────────────────────────────────── loading ──
def load_rows(csv_path: Path) -> list[dict[str, str]]:
    """Read the bibliography and refuse anything the page cannot render honestly."""
    with csv_path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        missing = [c for c in REQUIRED_COLS if c not in (reader.fieldnames or [])]
        if missing:
            sys.exit(f"::error::{csv_path.name} is missing columns: {missing}")
        rows = [{k: (v or "").strip() for k, v in r.items()} for r in reader]

    if not rows:
        sys.exit(f"::error::{csv_path.name} has no rows")

    known = {c for c, _, _ in CATEGORIES}
    problems: list[str] = []
    seen: set[str] = set()
    for r in rows:
        if not r["key"]:
            problems.append("a row has no key")
        elif r["key"] in seen:
            problems.append(f"duplicate key {r['key']!r}")
        seen.add(r["key"])
        if r["category"] not in known:
            problems.append(f"{r['key']}: unknown category {r['category']!r} "
                            f"(add it to CATEGORIES deliberately)")
        if r["verified"].startswith("verify_"):
            # RYA-854's hard requirement: no row may still be awaiting verification.
            problems.append(f"{r['key']}: still pending verification "
                            f"({r['verified']}) -- resolve it against the source "
                            f"before publishing")
        elif r["verified"] not in VERIFIED_LABELS:
            problems.append(f"{r['key']}: unknown verified state {r['verified']!r}")
        if not r["authors"] or not r["year"]:
            problems.append(f"{r['key']}: authors and year are mandatory")
        # A doi.org link must resolve the row's OWN doi -- a copy-paste that leaves a
        # neighbour's DOI in place is exactly the class of error this catches. A
        # non-doi.org url (arXiv abs page, NIST/VALD tool, project PDF) is unconstrained.
        if r["url"].startswith("https://doi.org/") and r["doi"] not in r["url"]:
            problems.append(f"{r['key']}: doi.org url does not carry its own doi "
                            f"({r['doi']!r} vs {r['url']!r})")
        if not r["url"] and not r["license_note"]:
            problems.append(f"{r['key']}: no url AND no note explaining why "
                            f"-- an absent link must say why it is absent")
        # `unconfirmed` claims "no stable public identifier for this source could be
        # found", so a row wearing that badge must not also be showing a link. The
        # implication runs ONE way only, and the converse is deliberately NOT asserted:
        # `verified` records how the METADATA was confirmed, while url/doi record
        # whether a PUBLIC LINK exists, and those are two independent axes. An
        # unpublished project document read straight off the PDF is `extracted` with no
        # url and is not "unconfirmed" -- conflating the axes would force a false badge
        # onto it. (The first version of this check asserted both directions and did
        # exactly that.)
        if r["verified"] == "unconfirmed" and (r["url"] or r["doi"]):
            problems.append(
                f"{r['key']}: verified='unconfirmed' but carries doi={r['doi']!r} / "
                f"url={r['url']!r} -- a row claiming nothing could be found must not "
                f"display a link")
    if problems:
        for p in problems:
            print(f"::error::bibliography: {p}", file=sys.stderr)
        sys.exit(1)
    return rows


def csv_digest(csv_path: Path) -> str:
    return hashlib.sha256(csv_path.read_bytes()).hexdigest()


# ────────────────────────────────────────────────────────────────── rendering ──
def esc(text: str) -> str:
    return html.escape(text, quote=True)


def citation(row: dict[str, str]) -> str:
    """Authors (Year). Title. Venue. -- title omitted when the CSV has none, which is
    itself a provenance statement: see the row's note."""
    bits = [f"{esc(row['authors'])} ({esc(row['year'])})."]
    if row["title"]:
        bits.append(f"<em>{esc(row['title'])}</em>.")
    if row["venue"]:
        bits.append(f"{esc(row['venue'])}.")
    return " ".join(bits)


def render_entry(row: dict[str, str]) -> list[str]:
    out: list[str] = []
    badge, badge_why = VERIFIED_LABELS[row["verified"]]
    gold = row["role_in_codex"].startswith(GOLD_MARKER)

    out.append(f'      <li class="src-entry" id="{esc(row["key"])}">')
    out.append('        <div class="src-head">')
    out.append(f'          <span class="src-cite">{citation(row)}</span>')
    if gold:
        out.append('          <span class="src-gold" title="Reference standard: the '
                   'Codex compares against this value, it does not re-derive it.">'
                   'GOLD</span>')
    out.append(f'          <span class="src-verified src-v-{esc(row["verified"])}" '
               f'title="{badge_why}">{esc(badge)}</span>')
    out.append('        </div>')

    if row["url"]:
        label = row["doi"] if row["doi"] else row["url"]
        out.append(f'        <div class="src-link"><a href="{esc(row["url"])}" '
                   f'target="_blank" rel="noopener">{esc(label)} &#8599;</a></div>')
    else:
        out.append('        <div class="src-link src-nolink">no verified link '
                   '&mdash; see note</div>')

    if row["role_in_codex"]:
        out.append(f'        <p class="src-role">{esc(row["role_in_codex"])}</p>')
    if row["license_note"]:
        out.append(f'        <p class="src-note">{esc(row["license_note"])}</p>')

    meta: list[str] = []
    if row["tickets"]:
        meta.append(f'<span class="src-tickets">{esc(row["tickets"])}</span>')
    if row["local_file"]:
        meta.append('<span class="src-local">local copy: '
                    f'{esc(row["local_file"])}</span>')
    meta.append(f'<span class="src-key">{esc(row["key"])}</span>')
    out.append(f'        <div class="src-meta">{" &middot; ".join(meta)}</div>')
    out.append('      </li>')
    return out


def render_nav() -> list[str]:
    out = ['  <nav class="codex-nav">',
           '    <a href="/" class="nav-logo"><span class="logo-text">EXOPLANET CODEX'
           '</span></a>',
           '    <ul class="nav-links">']
    for href, label in NAV:
        active = ' class="active"' if href == "/sources/" else ""
        out.append(f'      <li><a href="{href}"{active}>{label}</a></li>')
    out.append('      <li><a href="https://github.com/damienabraxas/exoplanetcodex" '
               'target="_blank">GitHub &#8599;</a></li>')
    out.append('    </ul>')
    out.append('    <div class="breadcrumb"><a href="/">Home</a> / '
               '<span>Sources</span></div>')
    out.append('  </nav>')
    return out


STYLE = """    .src-intro { padding:3rem 0 1rem; }
    .src-intro p { font-size:1.02rem; color:var(--text-mid); max-width:760px; line-height:1.85; margin-bottom:1.1rem; }
    .src-intro p strong { color:var(--text-bright); font-weight:400; }
    .src-intro a { color:var(--accent); text-decoration:none; }
    .src-intro a:hover { opacity:0.75; }
    .src-stats { display:flex; flex-wrap:wrap; gap:1px; background:var(--rim); border:1px solid var(--rim); margin:2rem 0 1rem; }
    .src-stat { flex:1 1 140px; background:var(--deep); padding:1rem 1.2rem; }
    .src-stat b { display:block; font-family:var(--mono); font-size:1.5rem; color:var(--accent); font-weight:400; line-height:1.2; }
    .src-stat span { font-family:var(--mono); font-size:0.58rem; letter-spacing:0.18em; text-transform:uppercase; color:var(--text-dim); }
    .src-toc { display:flex; flex-wrap:wrap; gap:0.5rem 1.2rem; margin:0 0 3rem; padding:1.2rem 0 0; border-top:1px solid var(--rim); }
    .src-toc a { font-family:var(--mono); font-size:0.68rem; letter-spacing:0.05em; color:var(--text-dim); text-decoration:none; }
    .src-toc a:hover { color:var(--accent); }
    .src-section { margin:0 0 3.5rem; }
    .src-section h2 { font-family:var(--serif); font-size:clamp(1.4rem,2.2vw,1.9rem); font-weight:300; color:var(--text-bright); line-height:1.25; margin-bottom:0.5rem; padding-top:1.6rem; border-top:1px solid var(--rim); }
    .src-section h2 .src-count { font-family:var(--mono); font-size:0.6rem; letter-spacing:0.15em; color:var(--text-dim); vertical-align:super; margin-left:0.6rem; }
    .src-blurb { font-size:0.95rem; color:var(--text-dim); max-width:760px; line-height:1.8; margin-bottom:1.8rem; }
    .src-list { list-style:none; margin:0; padding:0; }
    .src-entry { border-left:2px solid var(--rim); padding:0 0 0 1.4rem; margin:0 0 2rem; }
    .src-entry:target { border-left-color:var(--accent); }
    .src-head { display:flex; flex-wrap:wrap; align-items:baseline; gap:0.6rem; }
    .src-cite { font-size:1rem; color:var(--text); line-height:1.6; }
    .src-cite em { color:var(--text-bright); font-style:italic; }
    .src-gold { font-family:var(--mono); font-size:0.55rem; letter-spacing:0.18em; color:var(--accent-warm); border:1px solid var(--accent-warm); padding:0.12rem 0.4rem; white-space:nowrap; }
    .src-verified { font-family:var(--mono); font-size:0.55rem; letter-spacing:0.16em; color:var(--text-dim); border:1px solid var(--rim); padding:0.12rem 0.4rem; white-space:nowrap; }
    .src-v-extracted { color:var(--accent-green); border-color:var(--accent-green); }
    .src-v-crossref { color:var(--accent); border-color:var(--accent); }
    .src-v-unconfirmed { color:var(--danger); border-color:var(--danger); }
    .src-link { font-family:var(--mono); font-size:0.72rem; margin:0.5rem 0 0.7rem; word-break:break-all; }
    .src-link a { color:var(--accent); text-decoration:none; }
    .src-link a:hover { opacity:0.75; }
    .src-nolink { color:var(--danger); }
    .src-role { font-size:0.95rem; color:var(--text-mid); line-height:1.8; max-width:760px; margin:0 0 0.5rem; }
    .src-note { font-size:0.86rem; color:var(--text-dim); line-height:1.75; max-width:760px; margin:0 0 0.5rem; font-style:italic; }
    .src-meta { font-family:var(--mono); font-size:0.6rem; letter-spacing:0.08em; color:var(--muted); margin-top:0.4rem; }
    .src-tickets { color:var(--text-dim); }
    .src-local { color:var(--muted); }
    .src-key { color:var(--muted); }
    .src-v-publisher { color:var(--text-dim); border-color:var(--muted); }
    .src-provenance { margin:1rem 0 5rem; padding:1.4rem; border:1px solid var(--rim); background:var(--deep); }
    .src-provenance p { font-family:var(--mono); font-size:0.64rem; letter-spacing:0.04em; color:var(--text-dim); line-height:1.9; margin:0; word-break:break-all; }
    @media (max-width:600px) { .src-entry { padding-left:1rem; } }"""


def render(rows: list[dict[str, str]], digest: str) -> str:
    counts = {c: [r for r in rows if r["category"] == c] for c, _, _ in CATEGORIES}
    n_linked = sum(1 for r in rows if r["url"])
    n_local = sum(1 for r in rows if r["local_file"])
    n_gold = sum(1 for r in rows if r["role_in_codex"].startswith(GOLD_MARKER))

    L: list[str] = []
    a = L.append
    a("<!DOCTYPE html>")
    a('<html lang="en">')
    a("<head>")
    a('  <meta charset="UTF-8">')
    a('  <meta name="viewport" content="width=device-width, initial-scale=1.0">')
    a("  <title>Sources &middot; The Exoplanet Codex</title>")
    a('  <meta name="description" content="Every reference the Exoplanet Codex '
      'pipeline depends on, with a link to the original document. Generated from '
      'data/refs/bibliography.csv -- one source of truth, never hand-edited.">')
    a('  <meta name="author" content="Ryan Schmitt">')
    a('  <meta property="og:title" content="Sources &middot; The Exoplanet Codex">')
    a('  <meta property="og:description" content="Every reference behind the Codex '
      'measurements, linked to the original paper.">')
    a('  <meta property="og:url" content="https://exoplanetcodex.org/sources/">')
    a('  <meta property="og:type" content="website">')
    a('  <meta property="og:site_name" content="The Exoplanet Codex">')
    a('  <link rel="canonical" href="https://exoplanetcodex.org/sources/">')
    a('  <link rel="icon" href="/favicon.svg" type="image/svg+xml">')
    a('  <link rel="icon" href="/favicon.ico" sizes="any">')
    a('  <link rel="apple-touch-icon" href="/apple-touch-icon.png">')
    a('  <link rel="preconnect" href="https://fonts.googleapis.com">')
    # The other 22 site pages write these separators as a raw '&', which browsers
    # tolerate but which is not valid HTML. Emitting '&amp;' renders identically and
    # parses cleanly -- the only ampersand on the page that had to be escaped.
    a('  <link href="https://fonts.googleapis.com/css2?family=Space+Mono:ital,wght@0,'
      '400;0,700;1,400&amp;family=Crimson+Pro:ital,wght@0,300;0,400;0,600;1,300;1,400'
      '&amp;display=swap" rel="stylesheet">')
    a('  <link rel="stylesheet" href="/assets/css/codex.css">')
    a("  <style>")
    a(STYLE)
    a("  </style>")
    a("</head>")
    a("<body>")
    a("<!-- GENERATED FILE -- DO NOT EDIT BY HAND.")
    a("     generator : scripts/generate_sources_page.py (RYA-854)")
    a("     source    : data/refs/bibliography.csv")
    a(f"     sha256    : {digest}")
    a(f"     rows      : {len(rows)}")
    a("     Any hand edit is reverted by the next build and fails --check. Change the")
    a("     CSV, then regenerate. -->")
    a('<canvas id="starfield"></canvas>')
    a('<div class="page">')
    L.extend(render_nav())

    a('  <header class="page-hero">')
    a('    <div class="page-eyebrow">Sources // Every reference, linked to the '
      'original</div>')
    a('    <h1 class="page-title">The <em>bibliography</em></h1>')
    a('    <p class="page-subtitle">Every value the Codex publishes carries a cited '
      'source. This is that list &mdash; generated from one committed CSV, so the '
      'page cannot drift from the pipeline that reads it.</p>')
    a("  </header>")

    a('  <div class="src-intro">')
    a("    <p>An abundance is only as good as the atomic data underneath it. A log "
      "<em>gf</em> measured in a laboratory and a log <em>gf</em> fitted so that the "
      "Sun comes out right are <strong>not interchangeable</strong> &mdash; the second "
      "one makes a solar abundance a restatement of its own input. So the Codex tracks "
      "provenance per source, and this page publishes it.</p>")
    a("    <p><strong>Read the badges.</strong> "
      "<span class=\"src-verified src-v-extracted\">PDF</span> means the metadata was "
      "read off the paper itself; "
      "<span class=\"src-verified src-v-crossref\">DOI</span> means the DOI record was "
      "fetched and matched to this row; "
      "<span class=\"src-verified src-v-unconfirmed\">UNCONFIRMED</span> means no DOI "
      "or stable URL could be confirmed. Those rows are shown "
      "<strong>without a link and without a guessed DOI</strong>. Most are pre-DOI "
      "literature that no registry carries; a plausible-looking invented identifier "
      "would be a fabricated provenance record, so there isn&rsquo;t one.</p>")
    a("    <p><span class=\"src-gold\">GOLD</span> marks a reference standard: a scale "
      "the Codex compares against rather than re-derives.</p>")
    a("  </div>")

    a('  <div class="src-stats">')
    for value, label in ((len(rows), "references"),
                         (n_linked, "with verified link"),
                         (len(rows) - n_linked, "unconfirmed link"),
                         (n_local, "held locally"),
                         (n_gold, "reference standards")):
        a(f'    <div class="src-stat"><b>{value}</b><span>{label}</span></div>')
    a("  </div>")

    a('  <div class="src-toc">')
    for cat, heading, _ in CATEGORIES:
        a(f'    <a href="#{cat}">{heading} ({len(counts[cat])})</a>')
    a("  </div>")

    for cat, heading, blurb in CATEGORIES:
        group = sorted(counts[cat], key=lambda r: (r["authors"].lower(), r["year"]))
        a(f'  <section class="src-section" id="{cat}">')
        a(f'    <h2>{heading}<span class="src-count">{len(group)}</span></h2>')
        a(f'    <p class="src-blurb">{blurb}</p>')
        a('    <ul class="src-list">')
        for row in group:
            L.extend(render_entry(row))
        a("    </ul>")
        a("  </section>")

    a('  <div class="src-provenance">')
    a("    <p>GENERATED PAGE &middot; source: data/refs/bibliography.csv &middot; "
      f"sha256 {digest} &middot; {len(rows)} rows &middot; generator: "
      "scripts/generate_sources_page.py &middot; RYA-854. Never hand-edited: "
      "<code>--check</code> regenerates this file and fails the build on any "
      "difference.</p>")
    a("  </div>")

    a('  <footer class="codex-footer">')
    a('    <div class="footer-grid">')
    a('      <div><p class="footer-label">The Exoplanet Codex</p><p>Open science '
      '&middot; Stellar spectroscopy &middot; Montana, USA</p></div>')
    a('      <div><p class="footer-label">Navigate</p><a href="/glossary/">Glossary</a>'
      '<a href="/method/">Methodology</a><a href="/sources/">Sources</a>'
      '<a href="/mission-log/">Mission Log</a></div>')
    # Mirrors the footer the other 22 site pages carry, verbatim. Those link the v1.2
    # methodology / v1.0 glossary PDFs while the bibliography cites v2.0 -- a real site
    # staleness, but it belongs to the RYA-179 doc-sync pass, not to this page. Making
    # /sources the one page that disagrees would be worse than matching.
    a('      <div><p class="footer-label">Science</p>'
      '<a href="https://github.com/damienabraxas/exoplanetcodex">Pipeline Code</a>'
      '<a href="/assets/docs/exoplanet-codex-methodology-v1.2.pdf">Methodology PDF</a>'
      '<a href="/assets/docs/exoplanet-codex-glossary-v1.0.pdf">Glossary PDF</a>'
      '</div>')
    a('      <div><p class="footer-label">Contact</p>'
      '<a href="mailto:info@exoplanetcodex.org">info@exoplanetcodex.org</a></div>')
    a("    </div>")
    a('    <p class="footer-copy">&copy; 2026 The Exoplanet Codex &middot; Open '
      'Science &middot; All data public</p>')
    a("  </footer>")
    a("</div>")
    a('<script src="/assets/js/starfield.js"></script>')
    a('<script src="/assets/js/codex.js"></script>')
    a("</body>")
    a("</html>")
    return "\n".join(L) + "\n"


# ────────────────────────────────────────────────────────────── link checking ──
def verify_links(rows: list[dict[str, str]]) -> int:
    """Opt-in, NETWORKED: resolve every non-blank url and report the status.

    Deliberately NOT part of --check. A build gate that needs the internet fails for
    reasons that have nothing to do with the repo, and would eventually be disabled.
    Run this by hand when the bibliography changes.

    Publisher sites (aanda.org in particular) return 403 to automated agents. A 403 is
    a bot block, not a dead link, so it is reported separately from a 404 -- the DOI is
    still registered, which is what the row claims.
    """
    import urllib.error
    import urllib.request

    ua = "exoplanetcodex-RYA854 link check (mailto:ryan.damien.schmitt@gmail.com)"
    bad: list[str] = []
    blocked: list[str] = []
    for r in sorted(rows, key=lambda x: x["key"]):
        if not r["url"]:
            print(f"  ---- {r['key']:28s} (no link by design: {r['verified']})")
            continue
        req = urllib.request.Request(r["url"], headers={"User-Agent": ua},
                                     method="HEAD")
        try:
            with urllib.request.urlopen(req, timeout=45) as resp:  # noqa: S310
                code = resp.status
        except urllib.error.HTTPError as exc:
            code = exc.code
        except Exception as exc:                                   # noqa: BLE001
            code = f"ERR {exc}"
        print(f"  {str(code):4s} {r['key']:28s} {r['url']}")
        if code in (401, 403, 429):
            blocked.append(f"{r['key']} -> {code}")
        elif not (isinstance(code, int) and 200 <= code < 400):
            bad.append(f"{r['key']} -> {code}")

    print(f"\n{len(rows)} rows; {sum(1 for r in rows if not r['url'])} unlinked by "
          f"design; {len(blocked)} bot-blocked; {len(bad)} unreachable")
    for b in blocked:
        print(f"  bot-blocked (DOI still registered): {b}")
    for b in bad:
        print(f"::error::unreachable link: {b}")
    return 1 if bad else 0


# ───────────────────────────────────────────────────────────────────── driver ──
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="regenerate in memory and diff against the committed page; "
                         "exit 1 on any drift (a hand edit is a build break)")
    ap.add_argument("--out", type=Path, default=PAGE,
                    help=f"page to write (default {PAGE.relative_to(ROOT)})")
    ap.add_argument("--csv", type=Path, default=BIB_CSV,
                    help=f"bibliography (default {BIB_CSV.relative_to(ROOT)})")
    ap.add_argument("--verify-links", action="store_true",
                    help="NETWORKED: resolve every url in the bibliography and report "
                         "its status. Not part of --check, on purpose.")
    ap.add_argument("--site-root", type=Path, default=None,
                    help="also write <site-root>/sources/index.html, for deploying "
                         f"into a website checkout (sibling default: "
                         f"{SITE_ROOT_DEFAULT.name})")
    args = ap.parse_args(argv)

    if not args.csv.exists():
        sys.exit(f"::error::bibliography not found: {args.csv}")
    rows = load_rows(args.csv)

    if args.verify_links:
        return verify_links(rows)

    page = render(rows, csv_digest(args.csv))

    if args.check:
        if not args.out.exists():
            print(f"::error::{args.out} does not exist -- run without --check first",
                  file=sys.stderr)
            return 1
        current = args.out.read_text(encoding="utf-8")
        if current == page:
            print(f"OK: {args.out.relative_to(ROOT)} matches "
                  f"{args.csv.relative_to(ROOT)} ({len(rows)} rows)")
            return 0
        diff = list(difflib.unified_diff(
            current.splitlines(keepends=True), page.splitlines(keepends=True),
            fromfile="committed page", tofile="regenerated from CSV", n=2))
        sys.stdout.writelines(diff[:200])
        print(f"::error::{args.out.relative_to(ROOT)} does not match the "
              f"bibliography. The page is GENERATED: edit "
              f"{args.csv.relative_to(ROOT)} and regenerate, never the HTML.",
              file=sys.stderr)
        return 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(page, encoding="utf-8")
    print(f"wrote {args.out.relative_to(ROOT)} ({len(rows)} rows, "
          f"{sum(1 for r in rows if not r['url'])} without a verified link)")

    if args.site_root is not None:
        site_root = args.site_root if args.site_root.is_absolute() else \
            (Path.cwd() / args.site_root).resolve()
        if not (site_root / "assets" / "css" / "codex.css").exists():
            sys.exit(f"::error::{site_root} does not look like the website checkout "
                     f"(no assets/css/codex.css) -- refusing to write into it")
        target = site_root / SITE_PAGE_REL
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(page, encoding="utf-8")
        print(f"wrote {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
