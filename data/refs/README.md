# `data/refs/` — the bibliography (RYA-854)

`bibliography.csv` is the **single source of truth for every reference the Codex
depends on**. The public `/sources` page is GENERATED from it and nothing else:

```bash
python3 scripts/generate_sources_page.py            # write docs/site/sources/index.html
python3 scripts/generate_sources_page.py --check     # regenerate + diff; exit 1 on drift
python3 scripts/generate_sources_page.py --verify-links   # NETWORKED, opt-in
python3 scripts/generate_sources_page.py --site-root ../exoplanetcodex-site
python3 scripts/generate_sources_page.py --audit-library "<library dir>"   # coverage
```

**Never hand-edit `docs/site/sources/index.html`.** `--check` regenerates the page and
fails on any difference, so a hand edit is a build break, not a silent divergence.
Edit the CSV and regenerate. This is the RYA-775 pattern: hand-authored reference
pages rot, and nothing reports that they have.

## Schema

| column | meaning |
|---|---|
| `key` | stable citation handle, unique. Cite this in code comments and tickets. |
| `authors` | `Family, I.; Family, I.` — up to ~6, no `et al.` padding |
| `year` | the year a citation would use (the **volume** year, not Crossref's online-first date; where they differ the row says so in `license_note`) |
| `title` | the paper's title. **Blank where no authoritative record of the title could be reached** — see "unconfirmed" below |
| `venue` | journal/volume/page, or the imprint for a monograph |
| `doi` | bare DOI, no resolver prefix. Blank when none is registered |
| `url` | stable link to the ORIGINAL. `https://doi.org/…` preferred, `https://arxiv.org/abs/…` acceptable, database/tool URLs for NIST/VALD/Kurucz. **Blank rather than guessed** |
| `category` | one of the twelve values below |
| `role_in_codex` | what this source actually does for the project — the sentence that stops the row being decoration |
| `tickets` | `RYA-…; RYA-…` — the tickets that touched or depend on it |
| `local_file` | path inside `~/Documents/Exoplanet Codex/` for the copy on disk. Working papers, **not** redistributable repo assets |
| `license_note` | licensing note, and the place every verification caveat and correction is recorded |
| `verified` | how the row was confirmed — see below |

## `category` vocabulary

`solar_reference`, `atomic_data_gf`, `atomic_data_broadening`, `nlte`,
`nlte_synthesis`, `model_atmosphere`, `benchmark_params`,
`instrument_line_selection`, `cno`, `motivation`, `method_tool`, `project_doc`.

The generator holds this list with each category's display heading and ordering. A
category not in that list **fails the build** — adding one is a deliberate edit to
`scripts/generate_sources_page.py`, never a typo that silently drops rows off the page.

## `verified` vocabulary

| value | what was actually done |
|---|---|
| `extracted` | metadata read off the PDF held in the local reference library — the strongest state, because the paper itself was inspected |
| `crossref` | the DOI was fetched from `api.crossref.org` and its title/authors/venue matched this row |
| `publisher` | a publisher or database landing page was resolved and matched (NIST ASD, VALD, kurucz.harvard.edu, exoplanetcodex.org) |
| `unconfirmed` | **no DOI or stable URL could be confirmed.** `doi` and `url` are BLANK and `license_note` says why |

RYA-854 spec: `verify_doi` / `verify_ref` are *pending* states and **must not appear**.
The generator fails on them.

### `verified` and `url` are two axes, not one

`verified` records how the **metadata** was confirmed. `url`/`doi` record whether a
**public link exists**. They are independent, and the generator's invariant runs one
way only: `unconfirmed` ⇒ no doi and no url. The converse is *not* asserted, because an
unpublished project document read straight off its PDF is `extracted` and legitimately
has no link. (The first version of that check asserted both directions and forced a
false `unconfirmed` badge onto exactly such a row.)

### Where a pre-DOI source does have a stable link

Most pre-DOI literature *can* be linked without inventing an identifier — it just isn't
through `doi.org`:

* **ADS scanned articles** — `https://articles.adsabs.harvard.edu/pdf/<bibcode>`. Use
  this, **not** `ui.adsabs.harvard.edu/abs/<bibcode>`: the abstract page is a
  single-page app that returns HTTP 202 with an empty body for a *fabricated* bibcode
  just as readily as for a real one, so it cannot verify anything. The scan route
  returns 200 for a real bibcode and 404 for a fabricated one — it discriminates.
* **SIMBAD** — `https://simbad.cds.unistra.fr/simbad/sim-ref?bibcode=<bibcode>` returns
  server-rendered title/authors/pages. This is how the pre-DOI titles here were
  confirmed, and it corrected two that had been drafted from memory.
* **CDS VizieR** for the machine-readable table (e.g. `VI/10` for Kurucz &
  Peytremann 1975), **Internet Archive** for NBS monographs, and the **proceedings
  host** for workshop papers.
* **DOE OSTI** (`osti.gov/biblio/<id>`, JSON at `osti.gov/api/v1/records/<id>`) for US
  government-funded reports and monographs. This is what finally supplied the title,
  the Suppl. 3 number and a link for `martin1988` after Crossref, OpenAlex, NIST's own
  publications database, ADS, LoC, HathiTrust and Google Books all came up empty.

### Why ONE row still carries no link, and why that is not fixed by guessing

`fuhr1988` — *Atomic Transition Probabilities: Iron through Nickel*, JPCRD 17 Suppl. 4.
Title and volume are confirmed twice over, but no stable URL exists in any catalogue
checked. Its companion `martin1988` (Suppl. 3) *was* found in DOE OSTI, so the gap here
is a hole in that catalogue rather than a failed search.

**A blank `url` here is probably the correct permanent state, not an open task.** The
volume appears to be print-only — consistent with it being absent from every digital
catalogue tried, and with AIP's legacy supplement links returning 503. A row that says
"this exists in print and has no online copy" is complete; a row pointing at a *different*
document to fill the field would not be.

**JPCRD supplements are standalone monographs outside the article stream and were never
Crossref-registered.** Proved by enumerating all 103 JPCRD items 1987–1990 from
Crossref: neither 1988 supplement is among them. Don't re-derive that.

Two cautionary notes kept here on purpose, because both were **my** errors and both are
the same shape — trusting a *name* instead of opening the *document*:

* `schmitt_science_architecture` was recorded as *missing* because a search for the
  inventory's filename (`Exoplanet_Codex_Science_Architecture.docx`) came up empty. The
  document exists, published on the site as v3.0. **The absence was scoped to a
  filename, not to the document** — `feedback_absence_is_a_hypothesis`.
* The 2010 senior thesis was carried as **two rows**, `schmitt_beyond_metallicity` and
  `schmitt2010_thesis`, because the same file sits on disk under three names. The
  thesis row also carried a title I had paraphrased from the About page's prose rather
  than read off the PDF. `--audit-library` now hashes every `local_file` and fails when
  two rows name byte-identical content.

`martin1988` also still has a blank `title`: unlike its Fe–Ni companion, no primary
source held here gives one. Secondary sources agree on "Atomic transition
probabilities — Scandium through manganese", and the note says so, but that is not
promoted into the `title` field.

A fabricated DOI is a fabricated provenance record, and it is worse than an absent one
because it *looks* checkable. The page renders these rows without a link and shows the
reason. That is the data-stewardship position, not an unfinished task.

## Every document in the reference library must have a row

Ryan's standing rule: **an unlisted document is an uncited one.** That is checkable, not
just asserted:

```bash
python3 scripts/generate_sources_page.py --audit-library "<reference library dir>"
```

It fails **in both directions**: if any *document* in the library is not named by some
row, and if any two rows name **byte-identical** files.

🔴 **Run it EVERY SESSION.** The library gains documents between sessions — three
arrived mid-session on 2026-08-17 while this file was being written, and the audit is
the only reason they were noticed. The step is listed at the top of
[`LEDGERS.md`](../../LEDGERS.md) (session start) and in `DEV_CYCLE.md` §7 (session
close). There is no CI guard, because the library lives outside the repo and differs
per machine — so this is discipline, and the audit is what makes the discipline cheap. The second direction exists because a filename check cannot see that one
document is sitting on disk under three different names — only content can.
The path is an argument, never hardcoded — the library lives outside the repo on a
per-machine path, and a literal would both break on Sirius and trip the RYA-810 gate.

Non-reference assets (`.jpg`, `.png`, `.heic`, `.xml`, …), including project design
diagrams, are skipped **by extension**, never by a filename allow-list. They remain
project evidence but are not literature citations; an allow-list would rot silently
the moment a file is renamed.

## How agents use the local reference library

The library is a durable evidence holding, not a source of project law. Check existing
holdings before reacquiring a paper. For each scientific source record title, authors,
year, DOI or equivalent stable identifier, local filename, and relevant Codex role where
practical. Verify metadata and any page-specific interpretation by inspecting the source;
never invent a citation, DOI, page number, or reading from a filename. Ratified method
decisions belong in `docs/SCIENCE_STANDARDS.md` or a decision record and merely cite the
paper. Project diagrams and legacy instructions are governed as design/context artifacts,
not forced into the bibliography.

**Duplicates are settled by content, never by name.** A row's `local_file` names one
file and its `license_note` records any twin. This matters: `1002.4268v1.pdf` and
`BRUNTT2010.pdf` are the *same bytes twice* (md5 match), while `2602.14294v1.pdf` and
`aa59148-26.pdf` are the *same paper in two genuinely different files* (arXiv preprint
vs published A&A). Only a content check tells those apart, and they need opposite
treatment.

## Verification state as of RYA-854 (2026-08-17)

* **104 rows**; 0 in a `verify_*` state.
* 82 DOIs, **all present in the Crossref registry** (excluding the NIST `10.18434`
  and arXiv `10.48550` DOIs, which are not Crossref-registered and were confirmed
  against their own resolvers).
* **103 of 104 rows carry a verified link**; 1 carries none, for the reason above.
* `--verify-links`: **0 unreachable**. ~41 rows return HTTP 403 — the publisher's
  anti-bot response (aanda.org, OUP, Wiley), not a dead link; the DOI is registered,
  which is what the row claims. A 403 and a 404 are reported separately for exactly
  this reason.
* `--audit-library`: **57 documents in the library, 57 accounted for**; 7 images
  skipped.

## Relationship to `docs/references.md`

`docs/references.md` is the narrative, role-grouped reading guide. **This CSV is the
authority for author/year/venue/DOI.** Where the prose and the CSV disagree, the CSV
wins and the prose is the bug. Instrument and archive citations (HARPS, UVES,
ESPRESSO, SPIRou, CRIRES+, MAST/ESO archives) still live only in `references.md` and
are a flagged follow-on batch for this CSV.
