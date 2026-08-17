# `data/refs/` — the bibliography (RYA-854)

`bibliography.csv` is the **single source of truth for every reference the Codex
depends on**. The public `/sources` page is GENERATED from it and nothing else:

```bash
python3 scripts/generate_sources_page.py            # write docs/site/sources/index.html
python3 scripts/generate_sources_page.py --check     # regenerate + diff; exit 1 on drift
python3 scripts/generate_sources_page.py --verify-links   # NETWORKED, opt-in
python3 scripts/generate_sources_page.py --site-root ../exoplanetcodex-site
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

### Why `unconfirmed` rows exist, and why they are not fixed by guessing

Twelve rows carry no link. Most are pre-DOI literature — 1970s–1990s A&A volumes, NBS
monographs, JPCRD supplements, SAO Special Reports — that no registry indexes. For
seven of them the *title* is blank too: the Gaia-ESO v5 line-list report substantiates
authors/year/journal/volume/page and nothing more, so a title would have been invented.

A fabricated DOI is a fabricated provenance record, and it is worse than an absent
one because it *looks* checkable. The page renders these rows without a link and shows
the reason. That is the data-stewardship position, not an unfinished task.

## Verification state as of RYA-854 (2026-08-17)

* 99 rows; 0 in a `verify_*` state.
* 81 DOIs, **all 81 present in the Crossref registry** (excluding the NIST `10.18434`
  and arXiv `10.48550` DOIs, which are not Crossref-registered and were confirmed
  against their own resolvers).
* `--verify-links`: **0 unreachable**. 41 rows return HTTP 403 — the publisher's
  anti-bot response (aanda.org, OUP, Wiley), not a dead link; the DOI is registered,
  which is what the row claims. A 403 and a 404 are reported separately for exactly
  this reason.

## Relationship to `docs/references.md`

`docs/references.md` is the narrative, role-grouped reading guide. **This CSV is the
authority for author/year/venue/DOI.** Where the prose and the CSV disagree, the CSV
wins and the prose is the bug. Instrument and archive citations (HARPS, UVES,
ESPRESSO, SPIRou, CRIRES+, MAST/ESO archives) still live only in `references.md` and
are a flagged follow-on batch for this CSV.
