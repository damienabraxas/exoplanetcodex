# The Exoplanet Codex Science Product Package (SPP)

## What this document is

The Exoplanet Codex ships one Science Product Package (SPP) per released star. This document defines what an SPP contains, how it is versioned, and how it reads. It is the standing framework — apply it to every star, from the Sun through 55 Cancri A and beyond.

The SPP has two audiences at once: astronomers who need to see the methods, the diagnostics, and the caveats; and non-astronomers who need to understand what we measured and how confident we are. It is written for a smart reader who is new to the topic — not a field specialist, not a lay audience needing everything explained. Halfway between, with clear explanations, load-bearing numbers, and no hiding.

## What an SPP contains

Every SPP has the same twelve sections. The order does not change; a reader who knows one SPP can navigate any other.

1. **Star identity** — HD/HIP/common name, coordinates (RA/Dec, J2000), spectral type, distance, notable planets or features.
2. **Provenance header** — Underlying frozen reference artifact (path + sha256), register version, SPP version, generator commit, release date. Every value in the report traces to a specific version of a specific artifact.
3. **Stellar parameters** — Teff, log g, [Fe/H], microturbulence, cited to primary literature. No project-guessed values.
4. **The 27-element abundance table** — For every element: measured value, statistical uncertainty, method (EW / synthesis), engine (Engine A / Engine B / arbiter), NLTE grid source and vintage, tier (`gold` / `owed` / `data-gap`), verdict (PASS / CURATION-OWED / NLTE-OWED / DATA-GAP). Held-at-owed values are shown with their caveats — not hidden.
5. **Flagship science products** — C/O ratio, alpha/Fe, s-process/r-process ratios where measured, connection to planetary composition where applicable.
6. **Cross-region abundance panel** — For every element measured across multiple wavelength regions, a small-multiples panel showing UV / VIS / IR values against the combined estimate, with an agreement badge.
7. **Anomalies and discrepancies** — Its own section, equal billing with the wins. Every blend inflation, every held value, every element where two engines disagree, every line that gave us trouble. Peer-review-ready means peer-reviewable — a report that buries the mess is not science, it is marketing.
8. **Line-level diagnostics** — Plots of every problem line, every blend we encountered, every line of interest worth documenting. Curve-of-growth plots for elements with visible scatter. HFS synthetic fits for the odd-isotope species (Mn, Co, Cu, V, Sc, Ba II, Eu II). Blend traces for anything under investigation.
9. **Uncertainty budget** — For each element, what dominates the reported error: line count, continuum, NLTE grid interpolation, gf-value floor, 1D-vs-3D scale offset. If the budget is dominated by a defect we know about, name the defect.
10. **Peer-review caveats** — Where our method is limited. 1D-NLTE-vs-3D scale offsets. Scaled-Drawin atom annotations where they apply. NLTE grid vintages. Regions we chose not to trust.
11. **Data provenance** — Instruments, dates, archive IDs (ESO / MAST / KOA / CADC / TNG). Telluric correction status for any IR arm.
12. **Method provenance** — Full engine + NLTE grid + atom + line list stack, per element. Every choice is discoverable in one table.

## Voice — the standing rules

The voice guide applies to every SPP, every mission log, every newsletter, every outreach post.

**Plain words over jargon.** If a concept has a plain-English name, use it; introduce the technical term in parentheses on first use. "How much iron the Sun has, compared to what the field has settled on" beats "the Fe I ionization-corrected 1D-NLTE anchor."

**Lead with what changed and why.** Not "recent advancements in the field of stellar spectroscopy have enabled..." — instead "We measured 27 elements in the Sun. Here is what we got, where we agree with the literature, and where we don't."

**Show the mess.** Discrepancies, held values, blend inflations, and elements we do not yet trust get equal billing with the wins. A report that hides the trouble spots is not science, it is marketing.

**No hedging as decoration.** "May be" and "could potentially" only when the uncertainty is real. If we are confident, say so.

**No AI-slop tells.** No "it is worth noting that." No "in the ever-evolving landscape of." No "delve." No adjective stacks. No three-item lists where two would do. If a sentence could be cut by half and lose nothing, cut it.

**One idea per sentence.** Long sentences hide vagueness.

**Numbers are load-bearing.** Every value carries its source. No naked "roughly" or "approximately" without the actual number.

**Assume a smart reader who is new to the topic** — not a field specialist, not a lay reader who needs everything explained. Halfway between: someone motivated to understand, who deserves clear explanations.

## Document control

SPPs are **write-once, hashed, and immutable**, on the same discipline as `data/reference/solar/`. Once an SPP releases at version N, it cannot be edited — a correction produces version N+1 with a written change note. This protects the peer-review record: what we said on release date is exactly what we said.

**Versioning rule:** The SPP version tracks the underlying frozen artifact version. Sun SPP v3 is generated from `solar_abundances_v3.csv`. When gold v3 becomes v4 from a future re-run, the Sun SPP re-emits as v4. This coupling is machine-checkable: the SPP header carries the artifact sha256, and the generator refuses to build an SPP against a hash that no longer resolves.

**Storage:** `docs/products/<star>/spp_<star>_v<N>.pdf` (the PDF itself), `docs/products/<star>/plots/` (supporting figures), `docs/products/<star>/tables/` (abundance table + method-provenance table in CSV, for downstream tooling). A `docs/products/<star>/CHANGELOG.md` file records every version bump.

**Hash manifest:** `docs/products/<star>/hash_manifest.json`. Every artifact in the star's product folder is hashed at release. A future integrity check can verify the PDF you downloaded is the PDF we released.

**Immutability rule:** Editing a released SPP in place is a hard error — the generator refuses to overwrite an existing versioned file. Corrections produce a new version. If the correction is a typo or cosmetic fix, note it in the CHANGELOG; if it is a scientific correction, the changelog entry describes what changed and why.

## Distribution

For each released SPP, the following artifacts land together:

- **The PDF** — in `docs/products/<star>/spp_<star>_v<N>.pdf`, and attached to the release Linear ticket for archival redundancy.
- **The web bundle** — plots and tables in `docs/products/<star>/plots/` and `docs/products/<star>/tables/`, ready for ChatGPT Codex to pull into the star's website page. Pull directly from main, or hand-deliver on request.
- **A mission log entry** — blog-post format at `docs/products/<star>/mission_log_v<N>.md`. The story of the run: what worked, what surprised us, what took the longest. Written in the SPP voice.
- **Social posts** — templated per platform, drafted alongside the SPP release. Link back to the star's website page and the PDF.
- **A newsletter section** — Buttondown post covering the release, drafted alongside the mission log.
- **Future outreach materials** — as they are added, they attach to the same release folder.

## The release cadence

A star's SPP release follows the pipeline landing. For the Sun, the sequence is:

1. RYA-527 Phase 2 emits fresh two-engine artifact + verdict + disposition report.
2. Gold v4 freeze if material changes surface, otherwise v3 stands.
3. Sun SPP v3 or v4 generated from the frozen artifact.
4. RYA-179 doc-sync (Glossary, Method, Science Architecture) released in parallel.
5. Mission log + social + newsletter drafted from the SPP.
6. Web bundle handed to ChatGPT Codex.
7. Public release — one moment, all artifacts landing together.

The framework doc itself (this file) is peer to the three docs in RYA-179 (Glossary, Method, Science Architecture). The four together form the Codex's public documentation surface.

## What lives outside the SPP

Some things belong in other files, not the SPP:

- **The Codex Glossary** carries term definitions and acronyms.
- **The Method page** carries the pipeline architecture (VALD, normalize, EW / synth, NLTE, gold reference), general enough to apply to every star.
- **The Science Architecture doc** carries the two-engine design, the three-gate promotion, the disposition report structure — how the Codex thinks.
- **`CODEX_STATE_REGISTER.md`** carries the internal engineering state — never a public artifact.

The SPP is the star-specific narrative and data. Everything general lives in the peer docs above.

## Change control for this framework doc

This document itself is subject to normal RYA-179-style discipline. Changes to the framework require a ticket, a PR, and a register bump. When the framework evolves, existing released SPPs are NOT retrofitted — they released under the framework in effect at their release date. A new SPP framework version applies only to SPPs released after it. This preserves the peer-review record.
