---
name: codex-vald-extraction
description: The complete VALD3 line-list extraction and intake procedure for the Exoplanet Codex — the Extract Stellar form recipe, the mandatory HFS-splitting convention, the 100k truncation trap, per-system file naming, and the intake verification (truncation / coverage / HFS / threshold / ACCEPT-REJECT) that every raw delivery must pass before it is merged or built into a linelist. Use this skill whenever Ryan is preparing a VALD extraction, asks for "VALD parameters", "VALD inputs", "the extraction settings", "what do I put in the form", mentions submitting or downloading from vald.astro.uu.se, OR whenever a raw VALD delivery arrives and needs to be checked/unpacked/verified before use. Always use this skill before writing any VALD submission table or any intake/merge brief — the extraction recipe and the intake gate both live here so they cannot drift apart again.
---

# Codex VALD Extraction & Intake Skill

## Authority boundary

This skill owns only VALD extraction and raw-delivery intake. Resolve element execution
from `docs/ELEMENT_PROTOCOL.md`, scientific rules from `docs/SCIENCE_STANDARDS.md`, naming
and schemas from `docs/CONVENTIONS.md`, and current state from `LEDGERS.md`. Those sources
override stale examples here. Never turn a VALD delivery or a plausible `log_gf` into a
scientific decision, and never silently reconcile divergent oscillator strengths.

## Why this skill exists

The VALD procedure was previously scattered across Linear comments (RYA-64 how-to, RYA-269 HFS ruling, the intake check). That fragmentation is exactly how the **HFS-splitting setting went unnoticed for months** — the standard recipe never mentioned the field, so it was never checked, and mixed HFS-ON/OFF extractions could not be merged (RYA-269, 204 cross-extraction log_gf conflicts). This skill is the single source of truth. If the recipe changes, it changes here.

**Hard rule reminder:** VALD extraction is **always manual** via the web interface. There is no API path. The implementing agent never downloads from VALD — Ryan submits, Ryan drops the raw files, and the agent verifies and processes them.

---

## Part A — Submitting an extraction (Ryan, at vald.astro.uu.se → Extract Stellar)

### The standard form recipe

| Form field | Value | Notes |
|---|---|---|
| Extraction mode | **Extract Stellar** | not Extract Element / Extract All |
| Teff | per star — **from `config/stars.yaml`** | never hardcoded; see "Stellar parameters" below |
| log g | per star — **from `config/stars.yaml`** | |
| [Fe/H] (metallicity) | per star — **from `config/stars.yaml`** | |
| Microturbulence | per star — **from `config/stars.yaml`** | |
| Wavelength start / end | per request | split to dodge the 100k cap — see below |
| **Detection threshold** | **0.001** | *This is `central_depth` in the output.* The form labels it "Detection threshold." 0.001 = synthesis-grade: include every line absorbing ≥0.1% of continuum at line centre. The blend-aware Turbospectrum engine (RYA-285) needs these weak lines; the old EW-era 0.05 is **under-deep** and silently drops blends + trace species (RYA-381). See the tuning note below. |
| **Hyperfine splitting (HFS)** | **ON — ALWAYS** | The field that bit us. See Part B. Verify it is checked before *every* submit. |
| Line list configuration | default (VALD3) | |
| Format | **Long format** | |
| Wavelength unit | **Angstrom** | |
| Delivery | **FTP** where offered | higher capacity than browser; parser handles FTP format |

### Detection threshold (`central_depth`) tuning

The canonical extraction depth is **0.001** (ratified RYA-387). This is a hard convention, not a tuning knob — the synthesis engine is blend-aware and the weak lines are load-bearing.

- `0.001` — **canonical / synthesis-grade.** Matches the optical core; captures blends and trace n-capture / light-element lines the EW era discarded. Dense ranges will truncate at the 100k cap — **the answer is finer wavelength chunks, never a shallower threshold** (see the trap below).
- `0.05` — **EW-era vestige; under-deep. DO NOT USE.** Drops blends and trace species (Zr, P, S, neutron-capture). Any delivery measured at 0.05 is a heterogeneous-list defect (RYA-381) and fails the intake threshold check (Part C). Re-extract at 0.001 with finer chunks.
- `0.01` — historical intermediate; also under-deep relative to 0.001. Not used for new extractions.

The de-facto detection depth of a delivered file = `min(central_depth)` across its lines (VALD does not output anything below the submitted cut). The intake gate measures this directly — see Part C, check 3.

### Stellar parameters

**Do NOT hardcode stellar parameters in this skill.** They have drifted before (RYA-332/388: stale Teff/log g/[Fe/H] copied by hand into the recipe diverged from the calibrated values). The single source is `config/stars.yaml`, read via `config.constants.get_star_params(star_id)` (RYA-298). The submission table is **generated from constants, never typed.**

```python
from config.constants import get_star_params
p = get_star_params('solar')        # 'solar' | 'procyon' | '55cnc_a' | 'alpha_cen_a' | 'alpha_cen_b'
#   returns a dict → p['teff'], p['logg'], p['feh_ref'], p['xi']   (paste these into the form)
```

When a new target is scoped, add it to `config/stars.yaml` with the literature source noted in the issue that introduces the star — not to a table here.

### Species list (paste into the element selection field)

```
Li 1, C 1, N 1, O 1, Na 1, Mg 1, Mg 2, Al 1, Si 1, Si 2, P 1, S 1, Ca 1, Ca 2, Sc 1, Sc 2, Ti 1, Ti 2, V 1, V 2, Cr 1, Cr 2, Mn 1, Mn 2, Fe 1, Fe 2, Co 1, Co 2, Ni 1, Ni 2, Cu 1, Zr 1, Zr 2, Sr 1, Sr 2, Y 1, Y 2, Ba 1, Ba 2, Eu 2, CH 1, CN 1, C2 1, NH 1, OH 1, CO 1
```

(40 atomic species/ions for the 27-element list + 6 molecular species for C/N/O synthesis context. Molecular lines from VALD are non-authoritative — flagged for ExoMol/Brooke replacement, RYA-197.)

### The 100k truncation trap

VALD's web output caps at **100,000 transitions** and silently truncates beyond it. A truncated file has `WARNING: Output was truncated to 100000 lines` as **line 1** instead of the metadata header.

At the canonical 0.001 depth, dense ranges (near-UV forests, cool metal-rich blue regions) **will** hit the cap — that is expected. **Mitigation is always finer chunking, never a shallower threshold.** Split wide ranges into multiple narrower requests. Useful split points:

- **2000 Å** — also the vacuum/air convention boundary (λ < 2000 Å is vacuum in VALD), so splitting here keeps each file single-convention.
- Skip already-covered dense blue regions by starting a supplement higher (e.g., a red supplement starting at 5000 Å).
- For a forest that still truncates, halve the sub-range again. Many small files at 0.001 are correct; one shallow file is not.

Canonical per-star chunk plan (RYA-383/386): **1150–2000 / 2000–3780 / optical / 6910–9500 / 9500–25000**, each chunked finer wherever it still truncates at 0.001. Cool metal-rich stars (e.g., 55 Cnc) have especially dense blue/UV forests and need the most sub-splits.

### Naming convention

Linelists are **per-system**, not master-global (RYA-64 decision, 2026-06-09): `linelist_<star>.csv` (e.g., `linelist_solar.csv`, `linelist_55cnc.csv`, `linelist_procyon.csv`). Raw deliveries live in `data/linelists/` (confirm with Ryan whether a `raw/` subfolder is preferred). Do not let download-tool default folders (e.g. `VALD/<star>/`) create a parallel convention.

---

## Part B — The HFS convention (PERMANENT)

**Hyperfine splitting ON is the permanent Codex extraction convention — every star, every range, every submit.**

### Why
For HFS-sensitive species (Mn, Co, V, Sc, Cu, Ba II, Eu II, Na, Al, K, La — odd-isotope nuclei), an unsplit total-gf record saturates faster than the real distributed hyperfine components, biasing EW-derived abundances for any line of meaningful strength. Split components are the correct radiative-transfer input. This matters directly for the s/r-process anchors (Ba=s, Eu=r) and several 55 Cnc problem lines (Cu I 4767, etc.).

### Consistency
Solar and 55 Cnc optical extractions are HFS-ON; the solar calibration that passed gates (RYA-238) was computed on this convention. Mixed HFS-ON/OFF extractions **cannot be merged** — split components and unsplit total-gf records double-count opacity where they overlap. This is a hard merge gate.

### Fe is the exception that hides the bug
Fe is effectively HFS-immune (its dominant isotope ⁵⁶Fe has nuclear spin I = 0 — no hyperfine structure). So Fe-only validation (solar Fe, Procyon Fe) is numerically unaffected by HFS state — which is exactly why the convention error went undetected. **Do not let an Fe-only "it passed" lull you** — the multi-element run is where HFS bites.

### Trust nothing about the default
The HFS box state has historically varied across sessions for reasons never pinned down. **Physically confirm the box is checked before every submit. Trust the form default in neither direction.** Verify delivered files at intake (Part C) regardless.

---

## Part C — Intake verification (implementing agent, on every raw delivery)

Run this on each file **as it is unpacked, BEFORE any merge or build consumes it.** The merge conflict gate is the backstop; intake is where wrong files get caught and bounced. No partial merges — a merge runs only when all required files for that star are ACCEPT.

### Per-file checks

**1. Truncation.** Line 1 must be the VALD metadata header, not `WARNING: Output was truncated to 100000 lines`. Truncated → REJECT.

**2. Coverage.** Report delivered min–max wavelength vs requested. Use the shared coverage helper:

```python
def vald_coverage(path):
    wls = []
    with open(path) as f:
        for l in f:
            if l.startswith("'") and l[1] not in (' ', '_'):
                try: wls.append(float(l.split(',')[1]))
                except (ValueError, IndexError): pass
    if not wls:
        raise RuntimeError(f"{path}: no transitions parsed — wrong format or empty file")
    return min(wls), max(wls), len(wls)
```

**3. Threshold consistency.** The delivery's effective detection depth = `min(central_depth)` of its lines (field index 13). It must match the canonical **0.001**. Use the shared check (RYA-389):

```python
from vald_parse import verify_extraction_threshold
verdict, msg, eff = verify_extraction_threshold(path)   # 'ACCEPT' | 'FLAG' | 'REJECT'
```

- `ACCEPT` — effective depth ≤ 0.001 × tolerance (synthesis-grade). Good.
- `FLAG` — under-deep (e.g. 0.05 EW-era). The blend/trace lines are missing → heterogeneous-list defect (RYA-381). Re-extract at 0.001 with finer chunks; **do not** merge an under-deep file into a synthesis-grade linelist.
- `REJECT` — no parseable transitions (wrong format / empty).

The CI invariant `check_vald_threshold` (scripts/check_stewardship.py) enforces this on every committed `vald_*.txt`: any under-deep delivery must be tracked against its re-extraction ticket, or CI fails loud.

**4. HFS verification.** The **split-group signature is the primary, universal test**; Li 6707 is a fast confirmatory check valid ONLY where lithium survives the central-depth cut (see warning below).

- **PRIMARY — split-group signature (all files):** for each HFS-capable species present (Mn I/II, Co I/II, Cu I, V I/II, Sc I/II, Ba II, Eu II, Na I, Al I/II), measure the fraction of that species' lines sitting in split groups (≥2 same-species components within a tight wavelength distance, e.g. ±0.05 Å). **HFS-ON reads high (typically ~90%+ for odd-isotope-dominated species); HFS-OFF reads ~2%.** This is decisive and works regardless of wavelength range or stellar type. Report the per-species split fractions.
- **CROSS-CHECK (where an HFS-OFF reference overlaps):** a new HFS-ON file must have strictly MORE records than the quarantined HFS-OFF file in the overlap, the excess concentrated in odd-isotope species. Report the record deltas.
- **CONFIRMATORY — Li 6707 test (only where Li survives):** count records in 6707.5–6708.0 Å. ~15 components → ON; 2–4 → OFF. **⚠️ FAILS SILENTLY on Li-depleted stars.** Procyon (RYA-269) returned ZERO Li records even in a central_depth 0.01 extract — its depleted lithium never clears the depth cut, so a "0 records → OFF" reading is a **test misfire, not a bad file**. Hot stars, evolved stars, and any Li-poor target will do this. Use Li 6707 only to confirm an already-positive split-signature reading on a star known to retain lithium; never use it as the sole or primary test, and never read "0 records" as HFS-OFF.

**Decision:** trust the split-group signature. ON if split fractions are high across HFS-capable species AND (where checkable) the overlap delta is positive and odd-isotope-concentrated. A low split fraction on species that should split → OFF → REJECT.

**5. Verdict: ACCEPT / REJECT.** REJECT → quarantine the file (do not delete — provenance), report on the issue, Ryan resubmits with the box checked / threshold corrected. Report an intake table: file / coverage / truncation / threshold (effective depth) / HFS test used / result / verdict.

### Parsing & merge (when all files ACCEPT)

- Use the **shared parser** `data/linelists/vald_parse.py` (one parser for all stars — refactored from the RYA-223 Procyon logic in RYA-269). Never fork a second parser.
- Schema: `linelist_solar.csv` columns + `extraction_source` (provenance per line) + `wavelength_convention` (`vacuum` for λ < 2000 Å, `air` for λ ≥ 2000 Å — **do not convert vacuum→air at merge**; conversion happens downstream at the synthesis stage, not at the air/vacuum boundary).
- `blend_flag` = False on ingest (vetted exclusions only, RYA-209); proximity flagging is a separate downstream step (`vald_proximity_flag`).
- **Overlap dedup**: key on `(element, ion, round(λ,3), round(EP,3))`; keep the lower-threshold (more inclusive) source. **log_gf disagreement > 0.001 on a "duplicate" pair → CRITICAL, report, do not silently pick.** (This gate is what caught the HFS mismatch.)
- **Molecular species** (CH, CN, C2, NH, OH, CO): include but tag `notes = 'MOLECULAR — VALD non-authoritative, pending RYA-197'`.
- Parse-failure rate > 0.1% → CRITICAL.
- **Add gate**: assert zero HFS-convention conflicts after an HFS fix (the conflict count must come back 0).
- **Quarantine-source guard** (RYA-378): the linelist builder refuses any `*_quarantine` file as an input (`assert_not_quarantine_source` in `scripts/build_linelist.py`). Quarantined deliveries (HFS-off, truncated, under-deep) are provenance only — never an assembly source.
- Per-star canary gates where applicable: O I 7771/7774/7775 triplet present (NIR), Li I 6707 present.

### Quarantine / supersession

Superseded deliveries (wrong convention, truncated, under-deep, replaced by a re-extraction) are quarantined alongside the linelist, **never deleted** (provenance). Note the superseding file issue.

---

## Quick reference — building a submission table

When Ryan asks for "the VALD inputs," produce a table with: Teff / log g / [Fe/H] / Microturbulence (**from `config.constants.get_star_params(star_id)`, never typed by hand**), Wavelength start/end (from the split plan), Detection threshold **0.001**, and a reminder block: **HFS ON, FTP delivery, Long format, Angstrom, standard species list.** Always restate the HFS reminder — it is the field most likely to be skipped — and the 0.001 / finer-chunking rule, since the EW-era 0.05 is the most likely wrong value to creep back in.
