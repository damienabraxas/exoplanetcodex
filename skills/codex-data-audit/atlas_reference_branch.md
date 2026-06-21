# codex-data-audit — Atlas / Reference Data branch

> **Status (RYA-392):** the canonical `codex-data-audit` `SKILL.md` is **not on disk**
> in this repo or Drive (same situation flagged by RYA-388/386/370 — the skill lives
> outside the repo). This file is the **new "Atlas / Reference Data" branch**, written
> to be dropped into that `SKILL.md` verbatim. Until the canonical file is on disk, this
> is the source of truth for the branch. Do not fork a second copy of the whole skill.

The existing `codex-data-audit` branches target **instrument spectra** (product type,
BERV/`SPECSYS`, units trap, SNR floor, product-suffix family). **Reference atlases**
(ACE-FTS, NSO photatl, Wallace telluric; Kitt Peak optical, BASS2000, UV composites)
are a different animal: there is no BERV or exposure, but there *is* a unit convention,
a provenance pedigree, and a disk/telluric character that silently corrupt downstream
science if mis-tagged. This branch is reusable for everything RYA-162 still needs.

Reference implementation: `pipeline/audits/audit_reference_atlases.py` (RYA-392).

## When to use this branch
The dataset is a **published reference atlas / template**, not an observation of a
science target: a solar atlas, a telluric transmission atlas, a model template, a
composite. Telltales: native axis in **wavenumber (cm⁻¹)**, a citation/DOI, a "disk"
character, FTS provenance.

## Checks (every one asserted, never inferred-and-forgotten)

1. **Provenance + citation.** Source, primary citation, version/date, SHA-256, and the
   retrieval URL must be present — in a provenance sidecar (e.g.
   `*_provenance_*.json`) keyed by the stored filename, or an inline header.
   **Loud-fail** if an atlas file has neither. Nothing from memory: every value cited
   to its source.

2. **Axis convention — wavenumber vs wavelength, and air vs vacuum.** Classify the axis
   from its value range *and* its column label; `UNKNOWN` is a loud failure, never a
   silent guess. FTS atlases are **vacuum wavenumber (cm⁻¹)** natively. If the file
   also carries wavelength columns, verify them rather than trust the tag:
   - `wavelength_vac == 1e8 / wavenumber` to < 0.01 Å, and
   - `wavelength_air < wavelength_vac` by the expected offset (~6.3 Å near 2.3 µm,
     Birch & Downs 1994 / Edlén — the VALD3/iSpec convention). A swapped air/vac tag
     shows up here as the wrong sign or magnitude.
   Verify the cm⁻¹→wavelength conversion against a **known feature** (e.g. the CO (2-0)
   bandhead ≈ 4360 cm⁻¹ ≈ 2.293 µm). A missing or wrong convention tag is a loud fail.

3. **Coverage over the science segment.** The atlas must span the segment the science
   needs (e.g. 4255–4367 cm⁻¹ for the K-band CO arm). **Role-aware:** a *documented*
   partial product (e.g. the Wallace ASCII telluric ratio covers only the band middle
   4299.8–4338.6 cm⁻¹) is checked against its stated sub-range and the redirect for the
   wings (here: photatl `atmospheric` column) — it is **not** false-failed for missing
   coverage that is documented and covered elsewhere. Cross-check the measured range
   against the provenance-stated range (catches a stale sidecar). Allow an endpoint
   tolerance (~1 cm⁻¹) to absorb the FTS sampling grid.

4. **Telluric classification — free / residual / pure.** Read from provenance, not
   assumed:
   - **free** — space/telluric-corrected solar truth (ACE-FTS). For ACE specifically,
     confirm it is the **Hase+2010 WSpectra "complete solar spectrum" derived product**,
     NOT a raw ACE **occultation** transmission (occultation looks *through* the
     atmosphere → carries tellurics → silently poisons the telluric-free truth).
     Confirm from the source URL + citation + filename, not the value range.
   - **residual** — ground-based solar atlas with leftover tellurics (photatl: `solar`
     telluric-corrected / `atmospheric` residual / `total` observed; note the `solar`
     column is linearly interpolated across strong-telluric gaps).
   - **pure** — a telluric-only transmission atlas (Wallace ratio).

5. **Disk character — center vs integrated.** Record disk-center vs disk-integrated and
   **flag any mismatch with the science target** as a known systematic — do not
   silently treat them as equivalent. Example: NSO photatl is **disk-center**;
   reflected-Vesta is **integrated-disk**; center-to-limb variation changes CO line
   depths, so photatl is a caveated cross-reference, not an equal of the integrated
   truth.

## Permanent rules
- No silent unit assumptions — `UNKNOWN` axis or a missing convention tag is a **loud
  failure** (do not warn-and-continue).
- Loud-fail if any atlas lacks a provenance header/sidecar entry.
- Audit data **before** any downstream result leans on it; a downstream correlation is
  only as trustworthy as the atlas intake behind it.
