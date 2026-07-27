# Exoplanet Codex — engineering conventions

## Loaders declare + verify OBJECT, velocity frame, and wavelength scale (RYA-481 + RYA-264)

**Every loader must explicitly declare and verify three things about every frame
before that frame's photons reach any fit, EW, or synthesis: (1) its OBJECT
attribution, (2) its velocity reference frame, and (3) its wavelength unit+scale.
Selection is by authoritative FITS header — never by folder, filename, or glob.**

**The recurring failure this prevents** — one disease, six disguises, each a
plausible *wrong* answer that passed every check except an explicit one:

| # | incident | wrong… | how it bit |
|---|----------|--------|------------|
| 1 | Vesta-as-solar           | source   | reflected sunlight substituted for direct solar |
| 2 | glob-loading             | star     | globbed a tree, picked another star's frames |
| 3 | arm-registry default     | source   | a silent default resolved an arm to the wrong spectrum |
| 4 | Procyon tree held 55 Cnc | star     | `srho01cnc` frames in the "Procyon HST" tree (RYA-471) |
| 5 | α Cen folders mislabeled | star     | "α Cen A/HARPS" was 63 % α Cen B by OBJECT (RYA-301/479) |
| 6 | UVES O I 777 no sys-RV    | velocity | BERV applied, systemic RV not → −0.12 Å, χ²ᵣ 147 (RYA-478) |

Two root classes: **wrong-source/wrong-star** (attribution) and **wrong-wavelength**
(velocity frame). Both put good-looking photons in the wrong place.

### The standing rule (made executable)

`pipeline/frame_object_contract.py` is the single source of truth. Use it in every
loader:

```python
from pipeline.frame_object_contract import (
    assert_object, corrections_for_specsys, VelocityFrame, verify_line_position)

# §A — OBJECT attribution by authoritative header, never folder/filename/glob.
star = assert_object(h0['OBJECT'], expected='alpha_cen_a', context=f"{fn}: ")
#   • central canonical alias map (HD numbers, Bayer/proper names, composite
#     'SUN,FP,G2V'); case/separator-insensitive; one-to-one.
#   • refuses on ambiguity / unknown / wrong-star → raises (pass exc= for the
#     loader's own error type). A loader that can't confirm the star fails loud.

# §B — velocity frame, declared and corrected EXACTLY, then verified.
corrections_for_specsys(h0['SPECSYS'])              # raises on an unknown frame
VelocityFrame(specsys='TOPOCENT', berv_applied=True, systemic_rv_applied=False,
              wave_units='air', wave_scale='angstrom').validate()   # double-BERV trap
verify_line_position(wave_A, flux, 6562.79, tol_kms=8.0)   # BERV sign-check

# §C — wavelength unit+scale (nm/µm→Å, vacuum→air) from the CITED registry.
wave_air = WavelengthScale('SPIRou', native_unit='nm',
                           native_scale='vacuum').validate().to_air_angstrom(wave_nm)
#   • per-instrument convention is CITED to its DRS doc (never from memory);
#     unknown instrument → raises. Unit-sanity band gate catches a ×10/×10000
#     mismatch (the RYA-263 zero class). verify_vac_to_air = the scale sign-check.
```

Per-SPECSYS policy (loaders must apply exactly this, no more, no less):
`TOPOCENT` → loader **applies** BERV (UVES/CRIRES IDP); `BARYCENT`/`HELIOCEN` →
flux already barycentric, loader must **not** re-apply (HARPS S1D — the
double-BERV trap). Stellar **systemic RV** is separate from BERV: for any
rest-frame fit, shift to rest with systemic RV too (RYA-478).

**Wavelength scale (§C, RYA-264):** there is exactly **one** vacuum↔air converter
in the codebase — `pipeline/wavelength_util.py` (Birch & Downs 1994, VALD/NIST),
imported by both the RYA-426 UV path and this axis. The per-instrument unit/scale
lives in the cited `WAVELENGTH_CONVENTION` registry (verified against each DRS doc;
ESPRESSO and NIRPS are **vacuum**, not the table's guessed "air"). **Operation
order** is canonical and documented: convert unit+scale to air Å in the **observed**
frame **first** (n(λ) is evaluated at the observed wavelength), **then** apply the
§B velocity shift.

### Fail-loud (the principle shared by §A, §B, §C)

The correct response to *any* unverifiable attribution, velocity frame, or
wavelength scale is **raise/flag, never silently proceed on a default** — same
spirit as the RYA-288 broadening guard (a missing entry errors, it does not fall
back to solar). All three axes raise a common base, `FrameContractError`.

### How to apply

* **New loaders:** implement §A + §B + §C as part of the loader contract; the smoke
  test asserts OBJECT-selection (not glob), a post-correction known-line check, and
  the air-Å unit-sanity band.
* **Existing loaders:** audit against this checklist when next touched. Already
  wired: `hst_uv_loader` (target guard → central map), `uves_loader`
  (`expected_star` + `VelocityFrame.validate` + `WavelengthScale` air no-op),
  `spirou_loader` (`expected_star`; **vacuum-nm → air-Å now converted** — the
  RYA-481 OWED is closed, clearing SPIRou IR for production).
* **Authority ranking (RYA-495):** where headers are proven fallible, RV star-ID
  outranks OBJECT outranks folder — register that layer via `register_aliases`.

## Per-star output namespacing + frozen gold solar reference (RYA-469)

Every per-star pipeline product carries the star in its **path**
(`data/outputs/{star}/{star}_*`), so two stars cannot collide on a filename; the
gold-standard solar differential denominator is **frozen + versioned + immutable**
(`data/reference/solar/solar_abundances_v{N}.csv`, `CURRENT` pointer, hash-guarded).
Use `pipeline/data_namespace.py` for all of it; re-baseline the Sun only via
`scripts/promote_solar_reference.py` (bump, never overwrite). Full rule:
[`docs/design/adr_data_namespacing_and_gold_reference.md`](design/adr_data_namespacing_and_gold_reference.md).

## Artifact preservation: save-before-clean (RYA-461)

**The recurring failure this prevents:** gitignored artifacts produced inside a git
worktree — diagnostic plots, downloaded atlases, NLTE `.grd` grids, normalized-spectrum
intermediates — live **only** in that worktree. When the worktree is cleaned or removed,
they are **lost**. For the large NLTE grids that turns "we have the grid" into
"re-download 26 GB."

### The standing rule

> **Any gitignored artifact a pipeline or diagnostic run produces must be copied to the
> canonical local store as part of the run — never left only in a worktree.**

Concretely, right after a run writes a gitignored file it intends to keep (a plot, a
diagnostic CSV/JSON, a downloaded atlas, an NLTE grid, a normalized spectrum), call the
one-line helper:

```python
from pipeline.artifact_store import save_artifact
save_artifact("results/plots/solar_oi6300_diagnostic.png", kind="plots")
save_artifact("data/processed/solar_normalized.csv",        kind="data")
```

`save_artifact(path, kind)` copies the file into the store, de-duplicates by md5, and
records its provenance (source worktree/branch + date) in `ARTIFACT_MANIFEST.csv`. It
never deletes or mutates the source.

### The canonical store

Lives **outside** any worktree, in the directory that contains the worktrees
(default `~/Documents/Exoplanet Codex/`, override with `$CODEX_ARTIFACT_STORE`):

| subfolder      | holds                                                              |
|----------------|--------------------------------------------------------------------|
| `plots/`       | diagnostic plots (`*.png` / `*.pdf` / `*.svg`)                      |
| `data/`        | normalized-spectrum and other reusable data intermediates          |
| `diagnostics/` | diagnostic CSV/JSON outputs (audit / proof / verdict tables)       |
| `grids/`       | NLTE `.grd` / model grids (large, external, expensive to re-fetch)  |
| `atlases/`     | downloaded reference atlases (Kitt Peak, CALSPEC, IRTF, …)          |

`ARTIFACT_MANIFEST.csv` at the store root records every preserved file with its md5 and
provenance.

### Cleaning worktrees (the destructive half)

1. **Rescue first.** Before removing any worktree, run the save-before-clean rescue so
   every at-risk gitignored artifact is in the store. Rescue is **non-destructive** —
   it copies, never deletes.
2. **Only remove confirmed-merged worktrees.** Verify the branch's work is on `main`
   (a merge, or a patch-equivalent cherry-pick) before `git worktree remove`. When in
   doubt, **keep** the tree.
3. **Never delete an artifact that is not yet confirmed in the store.**

This rule is also the reason the `codex-artifact-preservation` step belongs in every
Mr. Code brief that produces gitignored output.

## Element status tracker must be updated with the element (RYA-594)

`data/audit/element_status_tracker.csv` is the git-tracked single source of truth for
"where does each of the 27 elements stand". It exists because per-element status has
drifted before: RYA-524's master 27×2 audit was needed precisely because status was
scattered across the verdict file, the RYA-463 registry and ticket history, and none of
them agreed.

**The standing rule:** a ticket that changes an element's

- **classification** (`WIRED-OK` / `GENUINELY-OWED` / `DONE-BUT-STALE` /
  `VINTAGE-INFLATED` / `WRONG-SPECIES` / `NLTE-VOID`), or
- **`verdict_value` / `tier`**, or
- **either engine's model-atom vintage** (`AB-INITIO` / `SCALED-DRAWIN` / `LTE` /
  `NLTE-VOID`)

**must update that row in the same change** — bump its `snapshot_date` and cite the ticket
in `source_tickets`. **A ticket that changes an element's status without updating this file
is INCOMPLETE**, the same discipline as the end-of-session Linear comment.

Two corollaries:

1. **Never silently reconcile a disagreement.** If the tracker and the live verdict
   (`docs/audit/solar_phase_c_verdict_rya371.md`) disagree, record it in
   `data/audit/element_status_tracker_drift.md` and flag it — do not pick one and
   overwrite.
2. **This file is the tracker, not a second pipeline output.** It is maintained outside
   the pipeline and reviewed by eye. A read-only one-way mirror generated *from* it is
   fine; a second writable copy is a single-source-of-truth violation.

The file carries the same rule as `#` comment lines in its own header, so it travels with
the data. Read it with `pandas.read_csv(path, comment='#')`.
