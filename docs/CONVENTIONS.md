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

**Since RYA-654 it is GENERATED — do not hand-edit it.** The hand-maintained era ended
because hand-maintenance is what let it drift: it was still carrying Co's demoted
blue-edge +1.188 and N's cleared NLTE debt long after the verdict channel had moved.
This is RYA-436's rule applied to the results side — *generate, do not hand-sync*.

```
python scripts/generate_element_status_tracker_rya654.py            # rebuild
python scripts/generate_element_status_tracker_rya654.py --check    # CI: committed == regenerated
```

**The standing rule is now about WHERE you make the change**, not whether you remember to:

- a **status** column (`verdict` / `verdict_value` / `sigma` / `n_lines` /
  `delta_vs_asplund` / `method`) is derived from the phase_c verdict channel — it is
  ratified canonical for status (RYA-654 §1). Fix it at the source and regenerate;
- `tier` is the RYA-522 ratified freeze tier, and `ion` / `regime_verdict` come from
  `config/physics_regime_rya400.yaml`. Same rule: fix the source;
- the **analyst** columns — **classification** (`WIRED-OK` / `GENUINELY-OWED` /
  `DONE-BUT-STALE` / `VINTAGE-INFLATED` / `WRONG-SPECIES` / `NLTE-VOID`),
  `action_needed`, `source_tickets`, **either engine's model-atom vintage**
  (`AB-INITIO` / `SCALED-DRAWIN` / `LTE` / `NLTE-VOID`) and `notes` — have no machine
  source. Edit **`data/audit/element_status_tracker_editorial.yaml`**, bump that entry's
  `editorial_updated`, cite the ticket in `source_tickets`, and regenerate.

**A ticket that changes an element's status without updating the relevant source and
regenerating is INCOMPLETE**, the same discipline as the end-of-session Linear comment.
`--check` is what enforces it: it diffs the committed CSV against a fresh regeneration,
so a hand-edit is a build break rather than a silent divergence.

Note `verdict` and `tier` are **different things** and have separate columns: a verdict is
live status, a tier is a freeze decision. Before RYA-654 one column carried both
vocabularies and the guard had to infer.

Three corollaries:

1. **Never silently reconcile a disagreement.** Cross-artifact contradictions are checked
   by `python -m pipeline.ledger_consistency_guard` (RYA-632). A disagreement is either
   fixed at its source or **ratified** in `data/audit/known_verdict_divergences.yaml`
   with a ticket — never overwritten, and never annotated just to reach green.
2. **This file is the tracker, not a second pipeline output.** A read-only one-way mirror
   generated *from* it is fine; a second writable copy is a single-source-of-truth
   violation.
3. **The generator will not run off an uncommitted verdict artifact.** The tracker traces
   to the ratified, committed phase_c run, never to a local Mac canary re-run.

The file carries the same rule as `#` comment lines in its own header, so it travels with
the data. Read it with `pandas.read_csv(path, comment='#')`.

## Isotope fractions live in the ENGINE or in the gf — never both (RYA-684)

Turbospectrum multiplies an isotope-coded species' population by `isotopfrac(Z, A)`
before forming the line opacity (`bsyn.f:1350`), and `makeabund.f` sets
`isotopfrac(Z, 0) == 1.0` precisely so a list can opt out — its own comment says the
zero code is "used in the case of no isotopes wanted in calculation … if for example
isotopic factors were included in gf-values".

So a Turbospectrum line list must commit to exactly one of two forms:

* **(A) isotope-coded `Z.AAA` + fraction-free log gf.** Each isotope block carries the
  FULL oscillator strength and the engine applies the fraction. The GES v6 HFS/ISO list
  (`GESv6_atom_hfs_iso.420_920nm`) and the Gerber NLTE deck are form (A).
* **(B) uncoded `Z.000` + folded log gf.** The split is already in the gf and the engine
  applies 1.0. `scripts/rya581_ba2_deblend_sirius.py` is form (B) — it writes its own
  Ba II 5853 HFS block as `56.000`.

**Isotope-coded AND folded is the error.** The fraction lands twice, the feature comes out
`sum_i f_i^2` too weak, and the fitted abundance absorbs `-log10(sum_i f_i^2)`: +0.3002 dex
for Eu II, +0.2694 for Ba II, +0.2415 for Cu I. The shipped TSFitPy `linelist_vald`
"for-grid" lists are form-(B) gf values written with form-(A) headers, which is what
RYA-684 measured and RYA-565 saw as a +0.300 VALD-vs-GES leg offset.

Two things follow, and both are guarded:

1. **HFS is not isotope structure.** Mn, Co, Sc and V are hyperfine-split but effectively
   mono-isotopic, so `sum f^2 == 1` and they are structurally immune. Do not reason from
   "this element has HFS" to "this element is exposed".
2. **A harness must never fit a TARGET species against a folded, isotope-coded block.**
   Call `pipeline.isotope_gf_convention.assert_target_convention(linelist, Z, ion)`
   immediately before `bsyn`. It reads which species are folded from the committed
   RYA-684 audit record rather than a hardcoded list, so re-vendoring a line list and
   re-running `scripts/rya684_isotope_gf_audit.py` keeps the guard honest. Exposed blocks
   in the BLEND model are recorded, not fatal — RYA-684 measured those at <0.01 % of window
   absorption in every window feeding a live value.
## A result artifact must not land without its generating harness (RYA-686)

**Every file committed under `data/results/` must be accompanied by the committed code
that generates it, and the two must be linked in `data/results/GENERATORS.yaml`.**

**The recurring failure this prevents:** Sirius computes, only the result JSON comes
back, and the harness that produced it stays in a scratch directory and is never
committed. RYA-567 makes compute-on-Sirius the rule; this is the missing preservation
half. Two verified instances, both of which cost real time in the RYA-672 batch:

- **The harness that never existed.** RYA-559's merge commit `564824a` shipped
  `data/results/solar_ba_synthesis_rya559.json`, a phase_c hook and a test — but not the
  synthesis harness. `scripts/rya559_ba2_synth_sirius.py` has never been committed on any
  branch, and is not on Sirius either. **A(Ba) = 2.410 is not independently checkable**,
  and RYA-581 had to rebuild the harness from scratch off the RYA-551 pattern plus a
  commit message.
- **The harness that shipped incomplete.** `scripts/rya560_zr2_synth_sirius.py` *is*
  committed, but declares no `argparse` at all — so the `--deblend` entrypoint RYA-585's
  brief called for did not exist and had to be added.

The cost compounds: every re-measurement of the same species pays the rebuild again, and
an unreproducible measurement cannot be independently checked at all.

### The standing rule

> **A result artifact and the code that produced it land in the same PR. If a
> measurement genuinely cannot be reproduced, say so — never point an artifact at a
> script that did not produce it.**

### The linkage: a manifest, and why not the alternatives

`data/results/GENERATORS.yaml` maps each artifact to its generator. Three mechanisms were
weighed; the manifest wins on four counts the others cannot cover.

| | naming convention | `_meta.generator` in the JSON | **manifest** |
|---|---|---|---|
| covers `.csv` / `.txt` (16 of 34 artifacts) | yes | **no** | yes |
| mutates historical artifacts | no | **yes** | no |
| can record "no generator exists" | **no** | yes | yes |
| carries the *invocation* | **no** | awkward | yes |
| can drift | no | yes | yes — **but the guard makes drift a build break** |

- **A naming convention** (`*_ryaNNN.*` → `scripts/*ryaNNN*.py`) is the simplest, and it
  *would* have caught the Ba case. It was rejected because it **proves the wrong thing**:
  RYA-485 shipped two artifacts and two scripts, so a shared ticket token cannot say which
  script generates which artifact — it would report green if either one went missing. It
  also cannot classify the seven artifacts carrying no ticket token at all
  (`procyon_co_*.csv`, `procyon_fe_spread*.csv`, `procyon_uves_oi777_phase2.*`); nor
  `rejection_ledger_solar_rya429.json`, whose real generator is `pipeline/lines_fit.py`;
  nor `zr2_deblend_rya585.json`, generated by `scripts/rya560_zr2_synth_sirius.py --deblend`
  — a *different ticket's* script, which is the normal shape of a follow-on re-measurement.
  It would need an exemption file for nine of thirty-four — i.e. a manifest, arrived at by
  a longer road.
- **A `_meta.generator` field** is genuinely self-describing and only three result JSONs
  carry a `_meta` block today, so it is not the established shape the count suggests. It
  cannot serve the 16 CSV/TXT artifacts without a second mechanism, and adding a top-level
  key to a bare wavelength-keyed dict like `sr2_synthesis_rya551.json` changes what every
  `for w, d in obj.items()` consumer iterates over.
- **The manifest's** one real cost is the "another thing to keep in sync" objection. That
  is answered by making the guard **bidirectional**: an artifact with no entry fails, an
  entry naming a missing artifact fails, and an entry naming a missing generator fails.
  Sync is not a discipline anyone has to remember — it is a build break, the same shape as
  RYA-654's `--check` on the element status tracker.

### The three statuses

The audit found three genuinely different things sitting in `data/results/`. Collapsing
them would lie, so the manifest names all three:

| status | means | requires |
|---|---|---|
| `COMMITTED` (default) | machine output, harness in the repo | `generator:` that exists |
| `HAND_AUTHORED` | a human wrote it — a decision record, not a program's output | `generator: null`, `sources:`, `note:` |
| `UNREPRODUCIBLE` | machine output whose harness was **never committed** | `generator: null`, `note:` |

`UNREPRODUCIBLE` is **the honest record of a defect, not an escape hatch.** Its membership
is frozen in `tests/test_result_generators_rya686.py`, so adding one is a deliberate,
reviewed edit rather than a quiet green.

### Running the guard

```
python scripts/check_result_generators.py                      # the landing gate
python scripts/check_result_generators.py --check-invocations  # + the RYA-560 check
```

Both run in CI inside the `test` job.

### What the invocation check does and does not catch

Recording `invocation:` lets the guard AST-read the generator's `argparse` and verify every
long flag in the recorded command is actually declared. That **does** catch the RYA-560
class — "committed, but not invocable as documented" — at the point where it matters, which
is when a brief names a flag.

It does **not** catch a flag that is declared but does nothing, a harness that runs and
produces different numbers, or a harness that no longer matches the artifact it once wrote.
Those remain review's job. The check is deliberately static: it never imports or executes a
harness, so it cannot be defeated by an import-time failure and cannot itself run compute.

## Two-engine driver inputs: one generated, seven committed (RYA-682)

`scripts/rya527_two_engine_run.py` is the Beta gate's driver. Its inputs split cleanly,
and the split decides how a missing one is repaired:

* **GENERATED — one.** `data/outputs/{star}/{star}_per_line_synth_v2.csv`, the Engine-B
  synthesis-v2 per-line table. `data/outputs/` is **gitignored**, so this is regenerable
  by design and never committed. A clean checkout will not have it, and that is correct.
  Regenerate on Sirius (RYA-567 — computation is Sirius-only):

  ```
  python -m pipeline.abundances_derive solar ATLAS9.Castelli synthesis-v2 --pin
  ```

* **COMMITTED — seven.** The dedicated Engine-B synthesis-harness measurements for the
  synthesis-required elements (CNO cross-arm, Mn RYA-473, Cu/V RYA-466, Sr II RYA-551,
  Zr II RYA-560 + RYA-585, Mg 5528 RYA-592). These are tracked; a missing one is a broken
  checkout, so the repair is `git checkout -- <path>`, never re-running a harness.

**Generate synthesis products on `venv312`, not `venv_ci`.** These are different
environments on purpose: `venv312` is the RYA-517 **reference** venv (py3.12.13 / numpy
**2.2.6**, exact pins, do not install into it), while `venv_ci` is built from
`requirements.txt`, which floors `numpy>=1.26.0` with **no ceiling** and has floated to
numpy **2.5.1**. The standing "use `venv_ci`, never `venv312`" rule is about running the
**test suite** — it does not apply to generating science products, and following it there
silently produces an empty artifact:

> `ispec/abundances.py:132` assigns a size-1 array into a scalar recarray slot. NumPy
> deprecated that in 1.25 and made it an **error in 2.3**. On numpy ≥ 2.3 every element
> loses its atom code, every line is written `status='failed'`, and the run **exits 0**
> with a full-length per-line table in which no row is usable. The frame is not empty, so
> the RYA-342 empty-set guard does not fire.

This is the same class as the RYA-313 `np.trapz` finding: an unpinned dependency floating
past what a vendored engine tolerates, invisible until something runs on the declared
stack. Three guards now make it loud instead of silent, all in
`pipeline/two_engine_inputs.py`:

1. `assert_engine_b_artifact()` — the generated input exists **and has ≥1 usable row**.
2. `assert_committed_inputs()` — every tracked input is present, so the driver can never
   emit a quietly smaller record set.
3. `assert_synthesis_stack()` — the running interpreter can actually build iSpec atom
   codes, naming the numpy cause and the RYA-517 remedy.

The driver runs (1) and (2) **before any compute**. Previously the check sat after the
whole Engine-A leg, so a clean checkout paid minutes of GES-linelist load, EW triage and
MOOG baseline to be told its input was missing; it now fails in ~2 s.
`_run_synthesis_v2_mode` additionally **refuses to write** a per-line table with zero
usable rows — never emit a canonical input that looks like a successful run.
