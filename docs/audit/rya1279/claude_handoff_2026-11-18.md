# Bride / RYA-1371 handoff to Claude — 2026-11-18

Ryan requested stopping Codex science work at the weekly limit, updating the PR
and leaving a complete continuation record. No new science run in this handoff.

## Start here

- Worktree: `/Users/ryanschmitt/codex/rya1371`.
- Branch: `ryandamienschmitt/rya-1371-oxygen-reference-provider`.
- PR: https://github.com/damienabraxas/exoplanetcodex/pull/620 — OPEN DRAFT.
- Latest frozen science checkpoint: `0a5f37776d656d7400c536f14ae19fd50a441bb9`.
- Main integrated for this handoff: `78ab9eb06ee00fb359611e71edb3d0e4aa5ba093`
  (PR626 / RYA-1388, includes prior main Cr landings). The commit containing
  this handoff is the continuation HEAD; use `git log -1` for its exact SHA.
- Register v299 is integration/handoff bookkeeping. Latest science remains v298.
- RYA-1371 remains In Progress / scientific HOLD. PR stays draft; Ryan merges
  and approves scientific adoption/publication. No validated solar 1D/3D claim.

Read `AGENTS.md` -> `LEDGERS.md` and its ordered ledgers; apply current science,
protocol, data-safety and development-cycle authorities. Then read this file,
`oxygen_interior_radiation_probe.md`, `oxygen_thermodynamic_contract_next.md`,
latest RYA-1371 comments and relevant code/tests. Full chronological evidence:
`independent_oxygen_reference_review.md`. Preserve old journal/source identities;
older “next” paragraphs are superseded by the latest sections, not erased.

## What is established

The existing held 5750K/logg4.5/A(O)8.70 reference state is an experimental
benchmark; canonical-solar observations remain unqualified. Full O graph:
41 levels, 256 BB, 22 BF, 56 depths. Distinguish fixed actual populations n,
reference LTE nStar, independently reconstructed rates, fresh RT and coupled SE.
Recent milestones are fixed-population diagnostics, not coupled closure.

| Milestone | Evidence / result |
|---|---|
| v294 emission/bottom profile axis | `emission_profile_v1`; both profiles retained, exact default reproduction |
| v295 full LTE shape | `full_lte_v1`; both OI/OII stages reconstruct ~1.08e-14; absolute EOS/partition normalization open |
| v296 joint inverse-input diagnostic | `joint_lte_v1`; actual n fixed, held-LW nStar shared by compiled opacity/emissivity/rates; source-forward collision directions preserved |
| v297 saved-field quadrature | `saved_field_quadrature_v1`; G16/G8 <=3.49e-9 integrates the assumed linear-J representation, not unknown true interior radiation |
| v298 actual interior probe | `interior_radiation_probe_v1`; 32 actual samples expose up to2.545323% J/interpolation difference, with all inherited J/coefficient/background values exact |

Source collision directions: CE/CH preserve downward, CI/CH0 upward, per record
before aggregation (741/302/39/24 records). The earlier all-upward algebraic
projection is only a canary. Do not substitute it as source collision physics.
Level41 is excited OII, not OIII; it has no collision edges. Invent none.

Latest grid: 19963 -> 19995, all ancestors/endpoints retained, 32 Gauss4 points
in 8 panels selected from ALL22BF/56depth/up-down rate-normalized field indicators.
Selection precedes new RT and never uses native residual/abundance agreement.
Worst new J/interpolation difference: +0.0254532324906587 at130.231096151674nm
VACUUM/top-depth0. Samples span105.6409007--130.4697080nm; they are selected
panels, not a global UV/observational validation. All32x56cells retained.

Fresh baseline reproduces all8prior arrays exactly. Both actual RT cases take
9 radiation iterations, ~6.6s excluding preflight. Max compiled rate refinement
1.99807917e-4; global balance residual barely changes2.72286197895160e-4 ->
2.72292163791942e-4. Stable global maximum does not certify local-field accuracy.
Independent all31136positive directional/depth rates match3.07532e-14 and total
chi/thermal eta4.44089e-16 under existing arithmetic screens.

Selected-panel Gauss4(actualJ-linearJ) field-response shadow: max coefficient
response2.14850265e-4 and residual2.66945671826650e-4. BB/C unchanged in that
shadow. It is not a global integral/error bound, replacement operator or fit.
All outcomes retained even when agreement worsens. No retuning/default adoption.

## Exact local assets and environments

```sh
cd /Users/ryanschmitt/codex/rya1371
git status --short
git log -1 --oneline
source config/hosts/mac.env
LW="$CODEX_BRIDE_EVAL/.venv/bin/python"
PROJECT="$CODEX_WORK_ROOT/venv/bin/python"
RAW="$CODEX_WORK_ROOT/bride_mac_checks"
"$LW" -m pytest -q tests/test_*rya1371.py
"$PROJECT" -m pytest -q tests/test_line_patch_rya1388.py
```

Held LW Python resolves to `/Users/ryanschmitt/bride_rya1118/.venv/bin/python`.
Keep environments separate: held LW lacks pandas. An initial combined test
invocation failed collection for that reason; separate reruns pass 143 Bride
in3.00s and13line-patch tests in0.44s. No packages installed/upgraded. No owned
science job left running. New integration commits use `[skip ci]`; local checks
are not hosted-CI proof. Do not start Sirius/hosted science jobs as a side effect.

Registered source/engine/build paths are in `config/path_register.yaml`; the
source engine pin remains32085d26c46e9919a986397b5514a9728abe2074. Mac paths are
resolved by `config/hosts/mac.env`. Private oracle/atomic raw assets remain local,
receipt-only; never reconstruct missing source data from metadata or redownload
without checking holdings. Preserve binaries/patch receipts/default constants.

Latest raw directories under `$RAW`:

- `rya1371-radiation-probe-plan-v1` — successful source-bound packet/generated plan.
- `rya1371-radiation-probe-baseline-v1` — fresh exact19963-node reproduction.
- `rya1371-interior-radiation-probe-v1` — actual19995-node RT/rates/coefficients.
- `rya1371-interior-radiation-probe-audit-v1` — independent full rate/coefficient audit.
- `rya1371-interior-radiation-comparison-v1` — all actual samples and shadow response.
- `rya1371-radiation-probe-preflight-failed-v1` — preserved failed pre-fix driver/note.

Compact copies and all raw receipts: `data/results/rya1371/interior_radiation_probe_v1`.
Last full verification:25834raw receipts/1181compact copies, exact length/hash/bytes.
These science artifacts are unchanged by main integration. Reproduction commands
for every stage and all failure details are in `oxygen_interior_radiation_probe.md`.
Existing output directories must never be overwritten or reused.

## Next concrete coding task — not yet implemented

Only5slots remain under the ORIGINAL20000point cap. A final nested4point probe
can retain ALL19995ancestors and finish at19999. Select its panel from actual
rate-normalized represented-field indicators over all BF/depth/directions;
never choose by desired native residual or a single favorable depth.

Current producer and loader explicitly require32points. Do NOT run that producer
on19995nodes: it must refuse the cap. Implement an explicitly versioned4point
protocol/producer/loader, with fresh registered output resource, source receipts
and negative tests. `gaussian_probe_grid` already supports4-point budgets, but
that alone does not qualify the end-to-end protocol. Preserve old reproductions:
changing pinned helper files intentionally invalidates old v1 receipts. Use
frozen versions or new source-bound controls; never edit receipts or weaken checks.

Required before another RT: unit tests for budget/retained nodes/state/physics/
changed-producer/overwrite refusal; fresh baseline reproduction; same one-thread,
600s,20000-point guards and native2GBcontrols. No silent tiling or point removal.
Then new actual J vs interpolation, all rate/coefficient audits, all samples,
receipt verification, self-review and register/Linear handoff. Even success is a
local numerical probe, not global convergence or physical adoption.

Beyond this cap a new mesh/compute strategy requires explicit scope and its own
normalization/coupling validation. Absolute oxygen/EOS/partition qualification
can proceed independently. Current finite-level LTE shape holds the measured
total; external partition normalization (LEVEL_PARTF=F), chemistry/EOS and full
thermal consistency are still unqualified. Do not silently replace totals.

Coupled action/preconditioner/conservation/convergence/restart, BB/finite-CRD
numerical convergence and canonical-solar observed AIR/intensity-vs-flux/
continuum/telluric/instrument/fullband/3D gates remain OPEN. Current coupled runner
refuses joint/emission/alpha controls. Fe is next after Oxygen foundation;
Fe atom/resource/FeII gaps persist, C/N queued, RYA-1374 gated. No gold, original
spectra/FITS, adopted abundance, public product, raw atomic table or source default
change without its existing scientific/user gate.

## Integration and debt Claude must not mistake for new science

Main has been integrated into this feature branch to clear PR conflicts; the PR
has NOT been merged into main. All1113main-only changed paths are retained exactly (three overlapping
metadata files resolved explicitly), and solar production feeds match main.
The pre-existing branch Frontier coverage source-pin difference is retained
unchanged; it is not a new handoff science/product change. Fe feed1.612 and ORCH subset-patch code belong to
RYA-1388/PR626, not an experimental Bride rerun. No production outputs regenerated.

Register conflict: retain branch scientific v253..v298 and main's distinct Cr
v253 landing as a main-scoped historical entry; v299 records integration. No
scientific history deleted or renumbered, no new scientific gate signed off.

Bibliography conflict: preserve distinct Gurell preprint/final-paper entries and
all main rows. Ward same-DOI/same-file intake alias `ward2023_crii` is consolidated
under current main key `ward2023_cr2`; old record is copied below and preserved in
Git0a5f3777. No active code/result references to old key were found; an old generated
sources-page anchor still records that history. No site regeneration/publication.
Main's adopted-data role remains intact; Raassen author-initial variants remain
explicitly unresolved in notes, no author verification invented.

Current bibliography audit stops before library inventory on six unsupported
`verified` values imported unchanged from main: bergemann_cescutti2010,
blackwell1984_cr1, nilsson2006_cr2, pinnington1993_cr2, ward2023_cr2, tozzi1985_cr1.
Do not promote those labels without source verification. Last completed inventory
at v298:215documents/132indexed/83unindexed/12missing local paths. Global path
checker had245present/18pre-existing broken. No libraries acquired here.

Integrated result-generator errors:5700. A fresh metadata-only current-main
snapshot/index baseline gives the EXACT same5700diagnostics, zero added/removed.
The increase from old Bride5120 is580inherited main diagnostics. Baseline only
materialized manifest/generator sources, never scientific output or physics.
Logs: `/private/tmp/rya1371-handoff-generators.log` and
`/private/tmp/rya1371-handoff-main-generator-baseline.log`.
Existing tracker/gold regeneration notice and historical immutable native-log
whitespace remain. No global-green or hosted-CI claim. Do not trim immutable logs
or manufacture generator registrations to silence these debts.

Run-bug ledger113rows; latest fixed output-resolver defect uses codex_path ONLY
for fresh write destination and require_codex_path for READ inputs; fixed at
4f35bbb1. Regression tests cover absent/existing/wrong-location destinations.
Old failed attempt preserved; no native/RT call on its failure. No new science bug
claimed in this documentation/main-integration handoff.

## Historical Ward intake record (preserved, not re-verified here)

```json
{
  "key": "ward2023_crii",
  "authors": "Ward, J. W.; Li, J. J.; Schwartz, J.; Nave, G.; Raassen, T. A. J. J.; Uylings, P. H. M.",
  "year": "2023",
  "title": "Branching Fractions and Transition Probabilities for UV Transitions in the Spectrum of Cr II",
  "venue": "The Astrophysical Journal 959, 8",
  "doi": "10.3847/1538-4357/acfafe",
  "url": "https://doi.org/10.3847/1538-4357/acfafe",
  "category": "atomic_data_gf",
  "role_in_codex": "Reference-library metadata intake only: candidate laboratory Cr II source. No numerical table audit, canonical_gf ingestion, grade, abundance or Oxygen-model adoption.",
  "tickets": "RYA-1371",
  "local_file": "Reference documents/Ward_2023_ApJ_959_8.pdf",
  "license_note": "Printed volume date 10 December 2023 and online publication 28 November 2023; both citation year2023. DOI/title/authors read from held PDF first page. First two PDF pages inspected during the standing library refresh on2026-11-09. No redistribution/data-reuse rights inferred; numerical tables remain unqualified. SHA256 c8124b2a7a4e03e7b3a2bb61c3e31b62ba754a27fef074d24dd9ffb4c4f1b05d.",
  "verified": "extracted"
}
```

## Handoff verification

No unmerged index entries. Register pointer, mirror, reconciler and freshness
pass. Resolved/new handoff files pass whitespace checks; immutable historical
log spacing is retained. All frozen RYA-1371 result files remain byte-identical
to0a5f3777. No science jobs or owned tool sessions remain after closure.
