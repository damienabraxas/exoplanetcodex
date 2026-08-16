# SCIENCE STANDARDS — authoritative decision records

Single source of truth for cross-cutting science-method decisions. New file (RYA-505);
flag for the RYA-179 docs-sync revision. Add sections as decisions are recorded; each
section cites its evidence and the ticket that set it.

---

## Corroboration-accept — registering an ab-initio value past a failed acceptance gate — RYA-544/545

**Decision date:** 2026-07-13 · **Set by:** Ryan (RYA-544 decision comment) · **First invocation:** Ti I (RYA-545).

A NLTE (or other physics) value MAY be registered **despite a FAILED acceptance gate** — but ONLY as a
**bounded exception** when **all five** of the following hold. This is *not* "accept when the gate is
inconvenient": a gate that fails because the point estimate is off does **not** qualify.

1. **Derived ab-initio, not tuned** to the gate — the validate-don't-tune firewall is intact (the value
   was produced by the physics, never fitted to pass).
2. **Reproduces an external published** result (independent of us).
3. **≥ 2 independent instruments** corroborate the gate's *point estimate*.
4. The gate failed on **precision** (wide SEM / measurement floor), **not** on the point estimate being off.
5. The **incumbent it replaces is confirmed wrong** (independently established).

The failed gate is then recorded as a **documented diagnostic** (why it is precision-limited), not erased and
not counted as a value defect. The threshold itself is **never lowered** to manufacture a pass.

**Directionality caveat (RYA-546):** this rule licenses *adopting* an externally-corroborated ab-initio value.
It does **NOT** license *retiring* a correction (e.g. dropping a NLTE grid to LTE) on a weak/precision-limited
balance — that is the opposite direction and needs its own affirmative evidence.

**First invocation — Ti I (RYA-545):** registered on the Mallinson-2024 ab-initio grid (δ=+0.0506, PySME/MARCS)
though the ionization-balance gate STOPPED on precision (SEM ~0.12–0.16 on the thin solar Ti pools, 6–7 Ti I /
3 Ti II). All five held: (1) derived from the grid, firewall intact; (2) reproduces Mallinson-2024 +0.052;
(3) two blend-aware instruments (PySME full-window blended + production TS synth-EW) agree LTE balance ≈0,
NLTE ≈+0.05; (4) the gate failed on SEM, not on the estimate (LTE re-graded to ≈0, NLTE dead on the grid);
(5) the incumbent (Bergemann-2011 scaled-Drawin +0.108) is confirmed inflated ~2× by the RYA-546 vintage audit.

---

## Hot-Teff NLTE grid coverage (F-star benchmarks: Procyon 6554 K, τ Boo 6400 K) — RYA-505

**Decision date:** 2026-07-02 · **Evidence:** `scripts/nlte_fstar_ceiling_rya505.py`
→ `data/results/nlte_fstar_ceiling_rya505.csv`. **Scope:** hot-benchmark fidelity upgrade;
NOT on the α Cen / 55 Cnc critical path (those are 5200–5800 K, already in-grid).

### Headline
The RYA-349 "6200–6500 K wall" is almost entirely an **under-loading artifact, not a grid
limit.** Our on-disk non-Fe grids are the benchmark-node **subsets** synthesised in RYA-410
(cool FGK dwarfs + 55 Cnc), not the full published grids. **9 of 12 non-Fe elements have a
published grid that already covers Procyon** — the fix is to re-load/re-synthesise the
Procyon node in the *same* family (self-consistent), not to mix codes. Only **Ba** is a
genuine real-limit clamp.

### Grid-selection hierarchy (mandatory order)
1. **Self-consistent extend** — re-synthesise/re-load the hot node from the *same* code +
   model-atom + model-atmosphere family we already wire. No new systematic. **Preferred.**
2. **MPIA-with-cross-check** — only where no self-consistent extension exists. MPIA/Bergemann
   grids are DETAIL/SIU on MAFAGS-OS (≠ our Amarsi/PySME + MARCS), so require a
   validate-don't-tune cross-check at overlap (the Sun + one cool star in both grids): the
   two grids' corrections must agree within **±0.05 dex** at overlap before the MPIA grid is
   trusted at hot Teff. Disagreement is a recorded finding, never a silent choice.
3. **Bounded clamp** — only for elements with no higher-ceiling grid at all: a Teff clamp at
   Procyon's **54 K** overshoot ONLY, monotonic-in-Teff elements only, with the flagged
   systematic propagated into the σ budget (RYA-282). Never extrapolate a 3D-NLTE correction
   past its grid; never a silent clamp (Fe/RYA-319 precedent).

### Per-element hot-Teff coverage map
| element | current family | on-disk ceil | published ceil (cited) | verdict |
|---|---|---:|---:|---|
| Na | Amarsi2020 GALAH · PySME · MARCS | 6200 | **8000** (Zenodo 3982506 v3 prov, 2500–8000) | self-consistent-extend |
| Mg | Amarsi2020 GALAH · PySME · MARCS | 6200 | **8000** (same release) | self-consistent-extend |
| Si | Amarsi2020 GALAH · PySME · MARCS | 6200 | **8000** | self-consistent-extend (near-LTE, low pri) |
| Al | Amarsi2020 GALAH · PySME · MARCS | 6200 | **8000** (v3 prov) | self-consistent-extend (near-LTE, low pri) |
| K  | Amarsi2020 GALAH · PySME · MARCS | 6200 | **8000** (v3 prov) | self-consistent-extend |
| S  | Amarsi2025 · PySME · MARCS | 6200 | **8000** (prov, 3000–8000) | self-consistent-extend |
| Ca | Mashonkina2017 · DETAIL · MAFAGS-OS | 6500 | unrecorded (subset) | **self-consistent-extend via Amarsi Ca `.grd`** (staged in `amarsi_galah/`, Family-A) with validate-don't-tune cross-check vs the MAFAGS-OS value |
| Mn | Bergemann MPIA · DETAIL/SIU · MAFAGS-OS | 6500 | 7000 (nlte.mpia.de survey) | **self-consistent-extend via Amarsi Mn `.grd`** (staged) w/ cross-check |
| Ti | **Mallinson2024 · PySME · MARCS** (RYA-545) | 8000 | **8000** (Zenodo 10753497, hull 2500–8000) | **self-consistent-extend** (ab-initio Grumer-Barklem-2020; supersedes Bergemann2011 MPIA, RYA-546; near-LTE, low pri) |
| Cr | Bergemann2010 MPIA · MAFAGS-OS | 6500 | 7000 (survey) | MPIA-with-cross-check |
| **Ba** | Korotin2015 · MULTI · MARCS | 6500 | **6500 (REAL LIMIT)** — Korotin2015 prov coverage 4000–6500 | **bounded clamp** (54 K, monotonic) or find a higher-ceiling Ba grid |
| Sr | Bergemann2012 INSPECT · MARCS (metal-poor) | 6000 | unrecorded | defer (off-critical; INASAN primary pending RYA-433) |

Priority order for execution (NLTE-significant at F-star Teff first): **Ca, Ba, Na, Mn**,
then Mg/K/S; Si/Al/Ti near-LTE and low priority.

### CNO (Step 1) and banked Procyon O (Step 2) — already correct, confirmed
- **C/O NLTE:** `pipeline/nlte_cno.py` already splits legs at `TEFF_3D_CEILING = 6500 K`:
  3D leg (tables 2/3) ≤ 6500 K; **1D-NLTE leg (tables 5/6) above** → **Procyon 6554 K uses
  the 1D-NLTE leg, which spans Teff 4000–8000 K (17 nodes) → in-grid.** No change (RYA-359).
  3D refinement stays future v2 (RYA-444/445), matching the IR v1/v2 boundary (RYA-504).
- **RYA-483 banked Procyon O:** primary indicator `OI_777_1D_NLTE`, A(O)=8.82, [O/H] +0.085
  — **used the 1D-NLTE correction, in-grid at 6554 K → NLTE-valid.** Its PROVISIONAL status
  is the cross-instrument continuum zero-point (~0.18 dex) + the terminal [O I] 6300 leg,
  **not** an off-grid NLTE correction. No re-derivation owed.

### Coverage after the decision
- **Procyon 6554 K:** in-grid-once-loaded for Na/Mg/Si/Al/K/S (+ Ti/Cr/Mn where the family
  extends); Ca/Mn via the staged Amarsi `.grd`; **Ba clamp-flagged** (only true edge case);
  Sr deferred; **Fe already NLTE-live** (RYA-319). C/O 1D-NLTE in-grid.
- **τ Boo 6400 K:** inside every ≥6500 grid → fully in the published grid for all except the
  Sr 6000 metal-poor case.

### Execution status (Steps 3–5) — recipe recorded, physical re-synthesis pending
The self-consistent extension is a re-run of the existing machinery — `pipeline/pysme_nlte.py`
+ `scripts/resource_clamped_grids_rya409.py` — adding the Procyon node (6554/4.00/+0.03) to
each Amarsi element's node list, from the **same `.grd` binary** (gitignored/freed; refetch
from Zenodo 3982506, md5 per each `*_amarsi2020_v3.prov.json`, e.g. Mg `cdc4449e…`). No
cross-code mixing for the Amarsi set. Ca/Mn: synth from their staged Amarsi `.grd` and run
the ±0.05 solar cross-check vs the current MAFAGS-OS registered value before swapping. Ti/Cr:
re-scrape MPIA to 7000 K + solar overlap check. Ba: implement the bounded-clamp-with-σ. This
session delivered the recon gate + decision record + the CNO/O confirmations; the grid
re-synthesis is the bounded next execution step (needs the multi-GB `.grd` refetch + PySME).

---

## Reference compute stack: Python 3.12 + numpy 2.2.x — RYA-517

**Decision date:** 2026-07-04 · **Evidence:** `requirements.lock`,
`data/audit/rya517_compat_probe/{findings.md, mac_v2_revalidation.md}`. **Scope:** the pinned
interpreter + scientific-library stack the solar anchor is validated on; binds BOTH machines.

### Decision
The reference stack is **Python 3.12.13 + numpy 2.2.x** (scipy ≥1.15.3, astropy ≥7.1.0,
pandas ≥2.2.3; exact resolved versions in `requirements.lock`). A stack change is a **science
event** — the anchor is re-validated on the new stack, never a silent infra bump.

### Why (the below-floor finding)
The banked v1 anchor ran on **py3.9.6 + numpy 1.26.4 — *below* iSpec's own declared minimum**
(`numpy>=2.2.5, scipy>=1.15.3, astropy>=7.1.0, pandas>=2.2.3`), an EOL accident (py3.9 EOL
Oct 2025). py3.12+numpy2.2 is the newest stack with a **proven iSpec C-extension build** and is
the FIRST anchor-validated stack that meets iSpec's floor.

### Validation (the science event, done — RYA-517 steps 2–4)
Solar FULL re-run reproduces banked v1 **to 0.000 dex on all 27 verdict elements**, Fe scatter
gate PASS (0.1390 ≤ 0.1398), RYA-371 verdict PASS=5/NLTE-OWED=1/CURATION-OWED=20/DATA-GAP=0 —
**identical on BOTH machines**: Mac (darwin arm64, prebuilt iSpec `.so`) and Sirius (linux
x86_64, iSpec C ext + MOOGSILENT **compiled from source**, gcc 15 needs
`CFLAGS="-std=gnu17 -fcommon -w"`). v2 ≡ v1 on both → **no re-bank** (Ryan's call: re-bank is
not automatic; here nothing moved).

### Standing rules
1. **Anchor-identity, not byte-identity.** numpy-2.x + platform libm/BLAS give three distinct
   `solar_normalized.csv`/`solar_ew.csv` md5s (py3.9, mac-py3.12, sirius-py3.12) for one
   identical set of abundances. The true floor, measured on the UNROUNDED per-line
   `normal_abund` at MOOG's native 1e-3 dex granularity (the finest quantity the pipeline
   emits): **max|Δ| = 0 dex EXACTLY on both splits** — same-machine v1↔v2 (numpy/stack axis)
   and cross-machine Mac↔Sirius (platform axis), 0 of 473 lines differ in either. The
   byte-differing intermediates prove sub-quantum FP noise exists upstream, but it never
   reaches even the finest reported abundance digit. Gate cross-machine/stack agreement on
   A(X)-identity, never on byte-identical intermediates.
2. **Full committed input set required for parity.** `data/processed/solar_ew_ges_reference.csv`
   is a *git-tracked* input (the 62-line GES Fe I pool); its absence makes the run *silently*
   fall back to a different Fe I pool (n 62→19, +0.12 dex, gate FAIL). Stage all tracked inputs
   before a cross-machine run.
3. **Native engines compile per-platform** (iSpec `.so`, MOOGSILENT) — see `requirements.lock`
   for the gcc-15 flags. Turbospectrum/standalone-SPECTRUM binaries are not on the solar
   EW→abundance path.

---

## Two-engine floor — quality-based per-line selection — RYA-525

**Decision date:** 2026-07-05 (ratified), pre-declared criterion ratified 2026-07-10.
**Evidence / build:** `pipeline/engine_selection.py`, `config.constants.TWO_ENGINE`,
`tests/test_two_engine_floor_rya525.py`; coverage map `data/curation/nlte_two_engine_coverage.csv`
(RYA-526). **Scope:** governing law for every element on every star. **The full 27-element
both-engine RUN is gated on RYA-526 grid coverage and executes under RYA-527** (this ticket builds
the law + selector + guards).

### The law
1. **Both engines, every element.** Every element is measured on BOTH **Engine-A** (1D-NLTE =
   EW + grid delta) AND **Engine-B** (synthesis; Turbospectrum LTE + TS-native Gerber NLTE).
   Neither optional. A missing synthesis grid is a grid-**acquisition** task (RYA-526/540),
   **never** a license for EW-1D-only. "Both ran" is enforced even when one is reported.
2. **Report the best, selected by a PRE-DECLARED PHYSICAL criterion — per LINE.** The reported
   value is the quality-selected best, chosen line-by-line on line/measurement QUALITY (σ, REW/
   saturation, blend flag [RYA-463], COG regime), then aggregated to the element. The criterion
   keys on the LINE, **decided before seeing which answer is closer to any reference.**
   FORBIDDEN: selecting the engine closer to Asplund/a literature value — that is tuning wearing a
   selection label (the RYA-161 firewall). The selector has **no reference-value input by
   construction**; a smuggled one RAISES (`ReferenceProximityError`).
3. **The rejected engine is recorded + shown, EXCLUDED from the value and its uncertainty budget**
   — per LINE. Never average a rejected engine into σ (the C=10.26 disease). The element budget
   carries only the winning lines' real uncertainties (inverse-variance combined).
4. **Cross-engine spread = a SEPARATE diagnostic**, never folded into the reported error bar.

### The pre-declared per-line criterion (ratified — encode, don't tune)
- **Validity gates (clause 1):** Engine-A eligible iff its grid is in-hull; Engine-B eligible iff
  `med_red_chi2 ≤ TWO_ENGINE['synth_chi2_gate']` — an eligibility floor ("synth didn't
  catastrophically fail"), **NOT** a quality selector.
- **Clause 2** — exactly one eligible → report it.
- **Clause 3 (clear regimes)** — both eligible → line regime decides: **CLEAN-WEAK** (unsaturated
  below the `saturation_knee_mA` AND unblended AND not a problem-child/HFS) → **Engine-A** (cleanest
  for clean weak lines); **HARD** (blended OR saturated OR problem-child/HFS) → **Engine-B**
  (synthesis handles the blend/saturation an EW cannot).
- **Clause 4 (INDETERMINATE regime ONLY)** — a line that is neither clearly clean-weak nor clearly
  hard → lower line-scatter σ; **exact tie → 1D-NLTE**, the differential zero-point / anchor scale
  (solar Fe is anchored 1D-NLTE 7.516), so ties stay on ONE consistent scale. Clause 3 governs
  every clear regime; clause 4 governs only the border — they do not overlap.
- **Clause 5** — neither eligible → no value; a cited disposition is recorded (never a silent PASS).

### Element aggregation + the cross-engine-mixing guard
- Reported element value = **inverse-variance combine of the per-LINE winners** (each line
  contributes only its winning engine's value + error; the rejected engine stays diagnostic).
- **`CROSS_ENGINE_MIX_GATE` (the Ti lesson):** when an element's per-line winners span BOTH engines
  AND the mean cross-engine Δ exceeds `TWO_ENGINE['cross_engine_mix_gate']` (dex), **FLAG +
  adjudicate** — do not silently average two disagreeing scales (mixing a 1D-NLTE and a synthesis
  scale that systematically disagree injects a regime-correlated bias, e.g. Ti's ~0.11 same-atom
  split, RYA-535/542). A threshold breach is a recorded adjudication, not a silent mean.

### Loud-fail guards (RYA-525 §3; siblings of RYA-409/518)
- **Missing synthesis grid → RAISE** ("acquire it / RYA-526"): read the RYA-526 two-engine ledger
  as the pre-declared exception list — raise where `disposition ∈ {acquire-task, build-task}` and
  the grid is genuinely absent, or where `wired-both` produced only one engine at runtime. Never on
  `wired-one` / `LTE-only-by-design` (documented, owned single-engine states).
- **A reported single-engine value with no cross-engine record and no documented disposition →
  RAISE.**
- **Selection keyed on reference-proximity → RAISE** (the tuning firewall).
- **No silent LTE:** the `abundances_derive.run()` `except Exception → print → continue` swallows
  around the NLTE calls are removed — a missing wired grid / NLTE failure propagates.

### Standing rules
- All thresholds live in `config.constants.TWO_ENGINE` (no inline knobs); the saturation knee and
  synth eligibility floor reference their existing SSOT homes.
- Ti stays **CHECK / excluded from the reported value** until RYA-535/542 resolves the same-atom
  systematic — but both engines still run for Ti and the spread is recorded. Mn carries the same
  open provenance question. Neither blocks the floor.

---

## All computation on Sirius — grid/atmosphere/engine provenance is single-sourced, absence → loud-fail — RYA-567

Successor to the RYA-555 calibration-ladder / RYA-528 truth-sync standards work, and a sibling of
the RYA-409 out-of-hull "no silent LTE" guard applied to the **data-provenance** axis.

### The law
**All heavy computation runs on Sirius, and every GRID / MODEL ATMOSPHERE / departure-coefficient
grid / synthesis ENGINE input for a compute step is single-sourced from the Sirius data root
(`/mnt/codex-data`, env `SIRIUS_DATA_ROOT`) — NEVER a local-Mac copy.** If a required compute-input
grid is absent under the Sirius root, the code **RAISES** naming the expected path; it never
substitutes a repo-local (`data/nlte_grids/…`) or `~/`-relative copy. A silent local-grid read for a
computation is a correctness-**and-provenance** defect of the same class as the silent-LTE fallback
(RYA-409): a repo-local copy can be stale or partial, so the number it produces is untrustworthy and
unattributable. (Set 2026-07-18 after a local-Mac grid was caught being read in place of the
authoritative Sirius-staged grid — the exact failure this standard closes.)

### The two categories (the distinction the resolver encodes)
- **COMPUTE-INPUT (heavy, Sirius-only):** the multi-GB `.grd` departure source grids, MARCS/ATLAS9
  atmospheres consumed by synthesis, the Turbospectrum/PySME engines, Gerber `.bin` grids, dep-grid
  HDF5. Read live to PRODUCE a number. → `config.constants.sirius_grid_path()` for the path (no
  existence check at import) + `require_sirius_grid()` / a loader loud-check at access +
  `assert_on_sirius()` gating the leg. **No local fallback, ever.**
- **COMMITTED-ARTIFACT (small, in-repo by ratified convention):** the KB-scale pre-derived NLTE
  delta-CSVs (`data/nlte_grids/*.csv`), the `amarsi2019_cno` tables, `data/threed_grids`. These are
  version-controlled RESULTS the Mac verdict path folds in — NOT live compute inputs — so they
  resolve **in-repo** via `committed_grid_artifact()` and do not loud-fail. Physically evicting them
  to the Sirius root is a **separate repo-wide migration** (RYA-559 successor); **no NEW
  compute-input grid may be added under `data/`** — route new source grids through the Sirius resolver.

### The resolver (single source, `config/constants.py`)
`SIRIUS_DATA_ROOT` (env-overridable) + `sirius_data_root()`, `sirius_root_present()`,
`sirius_grid_path()`, `require_sirius_grid(..., context=)` (loud-fail, no fallback),
`assert_on_sirius(context, require_subdirs=)` (refuse a compute leg off Sirius), and
`committed_grid_artifact()` (in-repo committed-result path). **No ad-hoc grid/atmosphere/engine path
literals** — every read is constructed through these.

### On-Sirius assertion (heavy-compute entrypoints)
`pipeline/nlte_bfactor_synth.py` (validate/derive), `pipeline/pysme_nlte.py` (`nlte_delta`), and the
`*_synth_sirius.py` / `ts_gerber_*` babsma/bsyn harnesses call `assert_on_sirius()` before touching a
grid or engine — so running them off Sirius fails with a clear "run this on Sirius" message naming the
missing `SIRIUS_DATA_ROOT`, rather than silently computing against whatever is local.

### Documented scope boundary (open, RYA-567 → follow-up)
The **Mac EW-measurement leg** (`pipeline/abundances_derive.py`) reads iSpec-bundled ATLAS9/MARCS
atmospheres + the MOOG/Turbospectrum engine via `ISPEC_DIR` (per-machine env, Mac default
`../ispec`). This is the **ratified Mac-banked measurement** (RYA-509: the solar anchor reproduces
from raw EXACTLY on the Mac; RYA-517 verified bit-identical Mac↔Sirius). It is NOT the grid-provenance
defect this standard targets, and forcing it Sirius-only would reverse Mac-banking and break the 9/9
gate — so it is **deliberately left on the local iSpec install** here. Whether to also move the EW
measurement leg's atmospheres to the Sirius root is an explicit architecture question for a follow-up
ticket, not silently decided inside this one.

## Abundances are reported PER INSTRUMENT and PER BAND — RYA-708

Ryan, 2026-08-09: *"we will keep everything separate for now, so a KITT abundance, a HARPS abundance etc. We can do a general abundance with everything as a showcase, with uncertainty, but I think the cool science will be per instrument and Band."*

### The law

An element's value is carried as **(instrument × band)**, not as a single number. A combined value may be published alongside, with its uncertainty, but it is a **derived showcase** and never the primary record. Nothing may collapse the per-instrument values into the combined one and discard them.

### Why — the combined number destroys the diagnostic

Solar Al I, 1D-LTE, same lines and the same measurement method throughout:

| line | state | HARPS (R~115k) | Kitt Peak (R~500k) | Δ |
|---|---|---|---|---|
| 6698.673 | clean | 6.279 | 6.260 | **−0.019** |
| 6696.023 | Fe I 6696.315 in-window | 6.629 | 6.537 | **−0.092** |

The two instruments agree to **0.019 dex** on the clean line and disagree by **0.092 dex** on the blended one — and the higher-resolution atlas returns the *lower* value, which is the direction resolving power predicts, because it charges less foreign flux to aluminium. Kitt Peak's equivalent width is 5.87 mÅ lower on the blended line and 0.82 mÅ lower on the clean one.

**So the cross-instrument delta is a blend diagnostic**, arrived at with no model, no synthesis and no line list — and it independently confirms the 12.775 mÅ of Fe I the synthesis blend census found. A single averaged A(Al) contains none of that.

The band split carries its own signal. Same element, same instrument:

| instrument | band | lines | mean | spread |
|---|---|---|---|---|
| HARPS | VIS | 2 | 6.454 | 0.350 |
| Kitt Peak | VIS | 2 | 6.399 | 0.277 |
| Kitt Peak | NIR | 4 | **6.415** | **0.075** |

The near-infrared band is roughly four times tighter than either optical set. Reporting one number hides which band earned the precision.

### Why HARPS-first was right, and what actually went wrong

Building against a single well-characterised instrument first is correct method: it keeps instrument systematics from being confounded with method bugs while the method is still being built. That was not the error.

The error is that **the scope decision was never recorded**. "We start with HARPS" silently stopped being a starting point and became an invisible boundary — and was then reported as a property of the sky, when the unresolved-element appendix printed "NO DATA, zero pixels" for lines two atlases on disk could see. **A scope choice nobody wrote down became a physical claim.**

### Standing rules

1. **Per-(instrument × band) is the primary record.** A combined value is derived, labelled as such, and carries its uncertainty.
2. **Cross-instrument agreement on a line is evidence about that line**, not noise to be averaged away. Disagreement concentrated on blended or problem lines is the system working.
3. **Never report a scope limit as a data limit.** If a wavelength was not measured because a run used one arm, say that; "no data" means no instrument we hold reaches it.
4. **Problem lines and problem children get the cross-instrument test first.** It costs one more measurement and needs no model, and it is the cheapest independent check the project has.

## All wavelengths, all data, all instruments, all models — RYA-708

Ryan, 2026-08-09: *"all wavelengths, all data, all instruments, all the models — that is the codex way."*

### The law

A measurement uses **every instrument that covers the line**, and a line is attempted on **every model that can reach it**. An element is not limited by the first instrument someone happened to build a pool from, and a line is not written off before every applicable rung has actually been executed.

Absence is a claim, and it carries the same burden of proof as a number. Three states are distinct and must never be collapsed:

- **not covered** — no instrument we hold reaches this wavelength. The only state that justifies "no data".
- **covered, not reachable here** — an instrument covers it but the file is on the other machine.
- **covered and loadable**.

`pipeline/coverage.py` is the single source for that question; `data/catalog/instrument_coverage.csv` is the registry. Asking whether *one file* reaches a wavelength, and printing the answer as though it were about the Codex, is the specific error this rule exists to prevent.

### Why it is a standard and not a preference

It was violated silently and at scale, and nothing in the repo could see it.

**The committed solar EW pool holds 808 lines and not one of them falls beyond 6910 Å.** The pool spans 3924.4–6905.3 Å; HARPS spans 3782.6–6910.0. Every measured equivalent width in the Codex came from a single instrument, while the **IAG** atlas (4047–10650 Å) and the **Kitt Peak** flux atlas (2960–13000 Å) sat on disk unused by that channel. Filtered to a usable depth window, **1466 IR lines across 23 elements** had never been touched. **Oxygen and phosphorus have zero usable optical lines** — their entire usable line set lies outside the instrument the pool was built from.

Aluminium is the worked case. It reported no value at all; its two pool lines over-claimed absorption by 1.8× and 4.7×; and its best lines — 7835/7836 and 8772/8773, graded **B/B+** against the optical pair's **C+** — were reported as "NO DATA, zero pixels" by an appendix that had asked one file. Measured across two atlases agreeing to 1–2%, they give **A(Al) = 6.415 ± 0.037** in 1D-LTE, against the **0.339 dex** spread of the optical pair.

### Standing rules

1. **Coverage is asked of the registry, never of a loaded array.** `W.min() <= x <= W.max()` is a statement about one file.
2. **A line reaches the appendix as unresolved only after every applicable rung has been executed** — EW/1D-LTE, then NLTE, then synthesis — with each rung's own outcome recorded. See Appendix A of the Science Product Package.
3. **Reaching a new wavelength is authoring, not guessing.** A line region needs `loggf`, `Ei`, `Ek` and `J` from a graded source; `nlte` set honestly; fit columns zeroed rather than inherited; and any inherited constant (damping above all) named as inherited.
4. **Two instruments covering one line is corroboration a single arm cannot give**, and it is worth the extra measurement.
5. **A model is skipped only for a recorded, citable reason** — no grid, no coverage, saturated core. "Nobody built the pool that way" is not such a reason.

## Ratified Constraints — structural re-check on every emission — RYA-674

Any element value, species selection, line inclusion, correction application, or gate result that has
been ratified by Ryan is a ratified constraint. Ratified constraints are protected by structural
re-check at emission time, not by discipline or by trust in cached state. Every emission path
(phase_c verdict, gold reference builder, disposition report, two-engine record) invokes
`pipeline/ratified_constraints.py::assert_ratified_constraints_satisfied()` before writing. Violation
is loud-fail, not warn-and-continue.

The pattern is deliberate: a ratified constraint that a downstream module could silently violate is
not actually ratified — it is a suggestion. The RYA-596 blank-cause honesty tripwire is the template.
Make the contradiction unrepresentable, not merely detectable. The alternative (each module
remembering to re-check each ratified constraint) has failed at least three times in the RYA-527 arc
(Fe method_scale, Li 1.409, Cr II 5.676) and will continue failing.

Adding a ratified constraint requires:

1. The ratifying Ryan decision (in the ratifying ticket's comments or description)
2. A new entry in `pipeline/ratified_constraints.py` registry, with the ratifying ticket's RYA-# as
   provenance
3. A test in `tests/test_ratified_constraints.py` that exercises the check

Removing a ratified constraint requires an explicit Ryan-ratified reversal in a new ticket. The
registry is append-only + revoked (with a revocation ticket), never silent-delete.

### The registry as it stands (RYA-674)

| `constraint_id` | type | ratified by | check semantics |
|---|---|---|---|
| `Li_6707_veto_1_409` | `FORBIDDEN_VALUE` | RYA-563 (RYA-103/458) | An element whose registry `required_treatment` is `upper_limit` (membership read from `problem_children.csv`, never a hardcoded element list) may not have its reported value sourced from the two-engine synthesis floor, and may never be emitted with verdict `PASS`. A floor **record** for such an element must carry its value under `diagnostic_only` / `diagnostic_value`, not `reported`. The phase_c upper limit itself (A(Li) 0.727) is the ratified treatment and passes. |
| `Cr_II_species_exclusion` | `EXCLUDED_SPECIES` | RYA-240 / RYA-558 | A species in `engine_selection.RATIFIED_EXCLUDED_SPECIES` (today: Cr II) may never carry a value in an emission unless the row is marked `diagnostic_only`. Cr is reported as Cr I. Scope is the **explicit, cited** list, not the registry-ion rule that additionally excludes Ti II / Si II from the selector — those are not ratified emission-time constraints. |
| `Fe_1D_3D_correction_required_on_solar_report` | `REQUIRED_CORRECTION` | RYA-553 (hardened RYA-681/674) | For every element with a reported-layer 1D→3D correction registered in `config/corrections_registry.yaml` (today: solar Fe), an element-level emission must (a) sit on one of the two recognised scales — a value on neither is the doubled-correction signature (7.416); (b) agree with its own declaration — a 3D value under a 1D declaration is gold v3's shape and re-arms the correction; and (c) sit on the **post**-correction scale, since a reported solar anchor carries the correction. Species-level diagnostic records are exempt: the floor's raw Fe I leg is not a claim about the reported anchor. |

| `unreliable_value_must_not_be_emitted` | `FORBIDDEN_VALUE` | RYA-679 / RYA-691 (ratified as a constraint RYA-699) | A row that DECLARES a reliability basis (`engineB_reliability` / `reliability_basis` / `reliability`) must declare a readable one: either RYA-679-gated, or `UNGATED — <why this artifact has no flag>`. A raw `reliable=False`, an `UNGATED` with no reason, and any string from outside the vocabulary are all refused — an unreadable basis cannot be told apart from an emitter that never looked, which is the RYA-691 defect restated. A row that declares **nothing** is out of scope: reliability is not computable over every artifact shape (the RYA-491/237 CNO cross-arm is a multi-indicator reconciliation, not a profile fit, so it has neither `dEW_dA` nor `railed`), and RYA-691 §3A forbids fabricating a uniform check over artifacts with genuinely different semantics. `diagnostic_only` is the ratified demotion here as elsewhere. The vocabulary is single-sourced in `pipeline/reliability_contract.py`, which both builds the basis and classifies it, so the producer cannot write a string the gate would refuse. |

### The scope line is the constraint (RYA-699)

`unreliable_value_must_not_be_emitted` is scoped to rows that *declare* a basis, not to all
rows, and that is the whole design rather than a weakening of it. Two failure modes sit on
either side. Scoped to every row, the check becomes a schema change — the verdict and gold
rows carry no reliability key, so the gate would fail everything until every emitter grew a
field it has no way to populate honestly, and the pressure would be to populate it
dishonestly. Scoped to nothing, RYA-691's six ungated reads come back the moment a seventh
consumer is written by someone who has not read that ticket.

What makes the narrow scope hold is that **silence and a false claim are different
failures, and only one of them is this constraint's**. A row that says nothing about
reliability is caught upstream, at the read, by the loud-fail RYA-691 installed in the
loader. A row that says something unreadable is caught here, at emission, in any module.
Neither check subsumes the other, and the reason the pair works is that both now speak one
vocabulary — before RYA-699 the basis strings were private f-strings inside a single
script, so a second emitter could produce a basis the gate had no way to read.

### Two row kinds, and why the distinction is load-bearing

`RowKind.ELEMENT_VALUE` is an element-level assertion ("element X's reported / frozen / proposed value
is V"). `RowKind.SPECIES_RECORD` is a per-species record in a diagnostic table, where a table may
legitimately carry species we would never report (the two-engine artifact records Fe II, Ti II and
Si II beside the reported ions).

Collapsing the two would be wrong in both directions. Applied to species records, the Cr II exclusion
would forbid the very diagnostic RYA-558 ratified keeping. Applied only to element rows, the floor
could keep writing Li 1.409 into `reported` for something downstream to adopt — which is exactly what
happened. So a species record may carry a vetoed or excluded species **only if it is marked
`diagnostic_only`**, which is the demotion RYA-558 and RYA-563 each specify in their own words.

### Correction bookkeeping is DATA, and the declaration is a list — RYA-674 / RYA-681

A frozen gold row declares which tabulated corrections its number already carries, in a
`corrections_applied` column holding a JSON list of identifiers from
`config/corrections_registry.yaml`. That list is the **single stored fact**. RYA-681's `scale_state`
and the human-readable `method_scale` prose are **views derived from it** at write time, never
independently computed — two stored facts that must agree is precisely how RYA-669 happened, and
`declared_scale()` raises if any two of them disagree.

`[]` is a positive statement ("no correction applied") and is written on every row. An **absent**
column means "this row predates the schema" and is a different state; the old empty-string fallback
silently meant "apply", which is the silence the doubled correction fell into.

`config/corrections_registry.yaml` declares **where a number lives, never the number**. Every quantity
is a `*_source` pointer into `config/constants.py`, resolved at load by
`pipeline/corrections_registry.py`. The "does this value look already-corrected?" bands are derived
from the magnitude — the two published scale centres are separated by exactly `|magnitude|`, so the
two-hypothesis decision boundary is their midpoint and the half-width is `|magnitude| / 2`. Revise the
tabulated offset and every band moves with it, with no edit anywhere else. A tabulated
`magnitude: -0.05` beside a tabulated `post_range: [7.45, 7.48]` would be three hand-maintained copies
of one literature value, none of which fails if only one is edited — the defect class this standard
exists to end, one level up.

`value_check_required()` states when the value-side check is mandatory rather than optional, and
computes it: a correction whose magnitude is at least half the acceptance gate's half-width cannot be
caught by the gate alone. Solar Fe is the archetype — `|−0.05|` against `FE_GATE`'s `0.05` half-width,
so a doubled correction lands exactly on the window edge and stays green (RYA-669 measured nine gate
tests passing on A(Fe I) 7.416). **`FE_GATE` is not narrowed to compensate**: it is a physical
acceptance window, not a bookkeeping check.

### Named inputs, not implicit ones

`scripts/phase_c_verdict_rya371.py --gold-version` and
`scripts/build_solar_reference_v2_rya522.py --gold-version` name the frozen gold an emission is
computed against, defaulting to `CURRENT`. The resolved version is stamped into the emitted summary,
so an artifact always says which frozen input produced it — and the guards run at full strength
against whatever is named. This is the sanctioned way to regenerate while `CURRENT` carries a row the
scale guard refuses to load. **A flag that skips a guard, or an in-memory repair of a frozen row, is
not.** A named input file is auditable; a silent correction is the pattern RYA-681 removed.

## A grade must name its subject — NIST grades gf, `MQ-` grades our measurement — RYA-711

Two grading systems in this project used the glyphs A–D, and they grade **different objects**:

| | `nist_grade` | our grade, now `MQ-*` |
|---|---|---|
| grades | the **atomic data** — % uncertainty on the transition probability | **our measurement** of the line |
| source | external lab/theory, 11 tiers AAA→E | internal composite of six weighted sub-scores → `line_score` 0–1 |
| depends on our spectrum? | no | yes — it moves if we re-observe |
| a `B` asserts | "log gf good to ≤10 %" | "composite `line_score` 0.60–0.80" |

Same glyph, unrelated claim. Ours therefore carries an explicit `MQ-` (measurement quality)
prefix — `MQ-A` / `MQ-B` / `MQ-C` / `MQ-D` — in the values themselves, not merely in the column
name, because a value gets copied into tables its column header does not travel with. **NIST's
letters stay NIST's and are never prefixed**: it is someone else's published scale and renaming it
here would misquote the source. Any table showing both must label each column with the object
graded, never a bare "grade".

`MQ-` rather than the `Q1–Q4` alternative because `Q` is spoken for: lines are **quarantined**, never
culled (below), and escaping one collision by starting another is not a fix.

Two pre-existing audit artifacts still carry the old bare letters —
`data/audit/procyon_outlier_loggf_rya281.csv` (which puts `nist_grade` and the old `line_grade` in
**adjacent columns**, the collision in the wild) and
`data/audit/fe1_scatter/fe1_per_line_residuals_rya407.csv`. Their generators now emit `mq_grade`;
the artifacts are historical snapshots and are corrected when their own tickets next regenerate them.

### Why the graded-gf firewall cuts at >25 % — the %→dex derivation

NIST grades are **percent** uncertainties on the transition probability; our gates are **dex**. They
meet only through `dex = log10(1 + ε)`. Against the ±0.10 dex RYA-561 ratification gate:

| NIST | ≤ % | dex | share of the ±0.10 dex gate | handling |
|---|---:|---:|---:|---|
| AAA | 0.3 | 0.0013 | 1 % | HIGH |
| AA | 1 | 0.0043 | 4 % | HIGH |
| A+ | 2 | 0.0086 | 9 % | HIGH |
| A | 3 | 0.0128 | 13 % | HIGH |
| B+ | 7 | 0.0294 | 29 % | HIGH |
| B | 10 | 0.0414 | 41 % | HIGH |
| C+ | 18 | 0.0719 | 72 % | falls through → MED/LOW |
| C | 25 | **0.0969** | **97 %** | falls through → MED/LOW |
| D+ | 40 | 0.1461 | 146 % | quarantined |
| D | 50 | 0.1761 | 176 % | quarantined |
| E | >50 | 0.2430 (at 75 %) | 243 % | quarantined |

The cut is **not** "25 % is a lot". It is the last rung whose gf uncertainty still fits inside the
gate: a NIST C is 0.0969 dex against a 0.10 dex tolerance — 97 % of it, **at** the gate but not past
it — while D+ is the first rung that exceeds the entire tolerance on atomic data alone. A line whose
gf could move the answer further than the gate allows cannot contribute to a value that gate is
meant to judge. Recorded because the code previously said only `# >25 % — cull`, which reads as an
arbitrary round number.

**C+ and C are in neither membership set, deliberately.** They fall through to reference-based
tiering, so a NIST C line can still land MED on a non-Kurucz reference. At 0.0719–0.0969 dex they
are at the gate, not past it, so they are reported **with the caveat stated** rather than excluded.
The consequence worth knowing: this is the one boundary where the tier is decided by
`loggf_reference` rather than by the NIST grade.

This cut is a correspondence between two published numbers — the NIST ASD accuracy ladder and the
RYA-561 gate. **It is not tunable**; it moves only if one of those moves.


## Frontier-band uncertainty is a RESULT, not a defect — RYA-777

Ratified by Ryan, 2026-08-11, and extended 2026-08-12. Append-only (RYA-674): the decision
is given, the citation is this section, and the clauses below are load-bearing rather than
advisory.

### VIS is the validated baseline

Solar elemental values in the visible are validated against the reference scale (Asplund
et al. 2021). Anchoring the metallicity to a known value **is** the validation — it is what
demonstrates that the atmosphere, the radiative transfer and the line physics are
calibrated. That is the trusted backbone, and it is the only band where a reference number
plays that role.

### Beyond that validation we do NOT chase better uncertainty

Larger error bars in the IR and near-UV are the honest measurement of a harder regime:
line crowding, molecular opacity, thinner NLTE grids, a pseudo-continuum that is never
directly observed. That uncertainty is captured, flagged and kept in the statistics. It is
a **result** — the map of where the frontier actually is — not a defect to be minimised.

A worked example of what this protects, from the Fe frontier chain: the near-UV median
line gap is 0.146 Å, *smaller than a strong line's own wings*, and the median continuum
sits at 0.607 — the true continuum is never observed there. A small error bar in that band
would be a false claim, not an achievement.

### Lines are excluded ONLY for physics or provenance

Ghost, blend, saturation, gf-provenance — the RYA-161 firewall and the RYA-711 quarantine.
**NEVER** to reduce a band's error bar, tighten its scatter, or move a band toward the
anchor. Excluding a messy-but-real IR/UV measurement to make a band look cleaner is the
prohibited move — "chasing the dragon".

The operational form of this is the exclude-vs-flag rule that RYA-807 wired and RYA-808
ratified: a line is removed from an aggregate only when its cause is **established**
(`required_treatment=exclude` **and** `status=active`). A line whose cause is not
established carries `investigate`/`owed`, stays in the aggregate, and is flagged. Removing
a line while the reason is still a hypothesis is the same error as tuning, wearing
different clothes.

### IR / near-UV validate against SAME-REGIME literature, not the optical anchor

Ratified 2026-08-12. When an abundance is derived in the IR or near-UV, its validation
reference is **published same-regime measurements** — not the Codex optical anchor and not
Asplund's optical-dominated recommended value. Two things, kept apart:

1. **Pipeline validation** — *does the method work?* Compare like-to-like against IR
   literature. That tests whether our IR pipeline reproduces what other IR studies get.
2. **The IR-vs-optical difference** — a **reported result**: the systematic offset between
   windows (IR NLTE, 3D, gf provenance, continuum). RYA-780 measured it for Fe:
   primary-sourced IR 7.5508, **+0.085** against the 7.466 optical anchor. That offset is
   the science, not a validation failure and not a tuning target.

Compare-don't-tune (RYA-161) binds both directions: we characterise agreement or
disagreement with the IR literature, and we tune toward it no more than we tune toward
Asplund.

**Code corollary, verified and enforced:** the optical `FE_GATE` [7.41, 7.51] must NOT be
applied to IR / near-UV band products. An IR value of 7.55–7.63 is the frontier result, not
a gate failure. Verified 2026-08-13: `FE_GATE` appears only in the phase_c / validation
scripts and nowhere in the band-product path, whose only gate is fit quality (RYA-342), so
a frontier cell cannot fail an optical gate by construction. This is the same gate-scoping
class as RYA-786 — a gate correct for VIS being wrongly applied to the frontier.

### Engines are presented, never ranked

Per RYA-712 and Ryan's 2026-08-11 ratification — *"there is no Primary. All Engines, LTE
and NLTE, are products that get presented."* A band cell carries every engine that reaches
it, side by side, each with its own value, σ, line count and reach. The cross-engine spread
is a diagnostic (RYA-525), never folded into an error bar, and a higher reach rate makes an
engine broader, not better.

### The solar model is a CHARACTERISED INSTRUMENT

VIS validated; IR and near-UV characterised with honest bars. Its value is that it can be
pointed at targets with no ground truth — Alpha Cen, 55 Cnc, eps Eri — and its per-band
uncertainties propagated honestly. **The discipline on solar, where we CAN check against a
reference, is exactly what lets us trust the number on a star where we cannot.** Shrinking
a frontier bar on the Sun does not make the instrument better; it makes every downstream
number on a ground-truth-free target quietly wrong.

### NO NEW GATE

This standard governs the **disposition philosophy**, not new machinery. The existing
apparatus already enforces "exclude for physics, not for cleanliness": the RYA-161 firewall,
the RYA-522 tiers, the RYA-711 quarantine, the RYA-586 per-element error-budget band, and
the RYA-807 registry gate. Nothing here adds a threshold, and nothing here may be used to
justify one.

---

## Use every tool we can — an archived model gets used where it adds value — RYA-817

Ryan, 2026-08-14: *"we use every tool we can. A capable model archived for convenience gets used where it adds value."*

### The law

Archiving a model is a **routing** decision, not a verdict on the model. A model that was
set aside for single-methodology consistency, pipeline simplicity, or any other
convenience is still a capable model, and where it can say something the live path cannot,
it gets run — as its **own data product**, beside the others, under RYA-712.

This is the natural consequence of *"all wavelengths, all data, all instruments, all
models"* (RYA-708) applied to our own shelf: a model we already hold and do not run is
indistinguishable, in the output, from one we never had.

### The reciprocal obligation, which is the whole of the discipline

Reactivating a model imports its **limits** along with its capability, and those limits are
usually undocumented in the code that wraps it. So a reactivation is not complete until:

1. **The model's DOMAIN is established from its own provenance, not assumed.** Not the
   parameter box its wrapper happens to check — the actual span of what it was built on.
   For a trained model that means recovering the training set; for a grid, the node
   coverage; for a fit, the data it was fitted to. RYA-817's worked case: the vendored
   Amarsi+2022 Fe MLP guards Teff/log g/vmic/A(Fe) and **nothing about the line**, so its
   171 + 12 training lines had to be recovered from the paper chain
   (`scripts/rya817_recover_amarsi_training_set.py`) before a single prediction could be
   trusted. They are optical, 4787.83–6810.26 Å, and nothing in the shipped code says so.

2. **Out-of-domain input is refused LOUDLY and named by axis** — never extrapolated
   silently, never quietly dropped. The size of the refused extrapolation is recorded as a
   diagnostic so the reader can see what was declined, but it may not enter a product.

3. **The reactivation is proved against a PUBLISHED result, on the author's own inputs.**
   Not against our own pipeline — that is the RYA-785 wrong-referee failure. RYA-817
   reproduces Amarsi+2022 Table 6's solar row (Fe I 7.47 → 7.46, Fe II 7.41 → 7.47) to
   ≤ 0.005 dex, using the Allende Prieto+2002 line list Amarsi actually used. Running the
   same control on the *training* list instead misses Fe I by 0.04 dex for reasons that
   have nothing to do with the engine — which is exactly why the control has to name its
   line list.

4. **Whatever the reactivation contradicts gets corrected at the source.** RYA-817's run
   showed the solar Fe I 3D-NLTE correction is ≈ 0.00 dex, not the −0.127 that
   `pipeline/nlte_corrections.py` had asserted in a docstring for months — a number
   reverse-engineered from a target rather than computed. The docstring was fixed in the
   same change.

### Standing rules

1. **"Archived" is never a reason not to run something.** It is a reason to check its
   domain first.
2. **A reactivated model is a NEW product key** (RYA-712), with its own value, σ and line
   count. It never merges into, replaces, or re-labels an existing engine's number.
3. **A domain check that admits everything is not a domain check.** Show it can reject —
   the RYA-805 rule, applied to model inputs rather than to spectra.
4. **A capability claim needs its domain stated in the same breath.** "We have a 3D-NLTE
   Fe engine" is true; "we have a 3D-NLTE Fe engine for the near-IR" is not, and the
   difference is 0 of 94 lines.


## The absence-of-evidence rule — an absence is a hypothesis, never a conclusion — RYA-833

Ratified by Ryan, 2026-08-15, after **six absence-claim misses in a single session** — every
one the same shape, every one falsified by measurement, and every correction *improving* the
outcome rather than costing us one.

### The law, verbatim

> **THE ABSENCE-OF-EVIDENCE RULE.** A claim that something DOES NOT EXIST or CANNOT BE
> REACHED — no lab data, no grid, no engine reach, no signal — is a hypothesis to be
> checked, never a conclusion to be asserted. Absence claims require a **dated, cited
> verification**, and they carry an **expiry**: the literature grows, our tables get
> ingested, our tools get re-pointed, so "we checked and it was not there" is only true as
> of when we checked. Before any absence claim GATES a decision (a product marked
> un-gradeable, a band marked unreachable, an element marked no-grid), re-verify it. Prefer
> specs that say "measure the reach" over specs that STATE the reach — the domain is
> empirical. Two sources of false absence to watch: (1) **architectural inference** — "the
> linelist ends, so the physics ends" (it does not); (2) **stale artifacts** — a
> column/tracker/survey describing a past state as if current. Counterweight, so this is
> discipline not paralysis: you cannot re-verify everything every turn — trust nothing
> LOAD-BEARING without a dated check, and treat an absence claim as load-bearing the moment
> it gates a decision.

### The evidence base — the six, dated and ticketed

| # | the claim, asserted | what measurement found | ticket |
|---|---|---|---|
| 1 | "Engine B can't follow past 9199.9 Å" | **187 of 239** followed | RYA-762 |
| 2 | "3D-NLTE is VIS-only" | true for Fe; **O I 7771 reaches the near-IR** | RYA-817 |
| 3 | "the Fe pool is 100 % Kurucz-floored" | a **lying column** — the values were fine | RYA-799 / 825 |
| 4 | "no lab gf in the near-UV" | **Belmonte 2017 covers the whole band** | RYA-822 |
| 5 | "only Fe/C/O/Li have off-solar 3D grids" | **Mg** (Matsuno 2024) and **Na** (Canocchi 2024) too | RYA-817 |
| 6 | "+0.033 is a rail artifact" | a clean **Elo dependence, corr 0.955** | RYA-831 |

**The tell, and it is consistent across all six.** Not one false absence came from a
measurement. Every one came from either **architectural inference** — the line list stops at
9199.9 Å, therefore the physics stops; the network splits on excitation potential, therefore
it is wavelength-agnostic — or from a **stale artifact** describing a past state as if it
were current: a `log_gf` column carrying intake metadata rather than the value the inversion
used, a survey with a cutoff date, a tracker generated from a superseded run.

The **positive** claims in the same session held up. It was specifically the negatives that
fell, which is why this rule is asymmetric and not merely "check your work".

### What it changes in practice

1. **An absence claim carries its check.** Not "there is no lab gf below 3780 Å" but "on
   2026-08-15, `fe1_lab_loggf.csv` (465 rows, 2132–11316 Å) carried 105 Fe I lines below
   3780 Å" — a statement with a date, a source and a number, which the next reader can
   re-run instead of inherit.
2. **Specs measure the reach; they do not state it.** A brief that says "the MLP is
   wavelength-agnostic" has decided the answer. A brief that says "domain-check each line
   and report the in-domain fraction" gets 0 of 94 (RYA-817) — the same conclusion, but
   earned, and with the *reason* attached.
3. **Scope an absence to what was actually examined.** RYA-822 reported "Ruffoni and Den
   Hartog have zero lines below 3780 Å". True of the two **partial VizieR extracts** in
   `data/linelists/primary_gf/`; false of the **sources**. The fix is to name the artifact
   inspected, never the source, unless the source itself was exhausted.
4. **A refusal is a finding, not a dead end.** Every one of the six corrections opened
   something: 187 recoverable lines, a near-IR 3D-NLTE line, two acquirable grids, a whole
   band of lab gf. Treating absence as provisional is not bookkeeping hygiene — it is where
   the results were.

### Relationship to RYA-674 — the procedural twin

RYA-674's ratified-constraint re-check catches **stale artifacts mechanically**: a boolean
`corrections_applied` schema and a shared constraint registry, so a value cannot silently
carry a correction it no longer describes. This rule catches **stale assumptions
procedurally**.

Same failure class — *something describing a state that reality has moved past* — attacked
from the two ends it can be attacked from. Where 674 can be made to fail loudly in code, it
should be; where it cannot, this rule governs. A guard is always preferable to a habit, so
if an absence claim can be turned into a dated, checkable artifact, turn it into one
(RYA-822's `nist_asd_FeI_3000_3780.prov.json` is the shape: the pull, its date, its access
route, and the reason the previous route was believed dead).

### Propagation

**Flagged for RYA-179** (Glossary / Method page / Science-Architecture): this rule belongs in
the narrative docs, and it accumulates there for the sync pass rather than being applied
piecemeal — per RYA-777/817, which added their standards here and left the narrative
propagation to 179.


## Line selection dominates gf — gf buys a tighter BAR, not a different value — RYA-842

Ratified by Ryan, 2026-08-16, after **five Fe/Al frontier tickets independently converged on
the same shape**. Unlike most entries here this one is *descriptive of current evidence*
rather than a rule, and it carries its own expiry — see the caveat at the end.

### The finding, verbatim

> **LINE SELECTION DOMINATES GF.** Across five Fe/Al frontier results, which lines enter the
> pool is first-order on the abundance VALUE and its line-to-line SCATTER; oscillator-strength
> work (grading, lab-gf re-pool) is first-order only on the reported error BAR. gf work buys a
> defensible, tighter uncertainty — it does NOT change the value or the spread.
> **Corollary for effort allocation:** when the goal is a different value or a smaller
> scatter, scrutinize line selection FIRST (depth cuts, blend rejection, EP range,
> saturation); reach for gf only to justify the error bar. Do not expect grading or lab-gf to
> move a value it cannot move (RYA-799/822 null results). This is descriptive of the current
> Fe/Al evidence — treat it as a prior to be re-checked per element, not an absolute
> (absence-of-evidence rule, RYA-833).

### The evidence base — five results, read from the committed artifacts

| # | ticket | what was changed | effect on VALUE | effect on SCATTER | effect on the BAR |
|---|---|---|---|---|---|
| 1 | RYA-817 | Fe VIS 7.586→7.466 decomposed | line selection **−0.097 / −0.141** | — | gf only **−0.026 / 0.000** |
| 2 | RYA-799 | graded the IR pool | **null** | **null** | grading alone, no σ gain |
| 3 | RYA-824 | swapped in primary lab gf | −0.026 IR, **0.000 VIS** | −0.009 IR, +0.007 VIS | **σ 0.200 → 0.052 / 0.060** |
| 4 | RYA-822 | graded the near-UV pool | 7.487 → **7.488** | 0.354 → **0.355** | — |
| 5 | RYA-836 | swapped in primary lab gf | 7.594 → 7.577 | 0.6513 → **0.6523** | tighter systematic |

**RYA-836 is the cleanest single demonstration**, because it varies both things in one run:

* changing the **gf** on the same 60 lines moved the scatter by **+0.001** (0.6513 → 0.6523)
  and the median per-line abundance by **exactly 0.000**, on a median Δlog gf of −0.026;
* changing the **line set** — the 61-line lab sub-pool against RYA-832's 40-line full pool —
  moved the scatter from **0.651 to 0.413**, a **0.238** difference.

Two hundred times the effect, from the axis nobody was working on.

RYA-824 shows the same asymmetry from the other side: a median Δlog gf of **+0.170** moved
the IR value by **−0.026** and the VIS value by **0.000**, while the gf systematic collapsed
from the blanket **0.200** to **0.052** (IR) and **0.060** (VIS). That is the entire return on
lab-gf work, and it is a real return — it is simply not the one people expect.

### Why this is worth writing down

Because the intuition runs the other way. "Bad oscillator strengths" is the reflexive
explanation for a discrepant abundance or a wide scatter, and it is a *plausible* mechanism —
which is exactly why it needs measuring rather than assuming. Four separate tickets spent
real effort on gf before the pattern was visible, and two of them (RYA-799, RYA-822) returned
clean nulls that only read as *results* once this shape was named.

The corollary is about **where to look first**, not about gf being unimportant:

* a **different value** or a **smaller spread** → interrogate the pool. Depth cuts, blend
  rejection, EP range, saturation, and which band the lines come from.
* a **defensible error bar** → gf. It is the right tool, and RYA-824's 0.200 → 0.052 is what
  winning looks like.

### ⚠️ This is a PRIOR, not a law — and it expires

Five results, on two elements, both in regimes where the pool is Kurucz-floored. That is
enough to reorder effort; it is **not** enough to conclude gf can never move a value. A
species whose lines are few, strong and gf-discordant could behave completely differently,
and the honest expectation is that some will.

So it is re-checked per element, not assumed — the same discipline as the absence-of-evidence
rule (RYA-833), of which this is the empirical companion. Both are what this session learned
about its own reasoning: RYA-833 about how we *state* things, this about where the error
actually *lives*. Neither licenses skipping the measurement.

### Propagation

**Flagged for RYA-179** — specifically the **uncertainty-breakdown section** of the
methodology doc, which this is the guide for: gf-floor vs lab-sub-pool vs line-selection vs
continuum, and which of those is worth spending on for a given goal. It accumulates there for
the sync pass rather than being applied piecemeal, per the RYA-777/817/833 precedent.
