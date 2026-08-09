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

## Engines are separate data products and are NEVER combined — RYA-712

Ryan, 2026-08-09: *"the two engines should never be added together and presented, they are seperate data products… this model says this in VIS, this is what we got in IR, this is the Uncertainty for each product. It breaks down scientifically better. No Ambiguity."*

### The law

The reported product key is **(instrument × band × engine)**. LTE and each NLTE engine are **separate products**, each carrying its own uncertainty. A general combined abundance may be published alongside as a convenience, clearly labelled derived — it is never the primary record, and no per-product value may be discarded once it exists.

**Engines are never averaged.** Not across lines, not within an element, not "because they mostly agree".

### What this supersedes

This **narrows RYA-525's element aggregation**. The two-engine floor computes the reported value as an inverse-variance combine of the per-LINE winners, so one number could be part Engine A and part Engine B. That is exactly the mixing this rule forbids at the reporting layer.

### Where this rule came from — it is older than it looks

This was not invented on 2026-08-09. **RYA-489** (2026-06-30) already required per-arm σ reported by wavelength rather than collapsed, and carried the arm-boundary table that RYA-708 independently re-derived from `data/catalog/instrument_catalog.csv` two months later. **RYA-307** (2026-06-14) already made regime a first-class provenance layer. The conventions existed; the *EW pool* never honoured them — 808 lines, one instrument, nothing past 6910 Å.

What 2026-08-09 actually changed is narrower than a new rule, and it is worth stating precisely because all three prior issues share one defect: **they key the product on fewer axes than reality has.**

| issue | keys the product on | missing |
| --- | --- | --- |
| RYA-307 | star, with regime beneath | **engine** (postdates it — RYA-525) |
| RYA-306 | star × element × ion | **arm** and **engine** |
| RYA-489 | arm | **engine** |

Two substantive amendments follow, both recorded on the issues themselves:

1. **RYA-489 §3.2 is inverted.** It made the inverse-variance combined number the headline with the per-arm breakdown beneath. The per-(arm × engine) products *are* the product; a combined number is subordinate and derived. Combination stays legal **across arms within one engine** (RYA-237) and is illegal **across engines**.
2. **RYA-307 §2's premise is false as an operating assumption.** It reasoned from *"done right, UV and optical return the same number."* Measured: Mg 5711 disagrees across instruments by 24% while 5528 agrees to 2.9%; Al 6696 disagrees by 0.092 dex with the higher-resolution instrument reading lower; Gerber and Bergemann-MPIA differ by +0.044 dex on Fe and +0.11 on Ti. §2's own better instinct — *"the agreement IS the product"* — is the part that survives, once **disagreement is also treated as a product** rather than as residual to be absorbed by a recommended value.

The headline is therefore a **selection, not a blend**: the showcased (arm × engine) determination per element, named with its instrument, carrying the grade that justifies the pick. Every other determination is still shown. Nothing is suppressed.

The floor's machinery is not repealed and stays valuable: **reference-blind per-line selection remains a quality diagnostic**, and `CROSS_ENGINE_MIX_GATE` — which already flags an element whose winners span both engines with a large cross-engine Δ — was an early recognition of this same problem. What changes is that a mixed aggregate is no longer a reportable product. Where the floor previously emitted one blended value, it now emits one value **per engine**, and the cross-engine Δ becomes a published comparison rather than a threshold to be silenced.

### The infrared makes this concrete rather than theoretical

The floor was designed and ratified where two engines exist — the optical. In the infrared that assumption fails, and the failure is the normal case:

| IR capability | elements |
|---|---|
| both engines | Ti · Mn · Si · Ca · Mg · Na · Ba · O |
| Engine-B (Gerber TS-native) only | **Fe** · Ni · Co · Sr |
| Engine-A (b-factor grid) only | Al · C · Cu · K · Li · N |
| no NLTE model at all — LTE only | Cr · V · Zr · Sc · Y · S · Eu · P |

**Eight of twenty-six** support a two-engine comparison in the IR. **Eight have no NLTE at all** there, and their infrared lines are LTE however well measured. **Fe — the anchor — is Engine-B only**, on 4000 unmeasured lines, and its sole optical cross-check returned CHECK.

Under the old framing that reads as a degradation. Under this rule it is simply the honest product list: an element with one engine publishes one NLTE product, an element with none publishes LTE, and each says which it is. Ryan: *"That is ok! We show LTE and TS as the products."*

### Standing rules

1. **Never average engines.** A value is produced by one engine, or it is a labelled derived combination.
2. **Every product carries its own uncertainty.** A shared error bar across engines implies a shared systematic that has not been demonstrated.
3. **One engine is a complete answer**, provided the product says which engine and which band.
4. **LTE is a product, not a failure state.** For eight elements it is the only honest infrared answer available today.
5. **The cross-engine delta is published, not resolved away.** Where two engines exist, their disagreement is a measurement about the models and belongs in the record.

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

---

## Lines are QUARANTINED, never culled — RYA-711

Ryan, 2026-08-09: *"I dont like cull persay, every line still gets measured within its
ability, and logged if it is a C grade or lower basically. That way it doesnt saturate our
product measurement, but is still defined in a table in each elements appendix."* And:
*"Maybe Quarantined?"*

Same objection as *"nothing should be suppressed. i dont like that word"* — and the same
answer. The word has to describe what actually happens.

### The law

A line that fails a quality threshold is **QUARANTINED**: held out of the reported
aggregate, and otherwise fully retained — measured, derived, logged, and **printed in that
element's appendix table with its reason**. It is never deleted, never omitted from the
record, and never silently absent.

**Quarantine is reversible by evidence.** That is the whole reason the word fits. Aluminium
is the worked example: the RYA-398 graded-gf firewall held both solar Al lines, which was
recorded as the reason Engine B could not help Al — *"a departure correction with nothing to
correct."* Correct on its evidence. Then RYA-708 measured Al on Kitt Peak at 6.443 ± 0.068
over 6 lines and the premise died. A cull cannot be lifted; a quarantine can.

### What the code already gets right

`pipeline/curate_nonfe_pools.py::apply_cull` **deletes nothing.** It adds `cull_reason` and
`kept` and returns every row, and the per-element ledger `{element}_cull_rya395.csv` carries
KEEP/HOLD with reason. All **219** currently-held lines retain their measured EW. The
behaviour was always right; only the vocabulary overstated it — and vocabulary leaks into
how results get reported.

### Two gaps this exposes (RYA-711)

1. **17 of 219 held lines carry an EW but no derived abundance** (Ca 6/16, Cr 8/55, Mn 3/6).
   "Measured within its ability" means a held line still gets an A(X) **with its inflated
   uncertainty stated** — a NIST D at 50% is ±0.176 dex, which is a number, not a blank. A
   blank cell is indistinguishable from "never attempted", which is the defect condition in
   `docs/ELEMENT_PROTOCOL.md`.
2. **The ledger reaches no appendix.** Its only consumers are its own writer and one
   adjudication script. The per-element quarantine table is owed by Appendix A.

### Quarantine is the TOP tier, not the only one — three states, not two

Ryan: *"for the really super bad problem children lol."* Quarantine is reserved for the worst
offenders. The code already implements three tiers via `_gf_tier()`; only the middle one was
never named, which is why it looked like a binary keep/cull:

| tier | NIST | dex on gf | what happens |
|---|---|---|---|
| **clean** | AAA–B | ≤ 0.0414 | enters the reported product |
| **caveated** | C+ / C | 0.0719 / 0.0969 | **reported, with the caveat stated** — carries `GF_UNVERIFIED`, appears in the appendix. Falls through `_gf_tier` to MED/LOW by `loggf_reference` |
| **QUARANTINED** | D+ / D / E | > 0.0969 | the really super bad problem children — measured, derived, logged, tabled in the appendix, **never in the aggregate** |

The quarantine line sits where a line's *atomic data alone* could consume the entire ±0.10 dex
ratification gate (RYA-561). That is the derivation, and it is why the boundary is principled
rather than a round number. A C-grade line at 0.097 dex is *at* the gate, not past it — so it
is caveated, not quarantined. That distinction was already the code's behaviour and is now its
stated intent.

### The quarantine ledger and the problem-children registry are the same subject

`data/registry/problem_children.csv` (RYA-463) already holds **25** entries with `problem_class`,
`severity`, `status`, and `governing_tickets`. It is the appendix table Ryan is describing, and
**5 of its 25 are already `resolved`** — Mn 6013/16/21, Sr II 4077, Mg b-triplet, K I 7665/99,
Co I 3845 — which is the reversibility claim demonstrated rather than asserted.

**Open defect: two parallel vocabularies for one phenomenon, with no join.** `cull_reason` emits
`SAT / HIERR / BLEND / BADGF / GRADE`; the registry emits `SATURATION_COG / MOLECULAR_BLEND /
BAD_GF / …`. Same physics, different spellings, nothing links a quarantined line to its registry
entry. Single-source it on the registry's vocabulary — RYA-711.

### Naming

`cull` survives in identifiers and filenames until RYA-711 renames them together with their
consumers — a filename change breaks readers, so it lands as one deliberate change, not
opportunistically. **New code says `quarantine`.** No new artifact, column, or message may
introduce the word `cull`.

---

## A larger error bar is a result; a hidden one is a defect — RYA-713

Ryan, 2026-08-09: *"I am ok with IR having bigger Error Bars if the science backs it. we
will be entering new science territory with our measurements and finds of host stars."*

### The law

**Ungraded atomic data does not disqualify a measurement. It enlarges its uncertainty, and
that uncertainty must be stated, propagated, and attributed to its source.** Refusing to
publish a well-measured line because its oscillator strength is ungraded is not caution —
it discards a real observation and hides the fact that we made it.

What is forbidden is the *unstated* error bar, not the large one.

### The distinction that must not blur

There are two different things and only one of them is an error bar:

| | |
|---|---|
| **ungraded gf** | the central value is plausible, the uncertainty is large and quantifiable → **report with the wider bar** |
| **wrong gf (a ghost)** | the central value is wrong — the Sun disagrees with the line list by 5× in depth → **quarantine; no error bar rescues a wrong number** |

The Fe I IR pass already separates these: ghosts are quarantined with a named root cause;
the 271 survivors are lines whose gf is merely *unverified*. Those are publishable under a
stated budget.

### How the budget must be built

The gf contribution has two parts that behave completely differently, and collapsing them
is the way an honest-looking number goes wrong:

* **Random component** — line-to-line gf scatter. Averages down as `σ/√N`. With 271 lines
  and σ_gf ≈ 0.2 dex, this contributes ≈ 0.012 dex to the mean.
* **Systematic component** — a shared offset in the source's gf scale. **Does not average
  down at any N.** Kurucz semi-empirical values carry known systematics (RYA-161), so this
  is a floor, and reporting only `σ/√N` would understate the truth by an order of magnitude.

So an IR product reports **both**: the observed line-to-line scatter (which is measured, and
which also *tests* the random assumption), and the irreducible systematic floor from the gf
source (which is asserted from the source's own literature and cannot be beaten by more
lines).

If the observed scatter is far *smaller* than the assumed random gf error, that is a finding
in itself — either the gf values are better than their grade suggests, or something is
suppressing the scatter, and it must be investigated rather than pocketed.

### Why this matters for the mission

Host stars are new territory. A measurement with a wide, honest, attributed uncertainty
extends the record; a measurement withheld because its uncertainty is unfashionable
extends nothing. The IR and near-UV will carry wider bars than the optical for some time,
and saying so plainly is the science, not an apology for it.

**This does not weaken the ratification gates.** A promotion to a frozen value still needs
its evidence. What changes is that a *reported product* may carry a large stated
uncertainty, where previously the absence of a gf grade stopped it from existing at all.

### The optical is the control; the IR and near-UV are the frontier — and they need each other

Ryan, 2026-08-09: *"we want what we know to be close right? Fe in Visible should be pretty
dang close to Asplund."* And: *"but in IR and UV it is the wild wild west."*

Both statements are load-bearing, and together they define how a wide error bar earns trust:

**The optical can falsify the method.** A(Fe) is known there — Asplund 2021 gives 7.46, our
own banked 1D-NLTE anchor is 7.466 (RYA-553), and 808 HARPS lines stand behind it. So when a
new harness runs on optical Fe, **it must reproduce that**. There is no wide-error-bar
defence available: a method that cannot recover a known answer in the band where the answer
is known has not earned the right to report an unknown one anywhere else.

**The IR and near-UV cannot falsify it.** There is no reference value to miss. Whatever we
report is, for many of these lines, the first number anyone has published — which means the
uncertainty is doing all the work, and nothing external will catch us if it is wrong.

**Therefore the optical control is a precondition, not a courtesy.** Before any frontier band
is reported, the *same harness, same code path, same continuum policy* must be shown to
reproduce the known optical answer. That agreement is what licenses the frontier number; the
error bar alone does not.

Practically this makes every element's protocol two-legged:

| leg | band | test | what failure means |
|---|---|---|---|
| **control** | optical | reproduce the known A(X) | **the method is broken** — fix before proceeding |
| **frontier** | IR / near-UV | report with a stated, attributed budget | a wide bar is a result, not a failure |

The control leg also *calibrates* the frontier budget: the optical residual against a known
answer is a direct measurement of the harness's own systematic, and that term belongs in the
frontier error budget rather than being assumed to be zero.

---

## Band policy: each regime gets its own analytical approach, checked at intake — RYA-713

Ryan, 2026-08-09: *"we should have a check on intake maybe. Check instrument and check band.
UV gets treated different than Vis or IR due to crowding. IR has more errors, telluric lines
etc. That way the products are filtered differently. **Not tuned, but a different scientific
and analytical approach for each band.**"*

`pipeline/band_policy.py` declares, per regime, which method is physically valid and why.
`check_intake(wavelength, method)` resolves it and **fails loud** — a wrong-for-the-regime
method does not produce a worse number, it produces a **different quantity**, so silently
allowing it and widening the error bar afterwards would misrepresent what was measured.

### The measured basis (Kitt Peak + our line inventory, 2026-08-09)

| band | lines/Å | median gap | continuum p95 | continuum median |
|---|---|---|---|---|
| near-UV 3000–3800 | 4.62 | **0.146 Å** | 0.916 | **0.607** |
| VIS 3800–6910 | 1.87 | 0.277 Å | 0.963 | 0.811 |
| red-optical 6910–10000 | 0.34 | **1.872 Å** | **0.997** | 0.991 |
| NIR 10000–24000 | 0.14 | 3.989 Å | 0.956 | 0.862 |

### What each band gets, and the property that forces it

| band | synthesis | profile-fit | interval | telluric | forcing property |
|---|---|---|---|---|---|
| **near-UV** | ✓ | ✗ | ✗ | — | median gap **0.146 Å is smaller than a strong line's wings** — no interval contains one profile and excludes its neighbours, and there is no isolated profile to fit. Only simultaneous modelling is valid. |
| **VIS** | ✓ | ✓ | ✗ | — | gap 0.277 Å leaves lines resolvable. **The control band** — ground truth exists. |
| **red-optical** | ✓ | ✓ | ✗ | **required** | cleanest regime, and therefore the most dangerous: a broken method looks *healthier* here with no reference to catch it. O₂ A-band and H₂O are terrestrial. |
| **NIR** | ✓ | ✗ | ✗ | **required** | continuum p95 0.956 vs median 0.862 is telluric, not stellar — before correction the flux is not a stellar spectrum at all. |

**Interval integration is forbidden everywhere**, on evidence rather than principle: against
the HARPS pool over 146 optical lines it returned a median EW ratio of **0.773 (−0.112 dex)**
with a 5× spread, because a separation-derived window clips wings on crowded strong lines and
over-reaches on isolated weak ones.

### Why this is not tuning

Tuning is choosing a treatment for the answer it produces. Every field in `BandPolicy` is
keyed on **observable properties of the regime** — line density, separation, continuum
level, terrestrial contribution — all measurable without knowing any abundance, and all
measured before the policy was written.

Two structural guards, not just a convention:

1. `BandPolicy` has **no field capable of holding an abundance, target, reference value or
   tolerance**. It cannot express "use method X to get answer Y" because there is nowhere to
   put Y.
2. `assert_not_tuned()` runs **at import** and fails if such a field is ever added — verified
   by test: injecting `target_abundance` raises immediately.

The policy may be revised by re-measuring the regime and recording what changed. It may
never be revised because a different method gives a nicer number.

### Relationship to RYA-306

This is the **`arm` axis** RYA-306's method-selection matrix was missing. That matrix keys
method on `(star × element × ion)`; regime is independent of all three, because the failure
modes are local in *wavelength* as well as in stellar parameters. The two compose: RYA-306
says which method suits this element in this star, band policy says which methods the regime
can support at all, and the intersection is what may run.

### Handlers are keyed on METHOD; the band selects the method — RYA-713

Ryan, 2026-08-09: *"If I was architecting the code, I would have a handler class for each
case, UV, VIS, and IR. Different tools for different work. Am I wrong?"*

Right instinct, one refinement on where the seam goes. Keying handlers on **band** would
give near-UV and NIR each their own copy of synthesis, and VIS and red-optical each their
own copy of profile fitting — four handlers, two pairs of near-duplicates, which is exactly
the drift the Ba→Al copy produced. So:

* **`pipeline/band_policy.py`** routes: band → permitted method + band parameters.
* **`pipeline/measure/`** implements: one handler per **method**.
* Everything that differs *by band rather than by method* — telluric requirement, continuum
  treatment, systematic floor — is already in the policy and is passed to the handler.

| band | handler | telluric | continuum |
|---|---|---|---|
| near-UV | `SynthesisHandler` | — | pseudo-continuum only |
| VIS | `ProfileFitHandler` | — | true |
| red-optical | `ProfileFitHandler` | required | true |
| NIR | `SynthesisHandler` | required | post-correction |

Two handlers, four bands; a fifth band routes to an existing handler.

### The control is per-HANDLER and does not transfer

A **control run** takes a method to the band where the answer is already known — A(Fe) =
7.46 (Asplund 2021), 7.466 banked (RYA-553) — and requires it to reproduce that. It tests
the *method*, not the star, and the optical is the only band where any method can be
falsified.

**It is earned per method.** `ProfileFitHandler` passing at −0.0129 dex licenses nothing for
`SynthesisHandler`, which fails in entirely different ways: incomplete line lists, wrong
blend abundances, broadening, pseudo-continuum placement. `assert_controlled()` enforces
this — a handler may run un-controlled **only** in the control band (that *is* the control),
and is refused anywhere else until it has passed.

Two consequences that are easy to lose:

* **`systematic_dex()` raises for an uncontrolled handler** rather than returning zero.
  Assuming a zero systematic would understate every frontier error bar, and a handler that
  has never been checked has not measured its systematic — it has only avoided measuring it.
* **`CONTROL_TOLERANCE_DEX = 0.05`** is not a target to tune toward. RYA-561's ratification
  gate is ±0.10 dex; a harness that consumes half of that before any physics has not earned
  a frontier run.

### Method control vs cross-band science — two different comparisons — RYA-713

Ryan, 2026-08-09: *"It should not be forced to require the same exact value as VIS. We want
to audit what is known in IR… We want to compare to VIS, sure, and hopefully within error
bars it fits, but also check against what is found from other sources in the IR. And we
document regardless."*

These are different questions and only one of them may gate.

| | **method control** | **cross-band / external comparison** |
|---|---|---|
| compares | our EW vs a banked EW | our abundance vs another determination |
| same lines? | **yes** | no — different lines, sometimes different work |
| same band? | **yes** (optical) | no |
| asks | *does this tool measure what a validated tool measured?* | *does the Sun look the same from a different band or a different group?* |
| about | the **instrument + method** | the **star and the physics** |
| **gates?** | **yes** — `assert_controlled()` | **never** — `BandComparison` has no verdict field |

**Why the second must never gate.** Requiring the IR to reproduce the optical value is
circular: it tunes the frontier to the control and erases precisely what we are trying to
observe. Lines in different bands form at different atmospheric depths, carry different NLTE
departures, and sit on different gf scales. **A cross-band difference is a result.** If IR
Fe comes back 0.08 dex above optical Fe, that is either physics or a systematic we have not
yet named — and both are worth publishing. What is *not* acceptable is quietly forcing them
to agree and reporting the agreement as validation.

**Three comparands, all reported, none adjudicating:**

* `internal-cross-band` — our IR against our own optical
* `cross-instrument` — our Kitt Peak against our IAG on the *same* lines (all 103
  in-aggregate IR lines have IAG coverage; not yet run)
* `external-literature` — independent published IR determinations

Enforced structurally: `assert_not_a_science_gate()` raises if `BandComparison` ever gains a
`passed`, `verdict`, `tolerance` or `gate` field — verified by test. The moment one appears,
someone will gate a frontier band on reproducing the optical value and the pipeline will
begin filing real astrophysics as failures.

**And we document regardless of the outcome** — agreement and disagreement are both results,
and the appendix carries whichever occurs.

### The synthesis handler and the near-UV pseudo-continuum — RYA-713

`pipeline/measure/synthesis.py`. Measures an abundance by **fitting the observed flux** in a
window at a grid of trial A(X), taking the χ² minimum — the target element's abundance is
the only free parameter and every other species stays at the converged composition, so
blends are **modelled rather than avoided**. It wraps `abundances_derive._synth_flux_at_abund`
(the shared synth-v2 generator, RYA-287) rather than building a second synthesiser.

**The pseudo-continuum is the near-UV's real problem, and it is not small.** Measured on the
Kitt Peak atlas over a ±1 Å window:

| band | flux median | flux p95 | envelope median |
|---|---|---|---|
| **near-UV 3100** | **0.193** | 0.602 | **0.287** |
| near-UV 3400 | 0.868 | 0.991 | 0.991 |
| near-UV 3700 | 0.854 | 0.927 | 0.926 |
| VIS 5200 | 0.981 | 0.994 | 0.990 |
| red-optical 8000 | 0.971 | 0.990 | 0.984 |

At 3100 Å the flux never exceeds **0.70** anywhere in the window — the true continuum is
simply not present in the data. A synthetic spectrum, by contrast, arrives on the true
continuum. Comparing them directly compares two normalisations, and the fitted abundance
silently absorbs the difference.

**Treatment: apply the same envelope operation to both sides.** A percentile-based,
piecewise curve through the least-absorbed pixels, so it follows real curvature rather than
imposing a constant (at 3100 Å the envelope spans 0.194–0.680 where a flat 95th percentile
would have used 0.602 everywhere). The definition then divides out of the comparison.

**What does not divide out, and is therefore recorded on every near-UV measurement:** the
envelope depends on the synthetic spectrum, so a residual pseudo-continuum systematic
remains, and it **does not average down with more lines**. `PSEUDO_CONTINUUM_SYSTEMATIC_NOTE`
travels with each result rather than being folded into the scatter.

**Guards, all verified by test:** `prepare()` names every missing context key rather than
defaulting one (an invented atmosphere or line list silently changes what is being
measured); the NIR refuses to run unless the context declares `telluric_corrected`, because
before correction the observed flux is not a stellar spectrum; a χ² minimum landing on the
edge of the trial range is quarantined `FIT-RAILED` rather than reported; and every failure
**quarantines with a reason instead of raising**, since a raised exception loses the line
silently and a silent drop is the defect (RYA-429).

**Not yet controlled.** `assert_controlled()` refuses to let it run in the near-UV or NIR
until it reproduces the known optical answer in the VIS.

### Spectra live on both machines; the path is resolved, never hardcoded — RYA-713/567

Ryan, 2026-08-09: *"the spectra should be on Sirius as well."*

Compute runs on Sirius (RYA-567) and grids are Sirius-only, but **spectra are held on both**.
A path hardcoded to one host is a latent failure that surfaces as the wrong error message:
the first synthesis control reported *"no Kitt Peak segment covers 4065.381 Å"* for every
line — a coverage-shaped message for what was really a missing directory on the remote host.

`scripts/measure_band_ew.py` now resolves the atlas in order: `CODEX_KP_ATLAS` env var, the
Mac path, the Sirius path (`/mnt/codex-data/spectra/…`), and **raises naming every location
it tried** if none has segments. Kitt Peak is 44 MB / 252 segments — staged.

The window-extraction path (`--extract-windows` / `--windows`, ~0.05 MB for 20 lines)
remains available as a transport optimisation, but it is no longer the way a remote run gets
its data.
