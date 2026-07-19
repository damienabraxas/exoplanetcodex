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
