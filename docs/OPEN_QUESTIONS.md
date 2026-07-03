# OPEN QUESTIONS

## RYA-498 — What does reflex do differently on the edge order vs standalone `wave_THAR`?

**Answer (2026-07-02): nothing — reflex fails `wave_THAR` identically.** This was
the brief's Step 6 question, premised on reflex iterating the wave bootstrap to
converge the edge order that broke the standalone hand-driven `wave_THAR`. That
premise is **disproven by the run**:

- esoreflex 2.11.2 orchestrated the full cascade end-to-end (mdark → orderdef →
  mflat → cal_contam → wave_FP → **wave_THAR** → sci_red; dataset saved to
  `reflex_end_products/`).
- reflex invoked `espdr_wave_THAR` **~20 times** (`wave_THAR_1` ×10 + `wave_THAR_2`
  ×10) — every attempt failed with the **same** error as the standalone run:
  `espdr_find_first_FP_ll: No FP line before/after the ThAr line (order 1, …)` →
  `espdr_get_all_FP_ll_per_order failed: Input data do not match` (status 13).
- So `wave_THAR` is a **deterministic recipe failure on this data**, not a
  convergence-needs-iteration problem. Reflex does not re-seed the recipe with a
  refined wave prior between attempts; each run uses the same inputs (the
  CALIB_DATA_DIR static wave matrix, latest shipped = **2022-11-01**) and fails
  identically.

**Root cause (both methods agree):** the edge order's FP-line → ThAr-line
association can't be built against the **~6-month-stale 2022-11-01 static wave
first-guess** for our 2023-04-29 data. The ESO archive α Cen NIRPS ADPs (2023–25)
reduced fine because the operational pipeline carries a **contemporaneous**
wave matrix (updated regularly), not the shipped static. NIRPS ships no wave
matrix newer than 2022-11-01.

**Second, independent blocker at `sci_red`:** the solar frames are
`DPR.CATG=CALIB`, `OBJECT=SUN,FP,G2V`, and **lack** `ESO OCS TARG SPTYPE` and
`ESO TEL TARG RADVEL` → `espdr_sci_red: espdr_get_science_params failed`. These
are science-target keywords a normal `OBJECT` science frame carries; the HELIOS
CALIB monitoring frames don't. Even a converged wave wouldn't yield an S1D
without supplying these (Sun: SPTYPE≈G2V, RADVEL≈0 — but that is header injection
into raw frames, a deviation to decide, not do silently).

**Net:** no `S1D_FINAL` produced. The reflex *engine* is fully working
(orchestration, DO, cascade, bookkeeping) — the remaining blockers are (1) a
fresher wave prior and (2) science-target keywords, both needing a decision:
- (1) obtain/derive a 2023-epoch `WAVE_MATRIX` prior (format + 3.2.6→3.3.12
  provenance to manage), or open an ESO support ticket (reflex is ESO-supported)
  on wave_THAR edge-order convergence with the shipped static;
- (2) decide whether to inject `SPTYPE`/`RADVEL` for the Sun, or treat these
  CALIB frames via a non-science extraction.

No silent order drop occurred — the failure is a hard stop at `wave_THAR`,
reported loudly here.

### Isolation test with the ESO demo (nirps-demo-reflex-0.4.0) — partial, inconclusive on wave
Ran the ESO known-good NIRPS reflex demo (18 GB; 36 HA raw frames incl. a real
`SCIENCE OBJECT,FP` frame + 5 darks + full cals, dated **2022-12-06** = ~1 month
past the 2022-11-01 static — vs our ~6 months). Purpose: confirm the engine is
correct and pin our `wave_THAR` failure to prior-staleness.
- **Confirmed:** the demo's real `SCIENCE OBJECT,FP` frame **auto-targets and
  builds a dataset** (DO reaches the `SCIENCE_FP` action, virtual/non-virtual
  calib tree assembled) — validating that our CALIB `SUN,FP,G2V` frames genuinely
  needed the `REFLEX.TARGET` toggle (a real SCIENCE-vs-CALIB difference), and that
  the engine builds datasets for proper science frames.
- **Inconclusive on wave:** the demo run then failed at an **mdark recipe-config
  error** (`ovsc_sig_clip_method parameter not found in input sop`) — *before*
  the wave step, and a *different* error than our solar run (whose mdark ran
  fine). Likely an invocation/recipe-config quirk of the demo run (launched with
  `-ROOT_DATA_DIR`/`-RAW_DATA_DIR` overrides on a fresh bookkeeping tree), not a
  wave result. Chasing it is out of the review's time-box; the wave finding
  stands on the two independent solar runs (hand-driven + reflex), which fail
  `wave_THAR` identically.

### 2026-07-03 — refined diagnosis (order-exclusion route) + make-do exit
Ryan's reframe: exclude the broken edge order (don't chase convergence); hard
time-box; a reference tolerates a dropped order. Investigating that **overturned
the "order 1 edge" framing and found the real cause:**

- **Order counts match:** orderdef traced **71** orders; the static
  `WAVE_MATRIX_A` covers **71**. So order 1 is *covered*, not a spurious extra
  edge order → not the "constrain the range" case.
- **The FP-line starvation is GLOBAL, not edge-specific.** The wave_FP
  `FP_SEARCHED_LINE_TABLE` has a **median of 3 FP lines/order** across all 71
  orders (order 1=4, order 2=3, orders 68–71=2–3; only orders 3–4 have 14–27).
  Total ≈ 363 lines. **reflex's own run got the identical 363/median-3** — so it
  is inherent, not a hand-config artifact. Excluding order 1 cannot help:
  `wave_THAR` would just fail at the next starved order (order 2 has 3).
- **The FP frame is good; the *matching* is starved.** The extracted `S2D_FP_FP`
  spectrum shows a proper dense FP comb (~**745 peaks in a mid-order**). So the
  FP data is fine — but the FP-line search keeps only ~3/order because it matches
  against the **~6-month-stale 2022-11-01 static wave first-guess**, which
  mis-locates the comb for our 2023-04-29 data. This is the "covered-but-
  mis-placed → stale first-guess" case, now confirmed precisely.
- **`wave_THAR` / `orderdef` expose no order-range / FP-threshold / exclude
  parameter** (full `--params` checked) — so there's no clean recipe-level lever
  to constrain the order set anyway.

**Conclusion:** the from-raw solar wave solution can't be built with the shipped
static prior at this epoch gap; the only real fix is a fresher (~2023) HA
`WAVE_MATRIX` first-guess, which is **not trivially grabbable** (the demo run
failed at an mdark-config error; the archive α Cen ADPs are 3.2.6 — a
format-provenance project Ryan scoped out; no ESO ticket). Per the make-do
time-box → **STOP**.

**Make-do exit — fallback reference (parked state):** use the on-disk external
solar IR atlas as the YJH reference instead of a self-reduced NIRPS solar S1D:
`data/spectra/exoplanetcodex-data/Solar Calibration/IR Reference Atlases/ACE-FTS`
and `/NSO_photatl` (Kitt Peak), with the explicit caveat that these are **not
instrument-matched to NIRPS** (resolution/LSF/sampling differ). Telluric verify
still available via `/Wallace_telluric` (RYA-390). IR *science* is unaffected —
it runs on the reduced α Cen ADP in the repo (RYA-507).

**Reusable reduction diagnosis (feeds RYA-508):** a NIRPS from-raw reduction with
the shipped static wave prior fails `wave_THAR` when the science epoch is far
(~6 months) from the newest static (2022-11-01) — not via an edge order, but via
*global* FP-comb-match starvation (median ~3 lines/order) against the stale
first-guess. Any future from-raw NIRPS target near an up-to-date static will not
hit this; targets far from 2022-11-01 need a contemporaneous `WAVE_MATRIX`.

---

## iSpec: theoretical-EW / synthesis child silently returns zero-init on child death — RYA-506

**Where:** `ispec/synth/spectrum.py`, `calculate_theoretical_ew_and_depth` (and the sibling
`generate_spectrum`). The SPECTRUM work runs in a child `multiprocessing.Process`; the parent
initialises `output_ew = np.zeros(num_lines)` and polls a `JoinableQueue` with
`while p.is_alive() and num_seconds < timeout:`. **If the child dies without enqueueing a
result, the loop exits with `output_ew` still all zeros and iSpec returns them with NO
exception and NO log** — only the separate *timeout* branch logs; the dead-child branch is
silent.

**Impact for us (RYA-506):** on macOS the Python 3.8+ default start method is `spawn`, under
which the iSpec child (which imports the C `synthesizer` extension and expects fork semantics)
dies on spawn → silent all-zero theoretical EWs → 100% of the Procyon Fe pool failed the
`theo < 5 mÅ` quarantine → "MOOG: No abundances." A defensible calibration number
(A(Fe I; NLTE) ≈ 7.571) was made to look unreproducible by an unreported subprocess death.

**Our mitigation (already applied):** force `multiprocessing.set_start_method('fork')` at
pipeline import (`pipeline/abundances_derive.py`), plus a loud all-zero-batch guard in
`_fe2_theoretical_ew` that raises instead of caching/quarantining on a failed synthesis.

**Upstream ask (report to the iSpec project):** the dead-child branch of the poll loop should
**raise** (or at minimum log LOUDLY) rather than return the zero-initialised array — a silent
all-zero synthesis result is indistinguishable from "every line is a null line" and corrupts
any downstream quarantine/abundance step. Ideally iSpec should also let the caller select the
start method / use a `fork` context explicitly for the synthesizer child on macOS.
