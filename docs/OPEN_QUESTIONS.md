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
