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

---

## F NLTE grid availability — RYA-758

**Question.** Does any post-Lodders-2025 NLTE model atom or departure grid exist for
F I, or for the HF vibration–rotation lines used to measure solar fluorine?

**Why it matters.** F is the halogen extension of the CHNOPS story, it has a real
CRIRES+ IR diagnostic window at 2.3 µm, and Lodders, Bergemann & Palme 2025 punted on
it — their Table 2 carries A(F) = 4.4 ± 0.2 flagged "sunspot", quality index D, with
the substance deferred to Lodders & Fegley 2023 (Geochemistry 83, 125957).

**What the RYA-758 audit found (2026-08-09) — the question is answered, and the
premise was wrong.** The blocker is not the grid, it is the diagnostic. **No F line
is detectable in the quiet-Sun photospheric spectrum at all.** The accepted solar
value comes from HF features in *sunspot umbrae*: HF has a dissociation energy of
5.87 eV, so the molecule only survives in the cool umbral gas, and the measurement
line is the (1–0) R9 feature at 23358.3 Å (Maiorca et al. 2014, ApJ 788, 149,
DOI `10.1088/0004-637X/788/2/149`). Separately, the fluorine literature reports no
NLTE corrections for that feature and expects LTE to hold for a ground-electronic-
state vibration–rotation transition. So F is recorded `not_measurable` /
`blocked_spectroscopy` in `data/audit/nlte_grid_inventory_beyond28.csv`: **F is not a
deferred element, it is an out-of-scope one** for a pipeline that measures FGK dwarf
photospheres. Our targets have no umbrae to observe.

**What would answer it (i.e. reopen it).** Not an NLTE paper. Only one of:
(a) identification of a genuine F I or HF feature in a *quiet* solar/FGK photospheric
spectrum, or (b) a decision to extend the codex to cool giants/M dwarfs where the
2.3 µm HF band is measurable — at which point an HF NLTE treatment would become the
next question rather than the first one.

**Watch triggers.** Any Amarsi / Bergemann / Sitnova paper naming F; any CRIRES+ or
IGRINS solar-atlas work claiming a photospheric HF detection; an updated
Lodders & Fegley halogen review; a codex scope change to giants.

---

## Mo NLTE model atom watch — RYA-758

**Question.** Has anyone published an ab-initio Mo+H collisional data set and an
NLTE model atom for Mo I?

**Why it matters.** Molybdenum is the metal at the active site of nitrogenase (the
FeMo-cofactor) and of the molybdopterin enzyme family — biological N₂ fixation
happens on a Mo atom. In the RYA-758 audit Mo is the **only** outside-28 element
that clears two of the four science gates (bio-significant + neutron-capture tracer),
and its two diagnostic lines, Mo I 5506.49 Å and 5533.03 Å, are optical and sit
inside the HARPS arm that already carries our measured EW pool. Mo is blocked purely
by atomic physics: the day a model atom exists, Mo is immediately workable.

**State as of the RYA-758 audit (2026-08-09): still blocked, confirmed by the most
recent dedicated study.** Mishenina, Kurtukian-Nieto, Gorbaneva, Amarsi, Psaltis &
Pignatari 2026 (A&A 705, A38, arXiv:2511.21190, submitted 2025-11-26) derive Mo and
Ru for 154 disk giants **in LTE**, and close with: *"we emphasize that defining NLTE
corrections for Mo and Ru abundance measurements is essential, as these remain
currently unknown."* Amarsi is a co-author, so this is the grid-building community
itself recording the gap.

**One trap, recorded so it is not re-imported.** That paper quotes an MPIA NLTE
correction of ≈**+0.15 dex** in the same paragraph as Mo I 5506.49/5533.03. It is a
correction for the **analogous Cr I lines at 5208.41/5206.02** (matching
3d⁵(⁶S)4s a⁵S₂–3d⁵(⁶S)4p z⁵P°₃,₂ configurations), used as a qualitative proxy — **not
a Mo correction**. Secondary sources misreport it as one. MPIA Spectrum Tools was
enumerated live on 2026-08-09 and serves H, O, Mg, Si, Ca I/II, Ti I/II, Cr, Mn,
Fe I/II and Co only; there is no Mo grid there or anywhere else public.

**What would answer it.** A published Mo I model atom with quantum-mechanical Mo+H
inelastic collision rates, plus a departure grid or a per-line correction table whose
(Teff, log g, [Fe/H]) hull contains the Sun.

**Watch triggers.** Barklem / Belyaev / Yakovleva inelastic-collision output naming
Mo; any Amarsi-group Zenodo `Grid/NLTE` version (concept DOI
`10.5281/zenodo.3888393`) adding a Mo file — that family gained Cu in 2025-03, S in
2025-09 and Ag in 2026-05, so it is the live channel to watch; any paper titled
"Mo NLTE" or "molybdenum stellar abundance NLTE"; a follow-up from the Mishenina
group acting on their own closing sentence. Interim fallback if Mo were ever needed
before a model atom exists: the Cr I analogue systematic above must be propagated as
a named uncertainty term, per the draft clause (c) in
`docs/proposals/science_standards_grid_availability_amendment.md`.
