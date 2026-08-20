# The Element Protocol — what "we have done this element" means

**RYA-711.** Ryan, 2026-08-09: *"we create 27 tickets, for each element, and we go through
each one, and verify, rerun, give it the Al treatment. If it fails, why, and document."*

This is that protocol. It exists because "the Al treatment" is ambiguous, and the
ambiguity is dangerous: aluminium was worked hard for a session and reached **4 of its 26
reachable lines**. Calling that done, and stamping it on twenty-six more elements, would
industrialise a 15% pass and call it a sweep.

---

## The definition of done

An element is **DONE** when every one of its reachable usable lines has an outcome. Not a
value — an *outcome*. A line that cannot be measured is done when the reason is recorded
and checkable.

```
usable          predicted central depth 0.05-0.60
reachable       inside an instrument we hold, per pipeline.coverage
DONE            every reachable usable line has an outcome
```

`scripts/line_accounting_rya709.py` is the scoreboard. An element is done when its
`unmeasured_reachable` is zero — every remaining line having an explicit, recorded
disposition rather than silence.

## The eight steps

Each is a gate. A step that cannot be completed is recorded as a stop, with its reason,
and the element carries that reason into the SPP appendix.

1. **ACCOUNT.** Run the line accounting. Know the four numbers before touching anything:
   usable, reachable, measured, unmeasured. Never start from the pool.
2. **COVER.** For every reachable line, which instruments see it. Never ask a loaded
   array — ask `pipeline.coverage`. A wavelength no instrument reaches is an acquisition
   target and leaves the element's scope.
3. **MEASURE.** Local continuum, window half-width from the **line separation** (not the
   FWHM — pair members must never share flux). Measure on **every** covering instrument,
   not the first one.
4. **GRADE.** A NIST ASD pull for every line that will carry a value. Unresolved
   components at one wavelength are **summed**, exactly as HFS is. Ungraded gf HOLDS the
   promotion — it does not block the measurement.
5. **REACH.** If the abundance path cannot see the line, author the line region from the
   graded pull: `loggf` / `Ei` / `Ek` / `J` sourced, `nlte` set honestly, fit columns
   **zeroed**, every inherited constant named as inherited.
6. **LADDER.** EW/1D-LTE, then NLTE, then synthesis — **executed, not narrated**. Each
   rung records its own outcome. Two states get named because they hide: *stopped early*
   and *escalated without cause*.
7. **CROSS-CHECK.** Where two instruments cover a line, measure both. Agreement is
   corroboration a single arm cannot give; disagreement concentrated on blended lines is
   the system working, not noise.
8. **REPORT.** Per (instrument × band), never one collapsed number. Every unresolved line
   defends its blank in SPP Appendix A with plots and measured evidence.

## The three honest outcomes

| outcome | meaning |
|---|---|
| **RESOLVED** | a value, per instrument and band, with its gf grades and its caveats |
| **BOUNDED** | an upper limit or a held value, with the bound's basis |
| **UNRESOLVED** | every applicable rung executed and failed, each failure recorded |

A fourth state is a **defect, not an outcome**: unresolved where a rung was never
attempted. That is a finding against the pipeline, not against the star.

## Aluminium, scored honestly against this

The worked example, and it is a **partial**:

| step | state |
|---|---|
| 1 ACCOUNT | done — 55 usable, 26 reachable, 26 unmeasured |
| 2 COVER | done — HARPS / Kitt Peak / IAG registered and verified |
| 3 MEASURE | **4 of 26** lines, on Kitt Peak; IAG depths only |
| 4 GRADE | done for those 4 (B/B+) and the optical pair (C+) |
| 5 REACH | done — four line regions authored |
| 6 LADDER | rung 1 done; **rung 2 blocked, no NLTE grid**; rung 3 not run |
| 7 CROSS-CHECK | done for 2 optical lines; **not done for the 4 IR** |
| 8 REPORT | Appendix A prototype published |

**A(Al) = 6.415 ± 0.037, 1D-LTE, Kitt Peak, NIR band.** Not promoted, nothing in the pool.
**Al is IN PROGRESS, not done** — 22 reachable lines remain untouched.

## What this protocol cost, measured

Aluminium's four lines took roughly one working session and produced **thirteen defects in
one adapted harness**, a false NO-DATA claim that reached the state register, a duplicated
instrument catalog, and six broken consumers from a CSV header. That is the honest rate.
Budget for it rather than being surprised by it.

## Step 3a — an element's CURATION travels with the element, not with the route

**Added 2026-08-19 (RYA-896/911/913 session).** Step 3 says measure on *every* covering
instrument. It did not say what a new instrument inherits, and that omission cost a full
session.

**The rule: when an element already has a curated pool, a new instrument or a new route
measures into THAT curation. It does not re-decide membership from scratch.**

Fe II is the worked example. RYA-352 established the solar Fe II EW-pool cull on a
computed, cited basis — saturation (REW past the COG-derived ceiling), measurement quality
(`ew_err/EW`), and blend flag — and applying it moves **A(Fe II; EW) 7.660 → 7.466** with
ionisation −0.144 → +0.050. That curation lives in
`abundances_derive._apply_fe2_ew_quality_cull` and is called from **`_load_solar_ews`
only**.

The band-product route never calls it. It re-derives membership from its own gates
(`REW > −4.9`, UNDER-PHYSICAL-WIDTH on the total Voigt FWHM, FEATURE-VERIFICATION), which
are reasonable gates and are **not the same gates**. On the six lines both pools contain
they disagreed on five:

| line | RYA-352 | band harness (as it stood) | after RYA-906/911 |
|---|---|---|---|
| 6084.102 | KEEP (clean) | ~~FIT-PINNED → rejected~~ | **KEPT** — the width was in the Lorentzian |
| 6456.380 | KEEP (clean) | ~~FIT-PINNED → rejected~~ | **saturation**, REW −4.889 |
| 6432.676 | CULL (HIERR) | kept — and its EW is +40% vs the original measurement | kept |
| 6247.557 | CULL (BLEND) | kept | kept |
| 6149.258 | CULL (BLEND) | kept | kept |
| 6369.459 | KEEP | kept | kept |

⚠️ **The third column is the correction, and it matters for how you read this section.**
FIT-PINNED — `abs(sigma_fit − floor) < 1e-4` — is **retired** (RYA-906/911, PR #315). It
rejected on the Gaussian sigma of fits that were measurably always Voigt, where sigma and
gamma are degenerate, so it recorded where the optimiser resolved the degeneracy rather
than how wide the line was. Its replacement tests the **total** Voigt FWHM against an
instrumental ⊕ thermal ⊕ microturbulent floor and fires on zero lines in both pools
checked. The disagreement in column two was therefore an artifact of the gate, and the
three disputed lines now fall out — or stay in — on physics.

The uncurated pool returned **7.656** where the curated one returns **7.466** — and a whole
session was spent decomposing that gap as a continuum defect, a loader defect and an
"unexplained residual" before anyone checked whether the curation was applied at all.

**How to apply, before MEASURE on a new arm:**

* **Name the curation that governs this element's pool**, and the code that enforces it.
  If the answer is "none", say so explicitly — that is a finding, not a default.
* **Confirm the route you are about to run actually calls it.** A route with its own
  quality gates is not the same thing; two independent gate sets on one element produce
  two pools and therefore two numbers.
* **A gate that removes a line the element's curation KEEPS is a defect to investigate,
  not a stricter standard.** FIT-PINNED rejecting 6084/6456 — the two lines RYA-352
  explicitly reinstated as clean — was the signal, and it was read as robustness.
  **This one was followed and it paid: the gate was wrong, not the lines** (RYA-906/911).
  The disagreement is worth trusting as a defect report about the gate.
* **Benchmark against the CURATED value.** Comparing a new pool to an uncurated historical
  number (Fe II EW-path ~7.700, which RYA-352 and RYA-715 both label blend/saturation
  inflated and diagnostic-only) makes a wrong answer look close and a right answer look
  broken.

**And the reason blends belong on the synthesis rung:** the EW route culls blended lines
because an equivalent width cannot separate two overlapping features. Synthesis fits them.
That is why the Fe II synth arbiter sits near 7.486–7.500 on lines the EW pool must drop —
the two routes are not competing measurements of the same pool, and LADDER (step 6) is
where that distinction is supposed to be made.
