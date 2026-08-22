# RYA-968 — Empirical gf-grade framework

**DESIGN SPEC — for sign-off. No implementation.** Child of RYA-958.
Firewall: RYA-161. Extends: RYA-855/850 gf rung. Doc sync owed: RYA-179.

---

## 0. Read this first — the measurement changes the ticket's premise

The ticket says the blanket 0.17 dex "should not drive our error." Before designing anything
I measured what an empirical replacement would actually say, on the RYA-959 VIS re-measure
(HARPS molecfit-corrected, 1D-LTE, in-aggregate, 3780–6910 Å), grouping each line's own
A(Fe I) by the gf tier RYA-945 assigned it:

| gf tier | n | mean A(Fe I) | **line-to-line sd** | median EW (mÅ) | median EP (eV) |
|---|---|---|---|---|---|
| **LAB** (primary laboratory) | **7** | 7.498 | **0.157** | **55.5** | 4.21 |
| NIST-C+ | 34 | 7.658 | 0.342 | 47.4 | 3.48 |
| OTHER | 86 | 7.657 | 0.311 | 47.2 | 4.26 |
| KURUCZ | 90 | 7.744 | **0.511** | **31.2** | 4.47 |
| VALD3 | 30 | 7.848 | **0.541** | **32.0** | 4.30 |

**Three findings, and each one changes the design.**

### A. The scatter really does order by gf tier — the hypothesis holds

LAB 0.157 → NIST-C+ 0.342 → KURUCZ 0.511 → VALD3 0.541, monotonic. That ordering was not
imposed; the tiers come from RYA-945's provenance ingest and the abundances from RYA-959's
re-measure, and nothing connects them but the physics. **This is the evidence that the method
is worth building.**

### B. 🔴 The empirical number is ~3× LARGER than 0.17, not smaller

Quadrature-subtracting the LAB floor:

| tier | sd | √(sd² − floor²) = empirical gf σ |
|---|---|---|
| NIST-C+ | 0.342 | **0.304** |
| OTHER | 0.311 | 0.269 |
| KURUCZ | 0.511 | **0.487** |
| VALD3 | 0.541 | **0.518** |

**The blanket 0.17 is not conservative — it is optimistic by about a factor of three.**
Replacing it with a measurement **widens** every ungraded error bar; it does not tighten them.

That is the opposite of the ticket's framing, and it is the single most important thing to
sign off on. RYA-850 already hit the smaller version of this — the cited lab σ (0.052–0.060)
came out *above* the 0.041 graded bound, and nothing was clamped. **The same discipline
applies here and it costs more.** If the framework is built and its answer is rejected because
the bars grew, we will have built a tuning knob.

### C. 🔴 The confound is present, and it is large

The high-scatter tiers are **systematically weaker lines**: LAB median EW **55.5 mÅ** against
Kurucz **31.2** and VALD3 **32.0**, with REW tracking it (−5.03 vs −5.27/−5.28). Weak lines
scatter more for reasons that have nothing to do with gf — continuum placement has less
leverage, per-line SNR is lower, and the EW→A inversion is steeper.

**So a naive quadrature subtraction attributes to gf a difference that is partly line
strength.** The ticket predicted exactly this ("done lazily it is the 0.17 with a local
accent"). It is not hypothetical here — it is measured, and it is in the direction that
inflates the answer.

### D. 🔴 The anchor is SEVEN LINES

In the flagship VIS cell the LAB pool is **n = 7**. The floor's own uncertainty is
`sd/√(2(n−1))` ≈ **±0.045 dex** on a 0.157 estimate. Quadrature subtraction is at its most
unstable when the subtrahend is poorly known, and matching on EW/EP to kill the confound (§C)
can only shrink 7 further.

**This is where the framework lives or dies, and it is a resourcing question, not a coding
one.** RYA-824 found the same wall — the lab tables hold 250 VIS Fe I lines and our measured
pool reached 9. RYA-945 lifted the *line list* to 199 LAB lines in VIS; it did not lift the
*measured* pool, which is still 7. §7 makes this a hard gate.

---

## 1. What is being replaced

`pipeline/error_budget.py` today:

```python
UNGRADED_GF_SYSTEMATIC_DEX = 0.17    # Kurucz semi-empirical (RYA-161)
GRADED_GF_SYSTEMATIC_DEX   = 0.041   # NIST grade B bound
```

and `pipeline/gf_rung.py` picks between them **per pool, all-or-nothing**:

| rung | condition | gf term |
|---|---|---|
| 1 | anything ungraded in the pool | 0.17 blanket |
| 2 | every line GF-LAB, <90 % carry a cited σ | 0.041 bound |
| 3 | every line GF-LAB, ≥90 % carry σ | RMS of the cited σ |

🔴 **The all-or-nothing rule is why RYA-855 moved 0 of 36 bars.** `decide()` returns rung 1 the
moment one line is ungraded — *"A pool is graded only if every line in it is"* — and every real
Fe cell is mixed (VIS: 7 LAB among 257). The rule is correct as written: a pool cannot inherit
a laboratory pedigree from a subset. **The fix is not to relax it. The fix is to stop grading
pools and start grading lines**, so a mixed pool is a quadrature over per-line σ instead of a
category that collapses to its worst member.

**That is the architectural core of this proposal.**

---

## 2. Method

### 2.1 Inputs

Per-line abundance tables that already exist, one per engine × band × instrument:
`data/results/band_products/{SPECIES}_{band}_{instrument}_{method}_{ENGINE}_lines.csv`
— carrying `wavelength_air_A`, `abundance`, `ew_mA`, `rew`, `ep_eV`, `red_chi2`,
`observed_depth`, `in_aggregate`, `excluded_reason`, plus RYA-959's `implied_width_A`.

Joined to `data/linelists/canonical_gf.csv` for `gf_tier`, `gf_sigma_dex`, `lab_source_tag`
(RYA-945). Join on wavelength **and** EP, never wavelength alone (RYA-780).

### 2.2 The floor — what the pipeline costs when gf is KNOWN

Within a **cell** (species × band × instrument × engine × treatment), take the GF-LAB lines
and compute the line-to-line dispersion of their per-line abundance. That is everything the
pipeline contributes when gf is not in question: continuum, blends, model atmosphere, NLTE
treatment, EW extraction.

Estimator: **MAD×1.4826, not sd** — n is small and one bad line must not set the floor.
Report both, and the disagreement between them, since a large gap is itself a diagnostic.

### 2.3 Per-line empirical σ_gf — three sources, in strict precedence

For each line, the gf uncertainty is the **first** of these that is available:

1. **CITED** — a published per-line laboratory σ (`gf_sigma_dex`, RYA-945/850). A measurement
   by someone with an apparatus beats any inference we draw from scatter. Median 0.030 dex
   (LAB) / 0.041 (NIST-C+).
2. **SELF-REPORTED** — the line's own spread across independent measurements
   (engine × band × instrument, where the same line is measured more than once). This needs
   **no anchor and no confound model**, because the line is compared only to itself. It is the
   most defensible number in the scheme and the ticket is right to give it precedence.
   Requires ≥ 3 independent measurements to be quoted.
3. **INFERRED** — the confound-controlled excess over the floor (§2.4). Used only where 1 and
   2 are unavailable.
4. **FALLBACK** — see §6. Not 0.17, and not Kurucz's.

### 2.4 The inferred term, with the confound controlled

Do **not** subtract raw tier scatters. The excess must be estimated at matched line strength.

**Primary method — matched comparison.** Bin lines by REW (the strength variable that showed
the confound) and by EP. Within each bin containing ≥ N_min graded and ≥ N_min ungraded lines,
compute `σ_gf(bin) = √(max(sd_ungraded² − sd_graded², 0))`. A bin without enough graded lines
yields **nothing** — it does not borrow a neighbour's floor.

**Secondary method — regression, as a cross-check that must agree.** Model per-line squared
residual against REW, EP, wavelength and a tier indicator; the tier coefficient is the gf term.
If the two methods disagree by more than their own uncertainties, **the framework reports that
disagreement and falls back**; it does not average them.

🔴 **The floor must be estimated from graded lines that look like the ungraded ones.** With the
LAB pool at median EW 55.5 mÅ and Kurucz at 31.2, there may be **no overlap bin at all** in
some cells. An empty overlap is a legitimate and reportable outcome: *this cell cannot be
graded empirically*. It is not licence to use the unmatched number.

---

## 3. Tiers — output is a classification, not a keep/drop

| tier | definition | disposition |
|---|---|---|
| **GOLDEN** | primary laboratory gf with a cited per-line σ | showcase; **never moves off lab gf** |
| **UNGRADED-CONSISTENT** | σ_gf within the cell's consistency bound | supporting product, tight empirical error |
| **UNGRADED-SCATTERED** | σ_gf measurable but wide | **quarantine** — kept, flagged, error sized to its *own* scatter, excluded from showcase |
| **INVALID** | artifact, not a line | removed from that spectrum's pool **with a recorded physical reason** |

**INVALID requires a named physical cause**, drawn from checks that already exist: RYA-959's
`gaussian_sigma_above_physical_ceiling` and `implied_width_exceeds_ceiling`, zero-flux telluric
cores, ghost/misidentification. **"Scattered" is never a reason for INVALID** — that is the
firewall boundary, and conflating the two is how a quarantine tier becomes a delete button.

The consistency bound is **declared before the run**, per §5.

### Quarantine, never delete (RYA-711)

Failed lines stay on the row. The failure *pattern* is diagnostic — a cluster failing in one
telluric region or one EP range names the next bug; every line stays accounted for, so a gap
is never silent; and today's scattered line may tighten under NLTE or after telluric
correction. RYA-959 already implements exactly this for the width ceiling; this reuses it.

---

## 4. 🔴 The firewall — how precision-not-accuracy is enforced STRUCTURALLY

RYA-161 in one line: **grade a line on its own precision, never on whether it agrees with the
expected answer.** A tightly-measured line that disagrees is a *finding*. A wildly-scattered
line that happens to sit on the reference is a *quarantine*.

Good intentions are not a mechanism. Four structural controls:

**F1 — The reference abundance is not an input.** The grading functions take per-line
abundances, EW, REW, EP, wavelength and gf provenance. They **do not** receive the reference
value, the gold value, or any pool aggregate. A function that cannot see the answer cannot
grade toward it. Enforced by signature, and by a test that greps the grading module for
imports of the gold/reference tables and fails on any.

**F2 — Grading is invariant under a constant offset.** Add a constant δ to *every* per-line
abundance in a cell and every tier assignment and every σ_gf must be **bit-identical**. Scatter
is offset-invariant; agreement is not. This single property is what makes "we graded on
precision" checkable rather than asserted. A dedicated test sweeps δ over ±0.5 dex.

**F3 — Thresholds are declared before the data is seen.** The consistency bound, N_min and the
bin edges live in a config block with a ticket reference and a fixed value; the run reads them.
Changing one is a diff, reviewable as such. A threshold chosen after seeing which lines it
excludes is the RYA-931 lesson — acceptance declared in advance, only the starting point
retried.

**F4 — The Cr canary is a blocking regression test.** Cr carries a standing +0.402 dex
residual (`pipeline/threed_corrections.py`). **If any change to this framework shrinks the Cr
canary, the framework is tuning and the build stops.** The canary must be *unmoved* by grading,
because grading is not supposed to be able to see it. This is the integrity check the ticket
names, wired as a test rather than a habit.

**A fifth, weaker control worth having:** report the mean-offset-by-tier (§0 table: LAB 7.498 →
VALD3 7.848) as a **finding, never as an input**. It is scientifically interesting — it may be
the Kurucz zero-point RYA-819/831 chased — and it must never enter a grade.

---

## 5. Declared constants (fill in at sign-off, before any run)

| name | proposed | rationale |
|---|---|---|
| `EMPIRICAL_FLOOR_MIN_LINES` | **12** | below this the floor's own σ exceeds ~20 % and the subtraction is unstable; VIS today has **7**, so VIS-HARPS **fails this gate as it stands** |
| `MIN_INDEPENDENT_MEASUREMENTS` | 3 | for a self-reported per-line spread |
| `CONSISTENCY_BOUND_DEX` | tbd at sign-off | CONSISTENT vs SCATTERED; declared, not fitted |
| `REW_BIN_EDGES` / `EP_BIN_EDGES` | tbd | must produce ≥ `EMPIRICAL_FLOOR_MIN_LINES` graded lines per used bin |
| `FALLBACK_GF_SIGMA_DEX` | measured (§6) | replaces 0.17; **expected ≈ 0.5**, not 0.17 |

---

## 6. Error-budget wiring

`gf_term(graded: bool)` is a two-branch switch and cannot express a mixed pool. Proposed:

```python
def empirical_gf_term(sigmas_dex, *, provenance) -> Term   # per-line σ, RMS'd
```

combining **per-line** σ in quadrature over the pool, exactly as `cited_gf_term` already does
for cited σ — RMS, never median, because these enter in quadrature and the median discards the
tail (RYA-850's stated reason).

Resolution order per line: **CITED → SELF-REPORTED → INFERRED → FALLBACK**, with the chosen
source recorded per line so a budget can state its own composition.

`UNGRADED_GF_SYSTEMATIC_DEX = 0.17` **is not deleted.** It is renamed to something that says
what it is (a Kurucz literature bound from RYA-161) and demoted to the last resort, reached
only when a line has one measurement, one engine, and no graded neighbour in its REW/EP bin.
🔴 **And even the fallback should be OUR measured value** — the pooled INFERRED σ across all
cells that passed the gates — with the Kurucz 0.17 kept only as the documented prior it always
was. On today's evidence that fallback is **≈ 0.5 dex**, and §0B is the argument for accepting
that rather than flinching from it.

**Every product's budget must state which source each line used.** A budget that cannot say
where its gf term came from is the RYA-873 defect ("MEASURED" printed beside a zero nobody
measured).

---

## 7. Gates — what must be true before implementation starts

1. **The anchor must reach `EMPIRICAL_FLOOR_MIN_LINES` in at least one cell.** Today VIS-HARPS
   has **7 GF-LAB lines in-aggregate**. Either the measured pool grows to reach more of
   RYA-945's 199 VIS LAB lines, or the first pilot runs somewhere the anchor is larger. **This
   gate is currently RED and it is the reason this is a design ticket and not a build ticket.**
2. **A matched REW/EP bin must exist** containing enough graded *and* ungraded lines. If LAB
   lines are all strong and Kurucz lines all weak with no overlap, §2.4 yields nothing.
3. **F2 (offset invariance) must pass** on synthetic data before any real grading runs.
4. **The Cr canary must be unmoved.**

---

## 8. Element-agnostic by construction

Fe I is the pilot **only because it has a laboratory anchor**. Nothing in §2–§4 is Fe-specific:
the floor is "the dispersion of lines whose gf is known", which is defined for any species with
a graded subset. `gf_rung.LAB_GRADED_SPECIES` is `{("Fe","I")}` today and the framework must
read it rather than assume Fe.

**The extension to un-anchored elements (Al, and the rest of the backbone) is the real prize
and it is NOT free.** Without a lab anchor there is no floor to subtract, so those elements can
use routes 2 (self-reported cross-engine spread) and 4 (fallback) but **not** route 3. Whether
a floor measured on Fe I transfers to Al is an open question, and the honest default is **it
does not** — a floor is a property of a pipeline *and a line list*, and Al's differs. Any
cross-element transfer needs its own evidence, not an assumption. **Recommend: state the
Fe-only scope in v1 and open transfer as its own ticket.**

---

## 9. Open questions for sign-off

1. **§0B — do you accept a ~3× wider ungraded gf term?** This is the decision. If the answer is
   "only if it comes out smaller", the framework must not be built.
2. **§7.1 — how do we grow the 7-line anchor?** Extend the measured pool toward RYA-945's 199
   VIS LAB lines, or pilot in a different cell? This gates everything.
3. **§5 — the consistency bound.** Declared where, and by whom, before the first run?
4. **§8 — Fe-only in v1?** Recommended.
5. **RYA-179 doc sync** — this changes the published error-budget methodology on the Method
   page. Owed in the implementation ticket, not this one.

---

## 10. Evidence index

All numbers in §0 are reproducible from committed data:
`data/results/band_products/FeI_3780_6910_harps_solar_harps_molecfit_corrected_PROFILEFIT_1D-LTE_lines.csv`
(RYA-959, 257 in-aggregate lines) joined on wavelength ±0.02 Å to `data/linelists/canonical_gf.csv`
`gf_tier` (RYA-945). Floor uncertainty from `sd/√(2(n−1))`.
