# RYA-959 — Fe I VIS re-measured, and the width ceiling that was never there

**Status:** measurement + Kitt Peak control landed; HARPS Engine-B re-deriving.
**Artifacts:** `data/measured/band_ew/FeI_3780_6910_*`.
**Code:** `pipeline/line_width.py` (the ceiling), `scripts/measure_band_profilefit.py`,
`pipeline/measure/profile_fit.py`, `pipeline/band_products.py` (two new columns).
**Test:** `tests/test_line_width_ceiling_rya959.py`.

## Why

RYA-958 found that 251 of 444 Fe rows in `data/measured/sol_ew_results_v1.csv` imply
physically impossible line widths — median implied FWHM 0.70 Å, 30 % above 2 Å, and
Fe I 5602.541 carrying 269 mÅ of absorption on a core 0.008 deep, which needs a 31.8 Å
Gaussian to hold it. Every one of those rows passed the harness's χ² and error checks.
This ticket re-measures the band from scratch and asks whether the fault is stale data or
live code.

**It is live code.** The fresh run reproduces the defect, and the mechanism is now named.

## The RCA — a floor with no mirror

`pipeline/line_width.py` was created by RYA-906/911 to hold one guard: a line may not be
NARROWER than instrumental (+) thermal (+) microturbulent broadening. Nothing ever
convicted a fit for being too WIDE. Both profile fitters bound the optimiser at

```
sigma_max = 0.40 Å          # instrument_sigma() / ProfileFitHandler.widths()
```

described only as *"generous; rotation and unresolved blends legitimately broaden lines"*.
Measured against the star's own ratified parameters, the widest Gaussian σ solar Doppler
physics permits at 5500 Å is **0.083 Å** — instrumental ⊕ thermal ⊕ microturbulent ⊕
macroturbulent ⊕ rotational, with `vmac` and `vsini` both counted as Gaussian σ, which
overstates both. The bound is **~5× outside physics**, and a fit with nothing good to
converge on ran to it and integrated a pedestal.

The fitted-σ distribution over 2 342 HARPS Fe I fits is **bimodal, with the ceiling
sitting in the valley**:

| fitted σ (Å) | fits | |
|---|---|---|
| < 0.030 | 832 | physical cores |
| 0.030 – 0.100 | 195 | the valley — the ceiling range is 0.057–0.104 Å |
| 0.100 – 0.399 | 354 | runaway |
| **= 0.400** | **961** | **pinned exactly on the optimiser's bound** |

Only **3.2 %** of fits fall between 1× and 2× the ceiling. The threshold is not chosen —
it is derived from `STAR_PARAMS`, and it happens to land in the sparsest part of a
distribution that separates on its own. A pinned fit returns `red_chi2 ≈ 0.02`: it looks
converged, which is why nothing objected for so long.

## What was added

Two checks, in `pipeline/line_width.py` so **both** fitters import them. RYA-906/911's
floor reached only one of the two sites and the other kept a refuted test for a whole
release; the module exists so that cannot recur, and a test asserts both sites resolve to
the same functions.

1. **`gaussian_sigma_above_physical_ceiling`** — the Doppler bound. `gamma` is
   deliberately **not** bounded: pressure broadening is Lorentzian and real, and a ceiling
   on it would quarantine the strongest lines in the band, which is exactly where
   RYA-945's laboratory gf backbone lives.
2. **`implied_width_exceeds_ceiling`** — the integral bound, and the only check in this
   harness that referees the EW itself rather than the profile's parameters. It asks what
   FWHM a Gaussian of the line's **observed** core depth would need to hold the EW.

Both quarantine, never cull (RYA-711): the EW stays on the row so a reader can see how far
off it was, and `attribute_root_cause` still records why the fit had nothing to hold.

`observed_depth` and `implied_width_A` are now **columns** on every row, passing or
failing. RYA-958 had to reconstruct this quantity by hand because the measurement never
recorded it — the RYA-911 mistake, repeated.

## 🔴 Two thresholds that were wrong, and the measurements that fixed them

**The implied-width ceiling was tuned, and is not any more.** It began as
`IMPLIED_WIDTH_ALLOWANCE = 3.0` — a multiple of the Doppler FWHM chosen to look generous.
Measured on the σ-clean fits, the implied/Doppler ratio is a **smooth continuum**:

| percentile | 50 | 75 | 90 | 95 | 99 | max |
|---|---|---|---|---|---|---|
| implied ÷ Doppler FWHM | 1.22 | 2.06 | 2.90 | 3.34 | 4.04 | 4.81 |

No valley anywhere, so any multiple cuts a populated distribution — the tuning this
project forbids. It was replaced by a fact about the harness:
`pipeline.lines_fit._integrate_profile` integrates the model over the **fit window only**,
so an equivalent rectangle wider than that window claims more absorption than the interval
can hold at the observed depth. That is arithmetic, not taste, and it needs no allowance.

**The depth must be the OBSERVED one, and RYA-958's diagnosis used the model's.** Running
the same test with `central_depth` from `linelist_solar.csv` — a VALD *model* quantity —
rejects **434 of 999** physically sound fits:

| depth used | median implied width | max | > 0.5 Å |
|---|---|---|---|
| **observed** (this guard) | **0.234 Å** | 0.692 Å | 58 / 999 |
| predicted (RYA-958's diagnostic) | 0.398 Å | 10.8 Å | 434 / 999 |

So a substantial part of the parent umbrella's **251 of 444** is line-list depth error, not
EW error. The count that stands is the one this harness now measures per row.
`verify_feature` already owns the observed-vs-predicted comparison (GF-GHOST); this test is
measurement against measurement.

## The fresh pools

3780–6910 Å, current post-945 line list, 2 359 candidates on each arm, 5 skipped.

| verdict | HARPS `solar_harps_molecfit_corrected` | Kitt Peak (control) |
|---|---|---|
| **in aggregate** | **257** | **318** |
| OVER-PHYSICAL-WIDTH | 1 319 | 1 192 |
| FEATURE-VERIFICATION | 425 | 530 |
| REW-SATURATED | 303 | 302 |
| CONTINUUM-UNPHYSICAL | 38 | 0 |
| BAND-POLICY (3780–3800 Å is near-UV) | 17 | 17 |
| **OVER-IMPLIED-WIDTH** | **0** | **0** |

The integral backstop fires **zero** times: the Doppler ceiling already catches the whole
pathological population, which is the outcome that justifies keeping the backstop hard and
unallowanced rather than tuning it into a filter.

**No impossible width survives.** In-aggregate implied width:

| | median | p95 | max |
|---|---|---|---|
| HARPS | 0.175 Å | 0.402 Å | 0.551 Å |
| Kitt Peak | 0.152 Å | 0.341 Å | 0.426 Å |
| *stale `sol_ew_results_v1`* | *0.69 Å* | — | *31.8 Å* |

The ticket's stated expectation was "order 0.1–0.3 Å at HARPS R". That is what the guarded
pool returns, and it was not aimed at — nothing in the guard knows about it.

## The values — and they are TOO HIGH, which the graded tier explains

Kitt Peak, fresh 304-line guarded pool, aggregated by the median as `build_product` does:

| leg | A(Fe I) | n | stat | syst | vs gold 7.466 |
|---|---|---|---|---|---|
| 1D-LTE (profile-fit) | 7.586 | 304 | 0.0214 | 0.1705 | **+0.120** |
| Engine-A (Bergemann NLTE) | 7.596 | 222 | 0.0235 | 0.1705 | **+0.130** |
| Engine-B (Turbospectrum flux-fit) | 7.503 | 252 | 0.0236 | 0.1700 | +0.037 |

HARPS `solar_harps_molecfit_corrected`: **7.656** (n=247) and **7.666** (n=181) on the same
two EW legs — **+0.19 dex against gold, and ~0.07 dex above Kitt Peak**.

**These numbers are high, and the budget already names why**: `gf scale (UNGRADED)` is the
dominant term on every product, at 0.170 dex, with the verdict *"more lines will NOT help;
fix the source."* The pool sits on the ungraded Kurucz semi-empirical gf floor, and the
next section measures that offset at −0.14 dex — which is most of the excess.

## 🔴 KITT PEAK IS NOT A CONTROL FOR THE VALUE, AND CALLING IT ONE WAS WRONG

An earlier draft of this document said the Kitt Peak arm "reproduces the incumbent 7.586
exactly" and treated that as validation. **It is not, and the reasoning was circular.**
The incumbent VIS number was itself derived from the Kitt Peak atlas on this same ungraded
Kurucz gf scale. Reproducing it demonstrates that the harness is deterministic and that the
width guard did not move the aggregate — worth knowing, and no more than that. A control
has to be something the measurement can be *wrong against*, and this one cannot be.

The same limit applies to the arm comparison. Kitt Peak versus HARPS is a **cross-check of
the measurement chain** — two instruments, one Sun, one line list, one gf scale. Both arms
inherit the identical gf offset, so their agreement can never detect it. What that
comparison *can* see is instrumental and continuum error, and it currently sees **0.07 dex**
of it, which is large and unexplained (see Owed).

**The only externally discriminating comparison in this run is the graded subset**, below.

## 🔴 THE TWO TIERS SEPARATE, AND THE GRADED ONE LANDS ON THE LABORATORY ANCHOR

Split by whether the line's gf is primary-laboratory (`gf_tier == LAB` in `canonical_gf`),
never merged, per RYA-946:

| leg | showcase (GF-LAB) | document (ungraded) | delta |
|---|---|---|---|
| 1D-LTE | **7.463** (n=13) | 7.604 (n=291) | **−0.141** |
| Engine-A | **7.464** (n=12) | 7.608 (n=210) | **−0.145** |
| Engine-B | **7.448** (n=9) | 7.505 (n=243) | −0.057 |

All three showcase values sit within 0.02 dex of **Den Hartog 2014's 7.45** and of this
project's own gold **7.466** — reached by three different treatments on pools they do not
fully share. The document tier sits 0.14 dex above, in incumbent territory.

That is the result: **the pipeline is not producing a wrong answer, the ungraded gf scale
is.** Restrict the pool to lines with a primary laboratory gf and the VIS band lands on
the accepted solar iron abundance. Do not read the 7.586 / 7.656 aggregates as this
band's answer — read them as the ungraded floor, with the offset now measured.

This **independently reproduces RYA-819/831** (a 148-line pool reading 7.586 against a
9-line lab pool reading 7.445) on a *fresh* pool, with a *different* line count, behind a
*new* width guard. And it turns the budget's dominant term from a quoted blanket into a
measurement: `gf scale (UNGRADED)` is carried at 0.170 dex, and the graded-vs-ungraded
offset actually measured here is **−0.14 dex**.

⚠️ **The showcase tier is n=13, sem 0.062**, so the −0.141 offset is ~2.3σ on its own
statistics. What makes it convincing is not that sigma — it is that three treatments, two
arms, and an independent earlier ticket all land in the same place as the laboratory value.
Quote the corroboration, never the sigma alone.

## The laboratory backbone still cannot be measured by EW, and now we know why

57 of the 240 lab-graded Fe I VIS lines reach the fitter at all — the rest are cut upstream
by the `central_depth ∈ [0.05, 0.60]` triage in `scripts/line_accounting_rya709.py`. Of
those 57:

| fate | HARPS | Kitt Peak |
|---|---|---|
| REW-SATURATED | 30 | 24 |
| OVER-PHYSICAL-WIDTH | 15 | 12 |
| FEATURE-VERIFICATION | 6 | 8 |
| **in aggregate** | **6** | **13** |

The dominant loss is **saturation** — the physics, not the plumbing. Laboratory gf papers
target strong clean transitions, which in the Sun are deep and mostly saturated, so the EW
route genuinely cannot measure this population. That is measured confirmation of RYA-958's
own design decision to route the deep graded lines to synthesis, and it is the reason
RYA-955's "8 → ~200" metric was never reachable by any line-list reroute.

## Owed

* 🔴 **The 0.07 dex HARPS-minus-Kitt-Peak gap.** Same Sun, same line list, same gf, same
  fitter, two instruments — they should agree, and 0.07 dex is not scatter (the stat term
  on each is ~0.02). Because both arms carry the identical gf offset, this gap is
  attributable to the measurement chain: continuum placement, the telluric correction, or
  the holdings' differing normalisation contracts. It is the one thing the arm comparison
  is actually able to detect, and it is currently unexplained.
* **The ungraded gf floor is the band's dominant error and it is not fixable here.** The
  aggregates run 0.12-0.19 dex above gold because 290 of 304 lines carry a Kurucz
  semi-empirical gf. Only 57 of the 240 laboratory-graded VIS lines reach the fitter at
  all, and 30 of those saturate. Closing this needs the graded pool to grow, which is
  RYA-958's synthesis leg and RYA-946's sweep, not another EW re-measurement.

* **`sigma_max = 0.40 Å` is still the optimiser's bound on both fitters.** This ticket
  convicts the fits that reach it; it does not tighten it. Tightening would let the
  optimiser find a physical solution for some of the 1 319 instead of quarantining them,
  and would change every EW where σ landed between 0.083 and 0.40 — a re-measurement
  decision with its own control, not a side effect of a guard.
* The 3780–3800 Å lines are refused by band policy (near-UV bans profile-fit) and by the
  holding's own 3782.6 Å blue edge. They belong to RYA-960.
