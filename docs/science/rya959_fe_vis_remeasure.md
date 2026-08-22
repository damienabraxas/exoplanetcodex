# RYA-959 — Fe I VIS re-measured, and the width ceiling that was never there

**Status:** measurement landed; abundance leg in progress.
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

* **`sigma_max = 0.40 Å` is still the optimiser's bound on both fitters.** This ticket
  convicts the fits that reach it; it does not tighten it. Tightening would let the
  optimiser find a physical solution for some of the 1 319 instead of quarantining them,
  and would change every EW where σ landed between 0.083 and 0.40 — a re-measurement
  decision with its own control, not a side effect of a guard.
* The 3780–3800 Å lines are refused by band policy (near-UV bans profile-fit) and by the
  holding's own 3782.6 Å blue edge. They belong to RYA-960.
