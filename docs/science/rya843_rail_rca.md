# RYA-843 — why a synthesis flux-fit runs to its abundance bound

**Status:** RCA complete. **No product published.** The criterion and the threshold are
RYA-847's; this document is the evidence RYA-847 is built on.

Reproduce with `scripts/rya843_rail_rca.py` (`--step atlas | curves | report`).
Artifacts in `data/results/rya843/`.

---

## The question, and why the obvious fix was refused

RYA-837's first NIR run reported **A(Fe I) = 7.588 ± 1.977** (n=40). The 1.977 was not
scatter: several lines had fitted to ≈12.49, which is `a_hi = A_solar + 5` — the search
rail. `_fit_synth_flux`'s `edge_pinned` detector should have caught them and did not,
because its tolerance is `1e-2` and three sat 0.017 away.

Widening the tolerance would have cleaned the product in one line. RYA-843 refused it:
`edge_pinned` is a shared fitter used by every band, and choosing a cut because it makes
*this* product look tidy is RYA-161 tuning. The ticket ordered **RCA first, detector
second** — a line the fitter runs to its bound is more likely a line-selection failure
than a fitter failure (RYA-842: line selection dominates).

That ordering was right. The fitter turned out not to be the problem, and the detector
turned out to be the wrong instrument entirely.

---

## 1. Five of six were telluric — a selection defect

| λ (Å) | A_fit | observed core depth | verdict |
|---|---|---|---|
| 11149.262 | 12.488 | 0.993 | H₂O 11120–11560 |
| 11251.115 | 12.469 | 1.000 | H₂O 11120–11560 |
| 11298.859 | 12.486 | **1.001** | flux **below zero** |
| 11374.078 | 12.478 | 0.987 | H₂O 11120–11560 |
| 11439.124 | 10.324 | 0.752 | H₂O 11120–11560 |

An observed depth of 1.000 means the atlas flux is zero in the line core; 1.001 means it
is negative. There is no Fe profile there to fit, so the optimizer runs the abundance to
its bound.

The project already knew: `telluric_reason` says the flux there is not stellar, and the
**EW route calls it**, skipping 29 NIR lines. The synthesis route inherited
`select_lines`, which knows nothing about tellurics. The mask simply never reached this
path.

**Control** (so the mask is not merely excluding everything): 6 of the railed set are
inside a telluric window; **1 of 30** well-behaved lines is.

Fixed at selection, which makes the detector tolerance *moot* for these rather than
tuning a threshold to hide them.

    7.588 ± 1.977   n=40   as first run
    7.524 ± 1.561   n=36   + RYA-759's status=='ok' filter restored (RYA-837)
    7.496 ± 1.139   n=32   + telluric mask at selection

The value converges toward the solar anchor as non-measurements are removed. That is a
**consequence** of removing garbage, not a target (RYA-844 firewall).

---

## 2. It is not the gf

Linelist `loggf` equals `canonical_gf.log_gf` to **0.000** for all five survivors *and*
all four controls, every one `VALD3 / single_source / hfs1`. RYA-834's redward extension
substituted nothing here. Closed.

## 3. It is not saturation — the standing hypothesis, refuted

These are the band's strongest lines (theoretical depth 0.296–0.583 against a control
median of 0.288), so a saturated core that cannot be deepened was the natural
explanation. The χ² curves refute it: **the synthetic core reaches the observed depth at
a sane abundance and χ² keeps falling well past that point.**

11593.588 — observed core 0.588, synthetic core 0.595 at **A = 8.996** — and χ² goes on
dropping to A = 10.496. Same shape on all five.

## 4. It is the window baseline

At ±1.4 Å the fit window holds 150–270 pixels and the line core a handful, so χ² is
dominated by how well the **baseline** is reproduced, not the line. Every railing line
sits in a window whose Kitt Peak median flux is depressed against a synthesis normalised
to unity, and A(Fe) is the only free parameter available to close the gap:

| λ (Å) | KP window median | A_fit |
|---|---|---|
| 11593.588 | 0.732 | 10.559 |
| 12638.703 | 0.857 | 12.486 |
| 12648.741 | 0.919 | 12.489 |
| 11973.046 | 0.927 | 10.056 |
| 11689.972 | 0.947 | 12.486 |
| *clean controls* | 0.966–0.996 | 7.44–7.45 |

**The rail is a normalisation failure wearing a fitter's clothes.**

⚠️ **Depression correlates but does not separate.** 11119.795 (0.872) and 11572.523
(0.866) are *more* depressed than three of the five railers and fit to plausible values.
The tidy version of this finding is false and is recorded as false.

⚠️ **The IAG arm cannot referee this band.** `telluric_policy` justifies its band list on
KP-vs-IAG, but the staged IAG atlas covers **5001.1–11083.4 Å** — below every one of the
five. Scoped to the artifact inspected, not to Baker+2020 (RYA-833). Whether a NIR
extension exists upstream is unchecked.

---

## 5. 🔴 The defect is bigger than the rail

Constraint measured as the fractional χ² rise from the minimum out to the ends of the
8-dex bracket — *did the objective pin A(Fe) at all?*

| λ (Å) | weaker side | σ_A (dex) | red_χ² | A_fit | in aggregate? |
|---|---|---|---|---|---|
| 12638.703 | **0.0 %** (min on the bound) | 0.753 | 392 | 12.486 | **yes** |
| 11689.972 | 0.0 % (min on the bound) | 0.958 | 244 | 12.486 | no |
| 12648.741 | 0.0 % (min on the bound) | 0.841 | 423 | 12.489 | no |
| 11572.523 | **1.4 %** | 0.543 | 687 | 7.979 | **yes** |
| 11593.588 | 2.2 % | 0.443 | 1226 | 10.559 | **yes** |
| 11119.795 | **2.2 %** | 0.766 | 1173 | 7.833 | **yes** |
| 11973.046 | 8.0 % | 0.267 | 222 | 10.056 | **yes** |
| 10145.561 | 46.7 % | 0.102 | 72 | 7.452 | yes |
| 12342.916 | 557.8 % | 0.028 | 2.6 | 7.438 | yes |

**11119.795 and 11572.523 entered this list as well-behaved controls. They are not
measurements.** Their χ² moves 1.4 % and 2.2 % across *eight dex of iron*; they landed at
7.833 and 7.979 by luck and were accepted.

So the defect is not "railed fits inflate σ". It is:

> **Unconstrained fits are accepted, and only the ones that land somewhere implausible
> are visible.** The rail is the loud half; the quiet half passes silently and biases
> nothing visibly at all.

This is what RYA-847 was filed to fix, and it is why the fix has to be pipeline-wide
rather than a NIR patch.

---

## 6. What this says about the criterion — and what it refuses to decide

`edge_pinned` is the wrong instrument: three of these never turned over, and two turned
over meaningfully but in the wrong place. Distance-to-bound cannot express either. Worth
recording how sharp that is — `a_hi` is **12.496** exactly, and two lines both *printing*
12.486 fell on opposite sides of `< 1e-2`, one 0.0100 from the bound and one a hair
under. **The verdict was decided in the 4th decimal of a rounded number.**

Three candidate replacements were measured. **None is adopted here.**

- **Fractional χ² rise across the bracket.** Separates cleanly on these nine. But the
  bracket is `[A_solar−3, A_solar+5]`, so a minimum near A = 10 is close to the +5 end
  and has little room left to rise. The metric is part *sharpness of the minimum* and
  part *proximity to the bracket edge*, and the second is a plausibility prior tied to
  A_solar entering through the back door. It is also caller-dependent, so it does not
  transfer between handlers that bracket differently.
- **σ_A from the rescaled curvature.** Bracket-free and in **dex**, so it is arguable on
  physics and transfers across callers. Δχ² is a *difference*, so the depressed-baseline
  offset cancels — which is precisely why it crosses bands where absolute red_χ² cannot.
  But it measures **precision, not accuracy**: it accepts 11973.046 at A = 10.056, whose
  minimum is genuinely sharp and genuinely wrong.
- **Absolute red_χ².** `pipeline/measure/synthesis.py` already refuses a line on
  `red_chi2 >= SYNTH_CHI2_GATE` (10.0, ratified by RYA-342), *before* it checks
  `edge_pinned`. **The band-product synthesis route never calls that handler**,
  reimplements accept/reject as `status != 'ok'`, and throws `red_chi2` away —
  `LineMeasurement` has no χ² field. That bypass is real and RYA-847 closes it. But the
  **constant does not transfer**: the good line 10145.561 fits at red_χ² = **72**, so
  gating at 10.0 here would refuse good lines.

The likely answer is two conditions — σ_A for *was A(Fe) determined*, red_χ² for *does
the model describe the data* — but **nine lines cannot settle that, and RYA-161 forbids
setting a cut from one product.** RYA-847 sets it from a sweep across every synthesis
band.

### Owed, found on the way

- **red-optical and NIR carry no continuum systematic.** Only the near-UV policy says
  "pseudo-continuum", so `error_budget.build()` adds the term there and nowhere else;
  these two bands get only the 0.030 telluric residual. §4 is the evidence they need one.
  The honest size is not 0.100 inherited from the near-UV — it has to be derived for
  these bands (RYA-161).
- **11119.795** sits 0.205 Å below the `TELLURIC_BANDS` H₂O edge at 11120.0, in a window
  at median flux 0.872, with no atlas able to referee it. RYA-847 item 7 excludes it
  conservatively.
- **Three implementations of the same accept/reject** now exist:
  `abundances_derive._fit_synth_flux`, a duplicate in `cno_synthesis`, and the
  translation in `measure/synthesis.py`. `cno_synthesis` additionally computes a
  curvature σ and **clips it to [0, 1]**, so an unconstrained fit reports `sigma_fit =
  1.000` — a number that reads like a measured 1 dex. Nothing gates on it.
