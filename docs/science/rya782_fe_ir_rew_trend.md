# RYA-782 — the Fe IR "+0.384 REW slope" is outlier leverage, not a curve-of-growth systematic

Band: Fe I, 6910–9199.9 Å, Kitt Peak solar atlas, 1D-LTE per-line pool (n = 99 with a
known gf source; 35 laboratory-tier).
Generator: `scripts/rya782_rew_trend.py` · artifacts: `FeI_IR_rew_trend.csv`,
`FeI_IR_rew_trend_summary.json`.

## The filed coefficient reproduces exactly

`A ~ EP + REW + is_loose` over the pool returns REW **+0.3838 (t = 2.55)**, EP +0.0167
(t = 0.47), is_loose +0.2938 (t = 3.97), matching RYA-760 to four decimals, and every
per-source median matches its table. Nothing below rests on a different pool.

## The ticket's framing is refuted by its own pool

RYA-782 states the finding is *"a gentle systematic SLOPE across the whole Fe IR pool, not
a single catastrophic bad line."* It is the other way round.

| cut | n | REW coefficient | residual sd |
|---|---|---|---|
| as filed | 99 | **+0.384 (t = 2.55)** | 0.311 |
| 3σ clipped (7 lines) | 92 | +0.127 (t = 1.39) | 0.181 |
| K07 removed | 52 | +0.081 (t = 0.51) | 0.192 |
| K07 removed + 3σ clipped | 49 | −0.051 (t = −0.41) | 0.150 |
| **laboratory pool** (= the delivered product) | 35 | +0.232 (t = 1.16) | 0.163 |
| laboratory pool, 3σ clipped | 34 | **+0.008 (t = 0.05)** | 0.133 |

The **marginal** REW slope over the whole pool is +0.103 (t = 0.71) — not significant. The
+0.384 exists only *conditional on* `is_loose`, a binary that lumps four sources of three
accuracies (FMW/GESB82c 0.08, MRW 0.08, K07 0.20) into one level. K07 alone is 47 of the 99
lines. A gf error that correlates with line strength — what a semi-empirical source is
expected to have — lands in the REW term, not the source term.

Refit with per-source dummies: is_K07 **+0.289 (t = 3.62)**, is_FMW +0.428 (t = 3.78),
is_MRW +0.013 (t = 0.07 — MRW behaves like the laboratory sources, confirming RYA-760).
REW alone explains **0.5 %** of the pool variance, and **0.0 %** after clipping; the source
dummies explain 11.7 % → 14.3 %.

## The three proposed remedies, as tests

**(b) strong-line saturation cut — already in force.** `pipeline/band_products.py`
`REW_SATURATION_CEILING = -4.90` is enforced, and the pool's REW range is
[−5.942, **−4.900**] — it tops out exactly at the ceiling. There is no saturated tail left
to cut. Cutting further *raises* the coefficient (EW < 100/80/60 mÅ → +0.378/+0.402/+0.472),
because it strips weak-line ballast rather than a saturated excess.

**(c) blend re-check — negative.** A Boltzmann-weighted strongest-neighbour ratio added as a
regressor carries nothing: −0.0071 (t = −0.66), and REW barely moves (+0.317 → +0.334). Note
the high-excitation Cr II / Fe II neighbours that look alarming in a raw ±0.15 Å dump are at
EP 11–13 eV and are negligible in the solar photosphere; the weighted metric correctly
discounts them.

**(a) microturbulence — refuted by locality, no synthesis run required.** A vmic error acts
through the *saturated* end of the curve of growth, so its slope must concentrate near the
knee and vanish in the clearly-linear regime. The pool shows the **opposite**:

| REW window | n | slope | t |
|---|---|---|---|
| clearly linear (≤ −5.2) | 49 | **+0.624** | +2.13 |
| mid (−5.2 … −5.0) | 35 | +0.734 | +0.98 |
| near knee (−5.1 … −4.9) | 32 | **−0.283** | −0.30 |

The trend lives where microturbulence has least leverage and reverses sign where it has
most. A vmic re-solve against this would be fitting a gf artifact — and, because the pool
sits below the saturation ceiling, would have little to act on in any case.

## What is actually there: seven catastrophic lines

| λ (Å) | src | acc | A | robust z | EW (mÅ) | centroid mismatch | blend |
|---|---|---|---|---|---|---|---|
| 8024.543 | K07 | 0.20 | **9.678** | +11.3 | 49.3 | **46.6 mÅ** | −6.53 |
| 7052.715 | K07 | 0.20 | 8.711 | +6.0 | 72.8 | 1.0 mÅ | −2.37 |
| 7810.814 | K07 | 0.20 | 8.641 | +5.6 | 70.3 | 1.0 mÅ | none |
| 7190.122 | FMW | 0.08 | 8.289 | +3.7 | 42.5 | 0.8 mÅ | −2.12 |
| 7120.021 | K07 | 0.20 | 8.254 | +3.5 | 38.7 | 0.8 mÅ | −2.08 |
| 7330.137 | FMW | 0.08 | 8.184 | +3.1 | 37.9 | 0.1 mÅ | −1.72 |
| 7107.459 | BWL | 0.02 | **6.979** | −3.5 | 30.9 | 0.6 mÅ | none |

**8024.543 is misidentified**, not blended: the measured feature is 46.6 mÅ from the nearest
catalogued Fe I, and the transition credited to it (log gf −4.746, EP 5.879 eV) cannot
produce 49.3 mÅ in the Sun. An identification cut removes this line and only this line; the
slope survives it (+0.364, t = 3.01), so mis-ID is one line, not the mechanism — **K07 is**.

**7190.122 and 7330.137 are already `QUARANTINED-SCALE-EVIDENCE` under RYA-780.** They are
that ticket's signal, correctly identified; RYA-782 does not re-adjudicate them. The
remaining five are new and are now registered in `data/registry/problem_children.csv`
(quarantined, never culled): four K07 lines plus **7107.459**, the one laboratory-tier
casualty and the only negative outlier — independently flagged by RYA-712 for the only
negative Gerber NLTE delta of 14 lines (−0.1063 against +0.043…+0.104). Two independent
anomalies on one clean line point at the atomic model, not the EW or a blend. Cause **not**
established; carried `owed`, not re-sourced.

## The stated error-budget term (RYA-713)

On the laboratory pool — the gf-homogeneous sample, and the one the product is built from —
the REW slope is **+0.232 ± 0.199 (t = 1.16), consistent with zero**. At n = 35 that is *not
detected*, not *proven absent*, so it is carried as a bound rather than a correction:

- best-estimate weak-to-strong swing across the pool's 0.52 dex REW span: **0.121 dex**
- 95 % upper bound on a hidden slope +0.622 → swing **≤ 0.326 dex**
- contribution to per-line scatter: **0.034 dex**

**Recommendation: carry 0.034 dex as the stated REW-trend term in the Fe IR error budget,
and do not apply a slope correction.** Correcting a slope that is consistent with zero on
the gf-homogeneous pool, and demonstrably driven by an excluded gf tier on the full pool,
would import a K07 artifact into a laboratory-tier product.

No coefficient here was chosen to move A(Fe) toward the optical anchor; the anchor appears
only descriptively (laboratory pool median 7.551, **+0.085** vs 7.466). RYA-161 firewall
respected.

## Consequence for RYA-760's reading

RYA-760 noted the REW trend "runs in the direction that PARTLY MASKS the source offset (why
raw +0.158 < controlled +0.294)". That mechanism does not survive: the masking term is
K07-driven, and K07 is neither the FMW population nor in the product. The FMW source offset
itself is unaffected by this finding — it rests on RYA-780's primary-scale adjudication, not
on the REW term.
