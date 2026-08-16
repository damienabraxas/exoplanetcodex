# RYA-841 — the near-UV pseudo-continuum, and what actually limits the band

**Star:** Sun, Kitt Peak flux atlas. **Band:** 3000–3780 Å, synthesis-only.
**Status:** investigation. Nothing here changes a product. **Not merged — Ryan reviews.**

```
python3 scripts/rya841_nearuv_continuum_sensitivity.py     # the lever, dA/d(delta)
python3 scripts/rya841_nearuv_continuum_methods.py         # placement disagreement
python3 scripts/rya841_nearuv_scatter_driver.py            # what drives the scatter
```

---

## The headline

RYA-836 made the pseudo-continuum the dominant near-UV systematic. **It has never been
measured, the arithmetic that produced it contains a unit error, and it is added twice.**

---

## 0. 🔴 Where the 0.147 comes from: the term is counted twice

`error_budget.build()` **already adds** the pseudo-continuum term for this band — the
near-UV `BandPolicy.continuum_treatment` reads *"pseudo-continuum only"*, which fires the
`if "pseudo" in ...` branch at `error_budget.py:170`. Both near-UV routes then add it again
in quadrature:

```
scripts/derive_band_products.py:263        syst = np.hypot(syst, NEARUV_PSEUDO_CONTINUUM_DEX)
scripts/rya836_nearuv_lab_gf_subpool.py    syst = np.hypot(syst, PSEUDO_CONTINUUM_DEX)
```

| cell | `build()` syst | reported | correct |
|---|---|---|---|
| **RYA-832 near-UV (product of record)** | 0.1972 | **0.2211** | **0.1972** |
| **RYA-836 lab-gf sub-pool** | 0.1081 | **0.1472** | **0.1081** |

Both reconstruct to the fourth decimal, and the `stat` values (0.0653, 0.0849) match the
published cells too, so this is the arithmetic that ran and not a lookalike. Only the
near-UV policy fires the branch — VIS, red-optical and NIR do not — so the defect is
confined to these two cells.

**This is mine, from RYA-832, and RYA-836 inherited it.** The comment at the call site says
the term *"is NOT in the line scatter, so it is added in quadrature"* — which is true, and
is not the question. It was already in the **budget**. Worse, the RYA-832 test I wrote
(`test_the_route_carries_the_759_pseudo_continuum_systematic`) asserts the constant is
0.100 and so **pinned the double-count instead of catching it**.

**Direction of the correction:** RYA-836 claimed the systematic fell 0.221 → 0.147 when the
pool moved to primary lab gf. Corrected, it fell **0.197 → 0.108**. The improvement is
*larger* than reported, and the "dominant term flipped from gf to pseudo-continuum" claim is
unaffected (`dominant` is computed from the budget's own terms, before the stray `hypot`).
**RYA-836's qualitative conclusion survives; its numbers do not.**

The fix is a two-line deletion, but it changes published product uncertainties, so it is
flagged here rather than applied.

---

## 1. The 0.100 was never derived, and the lever it assumes is wrong

`pipeline/error_budget.py:170-176` reasons:

> *"near-UV median flux 0.283–0.805 means the normalisation itself is uncertain at the
> ~10% level in the worst windows"* → `Term("pseudo-continuum", 0.10, ...)`

🔴 **That step equates ~10% in flux with 0.10 in dex.** The conversion dA/dδ appears
**nowhere in the repo**, so the budget assumes, without saying so, that **dA/dδ = 1.0**.

There is a second problem in the same sentence: *how far below unity the spectrum sits* is
**blanketing depth**, not **placement uncertainty**. The depression is largely known and is
modelled by the synthesis; what is uncertain is the residual after modelling it.

### Measured

Each line re-fitted with the observed flux renormalised as if the continuum had been
misplaced by a known δ (`f_used = f_atlas / (1 + δ)`), on a symmetric grid, everything else
held fixed:

| | dex per unit fractional continuum error |
|---|---|
| **measured median \|dA/dδ\|** | **2.42** |
| 16–84 percentile | 0.91 – 4.68 |
| **assumed by the budget** | **1.00** |

**The budget understates the lever by about 2.4×.** What the term becomes, for a range of
σ_δ:

| σ_δ | implied term |
|---|---|
| 1% | 0.024 dex |
| 2% | 0.048 dex |
| 5% | 0.121 dex |
| 10% | 0.242 dex |

Inverting it: **the 0.100 dex term implicitly asserts the near-UV continuum is misplaced by
≈4.1%.** That is a checkable claim about the Kitt Peak normalisation, and nobody has checked
it.

⚠️ **The lever is measured on a ±4% grid and is quoted only there.** That brackets the ~4%
misplacement the term implies, which is what it is for. It does not licence extrapolation to
the ~11% atlas-vs-local gap in §3 — those are compared by re-fitting, not by multiplying a
slope.

### 🔴 The response is not well defined for a third of the band

**12 of 40 lines are non-linear** (r² ≤ 0.98), and several respond in the **physically
impossible direction** — a higher continuum giving a *lower* abundance. Those are the
windows with χ²_r in the hundreds. Two lines return slopes above +36 with r² of 0.53 and
0.82, which are not derivatives of anything.

A single continuum systematic applied uniformly to every near-UV line describes neither
group.

---

## 2. σ_δ is missing, and the one diagnostic in the repo does not supply it

The term is `dA/dδ × σ_δ`. The lever is now measured; **σ_δ is not.**

`cont_ratio` (per-window synth/obs median ratio, RYA-759) looks like the right quantity —
median 1.019, MAD 0.068 — but it fails its own consistency check:

| | slope |
|---|---|
| across-line, dA/d(cont_ratio − 1) | **+0.49 ± 0.38** (n=40, r=+0.21) |
| controlled per-line, dA/dδ | **+2.42** |

Those disagree at roughly 5σ. If `cont_ratio` tracked a real continuum misplacement the two
would be the same physical response and would match. They do not — so `cont_ratio` is
dominated by **model** error in the window (missing or wrong blends), not normalisation
error, and **its 6.8% MAD must not be used as σ_δ.**

---

## 3. Continuum placement moves the answer by ~0.29 dex, and χ² cannot referee it

Three placements applied to the same 40 lines, each reduced to a per-line δ so all three go
through the same tested fitter:

| placement | median A | scatter | median δ |
|---|---|---|---|
| **atlas** (Kitt Peak as delivered) | **7.487** | **0.412** | 0 |
| `local` (window's own 95th-pct envelope) | 7.200 | 0.739 | −0.112 |
| `synth` (continuum where the model says) | 7.366 | 0.474 | −0.018 |

*(one railed fit excluded from each alternative — `_fit_synth_flux` is bounded at
`a_solar + 5`, and 3617.318 under `synth` returned A = 12.445, which is the optimiser
leaving the room, not an abundance)*

| pair | median difference |
|---|---|
| atlas − local | **+0.264 dex** |
| atlas − synth | **+0.020 dex** |

**Spread of the median across placements: 0.287 dex — about 2.9× the assumed 0.100.**

### ⚠️ A correction to my own first read

From the first line to complete (3000.468: atlas χ²_r 88 → local 229) I concluded the local
placement was *refuted by fit quality*. **That is wrong.** Across the pool:

| placement | median χ²_r |
|---|---|
| atlas | 118.0 |
| local | 124.3 |
| synth | 108.3 |

**`local` is worse than `atlas` on only 16 of 40 lines** — i.e. better on 24. Fit quality
does not discriminate between placements at all. Whether a placement is rejected is a claim
about the pool, not about whichever line ran first.

### What does discriminate

Two things, and neither is χ²:

1. **Pool scatter mildly prefers the atlas** (0.412 against 0.739 and 0.474). This is a
   non-circular argument — line-to-line scatter does not know which normalisation is
   "correct".
2. **The atlas and model-based placements agree in the median to 0.020 dex**, while the
   naive envelope sits 0.264 dex away. Two of three routes converge and the outlier is the
   method that has no line-free pixels to work with.

So the honest statement is **not** "the systematic is 0.287 dex". It is: *if the atlas
normalisation is accepted — and its tightest-pool and median-agreement give reason to — then
the term should reflect **the atlas's own** normalisation uncertainty, which is unmeasured.
If the placements were treated as equally defensible, the systematic would be ~0.29 dex,
nearly 3× the current term.* Either way, **the 0.100 is not derived from anything and is not
demonstrably conservative.**

The local envelope sits a median **11% below** the atlas continuum (range 4%–44%), which is
the quantitative form of what `band_policy.py` says qualitatively.

---

## 4. What drives the line-to-line scatter: nothing that was measured

RYA-836 attributed 0.238 dex of scatter to *"line selection"*. That is the difference
restated, not a cause. Every candidate cause was tested and **all of them fail.**

### 🔴 The selector does not do what its name suggests

`select_lines(min_sep_A=4.0)` applies the separation test **only between lines it has
already selected** (`rya759_nearuv_fe_product.py:127-134`). Its own docstring says so —
*"keeps the set spread across the band"*. It is a **sampling** rule, not a blend rejection:

| | |
|---|---|
| minimum separation between *selected* lines | 4.168 Å |
| median distance to the nearest *actual* neighbour | **0.0095 Å** |
| median other lines inside the ±0.40 Å fit window | **53** |

### No near-UV window is dominated by its own line

Window-wide, the target line supplies a **median 8%** of the absorption in its own fit window
(max 31%). Over a ±0.05 Å core it is 67%, with a median of 4 other lines inside.
*"Just use the clean lines"* is not available in this band, because there are none.

### Every correlation is null, in two independent pools

Seven metrics — window dominance, core dominance, blend depth, core count, window count,
nearest-neighbour separation, own depth — against |A − median|: **|ρ| ≤ 0.24 everywhere, in
both pools.** This extends RYA-759's *"correlates with nothing"* from raw line properties to
blend-aware ones.

### ⚠️ The one cut that looks like a win is refuted by its control

| depth_frac_core ≥ | RYA-832 pool | RYA-836 lab pool |
|---|---|---|
| 0.00 | 0.412 ± 0.047 | 0.651 ± 0.060 |
| **0.50** | **0.321 ± 0.047** (1.4σ) | **0.744 ± 0.089** (0.9σ) |

Cutting the 832 pool on core dominance takes its scatter 0.412 → 0.321, which looks like a
tightening worth having. **The same physical cut makes the lab pool worse.** A real criterion
cannot do both, so the apparent gain is sampling noise — and RYA-161 forbids taking it.

### Saturation looks like the explanation and is not

The lab pool sits mostly above the 832 selector's depth ceiling (median depth 0.944 vs
`DEPTH_CEIL = 0.90`), and its saturated lines come out **+0.117 dex** high — tempting,
because it would explain why the two matrix cells differ.

**Tested: bootstrap 95% CI [−0.711, +0.439], permutation p = 0.28.** Not significant.

### The pools are effectively disjoint — RYA-836's owed item 2, answered

**Only 1 of 60** lab lines passes the 832 selection. The labs and the clean-line criterion
barely overlap, so the 0.238 dex is a difference between two nearly disjoint samples and
cannot be decomposed line by line.

### Wavelength: ruled out as distribution, unexplained as trend

The pools sit at the same median wavelength (3463 vs 3451 Å) with the same blue fraction
(0.35 vs 0.38). But their **binned scatter runs in opposite directions**:

| bin | RYA-832 | RYA-836 lab |
|---|---|---|
| 3000–3300 Å | 0.531 | 0.459 |
| 3300–3600 Å | 0.390 | 0.688 |
| 3600–3780 Å | **0.247** | **0.817** |

The depth-selected pool gets *tighter* redward; the lab pool gets *worse*. Flagged, not
explained.

---

## 5. The honest tightening ceiling

**From line selection: zero.** No measured per-line physical property predicts the residual,
no defensible subset exists, and the one cut that appears to tighten is contradicted by its
control. Any near-UV number that looks tighter because lines were dropped is an artifact of
the dropping.

**From the continuum: the term is wrong in construction, and its size is unknown.** It is
double-counted (§0), its lever is understated 2.4× (§1), and σ_δ has never been measured
(§2). Correcting only the double-count moves the product of record's systematic
0.2211 → 0.1972; correcting the lever needs σ_δ first.

**What would actually move the band**, in order:

1. **Remove the double-count.** Two lines, no new physics, and it is simply wrong today.
2. **Measure σ_δ** — the Kitt Peak atlas's own normalisation uncertainty in 3000–3780 Å.
   This is the single missing number; with the lever now known, the term follows from it.
3. **Treat the non-linear responders separately.** A term applied per line to lines that do
   not respond to the continuum is not a continuum term.
4. Not more gf work, and not more line selection. Both are measured dead ends
   (RYA-822, RYA-824, RYA-836, and §4 above).

---

## Owed / found on the way

1. ⚠️ **RYA-759's gf provenance string is stale, and 4 lines did move.** It records
   *"canonical_gf.csv starts at 3780.0 Å … NOT available"*; **RYA-822 (`7a84c77`) extended it
   to 3000.003 Å**, 21,279 rows below 3780, and the committed
   `nearuv_fe_product_FINAL.json` still carries the old claim.

   The δ=0 control quantifies the consequence: **36 of 40 lines reproduce exactly, 4 moved**
   (−0.051, −0.016, −0.012, +0.012 dex). The **median is unchanged** — 7.4870 → 7.4875,
   which is exactly RYA-832's published 7.488 — so the product did not move, but "nothing
   changed" is not quite the right statement and the provenance string is now false.
2. **The opposite wavelength trend between the two pools** (§4) is unexplained.
3. **Fe I 3026.056** still wants adjudication (carried from RYA-836).
