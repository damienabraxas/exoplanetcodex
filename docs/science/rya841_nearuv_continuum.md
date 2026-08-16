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
measured, and the arithmetic that produced it contains a unit error.**

`pipeline/error_budget.py:170-176` adds the term whenever a band policy mentions a
pseudo-continuum, with this reasoning:

> *"near-UV median flux 0.283–0.805 means the normalisation itself is uncertain at the
> ~10% level in the worst windows"* → `Term("pseudo-continuum", 0.10, ...)`

🔴 **That step equates ~10% in flux with 0.10 in dex.** They are different units, and the
conversion between them — dA/d(continuum) — appears **nowhere in the repo**. The budget
therefore assumes, without saying so, that **dA/dδ = 1.0 dex per unit fractional continuum
error**.

There is a second problem inside the same sentence. *How far below unity the spectrum
sits* (median flux 0.283–0.805) is **blanketing depth**, not **placement uncertainty**.
The depression is largely known and is modelled by the synthesis; what is uncertain is the
residual after modelling it. The 0.10 uses the first as if it were the second.

---

## 0. 🔴 Where the 0.147 comes from: the term is counted twice

The ticket asks where the 0.147 comes from. It comes from **adding the 0.100 dex
pseudo-continuum term twice**, and the same defect is in the near-UV product of record.

`error_budget.build()` already adds the term for this band — the near-UV
`BandPolicy.continuum_treatment` reads *"pseudo-continuum only"*, which fires the
`if "pseudo" in ...` branch at `error_budget.py:170`. Both near-UV routes then add it
again in quadrature:

```
scripts/derive_band_products.py:263        syst = np.hypot(syst, NEARUV_PSEUDO_CONTINUUM_DEX)
scripts/rya836_nearuv_lab_gf_subpool.py    syst = np.hypot(syst, PSEUDO_CONTINUUM_DEX)
```

Reconstructed exactly from the published cells:

| cell | `build()` syst | reported | correct |
|---|---|---|---|
| **RYA-832 near-UV (product of record)** | 0.1972 | **0.2211** | **0.1972** |
| **RYA-836 lab-gf sub-pool** | 0.1081 | **0.1472** | **0.1081** |

Both reconstruct to the fourth decimal, and the `stat` values (0.0653, 0.0849) match the
published cells too, so this is the arithmetic that ran and not a lookalike.

**This is mine, from RYA-832, and RYA-836 inherited it.** The comment at the call site
says the term *"is NOT in the line scatter, so it is added in quadrature"* — which is true,
and is not the question. It was already in the **budget**. Worse, the RYA-832 test I wrote
(`test_the_route_carries_the_759_pseudo_continuum_systematic`) asserts the constant is
0.100 and so **pinned the double-count instead of catching it**.

**Direction of the correction:** RYA-836 claimed the systematic fell 0.221 → 0.147 when the
pool moved to primary lab gf. Corrected, it fell **0.197 → 0.108**. The improvement is
*larger* than reported, and the "dominant term flipped from gf to pseudo-continuum" claim
is unaffected (`dominant` is computed from the budget's own terms, before the stray
`hypot`). **RYA-836's qualitative conclusion survives; its numbers do not.**

The fix is a two-line deletion, but it changes published product uncertainties, so it is
flagged here rather than applied — that call is Ryan's.

---

## 1. The lever, measured

The only way to get dA/dδ is a controlled perturbation, so each line is re-fit with the
observed flux renormalised as if the continuum had been misplaced by a known δ:

```
f_used = f_atlas / (1 + delta)      delta = fractional error in the adopted continuum
```

δ > 0 means the continuum was placed too high, lines look deeper, and A must rise.
Everything else — atmosphere, linelist, window, bounds, broadening — is held fixed, which
is what makes the response attributable to the continuum alone.

**Positive control:** the δ = 0 column reproduces RYA-759's published per-line values
exactly (8.307, 6.837, 7.458, … line for line). The harness differs from the product of
record only by the perturbation.

⚠️ **The lever is measured on a ±4% grid and is quoted only there.** That brackets the ~6%
misplacement the 0.100 dex term implies, which is what it is for. It does **not** licence
extrapolation to the ~16% gap between the atlas and local-envelope placements (§3) — those
are compared by re-fitting, not by multiplying a slope.

*(Final lever numbers are filled in when the 200-fit run completes; the partial result at
7 lines is a median |dA/dδ| of ~1.2, and ~1.6 restricted to lines that respond linearly.)*

### The part that was not expected

**The response is only well defined for some of the band.** Lines that respond linearly
(r² > 0.98) all respond *positively*, with slopes ≈ +1.3 to +1.9. The rest respond
non-linearly, and several respond in the **physically impossible direction** (dA/dδ < 0 —
placing the continuum higher making the abundance *lower*). Those are the windows with
χ²_r in the hundreds.

That matters more than the slope itself: for a substantial fraction of near-UV lines, the
fit is not tracking the continuum at all, and a single "continuum systematic" applied to
every line describes none of them.

---

## 2. σ_δ is missing, and the one diagnostic in the repo does not supply it

The term is `dA/dδ × σ_δ`. The lever is now measured; **σ_δ — how badly the near-UV
continuum is actually placed — is not, and cannot be taken from `cont_ratio`.**

`cont_ratio` (per-window synth/obs median ratio, RYA-759) looks like the right quantity:
median 1.019, MAD 0.068. But it fails its own consistency check:

| | slope |
|---|---|
| across-line, dA/d(cont_ratio − 1) | **+0.49 ± 0.38** (n=40, r=+0.21) |
| controlled per-line, dA/dδ | **≈ +1.7** |

Those disagree at ~3σ. If `cont_ratio` tracked a real continuum misplacement the two would
be the same physical response and would match. They do not — so `cont_ratio` is dominated
by **model** error in the window (missing or wrong blends), not normalisation error, and
**its 6.8% MAD must not be used as σ_δ.**

### What the 0.100 implicitly asserts

Inverting the measured lever: the 0.100 dex term asserts the near-UV continuum is
misplaced by **δ ≈ 6%**. That is a checkable claim about the Kitt Peak atlas normalisation,
and **nobody has checked it.** It is not obviously wrong — it is simply unverified, and it
was arrived at by a route that could not have verified it.

---

## 3. Continuum placement: one method is refuted, not merely less accurate

Three placements applied to the same lines:

| placement | δ | what it is |
|---|---|---|
| `atlas` | 0 | the Kitt Peak normalisation as delivered |
| `local` | E − 1 | the window's own upper envelope (95th pct over ±2 Å) |
| `synth` | 1/cont_ratio − 1 | put the continuum where the model says it is |

**The placements disagree by a measured amount, and it brackets the assumed 0.100:**

| pair | median difference |
|---|---|
| atlas − local | **+0.132 dex** |
| atlas − synth | **+0.075 dex** |

The local envelope sits a median **16% below** the atlas continuum (range 4%–44%), which is
the quantitative form of what `band_policy.py` says qualitatively — the true continuum is
never observed, and how far the envelope falls short varies by a factor of ten across the
band.

⚠️ **A correction to my own first read.** From the first line to complete (3000.468: atlas
χ²_r 88 → local 229) I concluded the local placement was simply *refuted* by fit quality.
Across more lines that is **wrong** — it is often *better* (3087.4: 891 → 676; 3123.2:
999 → 883). Whether a placement is rejected is a claim about the pool, not about whichever
line ran first, so the report now joins the atlas χ² in and counts how many lines each
placement is worse on.

That the two placements differ by 0.13 dex while neither is clearly rejected by fit quality
is the honest statement of the problem: **the near-UV continuum is genuinely ambiguous at
the tenth-dex level, and the fit cannot referee it.**

⚠️ `synth` is partly circular — normalise to the model, then fit the model — and is
reported for distance only, never as a recommended placement.

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

Window-wide, the target line supplies a **median 8%** of the absorption in its own fit
window (max 31%). Over a ±0.05 Å core it is 67%, with a median of 4 other lines inside.
*"Just use the clean lines"* is not available in this band, because there are none.

### Every correlation is null, in two independent pools

Seven metrics — window dominance, core dominance, blend depth, core count, window count,
nearest-neighbour separation, own depth — against |A − median|:

**|ρ| ≤ 0.24 everywhere, in both pools.** This extends RYA-759's *"correlates with
nothing"* from raw line properties to blend-aware ones.

### ⚠️ The one cut that looks like a win is refuted by its control

| depth_frac_core ≥ | RYA-832 pool | RYA-836 lab pool |
|---|---|---|
| 0.00 | 0.412 ± 0.047 | 0.651 ± 0.060 |
| **0.50** | **0.321 ± 0.047** (1.4σ) | **0.744 ± 0.089** (0.9σ) |

Cutting the 832 pool on core dominance takes its scatter 0.412 → 0.321, which looks like a
tightening worth having. **The same physical cut makes the lab pool worse.** A real
criterion cannot do both, so the apparent gain is sampling noise — and RYA-161 forbids
taking it.

### Saturation looks like the explanation and is not

The lab pool sits mostly above the 832 selector's depth ceiling (median depth 0.944 vs
`DEPTH_CEIL = 0.90`), and its saturated lines come out **+0.117 dex** high — tempting,
because it would explain why the two matrix cells differ.

**Tested: bootstrap 95% CI [−0.711, +0.439], permutation p = 0.28.** Not significant. The
offset is quoted only with its test attached.

### The pools are effectively disjoint — RYA-836's owed item 2, answered

**Only 1 of 60** lab lines passes the 832 selection. The labs and the clean-line criterion
barely overlap, so the 0.238 dex is a difference between two nearly disjoint samples and
cannot be decomposed line by line.

### Wavelength: ruled out as distribution, unexplained as trend

The pools sit at the same median wavelength (3463 vs 3451 Å) with the same blue fraction
(0.35 vs 0.38), so an uneven distribution is not the explanation. But their **binned
scatter runs in opposite directions**:

| bin | RYA-832 | RYA-836 lab |
|---|---|---|
| 3000–3300 Å | 0.531 | 0.459 |
| 3300–3600 Å | 0.390 | 0.688 |
| 3600–3780 Å | **0.247** | **0.817** |

The depth-selected pool gets *tighter* redward; the lab pool gets *worse*. This is flagged,
not explained — whatever separates them is not a line property measured here.

---

## 5. The honest tightening ceiling

**From line selection: zero.** No measured per-line physical property predicts the
residual, no defensible subset exists, and the one cut that appears to tighten is
contradicted by its control. Any near-UV number that looks tighter because lines were
dropped is an artifact of the dropping.

**From the continuum: unknown, and the current term is unvalidated.** The lever is
measured; σ_δ is not, and the diagnostic that exists does not measure it. The 0.100 dex may
turn out to be roughly the right size — but if so that is a coincidence, because the route
to it equated percent with dex.

**What would actually move the band**, in order:

1. **Measure σ_δ** — the Kitt Peak atlas's own normalisation uncertainty in 3000–3780 Å.
   This is the single missing number, and everything else about the term follows from it.
2. **Treat the non-linear responders separately.** A term applied per line to lines that do
   not respond to the continuum is not a continuum term.
3. Not more gf work, and not more line selection. Both are now measured dead ends
   (RYA-822, RYA-824, RYA-836, and §4 above).

---

## Owed / found on the way

1. 🔴 **RYA-759's gf provenance string is stale.** It records *"canonical_gf.csv starts at
   3780.0 Å … NOT available"*; **RYA-822 (`7a84c77`) extended it to 3000.003 Å**, 21,279
   rows below 3780. The committed `nearuv_fe_product_FINAL.json` still carries the old
   claim. The values agree either way (δ=0 reproduces 759 line for line), so nothing
   moved — but the recorded provenance is now false and should be corrected.
2. **The opposite wavelength trend between the two pools** (§4) is unexplained.
3. **Fe I 3026.056** still wants adjudication (carried from RYA-836).
