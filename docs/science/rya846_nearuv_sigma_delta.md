# RYA-846 — σ_δ, and the honest near-UV continuum term

**Star:** Sun. **Band:** 3000–3780 Å. **Status:** measurement. **Not merged — Ryan reviews.**

```
python3 scripts/rya846_nearuv_sigma_delta.py        # sigma_delta + the term + revised cells
python3 scripts/rya846_nonlinear_responders.py      # the 12 lines with no derivative
```

---

## The answer

RYA-841 measured the lever; this measures σ_δ. The honest near-UV continuum systematic is

| | σ_δ | term (= lever × σ_δ) |
|---|---|---|
| **net** (reduction difference, Wallace-internal scatter removed) | **2.79%** | **0.071 dex** |
| **raw** (Kurucz vs Wallace at the lines — an upper bound) | 4.94% | 0.125 dex |
| assumed by the budget today | (≈3.9% implied) | 0.100 dex |

**The assumed 0.100 dex lands inside the measured range.** It was never derived and its
derivation was wrong (RYA-841: it equated ~10% flux with 0.10 dex, implicitly assuming
dA/dδ = 1.0), but it is not far off. What changes is that the number now has a measurement
behind it and a range attached, instead of being a round number.

**The lever used is 2.54, not 2.42** — the median over the 28 lines that *have* a
derivative. Averaging a slope over lines that have no slope is not a measurement (§3).

---

## 🔴 The obvious route is dead, and it would have answered zero

The Kitt Peak atlas ships as 250 `lm####` segments, and **every seam overlaps its neighbour
by ~1.0 Å** — 20 overlapping seams inside 3000–3780 Å. Two independent scans of the same
wavelengths would measure the atlas's own normalisation repeatability directly, with no
model involved. That is exactly σ_δ.

**They are not independent.** The overlapping pixels sit on an *identical* wavelength grid
with a flux difference of **exactly 0.000e+00**. The atlas is one continuous spectrum split
into files with 1 Å of padding.

Had that route been used without checking, it would have returned **σ_δ ≈ 0** and hence a
continuum term of ≈ 0 — a spuriously perfect answer, in the direction that flatters the
product. The check cost one command.

---

## 1. What was measured instead

A second, independent reduction of the same band was already on Sirius: **Wallace, Hinkle,
Livingston & Davis 2011 (NSO)**, a telluric-removed KPNO solar flux atlas fetched under
RYA-485 — its own `PROVENANCE.json` calls it the *"tie-breaker for the continuum lever:
distinguishes a KPNO-family systematic from a Kurucz-reduction-specific one"*.

| region | coverage | quality |
|---|---|---|
| 6 | 2899–3525 Å | clean — 0 overflow markers |
| 5 | 2995–4106 Å | 2,227 overflow markers **all inside 3000–3200 Å**; clean above |

Region 5's defects are **measured, not assumed**, to be confined to 3000–3200 Å (2,227 NaN
plus 7,460 out-of-range pixels there, 6.4% of that sub-band; every sub-band above 3200 Å has
0 NaN, 0 out-of-range, max flux ≤ 1.02). So the composite uses region 6 where it exists and
region 5 only above 3524 Å — covering **all 40 product lines**, none left out.

### Three controls, each able to abort the run

1. **Wavelength convention.** Wallace ships **wavenumbers**, so `1e8/ν` is a **vacuum**
   wavelength while the Kitt Peak files are **air**. At 3300 Å that is a **0.94 Å** error —
   wider than a line — and it would have read as a normalisation difference. After the
   Edlén conversion the cross-correlation residual is **+0.010 Å at r = 0.9981**;
   uncorrected it is −0.940 Å at r = 0.928.
2. **The two atlases must agree on the lines**, or a ratio between them is meaningless.
   Median per-bin spectral correlation **0.991**.
3. **Wallace against itself** — see §2, and it is the control that changes the answer.

⚠️ The Wallace files contain `*******` where a Fortran field overflowed. Those become NaN
and are **counted**; `genfromtxt(invalid_raise=False)` would have dropped exactly the rows
where something went wrong, silently.

---

## 2. The control that halves the answer

Regions 5 and 6 overlap, so Wallace can be compared **with itself**. Over 3200–3524 Å where
both are clean, 32 bins:

| | rms deviation from unity |
|---|---|
| Kurucz vs Wallace, at the 40 product lines | **4.94%** |
| **Wallace region 5 vs region 6 (internal)** | **4.07%** |
| net, in quadrature | **2.79%** |

**Most of the apparent reduction difference is internal to Wallace.** Without this control
the honest term would have been quoted as 0.125 dex when the part attributable to a genuine
Kurucz-vs-Wallace normalisation difference is 0.071 dex. Both are reported; the raw figure
is an upper bound, and the quadrature subtraction assumes the two are independent, which is
stated rather than hidden.

### The disagreement is structured, not noise

MAD 2.07% against std 4.69% band-wide — heavy tails, and they are coherent blocks, not
scatter:

| region | Wallace / Kurucz |
|---|---|
| 3125–3205 Å | **≈ 0.87–0.91** (Wallace ~10% lower) |
| 3445–3485 Å | **≈ 1.09** (Wallace ~9% higher) |

Spectral correlation holds at 0.987–0.991 *inside those blocks* — the lines agree, the
continuum **level** does not. It is also worst blueward (MAD 3.5% over 3000–3200 Å against
1.3% over 3200–3400 Å), which is what heavier blanketing should do to continuum placement,
and matches RYA-759's blueward degradation.

Because it is structured, the band-average and the value **at the measured lines** are
different questions. The latter is what enters the budget, and it is what is quoted.

---

## 3. The 12 lines with no derivative are not a continuum problem

RYA-841 found 12 of 40 lines whose continuum response is not linear (r² ≤ 0.98), several
with the physically impossible sign. They differ from the other 28 in a way the continuum
cannot explain:

| property | linear (28) | non-linear (12) | permutation p |
|---|---|---|---|
| **median A(Fe I)** | **7.413** | **7.880** | **0.0029** |
| median reduced χ² | 81.6 | 136.6 | 0.059 |
| median theoretical depth | 0.898 | 0.899 | 0.74 |
| median excitation potential | 2.855 | 2.928 | 0.83 |
| median core dominance | 0.774 | 0.529 | 0.26 |

They sit **+0.467 dex high at p = 0.003, at identical line strength.** A continuum
misplacement cannot do that selectively to 12 of 40 lines while leaving their depths alone.

⚠️ **The mechanism is not established, and an earlier draft of this claimed it was.** The
natural story is unmodelled opacity in the window — it would explain the high abundance, the
worse fits and the absent derivative at once. But the blend metrics do **not** confirm it:
core dominance is p = 0.26, and depth, EP and window occupancy are flat. What is established
is that they are a distinct population the continuum does not explain; naming the cause is
per-line work this ticket does not do.

**They are not dropped.** Dropping them would move the product **7.488 → 7.413** and the
scatter **0.413 → 0.347** — precisely the RYA-777 shape, a tighter bar bought by removing
lines. RYA-844 requires a stated *per-line* physical reason, and *"its continuum response is
not a derivative"* is evidence toward that case, not a substitute for it. They are carried
with their evidence (RYA-711) as candidate problem-children.

What does change is the lever: **2.54** over the 28 lines that have a derivative, against
2.42 over all 40.

---

## 4. The revised near-UV cells

On top of RYA-845's double-count removal:

| cell | published (assumed 0.100) | net (0.071) | raw (0.125) |
|---|---|---|---|
| Fe I near-UV **1D-LTE** | 0.1972 | **0.1842** | 0.2113 |
| Fe I near-UV **1D-LTE-LABGF** | 0.1081 | **0.0820** | 0.1320 |

The abundances, line counts, scatters and `stat` values are untouched — only the
systematic moves, and only because one of its terms is now measured.

**Not applied.** Which figure to adopt — net or raw — is a judgement about whether Wallace's
internal scatter should be credited against the Kitt Peak normalisation, and that is Ryan's
call, not a measurement.

---

## Owed

1. **σ_δ is a lower bound.** It measures reduction-to-reduction disagreement, not distance
   from an unobservable true continuum. Nothing here can bound the part both reductions get
   wrong together.
2. **The structured blocks** (3125–3205 Å, 3445–3485 Å) are unexplained. If one reduction is
   simply wrong there, σ_δ outside those blocks is smaller than quoted.
3. **The 12 responders want per-line adjudication** under RYA-844.
