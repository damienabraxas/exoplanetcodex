# The Fe 1D→⟨3D⟩ climb is dimensionality, not model family — RYA-1032

**Status: answered. No solver was built; the ticket's founding fact was wrong.**
Date 2026-08-29. Generator `scripts/rya1032_fe_dimensionality_ladder.py`, guard
`tests/test_fe_dimensionality_ladder_rya1032.py`. Every number below is derived at
run time from the live feed `data/products/solar/Fe.json` (v1.86) joined to
`data/catalog/model_registry.csv` — none is typed into this document by hand.

---

## 1. The founding fact was wrong, so the ticket collapsed

RYA-1032 was written to **solve the Fe departures ourselves**, on the stated fact:

> *"There is no Fe ⟨3D⟩ departure deck, anywhere. Verified 2026-08-24 across all three
> Sirius drives… So this is not an acquisition. There is nothing to fetch and nobody
> to email."*

That was a **`find` across three drives — a disk scan reported as a fact about the
source.** `NLTEgrid4TS_Fe_STAGGERmean3D_May-21-2021.bin` (92,945,908 B) had been on
the same MPG Keeper share we already fetch every Gerber deck from **since 2021**, and
TSFitPy's own downloader lists it as `[Fe] 3d_bin_link`. We had pulled Al, Cr, Eu and
Y from that share and never opened the Fe folder.

* **RYA-1035** found it (Step 0 = HAVE) and staged it.
* **RYA-710 / PR #385** wired it — merged 2026-08-25.
* `model_registry.csv` row 6 (`synth-mean3D-NLTE-gerber-stagger`) has been `live` since.

**Consequences.** Deliverables 1–3 — build a solver, reduce the 1D deck as a control,
solve the ⟨3D⟩ departures — are **moot**. The brief's *"do NOT reuse `model=gerber`,
that names a deck we did not get from Gerber"* is **void**: the deck *is* Gerber's, so
`model=gerber` is correct. What remained was deliverables 4–5, and those were already
produced by the consume route.

> A disk scan never measures the source. The same error produced the availability
> matrix's 13 `T2_BUILD_OWED` cells that RYA-1035 moved to `T2_FETCH_OWED`.

---

## 2. The ladder

Solar Fe I, VIS, GRADED pool, holding `solar_kpno_molecfit_corrected`. All synthesis
rows share `gf=kurucz`, so they are directly comparable.

| A | n | model |
|---|---|---|
| 7.447 | 67 | Synth · 1D-LTE (atlas9) |
| 7.451 | 67 | Synth · 1D-LTE · Gerber (marcs-ges) |
| 7.454 | 61 | Synth · 1D-NLTE · Bergemann (atlas9) |
| 7.497 | 67 | Synth · 1D-NLTE · Gerber (marcs-ges) |
| **7.552** | 67 | **Synth · ⟨3D⟩-LTE · Gerber · stagger** |
| **7.552** | 67 | **Synth · ⟨3D⟩-NLTE · Gerber · stagger** |
| 7.5116 | 50 | EW · 3D-NLTE · Amarsi |

---

## 3. The answer

Each difference below **holds a named axis fixed and varies exactly one**, and the
generator raises rather than reporting a number whose contract does not hold.

| Δ (dex) | comparison |
|---|---|
| **+0.101** | **atmosphere step** — 1D-LTE → ⟨3D⟩-LTE, model family fixed at Gerber, LTE fixed |
| +0.046 | 1D NLTE step — Gerber, marcs-ges, dimensionality and atmosphere fixed |
| +0.043 | model-family spread — Gerber − Bergemann at 1D-NLTE ⚠️ pools differ (67 vs 61) |
| +0.004 | **de-confounding control** — atlas9 → marcs-ges in LTE |

**The de-confounding control is what makes the headline defensible.** The family
spread compares `bergemann@atlas9` against `gerber@marcs-ges`, so it carries an
atmosphere change as well as a family change. Measuring that atmosphere change on its
own — in LTE, where no NLTE physics can contribute — gives **+0.004**. So the family
spread is genuinely family: **+0.039 de-confounded**.

> **Verdict: the 1D→⟨3D⟩ Fe climb is DIMENSIONALITY, not model family.**
> The atmosphere step is **2.35× the raw family spread, 2.59× de-confounded.**

This is the question the ticket was written to answer: *"A ⟨3D⟩ point would say
whether that climb is dimensionality or model family."* It is dimensionality, by a
factor of ~2.6, and the two effects are independent — the atmosphere step is measured
with the family held constant.

---

## 4. 🔴 The ladder is NON-MONOTONIC — ⟨3D⟩ overshoots full 3D

The brief expected ⟨3D⟩ to be *"the missing rung between 1D-NLTE and full 3D-NLTE"*,
landing between 7.45–7.50 and Amarsi's 7.511. **It does not.**

```
1D-NLTE 7.497   <   full-3D Amarsi 7.5116   <   <3D> 7.552
                                    overshoot = +0.040
```

⟨3D⟩ lands **+0.040 above** the full-3D reference it was supposed to approach. Averaging
the atmosphere is not the same as averaging the emergent spectrum, and here the mean-3D
approximation **over-corrects past** the full-3D answer rather than converging on it.

This is consistent with — and sharpens — RYA-1099's independently-measured **+0.105 ⟨3D⟩
inflation**, whose structural cause is still open (leads: the τ500-native STAGGER model
shipping neither τ_Ross nor P_g; the deck's A(Fe)=7.50 baseline; the missing
line-strength axis; vmic is **refuted** — ξ=0 made it 0.137 dex worse).

**Consequence for the programme.** The full cube is not a nice-to-have: it is the
independent check that should land back down at ~7.5116. That is the Bride
(Lightweaver, RYA-1119 / M4).

---

## 5. What this deliberately does not do

The ⟨3D⟩ NLTE effect is **not** computed here. Both published cells are **medians**,
and on this holding they collided on exactly 7.552 — so differencing them returns
**0.000 for a real per-line shift of +0.032**. That is RYA-1099 Finding 3, and the
correct statistic already exists in `pipeline.paired_differential` (RYA-1083), emitted
as `{stem}_{treatment}_nlte_effect.json`. Writing a second one here would be a second
implementation of one thing.

The generator reports the collision and refuses the subtraction; a test asserts that
refusal, so the deriver cannot quietly start publishing that zero as physics.

⚠️ **Three statistics of one dataset have been quoted interchangeably across tickets:**

| statistic | value |
|---|---|
| difference of medians (what the cells publish) | 0.000 |
| **median of per-line differences — the RYA-542 quantity** | **+0.032** |
| mean of per-line differences | +0.0267 |

**Publish +0.032.** They are the same data, not a disagreement.

⚠️ **Open, carried from RYA-1099:** the RYA-1083 fix is wired to this pair
(`treatment_axes.py:432` maps model 6 → model 5) but **has never fired** — zero
`*_nlte_effect.json` exist anywhere. The fix landed `942929a` at 2026-08-27 12:31:53;
the live ⟨3D⟩ products were written at 12:03, 28 minutes earlier, and were never
regenerated. That is re-run debt, and it belongs to RYA-1099.
