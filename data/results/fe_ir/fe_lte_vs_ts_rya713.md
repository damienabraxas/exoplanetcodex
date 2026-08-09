# Fe I — 1D-LTE vs TS-NLTE, optical and IR — RYA-712/713

**The first Fe NLTE measurement outside the optical.** Reported, not gated: Fe's
Engine A (Bergemann MPIA per-line table) stops at 6843.7 Å, so past that there is no
comparand and nothing to validate against. Per RYA-712, LTE and TS are separate
products and the delta between them is the reportable quantity.

## Result

| regime | lines (Å) | EW (mÅ) | EP (eV) | median δ = A☉ − A* |
|---|---|---|---|---|
| optical, strong | 5194.9 / 5202.3 / 5217.4 | 155 / 78 / 130 | — | **+0.0734** |
| optical, weak | 6745.1 / 5491.8 / 6105.1 | 8.8 / 9.5 / 9.8 *(measured)* | 4.19–4.58 | **+0.0543** |
| **IR** | 7212.4 / 7751.1 / 8526.7 / 9024.4 | 43 / 48 / 71 / 45 *(synth)* | 4.91–4.99 | **+0.0506** |

`NLTE departure engaged: True` in the IR — the atom and grid work out there.

## Reading

**The two unsaturated sets agree to 0.004 dex across 3533 Å of spectrum.** They are
not the same levels sampled twice: the optical set sits at 4.19–4.58 eV and the IR set
at 4.91–4.99 eV, and departure coefficients are per-LEVEL. Two distinct level
populations, two bands, one answer — Fe I solar NLTE ≈ **+0.05 dex**.

**The optical-strong +0.0734 was saturation-inflated**, as suspected. All three of its
lines sit above the 100 mÅ knee and its two strongest gave the worst ratios. Both
unsaturated sets land ~0.02 dex below it.

**The MPIA anchor of +0.010 now disagrees with three independent measurements**, not
one. That shifts suspicion toward the anchor rather than toward Gerber — and it is the
same shape as Ti, where Gerber gives +0.2206 against an MPIA anchor of +0.107
(RYA-548, open). Both anchor on Bergemann-MPIA; **Gerber runs high in both**. Two
elements is a pattern. The families differ in code and atmosphere (DETAIL/SIU on
MAFAGS-OS vs Turbospectrum on MARCS), so this is a cross-family systematic, not a bug
in either.

## What this does NOT establish

* **Not a validation.** LTE-vs-TS measures the *size* of the correction, not its
  *correctness*. If TS is wrong in the IR, this shows a delta and cannot say which
  side is right. A true IR cross-check needs a second engine, which Fe does not have.
* **The IR set spans a narrow EP range** (4.91–4.99 eV) — one level neighbourhood.
  Wavelength-independence across the *full* Fe I term structure is not shown.
* **Four lines.** 454 unmeasured Fe I lines exist in the runnable window; these span
  it but do not populate it.

## Hard boundary — 9199.9 Å

The GES level-identified linelist ends there. Past it a line carries no level
identification, so NLTE **silently falls back to LTE** in the IDL-SME path. 454 Fe I
lines are runnable in 6910–9199 Å; the remainder of Fe's ~4000 IR lines need a
level-identified linelist before *either* engine can touch them. Not a grid problem.

One line to watch: **9024.369 gives +0.0726**, well above the other three
(+0.0437…+0.0513). It is the reddest and closest to the linelist edge. Unexplained.

---

# UPDATE — Engine A **does** reach the IR. The two-engine comparison is now real.

Ryan: *"double check, there is Engine A for Fe? in IR???"* — **Yes. I had it wrong, and
the cause was the same bug one level deeper.**

## Why we thought Fe's Engine A stopped at 6843.7 Å

`scripts/build_fe_nlte_grid_rya319.py::matched_waves()` scrapes MPIA only for lines that
already exist in **our EW pool**:

```python
our = pd.read_csv(PATHS['solar_ew'])
our = our[(our.element == 'Fe') & (our.ion == ion)]['wavelength_air_A']
```

The EW pool is HARPS-only, 3924–6905 Å. **So 6843.7 Å is our own pool's ceiling reflected
back through a scraper filter into an artifact that was then read as the model's limit.**
Fifth occurrence of one bug.

## What MPIA actually serves (live query, 2026-08-09)

| ion | lines served | span (Å) | > 6910 Å | < 3791.7 Å |
|---|---|---|---|---|
| Fe I | 530 | 3121.8 – **11973.0** | 28 | 155 |
| Fe II | 1795 | 2900.8 – 9113.0 | 290 | 710 |

**We requested 332 of 2325 — 14%.** 1183 Fe lines sit outside our HARPS window with
Engine A available and were never asked for. **780 of those are near-UV lines reachable
by the Kitt Peak atlas** (155 Fe I + 625 Fe II at ≥2960 Å).

Verified live at the solar node — real values, not placeholders:
Fe I 3121.8 `+0.0100` · 3565.4 `+0.0080` · 7491.7 `+0.0120` · 8046.1 `−0.0010` ·
11200.8 `+0.0120` · 11973.1 `+0.0070`.

**Fe II is NOT clean outside the optical**: 6913.69, 8446.36 and 9112.95 return `nan`,
and 7711.72 returns a placeholder zero. Fe II needs the RYA-417 placeholder audit before
any IR/near-UV use. Fe I is sound.

## An inversion worth recording

**Fe I has 12 MPIA lines at 11200–11973 Å — past where Engine B can go.** The GES
level-identified linelist ends at 9199.9 Å. So in the J-band the situation is the
*reverse* of what was assumed: **Engine A reaches where Engine B cannot.**

## The line-matched result

Four transitions served by both engines (in MPIA's `felines[]` *and* under GES 9199.9 Å),
so this compares the same lines rather than two ensembles:

| line (Å) | EW (mÅ) | Engine A (MPIA) | Engine B (Gerber) | B − A |
|---|---|---|---|---|
| 7491.650 | 71.3 | +0.0120 | +0.0478 | +0.0358 |
| 8046.070 | **165.4** | −0.0010 | +0.0168 | +0.0178 |
| 8293.530 | 58.4 | +0.0160 | +0.0917 | +0.0757 |
| 8576.500 | 3.9 | +0.0120 | +0.0736 | +0.0616 |
| **median, all 4** | | **+0.0120** | **+0.0607** | **+0.0487** → PASS |
| **median, unsaturated** | | **+0.0120** | **+0.0736** | **+0.0616** → **CHECK** |

**The PASS is an artifact and must not be recorded as one.** It holds only because
8046.07 sits at 165 mÅ, far above the 100 mÅ knee, where inverting EW through a flat
LTE curve-of-growth is ill-conditioned. Excluding it — which the saturation rule
requires — gives **+0.0616 and a CHECK.**

## What this settles

**The Gerber-minus-MPIA offset is wavelength-independent.** Optical weak lines gave
+0.044; IR line-matched unsaturated gives +0.0616; the IR ensemble comparison gave
+0.0506 − 0.0120 = +0.0386. Across **5491–8576 Å** the systematic sits at
**≈ +0.04 to +0.06 dex** and shows no trend with wavelength.

That makes it a **cross-family systematic, not an IR artifact** — DETAIL/SIU on MAFAGS-OS
versus Turbospectrum on MARCS. Combined with Ti (+0.2206 vs +0.107, RYA-548 open), the
pattern is now two elements, two bands, one direction: **Gerber runs high against
Bergemann-MPIA everywhere we have looked.** Resolving it is a physics question about the
two families, not something a wider tolerance should absorb.

## Fixed en route

`ges_lines()` matched at 0.02 Å and rejected MPIA 8046.070 against GES 8046.046
(d=0.024) — two catalogues quoting one transition. Widened to 0.05 Å **with an ambiguity
guard**: the runner-up must be ≥2× further, else the run aborts rather than guessing.
Widening alone would have been the RYA-704 error (a coarse key silently merging lines).
