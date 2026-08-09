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
