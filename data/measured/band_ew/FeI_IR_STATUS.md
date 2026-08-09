# Fe I in the IR — measurement done, products BLOCKED on gf grading — RYA-713

Ryan: *"We showcase Engine A in one plot/product in the IR, and the same treatment for
B, and LTE."* Correct framing, and the framing I had wrong — I was leading with B−A,
which is a diagnostic, not a product.

**The three products cannot be built yet, and the blocker is not what I expected.**

## What was done

`scripts/measure_band_ew.py` (parameterised: element, ion, band, instrument are all
arguments — no element symbol in the logic, per RYA-701) measured Fe I on the Kitt Peak
atlas across 6910–9199 Å. 30 candidates → 29 measured, 1 with no atlas segment.

## Why raw window-integrated EW is not a measurement here

The first pass looked fine and was not. Correlating intrinsic line strength
(`log gf − θ·EP`) against `log EW` gave **r = +0.399** — for real lines of one ionisation
stage that should be strongly positive. Three failure modes, from the spectrum itself:

| line | what is actually there |
|---|---|
| 7252.635 | deepest feature **0.268 Å away at depth 0.907** — an uncatalogued strong line. The 198.8 mÅ is *that* line |
| 7296.577 | flux at the catalogued position ≈ 0.99. **The line is not there**; 0.16 mÅ was honest |
| 8075.149 | feature dead on (+0.008 Å, depth 0.255) but `log gf = −5.062` says invisible |

Crucially, **no catalogued neighbour flagged any of these** — our IR line inventory is
too sparse to see its own blends. An empty neighbour list means the *catalogue* is empty
there, not the spectrum. Absence of evidence, read as evidence of absence.

`verify_feature()` now checks that the deepest feature is at the catalogued position,
that absorption exists at all, and that observed depth resembles predicted. Failures are
**quarantined, not culled** (RYA-711) — measured, kept, reported, barred from the
aggregate with the reason recorded. Yield: **17 of 29 verified**, and r rises
+0.399 → **+0.562**.

## The actual blocker

**Zero of the 17 verified lines carry a NIST gf grade.** `canonical_gf.csv` holds 31,927
Fe rows and **10** graded ones.

An EW is only as good as the gf it is inverted through. A NIST-B gf is ±0.041 dex; an
ungraded Kurucz semi-empirical value is 0.1–0.3 dex. **The NLTE effect we are trying to
measure is +0.05 dex.** The atomic-data uncertainty is two to six times the signal.

So a 1D-LTE product built on these would carry an honest σ larger than the difference
between all three treatments — and the Engine-A/Engine-B comparison run earlier today,
though real, rests on the same ungraded gf values.

**This is Al's lesson repeating.** Al's four IR lines became usable only because Ryan
pulled NIST grades for them (B/B+). Reaching a new wavelength is authoring, not guessing,
and the authoring step that matters most is the gf.

## What this does NOT invalidate

The engine work stands as a statement about **the models**, because both engines were run
on the *same* lines with the *same* gf — the gf error is common-mode and cancels in the
comparison. Gerber-minus-MPIA ≈ +0.04…+0.06 dex, wavelength-independent from 5491 to
8576 Å, is unaffected. What the gf blocks is the **absolute abundance**, which is exactly
what a product is.

## Owed, in order

1. **NIST gf grades for the 17 verified Fe I IR lines** — the Al treatment. Until then no
   Fe IR product is defensible.
2. Fe II placeholder audit before any Fe II IR/near-UV use (`nan` at 6913.69/8446.36/
   9112.95, placeholder zero at 7711.72).
3. Then, and only then: three products — 1D-LTE, Engine A, Engine B — each with its own
   value, σ, line count and plot. Never merged (RYA-712).
