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

---

# Continuum normalisation — yes, completely different, and it was biasing the EWs

Ryan: *"I presume we are doing a completely different normalization of the continuum for
KITT data correct?"* — Yes, and checking it found a real error in my own harness.

## Two instruments, two normalisation histories

| instrument | arrives as | who set the continuum |
|---|---|---|
| **Kitt Peak FTS** | column 1 is **residual flux** | Kurucz, upstream — already normalised |
| **HARPS** | un-normalised | our pipeline |
| IAG FTS | residual flux | upstream |

So any continuum *we* fit on Kitt Peak is a **second** normalisation stacked on Kurucz's.

## How good is the atlas continuum? Measured, not assumed

95th percentile of the flux over ±8 Å:

| band | 95th pct | median | reading |
|---|---|---|---|
| 3200 Å | 0.894 | 0.541 | **depressed** — line blanketing |
| 3900 Å | 0.885 | 0.610 | **depressed** — line blanketing |
| 4600–9000 Å | 0.986–0.997 | — | atlas continuum is **excellent** |
| 11500 Å | 0.894 | 0.661 | **telluric**, not blanketing — inside the 11120–11560 H₂O band |

*(I first labelled the 11500 Å depression as blanketing. It is the H₂O band already in
the telluric list.)*

## The error this found

My first pass re-normalised every window to a local 95th-percentile side-band continuum.
Across 6910–9199 Å that changed EWs by a **median of −11.7 %, worst −71.4 %** — roughly
0.05 dex, the entire size of the NLTE effect being measured.

The cause is not the continuum *level* (median 0.9947) but **which side-bands were used**.
The two worst-shifted lines had side-bands at **0.902** and **0.936** — the side-bands
were themselves absorbed, so dividing by them removed real line flux:

| line | side-band p95 | EW shift from re-normalising |
|---|---|---|
| 7363.922 | 0.9022 | −71.4 % |
| 9013.977 | 0.9364 | −60.2 % |
| 7093.080 | 0.9776 | −30.7 % |
| 7751.108 | 0.9919 | −3.2 % |
| 8075.149 | 0.9953 | −6.1 % |

Where the atlas continuum is already good, a second local normalisation **is not
correcting an error — it is introducing one.**

## Fixed

`equivalent_width()` now takes `pre_normalised`, set per instrument from a declared map
rather than guessed. On pre-normalised data the atlas continuum is trusted, and the
side-bands are still *measured*: if they sit below 0.99 the line carries a **CONCERN**
noting the window is crowded and the EW may include unresolved neighbours. **13 of 29**
Fe I IR windows carry that concern.

Correlation of intrinsic strength against log EW, as the method improved:

| method | r |
|---|---|
| raw window EW, local re-norm | +0.399 |
| + feature verification | +0.562 |
| + atlas continuum trusted | **+0.587** |

## Why this matters beyond Fe

**This is a live candidate for the Mg 5711 disagreement.** Mg I 5528 agrees across
instruments to 2.9 % while 5711 disagrees by 24 %, and continuum methodology was the
untested suspicion. Two instruments with two normalisation histories, compared line by
line, is exactly the mechanism — and at 5711 Å the Kitt Peak atlas continuum measures
0.9862, i.e. good, so a second re-normalisation there would have biased one arm only.
Worth testing directly now that the policy is explicit.

**The gf blocker above is unchanged.** Better EWs do not fix ungraded oscillator
strengths, and no Fe IR product is defensible until the 17 verified lines are graded.
