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

---

# Appendix panels + root-cause attribution — every failure proved, and diagnosed

Ryan: *"if that line failed, you would still have to prove it to me in the Appendix, which
is a check in itself. From a plot we can determine, hey do we have the right line? Is it
normalized correctly? Is there a blend? is this a gf ghost?"* and *"In QA, we want to find
root causes. Why did it fail? What is the mechanism? Is it our model? the data? Something
wonky?"*

`scripts/plot_band_appendix.py` draws one panel per line — **passed and quarantined
alike**, quarantined first. A quarantined line without a plot is an assertion; with one it
is a demonstration the reader can overrule. Each panel answers the four questions:

| question | what the panel shows |
|---|---|
| right line? | catalogued position (red dash) vs deepest feature found (purple marker + offset) |
| normalised? | continuum in use (teal rule), side-band regions shaded, their 95th pct (dotted) |
| blend? | integration window shaded; profile drawn wide enough that a double core shows |
| gf ghost? | predicted depth (green bar) against the depth **at the line** |

## Symptom is not cause — the fault domain

A symptom says what we see. It does not say who owns the fix. Every quarantined line now
carries a **fault domain**, the **mechanism**, the **discriminator** that separated it from
the alternatives, and the **fix**. Where the evidence does not separate two candidates the
answer is `UNKNOWN` and both are named — a confidently wrong root cause is worse than an
honest undetermined one.

### Fe I, 6910–9199 Å, Kitt Peak — 12 quarantined

| domain | mechanism | n |
|---|---|---|
| **ATOMIC-DATA** | a real solar line missing from our list dominates the window | **4** |
| **ATOMIC-DATA** | log gf far too strong, or transition assigned to the wrong species | **2** |
| **ATOMIC-DATA** | wavelength error in our line list | **2** |
| METHOD | EW→abundance inversion on the flat part of the curve of growth | 3 |
| METHOD | integration window swallowed a *known* neighbour | 1 |
| MODEL | — | **0** |
| OBSERVATION | — | **0** |

**8 of 12 failures are atomic-data faults. Zero are model faults. Zero are data faults.**

The discriminators that produced that split are cheap and decisive:

* *ghost* → is there a feature of about the right depth within 0.6 Å? **Yes** ⇒ the line is
  real and our **wavelength** is wrong. **No** ⇒ the **gf** is wrong or the species is
  misassigned. Both atomic data, entirely different fixes.
  – 6963.016: *"nothing within 0.6 Å has a depth resembling the predicted 0.368; observed
  at the position is 0.014"* ⇒ gf or species.
* *blend* → is the interloper in **our** catalogue? **Yes** ⇒ we knew, and our window was
  too wide — a METHOD fault we own. **No** ⇒ our line list is missing a real solar line.
  – 7194.901: *"the dominant feature at 7195.036 (depth 0.587) has NO catalogue entry
  within 0.05 Å"* ⇒ our IR list is incomplete.

## What this converges with

Independently, **0 of the 17 surviving lines carry a NIST gf grade**. The root-cause split
says the same thing from the other direction: **the IR bottleneck is atomic data, not
physics and not the spectrum.** Our models reach the IR (both engines run there); Kitt Peak
sees it cleanly (95th pct 0.986–0.997 across 4600–9000 Å). What we lack is a graded,
complete, correctly-positioned IR line list.

That is a far more tractable problem than a missing grid, and it is the same conclusion
RYA-708 reached for Al: **reaching a new wavelength is authoring, not guessing.**

## Two flaws the plots exposed in my own harness

1. Panels printed `obs` from the **deepest feature in the window** while the verdict used
   the depth **at the catalogued position** — so 9179.742 read "obs 0.137" beside a verdict
   of "Sun shows none". Both true, of different points. Panels now report depth at the line,
   with the deepest shown separately when they differ.
2. `SIDEBAND_CLEAN_MIN = 0.99` fired on **13 of 29** windows including side-bands at 0.986,
   which is ordinary IR continuum — the flag was diluting itself. Lowered to 0.97; the
   genuinely absorbed cases measured 0.902 / 0.936 / 0.949. Now **5 of 29**.
