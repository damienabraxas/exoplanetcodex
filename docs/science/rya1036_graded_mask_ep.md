# RYA-1036 — EP through the product path, then the widening

**A LAB line was published UNGRADED because `_graded_mask` keyed on wavelength alone.**
Threading EP through and widening to 30 mA *with EP required* recovers **3 Fe I lines across
every Fe artifact, loses none**, and refuses three λ-coincidences that a wavelength-only
widening would have swallowed.

## The defect

`_graded_mask(waves)` matched a measured line to the `gf_tier == LAB` rows of `canonical_gf`
by nearest wavelength within **5 mA**, with no EP. The two files write the same transition to
different precision, so a genuine lab line whose canonical row sits further away came back
UNGRADED — **tier disagreeing with provenance**.

This is the product-path twin of RYA-871 (which threaded EP into the EW-route artifact) and
the mirror of RYA-1034: there an ungraded line wore a GRADED tier by λ coincidence; here a
lab line wore an UNGRADED one. Same defect, opposite sign.

## The sequence, and why it is not negotiable

**EP first, widening second.** Measured on this repo's `canonical_gf`: **378 Fe I clusters sit
within 5 mA of each other, and 340 of them disagree on gf while being separable by EP.** The
ticket's own example reproduces exactly:

| λ (Å) | EP (eV) | log gf | source |
|---|---|---|---|
| 6065.4820 | 2.609 | −1.530 | NIST ASD grade A |
| 6065.4850 | 4.956 | −3.471 | K07 |

**1.94 dex apart, 3 mÅ apart.** Widening the window on wavelength alone would swallow those —
the coin flip RYA-1033 killed, coming back through a wider door. (This is also the same line
RYA-853's EP guard refused, from the other direction.)

## What recovered

Routed through `pipeline.line_match` — the one canonical matcher (RYA-1033/1037) — with
`require_ep=True`, `tol_A=0.030`, `ep_tol_eV=0.05`.

| measured λ | EP | canonical λ | Δλ | ΔEP | source |
|---|---|---|---|---|---|
| **6705.1169** | 4.6070 | 6705.1010 | 15.9 mÅ | **0.0000** | **PRIMARY LAB Ruffoni2014** |
| 8598.8187 | 4.3865 | 8598.8290 | 10.3 mÅ | 0.0005 | PRIMARY LAB DenHartog2014 |
| 8945.1747 | 5.0331 | 8945.1890 | 14.3 mÅ | 0.0001 | PRIMARY LAB Ruffoni2014 |

Every ΔEP ≤ 0.0005 eV — **same-transition evidence**, not a looser window. Per artifact:

| artifact | graded before → after |
|---|---|
| HARPS VIS (molecfit) | 6 → 7 |
| KPNO VIS (molecfit / kurucz2005) | 13 → 14 |
| KPNO red-optical | 20 → 22 |
| IAG VIS | 5 → 6 |

**Zero lines lost anywhere.**

## 🔴 What it must NOT recover — and did, until a matcher defect was found

Three lines sit inside 30 mÅ of a LAB row and would be "recovered" by a wavelength-only
widening, while their EP says a different transition:

| λ | ΔEP to the row it matched |
|---|---|
| 6858.1396 | **1.006 eV** |
| 8713.1976 | **2.039 eV** |
| 8876.0059 | **0.436 eV** |

They came through on the first run. The cause was in the canonical matcher, not here:
`match()` applied the EP filter **only when there was more than one candidate**, so a **lone**
row inside the wavelength window was accepted with **no EP check at all** — necessary but not
sufficient. That is the Fe I 6065.490 shape exactly: one candidate, wrong level, no ambiguity
flag to warn anyone.

Fixed in `pipeline/line_match.py` under `require_ep=True` only, so no existing caller's match
set moves; an all-EP-rejected window is now reported as `unresolved` rather than crashing.
**RYA-1037's strict mode was necessary-but-not-sufficient until this ticket exercised it.**

## ⚠️ The published-value delta is NOT reported here, and cannot be

The ticket asks for the per-line published-value delta. **The `band_ew` artifact stores EW,
not abundance** — `abundance` is finite on **0 of 1826 rows**; A is derived downstream in
`derive_band_products`. So the value delta needs a **product re-run**, not a re-read, and that
needs the Mac synthesis stack.

What is established here is the **pool composition** change: +1 line on VIS, +2 on
red-optical, −0 anywhere, each with ΔEP ≤ 0.0005 eV. The value moves as a consequence and
should be reported from the re-run — the justification (these are Ruffoni/Den Hartog lab
lines matching to 0.0000 eV) is fixed **before** anyone looks at where A(Fe) lands (RYA-161).

## On the ticket's counts

The ticket says *"20 of 22 unresolved HARPS VIS lines"*. Measured today: **21 of 23** resolve
uniquely at 30 mÅ with EP required. The claim holds; the count moved by one because
`canonical_gf` has grown since (RYA-945/1053/1075). Note those 23 are unresolved against **all
of `canonical_gf`**, a different and larger set than the graded-mask misses this ticket fixes
— the graded-tier recovery is the 3 lines above.
