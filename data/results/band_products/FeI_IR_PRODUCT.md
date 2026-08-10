# Fe I in the infrared — the first product outside the optical — RYA-713

Ryan: *"what about the VALD references in the IR? Not everything needs to be NIST graded.
Especially if it is a good line with good science."*

That was the unlock. NIST accuracy grades cover a small fraction of what is actually
well-determined, and treating "ungraded" as "unknown" understated the data badly.

## The headline

**A(Fe I) = 7.533 ± 0.029 (stat) ± 0.043 (syst), n = 29**, red-optical 6910–9199 Å,
Kitt Peak, 1D-LTE, on laboratory oscillator strengths only.

Sources: Blackwell/Wells/Lynas-Gray (Oxford furnace) 11 · Ruffoni/Heiter/Lawler 2014 FTS 6
· Ruffoni 2014 MNRAS 5 · Bard & Kock 4 · Bard/Kock/Kock 2 · combinations 1.

## Why NIST alone was the wrong question

Pulling NIST ASD across 6900–9300 Å returns **504 Fe I rows, of which only 74 carry any
accuracy grade**. Of our 103 measured lines, 81 match NIST and **28 carry a grade** — and
those grades are mostly poor: 1×A, 4×B, 6×C+, 7×C, 5×D+, 5×D.

So "get NIST grades for the IR lines" has no good answer: **NIST largely does not have
graded oscillator strengths for solar Fe I in the infrared.** Only 5 of 103 reach B or
better.

But VALD carries the *source*, and a laboratory measurement is excellent whether or not
NIST graded it:

| provenance | n | gf uncertainty |
|---|---|---|
| Blackwell/Wells/Lynas-Gray — Oxford Fe I furnace | 13 | 0.020 |
| Ruffoni/Heiter/Lawler 2014 — FTS | 6 | 0.030 |
| Ruffoni et al. 2014 MNRAS — FTS | 5 | 0.030 |
| Bard & Kock 1994 | 5 | 0.040 |
| Bard, Kock & Kock 1991 | 2 | 0.040 |
| Fuhr, Martin & Wiese — NIST critical compilation | 12 | 0.080 |
| May, Richter & Wichelmann 1974 | 3 | 0.080 |
| GES-adopted Blackwell 1982 | 2 | 0.080 |
| **Kurucz 2007 — semi-empirical** | **54** | **0.200** |

**46 of 103 have laboratory or critical-compilation provenance.**

## The internal validation

Splitting by gf provenance and re-aggregating — **the measurement is untouched, only which
lines enter**:

| tier | n | A(Fe I) | line-to-line scatter | gf dex |
|---|---|---|---|---|
| all lines, blanket assumption | 101 | 7.639 | 0.357 | 0.170 |
| laboratory ≤ 0.08 | 46 | 7.621 | 0.231 | 0.054 |
| **best laboratory ≤ 0.04** | **29** | **7.533** | **0.156** | **0.028** |
| Kurucz semi-empirical only | 54 | 7.661 | 0.392 | 0.200 |

Two things fall out, and neither was arranged:

1. **The observed scatter tracks the gf quality** — 0.156 / 0.231 / 0.392 as provenance
   degrades. Independent confirmation that oscillator strengths, not the measurement, were
   the dominant scatter source.
2. **The abundance converges toward the optical as gf improves** — 7.661 → 7.639 → 7.621 →
   **7.533**, against the banked optical **7.466**.

## Cross-band comparison — reported, never adjudicated (RYA-712)

| | A(Fe I) | n |
|---|---|---|
| IR, best-laboratory gf | 7.533 | 29 |
| VIS, banked | 7.466 | 444 |
| **difference** | **+0.067 dex** | **+1.2 σ combined** |

**Consistent within the combined uncertainty.**

The blanket-gf version sat **+0.173** above the optical. Restricting to laboratory
oscillator strengths moved the IR **0.106 dex** toward the optical **without touching the
measurement** — exactly what you expect if the offset was a gf artifact rather than physics.

That is a comparison, not a validation: the IR was never required to reproduce the optical,
and had it stayed discrepant that would have been a finding rather than a failure.

## What this supersedes

The repeated conclusion that "the IR bottleneck is atomic data, and no defensible abundance
comes out until the lines are graded" was **half right**. The atomic data *is* the limiting
term — but the fix was not to obtain NIST grades that largely do not exist. It was to use
the provenance already carried in the line list. The 0.17 dex blanket was a floor imposed by
the worst source in the set, applied to lines that never deserved it.

## Owed

* The 54 Kurucz-gf lines are not wrong, only imprecise — they are retained with their 0.20
  dex term and reported separately, never merged with the laboratory set.
* Engine A on the laboratory subset (the 44-line Engine-A product spans mixed provenance).
* The same provenance split for the near-UV, once anything can measure it.
