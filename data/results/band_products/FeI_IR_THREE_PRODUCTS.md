# Fe I in the infrared — the three products — RYA-712/713

Ryan: *"lets try and wrap up our work on Fe, so we can start the other elements as well."*

**Fe I · red-optical 6910–9199 Å · Kitt Peak · laboratory-gf lines only.**
Three separate products. Never combined.

| treatment | A(Fe I) | n | stat | syst | scatter | vs VIS 7.466 |
|---|---|---|---|---|---|---|
| **1D-LTE** | **7.533** | 29 | 0.0289 | 0.0432 | 0.156 | +0.067 (**+1.2σ**) |
| **ENGINE-A** | **7.571** | 14 | 0.0548 | 0.0422 | 0.205 | +0.105 (**+1.5σ**) |
| **ENGINE-B** | **7.602** | 14 | 0.0633 | 0.0430 | 0.237 | +0.136 (**+1.7σ**) |

**All three are consistent with the optical**, and none was required to be — a discrepancy
would have been a finding, not a failure (RYA-712).

## What each is

* **1D-LTE** — profile-fit EWs inverted through the project's Turbospectrum curve of growth.
* **ENGINE-A** — plus Bergemann/MPIA per-line departure corrections, live query, median
  **+0.0145**. MPIA serves 14 of the 29 laboratory lines.
* **ENGINE-B** — plus Gerber TS-native NLTE, median **+0.0663**, on the 14 laboratory lines
  under the GES level-identification limit of 9199.9 Å.

The **cross-engine difference is +0.0423 dex on 6 shared lines**. That is a **diagnostic**,
reported alongside; it is not a product and the two are never averaged.

## The trade the separation makes visible

**Engine A and B have WORSE uncertainties than 1D-LTE**, not better — 0.0548 and 0.0633
statistical against 0.0289 — because applying a departure correction halves the sample
from 29 lines to 14. The NLTE corrections are real (+0.015, +0.066) but modest against
what the reduced line count costs.

That is a genuine trade, and it is visible *only* because the products are kept apart. A
single merged number would have hidden it entirely.

## Why laboratory-gf only

Restricting to laboratory oscillator strengths (Blackwell Oxford, Ruffoni FTS, Bard & Kock;
≤0.04 dex) rather than a blanket 0.17 dex assumption moved the 1D-LTE value **0.106 dex**
toward the optical **without touching a single measurement**, and tightened the line-to-line
scatter from 0.357 to 0.156. Both are what a gf artifact looks like when it is removed.

The 54 Kurucz semi-empirical lines are **not wrong, only imprecise**. They are retained
with their 0.20 dex term and reported separately, never merged into these.

## Flagged, not removed

**Fe I 7107.459** returns a Gerber delta of **−0.1063** while every other line in the set
runs +0.043 to +0.104 — the only negative in 14. It stays in the aggregate with this note
rather than being quietly dropped; a single outlier that changes a median is a finding
about the line or the atom, not a nuisance.

## Fe status

| leg | state |
|---|---|
| VIS control (profile fit) | **PASSED**, −0.013 dex |
| VIS reference | banked 7.466 (EW), 7.520 (synth-v2, reproduced exactly) |
| **IR, three products** | **done** — above |
| near-UV | **unmeasurable by profile fit** — 0 of 1007 lines, both ions, one with laboratory gf |
| NIR > 9199 Å | not started; beyond the GES level-identification limit |

**Owed for Fe:** Engine A/B coverage is 14 lines because MPIA and GES each serve only part
of the laboratory set — extending either would tighten both products. The near-UV needs
synthesis, whose line list is converted but not engine-validated.
