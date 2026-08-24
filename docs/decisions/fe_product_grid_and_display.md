# The Fe Product: (instrument × band) grid, telluric-corrected display, the KP pre-normalised rule

**Ticket:** RYA-1026. **Ratified by Ryan, 2026-08-23; deliverable 5 widened to the whole
Kitt Peak class by Ryan's 2026-08-24 standing-rule comment.**
**Scope:** structure, display rule and guards. **No value moves under this decision.**

---

## 1. The product unit is (instrument × band), not (element × band)

The reported Fe product is a **per-(instrument × band)** object: `Fe · HARPS · VIS`,
`Fe · Kitt Peak · red-optical` — never `Fe · VIS`. Inside each sits a grid:

| | 1D-LTE | 1D-LTE Synth | 1D-NLTE | 3D-LTE | ⟨3D⟩-NLTE | 3D-NLTE |
|---|---|---|---|---|---|---|
| **gf-graded** | | | | | | |
| **consistent** | | | | | | |
| **bad → appendix only** | | | | | | |

Rows are line tiers, columns are treatments.

**Why the instrument is part of the identity, not metadata.** Two instruments reaching
the same band have different normalisation histories, different telluric states and
different graded-line counts. Collapsing them into one "Fe VIS" number makes a
cross-instrument difference look physical when it may be methodological — the exact
confusion `band_products` already warns about for the KP/HARPS continuum split. The
graded-line count is **a real per-cell n**, never a shared number.

## 2. Gaps are first-class and LOUD

* An instrument that does not reach a band → a **declared coverage gap**.
* Most 3D columns do not exist yet: no 3D-RT synthesis engine in hand (RYA-444), and
  the ⟨3D⟩/3D-NLTE availability survey is RYA-817. Some 3D grids we will not have and
  may have to build ourselves — **that is fine and expected**.
* A missing model or missing reach renders as a **declared, loud gap** (RYA-833 shape),
  **never a blank that reads as zero**.

This is the same failure mode as every silent defect found while building this: a blank
cell, an exit-0 with no output, a zero returned for "not measurable". Wrong information
that does not announce itself is worse than an error.

## 3. Store the axes, derive the display (RYA-906)

Each cell is **route · scale · model-family**. The NLTE and 3D columns each carry a
**model family** (Bergemann/MPIA · Amarsi/Balder · Gerber).

Per RYA-282 these families are **independent measurements**. Where two cover the same
line, **the spread is a budget term, not a pick-one.** The grid stores the family per
cell; the display may collapse or split by family — granularity is Ryan's call.

## 4. Telluric display rule (ratified)

**Displayed science is telluric-CORRECTED input only.** Non-telluric data is not
displayed.

* Most current Fe products were built on `kpno_solar_atlas`, the 1984
  telluric-**uncorrected** atlas. That was a mistake. Those products rebuild on the
  corrected siblings — `solar_kpno_molecfit_corrected` (RYA-940),
  `solar_harps_molecfit_corrected` (RYA-931).
* **EXCEPTION, whitelisted: the KP2005 vs KP1984 pair is retained as the telluric
  CONTROL.** 2005 is telluric-free, 1984 telluric-retaining, so the pair *is* the
  molecfit validation. Removing the uncorrected half would destroy the only thing that
  demonstrates the correction works. A control is not a science product and is labelled
  as one — `as_control=True` must be stated explicitly in code, because a waiver
  obtainable by accident would make the gate decorative.
* Raw/uncorrected data is permitted **only** inside correction R&D for other stars,
  never on a displayed science product.
* Molecfit is trusted-as-working; refinement is deferred to a post-audit and does not
  block.

**Enforced at render** by `pipeline/telluric_display_policy.py`. An uncorrected product
does not *look* wrong — it renders a number with a bar, and the contamination is
indistinguishable from physics unless you already know.

### 🔴 This is stricter than the measurement gate, on purpose

`telluric_policy.gate_holding` (RYA-806) governs **measurement** and deliberately lets
`not-applied` Kitt Peak through, because per-line clean-line selection (RYA-460/786) is a
stated method *defined on* uncorrected data. RYA-1026 governs **display**, and rules that
basis out of the shipped product. Two questions, two gates — **passing the measurement
gate is not permission to ship.**

### 🔴 The clean set is DERIVED, never hand-written

`telluric_applied` is already a first-class **per-holding** column of
`holdings_manifest_registry.csv` (RYA-806), read through `telluric_policy.applied_state`.
A frozenset in the display module would be a second declaration of one fact (the RYA-845
shape) and would go stale silently *in the permissive direction*. The first draft of this
module did exactly that and got two of five entries wrong — it listed `solar_delbouille`,
an id that does not exist (it is `solar_delbouille_liege`), and called it clean when the
registry records it `not-applied`.

Consequence, from the registry as it stands: `solar_harps` and `solar_delbouille_liege`
are **BLOCKED** for display; `solar_kpno_molecfit_corrected`,
`solar_harps_molecfit_corrected` and `solar_kpno_kurucz2005_corrected` are **CLEAN**;
`solar_kpno` is **CONTROL_ONLY**. Unregistered ids return **UNREGISTERED**, kept distinct
from BLOCKED — "we hold this uncorrected" and "we never wrote down what state this is in"
are different problems with different fixes.

`solar_iag` is **CLEAN_WITH_ANOMALY**: the registry says `applied`, but RYA-944 found the
`iag_fts_solar_atlas` manifest routing to the telluric-**retaining** Reiners+2016 file
(46.25% of the O2 A-band below 0.5). No value moves under this ticket, so IAG still
renders — but the doubt travels with it via `anomaly()` instead of being resolved
silently in our favour. It needs its own ticket.

## 5. Do not normalise ANY Kitt Peak atlas — write in stone

**`PRE_NORMALISED = True` for the ENTIRE KP class.** For KP2005 this **REVERSES RYA-929**.

Both KP products ship their own continuum:

| product | what it ships |
|---|---|
| KP 1984 (classic flux atlas) | column 2 is pseudo-residual flux — **unity IS the continuum** (RYA-940 README finding) |
| KP 2005 (Kurucz `irradthu`) | absolute irradiance, its own continuum baked in |

Fitting or applying a normalisation continuum to either **adds a spurious TILT that
follows the saturated bands down and corrupts every EW/synth measurement in that
window.**

**It has bitten twice** — RYA-940 on the 1984 atlas, and the 2005 double-normalise that
forced the VIS re-run. That is why the deliverable was widened from *the 2005 file* to
**the whole KP atlas class**, and to any pre-normalised reference arm: IAG and Delbouille
are the same shape.

**In reverse: the only thing we do to a KP atlas on the way in is telluric correction**
(KP1984 → the RYA-940 corrected sibling). Never a continuum refit.

### Why a guard and not a note

RYA-929 set the flag `False` with a plausible-sounding argument from the file header, and
the harness re-normalised the product for months with nothing objecting. A ratified rule
that lives only in prose gets re-argued by the next reader who reasons from units. So
`pipeline/prenormalised_guard.py` makes the violation **raise** — unrepresentable, not
remembered.

It checks **two independent signals** — the ratified registry and the holding spec's own
`pre_normalised` flag — and treats *disagreement between them* as the bug. Checking one
would simply re-bless whichever was written down. `tests/test_rya1026_product_policy.py`
keeps the two in step by parsing the specs statically (`ast`, never importing
`measure_band_ew`, which `SystemExit`s when the KP atlas is unmounted — that would be an
environment check wearing a policy check's label), and additionally asserts that **every**
HoldingSpec declares `pre_normalised` explicitly rather than inheriting a default.

Three operations are refused, as three separate arguments rather than one boolean,
because they are three different mistakes and the third does not look like a mistake at
all:

1. `fitting_continuum` — fitting one of our own.
2. `applying_continuum` — applying one from elsewhere.
3. `pinning_unity` — re-pinning a **locally** normalised arm to exactly 1.0. Delbouille
   is normalised per-window: no absolute flux, no broad-band continuum shape, and the
   maximum over 2.38 M points is **0.9959** (RYA-944). Tidying that to 1.0 rescales every
   point by a number we chose.

**What is still allowed:** reading the product, measuring its shipped continuum for a
report, or comparing ours against it (RYA-911 `reference_continuum`). Those are readings.
The line is whether the number we divide by is ours or the product's.

An **unregistered** holding is **not** assumed normalised: assuming it would apply unity
as a continuum and inflate every EW silently. Absence of a claim is not a claim (RYA-833).

## 6. Scope

The **"consistent"** tier is in the schema but its handling is **deferred** (Ryan).
No value moves here — this ratifies structure, display and guards. Production runs per
(instrument × band) are the children: RYA-959 (VIS, in flight), RYA-961, RYA-908,
RYA-953, RYA-1027, RYA-1028, plus the 3D-column umbrella RYA-1029.

Glossary / method / architecture doc touches accumulate into **RYA-179**, not piecemeal.
