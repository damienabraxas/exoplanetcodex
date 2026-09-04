# RYA-850 — the graded (lab-gf) pool as the reported value

**Status:** wiring + reporting. No new measurement. **Not merged — Ryan reviews.**

> **Updated after RYA-853 (`eb02eb4`), corrected on re-referee.** The graded pool this
> ticket promotes is **refereed against the source papers and CLEAR** — Ruffoni 142/142,
> Den Hartog 203/203 and Belmonte 119/120, **464 of 465 exact on value and on cited σ, zero
> mismatches**. ⚠️ The "one genuine bad row" this note used to cite (Belmonte 3935.307) is
> **WITHDRAWN**: it was the referee reading Belmonte's *Published* (May et al. 1974)
> comparison column instead of its *This Experiment* column. There is no bad row, so the
> question of whether it reaches a cell here is moot. The 70%-wrong NIST extracts are a
> *different* pool and do not touch these cells.

```
python3 scripts/rya850_graded_products.py
```

---

## What was actually missing

The lab-gf pools were built (RYA-824 VIS/IR, RYA-836 near-UV), and
`error_budget.gf_term(graded=...)` has carried both terms all along:

| | dex | basis |
|---|---|---|
| graded | **0.041** | NIST grade B = 10% on the transition probability |
| ungraded | **0.170** | Kurucz semi-empirical, 0.1–0.3 dex (RYA-161) |

🔴 **`derive_band_products.py:770` passes `gf_graded=False` unconditionally.** Every
EW-route cell has been charged the Kurucz 0.17 regardless of what its lines actually are,
and RYA-824's two pools were never emitted as cells at all.

---

## The result: every graded cell tightens its bar

| band | ion | ungraded total | **graded total** | |
|---|---|---|---|---|
| VIS | I | 0.1732 | **0.0730** | ✅ |
| red-optical | I | 0.1767 | **0.0658** | ✅ |
| near-UV | I | 0.2077 | **0.1412** | ✅ |

**Main-table view (± visible, per spec item 4):**

```
A(Fe I; VIS,         1D-LTE) = 7.445 +/- 0.073   [stat 0.039, sys 0.061, n=9]
A(Fe I; red-optical, 1D-LTE) = 7.516 +/- 0.066   [stat 0.022, sys 0.062, n=20]
A(Fe I; near-UV,     1D-LTE) = 7.577 +/- 0.141   [stat 0.085, sys 0.113, n=59]
```

### This is not in tension with RYA-836 or RYA-842, and the reason is arithmetic

RYA-836 measured the near-UV lab-gf pool scattering **worse** (0.652 vs 0.413), and RYA-842
ratified that line selection — not gf — drives value and spread. **Both still hold.** The
reported bar is

    total = sqrt(stat² + sys²),   stat = scatter/√N

so the worse scatter enters only through `stat`, which **averages down**, while the gf term
enters `sys`, which **does not**. A pool can therefore scatter more and still report a
smaller total.

That is a consequence, not a re-litigation — so it is **checked per cell**, and a graded
cell that failed to beat its twin would be flagged rather than promoted. All three pass.

---

## The graded bands agree with each other

If they disagreed by more than their bars, "the graded value" would not be well defined and
picking one would be a choice dressed as a measurement.

| pair | Δ | combined σ | |
|---|---|---|---|
| VIS vs near-UV | −0.132 | 0.159 | **0.8σ** |
| VIS vs red-optical | −0.071 | 0.098 | **0.7σ** |
| near-UV vs red-optical | +0.061 | 0.156 | **0.4σ** |

Band-to-band spread **0.132 dex**, all pairs consistent. Notably the red-optical graded
value, **7.516**, lands exactly on the 1D reference (RYA-669/783); gold's 7.466 is on the
3D scale and is not the comparison.

---

## 🔴 Only 3 of 18 cells can have a graded twin

This is a finding, not a shortfall of the wiring. A primary-lab-gf pool exists **only** for
Fe I in the near-UV (RYA-836, n=59), Fe I VIS (RYA-824, n=9) and Fe I red-optical (RYA-824,
n=20).

**There is no lab-gf pool for Fe II in any band, nor for any ENGINE-A/ENGINE-B product.**
Nobody has measured one. Those 14 cells fall back to UNGRADED and say so, rather than being
relabelled graded on the strength of a term they are not entitled to.

⚠️ **The ungraded product is kept, never dropped** — an ungraded line is usually good data
with a worse-known oscillator strength. It stays labelled, with its wider bar, as the BROAD
number over every measurable line (RYA-712).

---

## Three things the spec asked for that need a decision

**1. ✅ DECIDED (Ryan, 2026-08-17): the pool's own CITED laboratory σ sets the published
bar.** RYA-850's spec text named `graded_gf_term` (0.041) twice, so this is a deliberate
deviation and the ticket description is amended to match — code that contradicts its own
ticket is how a reader loses the ability to tell a decision from a drift.

**Why the bound loses.** `GRADED_GF_SYSTEMATIC_DEX` is the worst grade we *accept*, and a
bound is the right answer only while the actual σ are unknown. These pools publish per-line
uncertainties — a measurement of the same quantity — and RYA-853 refereed them line-by-line
against the source papers before they were allowed near a bar.

| band | n | cited σ | vs generic | total, bound → **published** |
|---|---|---|---|---|
| VIS | 9/9 | **0.0600** | 1.5× | 0.0583 → **0.0730** (+25.1%) |
| red-optical | 20/20 | **0.0524** | 1.3× | 0.0570 → **0.0658** (+15.3%) |
| near-UV | 58/60 | **0.0522** | 1.3× | 0.1375 → **0.1412** (+2.7%) |

The near-UV barely moves because RYA-841's 0.100 dex pseudo-continuum term swamps any gf
term there — so the choice only bites on the two EW-route bands.

**Two guards changed the numbers, and both were worth having:**

🔴 **The near-UV σ is 0.0522, not the 0.0518 I first reported.** A wavelength window is not
a unique line ID (RYA-853), and 2 of the 60 lines match two lab rows each. Resolving them
with `iloc[0]` is precisely how RYA-853 manufactured 12-dex "defects", so they are counted
**unmatched** — n=58/60.

⚠️ **A partly-covered pool refuses the cited term.** The RMS describes the pool only if it
covers the pool; otherwise the unmatched lines silently inherit the matched ones'
uncertainty. Below 90% coverage the term is refused and the bound stands, which is the
honest fallback rather than a tighter number.

`cited_gf_term` lives in `error_budget`, not in this script, because any element that
acquires a lab-gf pool inherits it. It **replaces** the generic term rather than joining it
(both describe the oscillator strengths, so carrying both double-counts), refuses an
ungraded pool, refuses an unsourced σ, and is **not clamped** to the bound — a grade-A pool
would legitimately fall below 0.041, and clamping would turn the measurement back into the
assumption it supersedes.

**2. No combined headline is computed.** The ticket asks for a "combined/headline Fe value".
**RYA-712 forbids a cross-band combined product** and `pipeline.band_products` deliberately
has no `combine()`; the per-band values are separate measurements whose spread *is* the
result. So this emits the per-band graded table and reports the spread. Collapsing three
consistent bands into one number is a reporting decision for a human — the natural
candidate is the red-optical cell (tightest bar, largest lab-gf pool of the two EW-route
bands, and on the 1D reference).

**3. ⚠️ Still pre-RYA-847 — survivable, not blocking.** 847 is 6 of 7 items done but a
deliberate **numerical no-op today**: `SYNTH_CONSTRAINT = None`, and its own near-UV control
run excludes **0 lines**, because the threshold is set by a sweep still in progress. So
nothing here contradicts it.

When the sweep ratifies a cut, the **near-UV cell is the exposed one** (synthesis route) and
must be regenerated. 847's scope note puts the EW-route products — naming `1D-LTE-LABGF in
VIS/IR` explicitly — outside its remit, so the VIS and red-optical cells stand.

⚠️ One file, `scripts/derive_band_products.py`, is in both tickets' scope. 847 has not
touched it, so landing this first means 847 rebases onto it rather than the reverse.

---

## The reporting layer is general

`pipeline/graded_reporting.py` keys on the pool marker (`-LABGF`), not on Fe, so Al and the
rest inherit it: `is_graded`, `base_treatment`, `pair_products`, `element_table`,
`format_value`. It performs exactly one combination — `sqrt(stat² + sys²)` — and adds no
statistics beyond stat/sys/total.
