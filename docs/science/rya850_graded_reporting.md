# RYA-850 — the graded (lab-gf) pool as the reported value

**Status:** wiring + reporting. No new measurement. **Not merged — Ryan reviews.**

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
| VIS | I | 0.1732 | **0.0583** | ✅ |
| red-optical | I | 0.1767 | **0.0570** | ✅ |
| near-UV | I | 0.2077 | **0.1375** | ✅ |

**Main-table view (± visible, per spec item 4):**

```
A(Fe I; VIS,         1D-LTE) = 7.445 +/- 0.058   [stat 0.039, sys 0.043, n=9]
A(Fe I; red-optical, 1D-LTE) = 7.516 +/- 0.057   [stat 0.022, sys 0.052, n=20]
A(Fe I; near-UV,     1D-LTE) = 7.577 +/- 0.137   [stat 0.085, sys 0.108, n=59]
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
| VIS vs near-UV | −0.132 | 0.149 | **0.9σ** |
| VIS vs red-optical | −0.071 | 0.082 | **0.9σ** |
| near-UV vs red-optical | +0.061 | 0.149 | **0.4σ** |

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

**1. The generic 0.041 understates these pools.** RYA-824 recorded what its lines' *cited*
laboratory σ actually are:

| band | generic term | cited (RYA-824) | total, generic → cited |
|---|---|---|---|
| red-optical | 0.041 | **0.0524** | 0.0570 → 0.0657 (**+15%**) |
| VIS | 0.041 | **0.0600** | 0.0583 → 0.0729 (**+25%**) |

The spec asks for the generic term, so that is what is wired — but the measured figure is
the more defensible one, and using the generic value publishes a bar 15–25% tighter than
the pool's own gf uncertainty supports.

**2. No combined headline is computed.** The ticket asks for a "combined/headline Fe value".
**RYA-712 forbids a cross-band combined product** and `pipeline.band_products` deliberately
has no `combine()`; the per-band values are separate measurements whose spread *is* the
result. So this emits the per-band graded table and reports the spread. Collapsing three
consistent bands into one number is a reporting decision for a human — the natural
candidate is the red-optical cell (tightest bar, largest lab-gf pool of the two EW-route
bands, and on the 1D reference).

**3. ⚠️ The line sets are pre-RYA-847.** RYA-847 removes unconstrained synthesis fits, after
which the graded pool is lab-gf **∩** constrained. It is unmerged and actively in flight, so
**every number here is provisional**. The near-UV cells are the exposed ones (synthesis
route); the VIS/red-optical lab-gf pools are EW-route and outside 847's scope, so those two
should survive unchanged.

---

## The reporting layer is general

`pipeline/graded_reporting.py` keys on the pool marker (`-LABGF`), not on Fe, so Al and the
rest inherit it: `is_graded`, `base_treatment`, `pair_products`, `element_table`,
`format_value`. It performs exactly one combination — `sqrt(stat² + sys²)` — and adds no
statistics beyond stat/sys/total.
