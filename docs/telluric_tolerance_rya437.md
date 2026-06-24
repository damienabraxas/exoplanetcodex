# RYA-437 — Calibrating the telluric verification gate from the 13C/CO precision requirement

RYA-424 shipped the telluric data-input gate at `TELLURIC_RESIDUAL_TOL = 0.05` — a
provisional round number with no provenance. The Vesta CRIRES+ smoke landed at
0.061 / 0.065 and was correctly *flagged*; the gate **behaviour** was right (loud, not
silent) but the **threshold** was unvalidated. RYA-437 derives the real tolerance from
the science (validate-don't-tune: the number comes from the physics, never from wanting
a specific frame to pass), adds a CO-region-local metric, and RCAs Vesta.

Code: `pipeline/telluric_tolerance_rya437.py` (derivation), `config/constants.py`
(`TELLURIC_RESIDUAL_TOL` + `…_PROVENANCE` + `TELLURIC_CO_LOCAL_WINDOW_A`),
`pipeline/telluric_stage.py` (CO-local metric + gate wiring).

## A. Residual → 13C bias propagation, and the derived tolerance

A telluric mis-correction of typical amplitude `r` (the residual metric = median
`|1 − corrected_flux/continuum|` at telluric-dominated pixels) perturbs the apparent
**depth** of a CO feature by `~r` in continuum-normalized flux. For a weak line on the
**linear** part of the curve of growth, `EW ∝ N` (the abundance), so a depth
perturbation `r` on a feature of central depth `d`:

```
ΔEW/EW = (r·W)/(d·W) = r/d            (feature width W cancels)
ΔA     = (ΔEW/EW)/ln10 = r/(d·ln10)   (abundance bias, dex)
```

Inverting at the precision target `σ` (dex):  **`r(feature) = σ · ln10 · d(feature)`**.

### Per-feature (the binding one sets the gate)

The 12C/13C ratio `R = N(12C)/N(13C)` is carried by the **13CO** leg: 12CO is the
well-determined numerator (strong, deep), 13CO is the weak/uncertain denominator. The
ratio error is dominated by `ΔN13/N13 = r/d13`, so the gate is set by the **13CO**
feature — the shallower depth makes it the tighter (binding) tolerance.

| feature | λ_vac | central depth d | dA/dr = 1/(d·ln10) | r = σ·ln10·d |
|---|---|---|---|---|
| 12CO(2-0) bandhead | 2.2935 µm | 0.140 | 3.11 dex/unit | 0.0161 (looser) |
| **13CO(2-0) bandhead** | **2.3448 µm** | **0.093** | **4.68 dex/unit** | **0.0107 (BINDING)** |

**`TELLURIC_RESIDUAL_TOL = 0.0107` (1.1%), set by 13CO(2-0).**

### Inputs (cited — RYA-436 discipline)

- **Precision target σ = 0.05 dex** — the project's carbon-abundance precision, the
  Asplund-2021 solar carbon gate `cno_synthesis.SOLAR_VIS_GATES['C'] = (8.46, 0.05)`.
  The 13C abundance is a carbon-isotopologue abundance, so the carbon precision is the
  natural target (and tighter than C/O's 0.08). Implied 12C/13C fractional precision =
  `ln10·σ ≈ 11.5%`, consistent with high-quality stellar CO isotopic work (~10–20%).
- **Feature depths** — measured from the **ACE-FTS solar atlas** (Hase, Wallace, McLeod,
  Harrison & Bernath 2010, JQSRT 111, 521; telluric-free, roughly disk-integrated →
  matches reflected-solar Vesta / the integrated-disk solar Phase-B target), the segment
  vendored at `data/solar_reference/ir_atlases/ace_fts_solar_co_4255_4367.csv` (RYA-390).
  Disk-integrated depths are the conservative (shallower → tighter-gate) choice vs the
  ground-based NSO photatl cross-check (which gives d13 ≈ 0.12 → a looser 0.0138).

### Why this matters — the bias the round 0.05 would have admitted

`ΔA(13C) = r / (d13·ln10)`:

| telluric residual | ΔA(13C) bias | note |
|---|---|---|
| 1.1% | +0.051 dex | the calibrated gate (≈ the 0.05 target, by construction) |
| 2.0% | +0.094 dex | |
| 5.0% (old round gate) | +0.234 dex | would admit frames that cannot do the isotopic science |
| 6.0% (Vesta) | +0.281 dex | 12C/13C wrong by ~1.9× |

The calibrated gate is **~5× tighter** than the round 0.05. Reproduce:
`python -m pipeline.telluric_tolerance_rya437`.

## B. CO-region-local residual metric

Telluric-correction quality is wavelength-dependent — a global median can pass while the
CO region is dirty, or vice versa. The residual metric now reports **two** numbers
(`pipeline/telluric_stage.telluric_residual_metric(..., window_A=)`):

- **CO-local (primary for CO science)** — scored only within
  `TELLURIC_CO_LOCAL_WINDOW_A = (22900, 23500)` Å, i.e. 2.2935 µm (12CO) through
  2.3448 µm (13CO).
- **global (secondary)** — the whole corrected segment, retained for context.

The gate's verdict uses the CO-local residual when the science is CO; both land in the
conditioning manifest (`residual_co_local`, `residual_global`, `metric_used`).

## C. Vesta correction-quality RCA

Outcome: **FLOOR / COVERAGE-LIMITED** — Vesta does not pass the calibrated gate, and the
limiting factors are named. Validate-don't-tune: the gate was not loosened. Reproduce:
`python -m scripts.vesta_telluric_rca_rya437`.

1. **Coverage (dominant).** The binding science feature, 13CO(2-0) at 2.3448 µm, is
   **not on-chip** in either on-chip CO frame — it falls in a CRIRES+ inter-order/detector
   gap. K2217's 12CO segment covers 22835–22975 Å; K2192's 22879–23012 Å; neither reaches
   23448 Å, and no other segment covers it either. So K2192/K2217 **cannot measure 13C at
   all**, independent of telluric quality. A frame with 13CO(2-0) on a detector is a
   prerequisite for solar Phase-B / alpha Cen 13C.

2. **Telluric residual floor.** The measurable 12CO-region CO-local residual is **0.0317**
   on the reliable frame (K2217, 121 px) — ~3× the calibrated 1.07% gate. (K2192 reads
   0.081 but on only 10 px, its 12CO segment barely overlapping the CO-local window, so it
   is statistically marginal; K2217 is representative.)

   - *Not solar contamination of the metric*: the telluric-only scored pixels (0.0317) and
     the solar-coincident pixels (0.0509) are comparable — the residual is genuine telluric
     misfit on shallow K-band telluric lines, not the solar CO forest leaking into the score.
   - *Not a continuum-order fix*: a CONTINUUM_N sweep (1/2/3) moves the residual by 0.0003 —
     continuum-order-INSENSITIVE. An aggressive variant (CONTINUUM_N=5 + FIT_WLC + tight
     FTOL/XTOL) made it *worse* (overfit). So ~3% is a molecfit line-fit / data floor at
     this configuration, not a setup knob away from passing.

The consequence of Part A made concrete: at Vesta's best 12CO-region residual (3.2%) the
13C bias would be 0.032/(0.093·ln10) ≈ **0.15 dex** (and 13CO is not even on-chip). So the
calibrated gate names the **real precision ceiling**: the CRIRES+ K-band CO arm needs both
(a) frames with 13CO(2-0) on a detector and (b) ~1% telluric residual before it can deliver
12C/13C to the 0.05-dex carbon precision target. The round 0.05 would have passed K2217 (at
~6% global / 3.2% CO-local) as "analysis-ready" for a measurement it cannot support.

(`_molecfit_segment` gained backward-compatible `continuum_n` / `extra_cmd` RCA hooks so
this sweep is reproducible; the RYA-373 defaults are unchanged.)
