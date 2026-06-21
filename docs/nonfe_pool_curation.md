# Non-Fe EW line-pool curation (RYA-395)

GES-style quality vetting of the **non-Fe** solar EW line pools, the piece that
gates RYA-371's non-Fe GO. Closes the gap left by RYA-239: every measured non-Fe
element read systematically high (+0.17 … +2.25 dex vs Asplund 2021) because the
non-Fe pools never got the curation Fe I/II received.

```bash
python -m pipeline.curate_nonfe_pools --phase 1 --verify   # Mg Si Ca Ni
python -m pipeline.curate_nonfe_pools --phase 2 --verify   # Ti Cr Na Al Mn
python -m pipeline.curate_nonfe_pools --phase all
```

## The cardinal rule — validate, don't tune

Asplund 2021 is the **validation** target, never the tuning target. Every cull is
made on spectroscopic-quality grounds with thresholds **fixed in advance**
(`pipeline/curate_nonfe_pools.py` module constants), uniform across elements, and
applied **blind** to the abundance. A(X) is read only afterwards, to validate.
Enforced structurally: `cull_reasons()` is only ever handed the quality columns
(`_CULL_INPUT_COLS`) and raises if shown an abundance column; `--verify` shuffles
A(X) and asserts the cull mask is unchanged.

## Fixed quality criteria (abundance-blind, uniform)

| reason | test | physical basis |
|--------|------|----------------|
| `WEAK`  | `ew_mA < 5` | below detection floor (`PIPELINE['ew_min_mA']`) |
| `SAT`   | `log(EW/λ) > −4.9` (linear-COG knee) or `ew_mA > 100` | off the linear COG → EW→A biased high & vmic-sensitive |
| `HIERR` | `ew_err/ew > 0.5` | measurement too noisy (`PIPELINE['fe2_ew_err_frac_max']`) |
| `BLEND` | measured `blend_flag`, or strength-weighted contaminant opacity > 0.5 | genuine spectroscopic blend |
| `BADGF` | `nist_grade ∈ {D,E,F}` | graded-unreliable lab gf |

**Blend vetting** is a strength-weighted opacity-overlap metric — the principled
fix to the old `vald_proximity_flag`, which flagged on *any* nearby VALD line
regardless of strength (tagging ~87 % of lines, the known silent bug). Each
in-window neighbour is weighted by its solar-reference line opacity
(`A_sun + loggf − EP·θ`) and a Gaussian profile kernel; we cull only when
contaminants supply > 50 % of the target line's in-profile opacity. It uses the
*reference* solar composition (an input), never the derived A(X).

## The cull-vs-gf-scale fork

After quality curation, a persisting offset is a **gf-scale systematic** (RYA-161
territory), not stray bad lines — so it is **flagged for escalation, never
tuned**. The evidence is `low_gf_frac`: the fraction of the pool whose gf is
ungraded / Kurucz-theoretical (K03–K10, KP). Verdicts: `CLEAN` (reached
Asplund − Δ_NLTE on quality grounds), `GF_SCALE_RYA161` (offset persists →
escalate), `LOW_CONFIDENCE` (too few clean lines).

**Cr / Mn canaries.** Their NLTE Δ is *positive* (Cr +0.072, Mn +0.106, queried
live from the RYA-396 grids), so curation must bring LTE to **(Asplund − Δ_NLTE)**
— not to Asplund directly, or the positive NLTE overshoots. `--verify` asserts
that contract (`target_lte + Δ_NLTE == Asplund`, `Δ_NLTE > 0`).

## Outputs

`data/curation/nonfe_pools/`
* `{element}_cull_rya395.csv` — per-line KEEP/CULL + reason, provenance header.
* `curation_diagnostics_rya395.csv` — per-element A(X) LTE, scatter, REW slope,
  N_clean, NLTE Δ, target, down-correction, residual, low-gf-fraction, verdict.

Data sources (single-source, nothing hardcoded): measured EW =
`data/measured/sol_ew_results_v1.csv`; atomic data (loggf/EP/nist_grade/
loggf_reference) = `data/linelists/canonical_gf.csv`; RT species codes/damping =
the GES synthesis linelist (loggf/EP overridden from canonical); NLTE Δ = the
vendored MPIA/INSPECT/CDS grids queried at solar params.
