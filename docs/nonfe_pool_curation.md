# Non-Fe EW line-pool curation (RYA-395)

GES-style quality vetting of the **non-Fe** solar EW line pools, the piece that
gates RYA-371's non-Fe GO. Closes the gap left by RYA-239: every measured non-Fe
element read systematically high (+0.17 … +2.25 dex vs Asplund 2021) because the
non-Fe pools never got the curation Fe I/II received.

```bash
python -m pipeline.curate_nonfe_pools --phase 1 --verify   # Mg Si Ca Ni
python -m pipeline.curate_nonfe_pools --phase 2 --verify   # Ti Cr Na Al Mn
python -m pipeline.curate_nonfe_pools --phase all
python -m pipeline.curate_nonfe_pools --grade-restrict --verify   # RYA-398 independent-gf pass
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

## The independent-gf (grade) cull — `--grade-restrict` (RYA-398)

The 395 cull removed SAT/blend/noise but **not gf-grade**, so the surviving pools
were 43–81 % ungraded/Kurucz gf — which is why the full pool reads +2.19 on Cr
while RYA-239's graded-region subset read +0.345. `--grade-restrict` adds the
missing cull: a `GRADE` reason that keeps only **independent quality-graded gf**
(`ACCEPTED_GF_TIERS = {HIGH, MED}` — a NIST grade, or a non-Kurucz literature/lab
reference) and culls Kurucz semi-empirical (K03–K10, KP) and unreferenced gf.

It is **abundance-blind** (gf grade is atomic-data metadata, not an abundance) —
`--verify` shuffles A(X) and asserts the grade-cull mask is byte-identical. On the
graded subset we recompute 1D-LTE A(X), apply the registered NLTE per element, and
read the verdict against the fixed Asplund band `max(2σ_Asplund, 0.10)`:

* **VALIDATED** — graded gf + NLTE recovers Asplund within tolerance.
* **RESIDUAL** — gross offset gone, a residual survives (Cr ≈ +0.40 after NLTE):
  a real open puzzle (1D-LTE vs 3D, residual line quality). **Reported, never tuned.**
* **LOW_CONFIDENCE** — too few independent (graded) lines for a stable mean.

## The validation↔survey firewall (RYA-398)

There are **two gf paths and they must not cross**:

| path | gf source | why |
|------|-----------|-----|
| **solar validation** (this module) | independent (graded) gf only | Asplund is the *check*; reproducing it must come from independent atomic data, or the validation is circular. |
| **differential survey** (55 Cnc / α Cen) | RYA-161 astrophysical gf | 161 inverts the COG with the *known* Asplund abundance (`log_gf_astro = invert_COG(EW, A=Asplund)`), so it reproduces Asplund **by construction** — fine for [X/H], where the common-mode Sun↔target gf error cancels, but it is tuning if it touches the solar validation. |

`assert_no_astrophysical_gf()` is a loud guard (raises, never warns) on the
solar-validation path: if it is ever handed a `log_gf_astro` / `delta_log_gf` /
`gf_astro` column it refuses to proceed. **A residual on independent gf is
information — the next puzzle — never something to erase by switching to
solar-calibrated gf.**

## Outputs

`data/curation/nonfe_pools/`
* `{element}_cull_rya395.csv` / `{element}_cull_graded_rya398.csv` — per-line
  KEEP/CULL + reason, provenance header (full pool / graded subset).
* `curation_diagnostics_rya395.csv` / `curation_diagnostics_graded_rya398.csv` —
  per-element A(X) LTE, scatter, REW slope, N_clean, NLTE Δ, residual, verdict.

Data sources (single-source, nothing hardcoded): measured EW =
`data/measured/sol_ew_results_v1.csv`; atomic data (loggf/EP/nist_grade/
loggf_reference) = `data/linelists/canonical_gf.csv`; RT species codes/damping =
the GES synthesis linelist (loggf/EP overridden from canonical); NLTE Δ = the
vendored MPIA/INSPECT/CDS grids queried at solar params.
