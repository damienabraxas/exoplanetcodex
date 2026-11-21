# RYA-1354 v68: Mg I spin–orbit calibrated on the ion, second dipole form, full error budget

Status: **CANDIDATE / HOLD** (RYA-1352). Declared on RYA-1354 (2026-11-20 06:12 UTC) before the run.

Script: `scripts/rya1354_twovalence_model.py 'Mg I'`. Variants: `--variant=NAME` (10 runs, written to `data/audit/rya1354/mg1_model_v68/variants/`). Outputs are in `data/audit/rya1354/mg1_model_v68/`. Every line in `lines.csv` carries its own error budget, one column per term.

## 1. Spin–orbit calibrated on Mg⁺ only (declared primary)

Per-l scale s_l = median NIST/model fine structure of Mg II levels with splitting > 0.1 cm⁻¹:

| l | Scale |
|---|---|
| p | 1.222 |
| d | 1.100 |
| f | 0.958 (g takes the f value) |

No Mg I information enters.

| Blind prediction | Ab initio | **Ion-calibrated** | Reference |
|---|---|---|---|
| 3s3p ³P / 3p² ³P / 3p3d ³D splitting | 77 / 83 / 82 % | **94 / 102 / 100 %** | NIST |
| **4571 Å** | −5.824 | **−5.650 (+0.036)** | **lab −5.686** (Kwong 1982); theory −5.623 (TFF) |
| 4f ³F splitting | 100 % | 118 % | NIST |
| 14878 / 12083 Å (lab, spin-changing) | −0.078 / −0.071 | +0.142 / +0.149 | Pehlivan Rhodin 2017 |

Physics:
- The ion calibration fixes the p shell, which sets 4571 Å.
- It over-mixes 4f. There the ³F splitting comes from 3p3d admixture in the CI state, so the raised p and d scales inflate it, and the mixing scales as its square (2 log 1.18 = 0.14 dex = the miss).

**Variants (reported in `lines.csv` as `loggf_so_*`):**
- `abinitio`: unscaled.
- `tuned`: every resolved triplet fine-tuned to NIST energies, the Brage & Hibbert method used by Pehlivan Rhodin 2017. This overshoots 4571 Å (+0.104): the measured splitting also contains two-body magnetic terms that do not mix singlets and triplets.
- `hybrid` (**post-hoc, not adopted**): ion-calibrated, with energy fine-tuning only on the LS-broken configurations (nf, ng). All three lab spin-changing lines come within 0.04 (4571 +0.036, 14878 −0.001, 12083 +0.006). To be declared and tested blind on Ca I intercombination lines (Den Hartog 2021).

## 2. Second dipole form (velocity)

- **Control:** with one shared potential, velocity = ΔE × length to 10⁻⁶ (also a pytest).
- With our l-dependent potential the non-local offset is 0.03–0.08 % for Mg II valence states.
- **Lab lines:** dT ≤ 0.009 on all 26. Spearman(dT, |miss|) = −0.08, so the test has no power on this set.
- **Line list** (spin-allowed, n ≤ 8): ρ = 0.30 overall, but only through CF (dT–CF ρ = −0.50; within CF ≥ 0.1, ρ = −0.11).
- **Result:** dT is a diagnostic, not used in σ (as declared). Median dT is 0.004 and p90 0.25 (cancellation lines).

## 3. Convergence (6 full-CI rebuilds) and 4. physics sensitivity (4)

|Δ log gf| at p68 / p90, for lines with CF ≥ 0.1:

| Variant | n ≤ 8 | n 9–10 | n 11–12 | Energy max err |
|---|---|---|---|---|
| Rydberg 14 → 18 | 0.000 / 0.001 | 0.000 / 0.011 | 0.001 / 0.004 | 0.45 % |
| compact 14 → 18 | 0.001 / 0.011 | 0.001 / 0.012 | 0.001 / 0.020 | 0.45 % |
| l_max 4 → 5 | 0.001 / 0.010 | 0.001 / 0.002 | 0.001 / 0.002 | 0.45 % |
| grid 3200 → 4800 | 0.000 / 0.001 | 0.000 / 0.004 | 0.000 / 0.013 | 0.45 % |
| box 900 → 1200 a₀ | 0.000 / 0.001 | 0.000 / 0.002 | 0.000 / 0.009 | 0.45 % |
| **compact box 25 → 35 a₀** | 0.003 / 0.066 | 0.006 / 0.180 | 0.008 / **0.339** | **0.87 %** |
| α_d ± 5 % (two-electron terms) | 0.002 / 0.007 | 0.002 / 0.006 | 0.001 / 0.005 | 0.44–0.48 % |
| dielectronic cutoff r_c ± 20 % | 0.007 / 0.022 | 0.006 / 0.020 | 0.004 / 0.022 | 0.43–0.55 % |

**The model is converged.** Every knob that refines the basis moves lines by ≤ 0.02 dex (p90).

The physics inputs α_d and r_c matter at the 0.02 dex level.

The one large mover, compact box 35 a₀, **degrades** the basis: the energy error doubles. It is therefore not a convergence step. As declared, σ_conv still takes it, so σ for high-n lines is conservative. Proposed for v69: replace it with a tighter compact variant (r = 20 a₀, 18 functions).

## Error budget and grade

σ = hypot(floor 0.050, ½ polarization operator, family/lab term, max(spin–orbit treatment, physics), convergence) ⊕ CF-bin term. CF < 0.02 gets no σ.

Spin-allowed, n ≤ 8, CF ≥ 0.1 lines (n = 747):

| Term | Median (dex) | p90 (dex) |
|---|---|---|
| convergence | 0.0005 | 0.036 |
| physics | 0.0016 | 0.014 |
| spin–orbit | 0.0002 | 0.025 |
| **σ** | **0.055** | 0.136 |

σ is set by the lab floor.

Spin-changing lines carry the spin–orbit treatment term (~0.17 dex), so most weak intercombination lines are class E. Line classes: C 1,318, D 96, E 2,803.

Lab (spin-allowed, 24 lines): median +0.020, abs68 0.022. Coverage (in-sample): 83 % within 1σ, 96 % within 2σ.

## Other findings this stage

- **NIST level composition** (`data/reference/onevalence_nist/MgI_levels_percentages.csv`, new):
  - Every 3s nd ¹D₂ and ³D₂ level from 4d to 8d is **100 % pure**. This is gf-free confirmation that VALD's spin-changing nd values (which need ~50 % mixing) are wrong.
  - CI 3p² admixture vs NIST: 3s² 5.3 / 7 %; 3s3d ¹D **16.6 / 22 %** (under-mixed); 4d 11.4 / 13; 5d 7.0 / 7; 6d 4.3 / 4; 7d 2.8 / 3 %. The 3p3d admixture in 3s3p ¹P is 4.8 / 6 %.
  - Percentages depend on the orbital basis, so the comparison is approximate.
- **5711 Å, the line that sets the floor:** our model −1.703, Pehlivan's own calculation −1.742 and KP08/Chang–Tang CI −1.724 agree; the lab value −1.842 is the outlier. It is kept as an anchor and flagged.
- **RYA-1393 Mg I audit** (`mg1_model_v68/rya1393_intercomb_audit/`), 52 spin-changing features:
  - 10 REFUTED in the Sun (VIS/Y/J/H/K; the model is weaker in all 10).
  - 4 CONFIRMED (incl. 4571).
  - 15 blended, 22 with no atlas (UV).
- **Tests:** `tests/test_rya1354_twovalence_model.py`, 11 pass (Landé interval rule, l·s, fine-structure sum rule, gauge identity, spin–orbit linearity, LS-broken set, parity labels).

## Open for Ryan

None blocking. The hybrid spin–orbit rule goes to a blind test on Ca I (declared first).

Next is v69 per the roadmap: solar multiplet line ratios as an abundance-free referee, sister ions Al II / Si III blind, the NIST composition test formalised, and a replacement for the compact-box variant.
