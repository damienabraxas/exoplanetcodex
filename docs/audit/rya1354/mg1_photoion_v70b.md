# RYA-1354 v70b: Mg I photoionization from our own wavefunctions

**Status: CANDIDATE / HOLD, report only.** No synthesis, NLTE or product input is changed. The milestone was declared on RYA-1354 on 2026-11-21 (comment `257df6eb`) before any code was written or run.

Code:
- `scripts/rya1354_photoion.py`: the box runs, the C0 control and the T4 polarizability.
- `scripts/rya1354_photoion_analyze.py`: pooling and the tests.
- `scripts/rya1354_twovalence_model.py`: generalised so it can take a box continuum. Without one, it is regression-identical to v70.

Outputs: `data/audit/rya1354/mg1_photoion_v70b/`.

## Verdict

| | Declared band | Result | |
|---|---|---|---|
| C0: bound model intact | \|Δlog gf\| ≤ 0.01 (24 lab lines); energy ≤ v70 + 0.05 % | 0.0001 dex; 0.450 % (v70 0.450 %) | **PASS** (attempt 2; see below) |
| C1: box independence | spread ≤ 10 % at the test energies | 3.2 % (T2), 0.4 % (T3) | **PASS** |
| C2: gauge | reported; enters σ_model | velocity +3 % (T2/T3) | reported |
| T2: 3s3p ¹P, 354.5 nm | \|z\| ≤ 2 against 90 ± 16 Mb | **84.5 ± 3.0 Mb**, z = −0.34 | **PASS** |
| T3: 3s3p ¹P, 285 nm | \|z\| ≤ 2 against 81 ± 23 Mb | **67.0 ± 1.1 Mb**, z = −0.61 | **PASS** |
| T1: 3p² ¹S₀, position | \|E_r − 68,270\| ≤ 300 cm⁻¹ | **68,124** (−146) | **PASS** |
| T1: width | factor 1.5 of 280 cm⁻¹ | **608 ± 7 cm⁻¹** (×2.2) | **FAIL** |
| T1: peak σ | 785 ± 2×141 Mb | **185 Mb** | **FAIL** |
| T4: α₀(3s²) | ±2 % of 71.33 | **74.78** (+4.8 %) | **FAIL** |
| T4: α₀(3s3p ³P) | ±2 % of 101.2 | **107.41** (+6.1 %) | **FAIL** |
| T5: bound–free continuity | ratio within 15 % | 16 / 61 channels | **FAIL as declared**; see the breakdown |

## Method

- **Basis.** The v70 basis is kept on a finer grid (6,400 points), and box-continuum orbitals are added for l = 0–4.
  - Continuum orbitals are positive-energy eigenstates of the screened outer-electron potential inside a wall at R_box, with 0 < ε ≤ 0.115 hartree. There are 48–49 per l at R = 400.
  - They are projected against the core states and the v70 span, then re-diagonalised in h_l.
- **Hamiltonian.** The one-body matrix is no longer diagonal, so `hamiltonian()` now takes `B.h1`, and Coulomb blocks are computed on index subsets.
  - Regression: identical to the v70 matrices to 1e-16 on four symmetries of the real Mg I basis. There is also a unit test (`test_v70b_general_hamiltonian_matches_v70_block_algebra`).
- **Configurations.** All v70 pairs, plus Mg⁺ (3s, 4s, 3p, 4p, 3d) × εl.
- **Cross-section.** Box normalisation: df/dε = f_n / spacing. σ = 4.0336 Mb × df/dε, with ε in hartree.
- **Boxes.** Eight radii, R = 360–430 a₀, pooled; the pooled spectrum samples the resonance about every 50 cm⁻¹.
- **Gauges.** Length gauge with the v70 core-polarization dipole; velocity gauge as the bare ∇.
- **Energy coordinate.** ε is measured above the model's Mg⁺ 3s threshold. Lab photon energies are converted with NIST, so thresholds are aligned as in NLTE practice (Lind & Amarsi 2024).

## C0: attempt 1 FAILED_RETAINED, then a harness fix

- Attempt 1 at R = 400 gave an energy max error of 5.9 %. Files: `c0_400_attempt1_FAILED.*`.
- **Cause.** NIST lists 3p² ¹S at +6,604 cm⁻¹ above the limit; it is the T1 resonance. The C0 harness kept only roots below the limit, so the aligner forced that term onto a box state.
- **Fix.** C0 compares bound NIST terms only, in both runs. Attempt 2 passes: 0.450 % with or without the continuum, and the gf are unchanged.

## T2/T3: absolute cross-sections from 3s3p ¹P (PASS)

| λ (nm) | ε above limit (cm⁻¹) | Model, length (Mb) | Model, velocity (Mb) | Box spread | Lab (Mb) | z |
|---|---|---|---|---|---|---|
| 354.5 | 1,581 | 84.5 | 87.1 | 3.2 % | 90 ± 16 (Rafiq 2007) | −0.34 |
| 285 | 8,458 | 67.0 | 69.2 | 0.4 % | 81 ± 23 (Madsen & Thomsen 2002) | −0.61 |

The OP table in the OpenMULTI3D atom gives 73.5 Mb at the OP threshold. Ours, extrapolated to the threshold, is 81.7 Mb.

## T1: the 3p² ¹S₀ autoionizing resonance

The Fano fit uses 74 pooled points within ±2,000 cm⁻¹ (rms 16 %):

| | E_r (cm⁻¹) | Γ (cm⁻¹) | q | Peak (Mb) |
|---|---|---|---|---|
| Model | 68,124 ± 3 (threshold-aligned); 67,300 model-absolute | 608 ± 7 | 21.8 ± 1.9 | 185 |
| Lab | 68,270 (5) | 280 (5) | 24.5 | 785 ± 141 |

- The model resonance sits in the right place and has about the right shape (q), but it is **2.2× too broad**.
- Its integrated strength is about 0.5–0.8 of the lab value, so the peak is 4× too low.
- 3p² ¹S therefore couples too strongly to the 3s εs continuum in our CI.

Post-hoc diagnostic (labelled as such; not a declared test): the same 8 boxes, with only the ¹S and ¹P° symmetries, were rerun with the dielectronic polarisation term scaled to 0 (`posthoc_diag_dielectronic/`).

| Dielectronic scale | E_r (cm⁻¹) | Γ (cm⁻¹) | q | Peak (Mb) |
|---|---|---|---|---|
| 1 (as declared) | 68,124 | 608 | 21.8 | 185 |
| 0 | 68,463 | 418 | 29.2 | 269 |

- The baseline reproduces the full run exactly.
- Removing the dielectronic term narrows the resonance by about a third, but it remains 1.5× the lab width, and E_r moves 190 cm⁻¹ the wrong way.
- So the term over-couples 3p² ¹S to the continuum but is not the whole story. The remaining excess points to missing 3p² correlation: the compact basis gives 3p² too little weight relative to 3s εs.
- This diagnostic is not adopted. A v70d declaration would test it with a fitted r_c for the k = 1 dielectronic term on 3p² and a larger compact set.

## T4: static polarizability sum rule (FAIL as declared)

| | 3s² | 3s3p ³P |
|---|---|---|
| v70 basis (complete pseudo-spectrum) | 74.47 | 106.80 |
| Enlarged basis, R = 400 (declared) | **74.78 (+4.8 %)** | **107.41 (+6.1 %)** |
| Velocity gauge | 77.02 | 108.26 |
| NIST energies for matched terms (**post hoc**) | 72.90 (+2.2 %) | 104.90 (+3.7 %) |
| Reference (Mitroy+2010 Table XI) | 71.33 (RCI+MBPT); sum-rule experiment 74.9 (2.7) | 101.2 (3) |

**Finding.** The model's 3s3p ¹P excitation energy is **34,149 cm⁻¹ against NIST 35,051 (−2.6 %)**.
- α goes as 1/ΔE², so the resonance-line term alone is inflated by about 5 %. Its f is 1.724 (length) and 1.779 (velocity).
- The v70 energy gate, 0.45 % of the total two-electron energy, hides this excitation-energy error.
- **Line gf are unaffected**, because gf use NIST ΔE.
- Anything summed over states with model energies is affected: polarizabilities, C₆ and vdW from sum rules.

The experimental sum-rule value 74.9(2.7) happens to match ours, but the declared referee was the theory consensus 71.33.

## T5: bound–free continuity at the limit

- Below the limit: Rydberg members n* 7–12 (the v70 validity limit BOX_N), f·n*³, extrapolated to the threshold.
- Above the limit: box points at 440–2,630 cm⁻¹, extrapolated to the threshold.

**Harness corrections made after the first dry run (2 boxes), before the full read:**
1. States beyond n* = 12 are box-compressed, and one of them collects all the missing high-n strength, so the bound window stops at BOX_N.
2. The threshold-adjacent box states are skipped.
3. Initial terms inside the Rydberg window (n* > 6) are excluded, because they have no series sample.

**Result.** 16 of 61 channels are within 15 %, so T5 fails as declared. The structure is clear:
- **Low initial levels (n ≤ 4), strong channels: pass, mostly within 1–3 %.**

  | Initial term | Final symmetry | Continuum / Rydberg ratio |
  |---|---|---|
  | 3s² | ¹P | 0.963 |
  | 3s3p ³P | ³D | 1.018 |
  | 3s3p ¹P | ¹D | 0.990 |
  | 3s4s ³S | ³P | 1.032 |
  | 3s4s ¹S | ¹P | 1.023 |
  | 3s3d ¹D | ¹P | 1.010 |
  | 3s3d ¹D | ¹F | 1.005 |
  | 3s3d ³D | ³P | 1.005 |
  | 3s3d ³D | ³F | 1.025 |
  | 3s4p ³P | ³D | 0.994 |
  | 3s4p ¹P | ¹D | 1.009 |

  The model is continuous across the limit: the Rydberg series and the box continuum describe the same physics.
- **Fails, weak channels at a cancellation or on the 3p² ¹S wing:** 3s3p ³P → ³S (2.2) and 3s3p ¹P → ¹S (12.8).
- **Fails, initial levels n ≥ 5:** the extrapolation is dominated by the initial level's near neighbours, and some values go negative. This weakens the test design and says nothing about the model.

## Report: against the Opacity Project (OP) tables in `atom.mg86g_nq_qmax`

The ratio is ours over OP at the threshold.

| Level | Ours σ_thr (Mb) | OP (Mb) | Ratio |
|---|---|---|---|
| 3s² ¹S | 1.47 | 2.93 | 0.50 |
| 3s3p ³P | 16.2 | 16.2 | 1.00 |
| 3s3p ¹P | 81.7 | 73.5 | 1.11 |
| 3s4s ³S | 0.93 | 0.57 | 1.64 |
| 3s4s ¹S | 0.49 | 0.20 | 2.48 |
| 3s3d ¹D | 29.6 | 34.8 | 0.85 |
| 3s4p ³P | 22.5 | 27.4 | 0.82 |
| 3s3d ³D | 20.9 | 25.9 | 0.81 |
| 3s4p ¹P | 18.2 | 17.8 | 1.02 |
| 3s5s ³S | 2.55 | 2.36 | 1.08 |
| 3s4d ¹D | 34.3 | 37.4 | 0.92 |
| 3s4d ³D | 38.1 | 52.2 | 0.73 |
| 3s5p ³P | 33.7 | 46.4 | 0.73 |
| 3s4f ¹F / ³F | 19.4 / 19.2 | 35.6 / 35.0 | 0.55 |

- **Strong channels agree with OP within ~20 %.**
- **The s-states that differ** (3s², 3s4s, 3s5s ¹S) have near-threshold σ that sits at a cancellation. 3s5s ¹S even extrapolates to slightly negative.
- **Ground state: ours 1.47 Mb, against OP 2.93.** Lab threshold values of about 2.1 ± 0.3 Mb (Fung et al.) and 2.36 Mb (hook method) are known only from a search snippet; they are not held and not verified, so they are not used.
  - T5 shows our value is continuous with our own 3s np ¹P series. Those lines are class E cancellation lines (CF 0.01–0.07), so a small cancellation error moves both.
  - No Mg I lab or MCHF referee for the 3s np ¹P series is held. Next step: fetch the NIST MCHF Mg I collection.
- **Our 4f is about 0.55× OP**, and OP's sparse table makes shape comparison weak. The median log ratio over OP points in range is +0.08 dex.

## Open for Ryan / next

1. **T1 width (×2.2).** This is the main physics failure. See the diagnostic above. The candidates are the strength of the dielectronic polarisation term in the 3p² ↔ 3s εs coupling, and missing 3p² correlation. A declared v70d would test them.
2. **Excitation-energy error (−2.6 % on 3s3p ¹P)**, hidden by the total-energy gate. Proposal: add an excitation-energy gate for every downstream sum-over-states quantity, or adopt NIST energies there (the hybrid practice).
3. **Ground-state threshold** σ is 0.5 × OP. It needs the NIST MCHF Mg I 3s np ¹P f-values as referee. Not held; fetch them.
4. T5 is valid only for low initial levels. A redesigned T5 (a quantum-defect fit of the series) would test high-n.
5. Nothing here is adopted into NLTE. OP remains the production source.

## Sources

- Held: Baig 2022 (*Atoms* 10, 39), Mitroy, Safronova & Clark 2010, Bergemann 2017, Peralta 2022, Osorio 2015, Lind & Amarsi 2024.
- NOT HELD (paywall; rows marked NOT HELD, values taken via Baig 2022): Rafiq 2007, Madsen & Thomsen 2002, Bradley 1976, Bonanno 1986.
- NOT HELD: Butler, Mendoza & Zeippen 1993. TOPbase is down, so the OM atom carries the OP data.
