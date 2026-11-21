# RYA-1354 v65: our own two-valence-electron model, Mg I

Status: **CANDIDATE / HOLD**. Nothing here writes a gf, grade or σ into a product (RYA-1352). Declared on RYA-1354 (2026-11-19 19:39 UTC) before the run.

Script: `scripts/rya1354_twovalence_model.py 'Mg I'` (Mac, about 4 min). Outputs are in `data/audit/rya1354/mg1_model_v65/`: `summary.json`, `lines.csv` (2,288 fine-structure lines) and `solar_fill.csv`.

## Model

Two valence electrons outside the Mg²⁺ core. The core is the one-valence potential already fitted to Mg II in v63.

- **CI basis:** full CI in LS coupling, l ≤ 4, about 35 orbitals per l (ion orbitals, screened Rydberg states and compact correlation functions).
- **Electron–electron repulsion:** 1/r₁₂ computed exactly from Slater integrals.
- **Core polarization:** a dielectronic term −α_d (r₁·r₂)/(r₁³r₂³)·cutoff, plus the polarization-corrected dipole operator.
- **No spin–orbit.**
- **gf:** fine-structure components use 6j fractions on NIST energies.

## Pre-registered energy gate: PASS

Every NIST term the box can hold (outer n ≤ 12) is reproduced within **0.45 %** of the total two-electron energy. The worst terms are 3s² ¹S (822 cm⁻¹ overbound) and 3s3p ³P (693 cm⁻¹).

Process note: the first run aligned CI states to *all* NIST terms. It paired box-continuum states with the n = 58/59 p and 24/25 d Rydberg terms, which gave 3.3 % and a FAIL. The box limit (`BOX_N = 12`) was already in the script, but the aligner did not apply it. Terms with n > 12 (70) are now excluded rather than force-matched.

Unmatched terms: 4s4p and 4s5–7p ¹P. These are autoionizing doubly-excited states.

## Lab score (Pehlivan Rhodin 2017, 26 lines)

| Set | n | Model − lab | Pehlivan's ATSP2K − lab |
|---|---|---|---|
| Independent (lab lifetime) | 6 | median 0.000, mean abs 0.046 | |
| All spin-allowed | 24 | median +0.020, abs68 0.031 | median +0.006, abs68 0.024 |

**Mg I b triplet, reported and not tuned** (lab lifetime is theoretical): 5167 +0.023, 5173 +0.008, 5184 +0.034 dex.

The largest misses are:
- 5711 Å (3p ¹P–5s ¹S): +0.139. ATSP2K is also +0.10 here.
- 31157 Å: −0.082.

These two set the floor at **0.050 dex**, so no Mg I line is class A or B.

The two spin-changing lab lines (14878, 12083 Å) are flagged "no model" (no spin–orbit). Their lab gf, 5–6 % of the allowed lines, is direct evidence that 4f is LS-mixed.

## Grade (declared rule)

σ = hypot(floor, ½ polarization-operator effect, family abs68 vs lab). Families with fewer than 5 lab lines take the pooled 0.031 dex. Lines with CF < 0.1 also add the pooled term. The class comes from Li 2023 Table 5.

- Coverage (in-sample, 24 lab lines): 83 % within 1σ, 96 % within 2σ.
- Lines by class: C 1,206, D 985, E 97.

## Solar line list fill (1,697 Mg I rows → 1,144 features)

| Status | Features |
|---|---|
| graded | 366 (335 measurable) |
| no_model_beyond_box (n > 12, the near-limit 3p ³P–nd ³D forest at 2529–2546 Å) | 616 |
| no_model_spin_change (4571 intercombination line, etc.) | 145 |
| unassigned | 11 |
| no_model_not_in_ci | 6 |

Two grouping rules apply:
- **Degenerate NIST J levels** (e.g. 7f ³F J = 2,3,4) are compared at the term level, as summed gf, never by guessing J.
- **LS-broken configurations**, where the singlet–triplet separation is ≤ 3× the triplet fine-structure width or < 1 cm⁻¹ (all nf and ng), are compared at the configuration level: gf summed over S and J, which is conserved under mixing.

**Model − line list** (a diagnostic only; Mg I VALD gf is itself theory): median 0.000, abs68 0.089 dex on 335 measurable features. 31 features differ by > 0.3 dex.

## Known limits, open for Ryan

1. **Strong cancellation is under-graded.** Example: 4s ¹S–7/8/9p ¹P (CF 0.003–0.005) are 1–2 dex weaker than VALD but class D. The declared rule demotes CF < 0.1 by only one class. Proposal: CF < 0.02 → class E or no grade.
2. **Spin-allowed J components at high n.** For 3p/4p/5p ³P → nd ³D at n ≥ 9 the three components show a fixed −0.08 / −0.56 / −2.06 dex pattern against VALD. This is consistent with nd ¹D–³D mixing, which this model lacks (no spin–orbit), or with a change of VALD's source. No referee can tell them apart. Proposal: grade those at the multiplet level, or add a Breit–Pauli spin–orbit term.
3. **Mg I 4571** (intercombination) needs spin–orbit and is not covered.

Next per the declaration: v66 Ca I (referee Den Hartog 2021, 78 lab lines). Ca I has 3d in the valence space, so the floor is expected to be worse (cf. Ca II, 0.065).
