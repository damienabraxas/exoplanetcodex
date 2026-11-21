# RYA-1354 v67: Mg I with spin–orbit (intermediate coupling) and a measured cancellation term

Status: **CANDIDATE / HOLD** (RYA-1352). Declared on RYA-1354 (2026-11-20 05:41 UTC) before the run. This resolves the two open calls from v65 (`mg1_model_v65.md`).

Script: `scripts/rya1354_twovalence_model.py 'Mg I'` (about 4 min; `--cache=PATH` is for development only). Outputs are in `data/audit/rya1354/mg1_model_v67/`. The v65 outputs are kept unchanged as the pure-LS record.

## What changed

- **Spin–orbit.** One-body Breit–Pauli H_so = Σ ξ(rᵢ) lᵢ·sᵢ, with ξ_l = (α²/2)(1/r) dV_l/dr from our own Mg²⁺ potential, **unscaled**.
  - For each (J, parity), the basis is the NIST-matched LS CI states.
  - Diagonal: the NIST term centroid. Off-diagonal: the Racah spin–orbit matrix elements.
  - States are assigned to NIST J levels by largest component (Hungarian assignment).
  - gf uses NIST J energies.
- **Control.** With ξ = 0 the IC code reproduces every v65 LS gf exactly (max deviation 0.0, none missing). This validates the new phases.
- **Cancellation term.** CF < 0.1 lines now add the abs68 of model − line list in their CF bin (n ≤ 8, measurable, spin-allowed). That replaces the pooled lab term. CF < 0.02 gets class E and no σ.
- **Forbidden lines.** 2154 Å (E2) and 2089 Å (M2) are labelled non-E1 (KP08 Table 42).

## Declared tests

| Test | Result |
|---|---|
| (a) Fine structure | p triplets 77–83 % of NIST (3s3p ³P 0.769; Mg II 3p 0.85); 4f/5f 1.00/0.92. nd ³D is inverted in NIST (−0.08 cm⁻¹) and the one-body ξ cannot produce that. |
| (b) Lab spin-changing lines (Pehlivan) | 14878 Å **−0.078**, 12083 Å **−0.071** dex. v65 had no value. |
| (c) 4571 Å vs Tachiev & Froese Fischer BP-MCHF (−5.623, theory D) | −5.824 (−0.20). With ξ scaled to 3s3p ³P (×1.30): −5.596 (+0.03). |
| (d) 24 spin-allowed lab lines | median +0.020 (unchanged); abs68 0.031 → 0.025 |

The ξ-scaled variant is the declared σ indicator. It improves 4571 Å but worsens the lab 4f lines (+0.15/+0.16). A single global scale is wrong: p is 23 % short and f is already right. A per-l scale is the obvious refinement, and it is a method change for Ryan.

## Measured CF bins (spin-allowed)

| CF | n | abs68 (dex) |
|---|---|---|
| 0.05–0.1 | 37 | 0.154 |
| 0.02–0.05 | 7 | 0.331 |
| < 0.02 | | no σ, class E |

## Finding: VALD's Mg I intercombination gf are refuted by the observed Sun

v67 now grades spin-changing lines. Against the line list they look much worse: abs68 0.143 overall vs 0.092 on spin-allowed lines. VALD gives several of them near-allowed strengths. I measured them in the observed atlases (`rya1354_candidate_spectrum_check.load_atlases/measure`), with allowed Mg I lines in the same bands as positive controls:

| Line | VALD log gf | Pred. depth | Observed | Atlas |
|---|---|---|---|---|
| CTRL 15025 4s³S–4p³P₂ | +0.36 | 0.56 | 0.566 | CRIRES+ H |
| CTRL 17109 4s¹S–4p¹P | +0.06 | 0.54 | 0.514 | CRIRES+ H |
| CTRL 11828 3p¹P–4s¹S | −0.33 | 0.53 | 0.600 | CRIRES+ J |
| CTRL 10811 3d³D–5f³F | +0.02 | 0.38 | 0.520 | IAG |
| 20632 4p¹P–4d³D | +0.13 | 0.42 | **0.008** | CRIRES+ K |
| 11812 4p³P₁–5d¹D₂ | −1.04 | 0.24 | **0.027** | CRIRES+ J |
| 9822 4p³P₁–6d¹D₂ | −1.41 | 0.15 | **0.004** | CRIRES+ Y |
| 20867 5p³P₁–7d¹D₂ | −1.36 | 0.13 | **0.003** | CRIRES+ K |
| 7746 3p¹P–3d³D | −2.14 | 0.36 | **0.010** | IAG |

The lines are absent, as our model predicts (it gives them 10⁴–10⁵ times weaker). The lab agrees too: 26399 Å 4p¹P–4d¹D matches pure LS to 0.007 dex, which rules out the ~50 % 4d singlet–triplet mixing VALD's values would need.

Blue cases (4380, 4561 Å) lie in crowded Kitt Peak regions and cannot be attributed. Spin-changing lines are excluded from the CF-bin measurement for this reason.

**Line-list problem (noted, not fixed):** these VALD rows inject absorption that is not in the Sun.

## Grade and fill

- Line classes: C 1,404, D 773, E 2,040. Many E are weak intercombination lines with low CF.
- Lab coverage (in-sample, spin-allowed): 83 % within 1σ, 96 % within 2σ.

| Fill status | Features |
|---|---|
| graded | 438 (389 measurable, v65: 335) |
| beyond box | 658 |
| ungraded_cancellation | 31 |
| non-E1 | 2 |
| unassigned | 9 |
| not in CI | 6 |
| spin-change, no model | **0** (v65: 145) |

Model − line list, spin-allowed measurable: median −0.001, abs68 0.092 dex.

- 4571 Å: −5.824, σ 0.127, class D.

## Open for Ryan

1. **Per-l ξ calibration** (scale the p spin–orbit to the NIST p-triplet splittings, leave f): would move 4571 Å to about −5.60, in line with Tachiev & Froese Fischer. Currently declared unscaled.
2. **VALD Mg I intercombination rows** in `linelist_solar.csv` are refuted in the IR. This is a line-list issue for another ticket.

Next: v66 Ca I (Den Hartog 2021), with this IC machinery from the start.
