# RYA-1354 v61/v62: our own Na I atomic model fills the gf gaps, with a grade

Status: **CANDIDATE / HOLD**. Nothing here writes a gf, grade or sigma into a product. That needs Ryan's ruling under RYA-1352. Both milestones were declared on RYA-1354 before they were run.

This is the first element where **our own atomic model** supplies gf values, graded against experiment and independent theory, for lines that have no lab or NIST value.

## Why Na

No experimental gf exists for any solar Na I line:
- Scott et al. 2015: "No experimental gf-values are available for our solar lines".
- NIST (Kelleher & Podobedova 2008, JPCRD 37, 267; now in the local reference folder) uses one calculation and nothing else for those lines: Froese Fischer's B-spline CI.

The only experiments are:
- the 3p ²P lifetimes, which give the D-line A-values to class AA (Volz 1996, Jones 1996, Oates 1996);
- 1929 relative f-values for the UV 3s–np series.

## Referees

| Kind | Source | Local file |
|---|---|---|
| Experiment | D-line A = 6.16e7 / 6.14e7 s⁻¹ (KP08 Table 4, class AA) | `scripts/rya1354_na_model.py: EXP_D` |
| Independent theory [9] | Froese Fischer 2002 non-orthogonal B-spline CI, 385 E1 lines, length and velocity | `data/reference/na_mchf_ff2002/raw/Na_11.9.U.lsj.db.gz` |
| Independent theory [16] | Froese Fischer MCHF core-valence, 73 lines | `.../Na_11.16.mchf_cv.lsj.db.gz` |
| Energies | NIST ASD Na I levels | `.../nist_asd_na1_levels.csv` |

All of these come from the NIST MCHF/MCDHF collection (nlte.nist.gov/MCHF), downloaded on 2026-11-19.

The two Froese Fischer calculations agree to a robust 0.007 dex (22 lines). Above n ≈ 8, the B-spline set stops being a referee:
- Every transition into an n = 10 level is offset by the same −0.47 dex, the signature of box-confined pseudo-states.
- At 3p–10d, VALD's independent value (−2.37) agrees with ours (−2.39), not with Froese Fischer's (−1.89).

## v61: FAC (the O I pilot code)

FAC 9ab6def, run in both gauges and rescaled to NIST energies. Two variants: v61a, valence only (n ≤ 9), and v61b, which adds 2p⁵ core-excited configurations.

| | v61a | v61b |
|---|---|---|
| D lines, Babushkin (length) | **+0.048 dex** | +0.046 |
| D lines, Coulomb (velocity) | −0.028 | −0.021 |
| Gauge split | 17% (class C) | |
| vs Froese Fischer, all lines, robust σ | 0.070 | 0.068 |

**Diagnosis:** core polarization, the valence electron distorting the Na⁺ core, is missing. FAC's spectroscopic-orbital CI cannot represent it, so adding the 2p⁵ configurations changes almost nothing. For alkalis, this is the dominant correlation effect.

## v62: model potential with core polarization (`scripts/rya1354_na_modelpot.py`)

**Method.**
- **Potential:** one valence electron in an l-dependent parametric potential (the form of Marinescu, Sadeghpour & Dalgarno 1994), plus the core-polarization term −α_d/(2r⁴)·(1 − e^{−(r/r_c)⁶}), with α_d(Na⁺) = 0.9448.
- **Fit:** the four parameters per l are re-fitted here to the NIST terms. The fitted potential then reproduces every NIST level to n = 15 to better than 3×10⁻⁵.
- **Dipole operator:** polarization-corrected, r·[1 − α_d/r³·(1 − e^{−(r/r_c)³})].
- **Fine structure:** split by LS 6j fractions on NIST energies.
- **Solver:** shift-invert sparse eigensolver on a log grid, so states are identified by node count. A plain tridiagonal solve is ill-conditioned at small r and was rejected. Tested on hydrogen.

**Score.**

| Test | Result |
|---|---|
| D lines vs experiment | **−0.005 dex** (operator off: +0.011); Froese Fischer is +0.006 |
| All lines vs FF [9], n = 379 | median −0.003, robust σ 0.025 |
| n ≤ 7 lines (FF converged), CF ≥ 0.05 | σ68 0.002–0.009 dex |
| Asplund's 5 Na lines vs FF | +0.012, +0.010, +0.007, +0.007, −0.004 |

## The grade: per-line σ from the model's own indicators

- **Radial cancellation factor:** CF = |∫u₁ d u₂| / ∫|u₁ d u₂|. It plays the same role as Li's CF for O I.
- **Half the core-polarization operator effect.** This is the model's analogue of a gauge difference. On the principal series (3s–np), Froese Fischer lands almost exactly halfway between our values with and without the operator. For 3s–4p: ours −1.802, without the operator −1.712, FF −1.754.

**Candidate σ (dex)** is built from three pieces:
- base: √(0.006² + (½·operator effect)²). The 0.006 is the experimental floor, the model's D-line offset.
- n = 8–9: add, in quadrature, the measured model-vs-FF σ68 in the line's CF bin. This is an upper bound, because FF is degrading there.
- n ≥ 10: same as n = 8–9, flagged `no_referee`.

The class comes from Li 2023 Table 5 bands applied to σ, with the CF < 0.1 demotion.

**Validation (n ≤ 7, 185 lines vs FF):** 85% lie within 1σ and 91% within 2σ.

| Class | n vs FF | Median | Robust σ | Within 0.05 dex |
|---|---|---|---|---|
| A | 148 | −0.001 | 0.002 | 100% |
| B | 94 | −0.002 | 0.054 | 56% |
| C | 64 | −0.094 | 0.135 | 14% |
| D | 35 | −0.085 | 0.130 | 31% |
| E | 38 | −0.123 | 0.078 | 26% |

## The fill (`data/audit/rya1354/na_model_v62/solar_fill.csv`)

- **Coverage:** every Na I feature in our solar line list (143 features; unresolved fine-structure and HFS components are summed into one feature) gets our log gf, σ, class and CF. That includes all 102 measurable features (depth ≥ 1%).
- **Gaps:** 26 measurable features have no Froese Fischer/NIST value.

| Features | Class | σ | vs line list |
|---|---|---|---|
| **VIS 4213–4325 Å:** 3p–9d, 3p–11d, 3p–12d, 3p–11s, 3p–12s (depths up to 0.19) | C | 0.075–0.084 | ≤ 0.03 |
| **IR:** 9153.9 (3d–9f), 10182.9 (4p–9d), 18027.9 (4d–9f, depth 0.06), 18186.0 (5p–11d, class B), 19796.9/19806.5 (5p–9d) | B/C | 0.03–0.08 | ≤ 0.01 |
| **UV 2433–2491 Å:** 3s–np principal series, n ≥ 10 | E | 0.14–0.16 | 0.25–0.52 below |

The line list's values for the IR lines come from an unrecorded VALD source; our model agrees with them to 0.01 dex or better.

On the UV principal series, our model runs 0.25–0.52 dex below VALD. That series is the most sensitive to core polarization, and the grade says so.

## What is next

1. **Observed Sun and the ORCH:** screen the VIS 4213–4325 Å and IR gap features in the solar atlases, then measure them in the ORCH against the Asplund-line control (the paired test Fe I H passed).
2. **Principal series:** score the 3s–np series against the 1929 relative f-values (KP08 source "6n"). That is the only experiment that tests the series where our σ is largest.
3. **Next elements:** the same model-potential method applies to K I and the other one-electron spectra. Mg I and Al I need a two- or three-electron model, or GRASP, which is where the O pilot's machinery comes back in.

## Provenance (sha256 of the stored reference files)

```
840071dc0e02f40bc5ba87e35aff4dd54e8c7d80d0d9ede7be1b0d51914a765f  Na_11.16.mchf_cv.lsj.db.gz
6bfb6771ef89eace3677ffa6ae551fb2fc75d2d1befdc50796bb403211161b17  Na_11.16.mchf_cv-lev.lt.db.gz
67548b5751549d56f07d555c83563f81dfa1b9cc70e663924ef6c06d32a21721  Na_11.9.U.lsj.db.gz
0b2a9b48d8e876aed567aa73c9a3c63c45b4d5ead5c2837b056411eecc3b3c39  Na_11.9.U-lev.lt.db.gz
ba5e58f55256eec58dc6aa56788c60eb318c46aa2615eefc065cb7e9d4c8ac5b  LS_Na_Z11.gz
```

FAC v61 ran on Sirius (`/mnt/codex-data/codex/rya1354/toolchains/fac-9ab6def/sfac/sfac`, 5 min); its inputs and gzipped level/transition tables are in `data/audit/rya1354/na_model_v61/native/`. v62 runs on the Mac in 45 s: `~/codex/venv/bin/python scripts/rya1354_na_modelpot.py`.
