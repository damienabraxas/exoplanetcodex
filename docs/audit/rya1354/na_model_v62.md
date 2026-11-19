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

## v63: the same model for Mg II, and the σ rule generalised

The model is now one generic module for one-valence-electron spectra, `scripts/rya1354_onevalence_model.py`, configured per species. Na I reproduces v62 to 0.0007 dex, with no class changes.

**Mg II inputs.**
- Core polarizability α_d(Mg²⁺) = 0.4814.
- Potential re-fitted to the NIST Mg II terms (all levels within 3.4×10⁻⁵).
- **Anchor:** KP08 Table 44, 3s–3p A = 2.60e8 / 2.57e8 s⁻¹. This is class A, combining Froese Fischer, Siegel, Johnson and Ansbacher's measured lifetimes.
- **Referees:** Froese Fischer [9] and [16] for Mg II, from the NIST MCHF collection (`data/reference/onevalence_nist/mchf_mg2/`).
- This milestone was not declared on the ticket before it ran; it is recorded here and on RYA-1354 afterwards.

**Score.**
- h & k: −0.011 / −0.010 dex against the anchor.
- Asplund's Mg II lines (4481, 7877/7896, 9218/9244, 10915/10952 Å): within 0.015 dex of FF.
- **Referee spread:** FF [9] and [16] agree to ~0.01 dex for Mg II, as for Na, except at the Cooper-minimum 3s–4p line, where they differ by 0.25 dex.

**What changed in the σ rule.** Mg II's series involving d states carry 0.03–0.06 dex errors relative to FF; Mg⁺ 3d penetrates the core. The half-operator indicator does not see these, so σ now has more terms:
- σ = √(floor² + (½ operator effect)² + family²).
- **Family term:** the 68th percentile of |model − FF| for the line's transition family (s–p, p–d, d–f, …). It is measured on n ≤ 7 lines with CF ≥ 0.1, and our model is not shifted toward FF, since FF is theory, not lab.
- **Strongly cancelling lines (CF < 0.1)** add their CF-bin term instead, so a Cooper-minimum line cannot inflate its family's σ.
- **Floor:** the larger of the model's anchor offset and the anchor's own uncertainty, rounded up to 0.001. That gives 0.006 for Na and 0.013 for Mg II.

**Coverage** is measured on the same lines used to calibrate σ, so it is in-sample:

| | Within 1σ | Within 2σ |
|---|---|---|
| Na I | 92% | 100% |
| Mg II | 78% | 95% |

**Mg II grades.**

| Lines | Class | σ (dex) |
|---|---|---|
| h & k | B | 0.024 |
| 4481 Å | B | 0.017 |
| 9218/9244 Å | B | 0.019 |
| 3p–3d (2791/2798 Å) | C | 0.052 |
| 7877/7896 Å (4p–4d) | C | 0.05 |

Mg II has no class A lines.

**Fill:** all 63 Mg II features in the solar line list (`data/audit/rya1354/mg2_model_v63/solar_fill.csv`). Only one measurable feature lacks a FF value, 2202.7 Å (3d–10f); for Mg II, FF already covers nearly every solar line, so the value here is the grade rather than new lines.

**Next:**
- K I, once NIST's line server is back (it returned errors on 2026-11-19).
- Ca II: referees are Hettrich 2015 and the Safronova Ca II paper in the local reference folder.

## v64: Ca II, where the one-electron model reaches its limit

Declared on RYA-1354 before the run.

**Inputs.**
- Core 3p⁶, α_d(Ca²⁺) = 3.26(3) (Safronova & Safronova 2011, Sect. VI).
- Lowest valence n per l taken from the NIST levels: 4s, 4p, 3d, 4f. The fitted potential reproduces the NIST terms to 4×10⁻⁵.

**Referees.**
- **Experiment:** Hettrich et al. 2015, a single trapped ion. The H line (4s–4p₁/₂) has A = 2π × 21.57(8) MHz; the 8662 Å IR-triplet line (3d₃/₂–4p₁/₂) has A = 2π × 1.482(8) MHz.
- **Theory:** Safronova 2011 Table V (all-order relativistic; 55 lines, each with an uncertainty of 0.1–2.6%), in place of a Froese Fischer collection, which does not exist for Ca II.

**Score.**

| Anchor | Our model | Without polarization operator | Safronova |
|---|---|---|---|
| H line | **−0.051 dex** | | |
| 8662 Å | **+0.078 dex** | +0.033 | +0.006 |

The experimental floor is therefore 0.065 dex. **No Ca II line reaches class A or B.**
- Against Safronova, our class C lines agree to a robust 0.013 dex.
- The IR triplet (3d–4p) sits +0.072 above Safronova and is class D (σ 0.10). H & K are class C (σ 0.083).
- Asplund's Ca II lines (5001, 8249, 8255 Å) are within 0.012 dex of Safronova.

**Interpretation.** Ca⁺ 3d is a compact orbital that penetrates the core, and the 4s–4p rate is sensitive to core-valence correlation that a single polarization term only approximates. This is where a multi-electron calculation (GRASP, from the O pilot) would have to take over. The grade reports the limit instead of hiding it.

**Fill:** all 90 Ca II features in the solar line list (`data/audit/rya1354/ca2_model_v64/solar_fill.csv`); 26 measurable features have no Safronova value.

| Features | Class | σ (dex) | vs line list |
|---|---|---|---|
| **H band 6p–6d:** 16561.1, 16649.9, 16662.7 Å (depths 0.10, 0.13, 0.03) | C | 0.070 | −0.016 |
| 5p–nd and 4d–nf series, 2898–4110 Å | C/D | | |
| UV 3d–5p (2129–2132 Å) | E | | +0.68; near-cancelling line |

**Sparse-referee rule:** with only 55 referee lines, empty CF bins fall back to the other n regime, then to the overall referee scatter (flagged `fallback`).
