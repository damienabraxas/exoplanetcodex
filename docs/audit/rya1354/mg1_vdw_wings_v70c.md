# RYA-1354 v70c: solar wing test of van der Waals broadening (Mg I)

Status: **REPORT ONLY / HOLD.** Declared on RYA-1354 (comment bdf5d990, 2026-11-21) before any run. Attempt 1 is FAILED_RETAINED (`data/audit/rya1354/mg1_vdw_wings_v70c/attempt1/`). The attempt-2 mask is a POST-HOC correction, declared in reply c3ecd398 before the rerun. No synthesis input or published value changed.

## Method
- gf is fixed to lab values: Pehlivan 2017 exp, NIST AA/A, Smith & Raggett.
- A(X) is fixed to AAG21: Mg 7.55, Na 6.22, Ca 6.30.
- The free parameters are d = Δlog γ₆ relative to ABO (the ABO σ is scaled and α kept) and one continuum scalar per window.
- Turbospectrum 1D LTE through iSpec: MARCS is primary, ATLAS9 is the atmosphere term.
- Atlases: IAG 2020, KP1984 (4703), CRIRES+ H, with telluric pixels at min_mtrans < 0.9 dropped.
- The literature method is Anstee & O'Mara 1995 and Bergemann 2017.

Code:
- `scripts/rya1354_vdw_wings.py`
- `scripts/rya1354_vdw_wings_analyze.py`

Outputs: `summary.csv` and `summary.json`.

Attempt-2 mask rule, applied to every window: own depth ≥ 3 × blend depth, and at least 30 pixels.

## Result (d in dex; σ = curvature ⊕ A ⊕ gf ⊕ atmosphere)

| Line | d (MARCS) | σ | z raw | z vs controls | Unsöld Δχ² |
|---|---|---|---|---|---|
| Na D (ctrl) | −0.031 | 0.037 | −0.85 | | |
| Ca 6122 (ctrl) | −0.033 | 0.048 | | | |
| Ca 6162 (ctrl) | +0.004 | 0.045 | | | |
| Ca 6102 (ctrl) | −0.188 | 0.075 | | | (blend, Fe I 6102.17) |
| Mg 5172 | −0.111 | 0.042 | −2.66 | −0.78 | 974 |
| Mg 5183 | −0.064 | 0.044 | −1.46 | −0.28 | 4438 |
| Mg 5528 | −0.114 | 0.032 | −3.56 | −0.84 | 2122 |
| Mg 8806 | −0.010 | 0.044 | −0.23 | +0.28 | 4263 |
| Mg 15024 | −0.081 | 0.034 | −2.34 | −0.47 | 2581 |
| Mg 15040 | −0.125 | 0.042 | −3.02 | −0.93 | 523 |
| Mg 15047 | −0.019 | 0.041 | −0.47 | +0.19 | 1903 |
| Mg 15740 | −0.049 | 0.062 | −0.79 | −0.12 | 398 |
| Mg 15748 | −0.193 | 0.052 | −3.74 | −1.56 | 54 |
| Mg 15765 | −0.226 | 0.059 | −3.84 | −1.82 | 0 (no power) |
| Mg 4703 (T3, list ABO) | −0.311 | 0.051 | −6.09 | −2.75 | |
| Mg 5167 | UNMEASURABLE (Fe I 5167.49 at 0.17 Å) | | | | |

Δ_ctrl = −0.037. The control scatter is 0.086, inflated by Ca 6102.

## Verdict
- **Harness gate: FAIL.**
  - Ca 6102 vs its sisters: z −2.3 / −2.9.
  - Two marginal H-band pairs: z −2.19 and +2.28.
  - The b-triplet and 4s–4p (15024) pairs otherwise pass.
- **T1/T4 by the declared AND rule:**
  - 4 of 12 pass raw (5183, 8806, 15047, 15740).
  - All 12 pass control-differenced (|z| ≤ 1.82).
- **Finding:**
  - Every Mg line wants d ≤ 0. The median is −0.10 dex, against −0.03 for Na D and the clean Ca lines.
  - So in 1D LTE, our ABO for Mg I is about 0.07 dex broader than the Sun wants, beyond the control offset.
  - This is within the controls' spread, but systematic in sign.
- **T2:** Unsöld (our ⟨r²⟩) is rejected wherever there is power (Δχ² 400–4800). ABO with n*_r is indistinguishable from ABO (|Δ| ≤ 0.02).
- **T3, 4703:** the Sun wants the list ABO (fdamp 2806.269) narrower by 0.31 ± 0.05.
  - This is independently consistent with Bergemann 2017's low 1D 4703 abundance (7.44 against a mean of 7.51).
  - The upper level 5d ¹D has n* outside the ABO tables.
- **T4, H-band fill:**
  - Our ABO fill (−7.148 / −6.923) is only 0.05–0.08 dex above the list's implicit Unsöld × 2.5 (Turbospectrum fdamp 2.5).
  - So the fill changes little for these lines. The Sun slightly prefers the narrower value.
- **Untested:** the 3p² ³P lines at 2776–2783 Å (no atlas).

## Not done
- No adoption.
- Ca 6102 and the marginal IR pairs are not tuned away.
- 1D LTE is the stated caveat: a 3D/NLTE wing test is outside v70c.
