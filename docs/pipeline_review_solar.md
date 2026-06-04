# Solar Calibration Pipeline — Code Review Notes

**Date:** 2026-05-31
**Scope:** Solar finalization sprint (RYA-96, RYA-45, RYA-101, RYA-104, RYA-105)
**Reviewer:** Mr. Code (Claude Code)
**Status:** EW measurements complete and publishable. Abundance derivation deferred (see RYA-129).

---

## What Was Reviewed

All pipeline modules touching the solar calibration path:
`spectra_normalize.py` → `lines_fit.py` → `abundances_derive.py`

Supporting files: `config/constants.py`, `data/linelists/linelist_solar.csv`

---

## Items Checked

### 1. BERV Correction
**Status: Not directly reviewable at this stage.**
BERV correction is applied upstream in ESO's pipeline before the FITS files are delivered. The HARPS S1D spectra downloaded from the ESO archive already have BERV applied. The pipeline reads normalized wavelengths directly from `solar_normalized.csv` and performs no additional Doppler correction. This is correct for solar spectra (solar RV = 0 by definition).
**Action:** When 55 Cnc A spectra are processed, verify the RV correction sign convention: we add RV to bring stellar lines to rest wavelengths.

### 2. Wavelength Scale
**Status: ✓ Air wavelengths confirmed throughout.**
HARPS delivers spectra in air wavelengths. `linelist_solar.csv` uses the column `wavelength_air_A`. All EW measurements and line identifications use air wavelengths consistently. No Å/nm confusion found.

### 3. EW Units
**Status: ✓ mÅ throughout, no unit confusion found.**
- `_integrate_profile` returns EW in mÅ (multiplies by 1000 after trapz in Å)
- `solar_ew.csv` `ew_mA` column is in mÅ
- All threshold comparisons (`ew_min_mA`, `ew_max_mA`) are in mÅ
- `abundances_derive.py` converts correctly: `ew_mA * 1e-3` → Å before taking log

### 4. HFS Handling (Ba II, Eu II)
**Status: Partial — single-profile fit used, result is total HFS EW.**
Ba II 5853.668 and Eu II 6645.127 are measured as single Voigt profiles. For Ba II (23 HFS components) and Eu II (30 HFS components), the fitted EW is the total of all components — this is the correct quantity to use with a single-profile COG analysis at HARPS resolution (components unresolved). Both lines are flagged `hfs_total_ew` in notes.
**Limitation:** When MOOG integration arrives (RYA-129), Ba II and Eu II will need special HFS synthesis treatment. Do not derive A(Ba) or A(Eu) from these EWs without HFS MOOG synthesis.

### 5. Vectorization
**Status: ✓ No nested Python loops over wavelength arrays.**
All critical hot paths in `lines_fit.py` operate on NumPy arrays:
- `_local_renorm`: array slicing + `np.polyfit`
- `_fit_profile`: `scipy.optimize.curve_fit` on array
- `_integrate_profile`: `np.trapz` on array
- `_reduce_ew` in `abundances_derive.py`: fully vectorized

### 6. Solar Reference A(Fe)
**Status: ✓ A(Fe)☉ = 7.46 everywhere.**
`SOLAR_ASPLUND2021['Fe'] = 7.46` in `constants.py`. Searched all pipeline files — no hardcoded 7.50 (old Lodders 2003 value). `constants.py` has an inline note: `# NOTE: was 7.50 in Lodders 2003 used in 2010 thesis`.

### 7. 3D Corrections
**Status: ✓ Constant defined and ready; application deferred.**
`CORRECTIONS_3D = {'O_6300_3d_dex': -0.07, 'O_6363_3d_dex': -0.07}` in `constants.py`. The correction is applied in `abundances_derive.py` within `_abundance_one`. However, since `abundances_derive.py` is a stub pending MOOG integration (RYA-129), the correction is effectively deferred. It will be applied automatically once MOOG produces A(O) values — no code change needed at that point.
**Reference:** Caffau et al. 2008, A&A 483, 591.

---

## Known Issues — Filed as RYA-129

### Fe I Excitation Slope = +0.10 dex/eV

The slope of A(Fe I) per line vs excitation potential (EP) from the linear COG analysis is +0.10 dex/eV. The expected value for the correct Teff is |slope| < 0.01 dex/eV.

**What this means:** The linear COG cannot reproduce the correct slope because it treats all lines as forming at the same atmospheric depth with the same effective temperature. In reality, low-EP lines and high-EP lines form at different depths with different temperatures. A model atmosphere code (MOOG + ATLAS9) resolves this by integrating the line formation through the full temperature-pressure structure.

**What it does NOT mean:** Teff = 5777 K is wrong. The solar Teff is established by interferometry and bolometric luminosity. The slope reflects the limitation of the method, not an error in the stellar parameters.

### Linear COG Does Not Transfer Across Elements

The C_zero calibrated from Fe I gives ±0.5–2 dex errors for other elements. Root cause: partition functions, ionization equilibria, and line formation depths vary across the periodic table and cannot be captured by a single empirical constant. MOOG integration (RYA-129) is required for accurate multi-element abundances.

---

## What Is Publishable Now (Tuesday Deliverable)

| File | Status | Notes |
|------|--------|-------|
| `results/Solar/2026-05-31/solar_ew.csv` | ✅ Publishable | 806 lines, 22 elements |
| `results/plots/solar_ew_diagnostic.png` | ✅ Publishable | Tier 1 Fe line fits |
| `results/plots/solar_oi6300_diagnostic.png` | ✅ Publishable | O I Ni-subtraction diagnostic |
| `results/plots/solar_ca6122_diagnostic.png` | ✅ Publishable | Ca I narrow-window diagnostic |
| `pipeline/lines_fit.py` | ✅ Complete | All special-case handling done |
| `abundances_derive.py` | ⏳ Stub | MOOG integration required (RYA-129) |
| Solar abundance plots | ⏳ Deferred | Requires MOOG (RYA-129) |

---

## EW Measurement Summary (Solar)

- **Total lines measured:** 806
- **Elements with EW measurements:** 22 (19 blend=False + 3 SPECIAL)
- **Elements not measurable with HARPS optical:** Co, K, N, Sc
  - K: main lines at 7665/7699 Å — outside HARPS range (3780–6910 Å)
  - N: optical lines are UV or IR only
  - Co: all available lines below `min_fit_depth = 0.008`
  - Sc: no suitable unblended lines in optical

### Special-case handling verified

| Line | Treatment | EW | Notes |
|------|-----------|----|-------|
| O I 6300.304 | Ni I COG subtraction (Allende Prieto+2001) | 4.45 mÅ | ✓ in QA gate |
| O I 6363.776 | CN flag, narrow ±0.08 Å fit window | 1.80 mÅ | ✓ in QA gate |
| Ca I 6122.217 | Wide continuum ±2.78 Å, narrow fit ±0.10 Å | 120.7 mÅ | Improved from 243 mÅ |
| Ba II 5853.668 | HFS total EW | 166.7 mÅ | Flagged hfs_total_ew |
| Eu II 6645.127 | HFS total EW, narrow ±0.35 Å window | 6.8 mÅ | Flagged hfs_total_ew |
| Li I 6707.840 | CN blend flagged | 1.8 mÅ | Flagged CN_blend_possible |

---

## Bugs Found and Fixed During Sprint

| Bug | Ticket | Fix |
|-----|--------|-----|
| O I 6300 EW = 0.0 mÅ | RYA-104 | Ni COG subtraction + wide/narrow window approach |
| O I 6363 not measured | RYA-105 | Added to SPECIAL_MEASURES, ±0.25 Å continuum / ±0.08 Å fit |
| Ca I 6122 EW = 243 mÅ (expected 140) | Pre-sprint | LINE_FIT_WINDOWS ±0.10 Å wired up |
| Si I 6155 measuring contaminated blend | Pre-sprint | blend_flag=True in linelist |
| Ni I 6300.336 COG sanity cap wrong fallback | RYA-104 | ni6300_ew_lit_mA: 0.55 → 1.0 (Allende Prieto+2001) |
| O I continuum anchor hitting Sc II HFS | RYA-104 | LINE_WINDOWS ±0.30 Å gaps between Si I and Sc II |
| O I 6363 continuum anchor hitting deep photospheric lines | RYA-105 | LINE_WINDOWS ±0.25 Å |

---

## New Linear Issues Filed

- **RYA-129** — MOOG/ATLAS9 integration required; Fe I excitation slope 0.10 dex/eV documented
