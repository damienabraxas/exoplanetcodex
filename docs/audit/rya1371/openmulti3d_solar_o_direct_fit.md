# Direct OpenMULTI3D fit of the solar O I 777 nm triplet (RYA-1371, Phase 0)

2026-11-18. This follows `openmulti3d_solar_o_nlte.md`. That file applied OpenMULTI3D NLTE
corrections to the ORCH's existing LTE fits. Here the observed spectra are fitted directly
with OpenMULTI3D spectra, in LTE and NLTE.

## Method

- **Synthesis:** OpenMULTI3D 1D on MARCS 5777/4.4 (ξ 1 km/s). Atoms `o41f` (Barklem +
  Kaulakys H collisions) and `o41j` (Belyaev); `o41f_ref` gives the same triplet as `o41f`.
  - Grid: A(O) 8.55–9.05 in steps of 0.05.
  - Window: 7700–7845 Å, including the engine's Bergemann 2026 background line list.
  - The triplet itself uses the atom's log gf 0.369 / 0.223 / 0.002, the same as the ORCH.
  - Spec: `config/experiments/rya1371_openmulti3d_solar_o_wide.json`.
- **Observations:** the IAG atlas and the Kitt Peak atlas in both telluric corrections
  (Kurucz 2005, molecfit). All are read through the ORCH loader
  (`measure_band_ew.load_window_ex`), unchanged.
- **Profile fit:** χ² over ±1.1 Å per line, the ORCH's half-width for this band. The free
  parameters are:
  - A(O), linear in the model flux between grid points;
  - one Gaussian broadening (macroturbulence + instrument, km/s);
  - a wavelength shift.
- **EW inversion:** the observed EW within ±0.6 Å of the line, minus the synthetic blend EW
  (background-line run minus O-only run). This is inverted on the O-only curve of growth,
  measured in the same window with the fitted broadening.
- **Continuum, three legs:**
  1. `synth_ref`: the RYA-1232 nominal rule (`fit_model_guided` with the neighbourhood
     check), with the OpenMULTI3D synthesis as reference.
  2. `iag_ref`: the RYA-587 `contref` leg, with the IAG atlas as reference.
  3. `envelope`: the local upper envelope (`local_continuum.fit`, the Amarsi+2020
     local-maxima analogue).

Script: `scripts/fit_openmulti3d_solar_lines.py`. Outputs:
`data/results/rya1371/openmulti3d_solar_o_fit_v1/` (`fits.csv`, `summary.csv`,
`profiles.png`).

## Result

Envelope leg, means of the three triplet lines and all three holdings (holdings agree within
0.003 dex):

| Treatment | Atom | A(O) profile | A(O) EW | Line-to-line sd | Fit rms | Broadening FWHM |
|---|---|---|---|---|---|---|
| 1D LTE | either | 8.935 | 8.911 | 0.023 | 0.0023 | 3.8 km/s |
| 1D NLTE | `o41f` (Barklem+Kaulakys) | **8.786** | 8.776 | **0.003** | 0.0016 | 4.7 km/s |
| 1D NLTE | `o41j` (Belyaev) | **8.717** | 8.711 | 0.006 | 0.0015 | 5.0 km/s |

**Asplund σ for the `o41f` NLTE value** has three terms:

| Term | Value |
|---|---|
| Statistical (SEM of 3 lines) | 0.002 |
| NLTE, ½ × 0.149 | 0.074 |
| 3D, ½ × the Engine A 3D−1D difference | 0.008 |

Combined: **A(O) = 8.79 ± 0.07** in 1D NLTE.

**Continuum legs** (`o41f`, 1D NLTE, mean of all holdings):

| Leg | A(O) | Fit rms | Placed level |
|---|---|---|---|
| envelope | 8.786 | 0.0016 | 0.995 |
| iag_ref | 8.854 | 0.0061 | 1.003 |
| synth_ref | 8.722 | 0.0075 | 0.987 |

## What it shows

1. **The direct fit confirms the correction route.**
   - The direct 1D LTE fit gives 8.935, the same as the ORCH's own 1D LTE triplet (8.934)
     from a different synthesis code.
   - Direct NLTE: 8.786 (`o41f`) and 8.717 (`o41j`). The correction route gave 8.803 and 8.739.
   - The 0.02 dex difference comes from fitting the line shape rather than only its strength.
2. **NLTE is needed to fit the three lines together.**
   - In LTE, the three lines disagree by 0.023 dex (sd), with the strongest line highest.
   - In NLTE they agree to 0.003 (`o41f`) and 0.006 (`o41j`).
   - NLTE also fits the profiles better (rms 0.0016 vs 0.0023).
   - This is the line-to-line test the literature uses. Both H-collision atoms pass it, so
     it does not choose between them; their abundances differ by 0.07 dex.
3. **The continuum is the largest systematic.** Around 7772 the atlases' shipped continuum
   sits about 0.7% (IAG) and 0.2% (KP) below their own upper envelope.
   - Taking the shipped level (`iag_ref`) gives +0.07 dex.
   - The nominal rule with the OpenMULTI3D synthesis as reference (`synth_ref`) places the
     continuum at 0.987, below the observed upper envelope. That synthesis's line list lacks
     weak lines the Sun has, so the rule mistakes them for continuum. The result is −0.06 dex.
   - The envelope leg has 4–5× lower residuals (0.0016 vs 0.0061 and 0.0075) and makes the lines agree with each other.
   - When OpenMULTI3D enters the ORCH, the continuum should be placed with the ORCH's own
     VALD synthesis, as for every other treatment, so this is not a defect of the rule.
4. **1D NLTE 8.72–8.79 sits just above the 3D NLTE literature** (A21 8.69; Steffen,
   Bergemann and Magg 8.75–8.77). The `o41j` 1D value is already close to A21.

## Limits

- 1D only. The 3D term is borrowed from Engine A's ½(3D−1D) and is not measured with
  OpenMULTI3D. A solar 3D OpenMULTI3D run is Phase 2.
- One atmosphere: MARCS 5777/4.4, not the ORCH's interpolated 5772/4.438.
- Broadening is a single Gaussian, not radial-tangential macroturbulence plus rotation.
  The EW route does not depend on it and agrees within 0.01 dex.
- The OpenMULTI3D background line list is not the ORCH's VALD list.
