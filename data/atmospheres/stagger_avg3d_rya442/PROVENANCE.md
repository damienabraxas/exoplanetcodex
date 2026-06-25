# ⟨3D⟩ STAGGER solar model — RYA-442

`sun_avg3d_stagger.mod` — the spatially + temporally averaged ⟨3D⟩ solar model from
the STAGGER grid (Magic et al. 2013, A&A 557, A26), averaged on surfaces of constant
optical depth at 500 nm, in Turbospectrum-v20 averaged format (header `TAU5000 SCALE`,
NDEP=101; columns: log τ₅₀₀₀, T[K], nₑ[cm⁻³], V, depth-dependent v_mic[km/s]).

- Parameters: Teff 5777, log g 4.44, [Fe/H] 0.00 (member `p5777_g+4.4_m0.0_t02_st_z+0.00...`).
- Conversion to the TS-v20 averaged format: Ekaterina Semenova, July 2020.
- Source archive: `average_stagger_grid_forTSv20.zip` from the Turbospectrum/NLTE
  data collection on MPG Keeper (TSFitPy reference set,
  https://keeper.mpdl.mpg.de/d/6eaecbf95b88448f98a4/), `/atmospheres/`.
  Original ⟨3D⟩ stratifications: staggergrid.wordpress.com/mean-3d/.
- Retrieved 2026-06-24 for the RYA-442 disk-center 1D→3D CO probe.

INGESTION CAVEAT (RYA-442): this averaged 5-column `TAU5000 SCALE` format does NOT
match either of iSpec's two babsma model-input modes cleanly (MARCS-FILE .true. =
full MARCS columns; .false. = Kurucz RHOX columns). The probe fed it via
`atmosphere_layers_file` with MARCS-FILE .false.; the ingestion is NOT validated and
the ⟨3D⟩ probe result is treated as INDICATIVE ONLY, not a measurement.
