# Amarsi et al. 2021 Solar molecular CNO reference set (RYA-1136)

This holding pins CDS catalogue `J/A+A/656/A113`, Table 2: the 408 molecular
lines actually used in the homogeneous 3D LTE Solar CNO analysis inherited by
the AGSS21 abundance lineage.

Run `python3 scripts/ingest_amarsi2021_cno_rya1136.py` to regenerate the
normalized CSV and manifest. The source table supplies molecule, electronic
system, vibrational band, vacuum wavelength, lower energy, log(gf), equivalent
width, and model-by-model abundance. It does **not** provide full rotational
quantum identities, so these rows remain `CROSSMATCH_REVIEW` until joined to
the exact upstream molecular releases on physical identity.

The 408 used rows comprise 80 `12C16O`, 39 `C2`, 54 `CH`, 59 `CN`, 31 `NH`,
and 145 `OH` transitions. Whole synthesis distributions must not be described
as “AGSS21-used” merely because they contain nearby wavelengths.
