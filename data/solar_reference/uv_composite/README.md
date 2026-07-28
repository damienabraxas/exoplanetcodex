# CALSPEC solar reference composite — UV anchor (RYA-459, under RYA-162)

**Provenance: CITED-COMPOSITE — NOT a measurement.** Hubble cannot observe the Sun
directly, so there is **no direct solar UV spectrum**. This is a literature
composite + model, used as the cited UV reference and an absolute-flux cross-check.
Per the RYA-455 discipline, every value carries `provenance=cited-composite` and must
never be presented as a direct solar measurement; downstream UV abundances inherit
the flag. The audit gate (`pipeline/audit_solar_reference`) fails loud if a
UV/composite source is ever tagged `measured`.

- Source: CALSPEC `sun_reference_stis_002.fits`
  (https://archive.stsci.edu/hlsps/reference-atlases/cdbs/current_calspec/)
- Citation: Colina, Bohlin & Castelli 1996 (AJ 112, 307); Bohlin, Dickinson &
  Calzetti 2001 (AJ 122, 2118).
- Coverage: 119.5–2695.7 nm, **vacuum** wavelengths, FLAM (erg/s/cm²/Å, absolute).
- Resolution: **low** (Δλ ≈ 20 Å, R ≈ 150–300) — a flux composite, **not a line
  atlas**. It anchors UV flux/continuum and extends below the 296 nm Kitt Peak floor;
  it cannot resolve the diagnostic lines (use Kitt Peak for those).

## Composite sub-sources (the honesty breakdown — see `uv_provenance_rya459.json`)

| range (Å) | source | kind |
|-----------|--------|------|
| 1195–4100 | Woods et al. 1996 | **cited UV composite** |
| 4100–8700 | Neckel & Labs 1984 | measured (ground) |
| 8700–9600 | Arvesen et al. 1969 | measured (ground) |
| 9600–26950 | Castelli | model |

The whole product is tagged `cited-composite` because it blends cited-UV + model with
measured-visible; the UV (<4100 Å) — where NH 3360 / CN 3883 / UV C-N-O live — is the
cited part. The Kitt Peak atlas (measured, this ticket) is the resolved-line anchor
for everything ≥296 nm; CALSPEC carries only the deep-UV (<296 nm) and the flux scale.

## Columns (`sun_calspec_composite.csv`)

`wavelength_vac_A`, `wavelength_air_A` (Birch & Downs 1994), `flux_erg_s_cm2_A`,
`syserror`, `fwhm_A`, `provenance` (= `cited-composite`, on every row).
