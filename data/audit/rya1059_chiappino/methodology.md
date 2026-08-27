# Chiappino et al. 2026 methodology memo

- CRIRES+ settings: J1226 (1116–1356 nm), H1582 (1484–1854 nm), K2166 (1921–2472 nm).
- Observations: 0.4 arcsec slit, R approximately 50,000, final S/N per resolution element at least 40.
- Reduction: CR2RES v1.4.1 (dark/flat correction, nod-pair sky subtraction, arc wavelength calibration, optimal 1D extraction).
- Spectral synthesis: TURBOSPECTRUM in LTE with MARCS atmospheres.
- Atomic transitions: VALD3, including HFS entries for odd elements; the published table does not give per-line underlying references.
- Molecules: B. Plez online compilation. C uses 12CO, N uses CN, O uses OH; 13CO is retained separately for 12C/13C.
- Continuum/photon noise: continuum-placement uncertainty is stated as 1–2%. Multi-line random uncertainty is standard deviation / sqrt(N); single-line species receive 0.10 dex.
- Systematics: perturbations of ±50–100 K Teff, ±0.2 dex log(g), ±0.3 km/s microturbulence; reported abundance responses are generally <0.15 dex.
- Transfer firewall: the Liller 1 stars are cool RGB stars. Line usability, blends, continuum, LTE behavior, and telluric exposure do not transfer automatically to FGK dwarfs.
- Telluric handling and normalization are not specified with enough operational detail in the source text to inherit; every Codex spectrum still requires its own verified correction/provenance.
