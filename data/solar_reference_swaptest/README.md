# Solar O I 777 swap-test reference extracts (RYA-485/486)

777-triplet windows (air ~7762–7784 Å) extracted from the four independent solar
references for the RYA-484 Lever-3 continuum-lever swap-test. Each FTS file (Reiners,
Baker, Wallace) is wavenumber (cm⁻¹, vacuum): column 1 = wavenumber, column 2 =
normalized/telluric-removed flux; the loader applies 1e8/ν → vacuum Å → air (Birch&Downs).

- `reiners_oi777_vac.txt` — Reiners 2016 IAG base (CDS J/A+A/587/A65, spvis.dat nFlux)
- `baker_oi777_vac.txt`   — Baker 2020 telluric IAG (Zenodo 3598136, iag_telfree s)
- `wallace_kpno_oi777_vac.txt` — Wallace 2011 telluric-removed KPNO (NSO sptr.reg1 col2)
- (Kurucz 1984 KPNO is data/solar_reference/kpno_flux_atlas/kpno_OI_777_triplet.csv)

CANONICAL full atlases live on Sirius /mnt/codex-data/solar_reference/ (RYA-481/485);
these are the small reproducible 777-window extracts for the committed swap-test result.
