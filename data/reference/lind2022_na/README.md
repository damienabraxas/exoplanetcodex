# Lind et al. (2022) Na curves of growth at the solar-like node

RYA-1251. Lind, Nordlander, Wehrhahn et al. 2022, A&A 665, A33, CDS J/A+A/665/A33 `na.dat`
(corrected version 03-Jul-2026), extracted through the CfA VizieR mirror (cdsarc's direct file
is behind a bot wall): Teff 5750, logg 4.5, [Fe/H] 0.00, Vturb 1.0, all 35 lines,
[Na/Fe] -1.0..+0.8, LTE and NLTE log EW (mA). PySME + MARCS 1D NLTE -- the same method as
`scripts/extend_nlte_grid.py` -- so it is the independent check on the ENGINE-A extension
(RYA-410: agree within 0.04 dex at the Sun).

Solar NLTE-LTE at [Na/Fe] = 0 (same-EW inversion): 5149 -0.036, 6154 -0.056, 6161 -0.072,
10746 -0.007 -- vs AGSS21 Table A.1 MARCS -0.033/-0.047/-0.058/-0.006 (0.003-0.026 apart).
The banked ENGINE-A rows 5682/5688 (-0.121/-0.137) are 0.04 more negative than Lind's
-0.078/-0.096 -- at the RYA-410 limit; recorded, the banked table is not changed here.
