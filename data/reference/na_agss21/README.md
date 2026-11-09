# Solar Na: Asplund, Amarsi & Grevesse (2021)

RYA-1251. Source: A&A 653, A141, DOI 10.1051/0004-6361/202140445, local
`aa40445-21.pdf`, Table A.1 on printed page 27 (Na I block, five lines),
transcribed with all eight model columns, as RYA-1250 did for Mg
(`data/reference/mg_agss21/`). Wavelengths are the printed nm x 10; EW in pm.
The Na discussion is Sect. 2.4, page 9: full 3D non-LTE with the Lind et al.
(2011) 23-level atom; "3D non-LTE effects are minor (~ -0.03 dex)";
recommended 6.22 +/- 0.01 (stat) +/- 0.03 (syst).

Canonical joins: species, wavelength within 0.1 A, lower energy within 0.001 eV;
one row each (gf_097218, gf_097221, gf_097227, gf_097228, gf_167446), and the
canonical log gf equals the printed value for all five.

Registered in `data/reference/solar_nlte_tables.yaml` as the cited solar
correction for these lines: 1D leg = MARCS non-LTE minus MARCS LTE (ENGINE-A),
3D leg = 3D non-LTE minus MARCS LTE (ENGINE-A-3DNLTE). Same solar EWs; additive
corrections to our 1D-LTE fit, not our own 3D transfer.
