# Solar Ca: Asplund, Amarsi & Grevesse (2021)

RYA-1311. Source: A&A 653, A141, DOI 10.1051/0004-6361/202140445, local `aa40445-21.pdf`,
Table A.1 (Ca block: 9 Ca I features + 5 Ca II), all eight model columns, taken from
`data/reference/agss21_cited_sets/agss21_tableA1_NaKCa.csv`. Same layout as `na_agss21/`.
Sect. 2.4: full 3D non-LTE, 141 Ca I + 46 Ca II level atom; 616.13 dropped; recommended
6.30 +/- 0.03 (total).

Canonical joins: species, wavelength within 0.1 A, lower energy within 0.01 eV.
- 616.93 is the 616.9042 + 616.9563 doublet (Table A.1 note a, combined gf -0.308 and EW).
  Both components (gf_012415, gf_012416; Den Hartog 2021 lab gf -0.797/-0.478, sum -0.308)
  carry the doublet's per-feature abundances.
- Ca II 732.39 (forbidden) has no canonical row and is not in iSpec's GES v6 red-optical
  list: not measured, left out here (declared gap, RYA-1311).

Registered in `data/reference/solar_nlte_tables.yaml`: 1D leg = MARCS NLTE - MARCS LTE,
3D leg = 3D NLTE - MARCS LTE. Disc-centre intensity EWs; additive corrections to our fit.
