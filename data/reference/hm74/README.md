# Holweger & Müller 1974 solar model (HM74)

`hm74_cowley_castelli2002.csv` -- Table A.1 of Cowley & Castelli 2002, "Some aspects of the
calculation of Balmer lines in the sun and stars" (A&A 387, 595; arXiv astro-ph/0201258),
Appendix A: the Holweger & Müller (1974, Sol. Phys. 39, 19) T-tau5000 relation interpolated to
50 depths and put on the Kurucz RHOX scale (kappa5000 from Holweger, Steffen & Steenbock's ATM
code). Parsed from the arXiv PDF text (local copy: Reference documents/
Cowley_Castelli_2002_HM74_model_astro-ph0201258.pdf), 50/50 rows.

Use (RYA-1315): the 1D reference of Asplund, Amarsi & Grevesse 2021 Sect. 2.1's
temperature-structure term, half the mean (<3D> - HM74), for species Scott+2015 does not tabulate
(atomic C I, N I, O I).

Note: the printed rho at log tau = -5.31 (4.085e-10) breaks the monotone run (3.772e-10 ->
6.072e-10) and is likely a typo in the source table; T, Pgas, Ne and RHOX at that depth are
monotone and are what a model uses. Kept as printed.
