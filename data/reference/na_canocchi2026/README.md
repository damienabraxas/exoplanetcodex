# Solar Na: Canocchi et al. (2026) 3D NLTE

RYA-1251. Canocchi, Wang, Amarsi, Lind & Racca 2026, "3D NLTE Sodium abundances in
late-type stars", arXiv 2603.22975; grid and networks Zenodo 10.5281/zenodo.19201830
(CC BY 4.0; md5 of every file checked against the record). Balder 3D NLTE on Stagger,
Lind+2011 atom, MARCS 1D baseline -- the same atom and code family as Asplund+2021.

`na_canocchi2026_solar.csv` is written by `scripts/rya1251_na_canocchi2026_3d.py` from the
authors' networks (path register `grids.canocchi2026_na`): 3D NLTE minus 1D LTE at the Sun
(5772/4.44/0.00, vmic 1.00 km/s, A_1DLTE 6.25), seven lines. 6154/6160/10746 are cited by
AGSS21 Table A.1 (`../na_agss21/`) and stay there; their Canocchi values are a cross-check in
`data/results/rya1251/na_canocchi2026_validation.json` (Canocchi -0.035/-0.061/-0.009 vs
AGSS21 -0.012/-0.015/+0.027: 0.02-0.05 dex, different solar 3D models and 1D baselines).

Network noise is ~0.02 dex (the correction is not monotonic in the input abundance; the
paper's own held-out-Sun check is <= 0.025). The d_a1d_* / d_vmic_* columns record it.
8183/8194 carry the large corrections (-0.21/-0.18), as the paper finds.
