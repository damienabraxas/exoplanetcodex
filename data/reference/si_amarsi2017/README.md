# Si: cited solar 3D NLTE per line (RYA-1350)

`si_amarsi2017_lines.csv` joins two published per-line tables on the same disk-centre
equivalent widths (Scott et al. 2015a, Kitt Peak + Jungfraujoch average):

- Amarsi & Asplund 2017, MNRAS 464, 264, **Table 2**: `A_3d_lte` (lg e3L, = Scott's 3D LTE),
  `delta_3nlte_3lte` and `A_3d_nlte` (lg e3N) for the nine Si I lines A21 adopts; weights and
  W as printed. Si II 637.137 nm: the authors assume its 3D non-LTE correction is zero
  (printed in parentheses); A21 adopts that, and so do we, with a note on the row.
- Scott et al. 2015a (`data/reference/scott2015/scott2015_line_models.csv`): `A_mean3d_lte`,
  `A_hm_lte`, `A_marcs_lte` for the same lines.

Checked on build: every `A_3d_lte` equals Scott's `A_3D` to 0.0015 dex; every row joins to
one `canonical_gf.csv` line_id of the same species. The registered correction
(`solar_nlte_tables.yaml`) is `A_3d_nlte - A_marcs_lte`, +0.041 to +0.071 dex: cited and
added to our LTE fit, never presented as an in-transfer 3D calculation.

Al has no counterpart: Nordlander & Lind 2017 (A&A 607, A75) give per-line solar 3D NLTE
results only in Fig. 8, so Al's ENGINE-A-3DNLTE stays HELD.
