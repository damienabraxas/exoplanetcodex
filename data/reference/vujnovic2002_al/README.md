# Vujnovic et al. 2002 Al I/Al II tables

Raw machine-readable tables from CDS catalog `J/A+A/388/704`, acquired
2026-08-30. Article DOI: `10.1051/0004-6361:20020560`.

Source: `https://cdsarc.cds.unistra.fr/ftp/J/A+A/388/704/`

The files under `raw/` are byte-identical to the CDS responses:

- `ReadMe`: `a423ac834f8e8dc2f76fa642ec13abfb301895bd60e1ae40a27d366de010c71e`
- `table2.dat`: `8b1e541937226668919b1a7108fe47e02ea5ef7af301294d34c8ee3fbeb7b4d4`
- `table3.dat`: `692472e8e6ed42823719d3a942988df65e010aed6f3e751c5f109d99d1fe57b9`
- `table4.dat`: `39617da6a9c0cf699cf6bc387b23650246b4e53e691c0d4a127e60a61b03eff8`
- `table5.dat`: `67b69d40f137699d98915547e222c1dccb6067f5f3abd2b137b60f7f2421b0b5`

`scripts/build_al_intake_rya1132.py` parses all 106 fixed-width records and emits
normalized and manifest-crossmatch ledgers. It derives log gf only for finite Aki
values with a parseable upper-level J. Source limits remain limits. Promotion is
restricted to the manually adjudicated finite Al I transitions with stated Aki
uncertainties; ratio-only rows and 3092.839 are not promoted.
