# Kurucz line calculations (RYA-1354)

`gfZZII.pos.gz` is gzip -9 of http://kurucz.harvard.edu/atoms/ZZII/gfZZII.pos, retrieved 2026-11-18. ZZ is the atomic number and II the ion stage, with 00 for neutral and 01 for singly ionised.

| Code | Species | Code | Species |
|---|---|---|---|
| 1400 | Si I | 2400 | Cr I |
| 2000 | Ca I | 2401 | Cr II |
| 2200 | Ti I | 2500 | Mn I |
| 2201 | Ti II | 2501 | Mn II |
| 2300 | V I | 2601 | Fe II |
| 2301 | V II | 2700 | Co I |
| | | 2701 | Co II |
| | | 2800 | Ni I |

Each file holds lines between experimentally known levels. Every row is a Kurucz calculation, with reference codes K07–K18; no laboratory values are substituted, and the loader raises on any other code. The format is the gfall 160-column format. Fe I is in `../kurucz_fe1/`.

Al I has no `gf1300.pos` (HTTP 404); its lab file `../al_gf_lab/al1_lab_loggf.csv` carries Kurucz 1995 values in a column. The calibration that uses these files is `scripts/rya1354_theory_calibration_kurucz.py`.
