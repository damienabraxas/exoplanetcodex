# Na I independent theory and levels (RYA-1354 v61/v62)

Downloaded 2026-11-19 from the NIST MCHF/MCDHF collection (https://nlte.nist.gov/MCHF/, Elements/Na) and NIST ASD.

| File | Content |
|---|---|
| raw/Na_11.9.U.lsj.db.gz | Froese Fischer 2002, non-orthogonal B-spline CI (ref [9]); 385 E1 lines, length and velocity rows. This is the source of NIST's adopted value for every solar Na I line (KP08 Table 4, source "2"). |
| raw/Na_11.9.U-lev.lt.db.gz | its levels and lifetimes |
| raw/Na_11.16.mchf_cv.lsj.db.gz | Froese Fischer MCHF with core-valence correlation (ref [16]), 73 lines |
| raw/Na_11.16.mchf_cv-lev.lt.db.gz | its levels and lifetimes |
| raw/LS_Na_Z11.gz | the LS convergence trends for [16] |
| raw/nist_asd_na1_levels.csv | NIST ASD Na I energy levels (CSV export, 2026-11-19) |

The .db files are gzipped byte for byte (`gzip -9 -n`).
sha256 of the .gz files is in `docs/audit/rya1354/na_model_v62.md` provenance and was recorded at download.

Companion paper (local reference folder): Kelleher & Podobedova 2008, J. Phys. Chem. Ref. Data 37, 267,
`KelleherPodobedova_2008_JPCRD37_267_NaMg.pdf` (from srd.nist.gov/jpcrdreprint/1.2735328.pdf).
