# One-valence-electron referees and levels (RYA-1354 v63/v64)

| File | Content |
|---|---|
| NaI_levels.csv, MgII_levels.csv, CaII_levels.csv, KI_levels.csv | NIST ASD energy levels, CSV export 2026-11-19 |
| mchf_mg2/Na_12.9.U.lsj.db.gz, Na_12.16.mchf_cv.lsj.db.gz | Froese Fischer B-spline CI [9] and MCHF-CV [16] for Mg II, NIST MCHF collection (nlte.nist.gov/MCHF/Elements/Mg), gzip -9 -n |
| safronova2011_ca2/table5_ca2_A_f.csv | Safronova & Safronova 2011 (arXiv 1011.5949, local PDF SafronovaSafronova_2011_CaII_oscillator_strengths_arXiv1011.5949.pdf) Table V, p. 9: 55 Ca II E1 lines, vacuum wavelength, A, f, uncertainty %. Transcribed from the PDF text layer; all 55 rows satisfy f = 1.49919e-16 lambda^2 (g_u/g_l) A to 1%, and the count matches the paper's "55 np-n's, np-n'd, and nd-n'f transitions". |

NIST ASD *lines* queries returned a server error on 2026-11-19 ("Can't use an undefined value as an ARRAY reference at lines1.pl line 701"), so K I has levels but no referee yet.
| MgI_levels_percentages.csv | NIST ASD Mg I levels with leading-component percentages (perc_out) and references, CSV export 2026-11-20 (the Landé g column is not served for Mg I: "Invalid Column Setting"). Referee for CI composition and singlet-triplet purity (RYA-1354 v68/v69). |
