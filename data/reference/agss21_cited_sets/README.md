# Asplund-cited solar line sets for the remaining elements (RYA-1268)

AGSS21 (Asplund, Amarsi & Grevesse 2021, A&A 653, A141) re-analysed only Li, CNO, Na, Mg,
K, Ca, Fe (and Al/Si via Nordlander & Lind 2017 / Amarsi & Asplund 2017). For every other
element it ADOPTS the line set of one of three papers, with the changes listed in its
Sect. 2.4-2.6. This folder holds those sets, transcribed exactly.

| file | what | source |
|---|---|---|
| `tex/scott2015_NatoCa_1405.0279.tex` | authors' LaTeX (arXiv e-print, latin-1 -> utf-8) | Scott et al. 2015, A&A 573, A25, doi:10.1051/0004-6361/201424109 |
| `tex/scott2015_FePeak_1405.0287.tex` | same | Scott et al. 2015, A&A 573, A26, doi:10.1051/0004-6361/201424110 |
| `tex/grevesse2015_Heavies_1405.0288.tex` | same | Grevesse et al. 2015, A&A 573, A27, doi:10.1051/0004-6361/201424111 |
| `agss21_cited_lines.csv` | Table 1 of all three papers, 363 lines | `scripts/rya1268_transcribe_agss21_cited.py` |
| `agss21_cited_gf_refs.csv` | each paper's numbered gf reference list (`gf_ref` key) | same |
| `agss21_tableA1_NaKCa.csv` | AGSS21 Table A.1, Na / K / Ca rows (Mg is RYA-1250's) | `scripts/rya1268_transcribe_agss21_tableA1.py` on local `aa40445-21.pdf` |
| `agss21_reference_lines.csv` | the Asplund set per element, AGSS21's changes applied, NIST-matched | `scripts/rya1268_coverage_audit.py` |

Note the naming: AGSS21 calls the Sc-Ni paper "Scott et al. (2015a)" and the Na-Ca paper
"(2015b)", the reverse of their A25/A26 publication order. Here the files are named by
content, never by letter.

## Transcription

The tables are parsed from the authors' LaTeX, not the PDF: PDF text extraction scrambles the
multi-column rows (level designations, phantom signs). Wavelengths are the printed air nm x 10;
EW stays in pm. The `*` weight marker (gf weight reduced for a large gf uncertainty) is kept as
`gf_weight_reduced`. Footnoted rows (multiplets computed as one line, forbidden lines,
synthesis-updated rows) carry `blend_note`.

AGSS21 Table A.1 prints wavelengths to 0.01 nm; the 4-decimal wavelength is carried from the
Scott+2015 transcription where the same line is within 0.1 A. Ca I 616.93 (a) is the
616.9042+616.9563 nm doublet with combined gf and EW (Table A.1 note a) and keeps its printed
wavelength. Cross-check: the parser reproduces RYA-1250's hand-transcribed Mg rows exactly.

## What AGSS21 changed (column `agss21_gf_status`)

- `AS_PRINTED`: AGSS21 used the printed gf.
- `UPDATE_PENDING`: AGSS21 replaced the printed gf with a newer lab source we do not yet hold
  (Sc: Lawler+2019; V I: Lawler+2014; V II: Wood+2014a; Cr II: Lawler+2017; Co I:
  Lawler+2015). The printed gf is NOT the Asplund gf for these lines; no sigma is attached.
- `NOT_TABULATED`: Li (Wang+2021, one resonance doublet) and Eu (Lawler+2001 / Sneden+2009
  via HM synthesis) have no line table in AGSS21 or Grevesse+2015.

## Per-line gf sigma (`loggf_sigma_dex`)

None of the three papers prints a per-line gf error. A sigma is attached only from the NIST ASD
accuracy class of the SAME transition (species, wavelength within 0.05 A, lower level within
0.02 eV) and only when NIST's log gf is within 0.03 dex of the published one. A NIST class on a
different gf is recorded (`nist_loggf`, `d_loggf_nist`) but is not this gf's error. The
percent-to-dex bridge is `pipeline.gf_grades.nist_sigma_dex`.
