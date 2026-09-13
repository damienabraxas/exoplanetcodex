# Reference-folder content check — 2026-09-13

Checked `/Users/ryanschmitt/Documents/Exoplanet Codex/Reference documents` by
PDF file signature, not suffix: **93 PDF files, 84 distinct SHA256 hashes**.
Extracted available text across their pages, reviewed relevant passages and
first-page identities, and visually verified the scanned Garz title page.
This is a holdings/relevance audit, not full scientific adjudication of every paper.
Scanned pages without text require visual review; lack of keyword hits is not
proof that a document is irrelevant. Original reference files are unchanged.

## Directly applicable Si sources already held

| Held filename | Verified source | Use for RYA-1218 |
|---|---|---|
| `1973A&A....26..471G` | Garz 1973, A&A 26, 471–477, *Absolute Oscillator Strengths of Si I Lines between 2500 Å and 8000 Å* | Primary Si I measurements; seven-page scanned PDF with no extension. **Remove from missing-paper list.** Tables still need transcription. |
| `aa24109-14.pdf`; `1405.0279v2.pdf` | Scott et al. 2015, A&A 573, A25 | Published paper plus preprint; source selection, gf lineage, Si I/Si II comparison. |
| `1609.07283v1.pdf` | Amarsi & Asplund, solar Si 3D-NLTE analysis | Source tables and model assumptions; duplicate-named copy has identical bytes. |
| `Apslund 2021.pdf`; `aa40445-21.pdf` | Asplund, Amarsi & Grevesse 2021, A&A 653, A141 | Same PDF under two names; Si section and older theoretical-gf comparison. |
| `2301.11391v1.pdf` | Den Hartog et al. 2023, Si I/Si II transition probabilities | Modern primary laboratory gf; all three named copies are byte-identical. |
| `aa36291-19.pdf`; `ReadMeGaia-ESO.txt`; `ges_refs.bib` | Heiter et al. 2021, *Atomic data for the Gaia-ESO Survey* | Appendix B.8 explicitly distinguishes Garz + lifetime normalization, Nahar/Opacity Project, Kurucz, and the Si II experimental mean. Useful for decoding canonical provenance and blend flags. |
| `aa59148-26.pdf`; `2602.14294v1.pdf` | Elgueta et al. 2026, A&A 710, A111 | NIR line-selection/robustness evidence. Published and preprint versions differ; each has a byte-identical duplicate. |
| `TurboSpectrum.pdf` | Gerber et al. 2023, A&A 669, A43 | NLTE engine/departure-grid methodology; Table 1 identifies the Si model source as Bergemann et al. (2013). This does not replace the original atom paper. |
| `AbundancesLodders.pdf` | Lodders, Bergemann & Palme 2025, *Solar System Elemental Abundances from the Solar Photosphere and CI-Chondrites* | Material discussion of competing Si analyses, gf choices, ionization balance, and line-selection effects. Contains a reference to the published Pehlivan Rhodin et al. **2024**, A&A **682**, A184 Si transition-data paper. That primary paper was not identified in this folder. |
| `aa51889-24.pdf` | Sharma et al. 2024, A&A 691, A160 | Independent stellar Si line-analysis methodology; differential analysis must not supply solar-fitted gf to our reference pool. |
| `0004190v1.pdf` | Griesmann & Kling, *Interferometric measurement of resonance transition wavelengths in C IV, Si IV, Al III, Al II, and Si II* | UV physical wavelength identity; not a gf replacement or Solar abundance analysis. |

## Supporting material

- `Baker_2020_ApJS_247_24.pdf` and `aa27530-15.pdf`: corrected/raw IAG atlas provenance.
- `solar-spec.pdf`: Hase et al., ACE-FTS infrared solar atlas; useful IR context, not proof of a staged pipeline holding.
- `MARCS Stellar Model.pdf`, `3D Hydrodynamical Model Atmospheres (Magic et al. 2013).pdf`, and `atlas9atlas12.pdf`: atmosphere/model treatment.
- `aa26961-15.pdf` and `ReadMe PArtion functions for molecules and atoms.txt`: partition functions and equilibrium constants.
- `Nandakumar_2024_ApJ_964_96.pdf` and `aa48462-23.pdf`: IR giant-star Si methods; applicability to the Sun needs explicit evaluation.
- `2505.13607v1.pdf`: Koutsouridou et al., metal-poor-star NLTE/3D correction compilation; contextual source map, not an automatic Solar correction.
- `2511.04254v2.pdf`: Bergemann & Hoppe radiation-transfer review; methodological context.
- `Microturbulence.pdf` is actually the Asplund et al. 2009 solar-composition review. Conversely, `Apslund2010.pdf` is Mucciarelli's 2011 microturbulence paper. Filenames alone are unreliable.

## Still not identified as standalone primary papers

O'Brian & Lawler 1991a/b; Schulz-Gulde 1969; Blanco et al. 1995;
Matheron et al. 2001; the Si NIST critical compilation;
Pehlivan Rhodin et al. 2024 Si transition data; and the original Si
Bergemann et al. 2013 atom paper. References to them in held documents do
not establish that their primary tables are held.

`jpcrd372008911p.pdf` is the **aluminum** compilation, not silicon.
`1808.09478v1.pdf` is Al I/Al II theory, not Si theory.
The Oxford HTML table exports link DOI `10.1093/mnras/stab214` (the held
carbon transition-data paper); they are not the missing Si supplement.

The inventory CSV records filename, PDF page count, hash, extraction status,
and identical-byte siblings. No duplicate was deleted and no paper was downloaded.
Garz's existing bibliography entry now records its actual local filename.
