# Scientific references

This bibliography was reconciled against the project's local reference library
and live code/provenance records. It is grouped by role so a reproduction cites
only resources actually used. Local PDFs are working papers—not redistributable
repository assets—and the upstream publication remains authoritative.

## Scientific methodology and benchmark analysis

- Gray, D. F. (2005), *The Observation and Analysis of Stellar Photospheres*.
- Sousa et al. (2011), A&A 533, A141: homogeneous equivalent-width stellar
  parameters and line selection.
- Jofré et al. (2014), A&A 564, A133; Heiter et al. (2015), A&A 582, A49:
  Gaia FGK benchmark metallicities and fundamental parameters.
- von Braun et al. (2011), ApJ 740, 49: interferometric 55 Cancri parameters.
- Teske et al. (2013), ApJ 778, 132: 55 Cancri abundance context.
- Schmitt (2010), thesis: historical analysis context; not a substitute for
  current pipeline provenance.

## Software and radiative transfer

- Blanco-Cuaresma et al. (2014), A&A 569, A111: iSpec.
- Alvarez & Plez (1998), A&A 330, 1109; Plez (2012), Astrophysics Source Code
  Library `ascl:1205.004`; Gerber et al. (2023): Turbospectrum lineage.
- Sneden (1973), ApJ 184, 839: MOOG.
- Astropy Collaboration (2013, 2018, 2022): Astropy.

## Model atmospheres and solar scale

- Castelli & Kurucz (2004), arXiv:astro-ph/0405087: ATLAS9 ODFNEW grids.
- Gustafsson et al. (2008), A&A 486, 951: MARCS atmospheres.
- Magic et al. (2013), A&A 557, A26: STAGGER 3D atmosphere grid.
- Asplund et al. (2021), A&A 653, A141: canonical solar abundances.

## NLTE and 3D resources

- Amarsi, Nissen & Skúladóttir (2019), A&A 630, A104: C I/O I 3D- and
  1D-NLTE correction tables.
- Amarsi et al. (2022), A&A 668, A68: Fe 3D-NLTE/1D error models.
- Lind et al. (2011), A&A 528, A103: Na NLTE/INSPECT.
- Bergemann & Cescutti (2010), A&A 522, A9; Bergemann (2011), MNRAS 413,
  2184: Cr and Ti NLTE.
- Mashonkina et al. (2017), A&A 606, A147: Ca NLTE.
- Korotin et al. (2015), A&A 581, A70: Ba II NLTE.

The grid-specific `.prov.json` files in `data/nlte_grids/` are authoritative for
the exact artifact. Planned or archived resources must not be cited as applied.

## Atomic and molecular data

- Ryabchikova et al. (2015), Phys. Scr. 90, 054005: VALD3.
- Kramida et al., NIST Atomic Spectra Database, version recorded at retrieval.
- Den Hartog et al. (2011), ApJS 194, 35, and element-specific laboratory
  studies recorded in the canonical oscillator-strength store.
- ExoMol line-list publications and release manifests for each isotopologue;
  cite the exact dataset named by the molecular manifest.

## Instruments and archives

- Mayor et al. (2003), *The Messenger* 114, 20: HARPS.
- Cosentino et al. (2012), SPIE 8446: HARPS-N.
- Dekker et al. (2000), SPIE 4008: UVES.
- Pepe et al. (2021), A&A 645, A96: ESPRESSO.
- Donati et al. (2020), MNRAS 498, 5684; Cook et al. (2022), PASP 134,
  114509: SPIRou and APERO.
- HST STIS and COS Instrument Handbooks, version recorded with the run.
- ESO CRIRES+, NIRPS, FEROS, and pipeline manuals, version recorded with the
  reduced product.

Archive citations/acknowledgements are distinct from instrument papers:

- ESO Science Archive and the program IDs/product identifiers used.
- MAST and the HST program/product identifiers used.
- TNG/IA2, Keck Observatory Archive, CADC, PolarBase, Gemini Observatory
  Archive, NOIRLab Astro Data Archive, and IRTF/IRSA as applicable.

## Solar/reference datasets

- Kurucz, Furenlid, Brault & Testerman (1984), NSO Atlas No. 1: Kitt Peak
  Solar Flux Atlas.
- Colina, Bohlin & Castelli (1996) and Bohlin et al. (2001): CALSPEC solar
  composite/absolute-flux lineage.
- Hase et al. (2010) and the ACE-FTS solar atlas; NSO/Wallace infrared atlases
  where named in `data/solar_reference/ir_atlases/`.

Always record whether a product is measured, a cited composite, or a model.
HST cannot directly observe the Sun; the CALSPEC solar spectrum is never
described as a direct HST solar measurement.
