# Nandakumar et al. 2024 IGRINS H+K package (RYA-1056)

This directory preserves and normalizes CDS catalogue `J/A+A/684/A15` for
Nandakumar et al. (2024), DOI `10.1051/0004-6361/202348462`.

Run `python3 scripts/rya1056_ingest_igrins.py --acquire` to retrieve the 23
official CDS files and regenerate the two normalized tables, element delta,
and checksummed provenance record. Without `--acquire`, the script regenerates
derivatives exclusively from the vendored raw package.

The catalogue contains 50-star line-by-line abundances at 76 wavelengths for
21 elements. It does **not** contain excitation energies, level identities,
log(gf), gf references, or HFS/isotopic component data. Consequently every
canonical wavelength coincidence is explicitly a discovery candidate/HOLD;
none is a transition-level cross-match or a gf-grade promotion. Those fields
remain blank rather than being guessed.
