# Kurucz Fe I calculation (RYA-1354)

`gf2600.pos.gz` is gzip -9 of http://kurucz.harvard.edu/atoms/2600/gf2600.pos, retrieved 2026-11-18. The uncompressed sha256 is f06b83d4f3afb68864734fe46a3515da0e878a28a4d85ef832d7a95f1ebbfc4f.

- **Content:** 126,288 Fe I lines between experimentally known levels.
- **gf:** from Kurucz's least-squares-fitted (semi-empirical) calculation of 22 Oct 2014.
- **Reference code:** K17 on every row. No laboratory values are substituted.
- **Format:** the gfall 160-column format. Wavelengths are in nm, vacuum below 200 nm and air above.

`2600.readme` is Kurucz's own note on the calculation.

The theory-vs-lab calibration that uses this file is `scripts/rya1354_theory_calibration_fe.py`.
