# AGSS21 Solar Al lineage (RYA-1134)

This directory records the provenance behind the Solar Al reference set. It is
not a new abundance result and does not overwrite the 505-row RYA-1132 intake.

AGSS21 adopts the `A(Al) = 6.43 +/- 0.03` 3D-NLTE result of Nordlander & Lind
(2017). Nordlander & Lind state that they use the Scott et al. (2015) Solar line
selection and weights, except for the telluric-contaminated 10891 A line. Scott's
published table contains seven Al I features, leaving six traceable features.

AGSS21 calls these “five” lines. Because neither AGSS21 nor Nordlander & Lind
names a second exclusion, the generated sidecar preserves all six implied lines
and labels the published count mismatch. The sidecar also carries the excluded
10891 A feature so the denominator and rejection are auditable.

Build with:

`python3 scripts/build_al_agss21_reference_rya1134.py`

Sources:

- Asplund, Amarsi & Grevesse (2021), DOI 10.1051/0004-6361/202140445,
  arXiv:2105.01661.
- Nordlander & Lind (2017), DOI 10.1051/0004-6361/201730427,
  arXiv:1708.01949.
- Scott et al. (2015), DOI 10.1051/0004-6361/201424109,
  arXiv:1405.0279.
