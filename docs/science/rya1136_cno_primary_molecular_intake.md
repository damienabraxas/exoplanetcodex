# RYA-1136 primary CNO molecular intake

## Result

The complete Amarsi et al. (2021) Table 2 abundance-feature census contains
408 rows. All 80 CO rows retain their exact Li et al. (2015) physical-tuple
join. The acquired Brooke/Masseron primary holdings provide transition-level
evidence for the other 328 rows:

- 278 unique primary-transition matches;
- 32 explicitly summed unresolved-component matches;
- 9 rounded component sets for which Table 2 lacks the rotational label;
- 5 secure wavelength/band/energy identities with a source-release strength
  difference;
- 3 OH features absent from the acquired Brooke release; and
- 1 OH feature whose published Table 2 energy is inconsistent with the only
  same-band transition at that position.

Thus 390/408 rows have exact or explicitly summed physical identities. The
remaining 18 are retained as review rows, not silently promoted and not
dropped from the published-used census.

## Normalization rules

- Matching uses molecular-system and vibrational-band identity, wavenumber,
  lower-state energy, and oscillator strength. Wavelength-only joins are
  prohibited.
- Brooke CN, NH, OH, and C2 files tabulate `f`; conversion to `gf` applies the
  lower-state statistical weight `(2 J_lower + 1)`.
- Masseron CH directly tabulates `gf`.
- Brooke C2 lower energies are relative to the Swan lower-state band origin.
  A 0.0753 eV term-origin conversion puts all 39 Amarsi energies onto their
  published excitation scale.
- Unresolved feature strengths are sums in linear `gf`, never sums in log
  space. The exact source components travel with every admitted sum.
- Positions are compared in the primary lists' native wavenumber coordinate.
  A fixed Angstrom tolerance is invalid over the 0.4--15 micron span.

## Reproduction

Run, in order:

```console
python3 scripts/ingest_cno_molecular_primary_rya1136.py
python3 scripts/build_cno_intake_rya1136.py
python3 -m pytest -q tests/test_cno_primary_ingest_rya1136.py \
  tests/test_cno_intake_rya1136.py tests/test_cno_closure_rya1136.py
```

The compact evidence ledger is
`data/audit/rya1136_cno_intake/primary_molecular_crossmatch.csv`. The source
archives remain under `data/reference/cno_molecular_primary/`; the build does
not create an undocumented second line-list release.

## Atomic nitrogen closure

AGSS21 adopts the five N I lines from Amarsi et al. (2020) Table 1. Their air
wavelengths are 744.229, 821.633, 862.923, 868.340, and 1010.890 nm. The audit
now carries all five rather than the earlier three-line control subset.
