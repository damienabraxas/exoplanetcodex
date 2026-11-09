# RYA-1281: full-column H2 source and chemistry check

Scientific HOLD. This adds a source-qualified neutral-H2 equilibrium interface,
not a full EOS, H2 opacity, new oxygen spectrum, or production ORCH admission.
All work ran on the Mac; previous oxygen spectra remain separate.

## Source and units

Inspected held Barklem & Collet 2016 paper, DOI10.1051/0004-6361/201526961,
Eq7 and Sec3, and the author repository ReadMe/correction notice. Existing
`cno_molecular_primary/constants_barklem2016` holds dissociation energies only;
that frozen partial holding was not overwritten. New source root resolves via
`experiments.bride_molecular_tables`. Eight source files are pinned to author
commit `cd2f54ed46687659017b343763e95478ed93653e`, checked against Git blob IDs,
accepted read-only, and individually SHA256-receipted. BSD3-Clause notice retained.

Source: https://github.com/barklem/public-data/tree/cd2f54ed46687659017b343763e95478ed93653e/partition-functions_and_equilibrium-constants

Use corrected November2022 equilibrium/atomic-partition tables. Original and
corrected H2 equilibrium rows are exactly equal; correction affects other species.
Charged-molecule conventions are outside this neutral-H2 adapter. Adopted source
D0=4.478007eV. Kp=p(HI)^2/p(H2), in Pa, hence
Phi=n(H2)/n(HI)^2=kT/Kp in m3. No extra statistical-spin factor is applied.

Paper Sec3 recommends a cubic spline of log10(Kp)+D0*5040/T against log10(T).
This adapter explicitly uses not-a-knot boundaries, qualified1000-10000K, and
refuses extrapolation including a one-ULP domain excursion. The source holds
lower temperatures too; those are not qualified by this implementation.
The paper's reported1e-4 accuracy concerns its adaptive31-point grid, not a
verified bound on this rounded42-point table or this solar experiment.

## Numerical results and limits

Final immutable output: `experiments.bride_mac_grid/rya1281-h2-tables-v2`.
Earlier v1 remains preserved. Compact byte-identical evidence lives at
`data/results/rya1281/h2_tables_v1/`.

All56 bound solar5772/g4.438 rows, T4062.4452-9930.0128K, are covered, including
the four previously excluded hot rows. Native H6/He9 LTE states at unchanged
T/ne/nH/composition supply atomic ionization and H-minus ratios. Solving the
positive quadratic depletes HI when forming H2 and conserves HI+HII+Hminus+2H2.
Ne is fixed; this does not solve charge neutrality or a metal/CNO network.

- Maximum H2/total hydrogen nuclei:0.00078353268 (0.07835%).
- Maximum nuclei closure residual:4.4409e-16.
- Independent bounded Brent root: maximum absolute HI-fraction difference2.554e-15.
- Natural versus not-a-knot boundary choice changes solar Phi by at most5.462e-5 relative.
- Interior leave-one-knot-out max relative Phi differences: not-a-knot0.00074741,
  natural0.00175331, PCHIP0.00727370. These measure sensitivity to removed source
  knots, not physical uncertainty or a universal interpolation error bound.
- Independent Eq7 reconstruction with modern SI constants, source partition
  functions and LW atomic-H mass differs from tabulated log10(Kp) by
  +7.4936e-5 to+5.4178e-4dex. This small finite difference is **unassigned**;
  no claim of an upstream bug, no fitted constants, and no correction to the
  published Kp. Author H2 input also confirms the adopted D0. Resolve the
  source calculation/constants before claiming equation-level exact parity.

Six tests pass: three new source/unit/domain/independent-network/tamper tests
and three existing Stark/source/context/domain regressions. Native environment
warnings remain unchanged. Synthetic network controls span1000-10000K,
1e12-1e28m^-3 and multiple ionization/Hminus ratios; they are algorithm checks.

## Reproduction and next work

Load `config/hosts/mac.env`, use held Bride Python, then run:

```
python scripts/audit_bride_h2_tables.py --output <new-immutable-directory>
python -m pytest -q tests/test_bride_h2_tables_rya1281.py tests/test_bride_fe_stark_rya1281.py
```

For source re-intake in a new registered root, acquire the eight exact URLs in
`source_intake.json` with `.download` suffix; obtain author directory and LICENSE
GitHub contents metadata at the recorded commit. Run
`intake_bride_h2_tables.py --metadata <directory-json> --license-metadata <license-json>`.
Existing sources/receipt refuse overwrite. No raw tables are vendored here.

Next: source reconstruction closure and a paired H2 collider-only pressure
control with HI depletion explicitly installed. Such a control must remain
labelled partial: O collision/continuum and full EOS consistency are separate
requirements. Then fresh matched xi, broader opacity/continuum anchor controls,
and an independent matched1D O reference. No abundance/continuum/damping fit,
new observational correction, merge, publication, full3D or all-band claim.
C/N precede Fe abundance and another star.
