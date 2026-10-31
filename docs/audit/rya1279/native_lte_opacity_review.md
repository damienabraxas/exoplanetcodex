# Native LTE normalization and MARCS opacity intake

RYA-1281 / parent RYA-1278, 2026-10-28. Mac-only diagnostic work; scientific
admission remains **HOLD**. No new SE/RT solve or production emission.

## What was tested

The actual LTEPOP and GETWRD routines from the held pymulti/MULTI source commit
4d6aa3c39d634b38da772e469a5692435eba1354 were compiled with their native BLOCK DATA
constants. Original source files are verified against Git. DP include conversion
follows the native make script; MK capacity increases from30 to64 to hold O41.
No physical level, collision, cross-section or abundance was changed. The wrapper
sets the oxygen total to1 and returns normalized fractions. It deliberately
does not build native opacity, EOS or radiation transfer. Molecular QPART is an
unreachable stop guard. Existing GNU Fortran16.1.0 was used with bounds checks;
this is not evidence of the departure-grid producer's compiler or source version.

The exact native5750/4.5 atmosphere and full held O41 atom were used:41 levels,
56 depths, identical temperature and electron density supplied to both routines.
Local LTE fractions come from the hash-verified converged sourced-MN cap case,
not a new solve. Both native and local arrays are normalized by their own total.

| Quantity | Result |
|---|---:|
| Maximum absolute native/local LTE-fraction ratio minus1, all levels/depths |0.00061373818348|
| Top O II ground native/local LTE fraction |0.99950332258151|
| Native fraction conservation maximum error |5.5511151231e-16|

Thus the largest relative difference is0.0614%; top O II differs by0.0497%.
At **matched inputs**, the held LTEPOP implementation cannot explain the roughly
twofold O II departure difference remaining in the finite-cap controls. This
does not establish the producer's actual ne, oxygen total, EOS or LTE denominator.
It does not establish an absolute ion-fraction error from departures alone.

An independent non-square three-level/four-depth Saha-Boltzmann test covers
neutral excitation, electron phase-space factor, source constants, Debye lowering,
cm^-3/m^-3 conversion, conservation and COMMON-block reset across calls. Eight
invalid-input cases reject before entering native code. The first analytic test
failed at2.53e-11 relative because it promoted each literal separately; the native
Fortran constant-only product4.98E-4*5040 is default REAL before DP division.
Preserving that inspected arithmetic fixes the discrepancy without relaxing the
3e-14 test tolerance. An initially mistyped existing test filename ran no tests;
the corrected focused selection passes44 tests.

## Independent MARCS continuum reference acquired

Holdings registry and local grid inspection found no opacity companion. A bounded
Sirius file-metadata query found no exact `.opa` companion in its MARCS directory;
no Sirius science was started. Its Gerber inventory also exposed no native
producer ne/nstar/J/rate dump. These bounded searches are not exhaustive proof
that no such output exists elsewhere.

The official [MARCS search](https://marcs.oreme.org/data/) returned the exact model
identifier. Requesting that `.mod.gz` identifier with file_type=[opa] yielded one
391918-byte gzipped member in a399360-byte tar. The initial request with a bare
stem returned an empty archive plus an explicit missing-file error; it was never
used as science data. The successful archive, response and source are hash-receipted.
The complete grid bundle was not downloaded.

The [official reader](https://marcs.astro.uu.se/documents/auxiliary/readopa.f)
specifies1071 wavelength points,56 depth rows, interleaved absorption/scattering
ratios to the standard opacity OPS. The parser multiplies each coefficient by
OPS before converting cm2/g to m2/kg. Tests cover that interleaving and unit
conversion, zero microturbulence, malformed/truncated/extra/nonfinite input and
negative opacity. No coefficient is substituted into the solver.

At every depth, temperature, electron pressure and density match the held `.mod`
**exactly at printed precision**. All92 abundances and microturbulence also match.
The maximum |log10(tau_opa)-logtau_mod| is5.44513e-5; both printed scales are
retained, with no remapping or forced equality. This establishes a useful matched
structure reference, not the native MULTI producer's background opacity.

The inspected reader does **not** state the wavelength medium. Official docs
explicitly call the separate `.flx` wavelength list vacuum, but that statement
cannot silently be transferred to `.opa`. OPA source coordinates remain marked
UNVERIFIED; no wavelength-aligned cross-engine opacity or rate claim is made yet.
This is separate from the already established Kitt Peak **air** comparison policy.

## Evidence and next work

Native oracle: logical path experiments.bride_multi_lte_oracle, external
bride_inputs/rya1281/native-lte-v1. Final audit: bride_mac_checks/rya1281-native-lte-v2;
v1 remains immutable. Review copies live in data/results/rya1281/native_lte_v1.
OPA source/reader: experiments.bride_marcs_opacity, external
bride_inputs/rya1281/marcs-opacity-v2; initial stagingv1 is retained. Its receipt
is copied to data/results/rya1281/marcs_opacity_v1. GENERATORS.yaml records commands.

Next establish the OPA wavelength medium from an actual source contract, then
compare independent continuum absorption/scattering against the local background
at the verified depths, with active oxygen accounted for explicitly. Keep this
distinct from missing MULTI ABSDAT/native ne/nstar/J/rates. Do not fit opacity,
electron density, abundance or cap multipliers to the reference. Validate native
column interpolation before capped solar-parameter reruns. Then finish solar
oxygen diagnostics and proceed C/N, Fe, then another star. Full3D/all-band and
native producer parity remain open; no adopted value, default boundary, gold,
production product or verdict changed.
