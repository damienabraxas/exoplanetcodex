# RYA-1220 ticket update — 2026-09-23

## Completed

- Acquired the missing MARCS 2012 grid on Sirius from the NADC mirror and configured PySME to use the local copy.
- Repaired the NumPy 2 PySME extension build and native SMElib loading.
- Found and fixed the native NLTE control bug: SMElib departure coefficients are global; the LTE branch was not resetting them, so LTE reused the previous NLTE coefficients.
- Validated the repair with Na and K controls:
  - Na I 5682.633: −0.113 dex
  - K I 7698.964: −0.304 dex
- Completed native solar N corrections:
  - N I 7468.31: −0.0096 dex
  - N I 8216.34: −0.0128 dex
  - N I 8683.40: −0.0136 dex
  - median: −0.0128 dex
- Confirmed Turbospectrum remains the blend-capable route for NH/CN molecular diagnostics.

## Coverage findings

- N I 7442.29 and 8629.23 are measured reference features but need canonical long-format VALD term/J records before native PySME registration.
- N I 10108.90 has no complete registered solar holding coverage.
- NH 3360 and CN 3883 are molecular near-UV routes and require the molecular Turbospectrum line lists; they are not atomic PySME inputs.
- No approved near-UV atomic N I line is currently registered.

## Remaining work

1. Acquire canonical long-format records for N I 7442.29 and 8629.23 and register their grid-resolved levels.
2. Complete the UV molecular NH/CN synthesis uncertainty products.
3. Propagate the validated atomic N corrections and molecular blend/continuum terms into the instrument-by-band uncertainty tables.
4. Run the final all-instrument covariance and promotion checks.

Artifacts:

- `data/results/rya1220/pysme_nitrogen/native_solar_sweep.json`
- `data/results/rya1220/pysme_nitrogen/native_nlte_control_crosscheck.json`
- `data/results/rya1220/pysme_nitrogen/line_coverage_next_step.json`
- `data/output/rya1220/nitrogen_full_synthesis_sweep.json`
