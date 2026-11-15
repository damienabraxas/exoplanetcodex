# RYA-1354 — source-mapped radial readiness, v31/v32

Scientific admission remains **HOLD / EXPOSED_NONBLIND**. Adoption is false; accepted lifetime, uncertainty, covariance and calibrated sigma remain null. This checkpoint contains orbital readiness only: no new CI, LS decomposition, rates, lifetime or grid-convergence result.

## Retained failure and correction

The v31 runner incorrectly treated the complete FAC SavePotential serialization as an immutable central potential. The first ground/coherent10 native run succeeded, but the integration guard rejected old_solved.pot. Its receipt remains FAILED_RETAINED, with all 138 native paths preserved in fac_radial_probe_failed_v31 (one 1,277,453-byte compressed part).

Pinned primary source maps the 348-byte header and 69 arrays of 1,500 doubles. Only arrays 22–24 (W, dW, dW2; byte offsets 264348–300347 inclusive) changed. SetPotentialW derives these work arrays from orbital energy and kappa. Every other byte is protected by the corrected comparator; mutated header, grid, central and QED fields are rejected in regression tests. Restore roundtrips still require whole-file identity. Raw files and original solver source/executable remain unchanged.

A separately declared v32 reran all four cases using the paired v28 potentials and original SFAC binary. All passed independent archive replay: nine old spinors retain exact physical r, rVT0, P, Q samples and header metadata after adding three bound spinors. Each case has 19 overlap pairs before and 36 after; all printed diagonal overlaps are unity. Off-diagonal overlaps and Hamiltonian asymmetry remain diagnostics, not accuracy tolerances. The WaveFuncTable cached effective-potential column is retained but is not confused with a newly evaluated orbital field.

## Fresh orbital energies (eV)

| Case | 4s | 4p kappa +1 | 4p kappa -2 |
|---|---:|---:|---:|
| ground_coherent10 | -1.66621773 | -1.25423725 | -1.25370754 |
| ground_coherent100 | -1.66621709 | -1.25423697 | -1.2537081 |
| equal_coherent10 | -1.70111416 | -1.28366694 | -1.28320287 |
| equal_coherent100 | -1.70111385 | -1.28366766 | -1.28320251 |

The v32 archive preserves 255 paths in one 4,080,245-byte part. The independent audit verifies strict membership, reconstructed hashes, paired v28 potential identity, executable provenance and exact receipt replay. The machine report is data/audit/rya1354/source_curation_v2/fac_radial_preflight_v32.json.

An ad-hoc metadata inspection initially selected an older root module via Python stdin import precedence; explicit snapshot sys.path ordering corrected it. No native calculation changed.

## Resume point and limits

The proposed 114-group / 6,600-level 4s/4p valence CI has **not started**. Next prepare its complete native evidence/LS export strategy and memory bounds on Sirius, then declare and run that controlled extension. Existing gauge and potential-policy sensitivity remains unresolved. This is neither full-n4 convergence nor permission to assign line grades.

Bibliography audit still stops before traversal on pre-existing civis2013 and mallinson2024 verified-state vocabulary errors; no clean-library verdict. The unrelated full bug-ledger schema debt remains unchanged. This draft checkpoint uses the previously documented push exception for unrelated merged SEQUENCE debt, with manual file-size checks; no merge or publication.
