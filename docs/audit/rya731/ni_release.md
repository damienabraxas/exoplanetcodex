# Solar nickel release evidence

The Sirius runs use separate products for each ion, holding, band, model, line set, and continuum treatment. The OpenMULTI3D comparison uses the solar MARCS atmosphere and the same pinned Ni atom as Gerber. The STAGGER calculations use a horizontally averaged atmosphere.

## Line preparation and fit recovery

- The Ni I optical laboratory sweep adds the Blackwell and Huber furnace measurements alongside Wood et al. The calculated infrared oscillator strengths are carried by the separately adopted `NI_ELGUETA2026` set.
- The `NI_CLEAR2026` set applies Clear et al.'s calculated oscillator strengths to the 26 matched Ni II transitions. The override is local to the set and preserves the blend opacity in the surrounding window.
- The GES visible sidecar maps 15 missing upper endpoints to atom level 108. Its statistical weight of 35 and energy of 49,351.829 cm⁻¹ reproduce the weighted centroid of all five components of the NIST Ni I 3d⁸4s(⁴F)5s ⁵F term. The target J=2 component is at 50,346.427 cm⁻¹.
- The near-UV sidecar adds ten transitions through the corresponding 4d superlevel centroids. Each pair below reproduces the pinned atom's computed energy within 0.002 cm⁻¹. The endpoints use wavelength, lower excitation, and both J values so coincident lines keep distinct identities.

| Atom level | NIST component energies (cm⁻¹) | Component J |
|---|---:|---:|
| 160 | 56,624.599 and 56,885.215 | 6 and 7 |
| 162 | 56,766.480 and 56,973.654 | 4 and 5 |
| 163 | 56,801.586 and 56,954.098 | 5 and 6 |
| 178 | 57,700.590 and 57,767.467 | 2 and 1 |
| 179 | 57,789.532 and 57,829.332 | 4 and 5 |

The level data come from the [NIST Atomic Spectra Database](https://physics.nist.gov/cgi-bin/ASD/energy1.pl?spectrum=Ni+I&conf_out=1&term_out=1&level_out=1&j_out=1). The Ni atom pin and stage coverage are in [atom.ni538qm.json](../../../data/reference/nlte_atom_reach/atom.ni538qm.json).

Five Kitt Peak visible profile fits failed during synthesis sensitivity probes because scratch files disappeared. An isolated replay used the original EWs and recovered 6327.5712, 6378.2254, 6532.8474, and 6586.2836 Å. The fit of 6767.7806 Å converged and remains excluded at reduced χ² = 10.222. The replay fits are diagnostic evidence.

## Cross-checks and coverage limits

- Both Kitt Peak reductions yield near-UV Ni II on the Clear set. 1D LTE gives 6.392 and 6.557 dex, while mean-3D LTE gives 6.483 and 6.679 dex with 20 retained lines in each reduction. The atom has no Ni II bound-bound transitions. The visible Ni II fits retain no lines.
- CRIRES+ H-band Ni I retains seven LTE lines and five Gerber NLTE lines. Fixed-continuum mean-3D NLTE gives 6.092–6.147 dex across the two reductions, and the free-continuum treatment gives 6.197–6.232 dex. These are separate treatments and the continuum sensitivity is material.
- The Elgueta set yields a two-line wide-Y CRIRES+ measurement and a two-line Kitt Peak near-infrared measurement. The narrow Y and atlas fits that use just one line remain diagnostics. The registered set yields no CRIRES+ J or K aggregate on the run's holdings.

## Feed admission

Ni v1.161 carries 223 admitted Ni I products across seven treatments and five bands, with one additional quarantined row. All admitted records were joined by the publication key and their staged source CSVs verified against the recorded SHA-256. The two Ni II near-UV holdings have viable LTE fits, but the budget assembler could not resolve the model-atmosphere component for the Clear set. The visible Ni II cells contain no retained lines. The release keeps both dispositions in the run and publication reports.

## Ni II near-UV patch (v1.167)

Scott et al. (2015) tabulate no Ni II, so the assembler had no route to the model-atmosphere term. The Gerber and OpenMULTI3D LTE comparands fit the same lines, with the same fitter and window, on MARCS.GES with the departures withheld. The 1D-LTE route uses ATLAS9.Castelli. The paired difference therefore isolates the atmosphere grid. The assembler now prices `model_atmosphere` as the absolute paired median of A(MARCS.GES) − A(ATLAS9) on each product's pool. A peer unit is used only when its ATLAS9 pool is identical.

| Kitt Peak reduction | Lines | MARCS.GES − ATLAS9 | Same-atmosphere null |
|---|---:|---:|---:|
| Kurucz-corrected | 18 | +0.016 | 0.000 |
| molecfit | 21 | +0.020 | 0.000 |

Ni v1.167 adds six Theoretical Codex products: 1D-LTE at 6.557 ± 0.077 and 6.392 ± 0.110 dex, Gerber 1D-LTE at 6.558 and 6.386 dex, and OpenMULTI3D 1D-LTE at 6.558 and 6.386 dex. Line scatter dominates the uncertainty. The two mean-3D LTE products (6.679 and 6.483 dex) stay held: no 1D grid applies to them, and no published Ni II atmosphere term exists. The assembly read the committed nominal artifacts above. Those are md5-identical to the Sirius run units.
