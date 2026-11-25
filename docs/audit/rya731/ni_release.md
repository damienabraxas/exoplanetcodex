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
