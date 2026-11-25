# Bride oxygen grid recovery and qualification

## O41f MARCS departure grid

The approved run covers 3,929 atmosphere files and 43,219 abundance nodes.
The selection contains 1,419 plane-parallel and 2,510 spherical MARCS input
models. It has 13 temperatures from 4000 to 7000 K, ten metallicities from
−1.5 to +1.0, gravities from −0.5 to 5.0, and the matching tabulated mass and
microturbulence choices. Bride applies one-dimensional column transfer to
these models. The spherical designation describes the input atmosphere.

All 3,929 original atmosphere files arrived. The first run had 103 failing
abundance nodes. The negative-opacity scattering fix recovered 19. Two
controlled reproductions located distinct convergence problems:

- At 5750 K, log g 4.0, [Fe/H] −1.5 and A(O) 6.1, radiation Ng stalled
  even with 2000 inner formal solutions. With radiation Ng disabled, the
  same fixed-point tolerances were met in 12 population steps. The final
  independent rate residual was 6.48e-8 and particle-count error 6.66e-16.
- At 6500 K, log g 4.5, [Fe/H] −0.75 and A(O) 6.85, disabling
  population Ng recovered the previously timed-out calculation in 18
  population steps. The final rate residual was 1.71e-7 and particle-count
  error 6.66e-16.

The strict retries run on the Mac using two isolated workers and their
own per-node checkpoints. The ordinary solver defaults are retained.
The assembled binary is conditional on a complete exact-manifest audit.
The follow-up checks every binary record against the source array and
checks aux offsets, depth grids and total record count.

A held-out abundance interpolation check on the first staged checkpoint
tested 35,259 interiors with converged adjacent abundance records.
The 99th percentile of the maximum four-level photospheric error was
0.00665 dex for plane-parallel inputs and 0.00796 dex for spherical inputs.
This measures only the abundance direction on a 0.4 dex bracket.

## Independent evidence and scientific status

The first staged checkpoint compared 3,666 matching model/abundance cases
against the held GALAH Balder oxygen grid. Selection uses plane-parallel
t01 and one-solar-mass spherical t02 atmospheres, compares the four
identified triplet levels on overlapping log tau500 −3 to 0, and refuses
abundance brackets that cross failed nodes. O41f and GALAH differ in atom,
background, EOS and transfer details. [Amarsi et al. 2020](https://doi.org/10.1051/0004-6361/202038650)

| MARCS atmosphere | Cases | Median case maximum, dex | 90th percentile, dex | Worst, dex |
| --- | ---: | ---: | ---: | ---: |
| Plane-parallel | 1,169 | 0.035 | 0.219 | 0.517 |
| Spherical | 2,497 | 0.077 | 0.435 | 0.946 |

The worst case is 6750 K / log g 1.0 / [Fe/H] +0.5 / [O/Fe] +0.5.
Its largest triplet discrepancy is near log tau500 −2.97 and declines to
0.067 dex near −1. The comparison identifies a high-temperature,
low-gravity physics qualification gap.

For O53, the full Kaulakys run validated 55 rate pairs, including one
existing control. The 54 new rates and the empirical electron rates
converged on the Sun, the 4500 K metal-poor dwarf, and the 6500 K giant in
9, 6 and 16 steps, respectively. The 4000 K metal-poor endpoint
also converged in 11 steps. Solar observational preflight confirms data
coverage at 369.2 nm and all three new H-band groups. The atlas windows
contain stronger neighboring features such as V I at 369.22 nm and
Si I and Fe I near 1.589 µm. The O53 profiles used for the convergence
checks include the active oxygen atom and continuous background.
Quantitative observed comparisons require matching blend opacities,
pressure widths, and wavelength-dependent telluric and instrument
treatment.

Grid admission remains HOLD. Numeric completion and binary verification
cannot clear the giant discrepancy or the observed O53 profile gate.
G5 convergence across representative snapshots and G6 horizontal transfer
remain separate acceptance gates.

Work artifacts are under
`~/codex/bride_mac_checks/rya1282-grid-completion`. The scripts and logs
there retain the retry strategy, exact input model list, referee extract,
observed preflight, and assembly checks.

The catalogue audit matched O53 line centres to canonical O I rows using
air wavelength and lower excitation energy. It matched 4/4 near-UV and
9/9 H-band rows, and 17/69 visible, 25/47 red, 3/3 Y/J and 8/13 K-band
rows under the stated tolerances. The unmatched rows include lines from
highly excited states outside the O53 representation. The reported
matches establish catalogue coverage, while unresolved rows may also
include unresolved components or convention differences. The machine
audit lists every unmatched line ID for further classification.
