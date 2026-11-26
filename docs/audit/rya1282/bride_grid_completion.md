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

The isolated retries recovered all 84 remaining nodes. The ordinary
solver defaults are retained. The completed audit found 3,929 expected
models and 43,219 solved nodes, with zero missing or failed nodes. The
835,164,956-byte TurboSpectrum binary passed full readback against the
source arrays. The check verified all per-record model IDs, aux metadata,
byte offsets, depth grids, and departure coefficients.

A held-out abundance interpolation check on the completed grid tested
35,361 interiors with converged adjacent abundance records.
The 99th percentile of the maximum four-level photospheric error was
0.00665 dex for plane-parallel inputs and 0.00797 dex for spherical inputs.
The corresponding all-level 99th percentiles were 0.01024 and 0.04680 dex.
This measures only the abundance direction on a 0.4 dex bracket.

## Independent evidence and scientific status

The completed grid compared 3,672 matching model/abundance cases
against the held GALAH Balder oxygen grid. Selection uses plane-parallel
t01 and one-solar-mass spherical t02 atmospheres, compares the four
identified triplet levels on overlapping log tau500 −3 to 0, and refuses
abundance brackets that cross failed nodes. O41f and GALAH differ in atom,
background, EOS and transfer details. [Amarsi et al. 2020](https://doi.org/10.1051/0004-6361/202038650)

| MARCS atmosphere | Cases | Median case maximum, dex | 90th percentile, dex | Worst, dex |
| --- | ---: | ---: | ---: | ---: |
| Plane-parallel | 1,170 | 0.035 | 0.219 | 0.517 |
| Spherical | 2,502 | 0.077 | 0.435 | 0.946 |

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

At the worst GALAH atmosphere, the 7-ray control at the nearest exact
Bride node, A(O) 9.7, changed maximum triplet log b by 0.00275 dex.
The corresponding Gerber control uses the same O41f atom and gives
maximum Bride differences of 0.0360 dex at A(O) 9.5 and 0.0391 dex
at A(O) 9.7. The two solar Gerber controls agree within 0.0021 dex.
The two reference comparisons use different atoms and background physics.
The angular test rules out the original angular sampling as the primary
source of the hot giant discrepancy in this control.

The catalogue audit identified 49 unmatched rows whose estimated lower
and upper energies are both below the O I ground ionization threshold.
Another 30 unmatched rows have at least one endpoint above that threshold
by the catalogue-energy estimate. Resolving the former requires an explicit
level/transition crossmatch. The latter require separate classification
before extending the statistical-equilibrium atom.

An atomic-blend transfer control replayed the solar O53 populations with
the held VALD LTE blend opacity and toggled O line opacity. Peak
oxygen-attributable flux differences were 0.25%, 0.21%, 0.20% and 0.39%
in wavelength order. Differential areas were 0.347, 0.716, 1.201 and
1.915 mÅ. This isolates oxygen under the same assumed atomic background.
Molecular blends and observed instrumental and telluric transfer are
still required for an observational qualification.

In the held solar STAGGER run, 255 NLTE columns have paired LTE spectra
on the same 256-column lattice. Equal-area intensity averages yield
777 nm triplet NLTE/LTE equivalent-width ratios of 1.120, 1.112 and
1.100 at disk center. For 128-column draws from the same lattice,
the 95th-percentile sampling shift of the strongest line ratio is
0.00367 dex. Recorded wall times were 6,399 s NLTE and 7,152 s LTE.
These measurements describe one snapshot at disk center and leave
time, angle, depth and frequency representativeness unqualified.

## Interpolation and experimental atom controls

The completed-grid atmosphere holdout varied one axis at a time while keeping
geometry, mass, microturbulence, and A(O) − [Fe/H] fixed. For the four
triplet levels, the 99th-percentile maximum photospheric errors were:

| Axis | Plane-parallel, dex | Spherical, dex |
| --- | ---: | ---: |
| [Fe/H] | 0.00863 | 0.02074 |
| log g | 0.00752 | 0.03972 |
| Effective temperature | 0.02705 | 0.06156 |

The largest temperature holdout error was 0.08536 dex for the 4250 K,
log g −0.5, solar metallicity, 15-solar-mass spherical input at A(O) 9.6.
An eight-corner test varying temperature, gravity and metallicity together
covered 5,632 plane-parallel and 5,907 spherical interiors at matched
abundance offsets. Its 99th-percentile triplet errors were 0.02537 and
0.04568 dex, respectively.

Three NIST [O I] branches at 5577, 6363 and 6391 Å were added between
existing O53 levels, increasing the visible catalogue match count to
20/69. Replaying the original solar populations produced differential
flux areas of 1.4241, 1.1103 and 0.000526 mÅ, respectively. The 6300 Å
and 777 nm control windows had zero change at fixed populations.
Of the remaining 76 unmatched catalogue rows, 30 have an estimated
endpoint above the O I ground ionization threshold, four bound
transitions match a missing 6s ³S level, and 42 bound rows still lack a
secure identity in the held transition table.

A separate retry of the missing STAGGER column (12, 77) passed the
observable stop after 50 population steps with the current solver.
Its triplet equivalent widths were 80.61, 69.89 and 55.16 mÅ.
The aggregate above uses the original paired run for consistent inputs
and solver revision.

At 7000 K, log g 2.0, [Fe/H] +0.5 and A(O) 9.20, the O53 solution
converged in 33 cycles after initializing the shared states from the
O41f abundance bracket for the same star and its added states from the
6500 K O53 solution. Both Ng accelerators were disabled. The unchanged
convergence checks returned population change 6.46e-7, radiation
change 8.35e-7, independent rate residual 1.16e-7 and particle error
9.99e-16. The first two cold starts exceeded their time budgets.
The extended atom with the added forbidden branches also solved
self-consistently for the solar control in nine population steps.

Fixed-population spectra passed a six-band smoke test on the Sun,
the 4500 K metal-poor dwarf, and the 6500 K giant. The line-on and
line-off spectra were finite and positive at every sample. The 95
matched component windows per atmosphere were split across near-UV
(6), visible (32), red (25), Y/J (3), H (16), and K (13). Every window
responded to the oxygen opacity toggle. The canonical list contains no
O I rows in the 1.4–1.5 µm band gap.

The radiative O54 variant adds the 6s ³S₁ level from Li, 18 E1
links to existing states, and its ground-target photoionization
continuum from NORAD. It preserves all O53 collisional rates. The
catalogue recheck adds two visible and two K-band matches, reaching
22/69 and 10/13 respectively. The new level has no dedicated
electron or hydrogen collision rates and is a radiative sensitivity
experiment. The solar model at A(O) 8.70 converged in nine population
steps and produced finite flux in all four diagnostic windows.

The compiled TurboSpectrum interpolator passed the same solar and
6500 K giant readback on a small fixture and the full Bride binary. At the exact A(O) 8.60
node it reproduced the departure arrays within 2.2e-6 dex after
text export. At the requested A(O) 8.67 it stamped that abundance on
the output and selected the A(O) 8.60 coefficients for both inputs.
The difference from linear abundance interpolation was 0.0105 dex
for the Sun and 0.0398 dex for the giant across all levels and depths.
Off-node synthesis through this executable needs an explicit abundance
interpolation path or a validated nearest-node error budget.

Fresh Bride solves at A(O) 8.67 tested the off-node behavior directly.
For the Sun, the maximum photospheric triplet departure error was
0.00734 dex using the nearest node and 0.000208 dex using a linear
interpolation of departure coefficients. For the giant, the errors
were 0.02327 and 0.000853 dex, respectively. The two controls favor
linear interpolation for arbitrary abundances at exact atmosphere
nodes in this grid.

One O41f hot giant at 7000 K, log g 2.0, [Fe/H] −1.0, 0.5 solar
masses, and A(O) 8.20 still hit the inner cap tolerance after 2,000
radiation iterations from a cold start. Seeding the populations from
the validated A(O) 8.0 and 8.4 nodes and disabling both Ng accelerators
gave strict convergence in 25 cycles with rate residual 1.64e-7 and
conservation error 6.66e-16. The seed sets the starting populations.
The checkpoint stores the independently solved array.

## Experimental off-node coefficient export

The departure-grid CLI's `export-departures` command takes an exact
MARCS model ID, a trial A(O), and a staged Bride NPZ directory. It
interpolates the 41 O41f departure coefficients linearly between
adjacent converged A(O) nodes and writes a TurboSpectrum ASCII
coefficient file with the trial abundance and native depth scale.
The command rejects extrapolation and failed brackets.

The TurboSpectrum Fortran `read_departure` routine read the solar,
4500 K dwarf, and 6500 K giant exports with maximum departure
differences below 2.7e-8 dex at its single-precision read boundary.
A new 4500 K dwarf solve at A(O) 7.67 independently reproduced the
interpolated photospheric triplet departures to 2.36e-5 dex.
