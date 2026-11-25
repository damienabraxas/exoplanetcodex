# Experimental O53 radiative extension

The targeted extension adds 12 fine-structure levels to O41f and 143 electric-dipole
links. It uses NIST energies reported in Li et al. Table A.1 and the Babushkin
oscillator strengths from Table A.2. The original 256 line records and all
original rate coefficients are retained. The two O II indices move from 40–41
to 52–53. The current grid uses the pinned O41f atom. [Li et al. 2023](https://doi.org/10.1051/0004-6361/202245645)

| Air wavelength group | Terms | Resolved components |
| --- | --- | ---: |
| 369.2 nm | 3s ³S° → 5p ³P | 3 |
| 1588.8–1589.3 nm | 4p ⁵P → 5d ⁵D° | 9 |
| 1686.9–1687.5 nm | 4p ⁵P → 6s ⁵S° | 3 |
| 1745.3–1745.4 nm | 4p ³P → 5d ³D° | 6 |

The component counts include unresolved blends represented by single catalog
rows. The 5d ³D upper energies in Table A.1 are identical to the published
precision. The builder preserves J identities through the degeneracy. It logs
147 Li links touching the four new terms whose other endpoint is unrepresented
or collapsed in O41f.

The second radiative candidate adds 12 continua to the O II ground target using
four identified NORAD LS partial cross sections. Each J component of a term
uses the same cross section. The photon grid is aligned to each observed
threshold while preserving excess photon energy. The source tables are kept
intact. The solver input uses the O41f shortwave limit of 90.0064 nm because
the unrestricted extension enters wavelengths near 4.5 nm and produces
underflow in the current photospheric background and radiation acceleration.
The unrestricted LTE control converged in 44 formal solutions. Trimming to the
baseline domain took 10 and changed the tested intensity ratios by less than
3e-15. The receipt records the omitted point counts per channel. [NORAD O I data](https://norad.astronomy.osu.edu/), [Nahar 1999 partial cross sections](https://www.astronomy.ohio-state.edu/nahar.1/papers/apjs-pirc-oions.pdf)

The paired 5750 K, log g 4.5 solar-metallicity MARCS run with fixed LTE
populations converged its radiation in 8 steps with O41f and 10 with O53 plus
the ground-target continua. It returned finite, positive disk-center
intensities in all four added wavelength groups. Across the 777 nm control
window, the greatest relative change was 4.95e-8. The four minimum
O53/O41f intensity ratios were 0.999404, 0.989284, 0.997301 and 0.994623 in
wavelength order.

A bounded NLTE diagnostic on the same MARCS atmosphere converged in nine
population steps in 75.4 s. The 53-by-56 departure array is finite and positive
with a full-depth range of 0.660–9.948. This checks numerical execution with
the incomplete collision network. Scientific line comparisons are pending.

The candidate remains on scientific hold for grids. The Barklem O+H source
table and the precomputed oxygen Kaulakys matrices do not contain the four new
terms. The trial `kaulakys_O53_input.tsv` in the artifact directory appends
their statistical-weighted term centroids and effective principal quantum
numbers from observed binding energies to Barklem's original oxygen deck
(commit `dafeac77e189f468d39dbe43ad6cbf442a6b404b`). Running that deck requires
an IDL-compatible runtime and Barklem's `mswavef` and `integ_rate` libraries.
Amarsi et al.'s expanded electron-impact table stops near 12.65 eV, below
all four additions. Electron collisions, pressure broadening, excited-ion
channels and sensitivity to the
excluded radiative links also remain open. The new line calculations used
natural damping from the Li lifetimes and an experimental quadrature of 31
samples with MULTI QMAX 15. [Barklem’s Kaulakys implementation](https://github.com/barklem/kaulakys), [Amarsi et al. 2018](https://arxiv.org/html/1803.10531v2)

Rebuild with `scripts/rya1282_extend_oxygen_atom.py`. Use the held O41f file
for `--base` and the Li A.1/A.2 source directory for `--li-dir`. Add
`--norad-dir` with the original `o1.en.ls.txt` and `o1.ptpx.txt` to
include the ground-target continua. The Mac artifact is under
`bride_mac_checks/rya1282-o53-radiative-v4`. The paired LTE driver, arrays,
summary and source checksums are retained with the candidate.

## Collision sensitivity on the solar MARCS case

The headless GDL runtime executes Barklem’s published impulse calculation
with the pinned momentum-space wavefunction and rate-integration libraries.
The experimental O I batch covers 54 new same-core LS term pairs and one
existing control from 2000 to 20000 K. The scattering-length approximation
used 100 collision-energy samples. A portable trapezoidal integration over
the solver's dense log-momentum grid substitutes for IDL's INT_TABULATED in
that approximation. The full method uses Markwardt QPINT1D directly.
The GDL input, per-pair source quantum numbers, programs, logs and output
rates are retained in `rya1282-o53-radiative-v4`.
[Barklem implementation](https://github.com/barklem/kaulakys),
[GDL release](https://github.com/gnudatalanguage/gdl/releases/tag/weekly-release),
[QPINT1D](https://www.nv5geospatialsoftware.com/docs/qpint1d.html)

At 5000 K, the approximate-to-full rate ratios for the four target transitions,
in wavelength order, were 0.269, 0.521, 0.418 and 0.558. The full comparison
used 30 energy samples. The existing 4p ⁵P → 4p ³P upward control was
1.358e-9 cm³/s with the local full calculation versus 1.208e-9 in the published
100-point output. The approximate calculation gave 3.814e-9 versus 3.226e-9 in
the corresponding published matrix. The full 100-point term-rate calculation
is queued with paired checks on the Sun, a 4500 K metal-poor dwarf and a 6500 K
spherical giant.

The radiative candidate was compared with three experimental perturbations.
The hydrogen-only atom uses statistical final-J excitation branching,
reconstructed downward CH coefficients, and 362 fine-structure pairs.
The electron-only atom adds 57 same-spin E1 CE coefficients using the
neutral-target prescription, with no arbitrary rates for forbidden or
spin-changing processes. The combined case uses both additions. [Van Regemorter 1962, Eq. 18 and Table 2](https://articles.adsabs.harvard.edu/pdf/1962ApJ...136..906V)

| Solar case | Population steps | Maximum photospheric shift in log b, dex | Differential flux area near 1.589 µm, mÅ |
| --- | ---: | ---: | ---: |
| Hydrogen approximation | 9 | 0.0758 | −1.084 |
| Empirical electron excitation | 9 | 0.0282 | −1.058 |
| Combined | 9 | 0.0692 | −1.310 |

The flux areas integrate one minus the variant-to-radiative-reference flux
ratio in a shared window. These calculations test sensitivity at fixed
atmosphere and abundance. The observational comparison requires the band's
blends and telluric mask, and the collision uncertainty and pressure
broadening remain open.

## Full-rate run recovery

On 2026-11-24 at 07:18:44 MST, the macOS kernel logged that it killed GDL
PID 52851 under memory pressure. Its logged compressed memory was 90,923 MB.
The original batch retained 33 complete rows. The resumed run uses the
same full 100-point calculation in a fresh GDL process for each remaining
pair. A result is checkpointed only after the worker exits successfully,
prints its completion marker, and supplies a finite, nonnegative row on the
expected temperature axis. The merge checks every pair against the original
manifest and keeps the original rows. The detached workflow starts the
three-atmosphere Bride validation after the complete table is written.
The runnable recovery and its log are in the experiment directory as
`resume_full_rates.py` and `full_rate_resume.log`.

The first missing pair, 4p ⁵P → 5p ³P, completed in isolation in 130.4 s.
It wrote a fresh 11-temperature checkpoint with an upward rate of
7.444e-13 cm³/s at 5000 K. An in-flight RSS sample was 1.13 GiB.
The saved 100-point control from the original run was 1.360e-9 cm³/s
at 5000 K, compared with 1.358e-9 from the 30-point test.

The 4500 K MARCS dwarf at [Fe/H] −1 converged in eight population steps
with radiative O53 and six with the combined approximate collisions.
The largest shift in photospheric log b was 0.121 dex. The largest of the
five differential flux areas was 0.0311 mÅ in magnitude at 1.589 µm.
The four automated oxygen-extension tests passed on rerun.

The held RYA-1053 VALD source supplies zero damping entries on all 13
H-band O I rows recovered within the three target groups. The near-UV
pull supplies one 3692.395 Å component with a van der Waals entry of
−7.931. A separate fixed-population solar test applied classical Unsöld
broadening to all 18 resolved H-band components. It changed the three
differential flux areas by −0.0490, −0.0240 and −0.0669 mÅ, respectively.
This measures sensitivity to the width prescription under the same
atmosphere and populations. It leaves transition-specific broadening
qualification open.
