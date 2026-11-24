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
