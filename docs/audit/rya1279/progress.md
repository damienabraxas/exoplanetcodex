# Bride continuation: RYA-1279 / RYA-1280

Base: 4635c7421a7dee8671ffefe05d771c68da1cd24e (unmerged parent PR #597).
Worktree: codex/rya1279. Branch: ryandamienschmitt/rya-1279-1280-bride-oxygen-gates.
Ryan authorized continuing through the physics gates toward oxygen results and
subsequently said to proceed without waiting for CI. That supersedes an additional
permission checkpoint for the bounded experiment; scientific validation still applies.
No merge/publication authorization was given.

## Verified work

- Held Gerber `atom.o41f` receipt and MD5 match. This candidate has 41 levels,
  256 lines, 22 continua and 1106 collision tables. All triplet components exist.
  No level/rate reduction was made. Full source remains on Sirius; local small-file
  inspection is evidence only, never a production compute-input fallback.
- Found and fixed a converter defect: only the first TEMP grid was retained;
  976 later rate sets were truncated/misassigned. Each collision now owns its active
  grid: CE 130x7 + 611x10, CH 302x10, CI 39x10, CH0 24x10. All source rates and
  their own grids compare exactly after the documented cm3-to-m3 conversion.
- Found and fixed discarded damping coefficients: triplet packed ABO 453.234 had
  become Unsold 1. Native conversion now keeps sigma=453, alpha=.234, radiative
  gamma=3.89e7. Supported subset is zero/standard Unsold/packed ABO; unvalidated
  negative encodings, nonzero Stark and unsupported profiles/rates fail closed.
- Parser rejects unknown/truncated/extra rate data, bad grids, missing END and
  fixed radiative rates rather than skipping them. No existing Fe result is rerun.
- Synthetic actual-backend checks: continuum slab/thin limit with explicit
  boundaries, 3 formal solvers, depth/ray cases; actual CE/CI+SE vs independently
  computed Saha-Boltzmann with matched constant vintage and Debye lowering;
  alternate initial populations and rate scales; coupled CRD line+SE thermal
  enclosure. Deliberate wrong directional rates converge/conserve but fail the
  independent population test. Physical R+C diagnostics are distinct from MALI
  Gamma, which must not be mistaken for the physical SE operator.

## Numerical experiments retained honestly

Initial roundoff-only slab criterion failed four cases (largest absolute error
1.17e-7 for cubic Bezier). Inspection identified deliberate polynomial branches,
not a source/grid-tuning opportunity. Subsequent acceptance bounds are derived
from alternating-series remainders and rational conditioning in inspected 0.16.1
`LwInternal.hpp`, `FormalScalar.cpp`, `Bezier.hpp`. The tests do not choose a
constant from the observed discrepancies. These are bounded homogeneous controls,
not a universal spectral/abundance error budget.

The first LTE comparison used modern CODATA and omitted default Debye lowering;
that was a mismatched test reference. Matching the declared physics/constants is
required for an independent comparison. The corrected analytic reference shares
only constants with LW, never its LTE/population implementation. A hot synthetic
ion limit has a conservation-augmented condition number about 1.39e8; componentwise
roundoff assertions are inappropriate. Normwise backward error and an explicit
64*eps*kappa forward bound are tested, while component residuals remain reported.
Equation row equilibration preserves the homogeneous system; it changes no rates.
The coupled radiative control needed 51–54 iterations, exceeding an initial 20
iteration diagnostic budget. The bounded cap is now 80; population and J changes
must both reach 1e-11, with independently tested 1e-9 thermal-population accuracy.
No tolerance or physical input was changed toward an abundance target.

Exact vacuum (zero chi, undefined eta/chi source) is not labelled a passed backend
case. Positive-opacity thin-limit cases reach 1e-10 optical depth. Full-3D coupling,
nonuniform formal-transfer accuracy and real-atom resolution convergence are still
separate work; these controls do not certify a full-3D solver or a solar abundance.

## Provenance and baseline still being resolved

The distribution's literature attribution points to Gerber2023 and Bergemann2021.
[Bergemann2021 §3.1/Table1](https://arxiv.org/html/2109.01143v1) describes 122 levels,
not the held 41, and documents ABO and two alternative gf scales. No equivalence
of the reduced artifact to that paper or to Amarsi2018 is asserted. The historical
Gerber gateway anchor is not an identical-input independent numerical benchmark.
The raw source is not redistributed; TSFitPy's GPL does not supply a data license.

[Nahar2022](https://doi.org/10.1093/mnrasl/slac015) disputes photoionization data shown
in Bergemann2021. This is a primary methodological warning, not proof that the held
41-level artifact contains those disputed data. It must be traced rather than
silently inheriting or rejecting a model by author label. No cross-section altered.

Next: extract the already-held independent MULTI departure reference at the exact
MARCS benchmark node; freeze its atmosphere/atom/abundance/opacity assumptions and
compare, then run the bounded STAGGER O experiment if applicable scientific gates
are satisfied. MARCS 5750/4.5 is a benchmark node, not the canonical Sun's stellar
parameters and not a substitute for the STAGGER cube.

## Tests / operational evidence

Held Bride interpreter, from checkout:

```
PYTHONPATH=. /Users/ryanschmitt/bride_rya1118/.venv/bin/python -m pytest -q -p no:cacheprovider tests/test_multi_atom_rya1279.py tests/test_bride_numerics_rya1280.py tests/test_multi_atom_rya1119.py tests/test_multi3d_read_rya1118.py tests/test_bride_foundation_rya1278.py --tb=short
```

138 passed in 2.94s, exit 0. LW emitted 228 known stat_equil deprecation warnings.
Project interpreter (without LW), new parser/numerical plus reader/foundation tests:
40 passed, 7 skipped, exit 0. No full-suite/CI claim. Whitespace check exit 0.

A dedicated Sirius environment is being installed under the ticket work root using
LW0.16.1/numpy2.5.2/scipy1.18.1/numba0.67.0. No pre-existing environment modified.
Host load was 10.51/11.82/10.09 with 21 users; no heavy science job started. Future
jobs must be bounded, low priority, single-threaded initially, and separately logged.
Original FITS, frozen gold, Fe/Al products and published values remain untouched.


Further checkpoint: the explicit source Q0=3 for the triplet is now preserved
instead of guessing qCore=qmax/4. Q0=0 retains a documented automatic numerical
prescription; native frequency-grid equivalence is not asserted. The O grid has
an abundance axis; reference intake requests the full atom's stored abundance
exactly and refuses nearest-node/nearest-abundance substitution. The new reader
checks binary record identity/dimensions, uses (level,depth) ordering and checks
MARCS row identities and SI conversions on independent non-square synthetic data.
Source values are never adjusted to improve a comparison.

Sirius checkpoint: isolated environment installed successfully; pinned LW0.16.1
synthetic tests at 4b3da23f passed 95/95 in 38.34s, with the same 228 deprecation
warnings. Exact MULTI reference intake at 00814f76 succeeded: model 5750/4.5,
A(O)=8.70, 56 depths, 41 levels; max logtau500 mismatch to printed MARCS structure
1.89e-8 (Rosseland differs by 0.226 dex, so those axes are not interchangeable).
No solver had run at this checkpoint. Existing STAGGER prov.json was located;
its acquisition identity is useful, but its raw-EOS/byte-count claims are older
than the corrected reader/header audit and cannot override that evidence.

The bounded MARCS model-under-test runner now prepares a manifest by default.
It validates all source hashes, uses all MARCS header abundances with the exact
reference trace O abundance, explicit source gravity/units/boundaries, native
passive atoms, fixed ne, no molecular opacity, one thread, physical R+C diagnostics,
a 900-second budget, population/J/rate convergence 1e-6 and particle conservation
1e-10. These are declared numerical experiment controls, not abundance-admission
thresholds. It saves independent-reference departures and absolute triplet
intensities in W m^-2 sr^-1 Hz^-1, plus an isolated LTE comparator. It must label
background/EOS/model differences honestly; a converged experiment is not a passed
paper replication or solar-abundance product. The user has authorized continuing
these validation experiments; no new merge/publish permission is inferred.

## First real-input attempt and resolved format difference

At 85820b0b the prepared exact-node manifest passed source-hash checks but the
first run failed during wavelength-grid construction, before RT/SE. Source
continuum 4→41 has an energy-derived edge 82.7947 nm while its explicit table
runs 94.3253–106.4179 nm. Standard LW ExplicitContinuum clipped all its points
and raised IndexError. The original failed outputs/log remain intact under
Sirius ticket root `benchmark-q1-rays3-v1`; no fabricated spectrum was emitted.

The original MULTI2.3 source was located in the public pymulti distribution at
commit `4d6aa3c39d634b38da772e469a5692435eba1354`,
`pymulti/files/multi/mul23/source/mul23_sub.f`. ATOM 540–572 and FREQC 2325–2376
explicitly consume the supplied absolute wavelength table for QMAX<0, without
filtering it to the energy-derived threshold. The radiative Milne factor at
1010–1019 uses the actual photon frequency and LTE population ratio. Thus a
faithful legacy-model conversion must preserve the explicit table separately
from energy normalization, not truncate/relabel/shift its source data.

`MultiExplicitContinuum` implements that named convention, retains the energy
lambda0 and exposes both energy/table edges for audit. All 22 source continua
and 3778 base frequency points now construct. This is NOT a physical validation
of the source atom or an explanation of the level/table discrepancy; that caveat
remains. The converter does not manufacture cross-sections or alter level/rate
values. A synthetic conflicting-edge fixture verifies retention and zero opacity
outside the source table. Dry-run now checks native wavelength construction and
cross-sections before PREPARED, and records solver_started separately from setup.

Original-source permalink:
https://github.com/MingjieJian/pymulti/blob/4d6aa3c39d634b38da772e469a5692435eba1354/pymulti/files/multi/mul23/source/mul23_sub.f
