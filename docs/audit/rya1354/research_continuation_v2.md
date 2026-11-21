# RYA-1354 — two-method correlation and primary-source validation

**Scientific HOLD; exposed, exploratory nonblind work.** This continuation starts
at `f759d2c0257599c6558febccb2ff186916240f17` on the existing RYA-1354 branch and
draft PR 616. It resolves the failed n=5 initializer, extends the independent
method through the same valence ladder, and reaches a primary experimental
lifetime comparison. It does not establish calibrated uncertainties or adoption.

## Two-method valence sequence

The original Thomas–Fermi n=5 initializer failed again in a controlled replay.
Hydrogenic initialization completed on 220-point, extended 260-point and
half-spacing 440-point grids. These are **initialization controls**, not SCF-grid
convergence tests. The original 220-point grid already extends to about 1892 Bohr;
the original extent warnings did not justify assuming that the box was too short.
The successful n=5 GRASP run changes the initializer, retains the previous n=4
radial functions, and optimizes only the new layer. All failed attempts survive.

Both methods fix 1s2/2s2 occupations and n<3 radial functions, use valence SD, and
optimize new layers through spd/spdf/spdfg at n=3/4/5. Native term/J/parity labels
identify the target 3s 5S(o) and 3p 5P states. MCHF compares an LS multiplet with
the **sum** of GRASP fine-structure gf, never its term-average A with individual
component rates. Different Hamiltonians, optimized radials and native Rydberg
conventions remain explicit. Both methods share development lineage and omissions.

| n | GRASP summed gf B | MCHF gf length | Relative difference | GRASP A gauge differences |
|---|---:|---:|---:|---:|
| 3 | 5.87315 | 5.87023 | +0.0497% | 40.42% |
| 4 | 5.08822 | 5.08589 | +0.0458% | 7.73–7.74% |
| 5 | 5.03566 | 5.03338 | +0.0453% | 5.09–5.11% |

The last gf step is still about 1%. GRASP n=5 native separations are
12698.90/12701.66/12705.83 cm^-1 for J=1/2/3, farther from the observed separations
than n=4. Neither the shrinking gauge discrepancy nor the small method difference
is an error bar. Shared omitted core correlation is demonstrably significant.
The native stopping rule accepts either orbital self-consistency or sufficient
relative EOL-energy stationarity. `grasp_SCF_diagnostics` records both the actual
iteration count and final printed orbital residuals; a native stop does not certify
all orbitals reached ACCY. Separate radial/SCF precision qualification remains owed.

The modern serial MCHF/nonh builds use the pinned ATSP2K source and existing
libraries. Their native Si 1D n=4 control reproduces the upstream eigenvector
file **byte-for-byte**, energy -288.906186807 Hartree. Build/compiler/COMMON warnings
are retained. No upstream source was patched.

GENCL's explicit SD replacements from both open shells are used; its different
“type 2” occupation-restricted complex is not substituted for SD. The first oxygen
transition attempt reports inconsistent 3s/3p ordering, exits zero and produces no
required `.ls` file. Its raw receipt says completed, but that conclusion is rejected
in the derived summary. The driver now rejects this refusal and requires the file.
Adding an explicitly empty 3s to the upper **generator reference** establishes
consistent ordering, then regenerates native couplings and solves MCHF afresh.
No coefficient or term/seniority label is rearranged by hand. Occupation multisets
including CSF multiplicities are exactly equal before/after for both n=3/4 states.

The native MCHF printed air wavelength is also rejected: pinned `probab.f` uses
206030 in the dispersion numerator, giving a roughly 1.35 A discrepancy with the
project wavelength SSOT near 7887 A vacuum. Native vacuum energy/gf remain the
comparison axes. No solver patch or canonical wavelength change hides this issue.

## Core and Hamiltonian sensitivities

The versioned `sensitivity_plan.v2.json` and native driver define four fixed-orbital
n=3 controls. RCI interpolates onto its own default grid (RNT=2e-6/Z, H=.05,
N=590), so its same-CSF Dirac–Coulomb calculation first checks numerical agreement
with MCDHF. Printed energies differ by at most 0.01 cm^-1 and summed gf is identical.
This is a useful control, not a general radial-grid-convergence demonstration.

| RCI change from same-grid DC | Relative summed-gf change |
|---|---:|
| Transverse interaction, frequency scale 1e-6 (static limit) | -0.00647% |
| Vacuum polarization and self energy through n=3 | -0.00868% |
| One 2s hole allowed, fixed n=3 orbitals, Dirac–Coulomb | -2.4118% |

The core expansion shifts the three separations by approximately -185 cm^-1.
These coarse sensitivities are neither independent Gaussian errors nor terms to
combine blindly into sigma. Core-valence/n/l convergence and 1s sensitivity remain
owed before an accuracy claim. Recoil was disabled; no isotope precision is claimed.

Stationary-nucleus RCI emits divide-by-zero, underflow and denormal flags. Source
inspection finds unconditional `1/EMN` at pinned `setham_gg.f90:110`, even with
EMN=0 and both recoil terms disabled; the later recoil uses are guarded. A separate
**synthetic** positive-mass control, with recoil still disabled, removes the
divide-by-zero flag and reproduces the entire DC E1 table byte-for-byte. It is not
an oxygen isotope candidate. Underflow/denormal flags remain reported. The upstream
code was not patched and the original stationary-nucleus artifacts are preserved.

## Forbidden-line convergence

**The initial table below is unqualified for convergence:** n=4 hit its 100-cycle
SCF limit despite exit zero; n=5 froze those unqualified n=4 functions. A strict
semantic refusal check now stops that path. A same-model 500-cycle continuation
also repeats the cycle and fails explicitly. A -0.5 damping continuation also exhausts 500 cycles. Its earlier missing-prompt
attempt fails before SCF and is retained. A -0.9 damping-only continuation reaches
native stopping at n=4 (44 iterations) and n=5 (129 iterations). No physical
parameter, state weight or convergence tolerance is fitted to the external rates.

The five 2p4 ground levels now have a separate common-orbital EOL valence-SD
sequence through n=3/4/5. Core occupations and n<3 radial functions are fixed;
new layers use hydrogenic initialization and are optimized without observed energies.
Leading LS/J/parity identities are checked. Native M1/E2/M3/E4 cover all 16
angularly permitted one-photon channels among these five levels. Every channel's
two identical f1/f2 directional records are checked and retained once, never summed.

| n | Model 1D2 lifetime B / C (s) | Model 1S0 lifetime B / C (s) |
|---|---:|---:|
| 3 | 77.9715 / 77.9215 | 0.297515 / 4.34364 |
| 4 | 80.7780 / 80.7941 | 0.298556 / 0.415873 |
| 5 | 81.8283 / 81.8494 | 0.307748 / 0.392920 |

The explicit resumed sequence is separately recorded in
`forbidden_convergence_v2.json`; the original table and failures stay intact.

| n | Resumed model 1D2 B / C (s) | Resumed model 1S0 B / C (s) |
|---|---:|---:|
| 3 (original) | 77.9715 / 77.9215 | 0.297515 / 4.34364 |
| 4 | 80.7516 / 80.7673 | 0.298581 / 0.415841 |
| 5 | 81.8410 / 81.8624 | 0.307742 / 0.392327 |

**Native stopping is not full orbital convergence.** Pinned `scf.f90` accepts
either its function criterion or a sufficiently small relative EOL-energy change.
Final maximum printed orbital residuals are 2.42e-6, 3.37e-5 and 2.61e-5,
respectively, versus native ACCY=5.96046448e-8. All three exceed that criterion;
`ladder_qualified_for_convergence` remains false. The new summary exposes the
residual tables and stopping rule rather than promoting execution completion to
precision convergence. Resumed n=5 gives M1 0.00915920 s^-1 for the 6300 group,
62.69% above evaluated NIST; its E2 residual is +77.40%. The 6364 M1/E2 residuals
are +65.60%/+101.92%. Stabilizing native stopping leaves the large discrepancy.

These are model one-photon lifetimes, not adopted physical values. At n=5 the
6300-group M1 rate is 0.00916064 s^-1, about 62.7% above the evaluated NIST
0.00563 s^-1; the corresponding E2 rate is about 77.3% higher. The 6364-group
M1 and E2 residuals are approximately +65.6% and +101.7%. These NIST comparisons
are evaluated/theoretical, not laboratory gf tests. Despite close B/C agreement
for the model 1D2 lifetime, the external rates disagree substantially. The 1S0
gauges also remain inconsistent. Gauge agreement cannot calibrate absolute accuracy.
Core, Hamiltonian, radial and independent relativistic-method convergence and
uncertainty calibration remain owed; higher-order multiphoton contributions are
not certified. No E1 conversion is used for a forbidden rate.

## Primary sources and an actual lifetime comparison

The held Amarsi2019 paper points to Amarsi2018's oxygen atom. The latter's printed
page 4 calls Hibbert1991 oscillator strengths experimental. Hibbert's own method
description is CIV3 theory with empirical diagonal-energy adjustments. Both
statements are preserved; RYA-1352's primary-origin rule keeps these theoretical /
evaluated, without inheriting the secondary experimental label.

The Li2023 author preprint is now accessible. Selected method, comparison and
reference sections and Tables 2/3 were inspected. It supplies leads, not new
independent experimental observations. Following those leads to Brooks1977's
Table 3 shows that its 1302 A entry is explicitly borrowed from Smith1971.
Counting the two citations as independent measurements would be wrong.

Smith1971's primary experimental section and Tables 1/2 were inspected, including
visual inspection of the rotated Table 1. Its **direct** beam-foil measurement for
O I 2p3(4S)3s 3S(o), uniquely J=1, is **1.79 ± 0.17 ns**. The source combines
statistical rms fitting errors with its approximately 4% post-foil velocity error,
but does not specify enough about combination/confidence for Gaussian calibration.
The cascade-corrected two-exponential/background fit and common velocity scale are
retained as dependencies. Its recommended f=0.048±0.003 combines earlier studies;
it is not counted as another direct observation.

Using the frozen coverage_v8 wavefunctions, the new inventory evaluates every
angularly permitted one-photon channel to the six lower levels in the NIST level
inventory. Ground levels permit E1/M2/E3, and metastable 3s5S2 permits M1/E2/M3;
J upper=1 and J lower<=2 limit k to three. All 13 channels, including computed
zero E3 rows, are retained. E1 rates exactly reproduce the previous UV table.
Model-space lifetimes are **5.51246 ns (B)** and **2.59095 ns (C)**, residuals
**+3.72246 ns** and **+0.80095 ns** relative to the direct experiment. Conditional
model one-photon branching fractions are retained; they are not adopted physical
lifetimes/BFs. Higher-order multiphoton effects and configuration/radial/Hamiltonian
convergence are not certified. Unknown theory uncertainty stays null; no z score,
chi-square, coverage claim or experimental grade is assigned.

This is an exploratory **nonblind** comparison. The value was already exposed in
Li2023 before primary acquisition. Implementing-agent source inspection is not a
completed distinct-curator/human review or an independently held blind experiment.
Copyrighted source PDFs are preserved privately in Sirius `sources_curation_v2`,
with URLs/byte hashes and inspection limits in `source_receipt.json`; they are not
redistributed in Git. NIST source bytes are retained unchanged. Bibliography entries
record these exact limits; no whole-paper read is claimed for selected excerpts.

## Reference coverage and prospective elements

All nine physically identified, source-qualified NIST 3p5P–4d5D components now have
evaluated/theoretical comparisons in `evaluated_4d_validation.json`. They are not
forced onto individual rounded 615x source-grid aliases. Together with the existing
15 permitted comparisons, this covers the 24 permitted physical components; the
two forbidden reference groups retain their separate operators. The nine individual
rounded-grid aliases remain unresolved. No evaluated-theory row becomes a laboratory
gf accuracy test.

The new canonical 28-target/census/current-feed screen corrects the earlier Na/Mg
methodological shortlist: both already have active products. Li, P and S have no
active feed products in this snapshot and are prospective candidates for separate
curation. **No active feed does not mean never processed**; the historical tracker
is stale and remains a separate field. No unfinished element or existing LAB label
is treated as blind. Primary experiments, isotope/HFS, common upper-level/multiplet,
lifetime/branching and calibration-scale dependencies still require separate review
and an actual prediction/experiment freeze before any prospective reveal.

## Reproduction and verification

All solver work used isolated pinned Sirius builds, single-thread environment and
nice 19; other workers and production installations were preserved. Mac Python
3.12.13 was used only for primary-source inspection, metadata and native-text replay.
Sirius orchestration used Python 3.14.4 (standard library); numerical solvers used
GNU Fortran 15.2.0. No new synthesis/grid input was routed to the Mac.
Drivers create new directories, snapshot themselves and hash all inputs/outputs.
Native stdout/stderr includes every command prompt and executable SHA256.

Sirius launch commands (task-root-relative paths; choose new output directories
for another attempt, and retain the pinned toolchains/seed directories):

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 19 python3 rya1354_active_space.py --grasp-bin toolchains/grasp/build/bin --seed solver_smoke_v10 --n-min 5 --n-max 5 --previous-ladder runs/valence_sd_v4 --new-orbitals hydrogenic --out runs/valence_sd_v6
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 19 python3 rya1354_mchf_ladder.py --atsp toolchains/atsp2k-current --seed solver_smoke_v10 --n-max 5 --out runs/mchf_ladder_v3
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 19 python3 rya1354_rci_sensitivity.py --grasp-bin toolchains/grasp/build/bin --ladder runs/valence_sd_v4 --out runs/rci_sensitivity_v1
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 nice -n 19 python3 rya1354_forbidden_ladder.py --grasp-bin toolchains/grasp/build/bin --seed runs/coverage_v8/ground --n-min 4 --n-max 5 --resume runs/forbidden_ladder_v1/n4_ground --scf-cycles 500 --fixed-damping --damping-factor 0.9 --out runs/forbidden_ladder_v5
```

New capsules and their full byte hashes are recorded in the derived JSON summaries:

- `radial_decay_controls_v1.tar.gz`: initializer controls and full 3s triplet decay inventory.
- `n5_and_mchf_control_v1.tar.gz`: successful n=5 GRASP and native MCHF software control/build logs.
- `mchf_ladder_v3_all_attempts.tar.gz`: all three independent-method attempts, including the rejected first transition path.
- `rci_sensitivities_v1.tar.gz`: all Hamiltonian/core cases and synthetic warning control.
- `forbidden_ladder_v1.tar.gz`: initial ground-level sequence, including the unqualified n=4/5 diagnostics.
- `forbidden_numerical_controls_v5.tar.gz`: all four numerical continuations, including failed attempts and completed native stopping. SHA256 `c85eaa349a722949ed5314f39ce87a606a50adf4d8b0c525f9832af66528233d`.

Only `mcp.*` angular caches are omitted, retained on Sirius. Every other receipt-listed
artifact is checked, including failed attempts. Exact replay commands:

```sh
python scripts/rya1354_evaluated_4d.py
python scripts/rya1354_lifetime_validation.py
python scripts/rya1354_two_method_summary.py
python scripts/rya1354_rci_summary.py
python scripts/rya1354_forbidden_summary.py
python scripts/rya1354_forbidden_summary.py --qualified
python scripts/rya1354_prospective_screen.py
python -m pytest tests/test_rya1354_decay_inventory.py tests/test_rya1354_nonblind.py tests/test_atomic_benchmark_rya1354.py tests/test_grade_governance_rya1352.py tests/test_line_census_rya1353.py tests/test_line_match_rya1033.py -q
```

The bug ledger records the orbital-order false completion, rejected native-air
coordinate, source-origin/reused-observation hazards and stationary-mass warning.
It also records the unpropagated SCF iteration-limit failure; later layers cannot
silently inherit such a failed frozen basis. No new numerical criterion is inferred
from hosted GitHub CI behavior.
The prior n=5 failure row now points to the successful controlled continuation.
Production code, canonical/reference line lists, products, catalogs and stellar
parameters remain unchanged. Hosted CI remains infrastructure metadata, not science
evidence; the known bibliography vocabulary failures remain separate from focused
tests. PR 616 stays draft, RYA-1354 In Progress, with calibrated uncertainty,
reviewed adoption criteria and independent human scientific review still owed.

Checkpoint validation: **63 focused tests pass** (final native replay and refusal checks). Register freshness
passes with zero state surfaces changed. Library audit still exits 1 on the existing
civis2013/mallinson2024 vocabulary errors before library traversal. Raw NIST HTML
source whitespace is preserved; authored files are checked separately. The known
merged-ticket SEQUENCE hook exception remains disclosed for draft checkpoint pushes.
