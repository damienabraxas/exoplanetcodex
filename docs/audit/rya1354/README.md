Latest continuation: [explicit relativistic component constraints](component_constraint_v8.md).
Previous: [Cr candidate intake and oxygen node control](candidate_intake_and_node_control_v7.md).
Previous: [terminal MR audit and numerical continuation](terminal_continuation_v6.md).
Previous: [n4 multireference and orbital-relaxation work](mr_continuation_v5.md).

# RYA-1354 — O I atomic-structure pilot checkpoint

**Status: exploratory candidates and benchmark infrastructure; scientific HOLD.**
Latest: [balanced lifetime, core and orbital-model continuation](lifetime_continuation_v4.md).
Previous: [numerical stopping and physical-model error investigation](error_resolution_v3.md).
Previous: [two-method, core/Hamiltonian, forbidden-line and primary-lifetime continuation](research_continuation_v2.md).
This is not a completed blind benchmark. Calibrated uncertainty, reviewed two-method convergence and adoption criteria remain owed.
See [the 2026-11-08 continuation](continuation_2026-11-08.md) for the new nonblind
reference matrix, correlation ladder and native operator attempts. A named external
curator is not a prerequisite for exploratory calculations or honest external validation.

Base `adefcaf073ff70035954cf45a9a2283774ac7d06`, branch
`ryandamienschmitt/rya-1354-oxygen-pilot-blind-mcdhf`, worktree `codex/rya1354`.
Heavy solver compilation/calculation ran only on Sirius under the existing
`work` root, in `codex/rya1354`; no shared engine installation or PATH changed.
Other workers' O Bridge and Al/Si jobs/worktrees were preserved.

## What actually ran

1. Built pinned GRASP and the official ATSP book MCHF package; evaluated the modern
   ATSP2K biorthogonal path in separate research copies. Exact pins, compiler,
   libraries, commands, observed failures and upstream validation are in
   [toolchain_provenance.json](toolchain_provenance.json). The initial checkpoint sources were not patched. The later numerical investigation
   uses a separately declared convergence-only variant; see the latest dossier.
2. Ran GRASP's quadrature test and preserved serial integration workflow. The latter
   is its small uranium fixture, not an oxygen validation. Quadrature returned
   0.50000000000001932 for 0.5 (difference 1.93e-14).
3. Ran O I configurations `1s2 2s2 2p3 3s` (J=2 odd) and
   `1s2 2s2 2p3 3p` (J=1,2,3 even), zero configuration substitutions, point static
   Z=8 nucleus, default numerical controls, at most 100 SCF iterations. The solver
   stopped after 5/7 iterations respectively; a separate numerical tolerance study
   is not established. Optimized all radial orbitals; selected lowest ASF per J
   block, standard statistical weights for the three upper ASFs.
4. Obtained native LS compositions identifying the configured lower 5S(o) and
   upper 5P terms at approximately 99.998 percent. Biorthogonal transformation of
   separately relaxed Dirac orbitals preceded native E1 transition evaluation.
   These labels identify the calculation; they do not assert a canonical adopted
   transition/component/isotope join or accurate experimental wavefunctions.
5. Ran an independent nonrelativistic ground-state HF sanity calculation and
   limited-complex MCHF 5S/5P calculations. Gencl retains one lower and three upper
   CSFs; it is not an identical configuration model to GRASP. The older angular
   path rejects the open-shell overlap. The modern path initially rejected book
   coupling/eigenvector records. A value-preserving native-layout exporter now
   makes the modern biorthogonal E1 calculation complete. This is a limited
   theoretical comparison, not a converged or experimental benchmark.

The [retained capsule](../../../data/audit/rya1354/solver_smoke_v10.tar.gz) contains
all native inputs, output wavefunctions, LS compositions, energies, gauge tables,
stdout/stderr, per-file hashes and executable hashes. Its SHA256 is
`fe096c8978a577409454a52c1f599458b1599ec59f8241fcadac81ac6f7322f2`.
It contains generated results, not upstream source code or experimental data.
Final runner exit is **1**, intentionally retaining the rejected older angular
path as a scientific failure even though the modern path now completes. The
superseded v6 capsule and candidate v1 remain preserved as failure evidence.

## Theoretical candidates — no admission or uncertainty claim

[exploratory_candidates_v3.json](../../../data/audit/rya1354/exploratory_candidates_v3.json)
is generated strictly from the native gauge table and LS labels. Vacuum wavelengths
come from the computed energy difference, not a measured wavelength substitution.

| Upper J | Native vacuum wavelength (A) | A, Babushkin (s^-1) | gf, Babushkin | Gauge difference in A |
| --- | ---: | ---: | ---: | ---: |
| 1 | 8322.9089 | 3.18281e7 | 0.991608 | 2.998% |
| 2 | 8321.2675 | 3.18466e7 | 1.65299 | 3.010% |
| 3 | 8318.7824 | 3.18754e7 | 2.31489 | 3.028% |

These are the configured states associated with the familiar nominal 777 nm
triplet, but the uncorrelated calculation's native energies place them elsewhere.
No shift to the known line positions was applied. The difference is a reason to
study correlation/energy convergence, not an invitation to tune gf or a Solar
abundance. Gauge difference is `abs(A_B-A_C)/max(A_B,A_C)`; **it is not a 3% error
bar** and does not establish absolute accuracy. Both gauges and line strengths
are retained, with unknown uncertainty/covariance rather than invented numbers.

The independent MCHF LS multiplet gives gf(length)=4.95623 and
gf(velocity)=4.80704. The sum of GRASP Babushkin component gf is 4.959488,
about 0.066% higher. These are different small configuration models and
Hamiltonians, with native Rydberg conventions retained (MCHF reports
109733.5501 cm^-1). The agreement is not a controlled relativistic correction
or a calibrated accuracy estimate. Common omitted correlation can affect both.
The three upper fine-structure components are compared to the LS multiplet
**only through the summed gf**, not by treating the term-average A as three
measured component rates.

No full downward-rate inventory exists, so lifetime and branching fraction are
null/HOLD. Intercombination and other multipoles cannot be silently omitted.
No Landé g, HFS/isotope constants, broadening, collision or photoionization data
are inferred from this E1 run. No NIST experimental grade is attached to theory.

## Failures are part of the result

**CRITICAL — zero process exit can hide invalid atomic data.** Book Mltpol reports
`TOO MANY ELECTRONS FOR THIS CODE` for the dominant matrix element. Lstr then emits
an extremely small gf, which is invalid as an independent result. Modern ATSP2K
`biotr`, fed the retained book-version coupling records, reports
`CONFIGURATION DATA WRONG, PROGRAM HALTED` and still exits zero. Do not repair this
by changing term/seniority labels until the files pass. A native Si E2 control
first reproduced the pinned upstream reference file **byte-for-byte**, demonstrating
the modern executable works on its own format. Inspection then located the
(5X,A3) versus (1X,A3) coupling columns, missing closed/total orbital inventories
and coefficient placement differences. The exporter reformats those records,
preserving energies, every term/seniority label and each coefficient exactly.
The final native oxygen result is now available. The
reporting runner now propagates these observed semantic refusals into a nonzero
overall exit and the final candidate summary rejects the old comparison and identifies the modern
result as an exploratory multiplet comparison. Initial v7/v8 coefficient-format
failures and v9 extra operator request are retained; final v10 requests only E1. The native
quantum labels, coefficients and every failed attempt remain unchanged.

**WARNING — unique labeling exceeds resources.** GRASP's optional unique-LS-label
postprocessing attempted 25.6 GB even for the tiny lower-state case and failed.
[The receipt](failed_unique_label_receipt.json) and Sirius `solver_smoke_v2`
stderr preserve it. Ordinary LS composition output succeeds and is what is used;
no unique-label result or solver-source fix is claimed.

**Build limitations.** Initial CMake failed BLAS discovery because the development
link names were absent; exact installed `.so.3` files resolved the build. GRASP
has 34 legacy compiler warning occurrences in the retained build log. Book and
modern MCHF builds warn about legacy argument/COMMON layout. Modern `src90`
installation also reports unavailable MPI utilities even though its top-level
make exits zero. Serial libraries and the serial biorthogonal executable were
then built explicitly; that successful build did not solve input compatibility.
No compiler warnings or numerical failures were treated as scientific validation.

## Blind benchmark gate and calculation plan

The [v0 plan](../../../config/atomic/rya1354/benchmark_plan.v0.json) is a **draft**,
not a preregistration. It covers physically identified O I UV/permitted/IR and
separate forbidden-line groups, active-space/core correlation, Hamiltonian and
gauge comparisons, state tracking, explicit experimental-energy substitution,
complete decay inventories and correlated multiplet uncertainty. O II is deferred
until a target/diagnostic rationale exists; census presence alone is insufficient.
There are no approved numerical acceptance thresholds yet.

Existing project gf, NIST, theoretical and published experimental-comparison tables
have already been exposed. Their values cannot become a genuinely withheld
benchmark simply by renaming files or splitting multiplet components. The earlier curator request is superseded by the 2026-11-08 issue amendment:
curation is a distinct role, potentially agent-assisted with human review, not a
requirement for an external employee. No genuinely withheld experimental set,
model freeze or prediction freeze is established; exploratory work continues.
The historical v0 plan is retained; [v1](../../../config/atomic/rya1354/benchmark_plan.v1.json)
records the amended exploratory and nonblind validation routes.

The curator must inspect actual primary experiments, source uncertainties,
branching/lifetime normalization, complete physical identity and overlap with
the calculation/calibration material. Freeze experiments by shared upper-level,
multiplet and common experimental scale; never leak related components across
calibration/test partitions. A source called NIST, an accuracy letter or agreement
with a theoretical compilation is not proof of an independent experiment.

Required sequence: curator commits sealed bytes/hash with an independent witness;
reviewer approves model choices, thresholds, uncertainty model and exclusions;
freeze all model inputs and predictions before reveal; score every preregistered
row and retain failures; scientific reviewer assesses residuals, coverage and
common-mode errors. If genuinely unexposed measurements cannot be established,
call the study preregistered external validation, not a blind experiment.

`scripts/rya1354_benchmark.py` supplies tested controls for necessary angular
selection rules, E1 A/gf conversion with vacuum A and upper degeneracy, complete
decay-inventory lifetime/BF, and covariance-aware log-residual scoring. E1 conversion
refuses M1/E2. Angular parity/triangle checks are necessary, not sufficient;
relativistic mixing means a spin-forbidden label does not imply zero amplitude.
No forbidden-line rate was computed in the original checkpoint; the dated continuation now retains native M1/E2 diagnostics on HOLD. No guessed E1 formula supplies a missing
M1/E2 matrix element.

The scoring CLI requires frozen prediction and benchmark byte hashes, exact ordered
transition identity, matching units/quantity, a frozen model reference, explicit
non-Solar-calibration metadata, curator attestation, primary experimental origin,
declared independent evidence and compatible uncertainty convention. Invalid,
unknown or singular covariance refuses calibration instead of reverting to
independent errors or adding numerical jitter. Full covariance is in log10-space,
including explicit theory/experiment cross-covariance assumptions. Shared multiplet
systematics do not shrink as independent line noise.

This records attestations and hash integrity; software cannot independently prove
that a person never saw a benchmark. A score always remains review/HOLD. It neither
validates the experiment's lineage nor assigns a scientific quality grade. No
real benchmark residual/coverage/chi-square has been generated in this checkpoint.

## Relation to existing O work and science policy

Mandatory preflight inherits the completed RYA-1352 historical policy audit and
the RYA-1353 committed-holdings census, now both on main. Root/bootstrap, science,
ledger, protocol, uncertainty and product authorities remain unchanged. The
previous full-ticket/comment archives are in `docs/audit/rya1352/`; the current
RYA-1351/1353/1354 and RYA-1365 descriptions/comments were checked. External
exact-title GitHub Project Register/Target Star List/Canonical Element List artifacts
remain unidentified; local register, catalogs, `config/stars.yaml` and
`data/config/elements_master.json` are the explicit working authority pointers.

RYA-1353 has tuple candidates, not complete upper/J/isotope/component joins. No
new canonical physical match or stellar observability is asserted here. The
primary-only versus broader gf-grading conflict is still not a scientific
adoption decision; this work assigns no new Codex/Deep membership.

Current O Bridge is separate: RYA-1279 remains unsigned; RYA-1281's latest retrieved
checkpoint is `07e5c41e` (2026-11-08), still experimental/HOLD, with native H/density
normalization and earlier collision/background/reference gates. Those RT/chemistry
diagnostics do not certify an atomic-structure solver or a new gf. Likewise this
pilot does not settle the Solar O abundance, NLTE collision/atom completeness,
atmospheres, blends, continuum or tellurics. No stellar spectrum was fitted here.
IR abundance publication still requires verified exact-holding telluric correction
and the existing RYA-587 budget and human approval.

RYA-1365 is a separate optional learning track. These exposed candidate outputs
and fictional unit-test fixtures can be used there; no seven labs or learner
sign-off are claimed, and sealed experiments must stay outside learner material.

## Reproduction and verification

Source pins/build commands are in `toolchain_provenance.json`. Resolve the work
root from `config/path_register.yaml` on Sirius. Build into task-owned directories;
retain source revisions, licenses and compiler warnings. The older held MCHF
archive was inventoried before external acquisition; it includes precompiled old
binaries that were not used. Fresh official source builds were used instead.

```sh
# Sirius only; choose a NEW output path. No existing run is overwritten.
python3 scripts/rya1354_solver_smoke.py \
  --grasp-bin <grasp>/build/bin --atsp-bin <atsp-book>/bin \
  --atsp2k-bin <atsp2k>/bin --out <new-task-run>
# Expected exit 1: rejected book angular output is intentionally retained.
# Modern upstream control in a separate NEW directory:
python3 scripts/rya1354_atsp_control.py --atsp <atsp2k> --out <new-control-dir>
# Expected exit 0 and byte-identical reference E2 output; this is a software control.

# Offline replay of the retained result capsule, using a NEW extraction directory.
tar -xzf data/audit/rya1354/solver_smoke_v10.tar.gz -C <new-extraction-directory>
python scripts/rya1354_summarize.py --run <new-extraction-directory>/solver_smoke_v10 --out <new-candidates.json>
python -m pytest tests/test_atomic_benchmark_rya1354.py tests/test_grade_governance_rya1352.py tests/test_line_census_rya1353.py tests/test_line_match_rya1033.py -q
```

Mac offline tests use the existing `codex/venv/bin/python` (3.12.13); no solver
computation runs on the Mac. The combined focused suite passes **48 tests**.
It checks retained-file checksums and exact candidate replay, A/gf unit/degeneracy
consistency against the actual native table, failed zero-exit solver controls,
identity/hash/exposure/unit/source/covariance rejection, distinct multipoles,
complete-decay requirements and correlated multiplet arithmetic. The modern native O I transition result and the old refused path are tested
separately; no test claims experimental or converged atomic-data validation.

Refinement debt report is informational (21 registry rows, 10 not-fired ticket
rows, 5 no-ticket rows); none discharged. Bibliography library audit remains
blocked by the existing civis2013/mallinson2024 verification vocabulary errors.
Line-key audit retains the known 31 findings across 12 existing files; none in
these new scripts. No full-suite or library-success claim is made. Scientific
product/canonical/production code comparisons against the base are unchanged.

Bibliography: held Tachiev/Froese Fischer 2002 and MCHF 1991 documents were
inspected. Three source records added (`froesefischer2019_grasp`,
`froesefischer2007_atsp2k`, `li2023_o_atomic`), with primary identifiers and explicit
inspection limits. Publisher dates GRASP volume 237 to 2019 despite its upstream
README's 2018 citation year. No whole Li 2023 paper acquisition/read is claimed:
publisher/indexed method/comparison excerpts and author metadata were inspected;
direct full-text requests failed. No experimental source was adopted from excerpts.
Sources: [GRASP](https://github.com/compas/grasp),
[MCHF package](https://github.com/compas/atsp),
[Li et al. O I work](https://doi.org/10.1051/0004-6361/202245645),
[NIST ASD conventions](https://physics.nist.gov/PhysRefData/ASD/Html/lineshelp.html).

Two characterized bug-ledger rows record semantic solver failures and the optional
label-memory failure. No source spectrum, frozen reference, canonical gf,
classification, published abundance or uncertainty changed. No production run,
scientific publication, merge or Done transition is requested or implied.

## Exact continuation

The native modern software control and value-preserving record exporter now
support the limited O I multiplet comparison. Keep the initial incompatible inputs
and old angular output as failures, not physical predictions. Continue the v1 correlation and UV/IR/forbidden-operator work without waiting
for an external curator. Maintain the exposed nonblind matrix and resolve primary
provenance, physical identities and common-scale dependencies. Establish a distinct
curation role and genuinely unexposed sources only if undertaking a prospective
blind test; freeze that model and predictions before reveal. Reviewed criteria and
calibrated uncertainty remain prerequisites for adoption. The current three uncalibrated E1 rows are a starting
diagnostic, not fulfillment of the full RYA-1354 blind-benchmark acceptance.
