# RYA-1354 — continuation after the amended notes

**Exploratory / nonblind external validation; scientific HOLD.** No adoption,
production, publication, merge, completed blind test or calibrated accuracy claim.
The curator search is no longer a calculation blocker. Curation is a distinct
role, potentially agent-assisted with human review; no external employee or named
experimentalist is required. Existing project, NIST and literature data are exposed.
The historical v0 plan remains intact; `benchmark_plan.v1.json` records this amendment.

Branch `ryandamienschmitt/rya-1354-oxygen-pilot-blind-mcdhf`, worktree
`/Users/ryanschmitt/codex/rya1354`, starting SHA
`03e28b9b77a6f5b69f97617054d295f246f0740b`; original base `adefcaf073ff70035954cf45a9a2283774ac7d06`.
Fetched main `7f023fd7` has independent Ca work; this continuation does not edit it.
Draft [PR 616](https://github.com/damienabraxas/exoplanetcodex/pull/616) remains the handoff.

## Reference accounting and comparisons

`data/audit/rya1354/nonblind_v1/reference_matrix.json` accounts for all **26 O I
source-analysis-grid components**. Their presence in the Amarsi/AGSS21 source grid
does not prove all 26 were finally adopted or fitted in every Reference Grade
product. Actual O product grade/selector/tier/band/holding/treatment scopes and
recorded indicator IDs are retained. Existing canonical tuple matches provide
traceability, not new complete physical joins. Reference Grade does not mean lab gf.

A new official ASD request retains raw bytes, URL, retrieval time and SHA256,
plus configurations, parent terms, term/J, level energies, operator, accuracy
letters and references. Query 3000–10000 A, **vacuum** (`show_av=3`), energies in
cm^-1. Grid wavelengths are air; discovery uses the census's existing vacuum
column. A 0.03 A discovery radius is not a physical matching uncertainty. All
nearby candidates survive; no nearest-hit tie breaking. Isotope/HFS remain unknown.

There are **15 explicitly configuration/parent/term/J/parity-checked compiled-gf
comparisons**: 777 nm (3), 8446 A (3), and 926x A (9). Energy residuals are separate
from loggf residuals. The 777 nm values use the original v10 native prediction;
the other groups use coverage_v8. These small models show loggf residuals roughly
-0.009 to +0.004 dex against the compiled grid while their energies miss by about
541–1046 cm^-1. This illustrates why numerical gf agreement is not validation.

The nine 615x components have overlapping proximity candidates. Their 4d native
rates exist, but assigning one to each rounded grid component requires original
source-level component identity; they stay unresolved and unscored. The two
forbidden grid rows retain separate NIST M1/E2 entries and native operator tables;
no E1 formula or arbitrary sum supplies a missing forbidden gf score.
Thus all rows are accounted for: 15 compiled-gf comparisons, 9 identity holds,
2 operator-review rows, **zero experimental-gf accuracy rows**. No calibrated
chi-square, coverage, uncertainty, lifetime or branching fraction is inferred.

NIST citation receipts resolve T5803 to Hibbert et al. 1991, DOI
10.1088/0953-4075/24/18/010. The publisher abstract describes CIV3 calculations
with empirical diagonal energy adjustments: these are theoretical, not laboratory
gf measurements. T4539 is the Fischer/Saha MCHF/Breit–Pauli paper,
10.1103/PhysRevA.28.3169; only its citation/title was inspected. T5081 resolves to
Baluja/Zeippen 1988, 10.1088/0953-4075/21/9/007; its full method/table remains
uninspected. These inspections do not certify that the source grid used identical
upstream values. Source-grid primary provenance and covariance remain unresolved.
Three bibliography entries state these exact inspection limits.

The canonical target screen covers all 28 targets. Na/Mg are prospective
methodological candidates (light, priority 1), not blind or experimentally approved
selections. Canonical Na I has 113 rows, no LAB-tier row; Mg I has 473, including
16 LAB-tier rows already exposed in project holdings. Neither a LAB label nor an
unfinished survey element creates a holdout. A curator must inspect primary
experiments and group upper levels, multiplets, shared lifetimes/branching and
common calibration scales before any prospective split. Processing status is not
used as a proxy for blindness; live survey-state ranking remains owed.

## Sirius calculations and failures

All solver computation ran under `/mnt/codex-data/codex/rya1354` on Sirius with
the previously pinned research builds, `nice -n 19`, `OMP_NUM_THREADS=1` and
`OPENBLAS_NUM_THREADS=1`. Starting observed load was 15.07/14.08/11.57 on four
logical CPUs, available RAM 8.7 GiB. These execution choices reflect measured host
conditions, not a GitHub-CI-derived science cap. Other workers were left intact.

The valence SD ladder fixes 1s2/2s2 occupations, starts from retained v10
spectroscopic orbitals, freezes n<3 radial functions, optimizes n=3 spd and then
only the new n=4 spdf layer. Native LS output verifies lower 5S(o) J2 and upper
5P J1/2/3. Explicit upper ASF ordinals are 2/3/1 because ground states enter the
even correlation space. No empirical energy substitution or gf fitting occurred.

| Model | Native separation J1/J2/J3 (cm^-1) | Summed gf B | A gauge differences |
|---|---|---:|---|
| n=3 VV-SD | 14501.62 / 14503.99 / 14507.58 | 5.87315 | 40.42% |
| n=4 VV-SD | 12757.66 / 12760.22 / 12764.09 | 5.08822 | 7.73–7.74% |

These large changes do not establish convergence or sigma. n=5 spdfg failed in
5g Thomas–Fermi initialization, with radial-grid-extent warnings. No n=5 value
exists. The unseeded n=3 all-orbital attempt also failed SCF. Both failures remain
in the capsule and bug ledger. Core-valence, 1s policy sensitivity, radial-grid,
Breit/QED and a converged independent MCHF/Breit–Pauli comparison remain owed.

Coverage_v8 completes reference-configuration ground 3P/1D/1S, 3s 3S, 3p 3P,
3d/4d 5D calculations, native UV E1, 8446 E1, IR 3d E1, optical 4d E1, and ground
**M1 and E2** operators. Both M1/E2 same-J restrictions are explicitly disabled.
Native LS labels establish the selected configured states. Fully relaxed 3p
triplet attempts failed SCF even from retained starting orbitals; v8's successful
3p triplet freezes n<3 at the quintet seed and is a distinct sensitivity model,
not a converged replacement. The ground M1/E2 identical f1/f2 state lists produce
duplicate directional records: never sum them. Some E2 length/velocity results
differ by many orders of magnitude; they are diagnostic failures of accuracy,
not useful adopted rates. UV includes weak mixed-spin rows, not just strong lines.

All attempts valence_sd_v1–v5 and coverage_v1–v8 are preserved. Failed generator
inventories, omitted interactive same-J prompts, SCF refusals and initialization
failures are not discarded. The compact capsule excludes only `mcp.*` angular
caches, which remain on Sirius. Every other receipt-listed file is retained and
hash-checked by the offline summarizer. Executable hashes, native stdin/stdout,
wavefunctions, CSFs, LS labels and transition tables are included.

Capsule `data/audit/rya1354/continuation_20261108_compact.tar.gz`, SHA256
`c2d9f81d377b0e2ca9d5390128cd3694715bb5731fee314f96e8e296a439eee4`.
`native_summary.json` preserves parsed native results and all attempt statuses.

## Verification and exact continuation

Mac Python 3.12.13, existing `codex/venv/bin/python`; no atomic solver ran here.

```sh
python scripts/rya1354_native_summary.py
python scripts/rya1354_reference_matrix.py
python -m pytest tests/test_rya1354_nonblind.py tests/test_atomic_benchmark_rya1354.py tests/test_grade_governance_rya1352.py tests/test_line_census_rya1353.py tests/test_line_match_rya1033.py -q
python scripts/check_register_freshness.py
python scripts/generate_sources_page.py --audit-library '/Users/ryanschmitt/Documents/Exoplanet Codex'
```

Native parsing/hash checks and matrix replay pass; **51 focused tests pass**.
`git diff --check` flags only original whitespace in the four byte-preserved NIST
HTML receipts; those source bytes were deliberately not reformatted. Authored
files pass the whitespace check with those exact receipts excluded.
Register freshness passes with the existing stale generated tracker notice.
Library audit still exits 1 on existing civis2013/mallinson2024 verification-state
vocabulary, before library traversal; no new bibliography validation errors.
Production pipeline, products, canonical/reference line lists, catalogs and stellar
config are unchanged against the original base. One additional bug-ledger row
characterizes numerical failures, without asserting an upstream software defect.

Prior GitHub Continuous Integration runs were **cancelled, not passed**. No hosted
CI success or waiting dependency is claimed. Atomic configuration interaction is
part of the physics and unrelated to that service status. The push hook's known
merged-ticket SEQUENCE completeness failure is separate from science checks.

Next: resolve the nine 615x component identities and original gf origins from the
primary source; inspect common scales/covariance; qualify an enlarged radial grid
and high-l initialization before another n=5 attempt; establish independently
reviewed core-valence/Hamiltonian/two-method and forbidden-operator convergence.
Use a **new** Sirius output directory with the committed drivers and the retained
`solver_smoke_v10` seed. Do not keep restarting the same model hoping for a better
gf. Curation and any prospective sealed-test registration can proceed separately.
Scientific adoption still requires reviewed criteria, calibrated uncertainty and
Ryan's review; this checkpoint does not discharge those gates.
