# RYA-1352 — grade governance and theoretical evidence proposal

Status: **DRAFT — scientific decision required; no adoption or publication**.
Audit date: 2026-11-08. Base: `3baabba95f728342f069f708dc815e729a1d0447`.
Branch: `ryandamienschmitt/rya-1352-grade-governance-audit-codexdeep`.
Worktree: `codex/rya1352` on Ryan's Mac. Other worktrees, including ongoing
RYA-1350 Al/Si and RYA-1311 Ca work, are preserved. Current-main audit only;
unmerged branches are not silently treated as adopted state.

## 1. Result and decision required

The requested primary-lab-only Codex/Deep contract **does not describe current
main**. It describes the August ratification and the new RYA-1351/1352 mission.
Commit `ae521da874b862f830ccc2d7e0152ec6b9f0f135` (2026-10-09 UTC) explicitly
broadens Codex/Deep to any stored published gf uncertainty or NIST AAA through C.
Its message and code attribute this to Ryan. The retrieved RYA-1233 discussion
does not independently contain that October ruling. This is an authority conflict,
not permission to narrow the implementation or relabel existing products.

**Ryan's decision:** ratify a versioned future primary-experimental Codex/Deep
policy versus the October broad policy, and decide how existing labels are to be
explained or migrated. Preserve all current classifications and scientific values
until that decision and a separately reviewed migration. This dossier does not
declare either historical scientific policy silently superseded.

The proposed `THEORETICAL_CODEX` is a separate candidate selection, initially held.
Its name, acceptance thresholds, review roles and production integration require
scientific approval. The proposal cannot assign primary experimental provenance to
theory, even after a successful benchmark.

## 2. Preflight, history and authoritative surfaces

Read the repository bootstrap, LEDGERS authority routing, DEV_CYCLE, DATA_SAFETY,
science standards, element protocol, product package, conventions, uncertainty
contract, the ordered ledger surfaces, current element/target configuration and
the canonical gf/line-set schemas. The repository equivalents of Project Register,
Target Star List and Canonical Element List are `CODEX_STATE_REGISTER.md`,
`data/catalog/system_catalog.csv` plus `config/stars.yaml`, and
`data/config/elements_master.json`. No separately named GitHub Project Register
file was found; this mapping is explicit rather than a claim to have inspected an
unidentified external board. Main's ledgers retain old current-as-of stamps;
executable behavior and measured artifact counts are therefore stated separately.

Descriptions and all returned comments for RYA-1351, 946, 1169, 709, 587, 161,
945, 822 and 1131 are archived in [linear_preflight.json](linear_preflight.json).
Follow-up RYA-1233, 1185, 1213, 1174, 1278, 1279 and 1281 topic excerpts are in
[linear_followup.json](linear_followup.json). Every response reports
`hasNextPage: false`; comments are assessed by creation time, not server order.
The task itself had no comments at intake. Searches covered theoretical grading,
related O/Fe work and Bride work; RYA-1353 owns the actual O/Fe census, RYA-1355
the Fe calculation stress test. They are not duplicated here.

| Date / evidence | Definition / implication |
| --- | --- |
| RYA-946, 2026-08-22 comment `498f6b34` | Headline graded pool admitted laboratory or NIST C+; diagnostic ungraded pool stayed separate. Historical policy, later amended. |
| RYA-946, 2026-08-30 comment `2b66928c` | Codex = primary lab, depth <=0.60; Deep = same provenance, >0.60, synthesis. Asplund = independent external set membership. <0.05 separately below observability; Consistent evaluated tier parked. Explicitly supersedes earlier NIST admission. |
| RYA-946, 2026-09-03 comment `62392940` | Reference construction may amalgamate evaluated theory for lab-poor species, explicitly never calls it lab and keeps the Codex/Deep firewall. No borrowed Solar-fitted gf for absolute validation. |
| RYA-1185, 2026-09-03 comment `1b74806e` | Asplund Grade renamed Reference Grade; specific line_set preserves reference identity. |
| RYA-1213 and Fe emitter | Reference also covers the full primary laboratory pool, on our gf, with no depth gate. This is distinct from AGSS21 replication although both display Reference Grade. |
| `ae521da8`, current gf predicate and tests | Any stored sigma, LAB substring, or NIST AAA..C admits a row; PR2024-CALC explicitly tested as graded. This conflicts with the new mission's primary-only requirement. |
| RYA-161, June corrections | Astrophysical gf is differential-survey-only and parked; it cannot erase the Solar absolute residual or validate its own calibration. |

RYA-1169 distinguishes ten source-used Si rows from two rejected violet rows;
no familiar wavelength alone proves Asplund membership. RYA-945 documents the
`RU` collision (theoretical Raassen/Uylings versus laboratory Ruffoni); source
labels cannot determine experimental origin. RYA-822 distinguishes GF-NIST from
GF-LAB and warns about attaching a source uncertainty to a different gf value
(`SCALE-MISMATCH`). RYA-1131 retains molecular ambiguity/fallback rows rather than
promoting them. Molecular evidence remains in its existing separate domain.
Bride G1/G4 require independent atom/reference evidence and gated numerical
validation; their O I 777 work does not validate a new atomic-structure gf dataset.

Refinement debt report succeeded with the project Python: 21 registry rows,
10 not-fired ticket rows and 5 no-ticket rows. These include N atomic-data/model
debt and other elements; this governance task clears none. O/Fe tracker entries
do not establish absence of atomic-data debt. No phase close or freeze is claimed.

## 3. Current implementation and measured baseline

| Surface | Verified behavior |
| --- | --- |
| `pipeline/gf_grades.py:512` | `is_gf_graded` admits LAB substring, any numeric non-null sigma, NIST AAA..C. No method/origin check or positivity/finiteness validation. |
| `pipeline/element_prepare.py:206`, `pipeline/run_matrix.py:438`, `pipeline/line_sets.py:144` | Shared predicate is consumed in preparation, matrix grading and generated line sets. |
| `scripts/derive_band_products.py:360`, `:443`, `:467` | Same predicate feeds Codex/Deep; depth <=DEPTH_HI versus >DEPTH_HI; unknown depth is reported and excluded. Old LAB wording remains in diagnostics/docstrings. |
| `scripts/line_accounting_rya709.py:55` | DEPTH_LO=0.05 and DEPTH_HI=0.60. Operational floor semantics differ by route/vintage; do not create a new universal floor from this proposal. |
| `pipeline/al_grade_verification.py:105` | Older primary-lab membership calculation coexists with widened predicate. |
| `scripts/rya1178_emit_fe_schema.py:154` | Fe label/variant tables still describe Codex/Deep as laboratory-only; Reference has published-set and full-lab variants. |
| `pipeline/cno_grade.py:39` | Non-Fe grade stamping is selector/tier based, including default Codex Grade. A display grade alone proves neither origin nor experimental gf. |
| `data/reference/line_sets/{REGISTRY.csv,SOURCES.yaml}` | Source memberships, gf_basis, citations and source-used sets; preserve published gf independently of canonical adopted gf. |
| `data/linelists/canonical_gf.csv` | One adopted row with gf_tier/tag/sigma/DOI plus some alternatives; lacks a full versioned evidence graph, source-specific covariance and separate evaluated-underlying origin. |

[baseline.json](baseline.json) records source SHA256 for protected code/config,
canonical/reference data, catalogs and product feeds. There are **169,723 canonical
rows; 7,809 pass the current predicate**, of which 4,713 are not exact `gf_tier=LAB`.
Accepted tiers: LAB 3,096; NIST-C+ 4,164; PR2024-CALC 391; C2013-QDT 33;
E2024-ASTRO 7; PUBLISHED-SET 81; OTHER 12; VALD3 25.
These are **classifier counts, not measured/eligible line counts** and are not a
verified origin census. No inferred lab-only reclassification is adopted.
The current canonical population has zero accepted nonpositive/infinite stored
sigmas, although synthetic inputs demonstrate the predicate accepts them.

Zero protected paths changed; all accepted IDs and their ordered hash match the
base. All committed `data/products/` bytes match, which proves no local product
classification, abundance or uncertainty changed. Product grade counts are
recorded per feed in the report. No depth was remeasured, generator rerun or
external live site audited; the proof does not certify the existing science.

## 4. NIST semantics and uncertainty

NIST accuracy describes an evaluated transition strength, not experimental
origin or our measurement quality. ASD critically assesses values and their
uncertainty; a theoretical calculation can receive a published evaluation.
Preserve that exact evaluation and underlying theoretical provenance together.
Never invent a NIST class for new Codex theory or reinterpret a cited class as
proof of laboratory measurement. [NIST ASD help](https://physics.nist.gov/PhysRefData/ASD/Html/lineshelp.html).

| Accuracy | Relative bound |
| --- | --- |
| AAA / AA / A+ / A | <=0.3 / 1 / 2 / 3 percent |
| B+ / B / C+ / C | <=7 / 10 / 18 / 25 percent |
| D+ / D / E | <=40 / 50 / >50 percent |

Retain a prime suffix: it indicates LS-coupling component derivation and may
degrade accuracy. These assessments predominantly concern systematic effects;
they are not automatically Gaussian one-sigma measurements. Current project
positive dex conventions (e.g. log10(1+p)) must be identified as conventions.
For a bounded fractional error p<1, log-space bounds are asymmetric:
log10(1-p), log10(1+p). E does not supply a finite upper uncertainty.
Unknown confidence and covariance remain explicit; never convert them into zero.

For theory, convergence and length/velocity gauge consistency are diagnostics,
not independent truth. Blind withheld laboratory comparisons, independent method
checks where physically appropriate, configuration/correlation/relativistic
convergence and uncertainty coverage review are required. Thresholds depend on
species/transition class and must be approved before opening withheld results.
[NIST uncertainty review](https://www.nist.gov/publications/evaluation-uncertainties-atomic-data-spectral-lines-and-transition-probabilities)
provides primary methodological context; this dossier proposes no numeric theory
acceptance cutoff.

## 5. Proposed evidence contract and approval policy

The machine-readable [v0 field contract](../../design/theoretical_transition_evidence_rya1352.json)
defines transition identity, evidence, covariance, adoption, line selection and
holding-specific eligibility separately. It is a **schema design**, not a runtime
validator. No production module imports it. It contains no physical predictions.

`origin=evaluated` must retain `underlying_origin` and original evidence IDs, so
NIST redistribution cannot hide theory. `theoretical_method` identifies method,
code/version, expansion, gauge and observed-level substitutions. Preserve primary
lab, evaluated, theoretical, VALD/Kurucz delivery and astrophysical records as
independent evidence; an adoption points to one record without deleting alternatives.
DOI absence has an explicit reason and stable alternative identifier, never a fake DOI.
Physical joins require species/ion, levels/EP/J, medium/frame and isotope/HFS
identity; ambiguous assignments are held. A source gf and a summed feature gf
cannot share an uncertainty merely because their wavelengths overlap.

Known covariance carries ordered evidence IDs, native convention, units, source,
artifact hash, PSD/symmetry/diagonal checks and common scale/branching/lifetime
terms. Unknown covariance remains held where material; independence requires an
explicit reviewed assumption. The line-data sidecar supplies the existing
`rya587.uncertainty.v1` transition-data component with product response weights;
it does not replace that contract. Sun/target common-mode covariance and matched
physical pools remain mandatory for differential abundance/ratio publication.

Proposed approvals are separate: (1) governance/schema and policy version;
(2) preregistered benchmark and uncertainty plan; (3) named scientific review of
blind results and rejected candidates; (4) per-transition adoption; (5) validated
measurement integration; (6) independent publication approval. Ryan's approval is
required for scientific adoption/publication. Passing a schema/test or reviewing
this PR conveys none of the later approvals.

Proposed website wording: **“Theoretical Codex — candidate, review pending”**;
after explicit scientific acceptance, **“Theoretical Codex — [method], reviewed”**.
Show underlying origin, evaluator/release, source, uncertainty convention,
alternative gf, benchmark status and limitations. Reference Grade also names its
specific published set/full-lab/evaluated variant. Keep the legacy grade visible
with policy version during migration; do not claim reviewed theory is experimental.
No website patch or deployment is part of this PR.

## 6. Staged migration and regression acceptance

1. Review this dossier and settle the August/October/new-ticket conflict. Record a
   ratified decision without editing historical evidence or current feed labels.
2. Implement and test a strict sidecar validator in a separate reviewed change.
   Backfill only verified identities and origin; unknowns remain unknown/HOLD.
   Preserve canonical hash + original IDs, gf and classifications exactly.
3. Build a shadow proposal and row-by-row legacy-versus-proposed diff, covering
   all bands/species and accepted/rejected/unknown rows. No production reader uses it.
4. Require negative tests: evaluated theory masquerading as lab, forged NIST class,
   changed gf with an old source sigma, wavelength-only/ambiguous joins, components
   versus sums, missing DOI alternative, negative/NaN/infinite uncertainties,
   invalid/unknown covariance, Solar-tuned gf, missing exact-holding IR correction,
   benchmark leakage and approval scope mismatch. Test depth at/below/above 0.60,
   <0.05 and unknown depth; preserve reference membership independently.
5. Only after separate adoption approval, wire an explicit opt-in policy version
   and new class into selection/uncertainty/reporting. Test legacy default byte-for-byte
   and compare complete product identity + A + sigma + memberships. A changed
   existing product requires reviewed migration; no overwrite-by-label.
6. Scientific review precedes any production calculation or publication. Schema
   migration rollback selects prior immutable policy/evidence versions; no deletion.

The audit verifies all baseline protected-file hashes against the working tree
and recomputes the accepted-ID hash on every canonical row. Regression tests
characterize widened theory/NIST admission and invalid-sigma admission, and prove
on a disposable Git fixture that product and classifier mutations are detected.
They do not permanently freeze future scientific changes to this baseline.
Existing grade tests are also run. Future migration tests above are explicitly **owed**,
not presented as implemented validation of a deployed theoretical schema.

## 7. Findings, failures and reproduction

**CRITICAL — unresolved governance conflict:** primary-lab-only mission/ratifications
versus current predicate admitting theoretical/evaluated/astrophysical tagged values.
No current scientific result is declared invalid solely from a tier label. Resolve
policy before any affected scientific adoption or publication; preserve legacy results.

**WARNING — source/display ambiguity:** Fe laboratory-only prose and generic non-Fe
label stamping do not establish experimental origin. RYA-1174 already owns the
evaluated-versus-underlying-theory Al defect; this proposal addresses the schema
dimension and does not claim that ticket fixed. Full-lab Reference and source-set
Reference are distinct products despite one display label.

**WARNING — predicate uncertainty validation:** synthetic zero/negative/infinite
sigmas all pass current `notna()` admission. No such accepted stored sigma is found
in the current canonical baseline. Keep a separate future validator; narrowing
the production predicate here would violate zero-change scope.

**SUGGESTION:** scientific policy should have a ratified, versioned owner referenced
by classifiers and display text; neither a git message nor display name is enough.

Initial system `python3.12` lacks pandas; reran with existing project venv (Python
3.12.13, pandas 3.0.3). This is a read-only metadata audit, not a production compute
environment or numerical-engine validation. Initial fetch failed sandbox DNS;
authorized fetch succeeded. Linear `get_agent_context` returned tool-not-found;
repository and ticket evidence were used. Bibliography library audit exits 1 before
library traversal on pre-existing invalid `verified` values (`civis2013`,
`mallinson2024`). No clean-library claim is made and no unrelated entries changed.
NIST web sources above were directly inspected on the audit date; no primary gf
data was acquired. New source URLs are dossier evidence, not published product citations.

Reproduce from repository root (choose new report filenames; audit refuses overwrite):

```sh
git fetch origin
git show ae521da8 -- pipeline/gf_grades.py
python -m pipeline.refinement_debt_join --report
python scripts/generate_sources_page.py --audit-library "<resolved refs.library>"
python scripts/audit_grade_governance_rya1352.py --base 3baabba95f728342f069f708dc815e729a1d0447 --output /tmp/rya1352-new-audit.json
python -m pytest tests/test_grade_governance_rya1352.py tests/test_rya1233_gf_graded_definition.py tests/test_rya1233_grade_every_element.py -q
git diff --check
python scripts/check_register_freshness.py --since-main
```

Use the existing project interpreter `codex/venv/bin/python` for `python` above.
The report contains host/runtime and exact base SHA; no source spectrum, frozen
gold, classification, measured abundance or product bytes were edited. This task
does not change a scientific verdict, so no register/version bump is warranted.
Run-bug-ledger disposition: append characterized governance/display conflict and
latent invalid-sigma admission; RYA-1352 remains their review owner. Final check
results and PR/SHA are recorded in [validation.md](validation.md) and Linear.
