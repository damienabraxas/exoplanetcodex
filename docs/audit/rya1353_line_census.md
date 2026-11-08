# RYA-1353: O and Fe transition census and review dossier

Status: **review dossier; no transition admitted, grade changed, abundance derived,
or scientific product published**. The committed holdings census is reproducible.
A measurement-ready physical-identity census is **not established**: level/component
metadata and exact holding conditioning require further review. This limitation is
part of the deliverable, not a claim that the underlying transitions lack data.

## Baseline and authority

Worktree `codex/rya1353`, branch `rya1353-line-census`, baseline
`3baabba95f728342f069f708dc815e729a1d0447` after fetching origin/main.
The existing RYA-1352 worktree was inspected read-only; its changes are not included.
Current register, ledgers, element tracker, catalogs, product feeds, science standards,
element protocol, data safety, conventions, and project intake/brief/review skills
were consulted. Target and element definitions route to `config/stars.yaml` and
`config/constants.py`; holdings route to the catalog registry and its manifests.
GitHub issue searches for the exact titles Project Register, Target Star List, and
Canonical Element List returned no results through the connector. Those separately
named GitHub artifacts have **not** been positively identified; local authority
pointers do not prove their external identity.

The preflight reviewed RYA-1351/946/709/1131/945 and the inherited
1169/587/161/822 policies and discussion, plus current Bride O work
1278/1279/1137. RYA-1279 remains In Progress: its October 26 checkpoint on draft
PR 598 reports numerical convergence but an independent departure discrepancy and
open atom lineage/background/rate conventions. That unmerged branch is not an
admitted O model in this census. No Fe/O atom or grid calculation was run.

**Policy conflict, owned by RYA-1352:** RYA-1351 says Codex/Deep use primary
experimental gf. Current `pipeline/gf_grades.py:is_gf_graded` (RYA-1233) also
admits stored published uncertainty or NIST AAA through C; its comment cites
Ryan's October 9 decision. Historical RYA-946 comments first ratified primary-lab
Codex/Deep and then an independent evaluated/amalgamated reference tier. Both
claims are retained. The census records the current executable gf mask and
recorded experimental origin separately; it confers neither Codex nor Deep
membership. No competing science rule was resolved by this implementation.

The refinement-debt report exits 0 (informational): it lists 21 registry rows
over 18 elements; the O/Fe tracker rows have no Backlog refinement debt. The
global outstanding tickets are not discharged by this census. No verdict or state
register bump is warranted for an unadmitted audit.

## Measured inventory

Numbers below count **rows in the named snapshots**, not unique transitions in
nature, resolved observed features, or the number of usable abundance lines.

| Species | Canonical rows | Current code gf-graded | Recorded LAB tier |
|---|---:|---:|---:|
| O I | 145 | 93 | 0 |
| O II | 182 | 43 | 0 |
| Fe I | 25,671 | 1,610 | 428 |
| Fe II | 15,369 | 581 | 22 |

All 41,367 canonical rows and their original columns are retained unchanged.
The independent Solar intake has 38,302 O/Fe I/II rows, also retained. It is not
collapsed into canonical rows: its blend flags, depths and damping terms have
different source and feature semantics. Synthetic `solar:N` IDs enumerate the
filtered intake; they are not physical IDs or original CSV line numbers.

The source ledger contains 11,473 rows: held O I/II and Fe I/II ASD snapshots,
493 Fe I experimental rows (465 in the established table plus 28 Ladenburg
rows), 131 Fe II experimental rows, 26 O source-census rows, and 53 Fe reference
rows. Multiple source records of one transition remain separate alternatives.
Of the source rows, 5,424 lack gf in that **source row**. This does not establish
that no gf exists elsewhere. Canonical rows all carry finite gf in this snapshot.

Of the inspected experimental rows, 72 Fe I and 109 Fe II rows have no canonical
tuple candidate at the strict declared tolerances; two Fe I source rows encounter
duplicate canonical wavelengths. These are ingest/identity gaps despite an
experimental datum already being held. They are not targets requiring a new
theoretical gf because an experiment is supposedly missing.

The DOI-only bibliography join resolves 450 canonical rows and leaves 40,917 for
reference-code/missing-DOI review. This intentionally does not declare those
references absent: a code may resolve through VALD/NIST's reference block rather
than a DOI column. Shared ASD DOI keys also do not establish an ASD version.

## Reproducible method and artifact schema

Generator: `scripts/rya1353_line_census.py`. Dependencies are the existing Python
environment, NumPy and pandas. This is lightweight table accounting on the Mac,
not a synthesis, grid generation, fit, or production science run.

```sh
/Users/ryanschmitt/codex/venv/bin/python scripts/rya1353_line_census.py --out /private/tmp/rya1353-replay
/Users/ryanschmitt/codex/venv/bin/python -m pytest tests/test_line_census_rya1353.py tests/test_line_match_rya1033.py -q
```

The output directory must not exist. The reviewed snapshot is
`data/audit/rya1353_line_census/v4/`; the earlier local v1/v2/v3 probes are not part of
the review artifact. CSV payloads above one MB are deterministically gzip-compressed
(mtime zero). `input_manifest.csv` records SHA-256 and bytes for all 49 inputs,
including the generator, current grading/matcher implementations, bibliography,
source tables and available provenance sidecars. `summary.json` records baseline,
interpreter, platform, tolerance, and scope. No external live database is silently
substituted for a frozen pull.

| Artifact | Interpretation |
|---|---|
| `canonical_census.csv.gz` | Every canonical row; original gf, provenance, grades, IDs and EP; separate current gf mask, recorded LAB state, predicted-depth state, diagnostic candidates and holds |
| `solar_intake_rows.csv.gz` | Independent intake, including predicted depth, blend and damping fields; no inferred adopted gf or HFS merging |
| `source_records.csv.gz` | One record per inspected source row; source file and CSV row number, gf presence, source DOI, raw record, tuple candidate, delta gf and explicit identity gate; raw records retain upper/lower energies, J or TP codes where supplied |
| `candidate_review.csv.gz` | Complete source-row review queue with observability, weak/strong prediction, existing-experiment state, tractability, diagnostic value and M-dwarf HOLD; every candidate has an exclusion-from-admission reason |
| `coverage_counts.csv` | Canonical species × configured band × recorded gf tier × current gf mask counts |
| `solar_holdings.csv` | All 16 registered Solar holdings, conditioning declarations, catalog capability and explicit missing-span/mapping state |
| `holding_spans.csv` | Seven instrument-scope manifest spans plus eleven exact runtime HoldingSpec spans; source rows and declared/verified status |
| `holding_line_matrix.csv.gz` | Normalized species/manifest-span candidate line IDs; overlap is provisional until exact product frame, gaps and holding association are checked |
| `engine_scope.csv` | Existing declared element/band engine support, explicitly not per-transition applicability |
| `product_usage.csv` | Current/archived/quarantined published selection, grade, holding and model metadata, pointed to original JSON; no reverse inference of a line identity from aggregate n_lines |
| `exclusion_evidence.csv` | O/Fe curation exclusions and problem-child scopes, verbatim evidence and tickets; scoped text is not silently wavelength-matched into a new exclusion |
| `source_rejected_indicators.csv` | Existing CNO negative-selection ledger, retaining OH molecular context independently of the atomic census |
| `excluded_gf_alternatives.csv` | 72 raw Table 6 alternative cells from Ruffoni 2013, including bounds/missing entries; BF/effective-lifetime and recommended columns excluded from an independent Solar benchmark under RYA-161 |
| `canonical_provenance_resolution.csv.gz` | Exact DOI-to-bibliography resolution or explicit reference-code review |

The band labels within the existing synthesis scope come from
`config/synth_bands.yaml`. Below its floor the label is
`UV_below_synthesis_floor`; other uncovered regions are explicitly outside
configured synthesis bands. These extensions describe census bins, not new engine
support. No Fe IR reach is inferred from an instrument's 53,000 A catalog limit.

## Identity, uncertainty and observation boundaries

All tuple searches restrict species/ion and use `pipeline.line_match.match` with
`require_ep=True`, wavelength tolerance 0.005 A and EP tolerance 0.02 eV. There is
no nearest-row, rounded-key or wavelength-only fallback. The census adds a
conservative refusal for duplicate source wavelengths because the shared matcher
can collapse equal wavelengths without proving gf/component equality. Published
rounded source rows may therefore remain unmatched; tolerances are not widened to
manufacture coverage. Original source match decisions are preserved in raw records.

Results: 6,639 tuple candidates, 4,733 without a candidate at the strict tolerances,
72 missing-identity holds, 25 duplicate-wavelength holds and four ambiguous holds.
**Zero are asserted to be full physical-identity matches.** Canonical `physical_id`
is preserved, but cannot reconstruct missing upper levels, J, multiplet or
isotope/HFS component identity. HFS component counts are not component IDs.
Where the source has those fields, they remain in its raw record for adjudication.

ASD pulls declare mixed vacuum/air convention (vacuum below 2000 A, air above).
Canonical rows here start above that boundary. The source wavelength convention
and raw observed/Ritz values remain distinct. The canonical column name alone
cannot certify every upstream import's frame. Ambiguity/absence is kept as a hold,
never used to import a value or an experimental grade.

NIST accuracy is the evaluator's estimate for the transition probability; it does
not certify experimental origin. Observed and Ritz wavelengths are different
quantities. See the [NIST ASD line help](https://physics.nist.gov/PhysRefData/ASD/Html/lineshelp.html)
(bibliography family `nist_asd`). No new theoretical uncertainty is invented;
source sigma, NIST letters, TP reference codes, and unresolved covariance remain
separate. No RYA-587 abundance uncertainty budget is claimed from this census.

The Solar intake `central_depth` is **predicted**. For example, its strict tuple
candidates include O I 7771.944/7774.166/7775.388 A at EP 9.146 eV and predicted
depths 0.166/0.145/0.114. These are inventory values, not new measured feature
depths. They do not confer Codex membership or prove target/M-dwarf detectability.
Observed feature depth, SNR, local continuum, actual blend fraction, telluric
residual and transition-specific NLTE applicability remain unmeasured here.

Some manifests list several products for the same instrument. The census refuses
to transfer a registry holding's telluric state onto an instrument-only span:
in particular Baker and Reiners IAG products must not share conditioning by
accident. The exact registry declarations remain in `solar_holdings.csv`; span
conditioning remains UNKNOWN for instrument-only spans. The existing AST adapter
in `pipeline.si_protocol_audit.holding_specs` reads eleven exact declared spans
from `scripts/measure_band_ew.py` without importing the measurement driver or
opening an atlas; those rows carry their exact holding registry declarations.
Five holdings retain dynamic-loader or manifest-association holds, including the
raw IDPs, dynamic Kitt Peak readers, the literature-table holding and Delbouille.
These are not no-observation findings. Gaps, boundary frames,
and spectral resolution also require checking. No IR line is reportable on the
basis of this range-overlap matrix.

## Candidate priorities for scientific review

These are predeclared **review queues**, not numerical utility scores, approved
theory calculations, or ranked scientific results. There is no defensible
line-by-line observability/tractability score while the identity and conditioning
axes remain unknown. The machine queue retains that uncertainty.

1. **O I first:** review the admitted reference triplet's atomic identity and the
   other held O source indicators against independent primary atomic data. Existing
   evaluated gf is not a missing-any-gf case. Use the triplet as an independently
   benchmarked structure-method target only after the RYA-1279 model/physics gate.
   Separate lower-excitation/forbidden and high-excitation permitted indicators;
   Ni blending around [O I] is retained as problem-child evidence.
2. **Existing Fe experiments:** the 624 held experimental rows are benchmark or
   ingest work before being classified as missing-experiment targets. Freeze blind
   withheld comparisons before a theoretical calculation; retain Ruffoni's
   Ladenburg column independently of Solar-refined alternatives.
3. **Fe cross-band stress test:** prioritize ungraded lines with a predicted
   weak/strong candidate and a verified observing path, then UV/reference-only
   gaps. Weak/strong remains a predicted routing hint; Fe mixing/correlation,
   relativistic/gauge convergence and independent withheld benchmarks must establish
   tractability. A recorded ungraded/Kurucz/VALD row does not prove no experiment exists.
4. **O II and M-dwarf applicability:** require an explicit target/ionization and
   diagnostic justification before spending solver time. Every M-dwarf row remains
   HOLD for target-specific depth, molecular blends, pressure broadening, magnetic
   sensitivity/Lande factors and model applicability. Solar usability is not a
   transfer certificate. OH/CO molecular diagnostics need their own provenance and
   equilibrium/model contracts; atomic structure calculations do not supply them.

## Reference documents, checks and self-review

Applicable local library documents inspected with pypdf include AGSS21,
Wiese & Fuhr 2006, Den Hartog 2019, and Ruffoni 2013. Existing source tables,
ASD pull provenance and experimental-table sidecars were inspected before any
external search; no held table or paper was downloaded again. Their bibliography
keys already exist (`asplund2021`, `wiese_fuhr2006`, `denhartog2019`,
`ruffoni2013`, `nist_asd` and the Fe source keys). No bibliography rows were added.
The scan is applicable-source inspection, not a successful whole-library audit.

| Command / control | Result |
|---|---|
| Generator above, with output `data/audit/rya1353_line_census/v4` | Exit 0; Python 3.12.13, macOS arm64; 49 hashed inputs and 16 artifacts |
| pytest census + shared line-match tests | Exit 0; 20 passed; wrong EP, missing EP, ambiguous components, equal-wavelength/different gf refusal, canonical row/value conservation, all input checksums and distinct Baker/Reiners holding spans |
| Fresh replay to `/private/tmp/rya1353-census-replay-v4` and direct byte comparison | Exit 0; all 16 artifacts byte-identical |
| `python -m pipeline.refinement_debt_join --report` | Exit 0; informational outstanding debt preserved |
| `python scripts/generate_sources_page.py --audit-library '<resolved refs.library>'` | Exit 1 on both this branch and clean baseline; `civis2013` verified state `vizier` and `mallinson2024` free-text verified state are not accepted by the validator; no complete library-coverage verdict available |
| `python scripts/audit_line_keys_rya1037.py --check` | Exit 1 on both branch and clean baseline; identical 31 AST findings across 12 files and excess waiver counts in `pipeline/nlte_corrections.py`; no finding in this census script |
| `python scripts/check_register_freshness.py --since-main` after commit | Exit 0; zero state surfaces changed; register bump not required |
| Ordinary `git push` pre-push ledger gate | Exit 1: current origin/main contains merged RYA-1349 / PR 611 without a SEQUENCE landing entry. The hook explicitly documents `--no-verify` for ledger-incomplete draft work; used to push this review branch, with the baseline failure disclosed. No merge/publication gate bypassed. |

Self-review: **CRITICAL gate** — full physical identity and exact conditioning are
unproven; all new candidates remain inadmissible. **WARNING** — grading conflict,
incomplete reference-code bibliography resolution, source-only gaps, incomplete
holding-span association and the two baseline audit failures are explicit.
**SUGGESTION** — add upper-level/J/component lineage to canonical identity, use
exact holding dispatch for the next observational audit, and version the per-line
independent benchmark split. Existing products and canonical gf are byte-preserved.
Run-bug ledger rows record the grading-policy and baseline-audit hazards discovered
in this preflight. No production run or merger is authorized by this report.

## Review gate and exact continuation

This dossier supplies the reproducible current-holdings denominator, alternatives,
negative selections, and an evidence-based review order. It does **not** satisfy
the stronger claim that every identity, every observing path, every primary
reference code, or every M-dwarf diagnostic has been scientifically validated.
Ryan must review the dossier; RYA-1352 must settle the conflicting classification
contract before any new class is assigned. Next work is source-level upper/J/HFS/
isotope adjudication of `source_records.csv.gz` and exact HoldingSpec association
of `holding_spans.csv`, starting with O I, then the withheld Fe experiments.
Missing experimental gf, missing any gf, and unavailable observation must remain
separate states throughout. No automatic publication or Linear Done transition.
