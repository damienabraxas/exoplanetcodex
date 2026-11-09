# RYA-1373 review dossier

**Draft data/schema package; independent scientific review owed.** Public atlas HTML
is sequenced after RYA-1370 approval, as the issue requires. No production science,
new abundance/grade, full-coverage claim, model adoption or publication performed.

## Scope and provenance

Science base `078cf7dedfa961452c8bd9af43c4bf83c5463fdf`; branch
`rya1373-frontier`, worktree `codex/rya1373`. Website Roadmap implementation is
PR [#87](https://github.com/damienabraxas/exoplanetcodex-site/pull/87), head
`8d39d3c`, with hosted Roadmap CI passing. This data package is independent of its
public-page approval. Other science/website worktrees were left intact.

The repository intake, implementation-brief, scientific-code-review and
model-applicability procedures were used. The abundance-product-audit skill was
inspected for applicability; no completed abundance product is being signed off.
Ledgers/science/product/convention/development authorities and canonical catalogs
were inspected on the same main base as RYA-1370. RYA-1373 has no comments; linked
1314/751/970 scopes and comments were inspected, along with the 1356/1368/1370
context already gathered. The source census/policies and actual CNO/Fe/molecular
artifacts control the facts, rather than issue status. Existing refinement debt is
not resolved by this inventory (tracker references include 548/470/475/683/697).

Canonical JSON editorial selectors generate 30 entries with exact source receipts,
typed components, source-specific transitions and descriptive axis envelopes.
O/Fe entries are supported by primary retained release descriptions, tables and
measured atom/training receipts. CNO and held molecular lists are distinct from
cool-star discovery pointers. Scientific applicability, component compatibility,
input completeness and reuse gaps remain visible. There are no atlas-validated
entries or authorized execution combinations.

## Source conflict / defect disposition

Primary CDS `ReadMe` identifies Amarsi 2019 table4 as Fe II 3D LTE, while the old
RYA-1035 decision text says 3D NLTE and speculates it may already be held. The local
provenance says table4 was not vendored. The atlas follows the primary source and
records the contradictory historical claim, without modifying the historical
source or asserting an affected abundance. A proposed source correction is in
`proposed_feii_source_correction.md` for review. One WARNING provenance row is
appended to the canonical bug ledger, status characterized, owned by RYA-1373.
The atlas regression rejects changing table4 to FULL3D_NLTE.

## Validation and environment

Metadata-only generation and synthetic/negative checks run locally; no RT, SE,
atomic-structure engine, fitting or source-spectrum operations occur.

- Isolated Python 3.9.6 + jsonschema 4.25.1:
  `python scripts/generate_frontier_atlas.py --check`: exit 0.
  `python -m unittest tests.test_frontier_atlas`: exit 0, 12 tests.
  Actual source hashes are compared with pinned `git show <base>:<path>` bytes.
- Existing project Python 3.12 environment:
  `python -m pytest -q tests/test_run_bug_ledger_rya926.py tests/test_molecular_canonical_rya1130.py`:
  11 pass, 1 fails on pre-existing empty science_impact cells (nine legacy rows,
  tickets 587/1250). The same rows are present in the base commit; see the baseline
  evidence file. The new row has an explicit impact. No unrelated history filled in.
- `python scripts/reconcile_sequence_vs_git.py --since-main`: exit 0.
- `python scripts/check_register_version_pointer.py`: exit 0.
- `python scripts/check_register_freshness.py`: exit 0, with the existing stale
  generated-tracker notice preserved. No scientific state/SEQUENCE bump is warranted
  for a draft inventory without gate advancement.
- Python compile and `git diff --check`: pass.
- Existing CI discovers new unittest cases through pytest. Hosted runner execution
  and scientific review remain distinct from the focused local results. Full science
  suite and heavyweight integrations were not run on the Mac.

## Limits and review gates

The bibliography library audit was already attempted in this session on the same
science base and fails before enumeration on the existing civis2013/mallinson2024
verified-state vocabulary errors; no complete library audit is claimed. This atlas
joins existing citation metadata without changing it. Exact separately named GitHub
Project Register/Target Star List/Canonical Element List were not positively identified;
canonical committed counterparts and that limitation were documented during RYA-1370.

Live primary-source checks included the Fe paper/upstream surrogate repository,
Turbospectrum paper, ExoMol index, H2O and SiO landing pages. Other attempted ExoMol
molecule/licence pages were inaccessible through the web tool; two DOI landing pages
returned 403. Primary held CDS and molecular source documents remain inspectable.
Discovery entries therefore retain null release-audit dates and pending license and
applicability fields. No inaccessible page is treated as evidence of absence.

Review first: physical identity and treatment, exact hulls and unassessed constraints,
source/component independence and historical conflict. Then review version/reuse
terms, completeness and source acquisition. Finally approve a reviewed atlas release
and decide public HTML integration after Roadmap approval. Until then, continue only
bounded source curation and preserve the draft/no-execution contract.
