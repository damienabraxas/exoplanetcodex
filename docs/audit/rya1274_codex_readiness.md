# RYA-1274 — Codex readiness and standalone ORCH audit

## Scope, baseline and conclusion

Review branch: `ryandamienschmitt/rya-1274-codex-readiness`.
Base: `183310d8851629f43549e93269989032aa422bae` (fetched main, 2026-10-25).
The isolated `rya1274` worktree leaves other workers' worktrees untouched.
No scientific run, data download, abundance/product/gold change, production action,
publication or merge occurred. No engine or numerical algorithm was edited.

The thin bootstrap and local-reference discovery are validated. This is **not a
certification of complete ORCH/science readiness**: synthetic probes reveal completion
and routing defects, recorded below and in bounded child tickets. Existing baseline
test failures and library gaps remain explicit. See [inventory](rya1274/inventory.md)
for all skills, instruction paths, hooks, CI, runbooks and integration limits.

## Before and after

| Before, established from the base tree | This PR |
|---|---|
| No root or nested AGENTS/CLAUDE adapter in Git | Concise `AGENTS.md` routes to the existing authority/read-set and procedures. |
| DEV_CYCLE assigns mandatory strategy/implementation roles to Claude variants and copies historical numerical benchmarks | Agent-neutral discovery through review/PR/handoff; values and thresholds resolved from owning science sources; Ryan's scientific, destructive, production, publication and merge gates preserved. |
| Library described in bibliography prose but absent from path register | `refs.library` uses the existing resolver, an existing workstation holding and a host override; no library created/copied/downloaded. |
| LEDGERS prose says seven ledgers while table contains more, and points to sequence for gate blockers | Read in table order; gate-blocker routing names the disposition report. No state/count snapshot copied into the bootstrap. |
| Runbook calls prepare read-only, treats code as overriding policy and implies broad local-run/CI assumptions | States report writes, scientific authority boundaries, actual invocation limits and distinction between run/coverage/scientific completion. |

No new skill was created. RYA-926's nonduplication test is satisfied by reusing
LEDGERS, existing review/model/product/state skills, bibliography tooling and the
existing development/runbook surfaces. The run-bug ledger owns encountered defects.

## Findings with evidence

Line numbers below refer to baseline code unless the finding names a changed document.

| Severity | Evidence | Finding and disposition |
|---|---|---|
| CRITICAL | `pipeline/coverage_report.py:95` | `stale()` catches any resolver/hash exception and returns the empty current-state reason. Synthetic dependency failure produced PUBLISHED=1, STALE=0, done=true. **RYA-1275**, preserve diagnostic and fail closed before claiming currency. No real feed was touched. |
| WARNING | `pipeline/element_publish.py:53` | Budget unit expansion excludes the `reference` pool even though coverage enumerates it (`pipeline/coverage_report.py:104`). Mocked call confirmed set/codex/deep only. **RYA-1276**, verify complete routing with distinct fixture selectors. |
| WARNING | `run_pipeline.py:338` | Sweep branch accepts but does not forward publish/jobs; matrix+publisher composition exists only in single-element branch. **RYA-1276**, honor supported options or refuse explicitly. |
| WARNING | `scripts/generate_sources_page.py:531,583` | Actual held-library audit fails: 40 unindexed documents, 12 declared local_file paths unresolved under the existing audit rule. **RYA-1277**; missing paths captured in [resource report](rya1274/missing_references.json). Not every diagnostic is a new paper: duplicates/companions/internal notes require inspection. |
| WARNING | `tests/test_rya1235_post_run_report.py:21`, `data/catalog/linear_dossiers.csv:1` | Zn canonical target has no report dossier route. Fails identically on clean base. Existing RYA-1235 owns reporting; this PR does not invent a dossier assignment. |
| WARNING | `tests/test_run_bug_ledger_rya926.py:41` | Existing RYA-587 row has empty science_impact. Baseline and branch fail identically; appended RYA-1274 rows populate required fields. Existing historical row is not silently rewritten. |
| WARNING | `.github/workflows/ci.yml:69` | Committed CI still materialises LFS raw VALD files removed by RYA-1228. GitHub API reports workflow active; neither passing runner availability nor successful fresh checkout follows. No workflow/download triggered here. |
| WARNING | `docs/ORCHESTRATOR_RUNBOOK.md:7`, `docs/SCIENCE_STANDARDS.md:265` | Historical Mac production examples conflict with the general Sirius compute rule unless a task-specific authorization applies. Runbook now routes to current authority/ticket; this PR ratifies no compute exception. |
| SUGGESTION | `skills/codex-vald-extraction/SKILL.md:87`, `skills/codex-decision-log/SKILL.md:20` | Legacy raw-location and Claude-memory wording remains in focused skills. Resolve current resource/authority owners; memory cannot become canonical science. No broad skill rewrite needed for this adapter. |

The coverage probe replaces feed, ledger, descriptor and resolver in a temporary tree:
one synthetic published product, one recorded old hash, resolver raises
`RuntimeError("synthetic unavailable dependency")`. Observed counts and empty reason
are the implementation's output, not an agent's readiness self-rating. The pool probe
mocks `run_matrix.expand` and records `units_for`'s requested pools. Both are bounded
engineering evidence, not numerical/physical validation of a real element.

## ORCH understanding: source, lifecycle and integration

Read RYA-1232/1233/1250/1251 descriptions and all returned comments chronologically,
plus RYA-926 and its comments. Git history, not Done labels, established integration.
The [runbook architecture map](../ORCHESTRATOR_RUNBOOK.md) points to the current owners.

1. **Inputs and preparation.** `config/stars.yaml` via `get_star_params` owns pinned
   parameters. Catalogs/holdings and `HoldingSpec` define instrument/window/readers;
   raw or unwired holdings are reported separately. `element_prepare.prepare` audits
   literature, published selections/gf, cull candidates, ion-aware model coverage and
   rest-frame holdings. It writes reports. `apply` changes canonical inputs and culls;
   that remains an explicit human checkpoint. Missing NLTE rows and missing gf sigmas
   do not license fabricated corrections or complete budgets. Structured holds are
   scoped by species/window; unknown/malformed reports fail closed.
2. **Axes and route.** `run_matrix.expand` derives targets, ions, bands, corrected
   holdings, engine decks and pools from owning registries/policy. `RunDescriptor`
   and `resolve` carry species/ion/holding/window/deck/method/pool/star and produce
   preconditions plus exact stage arguments. Route is independent of treatment.
   Published-set and full-lab pools both belong to Reference Grade; Codex/Deep are
   distinct. Per-deck output directories prevent overwrites. Model/atom ion reach
   must support the requested physics, not merely contain an ion ground level.
3. **Execution and workload.** Planning is serial; stage execution uses a bounded
   thread pool. A step runs as a Python subprocess with explicit environment and
   interpreter checks. Expected holds/blockers, missing host assets, actual failures,
   dry intentions and built-unpublished work retain different statuses/reasons.
   There is no host-wide CPU scheduler or automatic Sirius routing. Observe current
   jobs/load before authorized work and never kill another session's processes.
4. **Fingerprints and resume.** `input_fingerprints`, `inputs_hash`, `ledger_key` and
   `is_current` record named inputs: descriptor/method, stage scripts, holding
   manifest/source bytes, selected-line recipe, scoped canonical gf/synthesis/NLTE
   rows, registry and cull inputs. The ledger records artifacts with hashes so
   unchanged unpublished work can SKIP. Published hand-campaign adoption is explicit
   in `run_matrix.run` (RYA-1232) and must not be mistaken for a fresh synthesis.
   Progress/latest reports and build fingerprints are evidence for a resume, not
   scientific sign-off. This audit does not prove a complete transitive dependency
   fingerprint or invalidate every possible cache; the synthetic tests cover their
   stated identities only.
5. **Budgets and admission.** Single-element `--publish` composes built synthesis
   units -> paired nominal/xi/profile/continuum/atmosphere/flux-response legs -> shared
   assembler -> deck emission and hold filtering -> `publish_product` gates.
   Source-correlated gf errors come through `gf_error_model`; paired responses must
   retain pool identity and recipe aggregation. Actual correction provenance
   distinguishes direct transfer, grid and cited solar-table deltas. Missing/failed
   legs are not zero uncertainty. Unit signatures and DONE stamps govern reuse;
   stage/superseded history keeps feed provenance inspectable. Publication may
   withdraw old SYNTH records; a narrowly filtered publish invocation can therefore
   have element-wide consequences. It is not a harmless reporting command.
6. **Verdict, coverage and reporting.** Literature verdicts diagnose outcomes and
   problem lines; they do not tune results. Coverage joins the feed, latest report,
   budget holds, quarantine and current input hash into PUBLISHED/HELD/NOT_RUN/STALE.
   DONE means no unrun/stale cells in that accounting and can contain holds. The
   exception-to-current defect above prevents unconditional trust in that claim.
   `post_run_report --dry-run` previews a dossier report; a real post is an external
   write. Science feed publication, generated website PR/deploy and Ryan merge are
   separate authorization boundaries.
7. **Sweep.** `run_sweep` derives the master target axis, prevents duplicate Fe-ion
   work, executes Fe first, gates dependent elements and reports omitted targets
   explicitly. It calls `run_matrix`, not `element_publish`; do not assume
   `--all-elements --publish` is a complete publication sweep on this baseline.
   The unflagged whole-star chain is explicitly legacy, not the current matrix.

### What the actual development/testing tickets established

- **RYA-1232** changed scientific inputs/conditioning and complete-budget campaigns;
  its description is historical and its only returned comment is an early RYA-1233
  handoff. Later behavior is established by merged code/artifacts and related comments,
  not the old fixed product counts or a Done pill. Current adoption, full-coverage
  corrected holdings and experimental model flags are inspected in their owners.
- **RYA-1233** evolved from MVP/front door through raw-holding exclusion, route/pool
  identity, per-deck paths, governing-process checks, unchanged-build reuse, budgets,
  superseded withdrawal and provenance preservation. Its early comment that budgets
  do not exist is superseded by `element_publish`. The latest comment still names
  manual NLTE extension, host routing/environment and workload gaps. Preparation and
  culling behavior must be read from the actual branch; no unconditional autonomous
  production/cull permission is granted here.
- **RYA-1250** exposed neutral/ionized grid confusion, departure engagement/reset,
  broadening and feature-component identity, cited-table vs full-3D distinctions,
  aggregation/provenance propagation, scoped holds and honest failed inversions.
  Its reopen comments correct premature completion/absent-feed claims and earlier
  erroneous gf hypotheses. Current Mg-fix branch adds canonical/cluster repairs,
  line-list/label sidecars and per-process Gerber workdirs. These are **unmerged
  changes**, not capabilities silently borrowed into this PR. RYA-1269 already owns
  the cross-element race audit; RYA-879/688 own concurrency/scratch coordination.
- **RYA-1251** adds subsequent shared routing/hold fixes, an NLTE abundance-axis
  bracket, failed-synthesis retries and HFS/telluric budget evidence. Its final reported
  holds are convergence responses, not absent continuum preparation. The branch is
  unmerged against this baseline; no result copied here and no scientific claim
  adjudicated. Mg/Na shared fitter changes must be reconciled by their owners.

The tests used here exercise orchestration/resolution/coverage/reporting and the merged
Mg ion/solar-table code without invoking scientific fits. No agent SDK/service is
imported by the CLI, matrix, descriptor, sweep, preparation, budget/publish or numerical
paths inspected. Direct CLI help executes successfully. This verifies a standalone
software interface and synthetic behavior, not complete production capability or
all model applicability on Sirius.

## Reference evidence and cold behavior

The path register now owns `refs.library`, with optional `CODEX_REFERENCE_LIBRARY`.
Existing home-root resolution is reused; another host must override or report missing
access. `data/refs/README.md` explains availability, parent-relative bibliography paths,
readable formats, missing-resource reporting and page/table provenance. It does not
copy documents or add a competing citation index.

The missing-path report also records filename candidates physically present in the
library. These candidates have not all been citation-matched: an unresolved declaration
is not proof that the paper is absent, and no download is warranted from this report.

A separate cold Codex agent started with no conversation history at AGENTS.md. Its
observed tool reads located the parameter owner, actual bibliography record and held
Asplund PDF, Table A.2, `FE_AGSS21` declaration and model capability accessor. The
primary session also ran deterministic missing-path and routing probes. The cold
probe read the bibliography's `Apslund 2021.pdf`: title/authors/year/DOI matched
`asplund2021` (10.1051/0004-6361/202140445). Table A.2 is PDF index 27, displayed page
28, printed A141 page 28 of 30. The primary session independently re-opened the title page and Table A.2 and verified
SHA-256 `bcd9e4527eec239df0f44ee3f2f64a7757896a8cf18eac1c343478489e959938`.
The runbook alias `aa40445-21.pdf` has identical bytes
in the actual holding; that alias remains an audit discrepancy, not a new citation.

Negative behavior:

- Missing library override raises FileNotFoundError naming the logical key and no
  fallback. Missing/unindexed synthetic evidence fails the existing library audit.
- A synthetic legacy instruction asking tuning is rejected under LEDGERS' authority
  protocol and SCIENCE_STANDARDS' reference-proximity prohibition. Preserve both
  claims and stop unresolved dependent science; no policy-reconciliation code invented.
- Actual `nlte_ion_capability("Fe", "II")` returns false with the atom's lack of
  bound-bound transitions, while Fe I returns true. This is a deck/ion-specific
  unsupported result, not a claim that all Fe II NLTE is unavailable.

The cold probe is targeted, not a claim of exhaustive rereading or a second historical
agent experiment on the old tree. Before/after discovery is established by Git absence
of AGENTS and the new observed routing. All scientific claims remain subject to their
own model/data protocols.

## Partial-task/stop classification

Inspected local session event metadata for the Mg session created 2026-10-18 and
read its Linear history. Task-complete events on 2026-10-23 coexist with the ticket's
subsequent reopening; event labels do not prove science completion. No explicit
permission-denied/timeout/termination output marker was found in that limited scan;
this is **not evidence that no such event occurred elsewhere**. Secret-bearing app
config and raw session logs are not committed.

| Evidence | Supported classification | Limit |
|---|---|---|
| RYA-1250 early unpublished-source/broadening/model-domain holds | Legitimate science hold | No guessed completion or replacement input. |
| RYA-1250 reopened after a missing-main-feed/incorrect-pool review | Context/completion-evidence gap | A task-complete event alone cannot explain the root cause. |
| RYA-1233 comments report timed-out Si cells and incomplete legs under load | Timeout/resource contention, as reported in ticket | Historical process/host state was not independently reconstructed. |
| RYA-1250 Gerber shared-workdir crash and per-process fix on unmerged branch | Resource contention with documented mechanism | No affected-element rerun authorized in this ticket; existing RYA-1269 owns it. |
| This session's get_agent_context call returns tool-not-found | Tool access limitation | Issue/comments remain accessible; discovery continued. |
| GitHub CLI initially failed network access, then approved retry succeeded | Sandbox/network access | Not classified as a scientific blocker. |

## Validation and acceptance limits

Interpreter: existing project `codex/venv/bin/python`, Python stack unchanged; Mac host.
`source config/hosts/mac.env` used for existing local holdings. Exact command groups and
outcomes are in [validation](rya1274/validation.md). No heavy Sirius run was scheduled.

Refinement report: exit 0, informational. It reports outstanding Solar Beta debt:
10 ticketed not-fired rows and 5 TBD rows, with snapshot ticket states; not a live
Linear claim. Backlog entries include RYA-523/683/475/548/470/697; no element-science
work is authorized here. The old registry's Mg debt predates current Mg work and is
not silently rewritten. This is not a freeze/phase-close ticket.

State update: a NATIVE readiness row records these validated findings as OPEN, v172
and a matching SEQUENCE review-branch line. No stellar mirror/scientific verdict
changed. Register freshness and target mirror checks run; inherited tracker-generation
notice remains informational. Four RYA-1274 defect rows appended with owning follow-ups.
The pre-existing malformed historical row and reporting failure are disclosed, not
laundered into a green result. Review must account for baseline failures; no full-suite,
Sirius or production validation is claimed.

Resume: review this PR; RYA-1275 fixes coverage currency, RYA-1276 repairs pool/CLI
routing, RYA-1277 reconciles library diagnostics. Reconcile with current main and active
Mg/Na branches before code work. Never restart production or merge from this handoff.
