# RYA-926 skills governance and workflow-gap audit

## Inventory and authority result

Reviewed all four pre-existing repo skills and the repository's canonical ledgers,
science/protocol/product/convention docs, specialized model/instrument/design/audit docs,
decision records, agent/development instructions, bibliography manifest, generator, and
local reference-library audit. `LEDGERS.md` is the smallest existing cold-start entry
surface, so it now owns the authority/read-routing, contradiction, and skills-governance
model rather than creating a competing entry document.

The legacy “Project Instructions for Claude” supplied in RYA-926 was treated as input
evidence. Durable scientific-method, provenance, skepticism, validate-don't-tune, and
review principles are already owned by science standards/protocol/product docs or are
routed there. Dated issue lists, file counts, paths, exact values, role names, and
instrument claims were not copied into skills. Legacy agent wording has no authority
over newer committed state.

## Existing skill dispositions

| Skill | Result |
|---|---|
| `codex-vald-extraction` | Preserved manual ownership, HFS-on, truncation/intake, split verification, quarantine/provenance, shared-parser, and no-silent-`log_gf` rules; added a tight authority boundary. |
| `codex-data-audit` | Re-scoped to scientific pipeline-entry safety. General provenance/frame/coverage/SNR methodology is separated from runtime-resolved instrument rules/config; BERV follows delivered-frame determination; coverage can yield DATA-GAP; checksums establish identity. |
| `codex-mr-code-brief` | Made agent-neutral and reconciliation-first; Linear mission no longer supersedes repo/science truth; code snippets are conditional; NaNs, hierarchy, milestones, instrument applicability, and complete handoff evidence are explicit. |
| `codex-state-register` | Preserves NATIVE/MIRROR/history discipline and now explicitly refuses to adjudicate unresolved science or compete with standards. |

## Candidate procedural-gap decisions

| Candidate | Disposition | Owner/rationale |
|---|---|---|
| Scientific cold start/context reconciliation | No new skill | `LEDGERS.md` is the canonical entry and now includes authority and contradiction routing; duplicating it would create two starts. |
| Scientific definition of done/evidence validation | No new skill | `ELEMENT_PROTOCOL.md` defines execution completeness; `SCIENCE_PRODUCT_PACKAGE.md` defines released evidence. |
| Model-domain/applicability audit | **New skill** | `codex-model-applicability` owns the reusable pre-run matrix across route/scale/model/atmosphere/gf and distinguishes unsupported, out-of-domain, and data gaps. Existing docs state facts but did not provide this procedure. |
| Reference literature/PDF intake | No new skill | `data/refs/README.md`, `bibliography.csv`, and `generate_sources_page.py --audit-library` already provide a maintained deterministic workflow; governance wording was strengthened. |
| Contradiction/escalation | No new skill | Cross-cutting behavior belongs at the `LEDGERS.md` authority entry, not in a competing catch-all skill. |
| Session handoff | No new skill | `codex-state-register` plus `DEV_CYCLE.md` owns state close; `codex-mr-code-brief` owns task evidence requirements. |
| General scientific-method review | No new skill | `SCIENCE_STANDARDS.md`, `ELEMENT_PROTOCOL.md`, and the SPP evidence rules already own it; another skill would duplicate law. |
| Scientific-software QA / code review | **New skill** | `codex-scientific-code-review` owns reusable wrong-data, provenance, cache-contamination, silent-fallback, and two-fixture review procedure. Existing science docs state law but do not provide this code-review method. |

## Committed defect-history result

No existing canonical committed bug/failure ledger was found. The nearest surfaces are
science-specific rejection/model-attempt ledgers and the local `debug/intake/*.jsonl`
trace; neither owns curated cross-ticket engineering/science defects. Therefore
`data/audit/run_bug_ledger.csv` is the single global ledger, with ownership/schema and
append/amend rules in `data/audit/RUN_BUG_LEDGER.md`. Session-close and brief workflows
require an explicit row disposition. The Kitt Peak→HARPS wrong-holding defect is seeded
from committed RYA-911/913 evidence; unrecoverable run/base/artifact facts are left
blank or `UNKNOWN`, not fabricated.

## Reference-library governance result

The existing manifest is `data/refs/bibliography.csv`; no parallel registry was created.
The audit initially reported `fe_recipe.drawio.xml` as uncited. That file is a project
workflow/design artifact, not literature, so `.xml` joins other non-reference asset
extensions in the audit rather than receiving fabricated citation metadata. Scientific
documents remain required to carry inspected provenance. Holding a paper never ratifies
its conclusions as project law.

## RYA-925 cold-start dry-read

Verdict: **STOP before the RYA-925 science run: cross-instrument behavioral canary owed.** A
cold agent can now identify authority, run data intake, enumerate every applicable model
path, use RYA-906 stored axes/derived display names, apply the element protocol, and
produce the required evidence/handoff without treating Linear or a skill as scientific
law. However, the exact Kitt Peak→HARPS regression class is not yet closed end-to-end.

RYA-906 is committed on main: route is resolved from handler/product evidence, not legacy
labels; `treatment` remains a permanent dual label; the near-UV Fe canary is synthesis;
and axis fields drive new reporting. RYA-925 must still resolve the live Al-specific line,
grid, holding, and implementation facts named in its own issue at execution time. Those
are required task inputs, not ambiguities in the governance stack.

Existing controls are real but incomplete: `tests/test_instrument_agnostic_rya913.py`
forbids instrument-specific flux readers in shared measurement routes, requires the
dispatch, forbids catalog instrument IDs as module-level constants, and carries positive
controls for both defect shapes. `scripts/derive_band_products.py` selects via
`instrument`, resolves a `holding_id`, records the holdings that served synthesis lines,
and fails if they differ from the holding selected for the requested instrument. The
distinguishing identity keys are therefore `instrument_id`, `holding_id`, band/window,
per-line wavelength/input set, route/model axes, and source provenance.

What is still missing is the comment's required behavioral proof: one test running two
deliberately distinguishable Kitt Peak/HARPS fixtures through the same public route and
asserting input identity → measurement identity → output identity, including that cached
measurements cannot collapse. Static AST checks and an in-run holding mismatch check do
not prove cross-run cache/state isolation. Until that canary exists and passes, a later
HARPS phase could regress to Kitt Peak data without the precise acceptance test RYA-926
now requires. File that implementation guard before executing RYA-925 or add it as a
bounded prerequisite; do not run Al science in this refactor.

RYA-925's eventual end-of-session report must also list every bug-ledger row created or
amended, or state `none discovered`, so the Al pilot's defects do not remain chat memory.

If later reconciliation finds conflicting model-domain evidence or an unratified
scientific choice, the authority protocol also makes it an `AMBIGUOUS-STOP` rather than
a guess.

No Al or Fe science was run or modified for this audit.
