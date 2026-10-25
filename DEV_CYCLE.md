# Exoplanet Codex — Development Cycle

## Scope and authority

This cycle applies to people and coding agents, including Codex, Claude and Sirius
controllers. Start with [LEDGERS.md](LEDGERS.md) and its ordered read-set; its authority
map governs reconciliation. Scientific rules remain in
[SCIENCE_STANDARDS.md](docs/SCIENCE_STANDARDS.md) and ratified decisions. This cycle
specifies execution and human gates, not scientific values or acceptance thresholds.
Use [DATA_SAFETY.md](DATA_SAFETY.md) for source protection.

## 1. Discover and plan

Work needs a Linear issue with scope, behavior, constraints and acceptance criteria.
Read its description and every comment chronologically. Reconcile current Git, ledgers,
relevant project skills, code/config/tests and science authority before implementation.
Record worktree, branch and base SHA; check other workers' changes and jobs. Use
[the briefing skill](skills/codex-mr-code-brief/SKILL.md) to check refinement debt and
report the result, even when this task does not authorize resolving it.

Ryan's instruction to execute a scoped task authorizes its bounded implementation,
validation, commits, PR and required Linear handoff. A separate strategy-agent session
is not a prerequisite. Continue nonblocking substeps without repetitive approval.
An element-run ticket can impose stricter code-change checkpoints; read the current
[orchestrator runbook](docs/ORCHESTRATOR_RUNBOOK.md) and that ticket before a run.

## 2. Implement

Stay within the issue. Keep code in its owning modules and resolve configuration,
paths and scientific inputs from their canonical owners. Never silently fall back,
suppress unexpected errors, tune toward reference values, or manufacture provenance.
Source spectra/FITS and frozen gold remain protected; preserve previous results and
other workers' files. Commit only intended files, with the issue ID in each message.
Do not stage the whole tree after tests: some tests rewrite tracked artifacts.

## 3. Validate

Test changed behavior on synthetic or disposable fixtures before production inputs.
Run applicable automated checks and required regression tests. For science changes,
also check numerical/physical validity, units/frames, model coverage, identity and
artifact/provenance completeness through the owning protocols and product contract.
Resolve benchmarks and thresholds from canonical configuration or inspected literature;
never copy them from historical development prose.

Report exact commands, exit codes, interpreter/host, results, skips and failures.
A passed process is not scientific validation. Compare pre-existing failures with a
clean current-main baseline when needed; do not hide new failures behind the baseline.
Check generated artifacts without rewriting unrelated products. CI definitions in Git
are not proof that the service or runner is enabled; verify actual check status.

## 4. Review

Self-review the final diff using
[scientific code review](skills/codex-scientific-code-review/SKILL.md), with
[model applicability](skills/codex-model-applicability/SKILL.md) and
[product audit](skills/codex-product-audit/SKILL.md) when applicable.
Ryan reviews scientific test outputs before production and signs off scientific
results. Preserve and cite unresolved conflicts; never adjudicate science by editing
an agent instruction or a state row.

## 5. Human gates and production

Stop the dependent action for unresolved ratified scientific conflicts, indispensable
missing data/access, destructive changes, scientific product changes, production,
publication or merge requiring Ryan's authorization. Explain the exact missing
permission or evidence; continue independent authorized work. Existing explicit
session authorization need not be requested again. Scope expansion needs a decision.

Production requires Ryan's explicit instruction to run on production data; a validation
command or a ticket's historical example is not that instruction. Resolve compute and
environment policy from the current science standards, register and host profile.
Protect original FITS, log input identity and generator version, and preserve outputs.
For an element run, `--publish` includes budget assembly, feed admission and withdrawal
of superseded products: it is a scientific product action, not a harmless report flag.
Ryan authorizes publication and merges PRs. Never mark Linear Done before his sign-off.

## 6. Evidence and handoff

Complete authorized work through discovery, implementation, validation, artifact checks,
self-review, PR and the issue's Linear handoff. Include branch/worktree, base and final
SHA, every changed file, exact commands/results, errors, blocked resources, open risks,
PR URL and exact resume instructions. Report the run-bug-ledger disposition under
[data/audit/RUN_BUG_LEDGER.md](data/audit/RUN_BUG_LEDGER.md).

For long runs checkpoint the branch/SHA, host/environment, input manifests, command,
log/report paths, owned job identifiers, completed and outstanding steps, blockers and
next command. Inspect current load before starting heavy work. Never kill another
worker's processes or use a host-wide pattern kill. A checkpoint or partial run is not
completion; the evidence must distinguish attempted, validated and still owed work.

## 7. Session close and state

Apply [the state-register skill](skills/codex-state-register/SKILL.md) when its update
triggers fire. Update affected NATIVE rows with cited evidence; regenerate MIRROR
blocks from source. Verify integration against Git before claiming settled state.
A register bump requires its matching SEQUENCE line in the same PR. Do not bump
scientific verdicts for a documentation-only change.

Run register freshness checks, including PR mode after committing. Run the bibliography
library audit each session per LEDGERS; use the resolver and citation procedure in
[data/refs/README.md](data/refs/README.md). Missing or unindexed resources are reported
precisely, never concealed by a passing claim.

| Linear state | Meaning |
|---|---|
| Backlog | Scoped, not started |
| In Progress | Implementation/validation underway |
| In Review | Evidence and PR ready for Ryan's review |
| Done | Ryan signed off |

## Standalone software invariant

The orchestrator and engines remain directly runnable from documented CLI commands,
with explicit environment, manifests and deterministic inputs/outputs. Agents may plan,
invoke, inspect and report; radiative transfer, abundance algorithms, uncertainty
propagation and publication gates must not depend on an agent or model service.
