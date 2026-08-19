---
name: codex-mr-code-brief
description: Write or review a Linear implementation brief for a coding or science agent. Use for any agent-executed ticket and always for pre-freeze, verification-gate, or phase-close work. Reconciles mission, current repo state, science authority, code, and refinement debt without duplicating mutable rules into the ticket.
---

# Agent implementation brief

## Why this exists

A brief is the durable mission and implementation contract in Linear. It must be
unambiguous enough for a cold agent, but it does **not** supersede current repository
state, ratified science standards, or committed implementation truth. Link canonical
sources instead of copying mutable values and instrument rules into every ticket.

Two failure classes motivate the rules below, and both have been paid for:

- **Orphaned refinement work.** RYA-524's science-refinement children (RYA-581 Ba
  deblend, RYA-585 Zr rescue, RYA-565 Eu adjudication) sat unfired through **eight**
  architecture tickets. Nothing in the briefing path ever asked "what is owed on the
  elements this brief touches?" The architectural children got executed because they were
  in the way; the science ones were not, so they were invisible. See RYA-676.
- **Stale premises.** RYA-672 was authored at 01:46Z asserting RYA-669 was "nearly
  complete"; 669 posted a STOP at 01:55Z. Four of 672's premises were wrong on arrival,
  including two command lines naming scripts that had never been committed. A brief is a
  claim about the repo, and claims decay.

## Pre-brief refinement debt check (mandatory)

Before drafting any pre-freeze, verification-gate, or phase-close agent brief:

1. Read `data/audit/element_status_tracker.csv` `refinement_debt` column.
2. If any element in the current phase has `Backlog:` in that column, list those tickets
   in the brief's Context section as **"orphaned refinement work to consider firing
   first"**.
3. If the parent ticket is Done and sibling Backlog tickets exist under it (check the
   parent's related issues in Linear), flag them explicitly.
4. If the brief is for a phase-close or freeze ticket AND any refinement debt is
   `Backlog:` for the current phase, add a section **"Refinement debt not yet resolved"**
   listing the tickets and the impact.

This check is not optional. Every brief must document either "no refinement debt in
current phase" (with the tracker read confirming) or list the debt.

```
python -m pipeline.refinement_debt_join --report        # informational, always exit 0
python -m pipeline.refinement_debt_join --phase-close   # exit 1 on open debt
```

Read the registry itself — `data/audit/element_refinement_registry.csv` — when you need
the *why*: it is the SSOT joining "this row is owed" to "this ticket resolves it", and its
header carries the admission rule. Two readings that are easy to get wrong:
`TBD - no resolving ticket` is the **loudest** row type (debt established, nobody has
filed a ticket — firing something cannot clear it), and an **empty cell means "no known
refinement path", not "nothing owed"**. Full discipline in `docs/CONVENTIONS.md`
§ "An owed row names the ticket that would fix it".

## Cold-start and verify the premises

Start at `LEDGERS.md` and reconcile applicable ledger/state rows, skills,
`SCIENCE_STANDARDS.md`, protocol/product/convention docs, current code/config/tests,
the issue description, and **all comments chronologically**. Check every factual claim
against fresh `origin/main` and record the base SHA. Cite `path:line` for code claims.
If a premise rests on an unmerged branch, name the branch and PR and state the contingency.

Grep **every** named precedent, not a subset. A partial sweep produces confident false
negatives: "no rchi2 ceiling exists on main" was written after grepping two of the five
harnesses that define one.

## Section contract

Number the sections; implementing agents reference them (`§3B`, `§4`).

1. **Context** — what surfaced this, with evidence. Include the debt check result.
2. **The defect / the goal** — stated as the *class*, not just the instance. If it is one
   of N, say so and enumerate.
3. **Spec** — complete behavior, interfaces, scope, constraints, and acceptance tests.
   Include code or pseudocode only when it materially constrains implementation. Each
   lettered sub-item is a decision to make or a thing to build. Where a
   judgement call is genuinely open, say *"establish X and report it"* rather than
   prescribing an answer you have not verified. **"Do not fabricate a uniform rule over
   things with genuinely different semantics; report the difference instead."**
4. **Do NOT** — the standing set, plus anything ticket-specific.
5. **Verification** — what evidence must exist for this to be done. Prefer
   before/after tables and a test that fails before the fix.
6. **End-of-Session Requirements** — branch/worktree, base SHA, final SHA/commit(s),
   every file changed, exact commands/results, errors, unresolved questions,
   collision/state checks, resume point, generated scientific artifacts where
   applicable, and PR link if authorized. State explicitly: no merge unless instructed.

Open the description with an agent-neutral launch instruction:

> **Agent brief — launch:** *"Read this issue's description and all comments
> chronologically, perform the repository cold start, reconcile current authorities,
> then implement."*

Use a true `parentId` only for actual work hierarchy. Use `relatedTo` for dependencies,
precedents, and history. Do not require a milestone when none applies. Instrument-specific
rules (APERO, BERV, tellurics, product columns, and similar) are applicability-gated and
linked from current instrument/science authority rather than pasted as universal rules.

## The standing Do NOTs

Copy these into every brief that touches results, and add ticket-specific ones:

- **Do NOT change any element's value or disposition** as a side effect of a contract or
  refactor ticket. If the correct fix would change what is emitted for any species,
  **STOP and report** — that is a science decision for Ryan, not a refactor.
- **Do NOT modify an existing frozen gold version or its hash manifest.** Resolve
  `data/reference/solar/CURRENT` and its current hash at authoring time. Gold is
  write-once; a ratified change creates a new version rather than editing the old one.
- **Do NOT tune toward literature.** Re-assessment is validate-don't-tune. A model
  changes because its physics or its provenance changed, never because the answer got
  closer to a reference value.
- **Do NOT merge.** Open the PR; Ryan merges.
- **Do NOT run grids on the Mac.** Resolve current compute and environment policy from
  the register and repo instructions; record the resolved runner and environment.
- **Do NOT `git add -A` after a pytest run** — the suite dirties a tracked file.

NaN behavior is semantic, never blanket masking: expected, scientifically valid missing
values may be represented and masked with an explicit disposition; unexpected or invalid
NaNs fail loudly with provenance. A brief identifies which class applies.

## Ratified constraints are not suggestions

If the brief asks an agent to enforce a rule Ryan has ratified, point it at
`pipeline/ratified_constraints.py` (RYA-674) rather than having it build a parallel check
in the module where the bug happened to bite. Adding a constraint needs a Ryan decision,
a registry entry citing the ratifying RYA-#, and a test — see
`docs/SCIENCE_STANDARDS.md` § "Ratified Constraints". An agent that *discovers* a rule of
that shape should **recommend** it and stop, not self-register it.

## Register and sequence discipline

A brief whose work changes a state surface must say so, and must require the register
bump and the `SEQUENCE.md` line **in the same PR** — see
`skills/codex-state-register/SKILL.md`. Run `scripts/check_register_freshness.py
--since-main` **after committing**; run before, and it compares an unstaged tree and
passes for the wrong reason.
