# Exoplanet Codex — Development Cycle

## Overview

All work on the Exoplanet Codex follows a formal development cycle. This applies to both Claude.ai (strategy/architecture) and Claude Code / Mr. Code (implementation). Neither AI skips stages or combines stages without explicit instruction from Ryan.

---

## The Cycle

### 1. PLAN
**Owner: Claude.ai**
- Scope is defined in Claude.ai before any code is written
- A Linear ticket (RYA-XXX) must exist before implementation begins
- Ticket includes: goal, acceptance criteria, files affected, data safety notes
- Mr. Code does not start work until Claude.ai has scoped it and Ryan has approved

### 2. IMPLEMENT
**Owner: Mr. Code (Claude Code)**
- Work from the Linear ticket spec exactly — no scope creep
- All new code goes in the correct pipeline module — never ad-hoc scripts in the root
- Pipeline outputs go to `data/results/` or `data/output/` only
- **`data/spectra/` is READ ONLY — see DATA_SAFETY.md**
- No deletions of any kind without explicit Ryan approval (see DATA_SAFETY.md)
- Commit messages reference the Linear ticket: `[RYA-XXX] description`

### 3. TEST
**Owner: Mr. Code, reviewed by Ryan**
- All new pipeline code is tested against **synthetic data or sample spectra first**
- Never run new/modified pipeline code directly against production FITS files on first run
- Test data lives in `data/test/` — small, synthetic, disposable
- Tests must produce a visible output (plot, CSV, printed summary) Ryan can inspect
- Solar calibration values are the benchmark: A(Fe)☉ = 7.46 ± 0.05 (Asplund et al. 2021)

### 4. QA
**Owner: Ryan (with Claude.ai support)**
- Ryan reviews test outputs before any production run
- Claude.ai checks results against known literature values where applicable
- Known validation benchmarks:
  - A(Fe)☉ = 7.46 ± 0.05
  - EW ranges: weak lines ~2–10 mÅ, strong lines ~50–200 mÅ
  - Gaussian fit residuals should be <5% of continuum
- QA failures go back to IMPLEMENT with a new or updated Linear ticket

### 5. PRODUCTION RUN
**Owner: Mr. Code, approved by Ryan**
- Ryan explicitly says "run on production data" — this is never assumed
- Outputs written to `data/results/[star]/[date]/`
- Original FITS files are never touched
- Run is logged with timestamp, pipeline version, and input file count

### 6. REVIEW & PUBLISH
**Owner: Ryan (with Claude.ai support)**
- Ryan reviews production results
- Claude.ai helps interpret, cross-check against literature, flag anomalies
- Only after Ryan approval does anything get committed to the repo or published to the site
- Linear ticket is marked Done after Ryan signs off

### 7. SESSION CLOSE — update the State Register
**Owner: Mr. Code, ratified by Ryan (RYA-516)**
- `CODEX_STATE_REGISTER.md` (repo root) is the **mutable current-truth ledger** — read it FIRST for "what is the current state of X", instead of reconstructing state from the ticket journal.
- Before posting any end-of-session Linear comment on a run that **settled / regressed / superseded** a component, **signed off a gate**, or hit a **milestone**, update the register per `skills/codex-state-register/SKILL.md`:
  - MIRROR rows (stellar params, abundances) are **regenerated from source** — `python scripts/gen_state_register_targets.py --write` — never hand-typed.
  - NATIVE rows (verdicts, statuses, gate states) are edited by hand with the establishing ticket cited; bump `Version:` + add a Changelog line.
  - Verify merge/integration state against **git**, not memory, before writing SETTLED.
- The end-of-session comment states which register rows changed (and flags any row a source contradicted).

---

## Hard Rules (Neither AI May Override)

1. No stage may be skipped without explicit Ryan instruction
2. `data/spectra/` is permanently read-only for all pipeline code
3. No deletions without Ryan's explicit "yes, delete X"
4. Test before production — always synthetic data first
5. One Linear ticket per piece of work — no untracked changes
6. Mr. Code commits reference the ticket ID in every commit message
7. **A gate cannot be signed off while its `CODEX_STATE_REGISTER.md` rows are stale** — the register is updated at every gate sign-off and the moment any component settles/regresses/is superseded (RYA-516; procedure in `skills/codex-state-register/SKILL.md`)

---

## Ticket States in Linear

| State | Meaning |
|-------|---------|
| Backlog | Scoped, not started |
| In Progress | Mr. Code is implementing |
| In Review | Ryan is reviewing output |
| Done | Ryan signed off |

---

## Sirius AI (Future — RYA-79)

When Sirius AI is online, it inherits all rules above. The autonomy model is:
- Sirius runs full pipeline measurement end-to-end
- Results are written to a review folder
- Ryan verifies results, then approves publication
- Sirius never deletes source data or overwrites previous results
