---
name: codex-state-register
description: The maintenance discipline for CODEX_STATE_REGISTER.md — the mutable current-truth ledger that sits alongside the immutable Linear journal. Use this skill at every session close, at every gate sign-off, and the moment any component settles / regresses / is superseded, to update the register correctly (NATIVE vs MIRROR rows, status vocab, citation rule, version + changelog bump). Also use it whenever asked "what is the current state of X" — read the register FIRST rather than reconstructing state from the ticket journal. A gate cannot be signed off while its register rows are stale.
---

# Codex State Register — maintenance skill

## Why this exists
Reconstructing "what is the current state of X" from the Linear journal means reading ~40 tickets backwards across Done → regressed → superseded transitions — the archaeology class of failure (RYA-516). `CODEX_STATE_REGISTER.md` (repo root) is the **mutable ledger** that answers that in one read. The tickets stay the immutable **journal** (history); the register wins for **state**.

This skill keeps the register TRUE. If the register and reality diverge, the register is the bug.

## When to update (triggers)
Update the register — and say so in the end-of-session Linear comment — whenever:
1. **A gate signs off** — MANDATORY. A gate cannot sign off while its register rows are stale.
2. **A verdict changes** — a component settles, regresses, or is superseded.
3. **A milestone is reached.**
4. **A validated "X works / doesn't" finding lands.**

Standing duty: when a register-worthy moment lands in conversation, proactively flag *"this belongs in the register."*

## The rules (do not break)
- **NATIVE vs MIRROR.** NATIVE rows (verdicts, statuses, gate states, grid/model/instrument selections, reopen triggers) are hand-maintained here — the register IS their source of truth. MIRROR rows (stellar params ← `config/constants.py`; abundances ← results tables) are **script-generated, never hand-typed** — regenerate, do not edit in place.
- **Cite every value.** No value from memory. A row without its establishing ticket / source file is `[confirm]` and is NOT settled until the citation is attached.
- **Correct loudly.** If a ticket-of-record or `constants.py` contradicts a row, fix the row and call it out in the close comment — never silently overwrite.
- **Do not adjudicate science here.** The register records current ratified state; it
  does not compete with `SCIENCE_STANDARDS.md` or decision records. Preserve an
  unresolved conflict as `OPEN`/`PENDING`, cite both claims, and stop for disposition.
- **Status vocab only:** `SETTLED` · `SETTLED-WITH-CAVEAT` · `REGRESSED` · `STALE` · `OPEN` · `NOT-SELF-SUFFICIENT` · `PENDING`.
- **Column contract (verdict tables):** `component | verdict | value | established by | status | reopen-only-if`. Do not drop or reorder columns.
- **`reopen-only-if` is the anti-archaeology rule:** a SETTLED row is not re-litigated unless its named trigger fires.

## Session-close checklist
Run this before posting the end-of-session Linear comment on any ticket that changed a state:

1. **Did this session settle / regress / supersede anything, sign off a gate, or hit a milestone?** If no → nothing to do. If yes → continue.
2. **MIRROR rows:** if `config/constants.py` / `stars.yaml` / a results table changed, regenerate:
   `python scripts/gen_state_register_targets.py --write` (then `--check` must print OK).
3. **NATIVE rows:** edit the affected row(s). New value + `established by` ticket + `status` (from the vocab) + `reopen-only-if`.
4. **Verify merge/integration state against git**, not memory (e.g. `git log origin/main | grep -i rya-NNN`) before writing SETTLED vs SETTLED-PENDING-MERGE.
5. **Bump the header `Version:` and add a one-line Changelog entry** (git holds the full diff).
6. **`SEQUENCE.md` bumps alongside.** Every PR that touches the register must also add one line to `SEQUENCE.md` under a `## YYYY-MM-DD` heading. If the date section already exists, append to it. One line per landing, newest first, format `- **RYA-XXX** — one-sentence summary; what it unblocks`. Keep under ~140 chars. No CI guard on this — it is discipline, not machinery.
7. **Commit** with the ticket ID: `[RYA-XXX] update state register: <what changed>`.
8. **In the Linear close comment,** state which register rows changed and why (and flag any row a source contradicted).

## Regenerating the MIRROR (params)
`scripts/gen_state_register_targets.py` reads `STAR_PARAMS` (from `config/stars.yaml`) and rewrites the block between
`<!-- BEGIN GENERATED: targets ... -->` and `<!-- END GENERATED: targets -->`.
- `--write` refreshes the block · `--check` exits 1 if it is stale (CI-friendly) · no flag prints to stdout.
Never hand-edit inside those markers.
