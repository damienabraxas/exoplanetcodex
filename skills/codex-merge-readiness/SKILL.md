---
name: codex-merge-readiness
description: Pre-merge checklist for the Exoplanet Codex, run before Ryan (sole merge authority) merges any branch to main. Use whenever a branch is "done"/ready, Ryan asks "is this mergeable" or "can I merge this", or Mr Code reports a finished ticket. Confirms the Sirius run (not just Mac), rebase-then-set-diff, clean failure-set delta, end-of-session comment, branch+SHA / paper-done check, and the PR-vs-stacked route -- so "is this actually ready to merge" is a standard gate, not re-derived each time.
---

# Codex Merge Readiness Skill

## Purpose

Ryan is the sole merge authority; Mr Code never merges. Nothing goes to main until it clears this gate. "Done" in Linear is not "mergeable" -- this is what turns one into the other.

## Checklist

- [ ] **Sirius run done, not just Mac.** A Mac run with skipped tests is not sufficient for a production merge (two-box: Sirius is the production runner, RYA-567). If the report says "Mac run, N skips -- Sirius owed," it is NOT ready.
- [ ] **Rebased onto current `origin/main`, THEN set-diff.** An un-rebased set-diff is blind to tests that exist on main but not the branch (a stack was once 153 tests short and the failure-diff could not see it). Rebase first, diff second.
- [ ] **Failure set-diff vs clean main is clean:** identical failure set, zero NEW failures. Record the +passes.
- [ ] **End-of-session comment posted** on the ticket: branch name, commit SHA, result, caveats.
- [ ] **Register + SEQUENCE bumped IN THIS PR, verified by the freshness check.** HARD, NON-SKIPPABLE (RYA-1184). A merge is not merge-ready without it. Run all three and require exit 0:
      `python3 scripts/reconcile_sequence_vs_git.py --since-main` (missing set EMPTY)
      `python3 scripts/check_register_version_pointer.py` (RYA-690 pointer names the newest row)
      `python3 scripts/check_register_freshness.py` (RYA-659)
      Each SEQUENCE line comes from that ticket's own end-of-session comment + its merged diff -- **never from memory**, and never guessed. A landing you cannot source goes in the ticket's adjudication block, not into the log.
      Install the gate once: `git config core.hooksPath hooks` (BLOCKING `hooks/pre-push`; CI is off per RYA-954, so the hook is the only thing standing here).
      🔴 This is the item that failed: RYA-1213 merged 2026-09-16 with the register still reading PENDING / "review branch only", caught only by the RYA-1221 audit days later.
- [ ] **Merge state confirmed from git, not the status pill.** `git log origin/main` -- a "Done" branch may be unmerged (paper-done trap).
- [ ] **PR-vs-stacked route decided.** A rebased whole branch can go straight to main; a stacked branch rides its parents -- know which, and that the parents are handled.
- [ ] **No tuning; provenance / labels checked** (RYA-161); values unchanged where the ticket said they would be.
- [ ] **Ryan explicitly merges.** Green-light rule; Mr Code never merges.

## Verdict

**MERGE / HOLD (name the unmet item).** If Sirius is owed, the rebase-then-diff was not done, or the register/SEQUENCE bump is missing, the answer is HOLD -- no matter how clean the Mac run looked.
