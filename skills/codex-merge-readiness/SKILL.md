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
- [ ] **Merge state confirmed from git, not the status pill.** `git log origin/main` -- a "Done" branch may be unmerged (paper-done trap).
- [ ] **PR-vs-stacked route decided.** A rebased whole branch can go straight to main; a stacked branch rides its parents -- know which, and that the parents are handled.
- [ ] **No tuning; provenance / labels checked** (RYA-161); values unchanged where the ticket said they would be.
- [ ] **Ryan explicitly merges.** Green-light rule; Mr Code never merges.

## Verdict

**MERGE / HOLD (name the unmet item).** If Sirius is owed or the rebase-then-diff was not done, the answer is HOLD -- no matter how clean the Mac run looked.
