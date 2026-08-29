---
name: codex-task-intake
description: The gate before starting any Exoplanet Codex task, answering any state/count/number/status question, or filing/updating any Linear ticket. Use this skill whenever Ryan says "start on X", "look into Y", "give Mr Code Z", "what did we decide", "double-check recent work", references "related tickets" or a value/count/status, or before creating any ticket. Enforces dedup-before-create, reading related tickets AND their end-of-session comments, loading live state from the authoritative surfaces (LEDGERS/register/SEQUENCE/live feed — NOT band_products or ticket-description tables), git-log-not-status merge truth, measure-from-artifacts-not-labels, and ticket labeling/linking discipline (title verb, parentId + relatedTo, priority, decisions-in-comments, ASCII, green-light rule). Its whole job is to stop repeated work and stale-number answers.
---

# Codex Task Intake & Ticket Discipline

## Purpose

We repeat work when we start blind — quoting a stale number, missing a ticket that already answered the question, or reconstructing state from memory. This is the gate before touching any task or filing any ticket. Two rules above all:

- **Never answer a state / count / number / status / "what did we decide" question from memory.** Check the authoritative surface first.
- **Never file a ticket without searching for the one that already exists.**

---

## Part 1 — Load context BEFORE you act (anti-repeat)

- [ ] **Dedup search.** Multi-term topic-cluster Linear queries (not one narrow term) for existing / related / duplicate tickets. Related work is usually 2-4 tickets, not one.
- [ ] **Read the related tickets AND their end-of-session comments.** *Done != merged.* The real result lives in the comment, not the status pill or the description. Scan the current recent cluster (the live 10xx range) — a sibling may already have done the work. (This session: RYA-311's measurement was already done and RYA-1093 was the open child; RYA-1083 had already fixed the "7.552" thing I almost re-investigated.)
- [ ] **Pull live state from the authoritative surfaces, in order:** LEDGERS.md -> CODEX_STATE_REGISTER -> SEQUENCE.md -> the **live product feed**. Do NOT read state from `band_products/` CSVs (gitignored, they lie) or from tables pasted into ticket descriptions (often stale by the time you read them).
- [ ] **Merge truth = git, not Linear.** `git log origin/main` after fetch. A "Done" branch/commit may be unmerged (paper-done trap).
- [ ] **Recall the task's memory + the last session's result** before proposing the next step.

---

## Part 2 — Ground-truth precedence (what to trust when sources disagree)

1. The **live feed / the actual committed artifact** — measured, this session.
2. **`git log origin/main`** — for merge state.
3. **Linear end-of-session comments** — for the real result of a ticket.
4. Ticket descriptions / status pills / memory — LAST, and only when nothing above contradicts them.

Measure from the artifact; never trust a label. Route, gf, atmosphere, and "departures applied" tags have all lied (RYA-1104).

---

## Part 3 — Ticket discipline (label + link so the next reader isn't lost)

- [ ] **Title verb:** BUILD / FIX / INVESTIGATE / SOURCE / AUDIT / REMOVE / RECONCILE / RUN — pick the one that names the work.
- [ ] **Link it:** `parentId` to the umbrella; `relatedTo` to every ticket the work reads from or affects. An unlinked ticket is a future duplicate.
- [ ] **Priority** set deliberately (Urgent only for a deadline or a hard block).
- [ ] **Decisions go in COMMENTS, not description rewrites** (RYA-270) — except ratified amendments to conventions.
- [ ] **ASCII-only in Linear comments** (Cloudflare WAF).
- [ ] **Firewall / provenance notes** where relevant: RYA-161 (never chosen to match a reference), single-source values with citations.
- [ ] **Green-light rule:** every Linear write needs Ryan's explicit go.

---

## Part 4 — Before you answer Ryan

- [ ] If the question touches a value, count, status, model, star, instrument, or a past decision — CHECK the source (feed / register / git / the ticket comment) before you answer. A confident wrong answer from memory costs far more than the 30 seconds to verify. If you catch yourself about to say a number from memory, stop and look.

---

## Output

No artifact — this is a gate, not a deliverable. If the intake surfaces a duplicate, a stale number, or an already-answered question, **say so up front and stop the wasteful path** before doing the work.
