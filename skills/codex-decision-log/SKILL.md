---
name: codex-decision-log
description: Where a ratified Exoplanet Codex decision gets recorded so it is not relitigated. Use whenever Ryan ratifies a decision, settles a science/architecture question, or changes a convention -- and before reopening any settled question. Routes each decision to the right home (owning-ticket comment per RYA-270, the register for conventions, memory for durable cross-session facts), flags supersessions, and requires checking the log before reopening. Kills the "we decided this already" relitigation loop.
---

# Codex Decision Log Skill

## Purpose

Half our friction is decisions that were *made* but not *lodged* anywhere findable, so they get relitigated sessions later (xi, the gate, the model count all took this hit). When a decision is ratified, record it immediately, in the right place.

## Where each kind of decision goes

- **Product / science / one-off decision** -> a COMMENT on the owning ticket (RYA-270), never a silent description rewrite.
- **Ratified convention / standard** -> the register (CODEX_STATE_REGISTER) or the relevant standards doc, AND a memory pin.
- **Durable cross-session fact** -> memory, so Claude carries it into the next session.
- **A decision that supersedes an earlier one** -> state "supersedes RYA-X / the earlier comment," and link it, so the old one is not re-followed.

## When a decision is made

- [ ] Which ticket owns it? Comment there (ASCII, green-light).
- [ ] Is it a convention? Register + memory pin.
- [ ] Does it supersede a prior decision? Say so explicitly and link the old one.
- [ ] If it involves a value, cite the source (single source of truth; no memory numbers).

## Before reopening a settled question

- [ ] Check the log first -- the owning-ticket comment, the register, memory. If it was decided, reopen ONLY with new evidence, and record WHY it is being reopened. "I don't remember deciding this" is not new evidence -- go read the log.

## Output

No artifact -- a discipline. The test: could a different session, or Mr Code, find this decision without asking Ryan? If not, it is not logged yet.
