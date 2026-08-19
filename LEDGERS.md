# Codex Ledgers — canonical read-set (read at session start, in order)

## Authority and read routing

This file is the canonical entry point for any agent working in this repository.
Read the seven ledgers below in order, then resolve task-specific sources before acting.

| Source | Owns | Does not own |
|---|---|---|
| `skills/` | Durable reusable procedure: **how** an agent performs work | Current values, task state, or new scientific law |
| Ledger/state chain below | Mutable current project, element, data, model, and workload state | Cross-cutting scientific method |
| `docs/SCIENCE_STANDARDS.md` and ratified `docs/decisions/` | Cross-cutting scientific rules and ratified decisions | Per-task mission or transient state |
| `docs/ELEMENT_PROTOCOL.md` | Element-analysis execution protocol and line-accounting ladder | Release-package format |
| `docs/SCIENCE_PRODUCT_PACKAGE.md` | Required scientific product and evidence package | Pipeline implementation state |
| `docs/CONVENTIONS.md` | Naming, schema, repository, and generated/native conventions | Scientific verdicts |
| Specialized model, instrument, design, NLTE, and audit docs | Conditional detail for the domain they state | Automatic authority outside their scope or vintage |
| `data/refs/bibliography.csv` and local reference documents | Inspected evidence and citation provenance | Project law merely because a paper is held locally |
| Linear description and comments in chronological order | Mission, acceptance criteria, holds/corrections, decisions, and history | Current repo state or authority to override ratified science |
| Current committed code, tests, and config | Executable implementation truth | Permission to contradict ratified science requirements |

For a task, read only conditional deep sources required by its instruments, models,
data, or products. A legacy agent instruction is contextual evidence, not authority;
resolve mutable values and workflows against the current canonical source before reuse.

### Contradictions: identify, preserve, and stop when necessary

Never silently reconcile conflicting sources. Record the exact claims and citations,
classify each source by the table above, check recency and committed/ratified status,
and preserve both in the audit or Linear report. Use the source that owns the kind of
fact only when the conflict is already dispositioned. If competing scientific rules,
provenance, units/frames, model applicability, or state cannot be resolved from a
ratified source, **STOP before measurement or publication and ask Ryan for a decision**.
A skill never overrides `SCIENCE_STANDARDS.md`, a ratified decision, or current state.

## Skills governance

All Exoplanet Codex agent skills are version-controlled project artifacts under
`skills/`. Review them like code/science infrastructure and check new or revised skills
for overlap and conflict with the authority map. Skills contain durable procedure and
safeguards. Resolve mutable project state, file counts, thresholds, paths, versions,
and numeric configuration from their canonical source at runtime rather than copying
them into a skill. Create a skill only for a reusable procedure, not as a second copy
of a document.

**If you are an agent, a Sirius-local model, or a collaborator opening this repo
cold: these seven files are what you read first.** They are the mutable state of
the project. Everything else in the repo is code, data, or history.

Until RYA-659 this list lived only in Claude.ai memory and chat, which is why the
State Register was able to drift eleven state-changing tickets without anyone
noticing. It is now a file.

| # | Ledger | Path | Source-of-truth for | Update trigger | Generated? |
|---|--------|------|---------------------|----------------|-----------|
| 1 | State Register | [`CODEX_STATE_REGISTER.md`](CODEX_STATE_REGISTER.md) | live state / verdicts / gates / grid + model + instrument selections / scope corrections | any gate sign-off, verdict change, milestone, or validated "we now know X" | **NATIVE** (hand) + MIRROR blocks (`scripts/gen_state_register_targets.py`) |
| 2 | Element status tracker | [`data/audit/element_status_tracker.csv`](data/audit/element_status_tracker.csv) | per-element status / tier / verdict / owed-work, 27 rows | any ticket changing an element's status | **GENERATED** (RYA-654) from phase_c + [`element_status_tracker_editorial.yaml`](data/audit/element_status_tracker_editorial.yaml) — `scripts/generate_element_status_tracker_rya654.py`; **do not hand-edit** |
| 3 | System catalog | [`data/catalog/system_catalog.csv`](data/catalog/system_catalog.csv) | star identity + pipeline lifecycle stage (13 systems) | a star advances a lifecycle stage | **NATIVE** — index-by-pointer; physics lives in `config/stars.yaml` |
| 4 | Instrument catalog (+ modes) | [`data/catalog/instrument_catalog.csv`](data/catalog/instrument_catalog.csv), [`instrument_modes.csv`](data/catalog/instrument_modes.csv) | instrument / mode capability + coverage (25 instruments, 11 modes) | a new instrument or mode use is verified | **NATIVE** |
| 5 | Holdings manifest registry | [`data/catalog/holdings_manifest_registry.csv`](data/catalog/holdings_manifest_registry.csv) | what data we already hold (anti-reinvent) | data acquired or verified for a system | **NATIVE** |
| 6 | Sequence log | [`SEQUENCE.md`](SEQUENCE.md) | narrative "what landed recently" overlay on the register | at every register bump, same PR | **NATIVE** — human-maintained, append-only, one line per landing |
| 7 | Element disposition report | [`docs/audit/element_disposition_rya663.md`](docs/audit/element_disposition_rya663.md) (+ machine twin [`data/audit/element_disposition_rya663.json`](data/audit/element_disposition_rya663.json)) | **can this element flip to PASS now, and if not what exactly is holding it** — the three ratified gates shown per element, plus stale-input evidence | any verdict / two-engine / gold change | **GENERATED** (RYA-663) — `scripts/gen_element_disposition.py`, `--check` in CI |

**Read #5 BEFORE proposing any download.** The holdings registry exists so we
stop re-acquiring data we already have.

## Standing session step — refresh the bibliography (RYA-854)

**The reference library gains documents between sessions, so check it every session:**

```bash
python3 scripts/generate_sources_page.py --audit-library "<reference library dir>"
```

It exits non-zero if any document in the library has no row in
[`data/refs/bibliography.csv`](data/refs/bibliography.csv), or if two rows name
byte-identical files. **An unlisted document is an uncited one.** Add the missing rows
(verify each against the document itself, never against its filename), then
regenerate the public page:

```bash
python3 scripts/generate_sources_page.py --site-root ../exoplanetcodex-site
```

This is *not* an eighth ledger — the bibliography is a citation catalogue, not mutable
project state, and it carries no register-bump obligation. It is on this page because
this page is what gets read first, and the check is worthless if nobody runs it. It
already caught three documents added mid-session on 2026-08-17.

**Read #7 BEFORE answering "why is element X still owed".** #2 tells you an
element's status; #7 tells you what is *holding* it, gate by gate. RYA-672 had to
reconstruct that by hand across three tickets because nothing routed a reader
here — which is precisely why it is now a read-set member (RYA-676).

## What each ledger is NOT

Naming convention (RYA-631) — these words are not interchangeable:

- **Catalog** = master enumeration, one row per entity (#3, #4, #5).
- **Register** = mutable current-state facts ledger (#1).
- **Tracker** = per-item work/progress status (#2).
- **Reference** = frozen, validated truth values — **"gold" lives ONLY there**,
  at `data/reference/solar/` (pointed at by `data/reference/solar/CURRENT`).

`SEQUENCE.md` is not a catalog/register/tracker/reference — it is a narrative log
that overlays the register; the four-noun convention still holds for state artifacts.

A value in the register or the tracker is *state*, not gold. Gold is frozen,
hashed and immutable; a measured-but-unfrozen value is recorded as such and
folds in at a ratified re-freeze (currently RYA-527, gold v3).

## Consistency enforcement

- **RYA-632** — results-ledger consistency guard: `pipeline/ledger_consistency_guard.py`
  cross-checks the tracker against the verdict artifacts (phase_c · gold ·
  physics_regime) and refuses SILENT disagreement. It is source-AGNOSTIC: it does
  not decide which artifact is right, so a stale tracker fails here — that is how
  "the tracker is updated on every merge" is enforced without a merge hook.
  Documented+ratified divergences live in `data/audit/known_verdict_divergences.yaml`;
  an entry without a `ratified_by` ticket is itself an error.
- **RYA-654** — tracker generator + `--check`: the committed tracker must equal a fresh
  regeneration, so hand-editing it is a build break. Status columns come from the phase_c
  verdict channel (ratified canonical); the analyst columns come from the editorial
  sidecar. The generator refuses to run off an uncommitted verdict artifact.

```
python scripts/generate_element_status_tracker_rya654.py --check
python -m pipeline.ledger_consistency_guard
```
- **RYA-676** — refinement-debt registry + join:
  [`data/audit/element_refinement_registry.csv`](data/audit/element_refinement_registry.csv)
  is the SSOT mapping **"this row is owed" → "this ticket resolves it"**, and
  `pipeline/refinement_debt_join.py` renders it into the tracker's generated
  `refinement_debt` column. It exists because the tracker said `owed`, Linear said
  `Backlog`, and nothing joined the two — which is how RYA-524's refinement
  children (RYA-581/585/565) sat unfired through eight architecture tickets.
  A cell reading **`TBD - no resolving ticket`** means the debt is established
  and **nobody has filed a ticket**; an **empty** cell means no known refinement
  path, which is *not* the same as "nothing owed". Adding a row requires a
  `provenance_ticket` — the admission rule lives in the registry's own header.

```
python -m pipeline.refinement_debt_join --report          # informational, always exit 0
python -m pipeline.refinement_debt_join --phase-close     # ESCALATED: exit 1 on open debt
```

- **RYA-659** — register-freshness guard: `scripts/check_register_freshness.py`
  loud-fails when a state surface is newer than the State Register, or when a PR
  changes state without touching it.

```
python scripts/check_register_freshness.py                    # history mode
python scripts/check_register_freshness.py --since-main       # PR mode
python -m pipeline.ledger_consistency_guard                   # ledger consistency
```

The authoritative list of state-changing surfaces (a superset of the ledgers
above — it also covers the verdict generator, the gold builder, `config/stars.yaml`
and `config/constants.py`) is single-sourced in
[`pipeline/state_surfaces.py`](pipeline/state_surfaces.py). Add a surface there,
not in a consumer.

**Both guards now run in CI** (`.github/workflows/ci.yml`, RYA-313/314 — a Sirius
self-hosted runner; all guard families run inside the `CI / test` job so the
required-check name never changes as guards are added). They remain runnable
standalone by the commands above. Enforcement differs, deliberately:

- register-freshness — **BLOCKING** (PR mode and push mode).
- refinement-debt — **INFORMATIONAL BY CONSTRUCTION** (RYA-676 §2C), not merely by
  workflow config: `--report` always exits 0, and the count carried in
  `pipeline.ledger_consistency_guard --json` cannot change that guard's exit code
  either. Seven registry rows are Engine-B `NO_MODEL_ATOM` gaps that need a model
  atom *acquired* — nobody clears those in a PR, and a permanently-red guard is a
  guard nobody reads. **Visibility is the deliverable, not gating.** A phase-close
  or freeze ticket (RYA-677) escalates it deliberately by running
  `--phase-close`, which exits 1 while any owed row in the phase still carries
  un-fired refinement debt. Do NOT flip the default.
- ledger-consistency — **`continue-on-error`, informational**, because it is
  PRE-DECLARED RED on six un-ratified `physics_regime` GET-DATA divergences
  (Co · N · P · Sc stale; K · Cu deliberate holds). Running it informationally
  surfaces any NEW contradiction without blocking every PR on the known red.
  It becomes blocking in the PR that closes **RYA-654**, by deleting one line
  of the workflow. Do NOT silence the red by annotating those elements into the
  exceptions file: an entry there requires a ratifying ticket, and annotating
  around an un-ratified contradiction is exactly the laundering this guard exists
  to stop.

## Maintenance discipline

Update a ledger row **the moment its component settles, regresses, or is
superseded** — and always at a gate sign-off: a gate cannot sign off while its
rows are stale. Every value cites its source; no value from memory. Full
procedure for the register: [`skills/codex-state-register/SKILL.md`](skills/codex-state-register/SKILL.md),
wired into session-close in [`DEV_CYCLE.md`](DEV_CYCLE.md). Tracker discipline:
[`docs/CONVENTIONS.md`](docs/CONVENTIONS.md).
