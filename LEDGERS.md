# Codex Ledgers — canonical read-set (read at session start, in order)

**If you are an agent, a Sirius-local model, or a collaborator opening this repo
cold: these five files are what you read first.** They are the mutable state of
the project. Everything else in the repo is code, data, or history.

Until RYA-659 this list lived only in Claude.ai memory and chat, which is why the
State Register was able to drift eleven state-changing tickets without anyone
noticing. It is now a file.

| # | Ledger | Path | Source-of-truth for | Update trigger | Generated? |
|---|--------|------|---------------------|----------------|-----------|
| 1 | State Register | [`CODEX_STATE_REGISTER.md`](CODEX_STATE_REGISTER.md) | live state / verdicts / gates / grid + model + instrument selections / scope corrections | any gate sign-off, verdict change, milestone, or validated "we now know X" | **NATIVE** (hand) + MIRROR blocks (`scripts/gen_state_register_targets.py`) |
| 2 | Element status tracker | [`data/audit/element_status_tracker.csv`](data/audit/element_status_tracker.csv) | per-element status / tier / verdict / owed-work, 27 rows | any ticket changing an element's status | NATIVE today; **GENERATED from phase_c is RYA-654** (not yet built) |
| 3 | System catalog | [`data/catalog/system_catalog.csv`](data/catalog/system_catalog.csv) | star identity + pipeline lifecycle stage (13 systems) | a star advances a lifecycle stage | **NATIVE** — index-by-pointer; physics lives in `config/stars.yaml` |
| 4 | Instrument catalog (+ modes) | [`data/catalog/instrument_catalog.csv`](data/catalog/instrument_catalog.csv), [`instrument_modes.csv`](data/catalog/instrument_modes.csv) | instrument / mode capability + coverage (25 instruments, 11 modes) | a new instrument or mode use is verified | **NATIVE** |
| 5 | Holdings manifest registry | [`data/catalog/holdings_manifest_registry.csv`](data/catalog/holdings_manifest_registry.csv) | what data we already hold (anti-reinvent) | data acquired or verified for a system | **NATIVE** |

**Read #5 BEFORE proposing any download.** The holdings registry exists so we
stop re-acquiring data we already have.

## What each ledger is NOT

Naming convention (RYA-631) — these words are not interchangeable:

- **Catalog** = master enumeration, one row per entity (#3, #4, #5).
- **Register** = mutable current-state facts ledger (#1).
- **Tracker** = per-item work/progress status (#2).
- **Reference** = frozen, validated truth values — **"gold" lives ONLY there**,
  at `data/reference/solar/` (pointed at by `data/reference/solar/CURRENT`).

A value in the register or the tracker is *state*, not gold. Gold is frozen,
hashed and immutable; a measured-but-unfrozen value is recorded as such and
folds in at a ratified re-freeze (currently RYA-527, gold v3).

## Consistency enforcement

- **RYA-632** — results-ledger consistency guard: tracker vs verdict artifacts
  cross-checked, contradictions unrepresentable at merge. *(built, not merged)*
- **RYA-659** — register-freshness guard: `scripts/check_register_freshness.py`
  loud-fails when a state surface is newer than the State Register, or when a PR
  changes state without touching it.

```
python scripts/check_register_freshness.py                    # history mode
python scripts/check_register_freshness.py --since-main       # PR mode
```

The authoritative list of state-changing surfaces (a superset of the ledgers
above — it also covers the verdict generator, the gold builder, `config/stars.yaml`
and `config/constants.py`) is single-sourced in
[`pipeline/state_surfaces.py`](pipeline/state_surfaces.py). Add a surface there,
not in a consumer.

**There is currently NO CI harness in this repo** (no `.github/workflows/`, no
Makefile, tox or pre-commit config), so both guards are runnable standalone and
are not yet enforced automatically. Wiring them in is a flagged follow-on shared
by RYA-632 and RYA-659.

## Maintenance discipline

Update a ledger row **the moment its component settles, regresses, or is
superseded** — and always at a gate sign-off: a gate cannot sign off while its
rows are stale. Every value cites its source; no value from memory. Full
procedure for the register: [`skills/codex-state-register/SKILL.md`](skills/codex-state-register/SKILL.md),
wired into session-close in [`DEV_CYCLE.md`](DEV_CYCLE.md). Tracker discipline:
[`docs/CONVENTIONS.md`](docs/CONVENTIONS.md).
