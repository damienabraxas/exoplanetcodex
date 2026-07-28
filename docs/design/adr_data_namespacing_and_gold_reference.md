# ADR — Per-star output namespacing + the frozen, versioned gold-standard solar reference

**Status:** Accepted (RYA-469, 2026-06-29). Standing convention; future briefs inherit it.
**Peer:** [`docs/CONVENTIONS.md` — artifact preservation (RYA-461)](../CONVENTIONS.md).
**Law it encodes:** *don't fuck with gold-standard data.*

## Context

`solar_abundances.csv` was a **bare, generic filename**, tracked-despite-gitignore under
`data/processed/`. Two failure modes:

1. Any incidental `abundances_derive` run showed it modified, so a **perturbed regen could
   be committed as "the baseline"** (the RYA-467 footgun: tests asserting a value against a
   file anyone could silently rewrite).
2. As Procyon / α Cen / 55 Cnc / Sirius begin producing abundances, they would **collide on
   the same generic names** (`solar_abundances.csv`, `solar_per_line.csv`, …).

The solar verdict is the **gold-standard differential denominator** every benchmark is
measured against. It must be immutable, versioned, and physically un-collidable.

## Decision

### 1. No bare per-star filenames — the star is in the PATH

Every per-star pipeline product is written under a per-star directory:

```
data/outputs/{star}/{star}_abundances.csv
data/outputs/{star}/{star}_per_line.csv
data/outputs/{star}/{star}_ew_integrity.csv
data/outputs/{star}/{star}_verdict.json
data/outputs/{star}/diagnostics/...          # e.g. fe1_triage_quarantine.csv
```

Because the star is in the path, **two stars physically cannot write the same file.**
Sirius float-drift reruns (RYA-420) and what-if branches each get their own namespaced
artifact and diff *named* files — they never fight over one generic name.

- `data/outputs/` is **gitignored** — these are regenerable working products.
- The single accessor is **`pipeline/data_namespace.py`** (`output_path(star, name)`,
  `outputs_dir(star)`, `diagnostics_dir(star)`). Do not hand-build these paths elsewhere;
  the `{star}_` prefix is applied centrally so a bare generic filename cannot be emitted.

### 2. The Sun is FROZEN and VERSIONED

```
data/reference/solar/solar_abundances_v{N}.csv   # write-once, immutable, COMMITTED
data/reference/solar/CURRENT                      # names the active version (e.g. "v1")
data/reference/solar/hash_manifest.json           # sha256 per version — the guard
```

- Each version embeds a **provenance header** (`#`-comment lines: source commit, frozen
  UTC, Fe anchor, frozen verdict counts, changelog vs the previous version, `supersedes`).
  Readers go through the comment-aware accessor; the header is part of the hashed artifact.
- `data/reference/` is **committed** (NOT gitignored) — it is the gold reference.
- **`v1` freezes the final solar verdict** `PASS=4 / NLTE-OWED=1 / CURATION-OWED=21 /
  DATA-GAP=0` (RYA-371/462/467b, Kitt Peak + K-NLTE applied), Fe I = 7.516.

### 3. Re-baselining is a version BUMP, never an overwrite

A working solar run writes to the **namespaced working path** (`data/outputs/solar/…`),
**not** the reference. Promotion is a deliberate, reviewed act:

```
python scripts/promote_solar_reference.py --apply --changelog "..."
```

It validates the Fe anchor (~7.516 — a perturbed regen is **refused, not frozen**), writes
the **next** version `v{N+1}`, records its hash, and repoints `CURRENT`. It **refuses to
overwrite an existing `v{N}`.** Old versions are kept forever.

### 4. The differential denominator pins a version

When a target computes `[X/H]` vs the Sun it reads a **pinned** `solar_abundances_v{N}`
(default `CURRENT`) via `data_namespace.read_solar_reference()` / `differential_denominator()`
and **records `solar_ref_version` in its own output provenance** (`abundances_derive.run`
stamps a `solar_ref_version` column for non-solar stars; the Phase-C verdict records it in
its summary). So re-baselining the Sun later **never silently changes** a target's
already-derived numbers — you re-run the target against the new version on purpose.

### 5. Immutability is CI-enforced

`data_namespace.assert_frozen_references()` re-hashes every committed `solar_abundances_v{N}.csv`
against `hash_manifest.json` and **fails loud** on any mismatch — or on a version present on
disk but absent from the manifest. A frozen version cannot be edited (even a comment) without
tripping the guard; a new baseline must go through `promote_solar_reference.py`.

## Promotion workflow (summary)

1. Run a solar pass → `data/outputs/solar/solar_abundances.csv` (gitignored, regenerable).
2. `python scripts/promote_solar_reference.py` (dry-run) — review Fe anchor + next version.
3. `--apply --changelog "..."` — writes `v{N+1}`, hashes it, repoints `CURRENT`.
4. Commit `data/reference/solar/` in a reviewed PR. Targets opt into the new version by
   re-running against it; their old outputs keep their pinned `solar_ref_version`.

## Consequences

- The footgun is closed: the committed gold is immutable and hash-guarded; working runs are
  namespaced and gitignored, so no incidental regen can become the baseline.
- Multi-target era is unblocked: every star's products are collision-free by construction.
- **Scope note / follow-on:** the production writers (`abundances_derive.run`,
  `phase_c_verdict`, the EW-integrity + Fe-triage diagnostic emitters) and the
  `problem_children` reader are migrated. Historical one-off diagnostic scripts that read
  `data/processed/{star}_*` directly (regenerable working files) are migrated lazily; they
  are not on the production path and not in CI. New code MUST use `pipeline/data_namespace.py`.
```
