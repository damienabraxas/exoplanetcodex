# Exoplanet Codex — engineering conventions

## Per-star output namespacing + frozen gold solar reference (RYA-469)

Every per-star pipeline product carries the star in its **path**
(`data/outputs/{star}/{star}_*`), so two stars cannot collide on a filename; the
gold-standard solar differential denominator is **frozen + versioned + immutable**
(`data/reference/solar/solar_abundances_v{N}.csv`, `CURRENT` pointer, hash-guarded).
Use `pipeline/data_namespace.py` for all of it; re-baseline the Sun only via
`scripts/promote_solar_reference.py` (bump, never overwrite). Full rule:
[`docs/design/adr_data_namespacing_and_gold_reference.md`](design/adr_data_namespacing_and_gold_reference.md).

## Artifact preservation: save-before-clean (RYA-461)

**The recurring failure this prevents:** gitignored artifacts produced inside a git
worktree — diagnostic plots, downloaded atlases, NLTE `.grd` grids, normalized-spectrum
intermediates — live **only** in that worktree. When the worktree is cleaned or removed,
they are **lost**. For the large NLTE grids that turns "we have the grid" into
"re-download 26 GB."

### The standing rule

> **Any gitignored artifact a pipeline or diagnostic run produces must be copied to the
> canonical local store as part of the run — never left only in a worktree.**

Concretely, right after a run writes a gitignored file it intends to keep (a plot, a
diagnostic CSV/JSON, a downloaded atlas, an NLTE grid, a normalized spectrum), call the
one-line helper:

```python
from pipeline.artifact_store import save_artifact
save_artifact("results/plots/solar_oi6300_diagnostic.png", kind="plots")
save_artifact("data/processed/solar_normalized.csv",        kind="data")
```

`save_artifact(path, kind)` copies the file into the store, de-duplicates by md5, and
records its provenance (source worktree/branch + date) in `ARTIFACT_MANIFEST.csv`. It
never deletes or mutates the source.

### The canonical store

Lives **outside** any worktree, in the directory that contains the worktrees
(default `~/Documents/Exoplanet Codex/`, override with `$CODEX_ARTIFACT_STORE`):

| subfolder      | holds                                                              |
|----------------|--------------------------------------------------------------------|
| `plots/`       | diagnostic plots (`*.png` / `*.pdf` / `*.svg`)                      |
| `data/`        | normalized-spectrum and other reusable data intermediates          |
| `diagnostics/` | diagnostic CSV/JSON outputs (audit / proof / verdict tables)       |
| `grids/`       | NLTE `.grd` / model grids (large, external, expensive to re-fetch)  |
| `atlases/`     | downloaded reference atlases (Kitt Peak, CALSPEC, IRTF, …)          |

`ARTIFACT_MANIFEST.csv` at the store root records every preserved file with its md5 and
provenance.

### Cleaning worktrees (the destructive half)

1. **Rescue first.** Before removing any worktree, run the save-before-clean rescue so
   every at-risk gitignored artifact is in the store. Rescue is **non-destructive** —
   it copies, never deletes.
2. **Only remove confirmed-merged worktrees.** Verify the branch's work is on `main`
   (a merge, or a patch-equivalent cherry-pick) before `git worktree remove`. When in
   doubt, **keep** the tree.
3. **Never delete an artifact that is not yet confirmed in the store.**

This rule is also the reason the `codex-artifact-preservation` step belongs in every
Mr. Code brief that produces gitignored output.
