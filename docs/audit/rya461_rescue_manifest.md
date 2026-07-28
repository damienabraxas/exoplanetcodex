# RYA-461 — artifact rescue + worktree cleanup manifest

_Executed 2026-06-28 on Mr. Code's local machine. Rescue is non-destructive and ran
BEFORE any cleanup; cleanup removed only confirmed-merged, clean worktrees._

## Canonical store (created)

`~/Documents/Exoplanet Codex/` (the directory that contains the worktrees; override
with `$CODEX_ARTIFACT_STORE`), subfolders:

```
plots/   data/   diagnostics/   grids/   atlases/
ARTIFACT_MANIFEST.csv   (every preserved file: md5 + source worktree/branch + date)
```

Convention codified in `docs/CONVENTIONS.md`; one-call helper in
`pipeline/artifact_store.py` (`save_artifact(path, kind)`), tested in
`tests/test_artifact_store_rya461.py`; rescue tool in
`scripts/rescue_worktree_artifacts_rya461.py`.

## Part B — rescue (non-destructive)

Scanned all 37 worktrees for gitignored, at-risk artifacts (diagnostic plots,
normalized-spectrum intermediates, NLTE `.grd` grids). Tracked files (committed test
fixtures `oi6300_synthetic_test.png`, committed `curation_diagnostics_*.csv`, committed
reference atlases) were correctly excluded.

- **found (gitignored, at-risk): 40**
- **newly rescued: 29** — 26 diagnostic plots + 3 normalized intermediates
- **already in store (dedup by md5): 11** — duplicate `solar_normalized.csv` copies
- **invariant: 29 + 11 == 40** ✓

### Plots rescued (26 → `plots/`)
engine_comparison_{solar,procyon}_{delta,ep_rew,nlines,xh} (8); fe_{histogram,
excitation_potential,reduced_ew,wavelength} (4); {solar,procyon}_continuum_diagnostic;
{solar,procyon}_ew_diagnostic; solar_{abundances,ca6122_diagnostic,ew_diagnostic,
fe_excitation,oi6300_diagnostic,residuals}; rya238_fe_validation (×2 versions);
rya248_fe_scatter; solar_continuum_diagnostic (×3 distinct versions, md5-suffixed).

### Normalized intermediates rescued (3 → `data/`)
- `solar_normalized.csv`  md5 `df4a49cf5d69876628ebca5f2ef8cdfc`  (22 MB; was identical in 11 worktrees)
- `procyon_normalized.csv` md5 `8ca4e8ac…`  (20 MB)
- `procyon_normalized__842b4d06.csv` md5 `842b4d06…`  (20 MB; a distinct earlier normalization)

### NLTE `.grd` grids — NONE FOUND (missing → re-download candidates)
No `*.grd` exists in any worktree (confirms the prior teardown reclaimed ~26 GB).
The Amarsi/GALAH NLTE source grids must be **re-fetched from Zenodo rec 3982506** if
needed again — this is the explicit grid-audit hand-off. (The small committed `.csv`
correction grids in `data/nlte_grids/` are unaffected — they are in git.)

## Part C — cleanup (confirmed-merged + clean only)

Merge status per worktree was established two ways: `git merge-base --is-ancestor`
(direct merge) and `git cherry origin/main` zero unmerged-patch count (cherry-pick
merge). Only worktrees that were merged AND had no uncommitted tracked changes were
removed; branches were left intact on the remote.

- **Removed (20):** rya-342, 371, 384, 423, 431, 433, 440, 441, 442, 443, 446, 447,
  448, 449, 452, 453, 454, 456, 458, 459.
- **Disk reclaimed:** ~529 MB.
- **Kept (17):**
  - core: `exoplanetcodex` (dev tree), `exoplanetcodex-main`, `exoplanetcodex-rya461`
  - unmerged (have un-landed commits): rya-341, 377, 378, 380, 381, 382, 416, 419,
    420, 424, 427, 437, 444, **460**.

**Note on rya-460:** the ticket listed it as a cleanup candidate, but it has 1 unmerged
commit (RYA-459 landed via PR#80; its follow-up RYA-460 has not) → **kept**.

## Guarantees

- Nothing tracked-and-committed was touched.
- No artifact was deleted before being confirmed in the store (rescue ran first; the
  rescue invariant 29+11==40 was asserted before any `git worktree remove`).
- Every removed worktree's at-risk gitignored artifacts are in the store
  (`solar_normalized.csv`, rya-371 `solar_continuum_diagnostic.png` →
  `__5e7b90da`, rya-446 `rya238_fe_validation.png` → `__111a3354`).
