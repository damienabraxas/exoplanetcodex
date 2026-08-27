# Worktree prune/archive plan — RYA-1087 (READ-ONLY)

**Reclaimable now (PRUNE): 20.3 GB** across 26 worktrees. 
A further **24.2 GB** in 30 worktrees needs archiving first. 
8 KEEP.


> Nothing here was deleted, pruned, or unprotected. The plan is the deliverable.


## Verdicts

| path | GB | branch status | prot | dirty | orphans | verdict | action |
|---|---|---|---|---|---|---|---|
| `exoplanetcodex` | 1.70 | merged | Y | 0 | 24 | **ARCHIVE** | save_artifact() 24 orphan file(s) first |
| `air` | 1.57 | detached |  | 3 | 37 | **ARCHIVE** | save_artifact() 37 orphan file(s) first |
| `rya1015` | 0.82 | UNMERGED(5/5 files differ) |  | 14 | 274 | **ARCHIVE** | git bundle/tag the branch before removal |
| `rya1070` | 0.81 | UNMERGED(11/11 files differ) |  | 0 | 0 | **ARCHIVE** | git bundle/tag the branch before removal |
| `rya1075` | 0.80 | UNMERGED(16/16 files differ) |  | 0 | 1 | **ARCHIVE** | git bundle/tag the branch before removal |
| `rya1069` | 0.80 | UNMERGED(8/8 files differ) |  | 1 | 0 | **ARCHIVE** | git bundle/tag the branch before removal |
| `rya1079` | 0.80 | UNMERGED(8/8 files differ) |  | 0 | 0 | **ARCHIVE** | git bundle/tag the branch before removal |
| `base1075` | 0.80 | detached |  | 1 | 0 | **KEEP** | uncommitted tracked changes present |
| `rya1072` | 0.80 | merged |  | 0 | 0 | **PRUNE** | merged/landed, no orphan artifacts, clean |
| `rya1078` | 0.80 | merged |  | 0 | 0 | **PRUNE** | merged/landed, no orphan artifacts, clean |
| `base1070` | 0.80 | detached |  | 1 | 0 | **KEEP** | uncommitted tracked changes present |
| `ci1059` | 0.80 | UNMERGED(5/5 files differ) |  | 0 | 0 | **ARCHIVE** | git bundle/tag the branch before removal |
| `rya821` | 0.79 | merged |  | 0 | 0 | **PRUNE** | merged/landed, no orphan artifacts, clean |
| `rya1081` | 0.79 | UNMERGED(37/37 files differ) |  | 0 | 0 | **ARCHIVE** | git bundle/tag the branch before removal |
| `main853` | 0.79 | detached |  | 0 | 0 | **PRUNE** | merged/landed, no orphan artifacts, clean |
| `rya1080` | 0.79 | merged |  | 0 | 0 | **PRUNE** | merged/landed, no orphan artifacts, clean |
| `rya1084` | 0.79 | UNMERGED(19/19 files differ) |  | 0 | 0 | **ARCHIVE** | git bundle/tag the branch before removal |
| `rya1071` | 0.79 | merged |  | 0 | 0 | **PRUNE** | merged/landed, no orphan artifacts, clean |
| `rya1077` | 0.79 | merged |  | 0 | 0 | **PRUNE** | merged/landed, no orphan artifacts, clean |
| `rya1083` | 0.79 | merged |  | 0 | 0 | **PRUNE** | merged/landed, no orphan artifacts, clean |
| `codex-tracker` | 0.79 | UNMERGED(7/15 files differ) |  | 0 | 0 | **ARCHIVE** | git bundle/tag the branch before removal |
| `rya853` | 0.79 | UNMERGED(13/13 files differ) |  | 0 | 0 | **ARCHIVE** | git bundle/tag the branch before removal |
| `maincheck` | 0.79 | detached |  | 0 | 0 | **PRUNE** | merged/landed, no orphan artifacts, clean |
| `rya1058-finish` | 0.78 | merged |  | 0 | 0 | **PRUNE** | merged/landed, no orphan artifacts, clean |
| `redopt` | 0.78 | merged |  | 2 | 4 | **ARCHIVE** | save_artifact() 4 orphan file(s) first |
| `rya1044r` | 0.78 | merged |  | 0 | 0 | **PRUNE** | merged/landed, no orphan artifacts, clean |
| `rya1050` | 0.78 | merged |  | 0 | 0 | **PRUNE** | merged/landed, no orphan artifacts, clean |
| `rya1051b` | 0.78 | merged |  | 0 | 0 | **PRUNE** | merged/landed, no orphan artifacts, clean |
| `rya1059` | 0.78 | merged |  | 0 | 0 | **PRUNE** | merged/landed, no orphan artifacts, clean |
| `rya1060-ci` | 0.78 | detached |  | 0 | 0 | **PRUNE** | merged/landed, no orphan artifacts, clean |
| `rya1063` | 0.78 | merged |  | 0 | 0 | **PRUNE** | merged/landed, no orphan artifacts, clean |
| `pr419` | 0.78 | detached |  | 0 | 0 | **PRUNE** | merged/landed, no orphan artifacts, clean |
| `rya1065` | 0.78 | merged |  | 0 | 0 | **PRUNE** | merged/landed, no orphan artifacts, clean |
| `rya1066` | 0.78 | merged |  | 0 | 0 | **PRUNE** | merged/landed, no orphan artifacts, clean |
| `rya987` | 0.78 | merged |  | 1 | 26 | **ARCHIVE** | save_artifact() 26 orphan file(s) first |
| `base1035` | 0.78 | detached |  | 0 | 0 | **PRUNE** | merged/landed, no orphan artifacts, clean |
| `rya1020` | 0.78 | merged |  | 4 | 0 | **KEEP** | uncommitted tracked changes present |
| `rya1040` | 0.78 | UNMERGED(9/11 files differ) |  | 0 | 0 | **ARCHIVE** | git bundle/tag the branch before removal |
| `mainchk` | 0.78 | detached |  | 0 | 0 | **PRUNE** | merged/landed, no orphan artifacts, clean |
| `tracker_merge_test` | 0.78 | detached |  | 211 | 0 | **KEEP** | uncommitted tracked changes present |
| `mainchk_wt` | 0.78 | detached |  | 0 | 0 | **PRUNE** | merged/landed, no orphan artifacts, clean |
| `rya992` | 0.78 | merged |  | 0 | 2 | **ARCHIVE** | save_artifact() 2 orphan file(s) first |
| `rya1006-b2` | 0.78 | detached |  | 0 | 0 | **PRUNE** | merged/landed, no orphan artifacts, clean |
| `rya1006-baseline` | 0.78 | detached |  | 0 | 0 | **PRUNE** | merged/landed, no orphan artifacts, clean |
| `rya1006` | 0.78 | merged |  | 0 | 0 | **PRUNE** | merged/landed, no orphan artifacts, clean |
| `rya1008` | 0.78 | UNMERGED(1/1 files differ) |  | 0 | 0 | **ARCHIVE** | git bundle/tag the branch before removal |
| `exoplanetcodex-rya1000` | 0.78 | merged | Y | 0 | 0 | **KEEP** | inside the artifact store root |
| `exoplanetcodex-rya935t` | 0.78 | UNMERGED(7/7 files differ) | Y | 0 | 0 | **ARCHIVE** | git bundle/tag the branch before removal |
| `exoplanetcodex-rya989` | 0.78 | merged | Y | 0 | 0 | **KEEP** | inside the artifact store root |
| `rya-951-graded-fe` | 0.77 | UNMERGED(1/1 files differ) |  | 1 | 0 | **ARCHIVE** | git bundle/tag the branch before removal |
| `exoplanetcodex-rya759` | 0.73 | UNMERGED(128/134 files differ) | Y | 0 | 0 | **ARCHIVE** | git bundle/tag the branch before removal |
| `exoplanetcodex-rya770` | 0.73 | UNMERGED(128/135 files differ) | Y | 0 | 0 | **ARCHIVE** | git bundle/tag the branch before removal |
| `rya-851-science-input` | 0.72 | detached |  | 0 | 0 | **PRUNE** | merged/landed, no orphan artifacts, clean |
| `exoplanetcodex-rya800` | 0.71 | UNMERGED(7/7 files differ) | Y | 0 | 0 | **ARCHIVE** | git bundle/tag the branch before removal |
| `exoplanetcodex-rya762` | 0.70 | UNMERGED(1/5 files differ) | Y | 0 | 0 | **ARCHIVE** | git bundle/tag the branch before removal |
| `exoplanetcodex-rya791` | 0.70 | UNMERGED(14/14 files differ) | Y | 0 | 0 | **ARCHIVE** | git bundle/tag the branch before removal |
| `exoplanetcodex-rya793` | 0.70 | UNMERGED(1/1 files differ) | Y | 0 | 0 | **ARCHIVE** | git bundle/tag the branch before removal |
| `exoplanetcodex-rya279` | 0.70 | UNMERGED(1/2 files differ) | Y | 0 | 8 | **ARCHIVE** | git bundle/tag the branch before removal |
| `exoplanetcodex-rya772` | 0.70 | detached | Y | 0 | 0 | **KEEP** | inside the artifact store root |
| `exoplanetcodex-rya758` | 0.70 | UNMERGED(8/8 files differ) | Y | 0 | 0 | **ARCHIVE** | git bundle/tag the branch before removal |
| `exoplanetcodex-rya546` | 0.70 | merged | Y | 0 | 0 | **KEEP** | inside the artifact store root |
| `exoplanetcodex-rya558` | 0.70 | UNMERGED(6/9 files differ) | Y | 3 | 0 | **ARCHIVE** | git bundle/tag the branch before removal |
| `exoplanetcodex-rya420` | 0.62 | UNMERGED(1/1 files differ) | Y | 0 | 0 | **ARCHIVE** | git bundle/tag the branch before removal |
| `exoplanetcodex-rya341` | 0.50 | UNMERGED(7/7 files differ) | Y | 0 | 6 | **ARCHIVE** | git bundle/tag the branch before removal |
