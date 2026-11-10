# Atomic Data Gap Register (RYA-1378)

This is the list of lines kept for future laboratory or theoretical study. It contains every line of an element that the grading keeps out of the reference (LAB) pool.

Ryan 2026-11-09: "all lines that are held or whatever due to grading should be added to a list for future study".

## How the files are made

The files are written by `pipeline/gap_register.py`, in two ways:

- `run_pipeline.py --star solar --element El --gap-register` writes one element's register.
- Every element run also writes the register at the end.

Each element gets two files:

- `solar_<El>.csv` holds one row per line.
- `solar_<El>.json` holds the summary counts.

## What the register contains

**Values come from `canonical_gf.csv`, never from here.** gf, its source, its tier and its uncertainty are not copied into the register; join on `line_id` to get them.

**`held_reason` says why a line is held:**

| Value | Meaning |
|---|---|
| `REGISTRY_EXCLUDED` | The line registry excludes the line. |
| `UNGRADED_GF` | The gf has no published uncertainty, so the line is in no graded pool. |
| `GRADED_NOT_LAB` | The line is graded but has no laboratory gf, so it is in the codex pool only. |
| `NOT_IN_CANONICAL` | A curated laboratory line that has no row in `canonical_gf.csv`. |

**`flags` describe the atomic data. One line can carry several:**

- `UNRATED_GF`
- `MISSING_GF`
- `PROVENANCE_GAP`
- `GF_DISCREPANCY`
- `SOLAR_DISCREPANCY`

Curated discrepancies are kept in `data/reference/atomic_gap_known.csv`.

**`priority` ranks lines for study:**

| Priority | Meaning |
|---|---|
| 1 | No graded line of that species exists anywhere in the band, so every cell in that band is held. |
| 2 | The line has a discrepancy flag. |
| 3 | The line is in the synthesis list. |
| 4 | The line is not in the synthesis list. |

**`codex_fit_evidence` holds the orchestrator's problem-line codes.** These describe how the line behaved when we fitted it, not its atomic data. They count as fit evidence only when the fits ran on or after commit 01f51846, the RYA-1350 continuum gate. Regenerate the register after a rerun.

## What the register does not say

- A missing grade is not evidence that the gf is poor.
- Being listed does not promote, exclude or change any line, gf, grade or product.
- The register does not establish whether a line can be observed (its depth or blending).
