# RYA-847 item 6 — the gated cells, and what the gate cost

`gated/` holds the four cells the NON-MINIMUM check touches, regenerated with it applied.
The other five synthesis cells are deliberately absent: regenerating a cell nothing
changed in churns its artifact and buries the real diff (RYA-845).

`rya847_gate_diff.csv` is the item-6 diff. `rya847_gate_caught_lines.csv` names every line
the gate excluded, with the metric that caught it.

## 🔴 The baseline is the SWEEP, not the published matrix

`rya847_pregate_control.csv` carries the pre-gate numbers, and they come from RYA-847's own
nine-cell sweep — the **same tree and the same code** as the gated run, with the gate not
yet applied. That is the control the diff uses.

**Do NOT diff these cells against `data/results/rya783/fe_product_matrix.csv`.** That
matrix was produced by the rya845 run and its line counts differ for reasons that have
nothing to do with this gate — four ATOMIC_BLEND registry exclusions, 11119.795, the
RYA-807 registry gate. Comparing to it measures this gate PLUS months of pool drift, which
is the mistake RYA-848 made and proved by control.

⚠️ The sweep's own per-line CSVs no longer exist: the working tree was re-synced after the
sweep and took the untracked output directory with it (the RYA-461 save-before-clean
convention exists for exactly this and was not followed). The control numbers are
therefore carried here as DATA with their source named, not recomputed. Every one of them
is checked against the gated artifact by the arithmetic in `scripts/rya847_gate_diff.py`:
`n_post + n_caught == n_pre` must hold per cell, and `scatter == stat * sqrt(n)` must hold
on both sides. A control row that does not belong to this run fails those checks.
