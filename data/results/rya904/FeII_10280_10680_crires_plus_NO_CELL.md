# RYA-904 — there is no CRIRES+ Fe II NIR cell, and here is why

The ticket asked for the Fe II cell **"if the 2 lines survive"**. They do not, and this
file exists because a missing product with no stated reason is indistinguishable from one
nobody tried (RYA-833).

## What was run

```
python3 scripts/derive_band_products.py --element Fe --ion II \
        --lo 10280 --hi 10680 --instrument crires_plus
```

It reached the corrected holding — `[holding] crires_plus ->
solar_crires_plus_y_rya794 (pre_normalised=True)` — and then refused, loudly:

```
no Fe 2 line in 10280.0-10680.0 A has theoretical depth in [0.15, 0.9] — check the column
```

Note that this command could not even have been asked before RYA-904 on two counts: the
holding was unaddressable, and `--ion` was decorative (`select_lines` was pinned to
Fe I, so `--ion II` would have returned **Fe I** lines labelled Fe II).

## Why

**The counts are from two different line lists, and the ticket quotes the other one.**

* The ticket's "48 Fe I + 2 Fe II in window" is from Elgueta+2026's `atomicy.dat`,
  which is their catalogue, not our production input.
* Our NIR synthesis list, `data/linelists/ispec_ir_9200_13000/atomic_lines.tsv`, holds
  **2 `Fe 2` rows in the whole 9200-13000 A extract** and exactly **one** inside this
  window:

  | wave_A | species | log gf | EP (eV) | theoretical_depth |
  | -- | -- | -- | -- | -- |
  | 9997.598 | Fe 2 | — | — | 0.087 |
  | **10501.503** | **Fe 2** | **-2.086** | **5.549** | **0.057** |

  (Our own line accounting agrees: one Fe II at 10501.503, predicted depth 0.057.)

* RYA-759's selection floor is `DEPTH_FLOOR = 0.15`, there to keep candidates visible
  above the band's crowding. At a theoretical central depth of **0.057** the single
  in-window Fe II line is **~2.6x below the floor** — it is not a marginal call.

## What this is NOT

It is **not** a telluric refusal, **not** a coverage gap, and **not** the dispatch defect
this ticket fixes — the loader served the window and the species selection was the thing
that had nothing to select. Lowering the floor to admit one 0.057-deep line in order to
report an n=1 Fe II abundance would be tuning a ratified selection rule to manufacture a
cell (RYA-161), and an n=1 product has no scatter to report either.

Fe II in this band is **owed coverage**, not a defect: it needs a line list that reaches
deeper, or a different window.
