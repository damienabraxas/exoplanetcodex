# RYA-783 matrix — stated capability limits on published cells

**RYA-1055. Annotation, not correction. No value in `fe_product_matrix.csv` is changed by
this file, and none should be.**

## The two annotated cells

| Fe | band | treatment | A(Fe) | n | n_excl | annotation |
| -- | -- | -- | -- | -- | -- | -- |
| II | VIS | `ENGINE-B-NLTE` | **7.470** | 8 | 3 | **cannot be shown to be NLTE — LTE-equivalent** |
| II | red-optical | `ENGINE-B-NLTE` | **7.461** | 2 | 2 | **cannot be shown to be NLTE — LTE-equivalent** |

Both are `kpno_solar_atlas`, and both come from
`gerber-nlte/FeII_*_kpno_solar_atlas_PROFILEFIT_products.csv` in the RYA-845 run.

## Why

`atom.fe607a` — the model atom BOTH registered Gerber decks load (`Fe` and `Fe@mean3D`,
`pipeline.gerber_nlte.DECKS`) — declares:

```
levels        548 Fe I   (idx 1-548,   0.000 -  7.505 eV)
               58 Fe II  (idx 549-606, 7.950 - 23.753 eV)
                1 Fe III (idx 607,     24.102 eV)
bound-bound  12,635 declared, 12,635 parsed, 12,635 Fe I -> Fe I
                  0 with BOTH levels Fe II
                  0 with EVEN ONE level Fe II
highest level index appearing in any radiative transition : 546   (< 549)
```

The 58 Fe II levels carry **no bound-bound transitions at all**. They are a pure
ionisation reservoir — the targets of Fe I photoionisation, which in turn ionise to
Fe III. bsyn applies departures **per line** and falls back to **departure = 1** for any
line whose levels are unidentified, so every Fe II line synthesises in LTE regardless of
which deck is loaded, and **no line list — VALD or otherwise — can change that.**

Measured by `scripts/rya1055_atom_ion_reach.py` on the staged bytes
(md5 `d08dc8232ed68eec65f9bb6631e82ea8`, the file `Fe_gerber2023.prov.json` pins);
artifact `data/results/rya1055/atom_ion_reach.json`.

## What the annotation does and does not claim

* It says **"cannot be shown to be NLTE"**, not "is wrong". Whether these cells were ever
  NLTE is **unrecorded and could not be determined**: the artifacts predate
  `provenance.txt`, predate the RYA-906 route/scale/model/atmos columns, and predate
  RYA-880's `nlte_delta_dex` / `nlte_source`.
* ⚠️ **The NLTE-vs-LTE per-line difference is NOT evidence of departures.** RYA-1055
  measured it at a median −0.0075 dex, nonzero on all 8 shared lines — and that is fully
  explained without any departure at all: `ENGINE-B-NLTE` ran on **MARCS.GES** while
  `ts-lte` `ENGINE-B` ran on the route's **own atmosphere**, and RYA-1045 measured that
  atmosphere step at **+0.004** for Fe I on the same instrument. The difference is
  evidence of a different **atmosphere**, not of applied departures.
* **The Fe I `ENGINE-B-NLTE` cells are NOT annotated and are not affected.** 12,635
  transitions is the whole Fe I term system; the limit is stage-specific.
* **`ENGINE-A` cells are NOT annotated.** That leg reads the MPIA/Bergemann per-line delta
  grid (`data/nlte_grids/Fe_Bergemann_MPIA.csv`, 6,400 Fe II rows over 80 lines,
  3805.5–6586.7 Å), never this atom. Its Fe II corrections are real and small: measured at
  the solar node they run **−0.002 to +0.016 dex, median +0.000**, with 146 of 160 entries
  inside ±0.005 — the tail is the three strong multiplet-42 lines
  4923.932 / 4924.921 / 5169.033 at +0.015 (Fe I control at the same node: median +0.011, up to +0.040). The limit is
  a property of the **deck**, not of Fe II.

## Annotate, do not delete

Ryan, 2026-08-30 and re-ratified 2026-09-03: these are the only Fe II Engine-B numbers we
hold, and *"a cell marked 'cannot be shown to be NLTE' is more useful than a missing one.
Deleting also destroys the only record that the question was ever asked."* The RYA-1050
pool guard is what prevents **new** ones: it refuses to emit a product when not one pooled
line carries an NLTE label, which is exactly this case (RYA-877's Fe II pool: 0 of 11).

## Where this annotation lives going forward

`scripts/rya783_fe_matrix_report.py._annotate_capability` emits a `capability_note` column
and prints a CAPABILITY LIMIT block, resolving the verdict through
`pipeline.gerber_nlte.nlte_ion_capability` — the same accessor the synthesis guard and the
product stamp read, so the report cannot state a reach the pipeline disagrees with. Run
today against the committed matrix it flags **exactly these 2 of 18 rows**, and leaves both
Fe I `ENGINE-B-NLTE` cells alone.

⚠️ **The committed `fe_product_matrix.csv` and `.png` are NOT hand-edited.** They predate
the generator in more ways than this one (absolute `_src` paths from before RYA-845, no
RYA-906 axis columns) and they cannot be regenerated from this checkout — their inputs are
the RYA-845 run under `/mnt/codex-data`. This is the same disposition RYA-869 took for the
four cells that still carry its pre-fix systematic bar: **flagged, not hand-edited.** The
plot likewise shows the two cells unannotated.

## The capability itself

Fe II NLTE is a **stated capability limit**, not a gap being papered over. A two-stage Fe
atom carrying Fe II bound-bound transitions is a fetch-or-build question of the same shape
as RYA-1035's ⟨3D⟩ deck hunt, and Ryan **deferred** it (2026-09-03) to the off-solar
programme, where the Fe I/Fe II balance does the log g work. At solar the missing effect is
literature-negligible: Amarsi's grid gives −0.0010 (vturb 0.75) and −0.0008 (1.50), and
Lind et al. 2017 §2.2.4 states outright *"For the Sun, NLTE effects on Fe II lines are
insignificant."*

**Fe II ⟨3D⟩-LTE is a different axis and IS available** — model 5
(`synth-mean3D-LTE-gerber-stagger`): the ⟨3D⟩-mean atmosphere is ion-agnostic. Never label
it 3D-NLTE. And model 5 == model 6 on Fe II **by construction**, which is itself the
phantom-departure check (RYA-1135).
