# Model registry — the bare `ENGINE-B` collision

Companion to [`data/catalog/model_registry.csv`](../../data/catalog/model_registry.csv).

## The hazard

`ENGINE-B` is a **live token** in `band_products.TREATMENTS`, and it does not mean what it
used to. Read from `pipeline/treatment_axes.py` (`LEGACY`), not from memory:

| token | scale | model | resolves to |
|---|---|---|---|
| `ENGINE-B` | `1D-LTE` | `none` | **model 1** — 1D-LTE |
| `ENGINE-B-NLTE` | `1D-NLTE` | `gerber` | **model 4** — 1D-NLTE Gerber |

Before RYA-906's physics-axis renaming, the bare string meant **the Gerber NLTE synthesis**,
i.e. model 4. After it, the bare string resolves to **model 1**.

**Same string, two physics — and an artifact written before RYA-906 does not record which
was meant.** That is the whole hazard: it is not a naming inconvenience, it is a silent
mis-attribution of physics across a rename boundary.

## How it is handled

`pipeline.model_registry.resolve("ENGINE-B")` **raises**, naming both candidates and the
canonical token for each. It does not pick one.

The guard lives in the registry resolver, **not** in `band_products.py`. RYA-1101 forbids
touching model code, and `ENGINE-B` is a legitimately live treatment with real published
products — a guard that raised wherever the string appears would break production runs. What
must not happen is a bare token being bound to a *model identity*, and that binding only ever
happens here.

## ⚠️ Not decided here

**The canonical binding for legacy `ENGINE-B` data is Ryan's call**, posted for sign-off on
RYA-1101 rather than settled in code. The options are:

1. **Bind to model 1** (what the code says today) — consistent with the current resolver,
   but silently re-labels pre-RYA-906 artifacts that meant model 4.
2. **Bind to model 4** for artifacts older than RYA-906, model 1 after — historically
   correct, needs a date or commit boundary that is itself a judgement.
3. **Refuse permanently** (today's behaviour) — every consumer must state which it means.

Until that is decided, option 3 stands, because it is the only one that cannot silently
publish the wrong physics.

## Two live code tokens that are NOT models in this roster

Found by reconciling the registry against `band_products.TREATMENTS`:

- **`ENGINE-B`** — the alias above.
- **`1D-LTE-LABGF`** (RYA-836) — same engine and route as model 1, differing in one input:
  the oscillator strength comes from a primary laboratory measurement instead of Kurucz.
  `band_products.py:60` documents it as *"a separate member rather than a variant"*. It is a
  real, live treatment and it is **not** one of the eight models in this roster. Whether it
  should be a ninth row is a scope question for RYA-1101, flagged rather than decided.

## Model 9 — "Frankenstein's Dog", and the premise that did not survive the code

Ryan ratified the Dog as a **distinct model** on RYA-1101 (2026-08-28) and RYA-1115 carries
that forward. It is in the roster as model 9. Its identity, in his words: *the 1D-LTE Gerber
leg computed on the mean-3D method's OWN reference atmosphere* — not the MARCS.GES
`synth-1D-LTE-gerber` (model 2).

**It has no `stored_token`, no `atmosphere` and no emitter, because none exists in the code.**
That is a measured result, not an omission:

- `treatment_axes.AXIS_NATIVE` registers exactly three tokens —
  `synth-mean3D-NLTE-gerber-stagger`, `synth-mean3D-LTE-gerber-stagger` and
  `synth-1D-LTE-gerber` (`atmos=marcs-ges`). None of them is a 1D-LTE leg on a mean-3D
  reference.
- `ATMOSPHERES` is a **closed** vocabulary — `("atlas9", "marcs-ges", "stagger-mean3d")`.
  There is no 1D reference atmosphere belonging to the STAGGER mean method for such a leg
  to run on.
- `derive_band_products.py` has no deck that would emit it: `--engine-b-deck` accepts
  `ts-lte | gerber-nlte | gerber-1d-lte | gerber-mean3d | gerber-mean3d-lte`, and every one
  of those maps to a token already in this roster.

So the row is written with the blanks left blank, under the same rule model 8 follows: **a
token that is not in the code is not written down here.**

### ⚠️ The stated proof is refuted by the artifact it cites

The RYA-1101 comment argues the Dog *"cannot exist without this 1D-LTE-on-mean-reference
leg"*, because RYA-1083's ladder computes the 1D→mean-3D atmosphere term at +0.092.

The ladder computes that term **without any such leg**. From
`scripts/rya1051_fe_nlte_ladder.py`:

```
"atmosphere_1D_to_mean3D": ("mean3d_lte", "lte_1d_atlas9", ...)
```

and from the committed `data/results/rya1051/fe_nlte_ladder.json`:

| term | minuend − subtrahend | median | n |
|---|---|---|---|
| `atmosphere_1D_to_mean3D` | `mean3d_lte` − `lte_1d_atlas9` | **+0.092** | 67 |
| `nlte_on_mean3D` | `mean3d_nlte` − `mean3d_lte` | +0.032 | 67 |
| `nlte_on_1D` | `nlte_1d` − `lte_1d_marcs` | +0.048 | 67 |

The +0.092 is **model 5 minus model 1** — the ⟨3D⟩-LTE product against the ATLAS9 1D-LTE
product. The ladder takes five per-line inputs and all five are already models 1, 2, 4, 5
and 6. The Dog is not among them.

That does not make the Dog wrong — it makes the *reason given for it* wrong, and the two
must not be confused. The open question, **for Ryan, not for this file**, is which of these
the Dog actually is:

1. **Model 5 under another name.** `synth-mean3D-LTE-gerber-stagger` already *is* "the
   Gerber deck's LTE limit on the mean-3D method's own atmosphere". If that is what the Dog
   means, model 9 is a duplicate of model 5 and should be withdrawn.
2. **A genuinely new leg** — 1D-LTE Gerber on the STAGGER method's own *1D reference*
   atmosphere, which would make today's +0.092 a slightly confounded number (it currently
   spans ATLAS9 → STAGGER-mean, mixing the ATLAS9→reference step into the atmosphere term).
   If so, model 9 is real, and building it needs a new entry in `ATMOSPHERES` and a new deck
   — neither of which this ticket may add (RYA-1101 forbids touching model code).

Reading 2 is the one that makes the Dog "missing", and it is the reading the row is written
under. **`status=not-emitted` records exactly that: the science calls for it, and nothing
produces it.** It is kept out of `pair_group=frankenstein` because it would be the
atmosphere-rung comparand, not a member of the mandatory 5/6 NLTE pair.

## The `line_set` hook (RYA-1115, for RYA-1111)

`line_set` is an open column with a **closed vocabulary**, declared in
`pipeline.model_registry.LINE_SETS`:

```
"-" | "asplund-graded" | "gbs" | "our-graded" | "our-deep-graded"
```

Which pool a measurement used is a property of the **product**, not of the model — any model
here can be measured on any of these sets — so every roster row carries `-` today. RYA-1111
populates it when it has products to key.

The column is opened now, with its guard already in place, so that the vocabulary exists
*before* the first value is written and cannot be typed from memory that day.
`verify_model_registry.py` check (f) refuses any value outside the tuple.

⚠️ **`consistent` is deliberately absent.** RYA-1105 removes the Consistent tier from the
active pipeline and the website; the going-forward set is Asplund Grade / Our Grade / Deep
Grade. But `--lines-tier consistent` is **still live** in `derive_band_products.py`. The code
and this vocabulary genuinely disagree today. Retiring the flag is RYA-1105's job, not this
ticket's, so the disagreement is recorded here rather than papered over.

## Correction: bare `ENGINE-B` is not deprecated

An earlier draft of `model_registry.csv` called the bare token a *"DEPRECATED alias"*. That is
false, and the row now says so. `--engine-b-deck` **defaults to `ts-lte`**, and that default
path sets `eb_treatment = "ENGINE-B"` with its own RYA-1044 provenance string. The bare token
is what the default deck emits today.

This makes the collision worse, not better: the ambiguous string is not a fossil to be read
out of old artifacts, it is being written by the default code path right now. Option 3
(refuse to bind) remains the only safe resolver behaviour until Ryan rules.
