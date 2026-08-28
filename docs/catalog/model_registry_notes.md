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
