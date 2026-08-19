# Decision record — what an n = 1 product publishes for its statistical uncertainty

**Ticket:** RYA-907. **Date:** 2026-08-19.
**Branch:** ryandamienschmitt/rya-907-bug-n1-products-publish-stat_dex-00-productsigma-or-00.
**Decision: OPTION (a) — the RYA-771 quantiser floor, with the line-to-line scatter term
carried as UNMEASURED and a `stat_basis` field naming which rule produced the number.**

---

## The question

`pipeline/band_products.build_product` refuses to invent a scatter from a single line:

```python
sigma = float(np.std(vals, ddof=1)) if len(vals) > 1 else None
```

`None` is the honest answer — the spread of one number is undefined, not small. So what
does the *published* bar say? RYA-907 §4.3 put three options:

**(a)** the quantiser floor, with the term marked unmeasured;
**(b)** a stated conservative proxy (e.g. the band's own multi-line scatter), clearly
labelled a proxy;
**(c)** refuse to publish a statistical figure at all and print `--`.

## The decision, and why

**(a).** The published statistical bar for an n = 1 product is `QUANTISER_FLOOR_DEX`
= **0.01 dex**, and `stat_basis` says in words that the floor is standing in for a
scatter nobody measured.

- **Against (c).** A blank is not free. It propagates as a NaN into
  `graded_reporting.total_sigma` and as an empty cell on the page, and the project has
  already been bitten by a row that arrives carrying no number: *"a row that arrives but
  carries no number is worse than one that does not arrive, because it looks like a
  measurement that failed rather than a wiring mismatch"* (`derive_band_products.py`,
  RYA-832). The measurement did not fail. One line was measured, and its abundance is
  real; what is missing is a *spread*, which is a different quantity.
- **Against (b).** A proxy from the band's multi-line pool is a number about a
  **different line set**. RYA-875 is the standing lesson: a residual computed as an
  18-line median against a scalar from a *different* 23-line set was an artifact of the
  mismatch, and vanished to 0.0000 once paired. Importing another pool's scatter would
  manufacture exactly that shape, and labelling it would not make it a measurement of
  this product.
- **For (a).** The floor is not an estimate of this product's scatter — it is a
  statement about the **instrument's resolution**. RYA-771 established that iSpec writes
  trial abundances with `%.2f` in both `bsyn` and `babsma`, so EW(A) is a staircase with
  0.01 dex treads. Two abundances closer together than one tread produce byte-identical
  synthetic spectra. So *no* product, at *any* line count, can honestly claim a bar finer
  than 0.01 dex; the n = 1 case is simply where nothing larger is known. It is the one
  number available that is a fact about the measurement chain rather than a guess about
  the line.

## What makes it non-silent

The floor alone would be a magic number, and a magic number is how the next reader
mistakes it for a measurement. Three things carry the claim instead:

1. **`stat_basis` is a published field**, written into every product artifact. Its three
   values are distinguishable: `measured`, `quantiser-floor … measured random RMS below
   one tread`, and `quantiser-floor … UNMEASURED`. The page can never render a bar whose
   origin is unstated.
2. **`ErrorBudget.describe()` lists the unmeasured term explicitly**, as
   `UNMEASURED  line-to-line scatter`, rather than dropping it. A term omitted from the
   printout reads as a term the budget does not have; this budget has it and could not
   measure it, which is the more useful thing to know.
3. **`Term.contribution()` raises** on an unmeasured term. There is no code path where it
   returns a number, so a future caller cannot re-acquire the defect by writing
   `or 0.0` — the value it would be reaching for does not exist.

## Scope

Applies to every element and band, not just Fe II red-optical. That cell is simply where
n = 1 landed first: across 39 Fe product rows on Sirius, `stat_dex == 0` occurred in
exactly two distinct cells, **both n = 1, and never once at n > 1**.

The n = 1 rows stay in the matrix. They are real measurements, and RYA-711 is explicit:
quarantined, never culled. The bar has to be honest, not the row deleted.

## Where the constant lives

`pipeline.error_budget.QUANTISER_FLOOR_DEX`, declared once and imported. Not re-typed at
any call site — RYA-845 is the standing lesson that a number written down twice drifts,
and the near-UV 0.100 double-count is what it cost.
