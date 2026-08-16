# RYA-825 — the `log_gf` column reports the gf the inversion actually used

**The defect:** `data/audit/line_accounting/per_line.csv` (and every pool built from it)
reported a `log_gf` copied from the linelist at intake, while the production synthesis
path resolved a **different** value out of `canonical_gf.csv` before inverting.
`_load_synth_resources` defaults `apply_canonical_gf=True`, which runs
`gf_resolver.apply_to_synth_array` (RYA-353). The column described an input nobody used.

Reproduce:

```
python3 scripts/line_accounting_rya709.py                      # regenerates the table
python3 scripts/rya825_reverify_scale_mismatch.py --pool <FeI_..._PASS.csv>
python3 -m pytest tests/test_gf_column_honesty_rya825.py
```

---

## Scale of the drift — and it was never Fe-only

| | |
|---|---|
| accounting rows | 21,828 |
| resolve against `canonical_gf` | 6,248 |
| **reported a gf the inversion did not use** | **3,197** |
| median \|Δ\| | 0.145 dex (max 9.643) |

Affected elements: **Fe 1773, Cr 326, Co 261, V 227, Ni 224, Mn 105, Ti 80, Sc 42, Mg 34,
Y 29, Na 28, Si 21** and more. Fe alone was 1,773 of its 7,541 rows. The ticket's
"generalizes" note is confirmed: the same `apply_canonical_gf` path serves every element,
so the stale column was pool-wide.

The remaining 15,580 rows sit **outside** canonical's 3780–9199.9 Å range, where the
resolver does not run and the intake value *is* the value used. Those are correct and
unchanged — flagged `gf_canonical = False` so the distinction is legible rather than
implied.

---

## What changed, and what deliberately did not

The table now carries four columns instead of one:

| column | meaning |
|---|---|
| `log_gf` | **the gf the inversion used** — canonical where the resolver reaches, intake otherwise |
| `log_gf_intake` | what the linelist said at intake. Reproduces the old column **exactly** |
| `gf_source_intake` | where that value came from |
| `gf_canonical` | whether the resolver supplied it |

That fourth column is what makes the other three readable: `log_gf == log_gf_intake` means
"they agree" only if you also know whether the resolver ran at all.

**A semantic change worth naming.** Where `gf_canonical` is true the value is now the
physical-line **total** that canonical stores. The accounting generator had aggregated HFS
components with `max`. Total is what the synthesis integrates, so it is the right meaning
for a column claiming to report the gf used — and the old convention survives intact in
`log_gf_intake`.

### No abundance moved, and that is checked rather than asserted

The accounting `log_gf` is **written and carried, never computed with**: a grep across
`scripts/` and `pipeline/` finds only the writer. The one consumer that *does* read it is
`pipeline/gf_grades.py` — RYA-799's grading — and its verdicts are precisely what the
re-verification below re-derives. Nothing else in the repo can move.

---

## Re-verifying RYA-799's SCALE-MISMATCH (spec item 4)

The same grading code, run twice, told only a different gf:

| verdict | on the stale column (RYA-799 as published) | on the corrected column |
|---|---|---|
| `GF-LAB` (graded) | 2 | **25** |
| `systematic:K07` | 221 | 235 |
| `systematic:K07/SCALE-MISMATCH` | **48** | **11** |

**Of RYA-799's 48 SCALE-MISMATCH lines: 11 are genuine, 37 were pure metadata artifacts.**

And the graded count goes from 2 to **25** — Ruffoni 2014 (15) and Den Hartog 2014 (10) —
a twelvefold recovery of lines that were being misclassified purely because the column
lied about what they were measured on.

Effect on the budget RYA-783 consumes: pool RMS σ **0.1993 → 0.1911 dex** against the
0.2000 blanket. Still modest, because 246 lines remain systematic-only — the σ story RYA-824
established is unchanged. What changed is that 23 more lines now carry a real bar instead
of a blanket.

---

## What this ticket did NOT fix

**The `*_PASS.csv` pool writer is not in this repository.** The pool carrying the stale
column is a gitignored artifact, and no committed script produces it — the `panel` column
points at appendix-era tooling that has since been replaced. So the fix is applied at the
**source**: the accounting table is honest, and any pool rebuilt from it inherits honest
columns. For pools that already exist,
`scripts/rya825_reverify_scale_mismatch.py` emits a corrected copy
(`rya825_pool_gf_corrected.csv`) without re-measuring anything.

That is a real gap and it is named rather than papered over: a measured pool whose
generator does not exist cannot be regenerated, which is the RYA-559 class the
`GENERATORS.yaml` manifest exists to prevent — here applied to a gitignored intermediate
rather than a committed result.

---

## The drift guard (spec item 3)

`gf_resolver.assert_gf_column_is_honest` re-derives the resolver's answer for every row a
table claims is canonical and raises on the first disagreement. It is asserted in
`tests/test_gf_column_honesty_rya825.py` **against the live artifact, not a fixture** —
a fixture cannot drift, and drift is the whole failure mode.

Nine tests, including the two that matter most: the guard **catches a column nudged by
0.05 dex**, and it **refuses a table with no `gf_canonical` flag** rather than passing
vacuously on a table it cannot check.

---

## Owed

1. **Apply the same column honesty to the non-Fe pools** as they are rebuilt — the
   accounting fix covers the inventory, but each element's measured pool carries its own
   copy. Feeds RYA-709.
2. **The missing `*_PASS.csv` generator** deserves its own ticket: an artifact that
   nothing can regenerate is the defect `GENERATORS.yaml` was written for.
3. **RYA-799's published numbers are superseded** by the re-verification above. Its
   document already carries the RYA-824 correction; the counts here supersede its
   `graded 2 / systematic-only 269`.
