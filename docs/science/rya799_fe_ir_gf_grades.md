# RYA-799 — grading, or honestly bounding, the 271 Fe I IR lines

**Pool:** `FeI_6910_9199_kpno_solar_atlas_PASS.csv`, 271 lines, 6911.5–9173.2 Å (Sirius;
three copies, all md5-identical `bb9d233e2223`).
**Ratified principle executed** (Ryan, 2026-08-12): *"we still do our best job, and if
there is no grade, no way to tie it to a graded line, then we note the systematic."*

Reproduce:

```
python3 scripts/rya799_fetch_fe_gf_lab.py
python3 scripts/rya799_grade_fe_ir_pool.py --pool <PASS pool>
```

---

## The headline, and it is not the one the ticket expected

**36 of the 271 lines have a primary laboratory gf in the literature. The pool was
measured on that value for 2 of them.**

The ticket's premise was that grades exist but were *flattened at intake* — VALD3's
per-line bibliographic provenance collapsed to the string `"VALD3"` — and that recovering
the provenance would recover the grades. The first half is true and is now undone. The
second half is not: un-flattening recovers a **source**, and for Fe I in this repo that
source is uniformly Kurucz. The grading gap is not a bookkeeping loss. It is that the
pool was built on semi-empirical values while better ones sat in the repo unused.

---

## What the provenance recovery actually found

**No re-extract was needed.** The per-line gf source was never lost: VALD3's long format
carries it on the **fourth physical line of every transition**, and
`data/linelists/vald_parse.py` reads only the first and discards the rest. The flattening
happened in the CSV builder, downstream of a file we already hold. Parsing that reference
line gives, for every Fe I transition in the solar delivery:

| band | Fe I transitions | gf source |
|---|---|---|
| 3780–6910 Å | 5912 | **100 % `gf:K14`** (Kurucz 2014) |
| 6910–9199 Å | 1696 | **100 % `gf:K14`** |

Not one Fe I line in either band carries a primary source in this delivery. That is not
an IR-specific defect — the optical band is in exactly the same state — and it is why
RYA-780 found no primary measurement behind any of its fourteen lines.

⚠️ **`gf:RU` in this file is Raassen & Uylings, not Ruffoni.** It appears 829 times, all
on Fe II and Cr II, and the file's own reference legend disambiguates it. Reading it as
Ruffoni would have manufactured 829 graded lines out of nothing.

### The surface the ticket did not know about

`data/linelists/canonical_gf.csv` **already carries a real per-line `loggf_reference`**
for this band — 5243 Fe I rows in 6910–9199 Å, of which 5049 are `K07` and **194 name a
compilation or primary source** (`FMW`, `BWL`, `GESHRL14`, `2014MNRAS.` — a truncated
Ruffoni bibcode — `MRW`, `BK`, `BKK`, `GESB82c` …). Its `nist_grade` column is entirely
empty for these rows.

So the provenance question was already answered in-repo. **The value question was not.**

---

## ⚠️ CORRECTION (RYA-824, 2026-08-15) — read this before the section below

**The grade verdicts below stand. The explanation of *why* does not, and it is corrected
here rather than quietly edited away.**

This document says the measured pool "was measured on VALD K14 gf" while better values sat
unused in `canonical_gf.csv`. RYA-824 tested that directly by re-inverting and found it is
**not what happened**. `abundances_derive._load_synth_resources` defaults to
`apply_canonical_gf=True`, which runs `gf_resolver.apply_to_synth_array` (RYA-353), so the
production synthesis line list **already carries whatever `canonical_gf` holds** — lab gf
included. Of 29 lab-covered lines, **18 were already on the lab value**; re-substituting it
changed the abundance by exactly 0.000.

What is K14 is the **`log_gf` column** in the measured pools and in
`data/audit/line_accounting/per_line.csv`. That column is stale metadata; it does not
describe the gf the inversion used, and this document compared against it.

Therefore:

* **The 48 `SCALE-MISMATCH` lines are a provenance-LABELLING defect, not 48 abundances
  derived from the wrong oscillator strength.** The count is right, the severity is not.
* **The proposed follow-up — "re-invert the 48 on the referenced gf" — is largely
  unnecessary.** Most of them are already there. The real fix is to correct the stale
  column.
* **The grade verdicts are unaffected**, and for the reason this document gives: a grade
  must describe the number that was used, and a column that does not describe the
  computation cannot certify it. Right verdict, wrong diagnosis.
* **The headline "36 of 271 have a primary lab gf; the pool used it for 2" is measured
  against that same stale column** and understates the truth. Against the gf the inversion
  actually used, the figure is 18 of 29 on the lines RYA-824 could test.

The corrected one-line summary of both tickets: **the Kurucz floor is in the σ we assign,
not uniformly in the values.**

---

## The trap: a grade describes the number that was USED

Matching the 271 PASS lines into `canonical_gf` on (λ ±0.01 Å, EP ±0.01 eV) ties 236 of
them, 110 to a non-`K07` reference. It is tempting to stop there and report 110 graded
lines. That would be wrong, and quietly so.

**The pool's `log_gf` and `canonical_gf`'s disagree by up to 0.96 dex.** Split by whether
the referenced value is the one the pool actually used:

| reference tag | value matches the pool | total |
|---|---|---|
| FMW | **39** | 39 |
| BWL | **16** | 16 |
| GESHRL14 | 0 | 13 |
| `2014MNRAS.` (Ruffoni 2014) | **0** | 13 |
| MRW | 0 | 6 |
| BWL+GESHRL | 0 | 3 |
| BK / BKK / GESB82c / … | 10 | 20 |

Every Gaia-ESO- and Ruffoni-updated value in the repo is one the pool did **not** use.
Attaching those references' accuracy to abundances derived from Kurucz values would
fabricate a pedigree: the number would look graded and would not be. So the harness
confirms every tie against the log gf the pool was measured with, and a tie that fails
that check becomes `SCALE-MISMATCH` — same bar, plus a flag, because *"a better gf exists
in-repo and was not used"* is actionable and must not be silently absorbed.

---

## The three terminal states, and the smoke test

| state | n | σ |
|---|---|---|
| `GF-LAB` — primary lab measurement, and it is the value the pool used | **2** | the paper's own per-line σ (0.02–0.03 dex) |
| `systematic:K07` — no graded tie | **221** | 0.20 dex (RYA-161) |
| `systematic:K07/SCALE-MISMATCH` — a better-referenced gf exists in-repo, but is not what the pool used | **48** | 0.20 dex + flag |

271 in → 271 out, **0 blank bars**, `graded 2 / systematic-only 269`.

Inside the plain-systematic bucket, two very different debts:

- **66 lines** whose gf *value* is confirmed against a named non-`K07` reference
  (`FMW` 39, `BWL` 16, …) but which cannot be **graded** because those compilations
  publish no per-line σ we can cite while NIST ASD is down. **Recoverable** the day ASD
  returns.
- **155 lines** with no reference tie at all. A genuine measurement gap.

Per-line σ for RYA-783's budget: graded lines mean **0.025 dex**, a factor 8 tighter than
the blanket. Pool RMS σ **0.1993 dex** against the 0.2000 blanket — i.e. **at 2 graded
lines the budget barely moves.** That is the honest answer: the grading work does not
shrink the Fe IR bar today, and saying so is the result.

---

## Axes that returned nothing, and the controls that prove they were asked

- **NIST ASD — externally blocked, re-verified.** `lines1.pl` returns HTTP 500,
  `Can't use an undefined value as an ARRAY reference`, on the exact RYA-760 recipe
  (tested 2026-08-15 for Fe I 6910–9199). Broken since at least 2026-08-11. This was the
  intended grade-letter axis. It is also the weaker one: RYA-760 established that `FMW`
  *is* a NIST compilation and VALD copies it, so ASD agreement proves only that nobody
  mistyped a transcription.
- **GES quality flags — zero coverage, exhaustively.** The ticket calls these "the
  cheapest in-hand grading axis". Every one of the **62** iSpec `42000_GES` and
  `42000_VALD` region files was swept: they carry 100–198 Fe I lines each, and **not one
  Fe I line above 6910 Å**. The GES "good for params" lists stop before this band. The
  positive control is in the same numbers — the files are full of Fe I, just all optical.
- **`data/linelists/nist_reference.csv`** — 40 rows, 10 of them Fe, spanning 5250–6648 Å.
  Zero in this band.

---

## RYA-711 dependency — satisfied, and deliberately not reused

RYA-711 items 1+2 are **merged** (`dc49658`, PR #243): our measurement grade is
`MQ-A/B/C/D`, with the prefix in the values, and NIST's letters stay unprefixed.

The ticket asks to "single-source the grade vocabulary through RYA-711's subject-named
grade". I followed the **principle** and not the literal value set: `MQ-` means
*measurement quality* — how well we measured a line — and a gf grade is a claim about the
**atomic datum**. Reusing `MQ-` here would be a category error wearing the costume of
single-sourcing, and it would put two different claims on one scale. So gf verdicts are
`GF-`/`systematic:` prefixed: they name their subject, as RYA-711 requires, and a test
asserts they are disjoint from both the `MQ-` set and NIST's unprefixed letters.

A pleasing corroboration of RYA-711 item 2: Belmonte et al. 2017 Eq. 1 is
Δlog(gf) = log(1 + ΔA/A) — the same percent→dex bridge RYA-711 derived independently.

---

## Owed / follow-ups

1. **Re-invert the 48 SCALE-MISMATCH lines on the referenced gf.** This is the actionable
   one and it is cheap: the values are already in `canonical_gf.csv`. 13 of them are
   Ruffoni 2014 measurements with published σ of 0.02–0.13 dex, so they would move from a
   0.20 blanket to a real bar. **This is a pool-rebuild, not a grading job** — it belongs
   with RYA-161/RYA-398 curation, and it should be a ticket rather than a silent fix here
   (changing the gf changes the abundance, and this ticket must not tune one).
2. **Re-run the ASD leg when the endpoint returns** — it converts 66 attributed lines
   into graded ones with no new measurement.
3. **RYA-379** (extend `canonical_gf` past 9199 Å) is the same gap one band redward, and
   **RYA-822** is its blueward twin.
