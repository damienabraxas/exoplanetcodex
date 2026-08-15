# RYA-824 — the primary-lab-gf Fe sub-pool

**Star:** Sun, Kitt Peak flux atlas, RYA-783 per-line pools.
**Rule executed:** only PRIMARY LABORATORY gf. Never astrophysical or solar-calibrated —
that would make the solar abundance a restatement of its own input (RYA-161 firewall).

Reproduce:

```
python3 scripts/rya824_lab_gf_subpool.py --pool-dir <RYA-783 band_products> --census-only
python3 scripts/rya824_lab_gf_subpool.py --pool-dir <RYA-783 band_products>       # Sirius
python3 scripts/rya824_lab_gf_subpool.py --pool-dir <dir> --from-per-line <saved csv>
```

---

## The headline: the gain is the ERROR BAR, not the value

| band | pool | A(Fe I) | line scatter | n | gf systematic |
|---|---|---|---|---|---|
| IR | PRODUCTION (canonical_gf) | 7.542 | 0.110 | 20 | 0.200 |
| IR | **PRIMARY-LAB-GF** | **7.516** | 0.101 | 20 | **0.052** |
| VIS | PRODUCTION (canonical_gf) | 7.445 | 0.111 | 9 | 0.200 |
| VIS | **PRIMARY-LAB-GF** | **7.445** | 0.118 | 9 | **0.060** |

Re-measuring on primary laboratory gf cuts the gf systematic by a factor **3.8 (IR)** and
**3.3 (VIS)** — 0.20 dex down to 0.05–0.06. It moves the *value* by **−0.026 dex in the
IR and 0.000 in the VIS**.

That is the ticket's lever, delivered, and it is smaller than the ticket expected. The
reason is worth more than the number.

---

## Why the value barely moved — and the premise that needed correcting

**18 of the 29 lines came back with ΔA exactly 0.000.** Not a failed substitution: a
**no-op**. `abundances_derive._load_synth_resources` defaults to
`apply_canonical_gf=True`, which runs `gf_resolver.apply_to_synth_array` (RYA-353 gf
single-sourcing), so **the production synthesis line list already carries the primary lab
gf wherever `canonical_gf.csv` holds one.** Those lines were never Kurucz-floored in
value. Substituting the lab gf into a line list that already had it changes nothing, and
0.000 is the correct answer.

What *is* Kurucz K14 is the **`log_gf` column** in the measured pools and in
`data/audit/line_accounting/per_line.csv`. That column does not describe the gf the
inversion used. It is metadata that drifted from the computation.

### This corrects RYA-799

RYA-799 compared `canonical_gf` against the pool's `log_gf` **column** and read the
disagreement as *"the pool was measured on Kurucz while a better gf sat unused."* The
disagreement is real; the conclusion was wrong. The inversion used `canonical_gf`. So:

- The 48 `SCALE-MISMATCH` lines are a **provenance-labelling** defect — a stale metadata
  column — not 48 abundances derived from the wrong oscillator strength.
- RYA-799's proposed follow-up, *"re-invert the 48 on the referenced gf"*, is largely
  unnecessary. On the lab-covered subset, **18 of 29 were already on it.**
- RYA-799's *grade* verdict survives intact and for the reason it gave: a grade must
  describe the number that was used, and a column that does not describe the computation
  cannot certify it. The verdict was right; my explanation of why was not.

The corrected statement: **the Kurucz floor is in the σ we assign, not uniformly in the
values.** RYA-799 applied a 0.20 blanket because no per-line uncertainty was citable. This
ticket makes one citable for 29 lines.

---

## What the 11 genuinely-changed lines show

| | |
|---|---|
| ΔA median | **−0.079 dex** (range −0.229 … +0.018) |
| Δ log gf median | +0.180 dex |
| sensitivity ΔA / Δlog gf | **−0.549** |
| EW range | 35–89 mÅ |

The weak-line limit is ΔA/Δlog gf = **−1.000**: A + log gf is conserved when the line is
optically thin. These lines return **−0.549**, absorbing roughly half the gf change,
because at 35–89 mÅ they sit on the flat part of the curve of growth where equivalent
width responds weakly to abundance.

**This is why the census's first-order estimate over-predicted.** The census projected
ΔA ≈ −Δlog gf ≈ −0.17 dex. The inversion returns roughly half of that on the lines that
move, and zero on the majority that were already on lab gf. A first-order estimate on a
saturated line is not a prediction, and the run is what decides — which is why both are
reported.

---

## The control

Every line was inverted **twice**: once with the line list untouched, once with the lab gf
substituted. The first inversion must reproduce the banked RYA-783 1D-LTE abundance, or
the harness is not the production path and the second number means nothing.

```
CONTROL — the untouched-line-list inversion against the banked RYA-783 value
  n=29  max |A_ctrl - A_banked| = 0.0000 dex  median 0.0000
  PASS — this IS the production path
```

Exact on all 29. The two inversions differ in exactly one input, so their difference is
attributable to the gf and to nothing else.

---

## Decomposing the sub-pool's proximity to the anchor

The sub-pool sits far closer to the gold anchor 7.466 than the full band does, and it
would be easy to read that as the lab gf working. It is mostly **selection**:

| | IR | VIS |
|---|---|---|
| full band, 1D-LTE (RYA-783) | 7.639 | 7.586 |
| the same 29 lines, production gf | 7.542 | 7.445 |
| the same 29 lines, lab gf | **7.516** | **7.445** |
| **attributable to line selection** | **−0.097** | **−0.141** |
| **attributable to the gf swap** | **−0.026** | **0.000** |

Line selection accounts for four to fourteen times as much as the oscillator strengths do.
The lab-covered lines are simply better lines. That is not nothing — it is a real,
defensible sub-pool — but it is not evidence that gf explains the band offset.

### Consequence for RYA-817

RYA-817 concluded that the 7.586 → 7.466 gap "is a gf zero-point matter, not a 3D one."
On the subset where that is testable, **switching to primary laboratory gf does not close
it**: the value moves ≤0.026 dex. That constrains the gf hypothesis on the lab-covered
lines only, and leaves the bulk untested, because the ~240 lines carrying a K07
semi-empirical gf in `canonical_gf` have no laboratory value to switch to. The hypothesis
is not refuted; it is un-testable on today's data for the lines that matter most.

---

## The census, including the zeros

| band | ion | lines | usable | derivable | median Δ log gf | lab σ |
|---|---|---|---|---|---|---|
| VIS | I | 9 | 9 | yes | +0.200 | 0.03–0.10 |
| IR | I | 20 | 20 | yes | +0.157 | 0.02–0.13 |
| VIS | II | 0 | 0 | yes | — | — |
| IR | II | 0 | 0 | yes | — | — |
| near-UV | I | 0 | 0 | census-only | — | — |
| near-UV | II | 0 | 0 | census-only | — | — |

Every pool is listed including the zeros: a band absent from the table would read as *not
looked at* rather than *looked at, nothing there*. Each zero was checked:

- **Fe II, all bands.** The vendored lab table is Fe I only. The Fe II candidate is
  Meléndez & Barbuy 2009 (A&A 497, 611; RYA-472) and it is a **different tier** —
  multiplets normalised on laboratory data, but relative gf *within* a multiplet from
  theory. It clears the RYA-161 firewall (not solar-calibrated) but it is not a primary
  laboratory measurement, and pooling it with one unlabelled would misrepresent both.
- **Fe I near-UV.** 901 measured lines and 64 lab lines share 3000–3800 Å and **not one
  pair falls within 0.02 Å**. Closest approach **0.0651 Å** — three times the tolerance
  with nothing in between, so these are disjoint line lists, not a tolerance artefact. The
  positive control is the IR in the same comparison: 68 pairs match at a minimum
  separation of **0.0000 Å**. RYA-759's 0.354 dex near-UV scatter cannot be attacked with
  the tables we hold.

**Scope, honestly: 29 lines, not 271.** The lab tables hold 250 Fe I lines in the VIS
against 73 in the IR, so the literature *is* optical-rich as the ticket expected — but our
measured VIS pool catches only 62 of those 250, and only 9 survive to the
abundance-producing set, against 20 of the IR's. **The binding constraint is our own
measurement coverage, not lab availability.**

---

## Errors hit

**The scratch-directory guard broke what it was guarding.** I routed the double inversion
to a private Turbospectrum scratch so the two calls could not read each other's working
files (the RYA-785 stale-workdir class) — and did not create the directory. Turbospectrum
writes into it and does not make it; the shared default only ever works because an earlier
run left it behind. All 58 inversions failed with *No such file or directory*. Visible only
because the products table came out **empty**: an inversion that returns NaN for everything
is loud, where a half-working one would have been quiet.

**My own product label was wrong.** The comparison pool was printed as `KURUCZ K14
(sub-pool)`. It is not K14 — it is `canonical_gf`, which is the whole point of this
ticket's finding. Relabelled `PRODUCTION (canonical_gf)`. RYA-711's lesson applied to my
own output: a label that names the wrong subject is worse than no label.

---

## Owed

1. **Fix the stale `log_gf` column** in `data/audit/line_accounting/per_line.csv` and in
   the measured band pools so it reports the gf the inversion actually used. That is the
   real defect RYA-799 found, correctly diagnosed only here.
2. **Fe II needs its own tier.** Meléndez & Barbuy 2009 is worth vendoring for the
   ionisation-balance anchor (RYA-472), under a label that says lab-normalised-with-
   theoretical-branching rather than primary-lab.
3. **The generalisation the ticket asks for** (RYA-709, 26 elements) should test the
   *σ* claim, not the *value* claim: the lever this ticket demonstrates is a factor ~3.5
   on the gf systematic for lab-covered lines, and near-zero on the value.
