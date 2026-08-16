# RYA-836 — the near-UV primary-lab-gf Fe sub-pool

**Star:** Sun, Kitt Peak flux atlas. **Band:** 3000–3780 Å, synthesis-only.
**Rule:** primary laboratory gf only — never a compilation, never solar-calibrated.

```
python3 scripts/rya836_nearuv_lab_gf_subpool.py --census-only
python3 scripts/rya836_nearuv_lab_gf_subpool.py                          # Sirius
python3 scripts/rya836_nearuv_lab_gf_subpool.py --from-per-line <csv>    # re-report
```

---

## The headline: the ticket's premise is refuted, and there is still a win

The ticket says *"the scatter IS the Kurucz floor"* and *"this is what actually tightens
the near-UV."* **The first is wrong and the second is wrong about which quantity.**

**Switching 60 near-UV lines from the production gf to a primary laboratory gf changes
their scatter by +0.002 dex.** Same lines, same route, one input different:

| | median A(Fe I) | scatter |
|---|---|---|
| control — production gf | 7.594 | **0.651** |
| primary lab gf | 7.607 | **0.654** |

Median ΔA is **+0.000 dex**. The gf is not what makes the near-UV scatter.

### What does make it: line selection

The RYA-832 cell scatters 0.413 over 40 lines chosen by theoretical depth with a 4.0 Å
minimum separation. This sub-pool scatters **0.651 over the same production gf** — because
its membership is *whatever the laboratories happened to measure*, not what is cleanly
measurable in a crowded band.

> **Line selection accounts for 0.238 dex of scatter. The oscillator strength accounts for
> 0.002.**

---

## The win that is real, and it changes what to do next

| | RYA-832 Kurucz cell (n=40) | RYA-836 lab-gf sub-pool (n=59) |
|---|---|---|
| line scatter | **0.413** | 0.652 |
| stat | 0.0652 | 0.0849 |
| **syst** | **0.1972** ~~0.2211~~ | **0.1081** ~~0.1472~~ |
| **dominant term** | **gf scale (UNGRADED)** | **pseudo-continuum** |
| total (stat ⊕ syst) | 0.208 ~~0.231~~ | **0.137** ~~0.170~~ |

Two things happened that are worth separating:

1. **The systematic fell 0.197 → 0.108** (corrected — see below), because the pool is now on gf with a cited
   per-line σ (0.020–0.120 dex) instead of an ungraded 0.20 blanket.
2. **The dominant term flipped from gf to the pseudo-continuum.** That is the actionable
   result: *once the near-UV is on primary laboratory gf, gf stops being the limiting
   systematic.* The next lever on this band is the 0.100 dex pseudo-continuum, not more gf
   work.

So: **does it beat ±0.413?** On *scatter*, no — 0.652. On *total reported uncertainty*,
yes — 0.170 against 0.231. The ticket's question conflates two quantities that moved in
opposite directions, and both answers are above.

---

## Relationship to RYA-824 — same operation, different outcome, same lesson

RYA-824 did this in the IR/VIS and the gf systematic fell 0.20 → 0.05. That is the same
effect seen here (0.197 → 0.108): **attaching a citable per-line σ**. What RYA-824 did *not*
do, and this ticket makes explicit, is reduce line-to-line scatter — 824's value moved
≤0.026 dex and its scatter was essentially unchanged too.

The consistent finding across both: **primary lab gf buys you a smaller error BAR, not a
smaller SPREAD and not a different value.** RYA-822 showed grading alone buys neither.

---

## The line-identification screen

One line was excluded from the aggregate and **carried, not dropped** (RYA-711):

**Fe I 3026.056** — production gf −5.038 against Belmonte 2017's −2.077, a **2.961 dex**
gap. That reproduces the maximum the ticket quotes.

A primary lab σ in this band is ≤0.12 dex and a Kurucz semi-empirical gf is good to perhaps
0.2–0.3. A whole-dex disagreement is outside any combination of the two, so the likelier
explanation is that the lab line and the line-list line are **different transitions** that
agree in wavelength and excitation potential to within the match tolerance. Substituting
there would put one line's gf on another line — the RYA-785 *"a same-species neighbour does
not cancel"* failure — rather than adopt a better oscillator strength.

The 1.0 dex threshold was fixed **before** looking at what it excluded, and the sub-pool is
reported both ways so the screen's leverage is visible rather than asserted:

| | median | scatter | n |
|---|---|---|---|
| screened (reported) | 7.577 | 0.652 | 59 |
| unscreened | 7.607 | 0.654 | 60 |

The screen moves the value by 0.030 dex and the scatter by 0.002. It is not doing heavy
lifting — which is itself worth knowing, and is why it is reported rather than trusted.

Only **1 of 61** in-band lab lines sits within 0.05 Å of another, so ambiguous matching is
not systemic here.

---

## Census

Reconciles the ticket exactly.

| | count | sources |
|---|---|---|
| primary-lab Fe I below 3780 Å | **105** | Belmonte 70, Den Hartog 31, Ruffoni 4 |
| **in band 3000–3780 Å** | **61** | Den Hartog 31, Belmonte 26, Ruffoni 4 |
| cited σ | 0.020–0.120 dex (median 0.030) | |

My RYA-824 census reported 64 for the same sources; that was a 3000–**3800** window and
the 3 extra sit in 3780–3800. No discrepancy.

**`GF-NIST` is excluded**, and the exclusion is load-bearing. RYA-822 assigned it to 604
near-UV lines, but `FMW` *is* a NIST compilation and VALD copies it (RYA-760) — admitting
it would let a compilation restate the headline it exists to test. Solar-calibrated gf is
excluded by the RYA-161 firewall.

---

## Spec item 4 — deliberately NOT done as written

The ticket asks to *"feed the sub-pool into the RYA-832 matrix cell as the TIGHT near-UV
product."* **It is not the tight product on the quantity that phrasing means.** Its line
scatter is 0.652 against the 832 cell's 0.413, so installing it as "the tight one" would
assert the opposite of what was measured.

What is done instead: it is emitted as its own matrix cell under the `1D-LTE-LABGF`
treatment, side by side with the Kurucz cell and never merged (RYA-712), with its own
value/σ/n. A reader can then see both — the broad-but-cleanly-selected 832 cell and the
lab-gf cell with the smaller systematic — which is the comparison that actually exists.

**The 832 cell remains the near-UV product of record.**

---

## Errors hit

`_report` crashed **after all 122 Turbospectrum fits had completed**: `~series.fillna(False)`
on an object-dtype column is bitwise negation on ints, turning `True` into `-1`, which was
then used as a column key. Nothing was lost, because the per-line CSV is written *before*
the report and `--from-per-line` rebuilds from it. That ordering is not luck — it is why
the expensive half and the cheap half are separated, and it is what let the
line-identification screen be added mid-flight without discarding the fits.

---

## Owed

1. **The pseudo-continuum is now the dominant near-UV systematic** (0.100 dex, and it does
   not average down). That is the next lever on this band, and it is a normalisation
   problem rather than an atomic-data one.
2. **A depth-selected lab-gf sub-pool** would separate the two effects properly: restrict
   to lab-covered lines that *also* pass the 832 selection rule. There are few of them,
   which is itself the point — the labs and the clean-line criterion barely overlap here.
3. **Fe I 3026.056 wants adjudication** — either VALD's gf is badly wrong or the lines are
   misidentified. It is carried and flagged, not resolved.

---

## ⚠️ CORRECTED BY RYA-845 — the systematic figures were double-counted

The `syst` values this document originally published were **inflated**. `error_budget.build()`
already carries the 0.100 dex pseudo-continuum term for the near-UV — the band's policy says
"pseudo-continuum only", which fires the branch in `error_budget.py` — and both product
routes then added it **again** in quadrature.

| cell | published here | correct |
|---|---|---|
| RYA-832 near-UV (product of record) | 0.2211 | **0.1972** |
| RYA-836 lab-gf sub-pool | 0.1472 | **0.1081** |

**Nothing about the conclusion changes, and the direction of the finding is unchanged.** The
fall on moving to primary laboratory gf is **larger** than was claimed here (0.197 → 0.108,
not 0.221 → 0.147), and the "dominant term flipped from gf to pseudo-continuum" statement is
untouched — `dominant` is computed from the budget's own terms, before the stray addition.

The abundances, line counts, scatters and `stat` values in this document are all unaffected;
only the systematic — and the totals derived from it — were wrong.

Found by RYA-841 while asking where the 0.147 came from. The defect originated in RYA-832's
route and was inherited here, and RYA-832's unit test asserted the constant *equals* 0.100,
which pinned the double-count rather than catching it.
