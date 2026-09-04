# RYA-853 — the Fe II arbiter lines claimed a NIST grade NIST never gave them

**Outcome: the claim is removed, the value is held, and the absence of a primary lab gf is
declared. No abundance moves.**

## The defect

`canonical_gf` labelled Fe II **6149.246** and **6247.557** `NIST ASD v5.11 grade B`. Live
NIST ASD — queried in **air**, **EP-matched on both sides**, uniquely matched — says:

| line | ours | NIST | Δ | NIST acc. |
|---|---|---|---|---|
| 6149.246 | −2.724 | **−2.854** | +0.130 | **E** |
| 6247.557 | −2.329 | **−2.444** | +0.115 | **D** |

Wrong on both axes, and **the grade is the load-bearing half**: `B` sits in
`NIST_GRADE_HIGH` (*"lab, ≤10 %, trusted"*) while E and D sit in `NIST_GRADE_CULL`. A
fabricated B publishes 0.041 dex on a line whose own source says 0.176–0.301.

## It lived in three files, propagating downhill

```
nist_reference.csv / nist_crosscheck.csv     value AND grade originate here — hand-maintained
          │                                   extracts whose own header warns rows "should not
          │                                   be assumed current", and whose grades RYA-853
          │                                   scope 1 measured 70 % wrong
          │  build_linelist.crosscheck_nist()
          ▼
linelist_solar.csv                           VALD3 value (−2.72) wearing the stamped B
          │  migrate_gf_single_source.py
          ▼
canonical_gf.csv                             adjudicated to the extract: −2.724 + grade B
```

Fixing `canonical_gf` alone would have been undone by the next build.

## Why the value is held rather than corrected

Every alternative source fails, and each was **checked, not assumed** (RYA-833):

- **Den Hartog 2019** — pure lab, and the referee scope 3 used. Its table **stops at 4584 Å**;
  these lines are redward. No coverage.
- **Meléndez & Barbuy 2009** — has both (−2.69 / −2.30), flagged **`S`**: solar-fitted.
  **Firewalled** by RYA-161 — adopting it would be the circularity this ticket exists to
  detect.
- **NIST ASD** — covers both, at **grade E and D**. And RYA-853 scope 3 measured our Fe II
  scale against pure-lab DH19 and found **NIST is the low one**. Adopting NIST's value would
  move the ionization arbiter ~0.13 dex onto the scale we just measured as low, on grade E/D
  data.

So there is no adjudicable source. The number stands as what it is — a VALD3-scale value —
the false NIST claim is removed, and `gf_sigma_dex` stays **uncited**, so the band budget
charges the ungraded blanket. That is the honest state for a line with no graded source, and
inventing a σ would be a second fabrication. **Owner for closing the gap: RYA-953** (Fe II has
no primary-lab gf table).

## What changed — exactly eight lines

| file | rows | fields |
|---|---|---|
| `canonical_gf.csv` | 2 | `loggf_reference`, `nist_grade` → blank, `adjudication_status` → `held_rya853` |
| `linelist_solar.csv` | 2 | `nist_grade` → blank |
| `nist_reference.csv` | 2 | `nist_grade` → blank, `notes` |
| `nist_crosscheck.csv` | 2 | `nist_grade` → blank, `notes` |

**`log_gf` is untouched in all four files.**

⚠️ Two near-misses worth recording, both caught by checking the blast radius rather than
trusting the edit:

1. The first cut used pandas to read and write each file. That **reformatted 46 unrelated
   `canonical_gf` rows** (`-1.2500` → `-1.25`) and rewrote both extracts wholesale. A
   provenance fix that silently reformats rows it was not asked to touch is the same class of
   defect it is fixing. Replaced with line-level surgery.
2. `nist_crosscheck.csv` has **mixed line endings** (21 CRLF, 5 LF). `Path.read_text()`
   normalises them on read, so the write-back rewrote 19 more untouched lines. Now read and
   written with `newline=''`.

## Scope 4 — both wavelength-only matchers closed

`build_linelist.crosscheck_nist()` matched on element + ion + **wavelength alone** within
0.010 Å and took the **first hit in file order**. That is how the fabricated `B` reached
`linelist_solar.csv`. It now requires the EP to agree within 0.05 eV and **refuses ambiguity**
rather than resolving it by proximity.

🔴 **This is not theoretical — it catches a live mis-stamp.** Measured over the real files:

| | stamps |
|---|---|
| old wavelength-only matcher | 74 |
| EP-guarded matcher | **73** |
| refused, EP disagrees | **1** — `Fe I 6065.490`, NIST-row EP **2.608** vs our line's **4.956** |

`rya347_fe2_atomic_data_audit.py` carried the same defect at ±0.1 Å. Guarded the same way —
and `step0` now emits `ep_eV`, because a guard fed NaN would refuse everything and read as
"no graded source", which is a silent false absence rather than a fix.

## Still open on this ticket

- **Scope 1 beyond Fe II**: 43 of the 45 tabulated corrections remain unapplied (the two Fe II
  grade rows are done here). The extract-vs-extract disagreements — including `Fe I 6065.49`,
  which this pass independently flagged from the other direction — are still live, and both
  strict xfails still xfail.
- **Scope 2**: ✅ **CLOSED, and the open item was withdrawn.** Belmonte Fe I **3935.3064**
  was carried here as −2.199 / σ 0.070 against "the paper's" −1.820 / 0.180. That −1.820 is
  the *Published* column — May et al. (1974), which Belmonte tabulates for comparison — not
  Belmonte's own measurement. Ours reproduces *This Experiment* exactly. Re-refereed with
  the column read positionally: **464 of 465 lines, zero mismatches on either axis.**
