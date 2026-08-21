# RYA-945 — the Fe I / Fe II primary-laboratory gf backbone in `canonical_gf`

**Status:** landed. **Artifacts:** `data/audit/rya945_fe_lab_gf/`.
**Scripts:** `scripts/rya945_ingest_fe_lab_gf.py`, `scripts/rya945_fetch_fe2_gf_dh19.py`.
**Test:** `tests/test_fe_lab_gf_ingest_rya945.py`.

## Why

The 148-line VIS Fe I pool reads 7.586 while the 9-line lab-gf pool reads 7.445
(RYA-819/831). That 0.14 dex is atomic data, not measurement — RYA-855 moved 0 of 36 bars
because every Fe cell is a mixed pool sitting on the Kurucz floor. RYA-824 proved the lever
is real *and* that it is bounded by **our coverage, not by lab availability**: the lab
tables hold 250 Fe I VIS lines and our measured pool reached nine of them.

This closes that gap in the **line list**. It re-derives nothing.

## What landed

`canonical_gf.csv` keeps all 167,739 rows — this is an **update**, not the RYA-834-style
extension — and rewrites 1,597 of them.

| tier | rows rewritten | median cited σ | vs the 0.20 blanket |
|---|---|---|---|
| LAB (DH14 / RU14 / BEL17 / DH19) | 395 | **0.030 dex** | 6.7× tighter |
| NIST grade C or better | 1,202 | **0.041 dex** | 4.8× tighter |

Fe I+II lines carrying a lab-or-graded gf at EP ≥ 1.2 eV: **475 → 1,826**.

### Per source

| tag | source | access | in table | matched | adopted |
|---|---|---|---|---|---|
| `DH14` | Den Hartog+ 2014, ApJS 215, 23 · `10.1088/0067-0049/215/2/23` | VizieR J/ApJS/215/23 | 203 | 199 | 199 |
| `RU14` | Ruffoni+ 2014, MNRAS 441, 3127 · `10.1093/mnras/stu780` | VizieR J/MNRAS/441/3127 | 142 | 129 | 129 |
| `BEL17` | Belmonte+ 2017, ApJ 848, 125 · `10.3847/1538-4357/aa8cd3` | arXiv:1710.07571 source (J/ApJ/848/125 is a 404) | 120 | 74 | 73 |
| `DH19` | Den Hartog+ 2019, ApJS 243, 33 · `10.3847/1538-4365/ab322e` | published PDF — **not on VizieR** | 131 | 22 | 22 |
| `NIST-C+` | NIST ASD, `astroquery.nist` | 3000–12935 Å Fe I, 3000–10502 Å Fe II | 2,180 graded C+ | — | 1,202 |

## The judgements

### Precedence is a one-way ratchet

`LAB > NIST grade C+ > Kurucz K## > VALD3`, written once in `_TIER_RANK`. A rewrite must
strictly improve the tier; 560 candidates were refused as downgrades or lateral moves. The
test re-derives this from the written file, not from the run's own report.

### The wavelength window is measured, not chosen

Binning the 2,023 graded NIST matches by wavelength residual and asking how often the two
sources then disagree by more than 0.2 dex:

| residual (Å) | n | frac \|Δ\| > 0.2 dex | median \|Δ\| |
|---|---|---|---|
| 0 – 0.001 | 1783 | 1.9 % | 0.0003 |
| 0.001 – 0.003 | 184 | 3.8 % | 0.0009 |
| 0.003 – 0.006 | 29 | 3.4 % | 0.0006 |
| 0.006 – 0.012 | 17 | **11.8 %** | 0.0025 |
| 0.012 – 0.020 | 10 | **30.0 %** | **0.1201** |

The last two bins are a different population — the chance-pairing signature RYA-822
diagnosed, where coincidences pile up at the tolerance *wall* while real matches sit at the
fourth decimal. The window is set inside that separation, at the ticket's 5 mÅ, which the
data independently vindicates. Zero ambiguous matches at that setting; randomised-null
match rate 0.04 % (Fe I) and 0.0 % (Fe II) against 86.5 % / 16.8 % real.

### RYA-834's band-specific ruling overrules the ticket's blanket precedence

In 9199.9–12935 Å, RYA-834 measured NIST ASD against the band's own line list: 34 of 56
matched, median |Δ| 0.0003 dex. Agreement at the fourth decimal is the same number arriving
twice (RYA-760: FMW *is* a NIST compilation and VALD copies it), so there NIST is
**recorded, never adopted**. Twenty-four rows take the grade and its σ without their value
moving. A band-specific measurement outranks a general rule.

### The result on values is very largely a null — and that is the finding

Per TP reference the NIST-minus-incumbent delta has median 0.000 with MAD 0.000–0.002. What
this ingest buys is not different numbers, it is **cited per-line uncertainties**: RYA-824's
"the Kurucz floor is in the σ we assign, not uniformly in the values", now actionable
across 1,600 lines instead of 29. Fifty-eight rows do move by more than the 0.20 dex Kurucz
systematic; every one is named in `outlier_ledger.csv`, largest +4.649 dex (Fe I 5487.738,
an intercombination line where Kurucz's semi-empirical value and NIST's graded experimental
one genuinely disagree).

### `RU` is Raassen & Uylings, not Ruffoni

The ticket calls the 5,113 `RU` rows a partial Ruffoni-2014 wiring. They are all **Fe II**
(and Cr II), 4200–9199 Å; Ruffoni 2014 is an Fe I paper. VALD's own reference block reads
`Fe 2: Raassen & Uylings … gf:RU` — a theoretical orthogonal-operator calculation. Ruffoni's
Fe I values reach `canonical_gf` under the GES codes `GESHRL14` and `2014MNRAS.` instead.
Treating `RU` as lab-already-wired would have credited 5,113 theoretical rows as primary
laboratory data. Pinned by test.

### Two pre-existing defects, refused rather than papered over

Nine rows claimed a NIST grade this run could not stand behind, and were demoted to
`OTHER` — value and note kept, claim withdrawn (`demoted_nist_rows.csv`):

* **five** sit on `hfs_n_components == 2` clusters (RYA-822). A single-transition ASD gf
  cannot describe a gf-summed two-component row — the same HFS error this script refuses to
  make on the way in, arriving from a previous ticket.
* **four** carry a stored `NIST ASD v5.11 grade A/B` letter with no verified σ (RYA-354).
  RYA-852 already showed two of them are stale: 6149.246 and 6247.557 are stored as B while
  live ASD reports E (0.301) and D (0.176) — understated 7.3× and 4.3×. Re-deriving a σ
  from the stored letter would launder that staleness into a number.

## Limits, stated

* **Fe II lab coverage is 22 lines and the binding constraint is our blue edge, not the
  source.** 109 of DH19's 131 lines lie blueward of `canonical_gf`'s 3000.06 Å edge. Of the
  22 that are reachable, 22 matched — a 100 % hit rate on what could be reached. Extending
  the table below 3000 Å would make ~109 more available.
* **63 Fe I lab lines found no canonical row** (BEL17 46, RU14 13, DH14 4) — mostly Belmonte's
  deep-UV lines below 3000 Å, same edge.
* **No canonical row is measured by more than one lab paper**, so the ticket's requested
  inter-source spread diagnostic has an empty population. The three Fe I sets are disjoint
  by construction: DH14 is high-lying *even*-parity, BEL17 high-lying *odd*-parity, RU14 the
  GES set.
* **The diagnostic flag carries two of the ticket's three criteria.** `blend_flag` is a
  column of the *measured EW table*, set per observation by the fitter and applied by
  `abundances_derive._prefilter`. A line list cannot carry it. The line-level criterion used
  here is the problem-children registry's standing `exclude` rows; the measurement-time
  filter still runs downstream.
* **Nothing was re-derived.** Whether the reported Fe value now lands in the literature band
  is the next ticket's measurement, not this one's claim.
