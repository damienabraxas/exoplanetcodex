# RYA-1110 — the GBS (Jofré+ 2014) solar Fe reference line set

**Status:** built, not merged. Sister of RYA-1109 (Asplund/AGSS21). Feeds RYA-1111, which
wires `line_set` to actual measurements.

| | |
|---|---|
| Source | Jofré et al. 2014, A&A **564**, A133 — `2014A&A...564A.133J`, bibliography key `jofre2014` (RYA-854) |
| Holding | `data/reference/jofre2014_gbs/` (VizieR `J/A+A/564/A133` + the paper's Tables 4/5, transcribed) |
| Line set | `data/linelists/reference_sets/gbs_solar_fe_rya1110.csv` — 159 rows, `line_set=gbs` |
| Coverage | `data/linelists/reference_sets/gbs_solar_fe_coverage_rya1110.csv` |
| Build | `scripts/rya1110_fetch_jofre2014.py`, `scripts/rya1110_build_gbs_fe_lineset.py` |
| Band | VIS only (every line resolves to VIS through `band_policy.resolve`, not asserted) |

## The set

The GBS solar replication set is `rew_class == "pass"`: **142 lines — 133 Fe I + 9 Fe II.**

| | |
|---|---|
| λ (air) | 4787.83 – 6820.37 Å |
| E<sub>low</sub> | 0.11 – 5.10 eV |
| log(EW/λ) | −5.9339 – −4.8201 |
| solar EW (mean over methods) | 8.2 – 90.7 mÅ |
| golden / not | 122 / 20 |
| published GBS log gf | 121 of 142 |
| joined to `canonical_gf` | 142 of 142 |
| our gf tier | 76 OTHER, 41 NIST-C+, 16 LAB, 9 VALD3 |

The file carries all **159** published solar lines, not just the 142; the other 17 are
retained with their class so the cut can be audited rather than taken on trust
(RYA-931 — quarantine, never cull).

## Two carriers, neither sufficient

Nothing in VizieR gives an oscillator strength. Nothing in the paper's line tables gives a
per-star selection. The line set this ticket asks for exists only as the join:

| | VizieR `J/A+A/564/A133` | paper Tables 4 / 5 |
|---|---|---|
| per-star selection | ✅ `table6.dat`, `ew.dat` | ❌ per stellar GROUP only |
| λ, E<sub>low</sub> | ✅ (λ 2 dp, EP 2 dp) | ✅ (λ 2 dp, E<sub>low</sub> 4 dp) |
| **log gf** | ❌ **absent entirely** | ✅ golden lines only |
| per-method EW / A(Fe) | ✅ six methods | ❌ |

The join is on λ (2 dp on both sides, exact) and is **controlled by EP**: VizieR's 2-dp EP
against the paper's 4-dp E<sub>low</sub> agrees to ≤ 6 meV on all 138 joined lines. Without
that control the gf column would rest on a wavelength alone, which is the RYA-853
`crosscheck_nist` failure exactly.

## 🔴 The paper's own −4.8 cut does not select the set it publishes

Sect. 6.1 states the first selection step verbatim: *"we selected those lines with
log(EW/λ) ≤ −4.8"*, to stay on the linear part of the curve of growth. Table 3 states the
Sun's selected-line counts: N(Fe I) = 150, N(Fe II) = 9. Both `table6.dat` and `ew.dat`
give exactly 150 + 9 for the Sun — and the identity holds for **all 34 benchmark stars**
(checked once against Table 3 read out of the PDF; the standing build guard
`_check_published_counts` pins the Sun's, which is the pair this artifact rests on). So the
published tables *are* the selected set.

Applied to the Sun's own published equivalent widths, **14 of those 159 lines exceed −4.8
on every one of the six methods that measured them.** Fe I 6393.60 sits at −4.63. It is
not a solar quirk:

| | all-method violations | of |
|---|---|---|
| Sun | 14 | 159 (8.8 %) |
| all 34 GBS stars | 713 | 4252 (16.8 %) |
| HD 140283, HD 84937 | 0 | 25, 21 (0.0 %) |
| μ Leo | 48 | 150 (32.0 %) |

The rate tracks stellar type exactly as it would if the cut had not been applied to what is
tabulated. **Neither published statement is overridden here.** The file carries both:
`gbs_selected_sun = 1` on all 159 rows (their selection, as published) and `rew_class`
(their stated rule, applied to their published EWs by this ticket).

## Why `rew_class` has three values

Each line carries up to six independent EW measurements and the paper does not say which
one the cut was applied to. For 156 lines it makes no difference. For 3 it does. Widening
the passing set to absorb them would launder an ambiguity into a decision (RYA-1072), so
there are two closed disjoint sets and an explicit remainder:

| class | rule | n |
|---|---|---|
| `pass` | every method ≤ −4.8 | 142 |
| `excluded` | every method > −4.8 | 14 |
| `ambiguous` | methods straddle −4.8 | 3 — Fe I 5389.48, 5930.18, 6335.33 |

## The join to our pool, and its measured null

λ+EP dual key through `pipeline.line_match` (RYA-1037), `require_ep=True`, tolerance
**0.015 Å**. Derived, not chosen:

* Jofré prints λ to 2 dp, so agreement cannot beat ±5 mÅ from rounding alone.
* Measured worst genuine disagreement: **5.0 mÅ** (Fe I, all 150), **14.0 mÅ** (Fe II, 2 of 9).
* 0.015 Å resolves 159/159 with **zero** ambiguous candidates; at 0.030 Å a genuine fork
  appears, so the window is not sitting on the ambiguity edge.
* **Null control, asserted at build time and in CI:** the same match with the GBS
  wavelengths displaced by ±0.2 / ±0.3 / ±0.5 Å resolves **0** lines at every displacement.
  **The EP half of the key is what makes it 0, and that is measured, not claimed:** drop
  EP and the same displacements find **12 – 22** chance "matches" while the undisplaced
  match falls to **149** of 159. Wavelength alone is both leakier and less complete.

Two positive controls, both from lines that are actually in this set:

* Fe I 6065.48 is `line_match`'s own documented trap — two `canonical_gf` rows 3 mÅ apart,
  EP 2.609 vs 4.956, **1.9 dex apart in log gf**. The join takes EP 2.609, the right one.
* Fe II 6149.26 has **three** candidates inside 0.08 Å (EP 3.889, 13.155, 13.436; log gf
  −2.724, −3.283, −4.983). Wavelength alone would be a 2.3 dex coin flip; EP settles it.

## 🔴 DECISION FLAG for Ryan — their gf or ours? (NOT resolved here)

Both columns ship (`log_gf_gbs`, `log_gf_ours`), plus `gf_synth_ges` for the version
question. Over the 142-line set:

| comparison | n | identical (≤ 0.0005) | median | mean | s.d. | max abs |
|---|---|---|---|---|---|---|
| GBS − ours | 121 | 67 | +0.0000 | −0.0033 ± 0.0040 | 0.0441 | 0.2170 |
| GBS − our GES seed | 120 | 86 | +0.0000 | −0.0023 ± 0.0029 | 0.0315 | 0.1800 |

Three things follow, and they point in different directions:

1. **The mean is not the issue.** Swapping to their gf moves the mean solar A(Fe) by
   +0.003 dex — inside its own standard error. The choice is not a systematic offset.
2. **The scatter is.** 54 of 121 lines disagree, up to 0.217 dex, so the choice
   redistributes ~0.044 dex of line-to-line scatter. Anything that reads per-line residuals
   (ξ, excitation balance, gf grading) sees a different pool, not a shifted one.
3. **The disagreements are concentrated where we are strongest, and that is the argument
   against the pure replication.** Every one of the five largest is a line where *our* value
   is a **primary laboratory** gf and theirs is GES-v3:

   | line | GBS (GES-v3) | ours | our source | Δ |
   |---|---|---|---|---|
   | Fe I 5775.08 | −1.297 | −1.080 | PRIMARY LAB DenHartog2014 | −0.217 |
   | Fe I 4950.11 | −1.670 | −1.500 | PRIMARY LAB DenHartog2014 | −0.170 |
   | Fe I 4907.73 | −1.840 | −1.700 | PRIMARY LAB DenHartog2014 | −0.140 |
   | Fe I 5784.66 | −2.532 | −2.670 | PRIMARY LAB Ruffoni2014 | +0.138 |
   | Fe I 5285.13 | −1.540 | −1.660 | PRIMARY LAB Ruffoni2014 | +0.120 |

   Taking "their gf" for the purest replication means **replacing primary laboratory
   measurements with the 2014-vintage GES-v3 compilation on exactly the lines where our
   pool has the best atomic data.** The GBS − GES column shows the same shape (max 0.180,
   also DenHartog lines), which says the drift is largely GES's own v3 → v5 adoption of
   post-2014 laboratory work, not a disagreement between us and Jofré.

The columns are both present so either can be selected; nothing downstream presumes the
answer. Recommended framing, not a decision: run the replication on **their line
selection** (which is the part the GBS method actually is) and our gf, and keep the
their-gf run as the control that isolates the atomic-data contribution.

## Coverage — the 142 selected lines against every VIS holding

| instrument | holding | span (Å) | span source | in span | telluric-excluded | reachable |
|---|---|---|---|---|---|---|
| harps | `solar_harps` | 3782.6 – 6910.0 | `HoldingSpec.span_A` | 142 | 0 | **142** |
| harps | `solar_harps_molecfit_corrected` | 3782.6 – 6910.0 | `HoldingSpec.span_A` | 142 | 0 | **142** |
| kpno | `solar_kpno` | 2960 – 13000 | RYA-708 registry (VERIFIED) | 142 | 0 | **142** |
| kpno | `solar_kpno_molecfit_corrected` | 2960 – 13000 | RYA-708 registry (VERIFIED) | 142 | 0 | **142** |
| kpno | `solar_kpno_kurucz2005_corrected` | 2990 – 10010 | `HoldingSpec.span_A` | 142 | 0 | **142** |
| iag | `solar_iag` (Baker+2020) | 5001.1 – 11083.5 | `HoldingSpec.span_A` | 132 | 0 | **132** |
| iag | `solar_iag_reiners2016` | 4047.5 – 5001.1 | `HoldingSpec.span_A` | 10 | 0 | **10** |
| crires_plus | all four | ≥ 9800 | `HoldingSpec.span_A` | 0 | 0 | 0 |

* **No line is telluric-excluded anywhere.** Not an assumption: the bluest registered
  telluric complex is O₂ B at 6867 Å and the reddest GBS line is 6820.37 Å.
* The two IAG holdings are complementary, not redundant — 132 + 10 = 142 exactly, which is
  the RYA-767 span declaration doing its job.
* `solar_kpno` and `solar_kpno_molecfit_corrected` declare `span_A=None` because their
  reader inventories its own segments, so `covers()` answers `True` for any wavelength.
  That is the *absence* of a coverage claim, not a coverage claim, so the report uses the
  VERIFIED instrument span from `data/catalog/solar_reference_holdings_rya708.csv` and
  names which source it used in `span_source`.
* `solar_vesta_crires_plus_idp` has no declared span and no registry row; it is reported as
  `UNDECLARED`, never as covered.

## Overlap with our pool (spec item 5)

All 159 lines join to `canonical_gf` (λ+EP, 0 unresolved, 0 ambiguous), so "overlap" is a
question about gf QUALITY, not about presence. Our graded pool is `gf_tier == "LAB"` —
the same key `derive_band_products` uses for `--lines-tier graded`.

| | Fe I | Fe II |
|---|---|---|
| our VIS `canonical_gf` rows | 10 127 | 8 427 |
| our VIS **graded** (`gf_tier=LAB`) pool | 240 | 10 |
| GBS selected set (142) | 133 | 9 |
| **GBS ∩ our graded pool** | **16** | **0** |

Tier breakdown of the 142: 76 OTHER, 41 NIST-C+, **16 LAB**, 9 VALD3.

🔴 **The Fe II overlap is zero for a structural reason, not by chance.** All ten of our
graded Fe II VIS lines are DenHartog+2019 and they span 4173.45 – 4583.83 Å. The GBS list
*begins* at 4787.83 Å. **The two sets are disjoint in wavelength** — a GBS replication
cannot draw on a single one of our graded Fe II oscillator strengths, and no line-by-line
choice would change that. Any Fe II leg of this replication runs on VALD3-scale and
NIST-C+ values, three of which (5991.38, 6084.11, 6456.38) are bare VALD3 and two of which
(6149.26, 6247.56) are the RYA-853 pair carrying a `NOT NIST` provenance note.

🔴 **On the 16-line Fe I overlap, the gf choice is not a rounding matter.** These are
exactly the lines where our value is a primary laboratory measurement, and they are where
Jofré's GES-v3 value disagrees most: mean −0.031, s.d. 0.098, and **5 of the 15 comparable
lines differ by more than 0.1 dex** (max 0.217 at Fe I 5775.08). Fe I 6232.64 is in our
graded pool and has no published GBS gf at all.

## 🔴 What the published record does not contain

Recorded so nobody re-runs these searches.

1. **The gf reference codes cannot be decoded from the copy we hold.** Tables 4/5 cite gf
   sources by integer code (102, 114, 129, 156, 158, 166, 167, 186, 187, 190). The decoder
   is those tables' own "References:" footnote, and in `1309.1099v2.pdf` it failed to
   typeset — it reads literally `References: 102: ????????. 114: ??. 129: ??. …`. The codes
   are transcribed anyway (they are real published values and group the lines correctly) but
   `loggf_ref_gbs_resolved` is `False` on every row. arXiv has only v1 and v2 (checked);
   `www.aanda.org` returns HTTP 403 to a scripted GET. Resolving them needs the publisher
   PDF, fetched by a human.
2. **Tables 4/5 do not cover every golden line.** `table6.dat` flags 193 lines golden; the
   paper tables list 183. The 10 with no row — Fe I 4939.7, 5083.3, 5166.3, 5506.8;
   Fe II 5256.9, 5316.6, 5316.8, 6113.3, 6149.3, 6369.5 — have no published gf here.
   Verified absent by searching the whole extracted PDF text for each wavelength, not by
   the table parse merely failing to find them.
3. **One of those 10 is in the solar selected set: Fe II 6149.26** (golden, `rew_class=pass`,
   EP 3.89). It is golden, it is used for the Sun, and the paper publishes no gf for it.
   Our own row for it is not strong either — `canonical_gf` carries `-2.724` under
   *"VALD3-scale value via nist_reference.csv; NOT NIST (ASD gives -2.854 acc E) — RYA-853"*.
   Both sides are weak on this line, and the file says so on both.
4. **Their line list is GES-v3** (Heiter et al. 2014, in prep., §3), where our
   `canonical_gf` GES seed is v5 (Heiter et al. 2021). The two are not the same
   compilation; the GBS − GES column above measures the difference rather than assuming it
   away.
5. **21 of the 142 selected lines have no published GBS gf at all** — 20 non-golden plus
   Fe II 6149.26. A pure their-gf replication cannot run on those lines; it runs on 121.
