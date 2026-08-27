# RYA-1075 — the HFS-per-isotope inflation in `canonical_gf`, corrected

**54 rows corrected. No published value moves. The guard RYA-684 lacked is installed.**

## What the defect actually is

`gf_resolver.cluster_physical_lines` groups components into physical lines on species + excitation potential + wavelength gap. It is **isotope-blind by construction** — and two isotopes of the same transition share the lower level and sit milli-Ångströms apart, so they land in one cluster. The GES v6 HFS/ISO list is RYA-684 *form (A)*: isotope-coded, with **each isotope's component set carrying the full oscillator strength**, because the engine applies `isotopfrac` afterwards and Σf = 1. So any consumer that writes `log10(Σ 10^gf)` over such a cluster gets `n_isotopes ×` the physical gf.

Demonstrated end to end on Eu II 6645, not inferred from the offset:

```
11 GES components, isotopes 151 (7 of them) and 153 (4)
  isotope 151 sums to   +0.1199      <- the whole transition
  isotope 153 sums to   +0.1198      <- the same total, as form (A) requires
  naive sum over all    +0.4208      == canonical_gf's published log_gf, exactly
  gf-weighted centroid   6645.0905 A == canonical_gf's wavelength_air_A, exactly
  offset                +0.3010      == log10(2)
```

## This is not RYA-684, and RYA-684's correction term would not have fixed it

RYA-684 measured an **engine-side** double application — `isotopfrac` applied on top of already-folded gf — whose offset is `−log10(Σ fᵢ²)` and which lives in the TSFitPy `linelist_vald` files. Its module docstring warns, correctly, not to reuse log10(2) as a correction anywhere. This is a **consumer-side aggregation** over a correctly-formed source, and its offset is `log10(n)` — a pure count, with no dependence on the abundances.

**La II settles which is which.** Parsed from `makeabund.f`:

| species | RYA-684 term −log10(Σfᵢ²) | count term log10(n) | measured |
|---|---:|---:|---:|
| La II | +0.0008 | +0.3010 | **+0.3010** |
| Nd II | +0.7258 | +0.8451 | **+0.8451** |
| Sm II | +0.7362 | +0.8451 | **+0.8451** |
| Eu II | +0.3002 | +0.3010 | **+0.3005** |

La is 99.911 % La-139, so the two hypotheses are 0.30 dex apart there and the data picks the count without ambiguity. Eu alone could never have separated them — which is why RYA-684, whose evidence was Eu, reasonably reached a different term.

## The correction could not ship alone

`apply_to_synth_array` shifts each cluster by `canon_total − cur_total`, and computed `cur_total` with the **same** isotope-blind sum. Measured, for Eu II 6645:

| | old resolver | new resolver |
|---|---:|---:|
| canonical **before** (+0.4208) | shift −0.0000 ✅ | shift **+0.3010** ❌ |
| canonical **after** (+0.1198) | shift **−0.3010** ❌ | shift −0.0000 ✅ |

Only the diagonal is safe. Today both numbers were equally inflated and the error **cancelled** — which is exactly why RYA-684 correctly concluded the GES surface reached no live value. De-inflating the table on its own would have scaled every component down by a factor n and turned a cancelling error into a live one. So `gf_resolver.physical_total` ships in the same change; neither half is correct alone.

## The consumer fix uncovered a second, larger defect on main

Making `apply_to_synth_array` isotope-aware changed the shift for 12 clusters — and those
are all rows this ticket did **not** correct. Measured across every multi-isotope cluster,
`main` vs this branch, as the effective gf the *engine* ends up synthesising:

| species | clusters | corrected here? | effective gf change |
|---|---:|---|---:|
| Eu II | 5 | yes | ≤ 7 × 10⁻⁵ |
| Nd II | 43 | yes | ≤ 2 × 10⁻⁶ |
| Sm II | 4 | yes | ≤ 2 × 10⁻⁶ |
| Ba II / Cu I | 2 | yes | ≤ 7 × 10⁻⁵ |
| **Nd II / Sm II** | **4** | **no — already correct** | **+0.8451** |
| **Ba II** | **3** | **no** | **+0.6990** |
| **Cu I** | **5** | **no** | **+0.3010** |
| Li I | 1 | no | 0.0000 |

On `main`, `apply_to_synth_array` normalised the isotope-**blind** cluster sum to the
canonical total. Turbospectrum then applies `isotopfrac` per isotope and synthesises the
**per-isotope** total, which is log10(n) below the number that was normalised. So the gf
that actually reached the engine was `log10(n)` **too weak** for every cluster whose
canonical value was a genuine physical total.

The irony is exact: **the rows whose canonical value was already correct are the ones that
were being mis-synthesised.** The 54 inflated rows were safe, because their inflation
cancelled the resolver's.

**One committed value is exposed.** `scripts/measure_cu_v_hfs_synthesis_rya466.py` calls
`_load_synth_resources()` with the default `apply_canonical_gf=True`, so its five Cu I
lines were synthesised 0.3010 dex too weak — a re-run would put **A(Cu) about 0.30 dex
lower**. It is **not** re-derived here: RYA-161 is report-before-re-deriving and this
ticket tunes no abundance. **A(Ba) 2.237 is not exposed** — RYA-581 builds its own VALD
blocks and never calls this path (as RYA-684 recorded) — and Nd/Sm/Eu carry no committed
abundance.

## What was corrected

| species | rows | isotopes coded | correction |
|---|---:|---:|---:|
| Ba II | 1 | 5 | -0.6990 dex |
| Cu I | 1 | 2 | -0.3010 dex |
| Eu II | 5 | 2 | -0.3010 dex |
| Nd II | 43 | 7 | -0.8451 dex |
| Sm II | 4 | 7 | -0.8451 dex |

**54 rows.** Full per-row provenance — original value, correction term, n_isotopes, per-isotope reconstruction, prior adjudication status — is in `data/linelists/canonical_gf_isotope_corrections_rya1075.csv`.

The write was surgical: of 178,680 × 25 cells, exactly 54 `log_gf` and the matching 54 `adjudication_status` changed. Nothing else in the table moved.

## What was NOT corrected, and why that matters

The selector is specific, not a species rule. 13 multi-isotope clusters were rejected:

* **`NOT_BUILT_BY_NAIVE_SUM` (12)** — the published value is not the naive cluster sum, so the row was not built by it. Four (3 Nd II, 1 Sm II) already carry the CORRECT value; the rest (5 Cu I, 3 Ba II) differ for unrelated catalogue reasons.
* **`NOT_FORM_A` (1)** — the isotopes do not each carry the full gf, so a count correction is undefined. This is **Li I 6707**, the positive control: a genuine multi-isotope cluster that must survive untouched, and does.

Measured separation: `|offset − log10(n)|` is at most **1.8 × 10⁻⁵** among corrected rows and at least **2.7 × 10⁻²** among rejected ones — a factor of ~1,500. The threshold sits nowhere near a boundary.

**La II 5971 is not this defect.** Its offset is also +log10(2), but its GES cluster is single-isotope: the source contains the row **twice**, byte-identical (same λ, gf, EP, reference, `iso=0`). Correcting it under an isotope rationale would have attached the wrong provenance to the right number. It is left alone and reported — one of 17 clusters in the GES v6 delivery containing duplicated rows.

## Corroboration

49 of the 54 corrected rows have a `gf_linelist_vald` sibling. Against it the corrected values agree to a **median 1.4 × 10⁻⁵ dex**, 95th percentile 7 × 10⁻⁴. Two exceptions, both reported rather than smoothed:

* **Cu I 5782.122** — corrected −1.7889 against a sibling of −1.9048, a **0.116 dex** gap. The inflation itself is proven from the source (offset exactly log10(2), form (A) confirmed); the residual is a GES-vs-VALD catalogue disagreement. RYA-684 already established Cu I's VALD surface is *folded*, so the sibling is not a trustworthy referee for this species specifically. Corrected on the source proof, flagged for adjudication.

* **5 Nd II rows have no sibling at all** and were corrected on the source reconstruction alone. Named individually in the classification output.

## The defect had propagated into an adjudication record

`RYA-354` stamped Eu II 6645 as *"STAMP graded (LWHS) — best-available cited, value frozen"* at **+0.4208**, citing Lawler, Wickliffe, den Hartog & Sneden 2001. The citation is right; **+0.4208 was never that paper's number**. The freeze worked exactly as designed — it froze an inflated input. Both the adjudication CSV and its pinned test are updated, with the reason recorded rather than the number quietly swapped.

`data/results/eu2_synthesis_rya565.json` had already recorded `canonical_loggf == ges_naive_sum_loggf == 0.4208` alongside `physical_total_loggf: 0.1198` and `naive_minus_physical_dex: 0.301`. **The fingerprint was written down in August and read as a checksum rather than a defect.** That artifact is a historical measurement record and is left as it stands.

## RYA-102 / A(Eu): the required validation

**A(Eu) does not move, and this is measured rather than argued.**

`scripts/rya565_eu2_synth_sirius.py` synthesised on `physical_loggf` = **+0.1198** — the per-isotope total — not on the canonical +0.4208. Its own result artifact records that. So the only path that has ever produced an A(Eu) already used the corrected number:

| | before | after |
|---|---:|---:|
| A(Eu) HARPS (RYA-565, LTE) | 0.702 | 0.702 |
| A(Eu) IAG (RYA-565, LTE) | 0.603 | 0.603 |
| effective gf reaching Turbospectrum | −0.9489…+0.4208 delivered, shift 0.0000 | identical |

**Effect on RYA-102: none, and RYA-102 is already closed.** It was resolved on 2026-06-29 by a fit-**window** change (`LINE_WINDOWS` 0.30 → 0.15 Å, EW 5.51 → 3.33 mÅ), and its 6–10 mÅ acceptance criterion was found to be 55-Cnc-calibrated rather than solar (Lawler+2001 solar is ~1.6–3 mÅ). A gf error cannot produce that symptom in any case: gf does not enter an EW **measurement**, only the abundance inversion afterwards — and there the correct value was already in use. The RYA-1075 ticket describes RYA-102 as open; it is not.

## Keyed on the stable id, not the row index

RYA-1077 landed on `main` mid-ticket and it applies directly here: `canonical_gf.line_id`
is `gf_NNNNNN`, assigned by **row position**, and it rots whenever a block of rows is
replaced — 1,739 committed references (25%) had already moved, some to a different
*species*.

The correction sidecar is exactly such a reference, so it is keyed on RYA-1077's
`physical_id` and carries `line_id` only as a human convenience. The apply, the verify and
the regression tests all resolve rows by `physical_id`; the guard itself never needed an id
at all, because it re-derives from the source. `--verify` checks the **table against the
committed sidecar**, not against a fresh classification — after a successful apply the
classifier correctly finds nothing, and "0/0 verified" would be a vacuous pass.

## The guard

`scripts/check_stewardship.py::check_isotope_inflation` (invariant 10), backed by `pipeline.isotope_gf_convention.isotope_inflated_rows`. It **re-detects from the source** rather than pinning the 54 corrected line ids — a pinned list passes forever while a new ingest reintroduces the defect on different lines, which is precisely how RYA-684 came to be closed with 54 live instances still in the table. It is registered **untracked**: there is no remediation ticket because the correct state is zero, so a hit is a real break.

```
corrected table          -> 0 violations
one row re-inflated      -> 1, named (gf_051798 Eu II 6645.0905, +0.3010 dex)
all 54 re-inflated       -> 54
Li I 6707 positive control -> never flagged, in any of the above
```

