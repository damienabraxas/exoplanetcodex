# RYA-517 steps 2–4 — dual-machine anchor re-validation (py3.12 + numpy 2.2)

**SCIENCE EVENT, not an infra bump.** The banked solar anchor (v1) ran on
Python 3.9.6 + numpy 1.26.4 — *below* iSpec's own declared minimum (numpy ≥ 2.2.5).
v2 re-runs solar FULL on Python **3.12.13 + numpy 2.2.6** (scipy 1.18.0, astropy 8.0.0,
pandas 3.0.3) — the FIRST anchor computed on a stack that meets iSpec's stated floor.
No re-bank performed (Ryan's call); this is the diff table only.

## Environment
- venv: `.venv312` (Homebrew python3.12.13). numpy pinned 2.2.x → 2.2.6; scipy/astropy/pandas
  at latest ≥ iSpec minimums.
- iSpec: the shipped `synthesizer.cpython-312-darwin.so` (v2023.08.04-20-g2126994) was
  installed into the package dir as `ispec/ispec/synthesizer.cpython-312-darwin.so`
  (ABI-tagged; python 3.12 prefers it over the untagged py3.9 `synthesizer.so`, so both
  interpreters coexist with no clobbering). After this, iSpec reports
  **synthesis=True, MOOG=True, Turbospectrum=True** under 3.12/numpy 2.2.
- RYA-514 force-fork + single-thread BLAS carry unchanged: start_method=fork on both stacks.

## Per-element delta — v1 (py3.9/numpy1.26) vs v2 (py3.12/numpy2.2)

| element | ion | A_X v1 | A_X v2 | ΔLTE | A_nlte v1 | A_nlte v2 | ΔNLTE | n v1 | n v2 |
|---|---|---|---|---|---|---|---|---|---|
| Al | I | 7.406 | 7.406 | 0.000 | 7.406 | 7.406 | 0.000 | 1 | 1 |
| C  | I | 10.260 | 10.260 | 0.000 | 10.260 | 10.260 | 0.000 | 1 | 1 |
| Ca | I | 6.308 | 6.308 | 0.000 | 6.324 | 6.324 | 0.000 | 2 | 2 |
| Cr | I | 5.950 | 5.950 | 0.000 | 6.022 | 6.022 | 0.000 | 7 | 7 |
| **Fe** | **I** | **7.510** | **7.510** | **0.000** | **7.516** | **7.516** | **0.000** | **62** | **62** |
| Fe | II | 7.657 | 7.657 | 0.000 | 7.657 | 7.657 | 0.000 | 3 | 3 |
| Li | I | 0.727 | 0.727 | 0.000 | 0.727 | 0.727 | 0.000 | 1 | 1 |
| Na | I | 6.370 | 6.370 | 0.000 | 6.264 | 6.264 | 0.000 | 2 | 2 |
| Ni | I | 6.946 | 6.946 | 0.000 | 6.946 | 6.946 | 0.000 | 2 | 2 |
| S  | I | 7.753 | 7.753 | 0.000 | 7.753 | 7.753 | 0.000 | 2 | 2 |
| Si | I | 7.892 | 7.892 | 0.000 | 7.888 | 7.888 | 0.000 | 7 | 7 |
| Sr | I | 4.961 | 4.961 | 0.000 | 4.961 | 4.961 | 0.000 | 1 | 1 |
| Ti | I | 5.364 | 5.364 | 0.000 | 5.471 | 5.471 | 0.000 | 10 | 10 |

**max |Δ| across all species (LTE & NLTE) = 0.000 dex. n_lines identical on every species.**

## Fe gate + verdict
- **Fe I anchor**: A(Fe I)_NLTE = 7.516 (n=62) — unchanged.
- **Fe scatter gate (RYA-446, profile G fe1_scatter_max=0.1398)**: Fe I std = 0.1390 → **PASS** (same as v1).
- **Solar Fe VERDICT (scale-robust)**: **PASS** [slope ✓ · ionization ✓ · scatter ✓].
- **RYA-371 27-element verdict** (regenerated under v2, classifier over the 27-element baseline
  vs Asplund 2021): **PASS=5 NLTE-OWED=1 CURATION-OWED=20 DATA-GAP=0** — per-element verdict map
  byte-identical to the committed 3.9 verdict (26 element keys, zero diffs).

## Honest finding — byte-level intermediates DIFFER, science outputs do NOT
The numpy 2.x FP path perturbs the byte-level intermediates:
- solar_normalized.csv md5: v1(3.9) `df4a49cf` → v2(3.12) `d2b92e9b`
- solar_ew.csv md5:         v1(3.9) `230c27c4` → v2(3.12) `bcfa16c6`

But the abundances are identical to reporting precision (0.000 dex on all 13 measured
species, incl. the iSpec-synthesis-arbitrated Fe II theo-EW leg). The numpy 2.x differences
are sub-millidex noise that washes out in the median aggregation + 3-dp rounding.

**Consequence for RYA-511 (dual-machine):** cross-machine agreement must be gated on
*abundance*-identity (Δ ≤ tol), NOT byte-identical intermediates — byte-identity will not
hold across numpy builds/platforms even when the science is invariant.

## Bottom line (Mac leg)
The solar anchor is **stack-robust**: moving 3.9.6/numpy1.26 → 3.12.13/numpy2.2 (onto and
above iSpec's declared floor) reproduces v1 exactly — every species 0.000 dex, Fe gate PASS,
27-element verdict identical. No scientific reason not to adopt py3.12; **no re-bank needed
or performed** (v2 ≡ v1).

---

# Steps 3–4 — Sirius (Linux x86_64) leg + cross-check

Ryan: "reference stack = py3.12 + numpy 2.2.x, migrate BOTH machines." Sirius is the
authoritative dual-machine leg (RYA-511). Unlike the Mac, Sirius ships **no** prebuilt iSpec
`.so` — the C extension (and MOOG) had to be **compiled from source on Linux**, which is the
single biggest RYA-517 unknown.

## Sirius environment built this session
- **Python 3.12.13** via pyenv (exact match to the Mac).
- venv `/mnt/codex-data/venv312`: numpy **2.2.6**, scipy **1.18.0**, astropy **8.0.0**,
  pandas **3.0.3**, Cython 3.2.8 — **identical resolved stack to the Mac venv312**.
- **iSpec C synthesizer COMPILED on Linux/gcc-15** → `ispec/synthesizer.so` (1.3 MB), imports
  clean, reports **synthesis=True** under py3.12/numpy2.2. **This resolves the probe's binding
  unknown: iSpec's C extension builds on Linux py3.12+numpy2.x.**
- **MOOGSILENT COMPILED** (gfortran 15.2.0) → EW→abundance baseline engine available
  (RYA-289 requires MOOG; it refuses to downgrade to SPECTRUM).

### gcc-15 build gotcha (documented for the migration)
The vendored SPECTRUM C code uses K&R empty-paren prototypes (`double gffactor();` then
`gffactor(m)`). gcc 15 defaults to **C23**, where `()` means *zero* args → hard error
("too many arguments to function"). Fix: build with **`CFLAGS="-std=gnu17 -fcommon -w"`**
(C17 keeps unspecified-arg semantics AND supports the C99 for-loop declarations that Cython
emits — `-std=gnu89` breaks the latter). The standalone `synthesizer/spectrum/spectrum` and
Turbospectrum binaries were NOT built (own Makefiles, same C23 issue) and are not needed by
the solar EW→abundance path (uses the `.so` + MOOGSILENT).

## Sirius per-element delta — v2 (linux) vs banked v1
After staging the same committed inputs as the Mac, **every measured species reproduces v1 to
0.000 dex** (LTE & NLTE), n_lines identical:

| element | Fe I | Fe II | Al | C | Ca | Cr | Li | Na | Ni | S | Si | Sr | Ti |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ΔNLTE (sirius−v1) | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |

- **Fe I anchor**: 7.516 (n=62) — unchanged. **Fe scatter gate**: std 0.1390 ≤ 0.1398 → **PASS**.
- **RYA-371 27-element verdict** regenerated on Sirius: **PASS=5 NLTE-OWED=1 CURATION-OWED=20
  DATA-GAP=0** — identical to Mac v2 and v1.
- **max |Δ(sirius v2 − v1)| = 0.000 dex.**

## Cross-check result
**Mac v2 (darwin arm64) ≡ Sirius v2 (linux x86_64) ≡ banked v1** — both platforms, both on
py3.12/numpy2.2 (above iSpec's floor), reproduce the anchor to reporting precision. The anchor
is platform- AND stack-robust.

## DUAL-MACHINE GOTCHA (critical for RYA-511)
`data/processed/solar_ew_ges_reference.csv` is a **git-tracked committed input** (the 62-line
GES Fe I pool), NOT a generated file. When it was absent on Sirius, the run **silently fell
back** to the canonical `sol_ew_results_v1.csv` for Fe I → a *different* 19-line pool →
A(Fe I) +0.12 dex, scatter 0.544, **gate FAIL**, and `aberr=nan` on many lines. Staging the
committed file restored exact reproduction. Lesson: cross-machine parity requires the full
committed input set (incl. tracked files under `data/processed/`), and a missing Fe I GES
reference degrades *silently* rather than erroring — a fallback worth making loud (follow-up).

## True numerical floor (full precision, split by axis)
The reported A(X) columns are `round(median, 3)`, and MOOG itself prints per-line abundances at
**1e-3 dex granularity** — that quantum is the finest quantity the pipeline produces, so it *is*
the floor; there is no sub-mdex float exposed to compare. Capturing the UNROUNDED per-line
`normal_abund` (iSpec `determine_abundances` return) on each stack and comparing all 473 lines
at MOOG's native granularity:

| split | axis isolated | max \|Δ\| | lines differing |
|---|---|---|---|
| **same-machine** (Mac): v1 py3.9/np1.26 ↔ v2 py3.12/np2.2 | numpy/stack | **0 dex (exact)** | 0 / 473 |
| **cross-machine** (both v2): Mac arm64 ↔ Sirius x86_64 | platform | **0 dex (exact)** | 0 / 473 |
| (both-axes) Mac py3.9 ↔ Sirius py3.12 | stack+platform | 0 dex (exact) | 0 / 473 |

Exact identity to MOOG's last printed digit on every line — **not** rounding-masked. The
aggregate A(X) inherit this (identical inputs → identical median → identical `round(·,3)`),
hence max|Δ| = 0 exactly for both splits.

## Byte-level intermediates: three distinct hashes, one science answer
| stage | v1 (mac py3.9) | v2 (mac py3.12) | v2 (sirius py3.12) |
|---|---|---|---|
| solar_normalized.csv | df4a49cf | d2b92e9b | 833554ad |
| solar_ew.csv | 230c27c4 | bcfa16c6 | 4e8e40e6 |

All three differ (numpy-2.x FP + platform libm/BLAS) — proving sub-quantum FP noise DOES exist
upstream — yet it is far below the 1e-3 dex MOOG quantum and never reaches even the finest
reported abundance digit (per-line max|Δ| = 0, above). Confirms: **gate cross-machine agreement
on abundance-identity, never on byte-identical intermediates.**
