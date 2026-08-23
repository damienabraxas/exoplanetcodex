# RYA-1001 — Al PHASE 0: full-band line census (UV / VIS / IR)

**Census only.** No synthesis was run, no abundance derived, no gf adopted. Everything
below is provenance and inventory. Executed on Sirius, worktree
`/mnt/codex-data/codex/rya1001`, off `origin/main` `4c031dd`.

Concurrency: the RYA-999 cap of 2 was respected — RYA-986's two Fe arms and a CI
`pytest` job were already on the box, so every step here ran `nice -n 10`,
single-threaded, and is I/O + catalog work.

---

## THE HEADLINE

The ticket's frame — *"Al is Fe's mirror: VIS empty-by-firewall, Al lives in the IR via
7835/8772 synthesis"* — is **half right and points at the wrong IR lines.**

Al does live in the IR. But **7835/7836/8772/8773 are not Al's best-graded lines and
never were.** Burheim 2023 publishes experimental log gf for **eight Al I lines inside
our wavelength reach**, and none of them is the 7835/8772 quartet. Two are the VIS pair
6696.02 / 6698.67 that the audit believed were firewalled; four more are in the NIR; and
the two best of all — **13123.42 Å at 1.5 % and 13150.75 Å at 3.1 %** — are deep solar
lines (central depth 0.474 / 0.444) that **we already hold real pixels for**, in the
CRIRES+ Vesta IDPs.

RYA-835 concluded `burheim_covers_any_target_line = false` and RYA-716 concluded "do not
queue Al for CRIRES+". Both conclusions rest on defects this census identified and both
should be reversed.

---

## OQ1 — Does Burheim 2023 cover 7835 / 8772?

**Answer: NO for those four lines — CONFIRMED. But the RYA-835 evidence for it was
invalid, and the same defect produced a false negative on six other lines.**

### The verdict, from the actual table

Pulled the arXiv e-print source for 2309.06273 (`45394corr.tex`) and read Table 3
(`tab:loggf_comp`), which is the table of derived log gf values. It has **exactly 12
rows**, spanning λ_vac 6697.864 → 41920.7 Å — i.e. 670–4192 nm, which reproduces the
paper's stated "12 lines, 670–4200 nm, 2–11 %" exactly. That is the control that the
right table was read.

None of the 12 is 7835 / 7836 / 8772 / 8773. Burheim's upper levels are 4p, 5p, 4f, 5s,
6s, 7s and 4d; the 7835/8772 quartet is 3d ²D – 5f ²F°, which he did not study. **The
absence is real, and it is proven with a positive control**: the same matcher, run over
the same 12 lines, produced 8 hits elsewhere at ≤ 0.005 Å. A matcher that finds nothing
anywhere proves nothing; this one discriminates.

### The defect behind the old answer

`scripts/rya835_al_gf_resolution.py` carries a 35-entry `BURHEIM_LINES` tuple. **It is
the union of two different columns.** Burheim's tables print `σ [cm⁻¹]` next to
`λ_vac [Å]`, and the tuple contains both — wavenumbers and wavelengths mixed into one
list of "wavelengths". Every "nearest Burheim line" computed from it is meaningless:

| RYA-835 claim | what that number actually is |
| --- | --- |
| `burheim_nearest_A = 7836.521`, 0.387 Å from our 7836.134 | a **wavenumber**, 7836.521 cm⁻¹ = λ_vac 12760.765 Å |
| `burheim_nearest_A = 6697.864`, "1.679 Å" from our 6696.185 | a **vacuum wavelength** — = 6696.02 Å in **air**, i.e. a direct hit that was read as a miss |

The second row is the costly one: it turned a covered line into an uncovered one.

### What Burheim actually covers, in air, in our linelist

| λ_air | λ_vac | transition | log gf | unc | central depth | band |
| --- | --- | --- | --- | --- | --- | --- |
| 6696.015 | 6697.864 | 4s ²S₁⁄₂ – 5p ²P₃⁄₂ | −1.46 | 9 % | 0.253 | VIS |
| 6698.673 | 6700.522 | 4s ²S₁⁄₂ – 5p ²P₁⁄₂ | −1.76 | 8 % | 0.167 | VIS |
| 11253.189 | 11256.270 | 3d ²D₃⁄₂ – 4f ²F₅⁄₂ | +0.167 | 5 % | 0.342 | NIR |
| 11254.924 | 11258.008 | 3d ²D₅⁄₂ – 4f ²F₇⁄₂ | +0.327 | 2 % | 0.361 | NIR |
| 12749.909 | 12753.397 | 3d ²D₅⁄₂ – 5p ²P₃⁄₂ | −2.29 | 9 % | 0.015 | NIR |
| 12757.275 | 12760.765 | 3d ²D₃⁄₂ – 5p ²P₁⁄₂ | −2.62 | 11 % | 0.008 | NIR |
| **13123.416** | 13127.005 | 4s ²S₁⁄₂ – 4p ²P₃⁄₂ | **+0.232** | **1.5 %** | **0.474** | NIR |
| **13150.753** | 13154.350 | 4s ²S₁⁄₂ – 4p ²P₁⁄₂ | **−0.098** | **3.1 %** | **0.444** | NIR |

The remaining four (38632 / 38721 / 41841 / 41920 Å) are mid-IR, outside every instrument
we hold or could hold.

### A second RYA-835 claim to retire

Its docstring says Burheim "DOES cover 10875.953 / 13127.005 / 13154.350 / 16723.541 /
16767.948 / 21098.84 / 21169.58". Only **13127.005 and 13154.350** are in Table 3. The
rest appear only in Table 2 (`tab:BRtable`), which tabulates **branching fractions with
no derived log gf**. Burheim does not grade them. In particular he does **not** grade
10891.7 Å (1.089 µm) — the line the litscan flagged as the telluric-excluded novelty
candidate.

### 🔴 The finding nobody asked for: the lab values are already in our linelist

`linelist_solar.csv`'s HFS components, summed, reproduce Burheim's experimental log gf
on **every line he measured in our range**:

| λ_air | Σ HFS components | Burheim | Δ |
| --- | --- | --- | --- |
| 6696.015 | −1.4598 | −1.46 | +0.0002 |
| 6698.673 | −1.7597 | −1.76 | +0.0003 |
| 11253.189 | +0.1670 | +0.167 | 0.0000 |
| 12749.909 | −2.2900 | −2.29 | 0.0000 |
| 12757.275 | −2.6198 | −2.62 | +0.0002 |
| 13123.416 | +0.2322 | +0.232 | +0.0002 |
| 13150.753 | −0.0977 | −0.098 | +0.0003 |

Seven exact agreements. **Negative control:** on the four lines Burheim did *not*
measure, our linelist carries K75/OP95 values that sit 0.03–0.15 dex off NIST — so the
VALD3 pull is not uniformly "the best number", it is Burheim's number exactly where
Burheim has one. The identification is not a coincidence.

**And `canonical_gf.csv` — the declared single source — is the stale one.** It carries
−1.347 / −1.647 for the VIS pair, which are *exactly* Burheim's `K95` column, i.e.
Kurucz 1995. So on Al's two VIS lines the adjudicated single source is **0.113 dex away
from the experiment, in the direction of a semi-empirical value the project's own cull
rejects when it is labelled honestly.**

---

## OQ2 — 8772/8773 and 7835/7836 provenance, line by line

**Answer: the "1995 J.Phys.B (lab)" attribution is REFUTED. It is Opacity Project
theory. 7835/7836 = Kurucz 1975 is CONFIRMED.**

The GES v6 atomic list (`GESv6_atom_hfs_iso.420_920nm/atomic_lines.tsv`) is the seed for
both, and Al I there carries exactly three reference codes: `K75` (76 lines),
`1995JPhB..` (5 lines), `WSM` (3 lines).

| line | GES `reference_code` | resolves to |
| --- | --- | --- |
| 7835.309, 7836.134 | `K75` | Kurucz 1975 semi-empirical — **CONFIRMED** |
| 8772.865, 8773.896 | `1995JPhB..` (truncated) | **Mendoza, Eissner, Le Dourneuf & Zeippen 1995, J. Phys. B **28**, 3485** — "Atomic data for opacity calculations. XXIII. The aluminium isoelectronic sequence". DOI `10.1088/0953-4075/28/16/006`, Crossref-verified. An **ab initio close-coupling Opacity Project calculation in LS-coupling. THEORY, not a lab measurement.** |

**How it was identified rather than guessed.** The truncated bibcode alone proves
nothing. Burheim's Table 3 independently transcribes the TOPbase/OP values as its `M00`
column, and the GES `1995JPhB..` values match them exactly on both lines where an
overlap exists:

| line | GES `1995JPhB..` | Burheim `M00` (TOPbase) |
| --- | --- | --- |
| 6696.023 | −1.569 | −1.569 |
| 6698.673 | −1.870 | −1.870 |

Two exact three-decimal matches, from a source that never saw our repo. Mendoza 1995 is
the OP aluminium calculation that TOPbase serves, so `1995JPhB..` ≡ TOPbase ≡ theory.

**This closes RYA-838 spec item 1 outright**, and it dissolves the puzzle 838 was filed
for: the "0.179 dex disagreement between NIST and the incumbent, 6× NIST's B+ accuracy"
is not a lab-versus-lab conflict at all. It is a graded compilation (NIST B+, 3 %)
against an *ab initio* opacity calculation with no published per-line bar. That is an
ordinary theory-minus-experiment offset, not an anomaly.

**It also retires a headline from RYA-716.** That comment's §4 concluded *"the three
lines that agree all carry the 1995 J.Phys.B **experimental** oscillator strengths"* and
built the 7836 outlier disposition on the experimental/semi-empirical split. There was
no experimental leg: the pool was three OP-theory lines and one Kurucz line.

---

## OQ3 — Does 6696 recover for free when the cull is re-run?

**Answer: NO — it was never culled. And the one line the cull keeps is graded off the
wrong transition, so the honest VIS graded pool is ZERO, not one.**

Re-ran `pipeline.curate_nonfe_pools --phase 2 --grade-restrict` on current code against
current `canonical_gf.csv`:

| line | EW | gf_tier | RYA-395 cull | RYA-398 grade-restricted |
| --- | --- | --- | --- | --- |
| 6631.218 | 34.53 mÅ | LOW (`K75`) | kept | **culled — GRADE** |
| 6696.185 | 59.65 mÅ | MED (NIST C+) | kept | kept |

`recovery_count = 0`, exactly reproducing RYA-835. No quality cut (WEAK / SAT / HIERR /
BLEND / BADGF) fires on either line.

### 🔴 But the surviving line's grade belongs to a different transition

`canonical_gf` row 6696.185 carries `nist_grade = C+`, `log_gf = −1.569`, and the
citation `NIST ASD v5.11; Kelleher & Podobedova 2008 … L4737`. Checked against the cached
NIST pull (`nist_asd_AlI_6600_8800.tsv`, 31 rows):

* **NIST ASD has no Al I line at 6696.185 anywhere in 6600–8800 Å.**
* The row those values came from is NIST **6696.015**, whose lower level is
  **3.1427 eV (4s ²S)**. Our 6696.185 row has **EP 4.0215 eV (3d ²D)**. A 0.17 Å
  wavelength gap and a **0.88 eV different lower level** — not the same transition.
* GES's own value for 6696.185 is `K75 −1.576`, i.e. **LOW tier**.

Under a match that tests the lower level and not just the wavelength, **6696.185 is
Kurucz semi-empirical → LOW → culled**, and the Al VIS EW pool contains **zero**
independently graded lines. The tracker's `A(Al) = 7.406` for this line — +0.98 dex above
the 6.43 anchor — is derived on that borrowed gf.

RYA-840 already flagged this line-ID discrepancy independently. This census confirms it
from the NIST side and adds the consequence: it is the sole survivor of the graded cull.

### The real reason the VIS pool is one line

Not the cull. **The measured EW pool has only two Al lines in it** — and neither of the
two Burheim-graded VIS lines (6696.015, 6698.673) is among them. The census finds
**7 GRADEABLE features in the VIS band**; the pool sees 2. The bottleneck is line
selection into the EW pool, not the gf firewall.

### And a laundering path in `_gf_tier`

`pipeline/curate_nonfe_pools._gf_tier` culls `K75` as LOW but tiers `VALD3` as MED,
because the rule is "not a Kurucz tag and not blank". For 6696.023 and 6698.673 the
`VALD3`-labelled values are **−1.347 / −1.647 — Burheim's `K95` column exactly**, i.e.
Kurucz 1995. The same semi-empirical source passes or fails purely on whether the
aggregator's name or the original source's name is in the string.

---

## OQ4 — Is Al's gf work homeless?

**Answer: NO. The premise is out of date — four open tickets already own it. What is
genuinely homeless is the *experimental* rung.**

RYA-161/162 are indeed Done and Al is indeed absent from RYA-697's successor cluster, and
`element_status_tracker.csv` still names `RYA-161/162 differential-gf` as Al's path under
`GENUINELY-OWED (gf)`. That tracker cell is stale. RYA-835 filed successors:

| ticket | state | owns |
| --- | --- | --- |
| **RYA-838** | Backlog | RCA the 8772/8773 gf disagreement; un-truncate `1995JPhB..` — **answered here** |
| **RYA-839** | Backlog | Al IR lab-gf pool rebuild (7835/36, + 8772/73 pending 838) |
| **RYA-840** | Backlog | the 6696 line-ID discrepancy — **corroborated here** |
| **RYA-778** | Backlog | the −0.075 dex 1D-LTE zero-point RCA |
| RYA-946 | Backlog | all-element lab-gf sweep umbrella |

**What no ticket owns:** ingesting **Burheim 2023** as Al's primary lab-gf source, the
way `data/reference/fe_gf_lab/fe1_lab_loggf.csv` serves Fe. There is no Al lab table in
the repo at all, and `pipeline/gf_grades.py`'s `GF-LAB` state is hard-wired to the Fe
file. Every one of 838/839/840 is scoped around *which theory value to prefer*, because
all of them were written believing no experiment covered these lines.

**Recommended home ticket (Claude PM to file):** *"BUILD: ingest Burheim 2023 as Al's
primary lab-gf table + generalise `gf_grades.GF-LAB` off the Fe-only file."* It
supersedes the theory-adjudication scope of 838 and unblocks 839.

Two corrections owed to existing surfaces, flagged not fixed (census only):
* `element_status_tracker.csv` Al row → point at 838/839/840, not the closed 161/162.
* `data/reference/litscan/Al.yaml` still does not exist (RYA-716 item 4, still open).

---

## The census — gradeable pool per band

`data/results/rya1001/rya1001_al_line_census.csv` (502 features collapsed from 1117 linelist rows).

| band | features | Al I / II / III | HFS carried | home | **GRADEABLE** | CAND-BLENDED | shallow | saturated | no home | Burheim-graded |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| FUV < 3000 Å (**no declared policy**) | 324 | 287 / 35 / 2 | 120 | 0 | 0 | 0 | – | – | **324** | 0 |
| near-UV 3000–3800 | 12 | 10 / 2 / 0 | 5 | 12 | **2** | 0 | 2 | 8 | 0 | 0 |
| VIS 3800–6910 | 26 | 24 / 2 / 0 | 13 | 26 | **7** | 4 | 13 | 2 | 0 | **2** |
| red-optical 6910–10000 | 41 | 41 / 0 / 0 | 12 | 41 | **6** | 2 | 33 | 0 | 0 | 0 |
| NIR 10000–24000 | 96 | 96 / 0 / 0 | 37 | 96 | **14** | 2 | 80 | 0 | 0 | **4** |
| > 24000 Å (no declared policy) | 3 | 3 / 0 / 0 | 1 | 2 | 0 | 0 | 2 | 0 | 1 | 0 |

**Gradeable pool = 29 lines** (+ 8 CANDIDATE-BLENDED). Tiering is on observable
properties only, thresholds fixed in advance and inherited from RYA-709
(`central_depth ∈ [0.05, 0.60]`), never from an abundance.

### Best gf source across the gradeable + candidate pool

| band | EXP-BURHEIM23 | NIST-B/B+ | NIST-C/C+/D | THEORY-OP95 | SEMIEMP-K75 | UNGRADED |
| --- | --- | --- | --- | --- | --- | --- |
| near-UV | 0 | 0 | 0 | 0 | 0 | 2 |
| VIS | **2** | 0 | 2 | 1 | 4 | 2 |
| red-optical | 0 | 4 | 4 | 0 | 0 | 0 |
| NIR | **4** | 0 | 0 | 0 | 0 | 12 |

* **The FUV is a hard zero.** 324 features — 287 Al I and 35 Al II, including the
  1670 Å Al II resonance line — and **no solar holding reaches below 2960 Å**. It is also
  a regime `band_policy` declines to declare a method for. This is the Al twin of
  RYA-909's near-UV Fe II result: an acquisition question, not a task.
* **The red-optical is NIST-graded but experiment-free.** All 8 graded lines there rest
  on the Kelleher & Podobedova 2008 compilation. That includes the entire 7835/8772 set.
* **12 of the 16 NIR gradeable lines are UNGRADED** — including the APOGEE H-band set
  (16718 / 16750 / 16763) and 21093 / 21163, which are deep (0.35–0.49) and heavily
  hyperfine-split. Burheim measures branching fractions for several but publishes no
  log gf.

### 🔴 A gf zero-point visible across the whole red-optical

Every red-optical Al I line for which NIST has a graded value sits **0.031–0.048 dex
above** NIST in our linelist, with no scatter to speak of:

| line | linelist | NIST (graded) | Δ |
| --- | --- | --- | --- |
| 7083.969 | −1.1110 | −1.1574 (C+) | −0.0464 |
| 7084.643 | −0.9348 | −0.9805 (C+) | −0.0457 |
| 7361.568 | −0.9030 | −0.9462 (C+) | −0.0432 |
| 7362.296 | −0.7268 | −0.7693 (C+) | −0.0425 |
| 7835.309 | −0.6490 | −0.6887 (B) | −0.0397 |
| 7836.134 | −0.4728 | −0.5131 (B+) | −0.0403 |
| 8772.865 | −0.3160 | −0.3487 (B+) | −0.0327 |
| 8773.896 | −0.1398 | −0.1714 (B+) | −0.0315 |

Eight lines, one sign, a 0.015 dex spread. That is a **scale offset, not scatter** — the
signature RYA-778 was filed to characterise, now measured. **Firewall (RYA-161): the fact
that closing it would move A(Al) upward is explicitly NOT the argument.** The argument is
provenance — a graded compilation against Kurucz-1975 and an OP calculation — and the
substitution is a pool rebuild (RYA-839), not a census action. Recorded, not adopted.

---

## HFS — mandatory, and it is carried exactly where physics says it should be

26 of the 37 gradeable/candidate features carry multiple components. The pattern is not
random and it is not a data gap:

| lower level | HFS in linelist? | example |
| --- | --- | --- |
| 4s ²S (3.143 eV) | **yes** | 6696.015 (6), 6698.673 (4), 13123.416 (6), 13150.753 (4) |
| 4p ²P (4.085/4.087) | **yes** | 10782 (12), 16750 (12), 16763 (10), 21163 (6) |
| 3d ²D (4.0215/4.0216) | **no** | 7835.309 (1), 8772.865 (1) |

²⁷Al is 100 % abundant with I = 5/2, so every level splits in principle — but the Fermi
contact term that dominates Al I hyperfine structure scales with s-electron density at
the nucleus. It is large for 4s, appreciable for 4p, and negligible for 3d. **So the
ticket's premise "the IR doublets have real hyperfine structure; the linelist carries
components" is false for 7835/8772 — and the linelist is right, not deficient.** Those
two are genuinely unsplit. 7836.134 and 8773.896 carry 2 components each, but those are
*fine*-structure J-partners, not hyperfine.

Consequence for Phase 1: HFS must be summed, not maxed — see the next section.

### 🔴 `line_accounting_rya709.features()` mis-handles Al's HFS in two ways

Verified by executing it, not inferred:

1. **It reports `gf = max` over the cluster — the strongest component, not the line.**
   For 6696.02 that is −1.886 where the six components sum to −1.460: a **0.43 dex**
   under-report on the very quantity being audited.
2. **`GROUP_A = 0.05 Å` splits real multiplets.** 13123.416 (internal gap 0.0577 Å)
   becomes two features, 13150.753 becomes two, 16750 (span 0.197 Å) becomes three.

Both are fixed in this census's own collapse, which groups on a **wavenumber** span —
because hyperfine splitting is bounded in energy, and no fixed-Å rule can serve both
1670 Å and 21163 Å. The threshold is measured, not chosen: on our own Al list the largest
gap *inside* a known multiplet is 0.070 cm⁻¹ and the smallest gap *between* distinct
lines is 1.345 cm⁻¹, so 0.30 cm⁻¹ sits in empty space with a factor > 4 on each side.

**Positive control:** the collapse reproduces `canonical_gf.csv`'s independently-recorded
summed `log_gf` on **20/20** overlapping lines to 0.005 dex, and its
`hfs_n_components` on 18/20. The two disagreements are a *canonical_gf* metadata bug —
3944.006 and 3961.520 are recorded as `hfs_n_components = 1` while the linelist splits
them into 4 and 6 components that sum to canonical's own value exactly.

---

## Per-instrument homes, and a silent hole in the coverage module

| instrument | span | host | Al lines it is the *only* cover for |
| --- | --- | --- | --- |
| HARPS | 3782.6–6910.0 | mac | – |
| Kitt Peak atlas | 2960–13000 | **mac** | – |
| IAG FTS | 4047.4–10649.9 | sirius | – |
| Delbouille/Liège | 3000–10000 | sirius | the two near-UV lines (3059.03, 3059.92) |
| **CRIRES+ (Vesta IDP)** | 9479–24855 | sirius | **13123.42, 13150.75 and everything > 13000 Å** |

### 🔴 `pipeline.coverage` cannot see 6 of the 10 solar holdings

`coverage.load_registry()` skips a holding whose `manifest_path` is not a
spectrum-location CSV — and the skip is a bare `continue`, silent:

| holding | why it is invisible |
| --- | --- |
| `solar_vesta_crires_plus_idp` | manifest CSV has no `loader` column |
| `solar_crires_plus_y_rya794` | manifest CSV has no `loader` column |
| `elgueta2026_vizier` | manifest is `MD5SUMS.txt` |
| `solar_harps_molecfit_corrected` | manifest is `.json` |
| `solar_kpno_molecfit_corrected` | manifest is `.md` |
| `solar_kpno_kurucz2005_corrected` | manifest is `.json` |

So the module that exists **specifically** to stop the project reporting "no data" for a
line we hold (RYA-708, filed because RYA-707 published exactly that about Al 7835/8772)
reports **zero instruments** for 13123.42 and 13150.75 — Al's two best-graded lines. It
also cannot see the corrected HARPS holding RYA-986 is running on right now.

The census works around it with an explicit `COVERAGE_BLIND_SPOT` table, but the fix
belongs in `coverage.py`: those skips should be loud, or the registry should carry the
loader.

### Real-pixel coverage, measured (header endpoints do not prove coverage)

Loaded all 18 CRIRES+ Vesta IDPs and tested actual quality-flagged pixels:

| line (air) | λ_vac recovered | setting | good px in ±1 Å |
| --- | --- | --- | --- |
| 10872.97 | 10875.951 | Y1028 | 59/59 |
| 10891.74 | 10894.720 | Y1028 | 60/60 |
| 11253.19 | 11256.271 | J1226 | 55/55 |
| 11254.92 | 11258.006 | J1226 | 56/56 |
| 12749.91 | 12753.397 | J1226 | 48/48 |
| 12757.27 | 12760.765 | J1228 | 48/48 |
| **13123.42** | **13127.006** | J1228 | **49/49** |
| **13150.75** | **13154.350** | J1226 | **50/50** |

All eight land on full, unflagged pixels. The recovered vacuum wavelengths agree with
Burheim's published λ_vac to **< 0.005 Å**, which closes the air↔vacuum chain end to end
and independently confirms the line identifications.

**This overturns RYA-716's "do not queue Al for CRIRES+".** That recommendation was
correct on its own evidence — Elgueta 2026's curated G-dwarf robust flag certifies 0 of
12 Al records — but it is *their* selection, and it predates knowing that Burheim grades
13123/13150 at 1.5 % and 3.1 %. The comment itself named "our own VALD selection with our
own quality gates" as the alternative needing justification. This is the justification.

### Telluric

Confirmed against the **live** gate (`telluric_policy.gate_holding`), which discriminates
— 2 refusals out of 8 holdings tested:

* **`solar_vesta_crires_plus_idp` → REFUSED.** `telluric_applied=not-applied`,
  `crires_plus` is `telluric_required=yes`. **Every NIR Al line above 10650 Å is gated
  behind a molecfit run on this holding.** This is where RYA-993 (well-mixed column) and
  RYA-998 (star-own-RV mask) get re-exercised.
* `solar_crires_plus_y_rya794` → passes, but covers only 10280–10680 Å and reaches no Al
  target line.
* `solar_delbouille_liege` → REFUSED (`not-applied`, `telluric_required=yes`), which
  removes the only Sirius-side cover for the two near-UV lines.
* HARPS / IAG / Kitt Peak → pass on the RYA-460/786 per-line clean-line basis. VIS and
  red-optical Al need no correction stage.

---

## Recommended Phase 1 scope and ordering

1. **Ingest Burheim 2023 as Al's primary lab-gf table** (`data/reference/al_gf_lab/`) and
   generalise `gf_grades.GF-LAB` off the Fe-only `fe1_lab_loggf.csv`. Everything else
   depends on this and nothing owns it today. Provenance case only — RYA-161.
2. **Molecfit the CRIRES+ Vesta J settings** (J1226 / J1228 / J1232). Gates every NIR Al
   line and re-exercises RYA-993 + RYA-998 on a new instrument/element mix, which is the
   cross-check this ticket exists for.
3. **IR-synth on 13123.42 + 13150.75.** Deep (0.474 / 0.444), 6- and 4-component HFS,
   experimental gf at 1.5 % / 3.1 %. The best-conditioned Al lines in the project, and
   entirely new ground.
4. **Then 11253.19 + 11254.92** (5 % / 2 %). Caveat to carry: 11254.891 and 11254.926 are
   0.035 Å apart and blend at R ≈ 86 000; Burheim's +0.327 is the **strong component
   alone** (he omits the weak partner as "more than an order of magnitude weaker"), while
   the observed feature sums to +0.354. Do not silently use one for the other.
5. **VIS re-measure on 6696.015 + 6698.673** with the summed-HFS Burheim gf — the two
   lines the firewall was believed to have emptied. Note this is *not* 6696.185, the line
   currently in the EW pool.
6. **Only then** 7835/7836 + 8772/8773, which is where the old work sat. They have no
   experimental gf and never will from Burheim; NIST B/B+ is the ceiling, and RYA-838's
   premise needs rewriting first given OQ2.
7. **External validation leg, not a tuning target:** Nordlander & Lind 2017 publish
   ⟨3D⟩-NLTE / ⟨3D⟩-LTE / 1D-LTE abundances on 7835/8772. Read the final value against
   them; never let them select lines or gf (RYA-161).

**Do not start Phase 1 on the VIS EW route.** The pool is 2 lines, one is graded off the
wrong transition, and both return A ≈ 7.2–7.4 against a 6.43 anchor. That is the RYA-523
saturation arm, and it is the documented poor arm for a reason.

---

## Lines orphaned, with reason

| lines | reason |
| --- | --- |
| 324 FUV features (< 3000 Å), incl. all 35 Al II and both Al III | no solar holding below 2960 Å, **and** no declared band policy. Acquisition question. |
| 1 feature > 24000 Å (24985.95) | beyond every held instrument |
| 2 near-UV (3059.03, 3059.92) | have a home only on Delbouille, which the live telluric gate **refuses** |
| 12749.91, 12757.27 | Burheim-graded but central depth 0.015 / 0.008 — far below the 0.05 floor. Graded and unusable. |
| 126 shallow features across VIS/red-opt/NIR | central depth < 0.05 |
| 10 saturated (8 near-UV + 2 VIS, incl. the 3944/3961 resonance doublet) | central depth > 0.60; not abundance lines for the Sun |

---

## Files

| path | what |
| --- | --- |
| `scripts/rya1001_al_census.py` | the census; wavenumber HFS collapse, gf ladder, tiering, positive control |
| `scripts/rya1001_crires_coverage.py` | real-pixel CRIRES+ coverage test |
| `data/results/rya1001/rya1001_al_line_census.csv` | 502 features × 47 columns |
| `data/results/rya1001/rya1001_band_summary.csv` | the per-band table above |
| `data/results/rya1001/rya1001_crires_coverage.csv` | the 8-line pixel test |
| `data/results/rya1001/rya1001_census_meta.json` | thresholds, ladder, positive-control result |
| `data/results/rya1001/RYA1001_AL_PHASE0_CENSUS.md` | this report |

Bibliography: the library audit is clean (59/59 documents rowed). Four rows added for
sources this census cites and the repo did not carry — Burheim 2023, Papoulia 2019,
Mendoza 1995, Nordlander & Lind 2017 — all DOIs Crossref-verified.

**NOT merged.**
