# RYA-350 — Canonical single-source `loggf` design

**Status:** DESIGN (Option A, ratified 2026-06-18). Implementation is **RYA-353** — do
not build the migration from this ticket. This document is the deliverable; RYA-350
closes when this design is signed off.

**Scope decision recap.** RYA-350 Step 0 found the synth-vs-`linelist_solar` `loggf`
divergence is **systemic** (23.6% of 13,984 cross-matched physical lines diverge
>0.05 dex, iron-peak dominated; Fe II only 16 of ~3,299). The ratified fix is
**Option A — structural single-source de-duplication**: one authoritative `loggf`
per physical line, both engines resolve from it, neither line-list keeps its own
copy. Per-line value *quality* (which catalog is right) is **out of scope here** —
that is RYA-354, run incrementally on top of the single source. This design only
removes the duplication mechanism.

---

## 1. The defect, precisely

`loggf` for one physical transition is stored independently in (at least) three
places, and they drift:

| # | Store | Column | Reader | Catalog | Consumed by |
|---|-------|--------|--------|---------|-------------|
| 1 | `ispec/.../GESv6_atom_hfs_iso.420_920nm/atomic_lines.tsv` (shared, outside repo) | `loggf` | `ispec.read_atomic_linelist()` — `abundances_derive.py:537` | GES v6 | flux synthesis (synth-EW RYA-285, synth-v2 flux fit) |
| 2 | `data/linelists/linelist_solar.csv` (+ per-star lists) | `log_gf` | `data/linelists/loader.py:63`; `lines_fit.py:363` (COG); `diagnostics_abundance.py:82` | VALD3 | Ni 6300 COG blend, Fe I diagnostic, line scoring |
| 3 | `ispec/.../regions/42000_VALD/limited_but_with_missing_elements_..._all_extended.txt` (shared, outside repo; 524 curated lines) | `loggf` | `ispec.read_line_regions()` → `_build_ispec_line_regions` `:198` → `determine_abundances` → `last_linemasks['loggf']` `:1588` | **VALD3** | **the live EW→abundance value** (every `*_per_line.csv` A(X)) |

> **Store-model correction (RYA-353 loader trace, 2026-06-18).** Store #3 is **not**
> GES "regions built from #1" (the `_build_ispec_line_regions` docstring saying "GES
> regions" is a misnomer — the file lives under `42000_VALD/`). Its `loggf` matches
> `linelist_solar` (VALD3) **exactly** on all 5 RYA-347 anchors (e.g. 6247.557 = −2.310,
> vs synth GES −2.435), confirmed by reading the file. So #3 is a **third independent
> VALD3 copy**, and it — not `linelist_solar.csv` (#2) — is what sets published EW
> abundances. Worse, **#2 and #3 have themselves drifted**: 201 of 524 region lines
> diverge >0.01 dex from `linelist_solar` (some HFS-artifact in a non-HFS-aware check,
> but real cases exist, e.g. Fe I 4924.301 Δ+0.132). Net: three genuinely independent
> copies — GES (#1, synth) vs VALD3-regions (#3, EW abundance) vs VALD3-linelist (#2,
> COG/diag) — and the synth-vs-EW cross-path divergence is **#1 vs #3**.

Store #1 (GES) feeds synth; stores #2 and #3 are both VALD3 but separate, drifted
copies, with #3 driving the live EW abundance (see correction box above). Because
`loggf` is degenerate with the floated abundance
(RYA-347: changing 6247.557's synth `loggf` moved A(Fe II) 7.552→7.467 at χ²ᵣ
unchanged 180→179), a cross-store Δgf is a **pure differential abundance bias**
between engines — invisible to χ²ᵣ, silently corrupting any cross-path comparison.

**What it does *not* bias.** Differential `[X/H] = A(star) − A(Sun)` on the same line
measured the same way cancels `loggf` exactly. So our published differential results
are sound *as long as each path is internally consistent*. The duplication bites at:
(a) cross-path comparison (RYA-349 — Δgf masquerades as engine/physics), (b) the
**EW↔synth escalation asymmetry** (a line escalating to synth for the star but not the
Sun breaks cancellation), (c) absolute-scale checks. This bounds urgency: real, worth
fixing structurally, **not** gating RYA-237.

---

## 2. Design: one `loggf` table, foreign-keyed

### 2.1 Canonical store

A single tracked, reviewable table **`data/linelists/canonical_gf.csv`**, one row per
**physical line**, keyed by the RYA-345 canonical identity:

| column | meaning |
|--------|---------|
| `line_id` | stable surrogate key (string, e.g. `gf_000123`), assigned once, never reused |
| `z`, `ion` | atomic number + ion stage (RYA-345 `species_key`); molecules → `z=NULL`, `mol_name` set |
| `mol_name` | molecule name for molecular lines, else empty |
| `wavelength_air_A` | air λ (centroid; see HFS §2.3) |
| `excitation_potential_eV` | lower-level EP — the secondary identity discriminant |
| `log_gf` | **the** authoritative value (total gf per physical line) |
| `loggf_reference` | the *underlying* source — the loggf paper/code (`K07`, `FMW`, `NIST ASD v5.11`), **not** "GES v6"/"VALD3" (stewardship rule) |
| `nist_grade` | NIST ASD grade (A/B/C/…) if the value is NIST-sourced, else empty |
| `hfs_n_components` | how many HFS components this physical line represents (1 = none) |
| `adjudication_status` | `seeded` (migration default) \| `nist` \| `pending` (divergent, no graded source — flagged for RYA-354) |
| `provenance_note` | free text: seed origin, Δgf at seed time, manual-pull flag |

`(z, ion, mol_name, wavelength_air_A, excitation_potential_eV)` is the natural key;
`line_id` is the surrogate both line-lists point at, so runtime resolution is an O(1)
join with **no tolerance ambiguity** (the fuzzy λ/EP match happens once, at migration,
not on every load).

### 2.2 Resolution: both paths read from the canonical store

Each source store **drops its own `loggf`/`log_gf` column** and instead resolves it
from `canonical_gf` (added once at migration). **All three** stores must be rerouted —
the trace shows three independent copies, not two:

- **Synth path** (store #1, `_load_synth_resources`): after `read_atomic_linelist`,
  overwrite the in-memory `loggf` by resolving against `canonical_gf` (species_key + λ +
  EP, or a pre-annotated `gf_line_id`).
- **EW→abundance path** (store #3, `_build_ispec_line_regions` `:198` → the 524-line
  `42000_VALD/...all_extended.txt` regions): overwrite the regions' `loggf` from
  `canonical_gf` **before** `determine_abundances`. **This is the load-bearing one** —
  it sets every published EW abundance. Once #1 and #3 both resolve from the same source,
  the synth-vs-EW gf divergence (RYA-349) is gone by construction.
- **COG/diagnostics path** (store #2, `loader.py:load_linelist`): populate `log_gf` from
  the same source; COG (`lines_fit.py`) and diagnostics inherit it.

Both shared `ispec/` files (#1, #3) live outside the repo, so the reroute happens in the
**reader** (overwrite-after-read), not by editing the shared files in place.

**No silent fallback** (project rule): a line that fails to resolve against
`canonical_gf` **raises** — it does not default. This reuses the RYA-345 0-match
loud-guard pattern.

### 2.3 HFS granularity — the one real design fork

The stores disagree on HFS representation: store #1 carries a single total-`gf`
line; store #2 carries N VALD3 components each at a fraction (this produced the fake
Na/Sc/La/K offsets in Step 0, removed by total-gf aggregation). The 524-line regions
file (#3) is curated and largely one-line-per-feature like #1. The abundance pipeline
already treats HFS features as **single absorption features** — "the EW from a single
Voigt fit = total HFS EW" (`docs/linelist_pipeline.md`, Ba II/Eu II/Li I). Therefore:

> **Canonical granularity = one row per physical line carrying the *total* gf**
> (`log_gf = log₁₀ Σ 10^{gf_component}`), which is the abundance-relevant quantity.

At migration, `linelist_solar` HFS components map to the single canonical `line_id`
for their feature (RYA-353 decides whether to collapse the rows or have all components
share one `gf_line_id`; no consumer needs per-component gf for abundance — to be
confirmed in the RYA-353 loader trace). `hfs_n_components` preserves the multiplicity.

### 2.4 Seed rule (migration default — value *quality* is RYA-354)

For each canonical line, seed `log_gf` to the **best current per-line value**, chosen
to **preserve published differential results by construction**:

1. **NIST-graded available** (from `linelist_solar.nist_grade` / `nist_crosscheck.csv`)
   → use NIST, `adjudication_status=nist`, cite grade. (Includes 6247.557, §4.)
2. **Stores agree** (|Δgf| ≤ 0.01) → that value; `adjudication_status=seeded`; cite the
   underlying reference.
3. **Stores diverge, no NIST grade** → seed to the value of the line's **primary
   measurement path** (synth-measured lines → GES store #1; EW-only lines → VALD3
   store #2), so the path that actually produces the number is unchanged. Record BOTH
   catalog values, the Δgf, and `adjudication_status=pending` → handed to RYA-354.
4. **Line in one store only** → that value with its provenance.

Consequence to state plainly: for a *divergent* line, the **non-primary** path's
absolute value changes (it now equals the primary). That is the unavoidable point of
de-duplication — you cannot keep two values. It does **not** move differential `[X/H]`
on path-consistent lines (gf cancels), and every such change is logged with
before/after + provenance. No value is invented, approximated, or pulled from memory;
a line needing a value with no citable source is flagged `pending`, never substituted.

---

## 3. Migration plan (for RYA-353 — not built here)

1. **Loader trace.** ✅ DONE (RYA-353, 2026-06-18) — see the §1 store-model correction
   box. Live EW→abundance gf = store #3 (the 524-line `42000_VALD/...all_extended.txt`
   regions, VALD3), **not** `linelist_solar` (#2); three independent copies confirmed.
   The migration reroutes all three readers (§2.2).
2. **Build `canonical_gf.csv`.** Reuse `scripts/audit_gf_duplication.py`'s match
   (species_key + λ + EP, HFS-aware total-gf aggregation) to union **all three** stores
   into physical lines, assign `line_id`, apply the §2.4 seed rule, emit provenance.
3. **Annotate source files** with `gf_line_id` (one-time, idempotent, scripted).
4. **Reroute loaders** (§2.2) behind the join; keep the old in-file gf column for **one
   transition commit** as a shadow with an assert-equal, then delete it.
5. **Validate** (§3.1) and land canonical-data changes on a review branch — Ryan merges.

### 3.1 Validation plan (acceptance for RYA-353)

- **Differential stability:** solar + per-star `[X/H]` on path-consistent lines
  unchanged within rounding (the de-dup must not move published differentials).
- **Cross-path agreement:** re-run `audit_gf_duplication.py` → **zero** duplicated-and-
  divergent gf remaining (every matched line resolves to one `line_id`/value).
- **Loud-guard:** an unresolved `gf_line_id` raises (negative test).
- **Provenance completeness:** every `canonical_gf` row has a non-empty
  `loggf_reference`; count `adjudication_status=pending` (the RYA-354 backlog).
- **6247.557 record:** before/after gf + abundance logged as provenance only (retired /
  blend-flagged — must not move a clean pool; do not anchor verification on it).

---

## 4. 6247.557 — catalog-correctness one-liner (specified; not applied here)

In scope for the gf-reconciliation track as a provenance fix, **independent of the
migration**. The synth store carries the more-wrong value:

| store | current `loggf` | correct | source |
|-------|-----------------|---------|--------|
| `atomic_lines.tsv` Fe II 6247.557 (EP 3.892) | **−2.435** | **−2.329** | NIST ASD v5.11, grade **B** |
| `linelist_solar.csv` Fe II 6247.557 | −2.31 | (already ≈NIST) | VALD3 |

Edit = change the single `atomic_lines.tsv` row `loggf −2.435 → −2.329` with an inline
NIST citation. **Not applied in this design ticket** because (a) RYA-350 is design-only
and (b) `atomic_lines.tsv` lives in the **shared `ispec/` tree outside the repo/worktree**
— editing it affects every worktree and is a canonical-data change needing explicit
review. It lands with the migration (RYA-353) or as a standalone reviewed edit, whichever
Ryan prefers. Effect (RYA-347 measured): A(Fe II) on this line 7.552→7.467 (−0.085 dex,
≈ −Δgf) at χ²ᵣ unchanged (180→179); **moves no pool** (6247 is retired from the synth
clean pool and blend-flagged in EW).

---

## 5. Interim guard (effective now, pre-migration)

Until the single source lands, to stop accreting onto the two-file mess:

> **Any new `loggf` consumer or new line added to either store must resolve gf through
> the planned canonical source — never add a fourth independent copy.** New lines get a
> `canonical_gf` entry (with provenance) first; consumers read the join, not a file
> column.

This is advisory until enforced. **RYA-355** makes it mechanical: a CI data-stewardship
invariant that fails loudly on any canonical value duplicated-and-divergent across files
(gf today; extensible to STAR_PARAMS, provenance). Run today it flags this gf defect
(proving it works); RYA-353 makes it pass.

---

## 6. Open questions for sign-off

1. ~~Store #3 vs #2 for the EW abundance~~ — **RESOLVED** by the RYA-353 loader trace:
   live EW gf = store #3 (regions file, VALD3); three independent copies; migration
   reroutes all three readers (§1 box, §2.2). **New decision this surfaces:** the canonical
   source must now also feed the shared regions file (#3) — confirm the overwrite-after-read
   approach (vs editing the shared `ispec/` regions file) is acceptable.
2. **HFS component rows** — collapse `linelist_solar` HFS components to one canonical
   line, or keep component rows all pointing at one `line_id`? (§2.3 — RYA-353 mechanic.)
3. **6247.557** — apply now as a standalone reviewed edit, or fold into RYA-353?
4. **`line_id` scheme** — opaque surrogate (recommended) vs composite natural key string.

---

### Spun-out tickets
- **RYA-353** — BUILD the migration (both paths → `canonical_gf.csv`). Blocked by this
  design + RYA-345. Decoupled from RYA-237.
- **RYA-354** — SCIENCE: per-line gf adjudication to best-available (NIST-graded
  preferred) *on* the single source; iron-peak first. Where `adjudication_status=pending`
  gets resolved.
- **RYA-355** — BUILD: CI stewardship invariant (no duplicated-divergent canonical
  values). Independently buildable; the "never silently again" guard.
