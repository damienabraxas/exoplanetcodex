# RYA-408 — Solar-EW canonical input source: Step-0 trace + fix

Root-cause fix for the RYA-406 incident: the solar Fe gate / abundance derivation read
its EW input from the **gitignored, regenerable** staging file `data/processed/solar_ew.csv`,
whose per-worktree content can silently diverge from the committed canonical. This is why
staging a stale runtime reddened the stewardship suite during RYA-406.

## Step 0 — read-path trace (resolves the RYA-397-vs-352 contradiction)

Every in-code consumer of a solar-EW table:

| Consumer | Path read | Gitignored? | Canonical? |
|----------|-----------|-------------|------------|
| `abundances_derive._load_solar_ews` (THE Fe gate / abundance EW pool) | `PATHS['solar_ew']` → `data/processed/solar_ew.csv` | **YES** | **NO** |
| `validate_fe_rya238` (gate; calls `abundances_derive.run()` → `_load_solar_ews`) | (same, via loader) | **YES** | **NO** |
| `diagnostics_abundance` | `PATHS['solar_ew']` | YES | NO |
| `check_stewardship.check_blend_flag` (propagation leg) | `PATHS['solar_ew']` | YES | NO (present-gated) |
| `scripts/build_fe2_ges_regions.py`, `build_fe_nlte_grid_rya319.py`, `rya305_fe2_triage.py`, `fe2_ew_quality_cull.py` | `PATHS['solar_ew']` | YES | NO |
| `pipeline/curate_nonfe_pools.py` (`_EW_POOL`) | `data/measured/sol_ew_results_v1.csv` | NO | **YES** |
| `scripts/build_nlte_grids_mpia.py`, `fetch_inspect_nlte.py` (RYA-396 scrapers) | `data/measured/sol_ew_results_v1.csv` | NO | **YES** |

**Verdict (CONFIRMED, not assumed):** the gate/abundance path read the **gitignored
runtime** — refuting RYA-397's "every live code path reads the canonical" and confirming
RYA-352/405/406. The RYA-396 scrapers and non-Fe pool curation already read the canonical;
the gate/abundance loader was never repointed. So Step 0 does **not** STOP — it proceeds.

## Step 1 — repoint to the committed canonical

`config.constants.PATHS['solar_ew_canonical'] = data/measured/sol_ew_results_v1.csv`
(new single-source constant). `_load_solar_ews` now reads it for the Fe II + non-Fe-I
pool (Fe I stays on the committed GES reference, RYA-330; unchanged). The read **loud-fails**
if the canonical is missing — it never falls back to the gitignored staging file.

Honest gate impact (validate-don't-tune; numbers **not** tuned):

| Quantity | runtime pool (pre-408) | canonical pool (post-408) |
|----------|------------------------|---------------------------|
| Fe II clean lines (post 352-cull + 305-triage) | 7 | **11 → 3 clean** (curated set) |
| Fe I slope (GES, unchanged) | −0.011 PASS | −0.011 PASS |
| Fe I−Fe II ionization (EW path) | −0.184 | −0.141 |
| A(Fe II) NLTE abs (EW path) | 7.700 | 7.657 |

The ionization / Fe II-absolute still read "high" on the EW path — **by design** (RYA-352;
the ratified ionization arbiter is **synthesis**, not the EW path — that gate-architecture
fix is RYA-406, coordinated separately). 408 owns only the **input source**: the gate now
runs entirely off committed inputs with **zero runtime staging**.

## Step 2 — demote runtime to staging + reviewed promotion

`data/processed/solar_ew.csv` is documented as **staging only** (regenerable `lines_fit`
output, never a gate/abundance input). Promotion staging→canonical is the reviewed
`scripts/promote_solar_ew.py`: dry-run by default; **STOPS writing nothing** on any
`blend_flag` conflict between staging and the canonical (the canonical's 11 vetted blends
are curation the raw staging does not own); promotes EW measurements; preserves canonical
coverage; reports — never auto-adds — staging-only lines. A genuine `lines_fit` flag bug is
a follow-up (fix `VETTED_BLENDS`), not a promotion-time overwrite.

## Step 3 — drift / identity guard

New stewardship invariant `check_solar_ew_canonical`:
1. **present + well-formed** — canonical exists under `data/measured`, has the schema, non-empty (else loud `StewardshipParseError`);
2. **IDENTITY** — `_load_solar_ews` source must read `solar_ew_canonical` and must **not** read `PATHS['solar_ew']` as its pool (a re-point to the runtime = UNTRACKED loud break);
3. **DRIFT** — if the regenerable staging file is present, every canonical line must match it on the **measured EW** (±0.5 mÅ); a divergence = stale / different-run staging (UNTRACKED). `blend_flag` is intentionally not compared here (it is curation owned by the canonical; its integrity is `check_blend_flag`'s job).

## Step 4 — pin the filename

Only `sol_ew_results_v1.csv` exists in-repo; there is no `sol_ew_results_v1_1.csv` (the
ambiguity was a stray reference, already negated in `data/measured/README.md`). The README
now states v1 is the single canonical and documents the RYA-408 contract + promotion.

## O I 6300.304 blend_flag — resolved (not muted)

The EW canonical carried `blend_flag=True` (correct — forbidden [O I] is blended with
Ni I 6300.339, RYA-104/208), but `linelist_solar.csv` + the vetted builder said `False`
(O I 6300.3 was **missing from `VETTED_BLENDS`**). Resolved by adding O I 6300.304 to
`VETTED_BLENDS` with the RYA-104/208 citation + a one-cell surgical patch of the line-list
flag (`False→True`, byte-preserving). Now linelist, vetted-builder, and EW canonical all
agree on **True** — 0 blend_flag definition mismatches. (The line-list's stale −9.776 gf is
out of scope — RYA-367 ruled the #2-store gf non-load-bearing.)
