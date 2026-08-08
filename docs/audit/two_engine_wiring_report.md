# Two-engine wiring audit — all canonical species (RYA-673)

**GENERATED — do not hand-edit.** Regenerate with `python scripts/rya673_two_engine_wiring_audit.py`.

**Discovery only.** No wiring was changed, no engine touched, no verdict regenerated. Every "wired" below means *the orchestrator's own coverage function returned a value for this species on real solar data* — not that a code path appears to exist.

`both` **10** · `A_only` **0** · `B_only` **6** · `neither` **11** — 27 species

## The headline

**11 of 27 species have no Engine B.** Of those, **7 are synthesis-required** — elements whose raw-EW leg is deliberately suppressed by RYA-520, so Engine B is not their second opinion, it is their *only* leg. Those rows are reporting on one engine or on nothing.

| species | treatment | why Engine B is missing |
|---|---|---|
| Mg I | `synthesis` | `DELIBERATELY_SKIPPED` |
| N I | `per_region_source` | `NO_MODEL_ATOM` |
| Co I | `synthesis` | `NO_HARNESS_INVOCATION` |
| Ba II | `synthesis` | `NO_HARNESS_INVOCATION` |
| Sc II | `HFS_sum` | `NO_MODEL_ATOM` |
| Eu II | `HFS_sum` | `NO_MODEL_ATOM` |
| Zr II | `synthesis` | `DELIBERATELY_SKIPPED` |

## ⚠ The RYA-525 loud-fail does not cover this

RYA-525 added a guard whose stated job is to refuse a synthesis-required element with no Engine-B value:

```python
species = sorted(set(a_pl) | set(b_pl) | set(ded_b), ...)
for (el, ion) in species:
    ...
    if synth_required and not has_B:
        loud.append(...)          # -> SystemExit
```

The guard iterates the **union of the three coverage sources**. An element absent from all three never enters the loop, so it is never tested — it is skipped in silence. The guard catches a *partially* covered species and is blind to a *completely* uncovered one, which is the more serious case.

Species that produce **no two-engine record at all**, and therefore never reach the guard: **Mg I, P I, N I, Co I, Al I, K I, Ba II, Y II, Sc II, Eu II, Zr II**.

## Per-species

| species | status | A | B | A reason | B reason | treatment | blocks Beta |
|---|---|---|---|---|---|---|---|
| Fe I | `both` | ✓ | ✓ | — | — | `astrophysical_gf_differential` | — |
| Fe II | `both` | ✓ | ✓ | — | — | `astrophysical_gf_differential` | — |
| C I | `B_only` | ✗ | ✓ | DELIBERATELY_SKIPPED | — | `exclude` | — |
| O I | `B_only` | ✗ | ✓ | DELIBERATELY_SKIPPED | — | `cited_substitution` | — |
| Mg I | `neither` | ✗ | ✗ | DELIBERATELY_SKIPPED | DELIBERATELY_SKIPPED | `synthesis` | **YES** |
| Si I | `both` | ✓ | ✓ | — | — | `astrophysical_gf_differential` | — |
| Ca I | `both` | ✓ | ✓ | — | — | `none` | — |
| Ti I | `both` | ✓ | ✓ | — | — | `astrophysical_gf_differential` | — |
| Ni I | `both` | ✓ | ✓ | — | — | `astrophysical_gf_differential` | — |
| Na I | `both` | ✓ | ✓ | — | — | `none` | — |
| P I | `neither` | ✗ | ✗ | NO_EW_POOL | NO_MODEL_ATOM | `per_region_source` | — |
| S I | `both` | ✓ | ✓ | — | — | `none` | — |
| N I | `neither` | ✗ | ✗ | DELIBERATELY_SKIPPED | NO_MODEL_ATOM | `per_region_source` | **YES** |
| Co I | `neither` | ✗ | ✗ | DELIBERATELY_SKIPPED | NO_HARNESS_INVOCATION | `synthesis` | **YES** |
| Cr I | `both` | ✓ | ✓ | — | — | `astrophysical_gf_differential` | — |
| Al I | `neither` | ✗ | ✗ | NO_EW_POOL | NO_MODEL_ATOM | `none` | — |
| K I | `neither` | ✗ | ✗ | NO_EW_POOL | NO_MODEL_ATOM | `NLTE_grid` | — |
| Ba II | `neither` | ✗ | ✗ | DELIBERATELY_SKIPPED | NO_HARNESS_INVOCATION | `synthesis` | **YES** |
| Y II | `neither` | ✗ | ✗ | NO_EW_POOL | NO_MODEL_ATOM | `per_region_source` | — |
| V I | `B_only` | ✗ | ✓ | NO_EW_POOL | — | `none` | — |
| Cu I | `B_only` | ✗ | ✓ | DELIBERATELY_SKIPPED | — | `synthesis` | — |
| Mn I | `B_only` | ✗ | ✓ | DELIBERATELY_SKIPPED | — | `synthesis` | — |
| Sc II | `neither` | ✗ | ✗ | DELIBERATELY_SKIPPED | NO_MODEL_ATOM | `HFS_sum` | **YES** |
| Li I | `both` | ✓ | ✓ | — | — | `upper_limit` | — |
| Eu II | `neither` | ✗ | ✗ | DELIBERATELY_SKIPPED | NO_MODEL_ATOM | `HFS_sum` | **YES** |
| Zr II | `neither` | ✗ | ✗ | DELIBERATELY_SKIPPED | DELIBERATELY_SKIPPED | `synthesis` | **YES** |
| Sr II | `B_only` | ✗ | ✓ | DELIBERATELY_SKIPPED | — | `synthesis` | — |

## What each reason class means, and the follow-on it implies

- **`NO_HARNESS_INVOCATION`** (Co I, Ba II) — the measurement **already exists** in the repo and the orchestrator simply never reads it. Cheapest possible class — a wiring ticket, no new science.
- **`NO_MODEL_ATOM`** (P I, N I, Al I, K I, Y II, Sc II, Eu II) — no validated Engine-B NLTE atom for the species. Needs atom sourcing + an RYA-534/548-style anchor validation BEFORE wiring is meaningful.
- **`NO_NLTE_GRID`** (none) — no 1D delta CSV for the Engine-A leg. Affects quality (the leg runs LTE), not wiring — no row is marked unwired for this alone.
- **`NO_EW_POOL`** (P I, Al I, K I, Y II, V I) — the curated EW pool has no surviving line for the species. A measurement gap, not a plumbing gap — line-pool / gf work, not a wiring ticket.
- **`DELIBERATELY_SKIPPED`** (C I, O I, Mg I, N I, Co I, Ba II, Cu I, Mn I, Sc II, Eu II, Zr II, Sr II) — a documented, ratified decision — the RYA-520 raw-EW suppression, or a harness that is wired and gated shut. **Not** a gap; do not file a ticket.
- **`UNKNOWN`** (none) — the audit could not determine the cause from the code. Deferred to Ryan by design rather than guessed.

## Recommended follow-on per non-`both` species

Recommendations only. **No tickets filed** — per §4, Ryan directs which get filed.

| species | recommended action |
|---|---|
| C I | None — ratified decision, working as designed. |
| O I | None — ratified decision, working as designed. |
| Mg I | **BLOCKS BETA** — synthesis-required with no Engine B (`DELIBERATELY_SKIPPED`). Needs the underlying gap closed first. |
| P I | Line-pool / gf work, not wiring. No wiring ticket. |
| N I | **BLOCKS BETA** — synthesis-required with no Engine B (`NO_MODEL_ATOM`). Needs the underlying gap closed first. |
| Co I | **FILE A WIRING TICKET** — result exists, orchestrator does not read it. Blocks Beta: this species currently reports on no engine. |
| Al I | Line-pool / gf work, not wiring. No wiring ticket. |
| K I | Line-pool / gf work, not wiring. No wiring ticket. |
| Ba II | **FILE A WIRING TICKET** — result exists, orchestrator does not read it. Blocks Beta: this species currently reports on no engine. |
| Y II | Line-pool / gf work, not wiring. No wiring ticket. |
| V I | Line-pool / gf work, not wiring. No wiring ticket. |
| Cu I | None — ratified decision, working as designed. |
| Mn I | None — ratified decision, working as designed. |
| Sc II | **BLOCKS BETA** — synthesis-required with no Engine B (`NO_MODEL_ATOM`). Needs the underlying gap closed first. |
| Eu II | **BLOCKS BETA** — synthesis-required with no Engine B (`NO_MODEL_ATOM`). Needs the underlying gap closed first. |
| Zr II | **BLOCKS BETA** — synthesis-required with no Engine B (`DELIBERATELY_SKIPPED`). Needs the underlying gap closed first. |
| Sr II | None — ratified decision, working as designed. |

## Per-species narrative

### Fe I — `both`

A: 81 EW-pool line(s). B: Engine B via synth-v2 per-line

### Fe II — `both`

A: 3 EW-pool line(s). B: Engine B via synth-v2 per-line

### C I — `B_only`

A: raw-EW leg suppressed by RYA-520 (required_treatment=exclude) — Engine B is this species' only valid leg. B: Engine B via dedicated harness (synth-v2 also covers it)

### O I — `B_only`

A: raw-EW leg suppressed by RYA-520 (required_treatment=cited_substitution) — Engine B is this species' only valid leg. B: Engine B via dedicated harness

### Mg I — `neither`

A: raw-EW leg suppressed by RYA-520 (required_treatment=synthesis) — Engine B is this species' only valid leg. B: harness IS wired (RYA-592) and is gated shut: it returns nothing until its reliability/concordance condition is met. Wiring is NOT the fix — the measurement is.

### Si I — `both`

A: 15 EW-pool line(s). B: Engine B via synth-v2 per-line

### Ca I — `both`

A: 4 EW-pool line(s). B: Engine B via synth-v2 per-line

### Ti I — `both`

A: 6 EW-pool line(s). B: Engine B via synth-v2 per-line

### Ni I — `both`

A: 3 EW-pool line(s). B: Engine B via synth-v2 per-line

### Na I — `both`

A: 2 EW-pool line(s). B: Engine B via synth-v2 per-line

### P I — `neither`

A: no line survives the curated EW pool for this species. B: no validated Engine-B NLTE atom: no RYA-534 Engine-B grid provenance on record for P

### S I — `both`

A: 1 EW-pool line(s). B: Engine B via synth-v2 per-line; ⚠ the ratified RYA-492 harness result (data/results/solar_s_costasilva_rya492.json) is NOT referenced by the orchestrator

### N I — `neither`

A: raw-EW leg suppressed by RYA-520 (required_treatment=per_region_source) — Engine B is this species' only valid leg. B: no validated Engine-B NLTE atom: no RYA-534 Engine-B grid provenance on record for N

### Co I — `neither`

A: raw-EW leg suppressed by RYA-520 (required_treatment=synthesis) — Engine B is this species' only valid leg. B: the RYA-564 synthesis result already exists at data/results/co_synthesis_rya564.json — the orchestrator never reads it

### Cr I — `both`

A: 3 EW-pool line(s). B: Engine B via synth-v2 per-line

### Al I — `neither`

A: no line survives the curated EW pool for this species. B: no validated Engine-B NLTE atom: no RYA-534 Engine-B grid provenance on record for Al

### K I — `neither`

A: no line survives the curated EW pool for this species. B: no validated Engine-B NLTE atom: no RYA-534 Engine-B grid provenance on record for K

### Ba II — `neither`

A: raw-EW leg suppressed by RYA-520 (required_treatment=synthesis) — Engine B is this species' only valid leg. B: the RYA-559 synthesis result already exists at data/results/solar_ba_synthesis_rya559.json — the orchestrator never reads it

### Y II — `neither`

A: no line survives the curated EW pool for this species. B: no validated Engine-B NLTE atom: no RYA-534 Engine-B grid provenance on record for Y

### V I — `B_only`

A: no line survives the curated EW pool for this species. B: Engine B via dedicated harness

### Cu I — `B_only`

A: raw-EW leg suppressed by RYA-520 (required_treatment=synthesis) — Engine B is this species' only valid leg. B: Engine B via dedicated harness

### Mn I — `B_only`

A: raw-EW leg suppressed by RYA-520 (required_treatment=synthesis) — Engine B is this species' only valid leg. B: Engine B via dedicated harness

### Sc II — `neither`

A: raw-EW leg suppressed by RYA-520 (required_treatment=HFS_sum) — Engine B is this species' only valid leg. B: no validated Engine-B NLTE atom: no RYA-534 Engine-B grid provenance on record for Sc

### Li I — `both`

A: 1 EW-pool line(s). B: Engine B via synth-v2 per-line; ratified UPPER_LIMIT (RYA-563) — may never carry a point value

### Eu II — `neither`

A: raw-EW leg suppressed by RYA-520 (required_treatment=HFS_sum) — Engine B is this species' only valid leg. B: no validated Engine-B NLTE atom: no RYA-534 Engine-B grid provenance on record for Eu

### Zr II — `neither`

A: raw-EW leg suppressed by RYA-520 (required_treatment=synthesis) — Engine B is this species' only valid leg. B: harness IS wired (RYA-560) and is gated shut: it returns nothing until its reliability/concordance condition is met. Wiring is NOT the fix — the measurement is.

### Sr II — `B_only`

A: raw-EW leg suppressed by RYA-520 (required_treatment=synthesis) — Engine B is this species' only valid leg. B: Engine B via dedicated harness

