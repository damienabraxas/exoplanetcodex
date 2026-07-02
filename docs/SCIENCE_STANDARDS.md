# SCIENCE STANDARDS — authoritative decision records

Single source of truth for cross-cutting science-method decisions. New file (RYA-505);
flag for the RYA-179 docs-sync revision. Add sections as decisions are recorded; each
section cites its evidence and the ticket that set it.

---

## Hot-Teff NLTE grid coverage (F-star benchmarks: Procyon 6554 K, τ Boo 6400 K) — RYA-505

**Decision date:** 2026-07-02 · **Evidence:** `scripts/nlte_fstar_ceiling_rya505.py`
→ `data/results/nlte_fstar_ceiling_rya505.csv`. **Scope:** hot-benchmark fidelity upgrade;
NOT on the α Cen / 55 Cnc critical path (those are 5200–5800 K, already in-grid).

### Headline
The RYA-349 "6200–6500 K wall" is almost entirely an **under-loading artifact, not a grid
limit.** Our on-disk non-Fe grids are the benchmark-node **subsets** synthesised in RYA-410
(cool FGK dwarfs + 55 Cnc), not the full published grids. **9 of 12 non-Fe elements have a
published grid that already covers Procyon** — the fix is to re-load/re-synthesise the
Procyon node in the *same* family (self-consistent), not to mix codes. Only **Ba** is a
genuine real-limit clamp.

### Grid-selection hierarchy (mandatory order)
1. **Self-consistent extend** — re-synthesise/re-load the hot node from the *same* code +
   model-atom + model-atmosphere family we already wire. No new systematic. **Preferred.**
2. **MPIA-with-cross-check** — only where no self-consistent extension exists. MPIA/Bergemann
   grids are DETAIL/SIU on MAFAGS-OS (≠ our Amarsi/PySME + MARCS), so require a
   validate-don't-tune cross-check at overlap (the Sun + one cool star in both grids): the
   two grids' corrections must agree within **±0.05 dex** at overlap before the MPIA grid is
   trusted at hot Teff. Disagreement is a recorded finding, never a silent choice.
3. **Bounded clamp** — only for elements with no higher-ceiling grid at all: a Teff clamp at
   Procyon's **54 K** overshoot ONLY, monotonic-in-Teff elements only, with the flagged
   systematic propagated into the σ budget (RYA-282). Never extrapolate a 3D-NLTE correction
   past its grid; never a silent clamp (Fe/RYA-319 precedent).

### Per-element hot-Teff coverage map
| element | current family | on-disk ceil | published ceil (cited) | verdict |
|---|---|---:|---:|---|
| Na | Amarsi2020 GALAH · PySME · MARCS | 6200 | **8000** (Zenodo 3982506 v3 prov, 2500–8000) | self-consistent-extend |
| Mg | Amarsi2020 GALAH · PySME · MARCS | 6200 | **8000** (same release) | self-consistent-extend |
| Si | Amarsi2020 GALAH · PySME · MARCS | 6200 | **8000** | self-consistent-extend (near-LTE, low pri) |
| Al | Amarsi2020 GALAH · PySME · MARCS | 6200 | **8000** (v3 prov) | self-consistent-extend (near-LTE, low pri) |
| K  | Amarsi2020 GALAH · PySME · MARCS | 6200 | **8000** (v3 prov) | self-consistent-extend |
| S  | Amarsi2025 · PySME · MARCS | 6200 | **8000** (prov, 3000–8000) | self-consistent-extend |
| Ca | Mashonkina2017 · DETAIL · MAFAGS-OS | 6500 | unrecorded (subset) | **self-consistent-extend via Amarsi Ca `.grd`** (staged in `amarsi_galah/`, Family-A) with validate-don't-tune cross-check vs the MAFAGS-OS value |
| Mn | Bergemann MPIA · DETAIL/SIU · MAFAGS-OS | 6500 | 7000 (nlte.mpia.de survey) | **self-consistent-extend via Amarsi Mn `.grd`** (staged) w/ cross-check |
| Ti | Bergemann2011 MPIA · MAFAGS-OS | 6500 | 7000 (survey) | MPIA-with-cross-check (near-LTE, low pri) |
| Cr | Bergemann2010 MPIA · MAFAGS-OS | 6500 | 7000 (survey) | MPIA-with-cross-check |
| **Ba** | Korotin2015 · MULTI · MARCS | 6500 | **6500 (REAL LIMIT)** — Korotin2015 prov coverage 4000–6500 | **bounded clamp** (54 K, monotonic) or find a higher-ceiling Ba grid |
| Sr | Bergemann2012 INSPECT · MARCS (metal-poor) | 6000 | unrecorded | defer (off-critical; INASAN primary pending RYA-433) |

Priority order for execution (NLTE-significant at F-star Teff first): **Ca, Ba, Na, Mn**,
then Mg/K/S; Si/Al/Ti near-LTE and low priority.

### CNO (Step 1) and banked Procyon O (Step 2) — already correct, confirmed
- **C/O NLTE:** `pipeline/nlte_cno.py` already splits legs at `TEFF_3D_CEILING = 6500 K`:
  3D leg (tables 2/3) ≤ 6500 K; **1D-NLTE leg (tables 5/6) above** → **Procyon 6554 K uses
  the 1D-NLTE leg, which spans Teff 4000–8000 K (17 nodes) → in-grid.** No change (RYA-359).
  3D refinement stays future v2 (RYA-444/445), matching the IR v1/v2 boundary (RYA-504).
- **RYA-483 banked Procyon O:** primary indicator `OI_777_1D_NLTE`, A(O)=8.82, [O/H] +0.085
  — **used the 1D-NLTE correction, in-grid at 6554 K → NLTE-valid.** Its PROVISIONAL status
  is the cross-instrument continuum zero-point (~0.18 dex) + the terminal [O I] 6300 leg,
  **not** an off-grid NLTE correction. No re-derivation owed.

### Coverage after the decision
- **Procyon 6554 K:** in-grid-once-loaded for Na/Mg/Si/Al/K/S (+ Ti/Cr/Mn where the family
  extends); Ca/Mn via the staged Amarsi `.grd`; **Ba clamp-flagged** (only true edge case);
  Sr deferred; **Fe already NLTE-live** (RYA-319). C/O 1D-NLTE in-grid.
- **τ Boo 6400 K:** inside every ≥6500 grid → fully in the published grid for all except the
  Sr 6000 metal-poor case.

### Execution status (Steps 3–5) — recipe recorded, physical re-synthesis pending
The self-consistent extension is a re-run of the existing machinery — `pipeline/pysme_nlte.py`
+ `scripts/resource_clamped_grids_rya409.py` — adding the Procyon node (6554/4.00/+0.03) to
each Amarsi element's node list, from the **same `.grd` binary** (gitignored/freed; refetch
from Zenodo 3982506, md5 per each `*_amarsi2020_v3.prov.json`, e.g. Mg `cdc4449e…`). No
cross-code mixing for the Amarsi set. Ca/Mn: synth from their staged Amarsi `.grd` and run
the ±0.05 solar cross-check vs the current MAFAGS-OS registered value before swapping. Ti/Cr:
re-scrape MPIA to 7000 K + solar overlap check. Ba: implement the bounded-clamp-with-σ. This
session delivered the recon gate + decision record + the CNO/O confirmations; the grid
re-synthesis is the bounded next execution step (needs the multi-GB `.grd` refetch + PySME).
