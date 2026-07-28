# RYA-401 — Beast round-2 NLTE-grid vendoring (Al/K/S/Cu/V/Sr): the dead-ends, recorded

RYA-400 ("the Beast") routed **6 elements to GET-GRID** (Al, K, S, Cu, V, Sr). This is
round-2 of RYA-396's vendoring (which did Na/Mg/Ba/Mn/Si via the INSPECT/MPIA scrapers).

**Outcome: 0 of 6 vendored — all six are documented gaps, each for a distinct reason.**
The GET-GRID label in round-1 was an *optimistic* audit call; round-2 probed the actual
sources live and the optimism did not survive contact. No grid was fabricated (the
ticket's explicit rule). Each element is flipped out of GET-GRID in
`config/physics_regime_rya400.yaml` with a `rya401:` block recording the probe, the gap
reason, and the recovery path; the RYA-400 audit script still **PASS**es.

This page is the "don't re-walk the dead ends" record (the ticket's words).

## What the sources actually host (probed live, 2026-06-21)

| Source | Hosts | Probe |
|--------|-------|-------|
| **MPIA SpectrumTools** (`nlte.mpia.de`, `build_nlte_grids_mpia.py`) | Ca, Co, Cr, Fe, H, Mg, Mn, O, Si, Ti | `…lines[]` form params — **none of Al/K/S/Cu/V/Sr** |
| **INSPECT** (`inspect-stars.com`, `fetch_inspect_nlte.py`) | Na, Mg (wired) + **Sr** (Sr II) | `A_from_e?element_name=X`: **Sr → form OK**; Al/K/S/Cu/V → **HTTP 500** |
| **CDS/VizieR** (Korotin-style deposits) | Ba (J/A+A/581/A70, vendored RYA-165) | no S/Cu correction-grid deposit found |
| **Amarsi 2020 GALAH** (A&A 642, A62) | Li/C/N/O/Na/Mg/Al/Si/K/Ca/Mn/Ba | **departure coefficients** (not delta-CSV) → needs RT line-formation conversion |

So the only scraper-reachable grid among the six is **Sr II via INSPECT** — and that one
hits an ion mismatch (below).

## Per-element verdict (round-2)

| El | Was | Now | Reason (cited) | Recovery path |
|----|-----|-----|----------------|---------------|
| **Sr** | GET-GRID | **GET-DATA** | INSPECT hosts **Sr II** (4077/4215); our solar set measures only **Sr I 6617** → ion mismatch. A Sr II grid won't apply to Sr I. `[Bergemann2012]` | measure Sr II 4077/4215 (blue/saturated) → INSPECT Sr II scrape is then quick |
| **Al** | GET-GRID | **HARD-carry-forward** | INSPECT 500; not on MPIA. Only Amarsi2020 **departure coefficients** (RT conversion, not a scrape). `[Nordlander&Lind2017; Amarsi2020]` | Amarsi2020 departure-coeff conversion (RYA-399-class) or carry-forward |
| **K** | GET-GRID | **HARD-carry-forward** | Amarsi2020 departure-coeffs **+ 0 measured lines + resonance in the O₂ telluric band**. `[Amarsi2020; Asplund2021]` | RYA-380 (telluric) + measurement + Amarsi2020 conversion — multi-blocker |
| **S** | GET-GRID | **HARD-carry-forward** | **No public machine-readable grid** (Korotin 2009 author/compute; tables only, 2≤logg≤4 partial). `[Korotin2009; Takeda2005]` | indicator decided (below); author-contact/digitise or carry as LTE-flagged |
| **Cu** | GET-GRID | **HARD-carry-forward** | Literature explicit: **no available NLTE grid for neutral Cu** (per-star corrections only). `[Korotin2018; Andrievsky2018]` | author-contact/compute or carry-forward (HFS-LTE) |
| **V** | GET-GRID (GAP) | **HARD-carry-forward** | V I NLTE **sparsely studied**, no public grid; HFS likely dominates over NLTE at solar. `[Asplund2021]` | carry-forward (HFS-LTE); HFS first |

## S — the indicator decision (the Beast flagged this matters more than the grid)

**Chosen indicator: the optical high-EP S I multiplet 6743/6757** (we measure 6748.7,
6757.0, 6757.1), **not** the strong NIR triplet 9212/9228/9237. Rationale (per the Beast's
`indicators` column + RYA-393 strategy logic):

- **Optical 6743/6757:** small solar NLTE (≈ −0.1 dex), telluric-free, already measured.
- **NIR 9212/9228/9237:** strong but **−1.1 dex NLTE**, Takeda 2005 vs Korotin 2009 disagree
  **~0.4 dex** at 9213, and sits among NIR tellurics — high-risk for a *solar* number.

Consequence: because the chosen optical indicator's solar NLTE is small, S can run
**LTE-flagged** at solar (a ~−0.1 dex known systematic, documented) until the Korotin
optical grid is obtained — rather than chase the −1.1 dex NIR correction with no grid.

## K ↔ RYA-380 (cross-reference, not re-solved here)

The K I 7665/7699 resonance doublet sits in the O₂ telluric A/B band. Telluric usability is
**RYA-380's** job (the per-night-GDAS red-optical/IR recipe) — not re-solved in RYA-401. K
stays blocked on three independent things (telluric + measurement + the Amarsi2020
departure-coeff conversion), so it carries forward.

## Net for RYA-371 Phase C

These 6 do not gain NLTE this round. They are **carried forward** honestly (validate-don't-tune):
Phase C reports them with their regime status (LTE-flagged + the documented gap), and a
benchmark in a different regime across the multi-star arc may unveil what the Sun won't.
The recovery paths above are the concrete follow-on work (Amarsi2020 departure-coeff
conversion for Al/K; Sr II / Sc / Co measurements; Korotin author-contact for S/Cu).

## Citations
Asplund/Amarsi/Grevesse 2021 (A&A 653, A141); Amarsi et al. 2020 (A&A 642, A62);
Nordlander & Lind 2017 (A&A 607, A75); Korotin 2009 (Astron. Rep. 53, 651); Takeda et al.
2005; Korotin et al. 2018 (MNRAS 480, 965); Andrievsky et al. 2018; Bergemann et al. 2012
(A&A 546, A90). Source probes: MPIA `nlte.mpia.de`, INSPECT `inspect-stars.com` (live, RYA-401).
