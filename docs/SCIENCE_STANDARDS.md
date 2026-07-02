# SCIENCE STANDARDS — authoritative decision records

Single source of truth for cross-cutting science-method decisions. New file (RYA-505);
flag for the RYA-179 docs-sync revision. Add sections as decisions are recorded; each
section cites its evidence and the ticket that set it.

---

## Multi-model NLTE at F-star Teff (Procyon 6554 K, τ Boo 6400 K) — RYA-505

**Decision date:** 2026-07-02 · **Evidence:** `scripts/nlte_fstar_ceiling_rya505.py`
→ `data/results/nlte_fstar_ceiling_rya505.csv`, `data/results/nlte_multimodel_rya505.csv`.
**Scope:** hot-benchmark fidelity upgrade; NOT on the α Cen / 55 Cnc critical path (those
are 5200–5800 K, already in-grid). The 55 Cnc metal-rich [Fe/H] edge is a separate axis
(RYA-410). Does not block RYA-349.

### Standing principle — multiple NLTE models == multiple instruments
We use every model that is useful, always — the same architecture as multi-instrument
abundances. The two independent NLTE model families we hold — **Amarsi/Balder/MULTI on
MARCS(/Stagger)** and **Bergemann/DETAIL/SIU on MAFAGS-OS** — are two independent
**measurements** of the NLTE correction. Policy:
1. **Run every model family that covers** an (element, Teff, logg, [Fe/H]) point. Emit a
   correction **per model** — the model analog of a per-instrument product.
2. **≥2 models cover** → report each, plus the **inter-model spread** as a measured
   model-systematic into the σ budget (RYA-282). Agreement within tolerance → a **combined
   / consolidated** value + agreement badge (never the word "blend"). Disagreement beyond
   tolerance → **flag + adjudicate**, never silently average, never silently pick one.
3. **Only 1 model covers** → use it, **flagged single-model** (no cross-check available
   there — a known, logged limitation, not a silent choice).
4. **No model covers** → bounded Teff clamp at the ≤54 K overshoot ONLY, monotonic-in-Teff
   elements only, systematic propagated into σ (RYA-282 / Fe RYA-319 precedent). Never
   extrapolate a 3D-NLTE correction past its grid; never a silent clamp.

**Atmosphere-baseline caveat (so the spread is meaningful, not confounded):** a correction
Δ = A(NLTE) − A(LTE) is defined against a baseline on its *own* model atmosphere. Amarsi
(MARCS/Stagger) vs Bergemann (MAFAGS-OS) differ in atmosphere as well as NLTE method, so the
**raw inter-model spread is an UPPER BOUND** on the model systematic — it carries both the
NLTE-method term and the atmosphere-baseline term. Report it as an upper bound; characterize
the two terms where feasible; never hide the atmosphere term inside "NLTE disagreement".

### Headline recon (the RYA-349 "6500 K wall")
The wall is almost entirely an **under-loading artifact, not a grid limit**: our on-disk
non-Fe grids are the RYA-410 benchmark-node **subsets**, not the full published grids. The
published grids reach F-star Teff (Amarsi GALAH 2500–8000 K; Bergemann/MPIA to 7000 K). Only
**Ba** (Korotin2015, real limit 4000–6500) is a genuine no-higher-grid case.

### Per-element, per-model coverage map (Procyon 6554 K)
On-disk ceilings read from the grid CSVs; published ceilings from cited provenance.

| element | Amarsi/MARCS family | Bergemann/MAFAGS-OS family | other | #models @ Procyon | posture |
|---|---|---|---|---:|---|
| Na | Amarsi2020, on-disk 6200 / **pub 8000** | — | Lind2011 INSPECT (MARCS) on-disk 6500 | **2** | run both; **solar spread 0.023 dex** |
| Mg | Amarsi2020, 6200 / **8000** | Bergemann MPIA, 6500 / **7000** | — | **2** | run both; **solar spread 0.035 dex** |
| Si | Amarsi2020, 6200 / **8000** | Bergemann MPIA, 6500 / **7000** | — | **2** | run both; **solar spread 0.007 dex** (near-LTE) |
| Ca | Amarsi Ca `.grd` staged (Family-A) | Mashonkina2017 (registered), 6500 | — | **2** (Amarsi via synth) | run both; overlap check owed |
| Mn | Amarsi Mn `.grd` staged | Bergemann MPIA (registered), 6500 / **7000** | — | **2** (Amarsi via synth) | run both; overlap check owed |
| Al | Amarsi2020, 6200 / **8000** | — | — | **1** | **single-model** (Amarsi), flagged |
| K  | Amarsi2020, 6200 / **8000** | — | — | **1** | single-model (Amarsi), flagged |
| S  | Amarsi2025, 6200 / **8000** | — | — | **1** | single-model (Amarsi), flagged |
| Ti | — | Bergemann2011 MPIA, 6500 / **7000** | — | **1** | single-model (MPIA), flagged; near-LTE |
| Cr | — | Bergemann2010 MPIA, 6500 / **7000** | — | **1** | single-model (MPIA), flagged |
| **Ba** | — | — | Korotin2015 (MULTI/MARCS), **6500 REAL LIMIT** | **1** | single-model + **bounded clamp** (54 K, monotonic, σ) |
| Sr | — | — | Bergemann2012 INSPECT (metal-poor), 6000 | **1** | defer (off-critical; INASAN primary pending RYA-433) |

Priority for chasing hot-Teff coverage (large F-star corrections first): **Ca, Ba, Na, Mn**,
then Mg; Si/Al/Ti near-LTE, low priority.

### Model cross-check at the solar overlap (Step 4 — measured now, from on-disk grids)
Where two families are already on disk, the solar (5772/4.44/0.0) inter-model spread — the
model-systematic meter — is:
- **Mg:** Amarsi −0.0216 vs Bergemann-MPIA +0.0132 → **spread 0.035 dex** (MARCS vs MAFAGS-OS
  → upper bound). Within 0.05 → **AGREE**.
- **Si:** Amarsi −0.0109 vs Bergemann-MPIA −0.0040 → **spread 0.007 dex**. AGREE (near-LTE).
- **Na:** Amarsi −0.1293 vs Lind2011-INSPECT −0.1068 → **spread 0.023 dex** (both MARCS-family
  → cleaner). AGREE.
All three pass a 0.05 dex tolerance at the Sun → the two families are consistent measurements
where they overlap; their spread is the σ contribution to carry to hot Teff. (Ca/Mn overlap
checks are owed once the staged Amarsi `.grd`s are synthesised.)

### CNO (Step 1) and banked Procyon O (Step 2) — already correct, confirmed
- **C/O NLTE:** `pipeline/nlte_cno.py` splits legs at `TEFF_3D_CEILING = 6500 K`: 3D leg
  (tables 2/3) ≤ 6500 K; **1D-NLTE leg (tables 5/6) above** → **Procyon 6554 K uses the
  1D-NLTE leg, Teff 4000–8000 K (17 nodes) → in-grid.** No change (RYA-359). Self-consistent
  3D-NLTE refinement is future v2 (RYA-444/445), matching the IR v1/v2 boundary (RYA-504).
- **RYA-483 banked Procyon O:** primary indicator `OI_777_1D_NLTE`, A(O)=8.82, [O/H] +0.085 —
  **used the 1D-NLTE correction, in-grid at 6554 K → NLTE-valid.** Its PROVISIONAL status is
  the cross-instrument continuum zero-point (~0.18 dex) + the terminal [O I] 6300 leg, **not**
  an off-grid NLTE correction. No re-derivation owed.

### Coverage after the decision
- **Procyon 6554 K:** 2-model (run both + spread) for Na/Mg/Si, and Ca/Mn once the staged
  Amarsi `.grd`s are synthesised; single-model-flagged for Al/K/S (Amarsi) and Ti/Cr (MPIA);
  **Ba single-model + bounded clamp** (the only no-second-model, real-limit case); Sr
  deferred; Fe already NLTE-live (RYA-319); C/O 1D-NLTE in-grid.
- **τ Boo 6400 K:** inside every ≥6500 K grid → in the published grid for all except the Sr
  6000 K metal-poor case.

### Execution status — recon + principle recorded; per-model wiring is the next step
Delivered this session: the standing principle, the per-element/per-model coverage map, the
CNO/O confirmations, and the solar-overlap inter-model spreads (Mg 0.035 / Si 0.007 /
Na 0.023 dex) computed from the on-disk grids. **Pending execution (per-model wiring at
Procyon):** (1) load/synthesise the Amarsi Procyon node (6554/4.00/+0.03) from the same
`.grd` (Zenodo 3982506, md5 per `*_amarsi2020_v3.prov.json`) for Na/Mg/Si/Al/K/S; (2)
synthesise Ca/Mn Amarsi from the staged `.grd` and run their solar overlap check; (3) load
the Bergemann/MPIA node to 7000 K for Mg/Si/Ti/Cr/Mn so both families report at Procyon;
(4) emit the per-model correction + inter-model spread per element into the RYA-282 σ layer;
(5) Ba bounded-clamp-with-σ. The registry must move from one grid per element to a
**per-model set** so both corrections are reported, never pre-combined.
