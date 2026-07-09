# RYA-530 — proposed CODEX_STATE_REGISTER section-edits

**Status: PROPOSED — for Ryan+Claude ratification** (register protocol: big/new rows decided together).
Base = whatever register version lands next (structure shown against v9; stack on top of the
534/535 register tip when those merge). Every value cites its ticket + merge commit — no value from memory.

---

## Edit A — `## Pipeline integrity` — add rows

```
| Single-source gf architecture | ONE authoritative log gf per physical line — `data/linelists/canonical_gf.csv` (145 887 lines) resolved at load by `pipeline/gf_resolver.py` for BOTH synth (branching-preserved HFS rescale) and EW paths; `pipeline/species.py` canonicalises every species/ion encoding (no silent mis-key); `scripts/check_stewardship.py` CI invariant fails the build on any UNTRACKED duplicated-and-divergent canonical value (gf/STAR_PARAMS/provenance/blend_flag) | RYA-345/350/353/354/355/358/368/408 (6eb2569/812d39c/bd881ad/aab6e2b/58a62e2) | SETTLED — wired in abundances_derive/synth/CNO; test-covered | a canonical value gains a second divergent source |
| HFS-resolved synthesis | flux-space synthesis fit over native GES hyperfine components (`cno_synthesis._fit_element`) — measures HFS-split lines the single-profile EW path SAT-culls; results feed problem_children/gold via RYA-519 | RYA-411/466/473/476 (eb5592d/723dec8/10a016c/77259c3) | SETTLED (method) — Mn measured this way is gold-tier (A=5.554 NLTE PASS, RYA-476 triplet-exact); Cu+V measured, held | new HFS element measured |
| EW-verification layer | per-line EW-integrity QA over the solar EW pool — flags BAD_FIT/ABUND_OUTLIER/COG_FLAG + literature LIT_DEVIATION cross-check; flags-only, `assert_no_ew_mutation` (validate-don't-tune); charter C I 5380 / Li 6707 / Eu 6645 | RYA-458 (e5321e6) | SETTLED — wired in abundances_derive (writes `{star}_ew_integrity.csv`) → aggregated by problem_children | — |
| Loader frame/OBJECT contract | executable anti-silent-substitution guard at every loader boundary — OBJECT-by-header (no glob), velocity-frame (no double-BERV), wave-scale (unit + vac/air, ×10/×10⁴ sanity gate); raise-never-default | RYA-481 + RYA-264 (b359189) | SETTLED — enforced; closes 6 documented incidents | a loader bypasses the contract |
| Multi-arm CNO arm-registry | per-star `STAR_ARMS` (region + diagnostics + loader + readiness); `run_cno`/`run_phase_a` resolve from it and loud-fail (`ArmNotWired`) instead of silently synthesising a non-solar star against Vesta sunlight | RYA-464 (9ebbfbf) | SETTLED — solar + Procyon declared (Procyon UVES ready; UV/IR deferred, cited) | new star needs its arms declared |
| Reflected-solar RV frame | asteroid-ephemeris conditioning (Vesta) — two-leg Horizons model empirically anchored to Fe-core bulk velocity, held-out verified, per-frame (no coadd); body-ID bug fixed | RYA-372 / RYA-394 (0874f6c/01d001b) | SETTLED (Vesta optical); reused by CRIRES+ IR (RYA-373) | — |
```

## Edit B — `## NLTE grids` — add 4 rows + correct 1

Add:
```
| Al I | Amarsi 2020 GALAH departure grid via PySME → `Al_Amarsi2020_PySME.csv`; atom Nordlander & Lind 2017 | subordinate 6696/6698 (not the +0.2 resonance 3961); solar median δ −0.022 | RYA-402 Fam-B (c10bd52) | SETTLED (selection + wiring) |
| K I  | Amarsi 2020 GALAH via PySME → `K_Amarsi2020_PySME.csv` | resonance 7665/7699, large negative (solar −0.27/−0.31); 7665 in O₂ A-band → clean 7699; measurement telluric-blocked (RYA-380) | RYA-402/462 (ea72a42) | SETTLED (selection + wiring); measurement GET-DATA-pending |
| S I  | Amarsi 2025 (A&A 703 A35) S departure via PySME → `S_Amarsi2025_PySME.csv` | optical 6748/6757 high-EP small negative (solar −0.016); NIR 9212/28/37 separate indicator | RYA-402 Fam-B (f1e35fa) | SETTLED (selection + wiring) |
| Sr II | WORKING = Bergemann 2012a/INSPECT → `Sr_Bergemann2012_INSPECT.csv` (applied). PRIMARY = Mashonkina 2022/INASAN NOT vendored (WAF-403, manual pull owed) | 4077 primary / 4215 Fe I-blended cross-check; [Fe/H] ceiling +0.0 → solar/55Cnc out-of-hull, loud | RYA-421/428/433 | PARTIAL — grid registered; primary pull + Sr II measurement owed (GET-DATA) |
```

Correct the existing **Na / Mg / Ba / Mn / Si** row — the Na/Mg/Si production source is now Amarsi-2020 PySME (RYA-410), not Lind-2011/Bergemann-MPIA:
```
| Na / Mg / Si | Amarsi 2020 via PySME (`*_Amarsi2020_PySME.csv`, RYA-410 re-source) — supersedes the legacy Lind-2011(Na)/Bergemann-MPIA(Mg/Si) CSVs (on disk, unwired) | feh axis is [Fe/H] RELATIVE | RYA-410/165/396 | SETTLED (selection) |
| Ba / Mn | Korotin 2015 (Ba, VizieR) · Bergemann MPIA (Mn) | — | RYA-165/396 | SETTLED (selection) |
```
Add a **note** below the table:
```
- Cu: departure inputs staged (`amarsi_galah/Cu_caliskan2024.*`) but grid NOT committed / NOT wired — blocker is measured-line quality (red_chi2 8–127), not the grid (RYA-466). V: genuine void, no public NLTE grid, HFS-LTE carry-forward (RYA-470/526).
```

## Edit C — `## 3D models & model atmospheres` — split applied-delta from RT engine

Add a row (keep line 112's STAGGER selection row):
```
| 1D→3D correction (applied delta) | `A(3D-NLTE)=A(1D-NLTE)+δ3d` from vendored `data/threed_grids/solar3d_metals_rya399.csv` (Si −0.01 Amarsi&Asplund 2017; Ti +0.06, Cr +0.03 Scott 2015); `3D_unavailable` flag if no node; wired via `nlte_corrections.apply_threed_corrections` | RYA-399 (a3c5796) | SETTLED (applied) — FINDING: 3D is NOT the lever for the Si/Ti/Cr residual (|δ|≤0.1); distinct from the still-absent 3D-RT *synthesis* engine (RYA-444) |
```
(Dead `CORRECTIONS_3D` scalar dict superseded by this grid path but not yet removed — RYA-474 owed.)

## Edit D — `## Instruments` — add UV conditioning + IR star-ID (siblings of registered rows)

```
| UV conditioning (no-telluric sibling of RYA-424 telluric) | STIS/COS FUV/NUV standing stage — vac→air @2000Å boundary, scattered-light + chromospheric-core masks, FUV synthesis-not-EW refusal, NLTE-grid-owed loud-flag; `analysis_ready` manifest (RYA-424 schema) + UV line-selection policy | RYA-426/190 | SETTLED (stage); UV arm production-deferred on RYA-359 |
| IR-native α Cen star-ID | RV-ephemeris orbit (Kervella 2016) predicts A/B heliocentric RV per epoch → PRIMARY discriminator for IR frames the optical classifier can't reach; `rv_bounds` rejects off-orbit frames | RYA-423/431 | BUILT — wired in star-ID scripts; RV>OBJECT authority ranking owed (RYA-495) |
```

## Edit E — `## Line lists` / instruments — loaders + converter + vac-air notes

Add rows/notes:
```
| Multi-instrument loaders | UVES (BERV/TOPOCENT + product guards + registry, wired), HST STIS/COS UV (wired ready=False pending RYA-359 grid), SPIRou; common `base_loader` SpectrumData contract | RYA-272/471 | SETTLED (UVES) / DEFERRED (HST UV) |
| Molecular linelist converter | ExoMol/MoLLIST `.states/.trans` → Turbospectrum `.bsyn`; CO round-trip hard gate; single path for all future molecular species | RYA-503 (PR#124) / RYA-360 | SETTLED |
| vac↔air SSOT | single Birch&Downs-1994 converter (`wavelength_util`), identity <2000Å; all loaders/conditioning + intake scripts routed through it (no second formula) | RYA-264/501 | SETTLED — no duplicate impl in-repo |
```

## Edit F — `## Scope / drift corrections` — dead-end guard note

```
| [O I] 6300 measured-O | continuum-renorm measurement = STOP (continuum-limited); production adopts Caffau-2015 8.73 literature, [O I] is a cross-check only. Do NOT re-attempt the renorm proof | RYA-447→455 (STOP proofs in pipeline/diagnostics/) |
```

---

## Placement summary
- **Pipeline integrity:** gf architecture, HFS synthesis, EW-verification, frame/OBJECT contract, multi-arm registry, reflected-solar RV (Edit A)
- **NLTE grids:** Al/K/S/Sr rows + Na/Mg/Si source fix + Cu/V note (Edit B)
- **3D models:** applied-delta row (Edit C)
- **Instruments:** UV conditioning + IR star-ID (Edit D)
- **Line lists:** loaders + molecular converter + vac-air SSOT (Edit E)
- **Scope/drift:** [O I] dead-end guard (Edit F)

Total: ~18 new rows / 1 corrected row / 2 notes. Header bump + one Changelog line (`RYA-530: capability-sweep reconciliation — registered N built-but-unregistered capabilities`).
