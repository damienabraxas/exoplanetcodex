# RYA-526 — Grid coverage program: every element on a synthesis grid

**Deliverable:** every one of the 26 canonical `TARGET_ELEMENTS` (config/constants.py:141)
has EITHER a wired NLTE synthesis grid OR a documented acquisition/build task with an owner.
This is the living coverage ledger (companion to `data/curation/nlte_grid_availability.csv`).

**Principle (RYA-525 two-engine floor):** "no grid → EW-1D cop-out" is dead — the *absence*
of a grid is a task, not a verdict. Note: **1D-LTE synthesis covers EVERY element** (the
Turbospectrum/iSpec linelist path); the tables below track the **NLTE departure grid** layer
on top of that. An element with "no NLTE grid" is still synthesizable in LTE — the task is the
NLTE enhancement, unless it is a genuine `NLTE_VOID` (V).

**Method:** ground-truth, no memory — dispositions verified against `NLTE_CORRECTION_ELEMENTS`
(config), `pipeline/pysme_nlte.py`, `pipeline/problem_children.py`, the Amarsi-2020 GALAH
Zenodo record `3982506` (element set verified live via the Zenodo API), and the on-disk
grid inventory on Sirius (`/mnt/codex-data/grids/nlte/amarsi_galah/`).

---

## A. WIRED — NLTE grid in production (15)

| El | Ion | Grid file | Reference | Ticket |
|----|-----|-----------|-----------|--------|
| Fe | I/II | Fe_Bergemann_MPIA.csv | Bergemann MPIA MAFAGS-OS 1D NLTE | 319/407 |
| C  | — | amarsi2019_cno/ (3D-NLTE synthesis) | Amarsi 2019 CH G-band + C I | — |
| O  | — | amarsi2019_cno/ (3D-NLTE synthesis) | Amarsi 2019 O I 777 + [O I] 6300 | — |
| Na | I | Na_Amarsi2020_PySME.csv | Amarsi 2020 (A&A 642 A62) via PySME; Lind 2011 | 410 |
| Mg | I | Mg_Amarsi2020_PySME.csv | Amarsi 2020 via PySME | 410 |
| Al | I | Al_Amarsi2020_PySME.csv | Amarsi 2020 via PySME; Nordlander & Lind 2017 | 402 |
| Si | I | Si_Amarsi2020_PySME.csv (+ 3D solar node) | Amarsi 2020 via PySME; Amarsi & Asplund 2017 (3D) | 410/399 |
| S  | I | S_Amarsi2025_PySME.csv | Amarsi 2025 (A&A 703 A35) via PySME | 402 |
| K  | I | K_Amarsi2020_PySME.csv | Amarsi 2020 via PySME; K I 7665/7699 | 462 |
| Ca | I | Ca_Mashonkina2017.csv | Mashonkina 2017 (A&A 606 A147) MPIA | 411/413 |
| Ti | I | Ti_Bergemann2011_MPIA.csv (+ 3D) | Bergemann 2011 (MNRAS 413 2184) MPIA | — |
| Cr | I | Cr_Bergemann2010_MPIA.csv (+ 3D) | Bergemann & Cescutti 2010 (A&A 522 A9) MPIA | — |
| Mn | I | Mn_Bergemann_MPIA.csv (+ Amarsi triplet grid, RYA-476) | Bergemann MPIA; Amarsi Mn (55Cnc triplet) | 411/476 |
| Ba | II | Ba_Korotin2015.csv | Korotin 2015 (A&A 581 A70) | 279/165 |
| Sr | II | Sr_Bergemann2012_INSPECT.csv | Bergemann 2012a via INSPECT (working); Mashonkina 2022/INASAN primary owed | 421/433 |

Sr carries a **primary-grid acquisition owed** (Mashonkina 2022/INASAN — INASAN nLTE.cgi
WAF-blocks programmatic clients, RYA-433; manual pull owed, owner Ryan). The applied
Bergemann-2012 working grid is production-valid; near-LTE at our metallicities.

---

## B. WIRED THIS SESSION — Part A of the ticket ✅ DONE

**N — wired into production (RYA-526).** Derived `N_Amarsi2020_PySME.csv` on Sirius via PySME
from the Amarsi-2020 N grid (`nlte_N_scatt_pysme.grd`, Zenodo 3982506, md5-verified), over the
standard 11-node box, and added N to `NLTE_CORRECTION_ELEMENTS` (ion 1, flag
`NLTE_Amarsi2020_PySME_1D`). Validation (validate-don't-tune): the **solar node reproduces the
committed RYA-369 load-test EXACTLY** — 7468 −0.0115 / 8216 −0.0145 / 8683 −0.0154, max |Δ| =
0.0000 dex. Production `nlte_corrections._mpia_element_delta('N', …)` returns those solar
deltas; N I red is near-LTE at the Sun (sane-negative, RYA-339). **4 weak-line COG rails** at
cool metal-rich nodes (7468 ×3, 8216 ×1, Teff ≤ 5172) were EXCLUDED → they fall out-of-hull and
loud-flag in production (RYA-409), never injected — consistent with N I red being a *warm*-star
indicator (RYA-369). Grid derivation ran on Sirius (`venv_pysme`, numpy 2.2.6); the `.grd` stays
on Sirius only (never the Mac).

---

## C. ACQUIRE / WIRE — grid available or route exists (owner: Ryan)

| El | NLTE grid situation | Route | Ticket |
|----|---------------------|-------|--------|
| **Li** | Grid staged on Sirius (`nlte_Li_scatt_pysme.grd` 7.6 GB). **RYA-540: PySME derivation ATTEMPTED** (`scripts/wire_li_cu_rya540.py`) → solar δ = **−0.030**, does **NOT** reproduce the Lind-2009 small-positive anchor (sign discrepancy; the grid J-label resolves oddly as J0→0). **NOT wired** (validate-don't-tune STOP). | Li I 6707 resonance+CN-blend derivation needs the dedicated ticket; not EW. | 103/458 |
| **Cu** | **RYA-540: DONE (grid).** Grid absent on disk but the prov pins a Zenodo DOI (15062813 v6) — re-acquired (md5 d6fce44e VERIFIED), PySME-derived `Cu_Caliskan2024_PySME.csv` (Cu I 5782 solar δ **+0.001** = Shi-2014), **REGISTERED** (PENDING-OK: registered grid, verdict stays GET-DATA). | Grid half complete + registered. Promotion to LOCKED still owed to **measured-line quality** (RYA-395: 5105/5218/5782). | 395/466 |

---

## D. NO NLTE GRID — but not the blocker; LTE-synthesis covers (owner: Ryan)

These have **no NLTE departure grid**, but NLTE is *not* the limiting factor — the blocker is
gf-scale / data-coverage / continuum, and 1D-LTE synthesis already covers the element. Each is a
documented data/curation task, not a grid void.

| El | Class | Blocker (not NLTE) | Ticket |
|----|-------|--------------------|--------|
| Ni | BAD_GF | gf-scale residual floor; differential-survey curation owed | 398/161 |
| P  | DATA_GAP | FUV needs HST/STIS; near-IR 10581/10596 ground-reachable, gf-limited | 119/460 |
| Co | CONTINUUM_LIMITED | blue-edge SNR-limited; extract cleaner red Co I lines + HFS | 460 |
| Sc | HFS_SUMMING | blue-edge HFS single line (4246), LOW_CONFIDENCE; HFS synth + cleaner line | 460 |
| Y  | DATA_GAP | Y II (dominant ion) absent from pool; acquire blue/UV Y II arm | 458 |
| Eu | HFS_SUMMING | Eu II 6645 HFS; HFS-resolved synthesis | 102/458 |
| Zr | (none needed) | Zr II is the **majority ion → LTE-robust**; route = synthesis, **no NLTE grid required** (measurable-owed, not a void) | 279/458 |

---

## E. BUILD — genuine NLTE void (owner: Ryan)

| El | Class | Situation | Route | Ticket |
|----|-------|-----------|-------|--------|
| **V** | NLTE_VOID | No Amarsi/GALAH grid, no usable neutral V I model atom (verified: not in Zenodo 3982506; `problem_children` NLTE_VOID). | **Build/bridge** (RYA-363 "can the Codex build it" is now mandatory) OR the interim **V II ionization-anchor** route (RYA-470 — the "Strontium move": measure around the unsolved model atom, V I−V II gap = empirical void depth). | 470/363/404 |

---

## Summary

| Disposition | Count | Elements |
|---|---|---|
| A. Wired (production NLTE) | 15 | Fe C O Na Mg Al Si S K Ca Ti Cr Mn Ba Sr |
| B. Wired THIS session (RYA-526) ✅ | 1 | N |
| C. Acquire/wire (grid exists) | 2 | Li Cu |
| D. No NLTE grid needed / LTE-synth covers | 7 | Ni P Co Sc Y Eu Zr |
| E. Build (NLTE void) | 1 | V |
| **Total** | **26** | all `TARGET_ELEMENTS` (now **16 wired**, 10 documented tasks) |

**Every element is accounted for:** wired, or a documented acquisition/build task with an owner
and governing ticket. The Amarsi-2020 GALAH Zenodo release (3982506, verified live) supplies
13 elements (Al Ba C Ca H K Li Mg Mn N Na O Si); of the canonical set that leaves Fe/Ti/Cr
(MPIA), S (Amarsi 2025), Sr (Bergemann/Mashonkina), Cu (Caliskan 2024) sourced elsewhere, and
Ni/P/Co/Sc/Y/Eu/Zr/V with no departure grid (V = the only true void).

## Ticket-premise corrections (confirm-don't-assume)

1. **Mn "8.4 GB gitignored grid" — NOT FOUND.** The largest file in the entire project tree is
   `ispec/input.tar.gz` (unrelated). Mn's production grid is the **committed 8.4 KB**
   `data/nlte_grids/Mn_Bergemann_MPIA.csv` — the PASS depends on a file that IS in the repo,
   not a missing 8.4 GB one. The only gitignored Mn artifact is the offline `.grd`
   (re-derivation input, on Sirius), exactly like every other PySME element. Stewardship guard
   added anyway (presence check on the committed CSV). The "8.4 GB" appears to be a KB/GB
   transcription error in the ticket.
2. **N "grid is Done" — half true.** N is registered + load-tested (RYA-369) but the production
   CSV was never derived and N was never added to `NLTE_CORRECTION_ELEMENTS`; the `.grd` is
   Amarsi-2020 (Zenodo 3982506), staged on Sirius. Wiring = the derivation in §B.

---

## Two-engine coverage refresh — RYA-526 (2026-07-10, additive)

**Why:** the tables above are the **Engine-A (1D-NLTE departure grid)** coverage map and predate
the RYA-534 Engine-B (TS-Gerber in-synthesis NLTE) rollout. RYA-525's loud-fail guard reads the
coverage ledger as its **pre-declared exception list** — it must see *both* engines per element or
it cannot tell a genuine gap ("acquire/build a grid") from a correctly-finished element ("LTE by
design"). This section adds the **Engine-B axis**; it does **not** alter the Engine-A ledger above
(the RYA-543 registry↔disk anti-drift test reads that unchanged).

**Method (ground truth, no memory):** disk-verified on Sirius, 2026-07-10 —
`/srv/codex/grids/nlte/gerber_ts` (Engine-B atoms/aux/provenance/md5 + `_cache_index.json`) and
`/srv/codex/grids/nlte/amarsi_galah/*.grd` (Engine-A PySME inputs), per the RYA-540 census
(`data/audit/rya540_disk_layout/`). Machine-readable companion:
`data/curation/nlte_two_engine_coverage.csv`. Gerber-2023 (A&A 669 A43) covers 14 target elements
(O Na Mg Al Si Ca Ti Mn Fe Co Ni Sr Ba Y); the `.bin` grids are freed-but-md5-pinned and
re-acquired on demand by `scripts/grid_cache.py` (RYA-540) — provisioned+deck-validated (RYA-534)
counts as **Engine-B wired**.

**State vocabulary (per engine):** `wired` (grid present + validated + pipeline-reachable) ·
`task` (grid genuinely absent, an acquisition/build task with an owner) · `LTE-only-by-design`
(no NLTE grid exists **and none is needed** — near-LTE or majority-ion species; a *finished* state,
NOT owed work — this is what stops 525 raising forever on a correctly-done element).

| El | Engine-A (1D-NLTE) | Engine-B (TS-Gerber) | Disposition | Owner · ticket |
|----|----|----|----|----|
| Fe | wired | task | wired-one | Ryan · 534/361 |
| C  | wired | task | wired-one | Ryan · 363 |
| O  | wired | wired | **wired-both** | — · 534 |
| Mg | wired | wired | **wired-both** | — · 534 |
| Si | wired | wired | **wired-both** | — · 534 |
| Ca | wired | wired | **wired-both** | — · 534 |
| Ti | wired | wired (CHECK) | **wired-both** | Ryan · 535 |
| Ni | task  | wired | wired-one | Ryan · 534 |
| Na | wired | wired | **wired-both** | — · 533/534 |
| P  | LTE-by-design | LTE-by-design | **LTE-only-by-design** | — · 460 |
| S  | wired | task | wired-one | Ryan · 361 |
| N  | wired | task | wired-one | Ryan · 369/363 |
| Co | task  | wired | wired-one | Ryan · 534 |
| Cr | wired | task | wired-one | Ryan · 361 |
| Al | wired | task | wired-one | Ryan · 534/361 |
| K  | wired | task | wired-one | Ryan · 361 |
| Ba | wired | wired | **wired-both** | — · 534 |
| Y  | task  | task | acquire-task | Ryan · 458 |
| V  | task | task | **build-task** (NLTE_VOID) | Ryan · 470/363 |
| Cu | task  | task | acquire-task | Ryan · 466 |
| Mn | wired | wired | **wired-both** | — · 534 |
| Sc | LTE-by-design | LTE-by-design | **LTE-only-by-design** | — · 460 |
| Li | task  | task | acquire-task | Ryan · 103/458 |
| Eu | LTE-by-design | LTE-by-design | **LTE-only-by-design** | — · 458 |
| Zr | LTE-by-design | LTE-by-design | **LTE-only-by-design** | — · 526 |
| Sr | wired | wired | **wired-both** | Ryan · 421/534 |

**Per-state counts (26 TARGET_ELEMENTS, each resolves to exactly one):**
- **wired-both = 9** — O Na Mg Si Ca Ti\* Mn Ba Sr  (\*Ti = CHECK, RYA-535: run both, record the spread, excluded from the reported value)
- **wired-one = 9** — A-wired/B-task: Fe C N S K Cr Al · B-wired/A-task: Co Ni
- **acquire-task = 3** — Li Cu Y  (grid exists/gettable; not yet derived+wired)
- **build-task = 1** — V  (the only true NLTE_VOID)
- **LTE-only-by-design = 4** — P Sc Eu Zr  (finished in LTE; NLTE grid neither exists nor needed)

**For RYA-525's guard:** raise ONLY where `disposition ∈ {acquire-task, build-task}` **and** the
grid is genuinely absent — never on `wired-*` (grid present) or `LTE-only-by-design` (finished).
`wired-one` elements run the one wired engine now and carry the other engine as a documented,
owned task (fed to RYA-540 for the Gerber-available set: Fe, Al, Y; and to RYA-363 for the
Gerber-absent set: C, N, S, K, Cr).
