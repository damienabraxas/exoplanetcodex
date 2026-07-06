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
| **Li** | Amarsi 2020 GALAH grid **AVAILABLE** (Zenodo 3982506 `nlte_Li_scatt_pysme.tar.gz`, 378 MB; on Sirius amarsi_galah). | Acquire → PySME-derive CSV → wire; measurement is synthesis (Li I 6707 is CN-molecular-blended, RYA-103) — HFS/synth, not EW. | 103/458 |
| **Cu** | Grid **EXISTS** (Caliskan 2024 `nlte_Cu_caliskan_Oct2024_pysme.grd`, registered in `pysme_nlte`); no production CSV yet. | Route = HFS-resolved synthesis + RYA-402 b-factor NLTE (RYA-466); Cu is NOT an EW-pool element (over-saturated HFS core). Derive CSV when the synthesis measurement lands. | 402/466 |

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
