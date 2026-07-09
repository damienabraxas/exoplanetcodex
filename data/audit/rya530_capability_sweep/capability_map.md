# RYA-530 — 300/400-series Element & Capability Sweep (reconciled capability map)

**Goal:** catch already-built-but-forgotten capabilities before the Beta sign-off (RYA-527)
and the Procyon redo (RYA-404) — the "we nearly rebuilt RYA-402 as LTE-only" failure mode.
This is the **capability axis** complement to RYA-524's verdict-axis reconciliation.

**Evidence base (verified, not from memory):**
- Live Linear inventory: 484 issues pulled (RYA-40…536); 232 in scope (≥300): 128 Done / 81 Backlog / 13 Todo / 8 In Progress / 2 Canceled.
- Code truth: worktree at **origin/main @ `fe0c8b6`** (PR #139 / RYA-533).
- Register: **`CODEX_STATE_REGISTER.md` v9** on origin/main (v10–v15 = the RYA-534 rollout, pushed but **not merged**; v16 = RYA-535 Ti revert, also unmerged).
- Method per capability: confirm merge (`git log origin/main`), stub-check (NotImplementedError = not built), wired-check (imported/invoked), register grep.

**Discipline:** register = the *state* ledger (verdicts/selections/gates), not a code inventory. A module without a
register row is only a defect when it is a **substantial, re-inventable capability** — the bar RYA-402 set. Register
edits below are **proposed** section-edits for ratification (per the register's own "big rows decided together" rule),
not unilaterally committed.

---

## 1 · BUILT ✗ NOT-REGISTERED — the dangerous class (register updates owed)

Every item below is real code (non-stub), merged to origin/main, and wired, yet has **no register footprint**.
These are the RYA-402 pattern: whole capabilities a future planner could believe don't exist.

| # | Capability | Evidence (files · merge) | Wired | Tickets | Proposed register action |
|---|---|---|---|---|---|
| 1 | **Single-source gf architecture + species normalizer + stewardship CI** | `pipeline/gf_resolver.py` (229L), `pipeline/species.py` (217L), `data/linelists/canonical_gf.csv` (145 887 lines), `scripts/check_stewardship.py`, `pipeline/audit/gf_store_consistency.py` · `6eb2569`/`812d39c`/`bd881ad`/`aab6e2b`/`58a62e2` | YES — `abundances_derive` L64/207/599 + cno_synthesis/nlte_cno/nlte_corrections/crires_telluric; tested (`test_gf_resolver`, `test_species`, `test_data_stewardship`, `test_gf_store_consistency_rya368`) | 345/350/353/354/355/358/365/368/408/450 | **NEW row** (Pipeline integrity) |
| 2 | **NLTE grids: Al, K, S, Sr** present + wired, absent from grid table | `Al_Amarsi2020_PySME.csv` (`c10bd52`), `K_Amarsi2020_PySME.csv` (`ea72a42`), `S_Amarsi2025_PySME.csv` (`f1e35fa`), `Sr_Bergemann2012_INSPECT.csv` (`0722ba3`) · all in `NLTE_CORRECTION_ELEMENTS` | YES (constants.py 704/738/746/759) | 402 Fam-B / 462 / 421 | **4 NEW grid rows** (Sr = PARTIAL) |
| 3 | **HFS-resolved synthesis** (measures HFS-split lines the EW path SAT-culls) | `cno_synthesis._fit_element` over GES `..._hfs_iso`; drivers `measure_mn_hfs_synthesis_rya473.py` (`eb5592d`), `measure_cu_v_hfs_synthesis_rya466.py` (`723dec8`), `resolve_camn_stops_rya411.py`, `stage_amarsi_mn_grid_rya476.py` (`77259c3`) | YES (output → problem_children/gold via RYA-519); **Mn gold-tier A=5.554 PASS** measured this way | 411/466/473/476 | **NEW row** (Pipeline integrity / Methods) |
| 4 | **EW-verification layer** (per-line EW-integrity QA, flags-only) | `pipeline/ew_integrity.py` (284L, `e5321e6`); `assert_no_ew_mutation`; charter C I 5380 / Li 6707 / Eu 6645 | YES — `abundances_derive.py:2265` writes `{star}_ew_integrity.csv` → problem_children | 458 | **NEW row** (Pipeline integrity) |
| 5 | **Unsaturated line-selection recovery (Na/Mg)** | `measure_unsaturated_namg_rya465.py` + `add_namg_recovery_to_pool_rya465.py` (`0c72a17`) | YES — Na I 6154/6160 live in canonical `sol_ew_results_v1.csv`; Mg all over-ceiling (documented gap) | 465 | **NEW row** (Line lists / scope) |
| 6 | **Data-input / loader / conditioning stages** (8 modules, no register footprint beyond the instruments-observed table) | see §1a below | mixed (see §1a) | 272/471/481/264/426/190/464/372/394/423/431/503/501 | **grouped rows** (Pipeline integrity + Instruments + Line lists) |

### 1a · Loader / conditioning cluster (breakdown of row 6)

| Cap | Module · merge | Wired | Tickets | Proposed |
|---|---|---|---|---|
| Loader frame/OBJECT contract (anti-silent-substitution: OBJECT-by-header, no double-BERV, wave-scale gate) | `frame_object_contract.py` (588L, `b359189`) | YES — enforced at every loader boundary | 481/264 | **NEW row** (HIGH — executable guard, sibling of the registered RYA-424 telluric row) |
| Multi-instrument loaders (UVES, HST STIS/COS UV, SPIRou) | `loaders/uves_loader.py` (`0fd75f2`), `loaders/hst_uv_loader.py` (`a77b6fa`), `spirou_loader.py` | UVES YES; HST UV wired `ready=False` (gated on RYA-359 C/O grid) | 272/471 | **NEW row** |
| UV conditioning stage (no-telluric sibling of RYA-424) | `uv_conditioning.py` (444L, `a97e745`) + `uv_line_selection.py` (`d5b6fdc`) | PARTIAL — standing stage, no production UV run yet | 426/190 | **NEW row** (HIGH — parity with registered RYA-424) |
| Multi-arm CNO arm-registry (loud-fail vs silent Vesta substitution) | `cno_synthesis.py STAR_ARMS` (`9ebbfbf`) | YES — `run_cno`/`run_phase_a` resolve per-star; solar+Procyon declared | 464 | **NEW row** (HIGH — wiring-integrity class) |
| Reflected-solar RV frame (Vesta asteroid-ephemeris, 2-leg Horizons) | `reflected_solar_rv.py` (665L, `0874f6c`/`01d001b`) | YES — `reflected_solar` loader kind | 372/394 | **NEW row** (compact) |
| IR-native α Cen star-ID (RV-ephemeris orbit discriminator) | `acen_orbit.py` (99L, `b5fa4a6`/`e98531d`) | script-wired (future target) | 423/431 | **NEW row** (compact; RV>OBJECT authority owed → RYA-495) |
| Molecular ExoMol/MoLLIST→Turbospectrum converter (CO round-trip gated) | `scripts/molecular_linelist_convert.py` (406L, PR#124) | YES (acquisition tool) | 503/360 | **NEW note** under Line lists (lists registered, TOOL is not) |
| vac↔air SSOT (single Birch&Downs-1994; intake scripts routed) | `wavelength_util.py` (`ba0b92d`/`2a1a3e8`) | YES — no duplicate impl in-repo | 264/501 | **NEW note** (SSOT invariant) |

---

## 2 · PARTIAL / under-registered

| Capability | State | Register gap | Tickets |
|---|---|---|---|
| **1D→3D applied correction (delta)** | `threed_corrections.py` (233L, `a3c5796`) applies `A(3D-NLTE)=A(1D-NLTE)+δ3d` from `solar3d_metals_rya399.csv` (Si/Ti/Cr), wired via `nlte_corrections.apply_threed_corrections`. FINDING: *3D is not the lever* for the Si/Ti/Cr residual. | Register line 112 **conflates** this applied-delta with the still-absent 3D-**RT synthesis engine** (RYA-444). Separate them. | 399/447 |
| **Sr II** | Working grid = Bergemann-2012a/INSPECT (registered as row 2 above). **Primary = Mashonkina-2022/INASAN not vendored** (WAF-403, manual pull owed); **no Sr II line in any EW pool** → A(Sr II) not yet measured (RYA-428 reverted a premature LOCKED). Metal-poor MARCS ceilings → solar/55Cnc out-of-hull, loud. | Grid registrable now; primary + measurement flagged GET-DATA. | 421/428/430/433/422 |
| **Cu** | Departure *inputs* staged (`amarsi_galah/Cu_caliskan2024.*`) but **no grid CSV committed, not wired**. Register prose already says "PySME-validated but unregistered" — accurate. Blocker = **measured-line quality** (5 Cu lines red_chi2 8–127), not the grid. | Keep as a note, not an NLTE-grid row. | 402/466 |

---

## 3 · Register source-drift (values stale, not missing)

- **Na / Mg / Si NLTE row** (register line 101) credits **Lind-2011 / Bergemann-MPIA**, but the wired grids are **`*_Amarsi2020_PySME.csv`** (RYA-410). Old CSVs remain on disk, unwired. → **correct the row's source**.
- **Ti "~2× model-atom" prose** in the v15/v16 Engine-B row was **reverted by RYA-535** (same Bergemann-2011 atom; cause narrowed to MAFAGS-OS-vs-MARCS atmosphere, not model-atom). Only relevant once 534/535 merge; do not cite v15's Ti framing as settled. (RYA-535's territory, not 530's.)

---

## 4 · Correctly-absent dead-ends (do NOT rebuild)

- **Measured [O I] 6300 via continuum renorm** (447/448/449/452/453/454) = two **STOP proofs** in `pipeline/diagnostics/` (`continuum_renorm_proof.py` `4c079a5` "STOP", `continuum_local_proof.py` `219de80` "STOP, continuum-limited"), **not wired**. Production treats [O I] 6300 as a continuum-limited cross-check and adopts **Caffau-2015 8.73** (literature). Add a one-line register note so this proof is not re-attempted.

---

## 5 · Owed-work dedup (already-built-but-still-open) + superseded

**Merged to origin/main but Linear pill still Backlog — recommend verify + close (NOT owed build):**
| Ticket | State | Merge evidence | Note |
|---|---|---|---|
| RYA-336 | Backlog | `9ccf064` scale-aware solar Fe gate | register already flags "Backlog pill is stale" |
| RYA-330 | Backlog | `5eb73db`/`04ec183` solar Fe I anchor cut + gate double-add | landed via RYA-331 step 5 |
| RYA-334 | Backlog | `5aecfd3` A_X_nlte double-add range guard | |
| RYA-317 | Backlog | `916f96a` Fe NLTE Teff-edge clamp (0-lines fix) | |

**Superseded / delivered by later work — recommend close + link:**
| Ticket | State | Superseded by | Evidence |
|---|---|---|---|
| RYA-361 (scope Gerber in-synthesis TS-NLTE) | Backlog | **RYA-533** (built the deck) | `gerber_ts/` deck + `ts_gerber_na_gate_rya533.py` on main; register v9 records "RYA-361 scoped it, never executed" → now built |
| RYA-362 (O I 777 NLTE post-hoc vs in-synth cross-check) | Backlog | **RYA-534** rollout | v10 changelog: "doubles as the RYA-362 post-hoc-vs-in-synthesis O cross-check (agree to 0.03 dex)" (on the unmerged 534 branch) |

**Genuine owed cleanup (confirmed still real):**
| Ticket | State | Confirmation |
|---|---|---|
| RYA-474 (remove dead `CORRECTIONS_3D` scalar dict, superseded by grid-based `threed_corrections`) | Backlog | dict **still present** `config/constants.py:1039`, referenced "legacy" in `uv_line_selection.py:24` → real owed removal |

**Genuine remaining owed-build (audit found NOT built — stays owed):** Sr II primary pull + measurement (421/428/430/433/422), Cu grid commit + line-quality fix (466), V ionization-anchor (470), FUV C I capability (487/348-P3), telluric infra installs (375/438/380/437/494/480), CI gate (313/314/436), UV production arm (359). The audit does **not** move these — it only removes the already-built items above from the owed list.

---

## 6 · Adequately BUILT & REGISTERED (spot-confirmed, no action)

Fork/runtime (514) · authoritative channel (521) · wiring-integrity/problem_children (463/519) · two-engine floor +
Engine-A/B substrate (525/529/531/533) · NLTE grids Fe/Ca/Ti/Cr/Na/Mg/Ba/Mn/Si/C/O/N · telluric routing (424) ·
N wiring (526) · solar gold v2 (522) · Sirius/compute-stack/venv-isolation (511/517/526) · NIRPS from-raw + α Cen ADP (498/500/494) · line lists OH/NH/CH/CO (503).

---

## Deliverable index
- This map: `data/audit/rya530_capability_sweep/capability_map.md`
- Proposed register section-edits: `data/audit/rya530_capability_sweep/proposed_register_edits.md`
