# Codex State Register

**The mutable current-truth ledger. Read this FIRST.**

The Linear tickets are an *immutable journal* — every step, dead end, and reversal, recorded and never rewritten. This register is the *mutable ledger* — the single authoritative statement of **what is true right now**, rewritten as things settle, regress, or supersede. When "current state of X" and "what the journal says happened to X" diverge, **this file wins for state**; the journal wins for history.

Built to be read fast by both Ryan/Claude.ai **and** local models on Sirius (Qwen etc.) — keep it clean-table, one fact per row, machine-parseable.

## How to read this — NATIVE vs MIRROR rows
- **NATIVE** — this register *is* the source of truth: verdicts, statuses, gate states, scope decisions, grid/model/instrument selections, `reopen-only-if` triggers. Hand-maintained here.
- **MIRROR** — reflects another single source (stellar params ← `config/constants.py` → `STAR_PARAMS`/`stars.yaml`; abundances ← results tables). Shown for fast access, but **regenerated from source, never hand-trusted**. A mirror value that disagrees with its source is a defect. Mirror tables are **script-generated** (`scripts/gen_state_register_targets.py`), not typed.

**Maintenance discipline:** update a row the moment its component settles, regresses, or is superseded — and always at a gate sign-off (a gate can't sign off until its rows are current). Every value cites its source (no value from memory). `[confirm]` = reconstructed from working memory, not yet source-cited — **not settled until the citation is attached.** Full procedure: `skills/codex-state-register/SKILL.md`, wired into session-close in `DEV_CYCLE.md`.

**Status vocab:** `SETTLED` · `SETTLED-WITH-CAVEAT` · `REGRESSED` · `STALE` · `OPEN` · `NOT-SELF-SUFFICIENT` · `PENDING`

**Version: v6** · _Last updated: 2026-07-06 · By: Mr. Code + Ryan (RYA-526 set-down + 525 architecture capture: two-engine floor ratified as governing law (build-pending); N wired into production NLTE + N I grid row; Sirius scope extended to all downloads/extractions + PySME venv-isolation; solar 5/1/20/0 flagged STALE pending the RYA-527 two-engine re-freeze)._

---

## Targets & stellar parameters — MIRROR ← config/constants.py

**Authoritative source is `config/constants.py` → `STAR_PARAMS` (loaded from `config/stars.yaml`, RYA-298). The table below is a script-generated view — do NOT hand-edit it; run `python scripts/gen_state_register_targets.py --write` to refresh, `--check` to detect drift.**

NATIVE caveats on the mirrored values (these are register annotations, not source data):
- **Sun Teff = 5772 K** is authoritative (`stars.yaml`, GBS Heiter+2015). The old `5777` (legacy `STAR_SOLAR`) and `5778` (`ASTRO['Teff_sun']`, a physical constant) were single-sourced away by RYA-298/355 — the seed's "5772 or 5777" discrepancy **resolves to 5772**.
- **Procyon** canonical = GBS (Heiter+2015 6554/4.0, Jofré+2014 [Fe/H]=+0.03, ξ=1.8) — **NOT** the former Allende-Prieto 2002 spectroscopic 6530/3.96, and **NOT** 6804. `feh` is SOLVED, so `feh_ref` is a reference seed only.
- **α Cen A / B [Fe/H] mirror = +0.20**, but this is **PENDING ratification (RYA-435, Backlog)**: `stars.yaml` stores 0.20/0.20 while citing Jofré 2014 Table 3 (NLTE), which actually gives **A=+0.26, B=+0.22** — a citation/value mismatch. The α Cen A abundance run (Aug-2026 deadline, RYA-116) must NOT anchor on the unreconciled value until RYA-435 lands.
- **55 Cnc A** is the LAST target; ξ is SOLVED (init 1.0, x-check 0.85–1.35), not pinned.
- `synthetic_no_logg` is a test fixture (no fundamental logg), not a science target — shown for completeness, flagged.

<!-- BEGIN GENERATED: targets (scripts/gen_state_register_targets.py) -->
<!-- generated from config/constants.py STAR_PARAMS (stars.yaml, __version__=2026-05-28); regenerate with: python scripts/gen_state_register_targets.py --write -->
| Star | Teff (K) | log g | [Fe/H] ref | ξ (km/s) | pin/solve policy | Source (verbatim from stars.yaml) |
|---|---|---|---|---|---|---|
| Sun (Star Zero) | 5772 ± 1 | 4.438 | +0.00 | 1.00 (pinned) | pin: teff,logg,xi · solve: feh | GBS: Heiter+2015 (Teff/logg), Jofré+2014 ([Fe/H], ξ) |
| Procyon | 6554 ± 84 | 4.0 | +0.03 | 1.80 (pinned) | pin: teff,logg,xi · solve: feh | GBS: Heiter+2015 (Teff/logg), Jofré+2014 ([Fe/H], ξ=1.8); broadening: Bruntt+2010 vsini=2.8 + adopted RT vmac=6.0 (Gray scale; FT estimator deferred to RYA-316) |
| α Cen A | 5792 ± 16 | 4.3 | +0.20 | 1.10 (pinned) | pin: teff,logg,xi · solve: feh | GBS: Heiter+2015 (Teff/logg), Jofré+2014 ([Fe/H], ξ); broadening: Bruntt+2010 (MNRAS 405,1907) Table B1 vsini=1.9 + vmac=2.3 (RT, Gray scale) |
| α Cen B | 5231 ± 20 | 4.53 | +0.20 | 1.00 (pinned) | pin: teff,logg,xi · solve: feh | GBS: Heiter+2015 (Teff/logg), Jofré+2014 ([Fe/H], ξ); broadening: Bruntt+2010 (MNRAS 405,1907) Table B1 vsini=1.0 + vmac=0.8 (RT, Gray scale) |
| 55 Cnc A | 5196 ± 24 | 4.45 | +0.31 | 1.00 (solved; x-check 0.85–1.35) | pin: teff,logg · solve: feh,xi | von Braun et al. 2011, ApJ 740, 49; ξ x-check Teske+2013 / Ecuvillon+2004 |
| synthetic_no_logg (test fixture) | 5800 ± 100 | 4.4 | +0.00 | — | pin: — · solve: teff,logg,feh,xi | synthetic test fixture — no fundamental logg available |
<!-- END GENERATED: targets -->

_Target order: Sun → Procyon → α Cen A/B → 55 Cnc A (**55 Cnc LAST**). 27-element canonical list per RYA-109 (`TARGET_ELEMENTS` in constants.py)._

**A(Fe)☉ reference (NATIVE, from constants.py):** `SOLAR_ASPLUND2021['Fe'] = 7.46` (Asplund 2021, 3D-true) is the differential zero-point. The pipeline currently outputs 1D-NLTE **7.516** (~+0.05 above 3D-true) — see the Solar Fe I row for how the scale-aware gate handles this (RYA-336).

---

## Solar anchor (Star Zero) — NATIVE

| Component | Verdict | Value | Established by | Status | Reopen only if |
|---|---|---|---|---|---|
| Solar Fe I (value) | 1D-NLTE anchor — PASS | 7.516 (~+0.05 above 3D-true 7.46) | RYA-247 (NLTE unit fix), RYA-336 (scale-aware gate) | SETTLED-WITH-CAVEAT — scale-aware gate PASSES; **RYA-336 IS on main (`9ccf064`); its Backlog pill is stale** | applied 3D-NLTE correction (RYA-285 / Magic 2013) |
| Solar Fe I (scatter) | honest floor, NOT a defect | 0.138 dex (RYA-407 verdict c) | RYA-407 (PR #48, MERGED `220f263`) | SETTLED as honest floor — no cited mechanism reaches 0.10; gate threshold too strict → recalibrated to a G-anchor acceptance profile by RYA-446 (MERGED #71) | principled new vetting genuinely lowers it |
| Solar Fe II (arbiter) | balanced — PASS | synth arbiter ≈7.486–7.500, ΔFe(I−II) −0.007…−0.015 | RYA-305/341 (arbiter ratified), RYA-406 (gate scores it) | SETTLED — the ratified ionization arbiter is **SYNTHESIS** Fe II; on main (RYA-405 confirmed the solved stack is integrated, NOT regressed) | 3D-NLTE, or new synth Fe II (RYA-338 flux-space) |
| Solar Fe II (EW-path) | high BY DESIGN | ~7.700 (blend-limited EW pool) | RYA-352 / RYA-406 | SETTLED as DIAGNOSTIC — demoted to a reported line (not a gate metric); EW-vs-synth ~0.21 = documented blend-bias, within band | — |
| Solar Fe II gate scoring | reconciled to the arbiter | Fe I−Fe II −0.007…−0.015 PASS (cited tolerance) | RYA-406 (PR #47, MERGED `2ecc2a2`) | SETTLED — merged into the validated main | — |
| Solar Fe II 509-report anomaly | 7.657 / n=3 | neither arbiter (7.486) nor EW-path (7.700) | RYA-509 | OPEN — RYA-509 reproduced the banked v1 baseline, which predates the RYA-406 gate reconciliation; per-line provenance owed → RYA-515 | — |
| Solar reproducibility | from-raw deterministic | Δ=0.000 on 13 species; EW md5 `230c27c4` | RYA-509 (Mac leg; pushed, NOT merged) | SETTLED for reproducibility ONLY — **reproducibility ≠ correctness**: it reproduces the *banked v1 baseline* (hence the 509 Fe II anomaly above), it does not re-derive the corrected arbiter | cross-env reproduction (RYA-511) |
| Solar 13-species verdict (v1) | banked baseline | PASS-4 (O, C, Fe, K) | pre-June v1 baseline | STALE — predates the June Fe stack; does not reflect the RYA-406/407 gate state | RYA-251 sign-off supersedes |
| NLTE unit convention | absolute A(Fe); `_A_FE_SOLAR = 7.46` offset | — | RYA-247 fix | SETTLED for Fe; **Ca/Ti/Cr grids not yet verified for the same convention** | verify Ca/Ti/Cr grid units (RYA-245/256) |
| Solar gold reference | **v2 — verdict-sourced, tiered** | C 10.26→**8.491**; tiers: gold=C/O/K/Mn/Fe/Sc, gf_floor=Cr/Si, upper_limit=Li, **owed(held, no value)**=S/Sr/N/Co/Ti/Ni/Al/Na/Ca/P/Mg/Ba/Y/Zr/Eu/Cu/V | RYA-522 (frozen+hashed; RYA-521 channel; RYA-520 C fix) | SETTLED — `CURRENT`→v2; v1 retained immutable + SUPERSEDED (C=10.26 RYA-520 artifact) | a new ratified verdict re-freeze (v3) |
| Solar 27-el verdict | **5 PASS / 1 NLTE-OWED / 20 CURATION-OWED / 0 DATA-GAP** | PASS: O, C, Fe, Mn, K | RYA-519 (wiring) / RYA-371 (phase_c); per-element detail → RYA-463 registry | **STALE** — the live 5/1/20/0 predates the N wiring (RYA-526: N no longer NLTE-OWED) and the 491/492/520 merges (RYA-524 audit). Re-freeze on the two-engine floor = **RYA-527** (verdict v-next, gold v3) — THE Beta gate | RYA-527 re-freeze lands |

---

## Procyon (calibration check) — NATIVE

| Component | Verdict | Value | Established by | Status | Reopen only if |
|---|---|---|---|---|---|
| Procyon Fe | reproducible anchor | Fe I 7.571 (NLTE) / Fe II 7.535 | RYA-506 (clean-from-raw, MERGED PR #126 `3425995`) | SETTLED — supersedes the un-reproducible 7.593 | — |
| Procyon Fe I/II NLTE gap | diagnostic | +0.036 (NLTE widens from +0.023 LTE) | RYA-506 | OPEN — carried to the ξ-pin review (RYA-322) | RYA-322 |

---

## Pipeline integrity — NATIVE

| Component | Verdict | Value | Established by | Status | Reopen only if |
|---|---|---|---|---|---|
| theo-EW spawn/fork bug | root-caused + fixed | force-fork + all-zero guard + tmp makedirs | RYA-506 (MERGED `3425995`) | SETTLED | — |
| fork coverage | all entry points + both machines | centralized `pipeline/_runtime.py` force-fork + BLAS pins, imported by every entry point | RYA-514 (MERGED #129) | SETTLED — fork confirmed on Sirius py3.12 | — |
| Authoritative channel | **the phase_c verdict is THE single source** | raw `run()` EW output = DIAGNOSTIC-ONLY (in-file labeled); divergence guard on the gold freeze | RYA-521 (MERGED) / RYA-520 | SETTLED — consumers read `read_solar_reference('CURRENT')`; `assert_authoritative_is_verdict` rejects the raw file | — |
| Wiring integrity ("merged ≠ wired") | a merged capability can be orphaned (never invoked) — a silent-omission class | RYA-463 problem-children registry was BUILT but ORPHANED; now wired as the per-element disposition source + completeness/loud-fail guards | RYA-519 (MERGED #131) | SETTLED — run('solar') emits a value or a cited disposition for every target; 0 silent | a new orphaned merged capability |
| grade application | `line_grade` applied to EW aggregation | — | RYA-329 | OPEN — computed, not applied (silent leak). NOTE: RYA-407 found the Fe I scatter is an honest floor, so 329 is NOT the scatter driver (earlier guess retracted) | — |
| per-line engine provenance | `{engine, correction, source}` stamp per line | — | RYA-250 / RYA-512 / RYA-515 | OPEN — no per-line stamp yet; verdicts not fully auditable | — |
| Engine (all 27 elements) | Turbospectrum synthesis + MOOG EW comparison | `RADIATIVE_TRANSFER_CODE='turbospectrum'`, `EW_BASELINE_CODE='moog'` | RYA-285 (Done) | SETTLED — reverses MOOG-for-speed; Turbospectrum-via-EW not viable (RYA-234) | — |
| Two-engine floor (governing law) | every element on BOTH 1D-NLTE + synthesis; report the single best engine by a PRE-DECLARED line-quality criterion (σ / REW / blend flag / COG — NEVER reference-proximity); missing synth grid → RAISE (acquire, never silent EW-1D); rejected engine recorded but excluded from value + budget; cross-engine spread = separate diagnostic, never in the error bar | RYA-525 (ratified 2026-07-05, triggered by RYA-524 audit) | RATIFIED — **build pending** (`SCIENCE_STANDARDS.md` + selection fn + loud-fail guard not yet wired; only RYA-526's grid-presence guard core landed) | build lands (RYA-525) → upgrade to SETTLED |

---

## NLTE grids — NATIVE (selection is the register's truth; coverage is per-target)

| Species | Grid / source | Coverage note | Cite | Status |
|---|---|---|---|---|
| Fe I | Amarsi 2022 neural-network NLTE (absolute A(Fe)) | confirm bounds cover 55 Cnc (5196/4.45/+0.31) | RYA-251 Phase-3 / RYA-247 (unit fix) | SETTLED (selection); coverage per-target [confirm] |
| Ca I | MPIA MAFAGS-OS / Mashonkina 2017 (`Ca_Mashonkina2017.csv`) | confirm coverage at [Fe/H]=+0.31 | constants.py `NLTE_CORRECTION_ELEMENTS` / RYA-235 | SETTLED (selection) |
| Ti I | MPIA MAFAGS-OS / Bergemann 2011 | — | constants.py / RYA-235 | SETTLED (selection) |
| Cr I | MPIA MAFAGS-OS / Bergemann & Cescutti 2010 | Cr II deliberately excluded (RYA-240 COG artifact) | constants.py / RYA-235/240 | SETTLED (selection) |
| Na / Mg / Ba / Mn / Si | Lind 2011 (Na, INSPECT) · Bergemann MPIA (Mg/Mn/Si) · Korotin 2015 (Ba, VizieR) | feh axis is [Fe/H] RELATIVE — no solar offset (unlike Fe) | constants.py `NLTE_CORRECTION_ELEMENTS` / RYA-165/396 | SETTLED (selection) |
| C I / O I | Amarsi 2019/2020 grids (owned by `nlte_cno`) | C/O flagship prerequisite | RYA-359 | OPEN — status per C/O track |
| N I | Amarsi 2020 N grid → `N_Amarsi2020_PySME.csv` (Zenodo 3982506, md5-verified; derived via PySME on Sirius, grids never on Mac) | N I red near-LTE at the Sun (solar δ −0.0115/−0.0145/−0.0154, max\|Δ\|=0.0000 vs the RYA-369 load-test); **warm-star indicator** — 4 cool metal-rich weak-line COG nodes excluded → out-of-hull loud-flag (genuine weak-line info-loss; bracket-widening worsened it), not shipped as spurious +0.2 | RYA-526 (wired into `NLTE_CORRECTION_ELEMENTS`) | SETTLED (selection + wiring) — production resolves the solar δ live | 3D-NLTE, or a new N I grid |
| — (principle) | Amarsi/Balder **and** Bergemann/MPIA are **independent measurements**; inter-model spread feeds the uncertainty budget | not "pick one" | RYA-282 / project | SETTLED (principle) |

---

## 3D models & model atmospheres — NATIVE

| Item | Selection | Use | Cite | Status |
|---|---|---|---|---|
| 3D hydro grid | STAGGER — Magic et al. 2013 | 1D→3D Fe correction trajectory; C/O 3D path | RYA-336 / project | SETTLED (selection); 3D-RT engine not yet in hand (RYA-444) |
| Model atm (FGK dwarfs) | ATLAS9 / Castelli–Kurucz (plane-parallel, LTE) | 1D LTE base | constants.py `MODEL` | SETTLED (selection) |
| Model atm (M dwarfs) | MARCS.GES | 1D LTE base | [confirm] | [confirm] |

---

## Line lists — atomic + molecular — NATIVE

| List | Source + version | DOI / cite | Status |
|---|---|---|---|
| Atomic | VALD3 (manual web extraction; see `skills/codex-vald-extraction`) | vald.astro.uu.se | SETTLED |
| OH (mid-IR) | ExoMol **MYTHOS** 20240526 | 10.1093/mnras/stae2803 (Mitev 2025) | SETTLED — verified RYA-503 |
| NH (mid-IR) | ExoMol **kNigHt** 20240301 | 10.1093/mnras/stae1340 (Perri & McKemmish 2024) | SETTLED — RYA-503; 2kNigHt (no DOI) deferred |
| CH (mid-IR) | **MoLLIST** (Masseron 2014) | 10.1051/0004-6361/201423956 | SETTLED — RYA-503 |
| CO (mid-IR) | CO_IR_Li2015 (vendored) | RYA-360 / RYA-503 | SETTLED — CO round-trip validated (117783/117783) |
| Molecular scope | mid-IR ro-vibrational fundamentals acquired; synthesis is **v2** (Orion-gated) | RYA-503/504 | SETTLED (v1/v2 boundary) |

---

## Instruments — NATIVE

| Band | Instruments | Archive | Notes |
|---|---|---|---|
| Optical | HARPS (R~115k), ESPRESSO, UVES, FEROS | archive.eso.org | HARPS-N (TNG, ia2.inaf.it/tng); Keck/HIRES (koa.ipac.caltech.edu) |
| UV | HST **STIS** (all science gratings) + **COS G130M/G160M FUV** | mast.stsci.edu | **WFC3 EXCLUDED (RYA-119)** — corrects the pinned doc |
| Near-IR / IR | CRIRES+, NIRPS | eso.org | SPIRou / ESPaDOnS (cadc.nrc.ca); NARVAL (see Scope/drift, source [OPEN]) |
| IR reference atlases | ACE-FTS (telluric-free), NSO Kitt Peak photatl, Wallace NIR | on disk | solar IR reference / atlas fallback (RYA-513) |
| Telluric routing | CRIRES+/UVES-red/ESPRESSO-red/FEROS-red/NIRPS = molecfit + per-night GDAS · SPIRou = APERO+Wapiti | — | RYA-424 (confirmed). **NIRPS nuance:** 424 lists NIRPS under molecfit, but the actual α Cen NIRPS work (RYA-494) used the native-DRS `FLUX_TELL` product + Wallace verify — reconcile before an IR NIRPS run |

---

## Infrastructure — NATIVE

| Component | Verdict | Value | Established by | Status | Reopen only if |
|---|---|---|---|---|---|
| Sirius | **authoritative runner — BUILT** | `/mnt/codex-data`: venv312 (py3.12.13/numpy2.2.6, byte-identical to Mac), iSpec compiled, 69 GB Amarsi grids, engines, RYA-514 fork confirmed | RYA-511 Phase 0 (done) | SETTLED — all runs execute on Sirius now; reproduced Fe I 7.516/n=62 | — |
| Compute stack | **py3.12 + numpy 2.2** (reference stack, off EOL 3.9.6) | exact-zero cross-machine floor (per-line bit-identical); Mac↔Sirius drift = **null** | RYA-517 (merged) | SETTLED — Sirius alone authoritative (no dual-machine cross-confirm required) | new engine ceiling |
| Trust model | Sirius = the authoritative runner (null drift) | forced-fork + single-thread BLAS both machines | RYA-506 / RYA-511 / RYA-517 | SETTLED — "Sirius = native fork" dead (py3.14 forkserver); 517 proved null drift → Sirius authoritative | — |
| Sirius scope (downloads too) | **all downloads + grid extractions run on Sirius, never the Mac** | a Mac-side `.grd` extraction filled the disk | RYA-526 | SETTLED — extends "runs on Sirius" to downloads/extractions | — |
| venv isolation | **PySME in its own dedicated venv**; the RYA-517 reference `venv312` (astropy/pandas exact pins) stays untouched — never install into it | a transient `pysme-astro` install downgraded venv312 (astropy/pandas); restored to exact pins | RYA-526 / RYA-517 | SETTLED | deliberate reference-stack bump |
| NIRPS solar from-raw | walled | wave_THAR global FP-comb starvation vs stale 2022-11 prior | RYA-498 (parked) | SETTLED (make-do: atlas fallback) | 2023-epoch 3.3.12 WAVE_MATRIX becomes trivially available |
| α Cen NIRPS | adopt existing reduced ADP | S1D_FINAL_A (HA, 28 frames) | RYA-500 (reframed) / RYA-494 | SETTLED — absolute-mode, version-caveat; NIRPS = α Cen **A** not B (RYA-494) | fresher prior (as above) |
| Orion | planned 2nd runner (3D-RT / molecular-synthesis era) | — | project | PENDING | — |

---

## Gates — NATIVE

| Gate | Detail | Status | Cite |
|---|---|---|---|
| Solar β sign-off | 4-phase gate. Phase-1 Fe: Fe II arbiter PASS (406, merged), Fe I value PASS (336, on main), Fe I scatter threshold recalibrated (446, merged); remaining = explain the 509 Fe II 7.657/n=3 anomaly (via 515) | OPEN — much closer than "fails everything"; nothing proceeds to a science target until all 4 phases sign off | RYA-251 (supersedes RYA-101) |

---

## Scope / drift corrections — NATIVE
_Supersede the pinned project-instructions doc (a May-2026 snapshot) wherever they conflict._

| Item | Current truth | Cite |
|---|---|---|
| Target order | Sun → Procyon → α Cen A/B → 55 Cnc A (**55 Cnc LAST, not primary**) | project / RYA-109 |
| HST | STIS + COS G130M/G160M only; **WFC3 EXCLUDED** | RYA-119 |
| Telluric | instrument-aware (see Instruments) — NOT cr2res+molecfit universally | RYA-424 |
| NARVAL source | claimed **tbl.omp.eu** (PolarBase = CFHT/ESPaDOnS only) | **[OPEN]** — no ticket citation found in-repo; do NOT treat as settled until sourced |
| Target engine | Turbospectrum synthesis (all 27) + MOOG EW comparison | RYA-285 / constants.py |
| VALD | manual web extraction (per skill); synthesis-era central-depth 0.001 | RYA-389 [confirm] |

---

## Versioning & edit protocol — NATIVE

**History:** the register lives in git — every edit is a commit, so `git log` / `git diff` IS the full historical version record. No separate versioning system.

**In-file:** the header carries `Version: vN · date · by · change`. The Changelog below logs each meaningful update in one human-readable line (git holds the full diff).

**When the register gets updated (triggers):**
1. A gate signs off — **mandatory**; a gate cannot sign off while its rows are stale.
2. A verdict changes — a component settles, regresses, or is superseded.
3. A milestone is reached.
4. "We now know X works / doesn't" — a validated finding.

**Who decides:**
- BIG moments (verdict changes, new SETTLED rows, milestones) — **Ryan + Claude decide together.**
- MIRROR rows (params, abundances) — script-regenerated from source (`scripts/gen_state_register_targets.py`), automatic, no decision needed.
- Mr. Code proposes register updates in end-of-session comments; Ryan + Claude ratify the significant ones.

**Claude's standing reminder duty:** when a register-worthy moment lands in conversation — a verdict, a confirmed fix, a milestone, a "this works now" — Claude proactively says *"this belongs in the register."* Ryan need not remember to ask.

## Changelog
- **v6** (2026-07-06) — RYA-526 set-down + RYA-525 architecture capture. **Two-engine floor** ratified as governing law (RYA-525, build-pending — SCIENCE_STANDARDS.md/selection-fn/loud-fail guard not yet wired). **N wired** into production NLTE (RYA-526: `N_Amarsi2020_PySME.csv`, near-LTE at the Sun, warm-star indicator; 4 cool-metal-rich COG nodes excluded → loud-flag) + N I NLTE-grid row. **Infra:** Sirius scope extended to all downloads/extractions (Mac disk-fill); PySME venv-isolation, venv312 pinned/untouched (transient downgrade, restored) — both RYA-526 incidents. Solar 5/1/20/0 verdict flagged **STALE** pending the RYA-527 two-engine re-freeze (verdict v-next, gold v3).
- **v5** (2026-07-05) — RYA-522 set-down (Ryan-ratified checkpoint, the picture stopped moving): solar gold reference re-frozen **v2** from the verdict channel (RYA-521), **C 10.26→8.491** (RYA-520 saturated-C I-5380 fix), **tiered confidence** (gold C/O/K/Mn/Fe/Sc; gf_floor Cr/Si; upper_limit Li; owed/held incl S/Sr/N/Co — no value immortalised), v1 immutable+SUPERSEDED. Solar verdict → **5/1/20/0**. Sirius = authoritative runner, py3.12/numpy2.2, null cross-machine drift (RYA-511/517). Verdict = single authoritative channel; raw EW diagnostic-only (RYA-521). "merged ≠ wired" wiring-integrity row (RYA-463 orphaned → wired, RYA-519). +2.1-tail Sr flagged for saturation-trace.
- **v4** (2026-07-04) — VENDORED (RYA-516, Mr. Code): confirmed every row against its ticket + `constants.py`; **Targets table now script-generated** from `STAR_PARAMS` (was hand-typed). Corrections from source: Sun Teff resolved to 5772; Procyon [Fe/H] +0.03 / ξ 1.8; α Cen A/B params emitted from source (0.20) with the RYA-435 ratification caveat; 55 Cnc [Fe/H] +0.31. Merge states verified vs `origin/main` @ `3425995`: RYA-406 (#47) + RYA-407 (#48) + RYA-336 + RYA-446 (#71) confirmed **merged** → upgraded from SETTLED-PENDING-MERGE to SETTLED. Telluric [confirm] resolved to RYA-424 (+ NIRPS nuance). NARVAL source flagged **[OPEN]** (no ticket). Wired maintenance into `DEV_CYCLE.md` + `skills/codex-state-register`.
- **v3** (2026-07-04) — Fe rows corrected after verifying RYA-406/407: no Fe II regression (arbiter 7.486 balanced, on main; 406 reconciled the gate); Fe I scale-gate 336 already on main; Fe I scatter 0.138 = honest floor (446). Retracted the "329 drives scatter" guess. Added the versioning & edit protocol.
- **v2** (2026-07-04) — expanded to full personal-codex: targets & params (MIRROR), NLTE grids, 3D models, line lists, instruments; added the NATIVE/MIRROR distinction.
- **v1** (2026-07-04) — seed: solar anchor, Procyon, pipeline integrity, infra, gates, scope/drift corrections.

---

_Rows tagged `[confirm]` are reconstructed from working memory and carry citation debt — a row is not truly settled until its authoritative source is attached (per the value-provenance rule). Maintenance procedure: `skills/codex-state-register/SKILL.md`._
