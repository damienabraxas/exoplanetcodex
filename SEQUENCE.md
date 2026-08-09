# Codex Ticket Sequence — chronological landing log

**Read this second, after `LEDGERS.md`, for a quick "what happened recently" catch-up.**

One line per merged/landed ticket. Newest first. This is a narrative overlay on
the register — the register carries the deep "why" in `CODEX_STATE_REGISTER.md`;
this file carries the sequence.

**Discipline:** bump this at the same time you bump `CODEX_STATE_REGISTER.md`,
in the same PR. Append-only (never rewrite history; if a landing was wrong or
superseded, add a new line noting the correction, don't edit the original).

**Format:** `- **RYA-XXX** — one-sentence summary; what it unblocks`. Keep it
under ~140 chars per line. If you need more, it belongs in the register, not here.

---

## 2026-08-08

- **RYA-700** — run order ratified: Sun → Alpha Cen A → Alpha Cen B → Procyon → 55 Cnc A; the second star measures the infrastructure, so it adds the least new physics
- **RYA-694** — systems catalog reconciled with the public site: 19 systems, 5 published-but-untracked added, spectral types backfilled, Copernicus recorded
- **RYA-695** — Kitt Peak wired as Engine B (N/K/Sc); P I refused as a RAILED fit; `none-published` splits impossible tasks from unpulled grids; tracker gains chosen_engine/selection_reason/models_tried; Phase 3 re-emit, Ba 2.237
- **RYA-690** — register header collapsed 8 lines → 1, three orphaned landings rehomed (v33/v37/v38); structure guard + SEQUENCE merge=union
- **RYA-680 + RYA-691** — Co I (4.960) and Ba II (2.237, the RYA-581 deblend, NOT 559's 2.410) wired into _dedicated_engine_B(); `reliable` honoured at every read, silent NLTE→LTE `or` killed (was live for V); no value moved
- **RYA-692** — LTE_ONLY_BY_DESIGN split from NO_MODEL_ATOM: the wiring audit no longer reports P/Sc/Eu's ratified LTE-only disposition as missing atoms; 3 phantom debt rows retired
- **RYA-676** — refinement debt architecture: element_refinement_registry.csv (SSOT) + refinement_debt tracker column + LEDGERS.md + CI guard extension + codex-mr-code-brief skill pre-check; structurally prevents the RYA-524 → 581/585/565 orphan class recurring
- **RYA-585** (backfill; landed PR#189 2026-08-08, no SEQUENCE line at the time) — Zr II deblend fixed the MODEL (rchi2 83→0.39); lines intrinsically insensitive, Zr stays owed, line set declared EXHAUSTED
- **RYA-675** — staleness detector narrowed: artifact_age_stale vs cross_channel_disagreement are now distinct signals with distinct remedies; unblocks honest Ca promotion decision at v4 freeze
- **RYA-674** — ratified constraints re-checked at EVERY emission (registry + gate + SCIENCE_STANDARDS); Fe/Li/CrII protected structurally; `--gold-version` unblocks re-emit
- **RYA-679** — ONE reliability rule; red_chi2 ceiling RETIRED (60.0's sigma_flux rationale refuted); 80 records re-adjudicated, 0 flips; Sr II 4077 kept
- **RYA-682** — two-engine driver inputs preflighted; numpy>=2.3 silently emptied the Engine-B artifact (generate on venv312, not venv_ci)
- **RYA-681** — Fe 1D→3D guard re-keyed on the VALUE + scale-identity gate check; 7.416 now fails; phase_c LOUD-FAILS on gold v3 pending a v4 (RYA-669)
- **RYA-684** — isotope fraction double-applied on 5 VALD-list species (Eu +0.3002); NO live value exposed; convention guarded
- **RYA-686** — a result cannot land without its generator (GENERATORS.yaml + CI guard); RYA-559's Ba harness recorded UNREPRODUCIBLE
- **RYA-581** — Ba II 5853 deblended by in-window blend fit: A(Ba) 2.410 → 2.237, verdict → PASS but HELD (one line, gate 3 unevaluable)
- **RYA-673** — Engine A/B wiring audit across all 27 canonical species; **10 `neither` / 7 `B_only` / only Fe of the 6 PASS elements confirmed on both engines**; 6 synthesis-required species have no Engine B; per-element wiring tickets owed

## 2026-08-07 — gold v3 freeze + pre-527 cleanup

- **RYA-565** — Eu II LTE HFS synthesis DONE → owed-no-value (dEW/dA 13.9 vs floor 40, linear-COG); "finished treatment" wording retired in both registries
- **RYA-668** — Science Product Package (SPP) framework doc landed at docs/SCIENCE_PRODUCT_PACKAGE.md; peer of Glossary/Method/Science-Architecture; unblocks per-star SPP tickets
- **RYA-664** — Na Gerber prov gate block populated (writeback of RYA-533 result); Engine-B gate 1 clears for Na
- **RYA-665** — gold v3 FROZEN (Fe I 7.466 / Mn I 5.466 gold; Co PASS at owed-HELD; N off NLTE-OWED; Ba phantom killed); CURRENT→v3; Phase 1 of split RYA-527 complete
- **RYA-663** — pre-527 per-element disposition report generated; Ca = provisional flip, Na = candidate second flip pending RYA-664
- **RYA-654** — element_status_tracker becomes GENERATED from phase_c + editorial sidecar; physics_regime GET-DATA divergences adjudicated on EW-vs-synthesis axis
- **RYA-653** — shared blank-cause honesty tripwire extracted to `pipeline/provenance_honesty.py`; corrected gold candidate for Ba phantom (promoted to v3 via RYA-665)

## 2026-08-06

- **RYA-660** — Sirius storage crisis: OS drive stalled + remounted read-only, recovered same night; 30→180s SCSI timeout udev rule as mitigation; host-side cause remains OPEN RISK; follow-ons RYA-661/662

## 2026-08-05

- **RYA-313 / RYA-314** — CI + merge gate real for first time in repo history; Sirius self-hosted runner; `CI/test` required; merge-commit style ratified; RYA-506 iSpec makedirs + NumPy 2.0 `np.trapz` regression fixed

## 2026-08-04

- **RYA-659** — register re-synced from 11-ticket drift (RYA-556..652 backlog); `LEDGERS.md` startup index created; register-freshness CI guard wired

## 2026-07-17

- **RYA-553** — solar Fe 1D→3D correction APPLIED (7.516 → 7.466 on true 3D scale); `FE_GATE [7.41,7.51]` restored as real solar gate; unblocks gold v3 freeze

## 2026-07-14

- **RYA-549** — Fe anchor vintage confirmed BENIGN (MPIA δ +0.010 solar, ionization-balance-gated); ab-initio migration deferred to RYA-550 low-pri; unblocks RYA-527

## 2026-07-13

- **RYA-545** — Ti I wired onto Mallinson-2024 ab-initio grid (retires Bergemann-2011 scaled-Drawin +0.108); solar δ +0.0506; corroboration-accept
- **RYA-546** — Mn NLTE re-based to ab-initio Amarsi-2020 +0.024 (was scaled-Drawin +0.107); reverses RYA-411; A(Mn)☉ 5.554 → 5.466 PASS

## 2026-07-10

- **RYA-525** — two-engine floor BUILT (per-line reference-blind selector + inverse-variance aggregation + `CROSS_ENGINE_MIX_GATE` + loud-fail guards)

## 2026-07-09

- **RYA-534** — Family-A TS-Gerber NLTE rollout completed (10/11 clean; Ti CHECK honest strict-xfail; atom swap owed as RYA-548)
- **RYA-530** — capability-sweep reconciliation: 14 BUILT-but-unregistered capabilities registered; Na/Mg/Si NLTE source drift corrected to Amarsi-2020 PySME
- **RYA-361** — closed as Duplicate of RYA-533 (Gerber TS-NLTE vision executed under new forensics-first framing)

## 2026-07-06

- **RYA-533** — TS-native Gerber NLTE deck BUILT + Na-validated on Sirius (median δ −0.068 vs −0.107, PASS); Turbospectrum_NLTE v20.1 compiled
- **RYA-531** — corrected RYA-529 PySME-wholesale row to two-family NLTE-derivation map + Engine-A/B distinction
- **RYA-526** — grid coverage: N wired (RYA-369), missing grids acquired, Mn grid vendored

## 2026-07-05

- **RYA-522** — solar gold reference v2 re-freeze from verdict channel (tiered confidence); C 10.26 → 8.491 (saturated C I fix, RYA-520)
- **RYA-517** — reference stack ratified as py3.12+numpy 2.2; null cross-machine drift confirmed
