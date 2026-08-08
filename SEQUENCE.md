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

- **RYA-682** — two-engine driver inputs preflighted; numpy>=2.3 silently emptied the Engine-B artifact (generate on venv312, not venv_ci)
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
