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

## 2026-08-19

- **RYA-925** — Near-UV Al retracted: 3057 Å is abundance-insensitive in a 146-transition blend; appendix evidence retained, no product.
- **RYA-925** — Calibration gate corrected: Kitt Peak VIS/red-optical replicate within combined 1σ; IR/model continuation remains open.
- **RYA-925** — Al Kitt Peak matrix banked across EW/synth LTE + Amarsi; non-Fe identity defects fixed; VIS fails literature validation, no tuning.
## 2026-08-20

- **RYA-929** — full Kitt Peak/Kurucz/IAG sweep proves broad correction agreement; line-level caveats remain diagnostic-only and no abundance is promoted.

## 2026-08-19

- **RYA-929** — Kurucz 2005 Kitt Peak irradiance authenticated, checksum-pinned, staged on Sirius, and registered for independent telluric comparison.
- **RYA-927** — shared telluric route contract covers all catalog instruments; HARPS/Kitt Peak clean-line paths coexist with molecfit/GDAS correction routes.
- **RYA-926** — authority and skill governance settled; RYA-925 dry-read exposes a missing cross-instrument canary.

## 2026-08-18

- **RYA-878** — ANGLE 1 made a MEASURED angle: the production path banks a synthetic EW and one definition serves both sides. Adapter reproduces the engine **1.000**; the **engine** carries the +0.1523. Mostly definitional (range explains 70–88%), residual + per-line spread OPEN.
- **RYA-875** — the SynthesisHandler residual was a LINE-SET artifact: an 18-line median vs a scalar from a different 23-line set. Paired, the offset is **0.0000** (17/18 within ±0.009). RYA-873's declaration deleted because the numbers agree; no bar moves. ANGLE 1 still open.
- **RYA-873** — the harness term's prose is derived from PROVENANCE, so an uncharged residual stops printing "MEASURED". Contamination refuted (15/18 uncontaminated lines still span 0.50–2.07) ⇒ 0.0100 unestablished; nothing charged, no bar moves. RCA = RYA-875.
- **RYA-847** *(part 2)* — the sweep found NO transferable threshold, so the gate is the zero-parameter non-minimum check; near-UV 7.488 → 7.498.
- **RYA-869** — the harness residual follows the HANDLER, not the treatment label; `ENGINE-B-NLTE` was charged the profile fitter's 0.0129 and labelled `ProfileFitHandler` in its own budget. 4 published Fe bars 0.1705→0.1700 / 0.1731→0.1726; no value moves.
- **RYA-855** *(follow-up)* — the two mirrored harness rules in the rung audit are deleted; one rule now, in `pipeline/harness_residual.py`.
- **RYA-871** — the EW per-line artifact carries `ep_eV`; the gf resolver keys on wavelength AND EP, and the tolerance travels with the key so a keyless line is never widened. 107→130 lines priced; 0 rungs, 0 bars, 0 values move.
- **RYA-873** *(filed)* — `SynthesisHandler` is charged 0.0000 while its own banked control measured 0.0100 and PASSED, under prose reading "MEASURED ... not assumed zero".

## 2026-08-17

- **RYA-837 + RYA-843** — IR synthesis-context wiring lands; the NIR rail RCA finds a NORMALISATION failure, not a fitter one, and the real defect is that UNCONSTRAINED fits are accepted (two lines at 7.833/7.979 with chi2 flat across 8 dex). No product published.
- **RYA-848** — the CNO curvature sigma IS the published C/N/O sigma_stat; rescaled by sqrt(red_chi2) and the railed-fit sigma=0.000 arithmetic fixed. **sigma-only, proved by a same-inputs control**; solar O sigma_stat 0.041 → 0.416. CNO product set flagged STALE (banked 2026-06-27, 11 gf commits since).
- **RYA-847** *(in progress)* — three copies of the synthesis accept/reject collapsed into `pipeline/fit_constraint.py`; measure and decide split (`constraint_gate.py`); `SYNTH_CONSTRAINT` left **deliberately None** until a cross-band sweep sets the metric and the cut (RYA-161). Appendix now names excluded lines with their physical cause.

## 2026-08-16

- **RYA-819/831** — gold **v5**: Fe I provenance corrected (Magic −0.05 is wrong in magnitude AND shape); **value 7.466 unchanged**
- **RYA-834** — `canonical_gf` red edge 9199.90 → **12934.67 Å**; 28 Fe I lines adjudicated on PRIMARY LAB (NIST disqualified — agrees to 0.0003 dex, a compilation echo). 762 unblocked; products owed

- **RYA-762** — Fe 9199–13000 Å inventory banked: Engine B **187/239** (the ticket said zero); unmodelled-Pa-beta systematic RETRACTED (Hlinedata carries it); products parked on RYA-379

## 2026-08-15

- **RYA-823** — model-atom levels keyed by the coordinate each level HAS; (J,energy) ∪ term-label union. Cr I/II → REACH-UNKNOWN (was about to land `SERVED, reach 3` of 5353); Fe/Ti/Mn gain 7 rows, all gains, zero Engine-A change

- **RYA-822** — canonical_gf reaches **3000 Å**; near-UV Fe I gf graded — and it changes NOTHING (7.487→7.488): the 0.354 scatter is the Kurucz floor. 61 lab lines exist, **6** were used
- **RYA-759** — near-UV Fe I **7.487 ± 0.120** (n=40, 3000–3780 Å, 1D-LTE); the Balmer "hole" was never there — TS's own `Hlinedata` covers it

## 2026-08-14

- **RYA-759** — CORRECTION: the near-UV route was merged 2026-08-11 (#218), not "unmerged" as v71 said; gate now OPEN but a COMPENSATING-ERROR pass

## 2026-08-13

- **RYA-815** — a reference self-contradiction now withholds THAT element (INDETERMINATE, cell named) while the other 27 proceed; RYA-681's refusal preserved
- **RYA-777** — Fe frontier matrix on the surfaces + frontier standard RATIFIED; tracker regen BLOCKED by RYA-669 (phase_c refuses) and the freshness guard blind spot instrumented
- **RYA-811** — gold **v4** frozen: Fe I label 1D-NLTE → **3D-NLTE**, value 7.466 UNCHANGED; unblocks phase_c + RYA-777 A.1
- **RYA-810** — batch 4 FINAL: 97→**0** literals; the audit is now a HARD GATE, not a ratchet
- **RYA-810** — batch 3: all 19 GRID literals onto the register (116→97); proven equivalent by inode, not by string
- **RYA-810** — batch 2: `pipeline/` is now literal-free; +6 register entries and a `work` root for the Sirius drivers
- **RYA-810** — batch 1: 16 path literals onto the register, retiring ALL 15 that carried a username; `{repo_parent}` derivation means no config and no personal path
- **RYA-810** — path REGISTER (`config/path_register.yaml` + `codex_path`/`require_codex_path`); repo carries structure, roots come from env, not commits
- **RYA-806** — `telluric_applied` per holding, determined from headers and gating the arm; NIRPS is corrected where CRIRES+ is not, so the two axes cannot merge
- **RYA-805** — the 18 Vesta IDPs are NOT telluric-corrected (headers + O₂ band + r=0.996 vs water vapour); no corrected variant exists; RYA-373 shrinks to a run in RYA-797
- **RYA-711 (1+2)** — our grade is now `MQ-A/B/C/D`, never NIST's bare letters; the >25 % gf cut is derived (C = 0.0969 dex vs the ±0.10 gate)
- **RYA-796** — `crires_plus` arm for `load_window`; it REFUSES the staged IDPs (TOPOCENT + moving reflector) — unblocks RYA-797 once RYA-372/373 conditioning runs
- **RYA-796** — co-add is a ROTATION question: 2 of 5 duplicate settings sit 166°/163° apart in sub-observer longitude, the opposite face of Vesta

## 2026-08-12

- **RYA-785** — Fe Gerber Engine-B NLTE deck PASSES (+0.0579 vs published +0.06); the CHECK was an MPIA anchor + 2 blended lines ⇒ unblocks RYA-798
- **RYA-794** — the 18 CRIRES+ Vesta IDPs re-pulled from ESO (lost from us, not the archive); Y certifies ZERO solar Fe I, but J+H certify 74 ⇒ unblocks RYA-797
- **RYA-784** — ENGINE-B wired into `derive_band_products.py`; the driver reserved it and produced no value ⇒ unblocks the RYA-783 Engine-B leg
- **RYA-786** — `telluric_basis` added to `instrument_catalog.csv` + `telluric_policy.py` is the single source; splits "corrected" from "line_selection"
- **RYA-783** — IAG arm added: the telluric split recovers 94 Fe I IR lines (89 KP vs 183 IAG in-aggregate, 2.1x)
- **RYA-789** — acquired Elgueta 2026 `J/A+A/710/A111`; **sp/ holds the Vesta-derived solar spectrum** ⇒ RYA-787 takes the reduced-spectrum route
- **RYA-789** — fetch trap: the mirror links subdirs with no trailing slash, so sp/ saved as an HTML page and reported VESTA=False with every guard green

## 2026-08-11

- **RYA-776/773** — coverage reference REFRESHED post-773: Al I red-optical REACHABLE-NOT-EXTRACTED → **SERVED** (4 lines); it went stale within one merge and no CI guard caught it
- **RYA-776** — generated `data/catalog/engine_coverage.csv`: engine × wavelength reach as a LOOKUP, not a re-derivation; unblocks RYA-306/775 and gives RYA-773 its answer (Al doublets are reachable-not-extracted)
- **RYA-776** — 4 decoding traps each faked an absence (38 UNCOVERED → 2): super-levels, cumulative atom energies, continuum stage, missing `ion` column dropping Ca+Cr

## 2026-08-08

- **RYA-708** — per-(instrument × band) abundances ratified; the cross-instrument delta is a **blend diagnostic** (Al: clean line agrees to 0.019, blended disagrees by 0.092)
- **RYA-708** — coverage service rebased onto the EXISTING instrument catalog (my duplicate deleted); IAG + solar holdings registered; near-UV sweep adds 3633 lines → **5099 untouched**
- **RYA-708** — "all wavelengths, all instruments, all models" written into SCIENCE_STANDARDS; the EW pool is **HARPS-only**, 0 of 808 lines beyond 6910 Å
- **RYA-708** — Kitt Peak found + registered (2960–13000 Å, widest arm); Al now has THREE arms, KP corroborates IAG to 1–2%
- **RYA-708** — coverage service + instrument registry; **corrects v48's false Al NO-DATA** (IAG covers 7835/8772; 8772 is the element's strongest line)
- **RYA-706** — promotion gate ratified: promoted lines get a gf check. Al pilot reaches 6.431 (Asplund 6.43) but **STOPS on 0.339 dex scatter**; all 654 candidates HOLD ungraded
- **RYA-706** — stage-2 (fit→pool) ledger built: 8149 drops classified, **913 UNEXPLAINED**, **229 recovery candidates** incl. 2 clean Al lines for an element reporting nothing
- **RYA-707** — SPP Appendix A mandatory: every unresolved element defends its blank with plots + measured evidence; tool ladder (EW→NLTE→synth) is the spine; solar Al re-diagnosed from the spectrum
- **RYA-705** — refinement-debt gate made true: `In Progress` rows were counted NOWHERE (a phase could close over mid-flight debt); Ca re-homed 562→561; discharged rows now printed, 17→15
- **RYA-694 CORRECTION** — RYA-694 was auto-closed by PR #205's branch name with none of its work done; reopened. Real tickets: RYA-700/701/702 (see register v46)
- **RYA-699** — RYA-691's reliability rule ratified as the 4th RYA-674 constraint (vocabulary single-sourced so the gate can read what the loader writes); codex-mr-code-brief skill homed in `skills/`, narrowing RYA-386
- **RYA-702** — run order ratified: Sun → Alpha Cen A → Alpha Cen B → Procyon → 55 Cnc A; the second star measures the infrastructure, so it adds the least new physics
- **RYA-701** — Al→RYA-523 (saturation, not gf), Y→RYA-683+523 (ion label + saturation), S/N/P/Cu gf rows→RYA-697 (they escalated to closed RYA-161)
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
- **RYA-711** — the ELEMENT PROTOCOL lands (authored, never merged; RYA-709 cited a file that did not exist). **+Step 3a: an element's curation travels with the element, not the route** — Fe II uncurated 7.656 vs curated 7.466
- **RYA-906** — physics-axis naming: store `route/scale/model/atmos/gf`, derive the display name, keep `treatment` as a permanent dual label. **Route from the HANDLER, never the label** — `1D-LTE` is used by both routes
- **RYA-905** — solar_harps `telluric_applied` measured from flux: **not-applied**. 🔴 kpno_solar_atlas scores WORSE (51.3 vs 22.06) on the same metric — the atlases are not uniformly corrected; per-line clean-line selection is the method, not correction

