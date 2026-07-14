# RYA-524 — Master 27×2 audit (all 27 × both engines) — READ + CLASSIFY

_Ground truth read directly at `origin/main` **b69f7a6** (post RYA-545 Ti wire + RYA-546 Mn + RYA-540 Li/Cu). Firewall: read + classify only — no gf/gold/wiring/tuning changes, no merge._

## 0. Freshness (the headline)

**The live verdict `docs/audit/solar_phase_c_verdict_rya371.md` is dated 2026-06-29; gold `solar_abundances_v2.csv` frozen 2026-07-05.** Both PRE-DATE a large stack of merged work: **RYA-491/492/520** (finish-solar O/N/S/C, S gf, C raw-EW), **RYA-526** (N NLTE grid wired), **RYA-540** (Li/Cu grids), **RYA-545** (Ti→Mallinson ab-initio), **RYA-546** (Mn→ab-initio; Cr vintage). The verdict is globally several re-runs behind → the whole table needs a re-emit/re-freeze (RYA-527, the Beta gate). Two supporting docs are also stale in-place: the wiring matrix `nlte_grid_availability.csv` still lists **`Ti_Bergemann2011_MPIA.csv`** (RYA-545 changed `constants.py` to Mallinson but not this CSV), and the RYA-463 registry is dated 2026-07-04.

## 1. The single most important output — displayed values resting on a SCALED-DRAWIN (known-inflated) leg

| rank | element | displayed value | leg it rests on | vintage | stakes / mitigant |
|---|---|---|---|---|---|
| **1** | **Fe I / Fe II** | **7.516 gold PASS (THE anchor)** | `Fe_Bergemann_MPIA.csv` (Engine-A live; Amarsi-2022 NN MLP present but ARCHIVED) | **SCALED-DRAWIN (suspect)** | HIGHEST — the whole solar scale hangs off Fe. Mitigant: ionization-balance-gated (RYA-407, arbiter 7.486, scatter 0.139 honest floor) constrains the correction; Δ vs Asplund only +0.056. Still: the anchor's NLTE leg is a pre-2020 MPIA grid, and the register's *selection* table claims Amarsi-2022-NN — code≠selection. **Confirm Fe MPIA H-collision recipe + reconcile register-vs-code.** |
| **2** | **Cr** | **6.022 gf_floor CURATION-OWED** | `Cr_Bergemann2010_MPIA.csv` + 3D | **SCALED-DRAWIN, SH=0** | The +NLTE(+3D) correction is applied to the displayed value; RYA-546 flagged SH=0 → likely spurious (→ drop-to-LTE gate, §3). Cr is also gf-floor-owed, so the value isn't gold, but the leg is inflated. |
| 3 | Ca | (owed, no gold value) | `Ca_Mashonkina2017.csv` | SCALED-DRAWIN (MPIA) | LOW — value is owed (gf floor); Engine-B `ca105b` (ab-initio) available; Ca NLTE small (~+0.02). Note, don't rush. |
| 4 | Ba / Sr | (owed / wrong-species) | `Ba_Korotin2015` / `Sr_Bergemann2012_INSPECT` | SCALED-DRAWIN (pre-2020) | LOW — both ionised-resonance, **near-LTE** (δ≈−0.01), so the H-collision recipe is immaterial; Sr's real problem is wrong-species (§2). |
| — | Mn | 5.470 gold PASS | **displayed = ab-initio Amarsi HFS δ+0.024 (RYA-546) ✓** | AB-INITIO | NOT on scaled-Drawin for the displayed value. BUT the EW-path grid *pointer* `Mn_Bergemann_MPIA.csv` is still scaled-Drawin (superseded, repoint = RYA-527) → a non-HFS re-run silently reverts to +0.108. FRAGILITY flag. |
| — | Ti | 5.471 (gf-owed) | **displayed EngA = Mallinson ab-initio ✓ (RYA-545)** | AB-INITIO | Displayed value is clean. Engine-B leg (`ti503b`) is still scaled-Drawin but NOT displayed (follow-on RYA-548). |

**Bottom line: Fe (the anchor) and Cr are the two displayed values resting on a scaled-Drawin leg that matter.** Ti and Mn — the two we just fixed — are clean on their displayed leg.

## 2. The 27×2 table

Cols: verdict (value / tier / method / Δ vs Asplund) · Engine-A grid / vintage · Engine-B atom / vintage / gate δ · RYA-463 registry · done-tickets · classification.
Vintage: **AI**=ab-initio · **SD**=scaled-Drawin · **LTE** · **VOID**=NLTE-void.

| El | Verdict (val / tier / method / Δ) | Engine-A / vintage | Engine-B / vintage / δ | Registry | Tickets | Classification |
|---|---|---|---|---|---|---|
| Fe I | 7.516 / gold PASS / EW+NLTE / +0.056 | Fe_Bergemann_MPIA / **SD** | not in deck / LTE | BAD_GF(Procyon-pred) | 319,406,408,446,506,534 | **WIRED-OK + VINTAGE-FLAG** (anchor on SD; balance-gated) |
| Fe II | ~7.486 arbiter (balance) / diagnostic | (ion 2, same grid) | not in deck / LTE | (arbiter) | 305,341,406 | WIRED-OK |
| C | 8.491 / gold PASS / synth / +0.031 | Amarsi2019 nlte_cno / **AI** | not in deck / LTE | C I 5380 SAT=exclude | 359,237,491,493 | WIRED-OK |
| O | 8.735 / gold PASS / synth / +0.045 | Amarsi2019 nlte_cno / **AI** | atom.o41f / **AI** / −0.105 PASS | [O I]6300 CONT_LIM | 359,534,449,447,365,483 | WIRED-OK |
| N | 8.202 / owed NLTE-OWED / atlas 1D / +0.372 | N_Amarsi2020_PySME / **AI** (wired RYA-526) | not in deck / LTE | DATA_GAP owed | 369,526,460,491 | **DONE-BUT-STALE-VERDICT** (grid wired post-verdict; apply −0.37) |
| Mg | — / owed / EW culled / — | Mg_Amarsi2020_PySME / **AI** | atom.mg86d / **AI** / −0.023 PASS | SAT_COG=synth | 410,534,465,354,165 | GENUINELY-OWED (gf cull + b-triplet sat) |
| Si | 7.888 / gf_floor / EW+NLTE+3D / +0.378 | Si_Amarsi2020_PySME / **AI** | atom.si340 / **AI** / −0.033 PASS | BAD_GF | 410,534,399,165,417 | GENUINELY-OWED (gf floor → 161/162) |
| S | 7.753 / owed / EW+NLTE / +0.633 | S_Amarsi2025_PySME / **AI** | not in deck / LTE | (low-conf) | 402,**492**,491,361 | **DONE-BUT-STALE-VERDICT** (RYA-492 gf → ~7.486; "expect 7.12" refuted) |
| Ca | 6.324 / owed / EW+NLTE / +0.024 | Ca_Mashonkina2017 / **SD** | atom.ca105b / AI(likely) / −0.009 PASS | (none) | 235,411,413,534 | GENUINELY-OWED + EngA-SD (owed anyway; EngB AI avail) |
| Ti | 5.471 / owed / EW+NLTE+3D / +0.501 | **Ti_Mallinson2024_PySME / AI (RYA-545)** | atom.ti503b / **SD** / +0.221 CHECK | BAD_GF | 545,544,542,535,546,534,399 | GENUINELY-OWED(gf) + **EngB-SD→RYA-548**; displayed EngA=AI ✓; matrix stale |
| Co | 6.128 / owed / atlas 1D / +1.188 | none / LTE | atom.co247qm / AI(likely) / +0.099 PASS | CONT_LIM | 534,460 | GENUINELY-OWED (blue-edge, not trusted) |
| Ni | 6.946 / owed / EW / +0.746 | none / LTE | atom.ni538qm / AI(likely) / +0.018 PASS | BAD_GF | 534,365,450,543 | GENUINELY-OWED (gf floor → 161/162) |
| Na | 6.264 / owed / EW+NLTE / +0.024 | Na_Amarsi2020_PySME / **AI** | Gerber Na (RYA-533 sep. gate) / **AI** / −0.068 PASS | (none) | 402,533,529,465,354,165 | GENUINELY-OWED (thin pool → 161/162) |
| Al | 7.406 / owed / EW+NLTE / +0.976 | Al_Amarsi2020_PySME / **AI** | not in deck / LTE | (none) | 402,354,534 | GENUINELY-OWED (gf floor → 161/162) |
| K | 5.099 / gold PASS / atlas+NLTE / +0.029 | K_Amarsi2020_PySME / **AI** | not in deck / LTE | NLTE_OWED RESOLVED | 402,462 | WIRED-OK |
| P | 6.610 / owed / atlas 1D / +1.200 | none / LTE | not in deck / LTE | DATA_GAP | 460 | GENUINELY-OWED (gf/data → 161/162) |
| Ba | — / owed / EW culled / — | Ba_Korotin2015 / **SD** (near-LTE) | atom.ba111 / SD-cand / −0.018 PASS | SAT_COG=synth | 534,396,354,165 | GENUINELY-OWED + EngA-SD(near-LTE, low-stakes) |
| Y | — / owed / EW culled / — | none / LTE | not in deck / LTE | DATA_GAP Y II | 458 | GENUINELY-OWED (data gap — Y II absent) |
| Eu | — / owed / EW culled / — | none / LTE | not in deck / LTE | HFS_SUM | 354,458,371 | GENUINELY-OWED (HFS-synth Eu II 6645) |
| Mn | 5.470 / gold PASS / HFS synth+NLTE / +0.050 | Amarsi2020 GALAH **AI** δ+0.024 (RYA-546); pointer Mn_Bergemann_MPIA=SD SUPERSEDED | atom.mn281kbc / **AI** / +0.043 PASS | HFS_SUM RESOLVED | 546,476,473,411,468,534,354 | **WIRED-OK (displayed AI ✓)** + FRAGILITY (EW-path pointer SD → RYA-527) |
| Cr | 6.022 / gf_floor / EW+NLTE+3D / +0.402 | Cr_Bergemann2010 / **SD (SH=0)** | not in deck / LTE | BAD_GF(canary) | 546,235,399,240 | **VINTAGE-INFLATED**; drop-to-LTE gate = **STOP-PRECISION-LIMITED** (§3; Cr II n=1) + gf-floor-owed |
| V | 3.917 / owed / synth LTE / +0.017 | none / **VOID** | not in deck / LTE | NLTE_VOID | 466,354 | **NLTE-VOID** (WIRED-OK, LTE-flagged; RYA-470 V II anchor) |
| Sc | 3.203 / owed(gold-conf) / atlas 1D / +0.063 | none / LTE | not in deck / LTE | HFS_SUM | 460 | GENUINELY-OWED (blue-edge HFS single line) |
| Cu | 4.345 / owed / synth+NLTE / +0.165 | Cu_Caliskan2024_PySME / **AI** | not in deck / LTE | HFS_SUM | 540,466,402,354 | WIRED-OK (fragile: b-factor .grd offline) + gf residual |
| Zr | — / owed / EW culled / — | none / LTE (Zr II majority) | not in deck / LTE | SAT_COG=synth | 354,526 | GENUINELY-OWED (measurable, synth) |
| Li | 0.727 / upper_limit / EW / −0.323 | grid present UNWIRED (RYA-540 STOP) / LTE | not in deck / LTE | MOLECULAR_BLEND=UL | 540,103,458,371 | WIRED-OK (terminal upper-limit by design) |
| Sr | 4.961 / owed / EW (**Sr I 6617 SUSPECT**) / **+2.131** | Sr_Bergemann2012_INSPECT / SD (near-LTE) | atom.sr191 / SD-cand / −0.013 PASS | NLTE_OWED Sr II owed | 421,534,433,428 | **WRONG-SPECIES** (verdict reads Sr I EW; wired = Sr II near-LTE) |

## 3. Cr I/II LTE ionization balance (folded-in RYA-546 Cr drop-to-LTE gate)

Ran on Sirius — `scripts/rya524_ion_balance.py --element Cr`, MARCS solar, weak/moderate lines (linear COG), Cr II = untouched validator, railed inversions dropped not clamped.

```
Cr I : n=10/16  A_LTE = 6.058  MAD 0.236  SEM 0.114
Cr II: n=1/2    A_LTE = 5.957  (only 4588.20 survives; 6053.47 railed)  SEM n/a
LTE balance A(Cr I) − A(Cr II) = +0.102   (gate |.| < 0.05)
```

**Verdict: STOP-PRECISION-LIMITED.** Cr I SEM 0.114 ≫ 0.05 and Cr II has only **n=1** usable weak line (the rest of our 5 measured Cr II lines are saturated, >80 mA — off the linear COG). Per the RYA-546 Cr-brief discipline ("SEM ≫ 0.05 → STOP + flag for Ryan; do NOT certify a retire on an uncertain balance"), **the drop-to-LTE cannot be certified on our solar Cr data** — the same thin-majority-ion precision wall that STOPped Ti (RYA-545). Notes: (a) both ions sit ~+0.35–0.44 above Asplund 5.62 = the Cr gf-floor (BAD_GF canary), consistent with the verdict; (b) the point-estimate imbalance is +0.10 with Cr I *above* Cr II — and the scaled-Drawin SH=0 NLTE correction is *positive* (pushes Cr I further up), so it would worsen balance if trusted — directionally supporting "spurious/inflated," but NOT certifiable. **Cr stays VINTAGE-INFLATED + gf-floor-owed; the retire is precision-blocked (flag for Ryan / needs more clean Cr II weak lines).**

## 4. Classification tally

- **WIRED-OK: 9** — Fe I, Fe II, C, O, K, Mn, V(void-flagged), Cu, Li(terminal). _(Fe carries the VINTAGE-FLAG; Mn the fragility flag.)_
- **GENUINELY-OWED: 13** — Mg, Si, Ca, Co, Ni, Na, Al, P, Ba, Y, Eu, Sc, Zr. All behind RYA-161(Done)/162(WIP) differential-gf + curation; honest, correctly reflected.
- **DONE-BUT-STALE-VERDICT: 2** — N (grid wired RYA-526, apply −0.37), S (RYA-492 gf → ~7.486).
- **WRONG-SPECIES: 1** — Sr (Sr I EW +2.13 vs wired Sr II).
- **VINTAGE-INFLATED: 1** — Cr (SH=0 scaled-Drawin; §3 gate = STOP-precision-limited, retire blocked). _(Ti Engine-B is vintage-inflated too but not displayed → RYA-548.)_
- **NLTE-VOID: 1** — V.
- **MERGED-NOT-WIRED: 0** (closed by RYA-519/462/466). **Unclassified: 0.**

## 5. Batched fix plan (fixes are separate approved batches — this audit only produces the list)

**Batch A — RE-DERIVE (Engine-A scaled-Drawin → ab-initio):**
- **Fe (anchor, HIGH)** — confirm the `Fe_Bergemann_MPIA` H-collision recipe; if scaled-Drawin, migrate to the ab-initio Amarsi-2022 NN grid the *register already names* (reconcile code-vs-selection), gated on Fe I/II ionization balance. High-stakes (the whole scale).
- **Ca** — Mashonkina2017(SD) → Amarsi Ca ab-initio `.grd` (staged, RYA-505 hierarchy), validate-don't-tune cross-check. LOW priority (owed, small NLTE).
- (Ba/Sr SD legs are near-LTE → no re-derive needed; note only.)

**Batch B — ATOM-SWAP (Engine-B old→ab-initio):** Ti `atom.ti503b` → ab-initio (**RYA-548 already open**; flips the RYA-534 xfail honestly).

**Batch C — DROP-TO-LTE:** Cr — **gate ran (§3) = STOP-PRECISION-LIMITED** (Cr II n=1, SEM 0.114). The retire is NOT certifiable on our solar Cr data → **flag for Ryan; hold Cr on the scaled-Drawin leg with a VINTAGE-INFLATED tag** until either more clean Cr II weak lines are measured (to break the precision wall) or a science call accepts the drop on the literature + directional evidence. Cr's dominant issue is the gf-floor anyway (Batch E).

**Batch D — WIRE / REROUTE / RE-EMIT:**
- **Sr WRONG-SPECIES (flagship)** — reroute the verdict to the wired Sr II; complete Sr II 4077/4215 measurement (RYA-421 Step 1) + unblock silent Sr II line-drop (RYA-429). Retires the bogus +2.13. Folds RYA-523.
- **N stale** — apply the wired RYA-526 NLTE (−0.37) or ratify data-channel-owed.
- **S stale** — re-emit with RYA-492 Costa-Silva gf (7.753 → ~7.486).
- **Doc sync** — update the stale wiring matrix (`nlte_grid_availability.csv` Ti row still says Bergemann2011); repoint the Mn EW-path grid pointer (RYA-527).
- **Verdict re-freeze (umbrella RYA-527)** — the verdict is dated 2026-06-29; one re-run folds N/S/Ti/Mn + O/C refinements. **THE Beta gate.**

**Batch E — OWED / VOID (honest, no action beyond the named survey):**
- V — NLTE-VOID → V II ionization-anchor (RYA-470).
- The 13 GENUINELY-OWED gf-floor/thin-pool elements → RYA-161/162 differential-gf survey.

## 6. Beta-gate read
25/27 of the owed-list is honest and correctly reflected. **Do NOT sign Beta until:** (1) **Sr** wrong-species is rerouted, (2) **N/S** stale values re-emitted, (3) the **Fe-anchor scaled-Drawin vintage** is confirmed/reconciled (highest-stakes new finding this pass), and (4) the verdict is re-frozen (RYA-527) so it stops being several merges behind. Cr per §3.
