# Codex Ticket Sequence — chronological landing log

> **Current as of `main` 36d5c3d5** (2026-09-17, RYA-1184) — the landings below were rebuilt from `git log --merges --first-parent origin/main`, not from the previous text of this file. The 2026-08-31 and 2026-09-03 blocks were missing entirely: this log had stopped at 2026-08-30 while sixteen PRs merged.

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

## 2026-09-17

- **RYA-1222** — orchestrator MVP: `run_pipeline.py --star X --element Y` drives the band×instrument×engine matrix idempotently, one run-report, loud-fail-continue; never writes the science feed
- **RYA-1184** — SEQUENCE-vs-git reconciler + register Version-pointer guard + BLOCKING pre-push hook; re-arms the RYA-659 check while CI is off (RYA-954)

## 2026-09-15

- **RYA-1218** — Si preflight checkpoint (campaign still In Progress); fixed a source-provenance defect where Si II 6371.370 A carried the Si I Garz scale from the RYA-1169 generator
- **RYA-1217** — Al restart Gate 0 CLOSED on main c64eccfe: the campaign cannot start from the required verified pool

## 2026-09-11

- **RYA-1212** — publish gate refuses the 0.17 UNGRADED_GF blanket (no override) and wires RYA-968 per-line gf sigma; 0 of 92 live products carry 0.17 and the 4 Reference bars drop 0.171 → ~0.049–0.052, A and n_lines unchanged
- **RYA-1211** — premise refuted: the 15 "K07" Reference lines were never Kurucz-sourced; the gf join missed because AGSS21 prints nm to 2 dp (0.1 Å) against a 0.02 Å tolerance. 14 of 21 upgraded, citable-σ coverage 5/21 → 17/21, no log gf changed
- **RYA-1210** — published the Fe I / Fe II appendix pages, tracker refresh and social forest from the live feed (v1.158 → v1.163, no abundance, membership or line-count change); the four Asplund-reference systematic terms moved 0.17 → 0.0475 dex
- **RYA-1208** — Gerber matrix v2, rebuilt on current main because merging the stale branch would have reverted RYA-1209's near-UV values; 17 of 19 cells published, and it carries the RYA-1206 `_unlabelled` UnboundLocalError fix main still lacked
- **RYA-1203** — post-RYA-1191 NIR re-ingest rebuilt on current main; clears exactly the two SHA_MISMATCH rows from the RYA-1080 guard and makes all three KP-molecfit NIR Fe I legs agree on one guarded pool (n=23/6/23)

## 2026-09-10

- **RYA-1209** — dropped the K07 line from the near-UV pool so the budget clears at gf rung 3; all four Fe I near-UV products carry post-molecular values below the 7.596/7.642 they replaced, no value calibrated

## 2026-09-09

- **RYA-1207** — near-UV molecular opacity into production through a SEPARATE `.bsyn` list plus a per-band `use_molecules` (the atomic list is unchanged); paired lever −0.050 dex on Fe I, −0.29 on Fe II; 4 Fe II published, 4 Fe I refused on a mixed-gf pool
- **RYA-1206** — CRIRES+ H-arm ENGINE-B-NLTE product emitted (A=7.549, n=9, n_excluded=14); the 12 unlabelled lines are TRANSITION-ABSENT from the Gerber deck, so the boundary is excitation not wavelength and Gerber/ENGINE-A are complementary in H
- **RYA-1204** — near-UV opacity physical fix: two housekeeping fixes plus diagnostics; both levers were measured through monkeypatched runs and never installed, so nothing published moves

## 2026-09-08

- **RYA-1183** — CNO closure test fix (test-only, `tests/test_cno_closure_rya1136.py`). ⚠️ *Sourced from the merged diff alone — this ticket carries NO end-of-session comment; see the RYA-1184 adjudication block.*

## 2026-09-07

- **RYA-1196** — telluric consumer wired to pass its instrument; latent-defect close (the EW route is dead — 0 of 70 live products use it) plus an AST guard over five production paths so a resolver call carrying only a wavelength fails the suite
- **RYA-1195** — retired the CNO manifest's `date.today()` for a static date anchored to the vendoring commit; the stamp was false, not merely non-deterministic, and one unrelated merged commit had already overwritten it *and* reverted RYA-1150's reconciliation
- **RYA-1194** — re-keyed the telluric consumer axis per holding via a new `holding_basis()` rather than re-keying `basis()`; verified measurement outranks the registry label, and exactly ONE decision of 26 holdings moves
- **RYA-1193** — pinned the IR telluric window policy: `ENUMERATION_COMPLETE_TO_A` ratified as a standing concept, so adding spot-flagged bands is a FLOOR and can never relax the readiness gate into looking like a completed survey
- **RYA-1192** — CRIRES+ compared against RAW for the first time (Elgueta's spectra and our Vesta IDPs are the same night, reduced twice): Y verified-corrected, H corrected in 15 of 18 windows; "12 lines in inter-order gaps" was a `len>=20` threshold artifact and the arm serves 25 of 29
- **RYA-1191** — telluric close-out, with its own O₂-γ finding WITHDRAWN: the difference template was mostly a reduction difference. A displaced null proves registration; only a clean-window null proves the template. The deep-band verdicts survive
- **RYA-1182** — separated 8,977 molecular rows out of `canonical_gf.csv` by line surgery (`git diff --numstat` exactly `0 8977`), atomic store intact; the strongest verification self-skips once committed, so the standing guards are the invariant check and the RYA-945 tripwire
- **RYA-1181** — carried the CNO identity fields the ingest already parsed and discarded (J″ survived only inside gf, unreadable; 405 of 408 rows now carry one) and staged 105,858 held UV transitions nothing had counted
- **RYA-1180** — typed all 11 CNO sources and re-graded Li2015_CO as a REDISTRIBUTION (twice derived) carrying the intake's entire 80-row match class; replaced six byte-identical constants rows with dissociation energies parsed positionally from table1.dat's FOURTH (adopted) column
- **RYA-1179** — made the line-key join guard scope-aware: the EP test was FUNCTION-scoped, so one `ep` laundered every wavelength-only comparison below it. Repo scan 0 → 19 findings across 15 files, 11 genuine joins including the live feed writer
- **RYA-1172** — typed the CNO atomic gf pedigree from the compilers' own words: C/N are a 2006 MCHF partial update, O rests on the 1996 Monograph 7 OPACITY Project — one "NIST grade" label spanning ten years and two methods; 403 higher-stage rows read NOT_ESTABLISHED
- **RYA-1171** — corrected the O I 777 triplet's NIST grade A+ → A (σ 0.0086 → 0.0128, 1.49×) with `log_gf` byte-identical on all 169,703 rows; the over-claim had three copies, and the sweep found 2 suspected others (Li I 6707 ×2)
- **RYA-1170** — made the intake DOI resolve from `bibliography.csv` at build time instead of being copied, so divergence is impossible rather than detectable; the copied DOI had drifted onto a DIFFERENT PAPER, and the SSOT holds 2 of the 11 sources this intake cites
- **RYA-935** — live status/tracker refresh (`live_status.json`, `live_tracker.html`, model-availability findings). ⚠️ *Sourced from the merged diff alone — this ticket carries NO end-of-session comment; see the RYA-1184 adjudication block.*
- **RYA-515** — Fe per-line provenance: 12 live products contradict their own committed evidence after RYA-1191's re-run, and the HARPS pair is the dangerous shape — the re-measurement went to a different wavelength stem, so stale evidence still vouches for the stale product; resolution re-keyed on `selector`, not `tier`

## 2026-09-04

- **RYA-1135** — first Fe II ⟨3D⟩-LTE product (A=7.640, n=9; +0.071 over its 1D-LTE sibling is an atmosphere shift, not NLTE physics), and `nlte_ion_capability` became a GATE inside `assert_linelist_supports_nlte` — 854 labelled Fe II lines in-window meant one such line would have emitted LTE under an NLTE label
- **RYA-853** — Fe I lab-gf integrity close-out: the "genuine bad row" WITHDRAWN (a parser had read Belmonte's comparison column), leaving 464 refereed rows at 0 log gf and 0 σ mismatches without a byte of `canonical_gf` changing; 37 corrections remain Ryan's call
- **RYA-715** — regenerated the Fe II dossier from the artifacts and measured its open decision MOOT as posed: one of the two gates is retired, all six contested lines are shallow and non-lab-gf, and no live Fe II product uses the EW route the cull governs

## 2026-09-03

- **RYA-1214** — first Solar CNO products: the merged diff creates `data/products/solar/{C,N,O}.json`. ⚠️ *Sourced from the merged diff alone — this ticket carries NO end-of-session comment; see the RYA-1184 adjudication block.*
- **RYA-1213** — Reference Grade across all five bands on the lab pool with no depth gate. ⚠️ *Sourced from the merged diff and the RYA-1221 post-merge audit — this ticket carries NO end-of-session comment; see the RYA-1184 adjudication block.*
- **RYA-1187** — built the (holding × band × engine) applicability matrix, 60 holding-cells + 32 engine-cells with 0 silent empties; "zero Deep Grade in red-optical/NIR" is NOT a gap, and "Reference appears only in VIS" is FALSE — 17 AGSS21 lines fall outside VIS
- **RYA-1190** — frontier near-UV opacity (diagnostic): Part A refuted as posed — our near-UV list already IS the VALD extract, 0 lines to add, payoff of a deeper re-extraction ≤0.028 dex; red-optical is uncorrected telluric, NIR-H undetermined
- **RYA-1189** — continuum root-cause analysis (diagnostic): the near-UV +27.9% shift is BLEND-DRIVEN, not a continuum error (59 lines/Å, 0 of 10 clean side-bands, no isolated Fe I line exists there); red-optical and NIR-H are the real candidates
- **RYA-1055** — Fe II NLTE capability limit established as a DECK limit: `atom.fe607a` carries 12,635 bound-bound transitions and not one involves an Fe II level; all 10 live Fe II products carry the stamp and the two ENGINE-B-NLTE cells are annotated, not deleted

### Corrections (RYA-1184, 2026-09-17) — append-only, prior lines left intact

- **Ordering defect:** a `## 2026-09-15` section sits at the END of this file (below `## 2026-07-05`), where the RYA-587 landings were appended to a newest-first log. Those lines are correct and are NOT moved (append-only); this note records that they are at the bottom, not the top.
- **Measured drift (git, not memory):** at `origin/main` 36d5c3d5 the reconciler reported **35 tickets across 42 merges** landed since the `f3688d5` stamp and unrecorded, spanning 2026-09-03 → 2026-09-17. All 35 are now recorded; each line above is sourced from that ticket's own end-of-session comment plus its merged diff.
- **Four landings carry NO end-of-session comment** — RYA-1213, RYA-1214, RYA-1183 and RYA-935. Their lines are marked and sourced from the merged diff alone, and they are listed in the RYA-1184 adjudication block for Ryan's confirmation rather than treated as settled.
- **RYA-850** (merged 2026-08-17, 953b98d0, PR #288) is recorded nowhere in this file or the register, and it promoted the graded lab-gf pool to the PRIMARY reported value. The register rows are now annotated; this landing predates the `f3688d5` stamp, so the reconciler does not surface it.

## 2026-09-13

- **RYA-1218** — Si preflight and source correction; exposes identity/grade holds before abundance work.
- **RYA-1134** — Provisional 505-row Al adjudication and holding/engine dispositions; RYA-1217 handoff remains held.

## 2026-09-03

- **RYA-1185** — reconciled register/SEQUENCE/LEDGERS to merged reality (stamped `main f3688d5`) and consolidated `Fe.json` v1.87 → v1.92
- **RYA-1185** (Part A) — folded the per-band dA/dxi into every product: **17 ALIASED + 13 KEY_AMBIGUOUS → 0**; MEASURED 19 → 51; no published A moved
- **RYA-1185** (Part A) — carried the 4 RYA-1106 Asplund products into the feed through `publish_product.py --line-set asplund` (the route RYA-1178 named); feed 66 → 70
- **RYA-1185** (Part A) — grades onto the ratified taxonomy: `Our Grade` → **Codex Grade** (42), `Asplund Grade` → **Reference Grade** (4), Deep Grade (24)
- **RYA-1168** — dA/dxi for the four near-UV Fe pools; the missing run behind RYA-1113 F7 — the arm was never skipped, the BAND was
- **RYA-1163** — dA/dxi for the three IR Fe pools; the xi campaign was band-pinned to VIS
- **red-optical xi** — dA/dxi for the last band, the one nobody had a ticket for; completes per-band xi coverage
- **RYA-1178** — re-emitted `Fe.json` on a complete publication schema (line_set, atlas/telluric state, wavelength range, xi terms, IR irreducible)
- **RYA-1173** — built the AGSS21-lineage solar Al reference set from its primaries and enforced RYA-946's census gate; **AGSS21 publishes no Al line list**
- **RYA-1142** (A6) — the A6 status now follows its own AST evidence instead of being hardcoded
- **RYA-1143/1144/1150** — fixed the three CNO intake defects the QA found (wavelength-only N I join, echoed "published" values, verdict drift)
- **RYA-1169** — Si gf-graded line sets + AGSS21 reference census; output deliberately **BLOCKED**, not frozen

## 2026-08-31

- **RYA-1160** — acquired NIST ASD C/N/O gf: 3,611 graded transitions, the first graded gf source CNO has ever had; O I 777 grade was overstated A+ vs A
- **RYA-1141** (E-series) — audited the Al intake against the recipes the ticket actually names

## 2026-08-30

- **RYA-1132** (source ingest) — normalized 106 Vujnovic rows; promoted six Al I plus Johnson Al II 2669 to GF-LAB; UV policy still blocked
- **RYA-1132** (web follow-up) — found UV lab sources; UV moves to crossmatch review, IR has no new source and remains pipeline-blocked
- **RYA-1132** — froze 505 Al I/II candidates/controls; VIS documented, UV atomic-data-blocked, IR pipeline-blocked; gate stays closed

## 2026-08-25

- **RYA-1044** — the Engine-B band-product leg was UNREACHABLE since it was "wired" (2 tickets Done, 0 products); now runs, and its first product shows +0.141 of the +0.187 was ATMOSPHERE
- **RYA-1040** — ⟨3D⟩ treatment pair as AXES + route wired; the run found the deck was read **TRANSPOSED** and no ⟨3D⟩ deck could synthesise at all
- **RYA-1013** — staged + verified the **O ⟨3D⟩ deck**: its CONSUME leg works (via Gerber, not the SIU route it names); build becomes a referee, not the path
- **RYA-1035** (follow-up) — swept the ⟨3D⟩ aux defect across all 17 published decks: **Mn has Fe's identical defect**; the other 15 are clean
- **RYA-1029** — the ⟨3D⟩ atmosphere loader; the "still owed PHYSICS step" is an interface contract, and iSpec's τ overwrite is an identity here (5e-5)

## 2026-08-24

- **RYA-710** — Fe ⟨3D⟩ **WIRED**: `Fe@mean3D` registered against the PLAIN aux (opposite of Al); cell `T2_FETCH_OWED` → `T2_CONSUME_WIRED`
- **RYA-1035** — Step 0 came back **HAVE**: the Fe ⟨3D⟩ deck was on our own Keeper share all along; 13 matrix cells were mis-stated as builds ⇒ new `T2_FETCH_OWED`
- **RYA-821** (follow-up) — 🔴 **Engine-B NLTE could not run Al AT ALL, on either deck**: the call site hoisted departures out of the χ² loop with no abundance, which is right for Fe (1 value) and REFUSED for Al (31 values, b differing by up to **10.28** between adjacent 0.1-dex nodes). The matrix blamed 'a registry line' — wrong; it was the call site. Now decided **per deck**, evaluated per trial, with `departures_at_abundance` interpolating linearly between brackets (exact at a node, verified byte-identical). **Unblocks Al 1D-NLTE as well as ⟨3D⟩.** Still owed: the ⟨3D⟩ band product needs a `gerber-mean3d` deck choice AND the ⟨3D⟩ atmosphere, not MARCS.GES

- **RYA-821** — **the Al ⟨3D⟩ consume route RUNS**: the deck is **read directly**, bypassing the vendor interpolator, which structurally cannot consume a ⟨3D⟩ deck (needs τ_Ross+P_g; the public STAGGER archives ship neither). Layout verified twice — field sizes and the aux's own pointer stride both give 287,348 B. n_dep=101, log τ −5…+5, median b→1.0015 at depth. 🔴 **Record and aux use DIFFERENT naming conventions — compared on PHYSICS, never strings.** 🔴 **The matrix's `_MEAN3D_WIRED` was hard-coded and LIED about main for days; now DERIVED from `gerber_nlte.DECKS`.** Al: BROKEN → T2_CONSUME_WIRED. Still owed: an end-to-end band product

- **IAG split** — `solar_iag` (Baker+2020, corrected) and `solar_iag_reiners2016` (the BLUE arm) are separate holdings. 🔴 **The 'IAG is uncorrected' doubt was a STALE CATALOGUE ROW, never the reader** — served flux measures 0.18%/0.00% in O₂/H₂O vs 46.25/51.63 raw, so RYA-1026's CLEAN_WITH_ANOMALY was WRONG and is retired. 🔴 **Not a data gap**: Baker starts at 5001.1 Å, the blue 954 Å is Reiners-only and holds **136 primary-lab Fe I lines (29%)**, and there is nothing to correct there (measured vs KP2005 residual: 10.48 vs 10.35% below 0.5). `solar_iag` had declared NO span (RYA-767 shape); both do now, in AIR. No posted value moves — zero Fe IAG products exist

- **RYA-1027** — Fe I **near-UV graded product** on telluric-corrected KP: **A(Fe I) = 7.596 (n=58 deep-graded), ENGINE-A 7.606** vs RYA-759/822's 7.577 ⇒ **+0.019, values HOLD**. 🔴 **The at-gate graded tier is n=1 (3026 Å, continuum-limited, 8.529) — a DECLARED GAP, not a value**; RYA-954 reproduced at band scale, lab Fe I is a deep/saturated population. 🔴 **Wide bar is a RESULT (RYA-777)**: MAD 0.272, and it is FIT-QUALITY limited (red_chi2 47–231 on every line), **not** rail-limited (edge_distance ≥ 1.90). No χ²ᵣ cut — RYA-847 measured that no threshold transfers

- **RYA-1030** — normalisation is **determined from the FLUX at intake**, never asserted from a label; `normalization_state` joins `telluric_applied` and `observed_conditioning` as the third conditioning axis. 🔴 **A declared flag and a MIS-ROUTED file agree perfectly and are both wrong** — the KP2005 class (0.022 dex, for months) and the still-open IAG one. Detector = upper envelope; **both level and slope tests are load-bearing, neither sufficient**. 🔴 **The first tolerance (0.05) was WRONG** — the known-normalised KP1984 atlas failed its own test below 4200 Å; the fix was to establish the DOMAIN (blue edge 4500 = the existing RYA-451/460 class, telluric bands, fill values) and only then re-derive **0.15** from the measured gap. 55/55 validated. Wired into preflight check 7 + the data-audit skill. LOUD STOP, never an auto-fix

## 2026-08-23

- **RYA-1026** — the Fe product is ratified as an **(instrument × band) GRID**; **displayed science is telluric-CORRECTED only** (KP2005-vs-KP1984 whitelisted as the molecfit CONTROL, `as_control=True` stated in code) and the **whole Kitt Peak class is PRE-NORMALISED**, reversing RYA-929 for KP2005. 🔴 **DO NOT NORMALISE ANY KP ATLAS** — a second continuum tilts the spectrum and it has bitten TWICE (RYA-940 on 1984; the 2005 double-normalise forced the VIS re-run). Both rules are GUARDS that raise (`prenormalised_guard`, `telluric_display_policy`), not notes — RYA-929's flag said False for months with nothing objecting. 🔴 **The display gate is STRICTER than `telluric_policy.gate_holding`**: passing the MEASUREMENT gate is not permission to SHIP. 🔴 **The clean set is DERIVED from `applied_state`** — a hand-written one got 2 of 5 wrong and invented the id `solar_delbouille`. No value moves. Unblocks RYA-959/961/908/953/1027/1028 and the 3D umbrella RYA-1029

- **RYA-1005** — Al's Gerber NLTE deck wired; **its grid has an abundance axis Fe's does not** (31 vs 1 A(X)), so `for_node` now interpolates at the synthesis abundance. Gate PASSES but the deck reports b~1 and a delta indistinguishable from LTE at the 0.01 dex quantiser floor — recorded, not tuned.

- **RYA-1002** — Al gets a **primary-lab gf table** (Burheim 2023, 12 lines, 2-11 %) and `gf_grades.GF-LAB` is generalised off the Fe-only file — proved inert on Fe (931 verdicts byte-identical). 7/8 in-range Al lines reach GF-LAB on `linelist_solar`; only 1 on `canonical_gf`.
- **RYA-1006** — the artifact stem gains the **observed-spectrum conditioning** axis (`--local-renorm`, `--degrade-to-R`). 🔴 **It had ALREADY overwritten the rya984_graded_163 anchor — KP 7.417→7.337, HARPS 7.535→7.339, provenance BYTE-IDENTICAL.** Summaries restored from git; per-line CSVs unrecoverable, re-run owed; anchor now loud-fails on an unverifiable part
- **RYA-992** (follow-up) — stage-4 synth GoF cut made **PER ARM** via a measured `ARM_SCALE` registry; `instrument` now required. 🔴 **RYA-995's premise is BACKWARDS — the HARPS bar goes UP (1.0 → 1.357), not down.** Unblocks RYA-968 stage 4
- **RYA-1001** — Al Phase-0 census: Burheim 2023 grades **8** of our Al lines (RYA-835 merged the paper's cm^-1 and A columns); `1995JPhB..` is Opacity-Project **theory**, not lab; 6696.185's NIST C+ is off a 0.88 eV-different transition ⇒ VIS graded pool is **zero**. Unblocks Phase 1 on 13123/13150 via CRIRES+.

## 2026-08-22

- **RYA-988** — cited RT **vmac** for tau Ceti (0.9, Bruntt+2010 B1) + eps Eri (1.8, Saar & Osten 1997 T3). 🔴 **eps Eri absent from Bruntt; VF05 publish NO vmac; Doyle+2014 out of range on BOTH axes.** 🔴 **Bruntt B1's 'Saar' columns are TRANSPOSED (4/4).** Adopted as PAIRS — vsini moves too
- **RYA-963** — α Cen A CRIRES+ telluric: **6/6 gated PASS** (0.21→0.02); new holding. 🔴 **`acen_orbit` A/B CONTESTED: CRIRES pair says inverted, NIRPS anchor disagrees**
- **RYA-964** — `resolve_star()` at intake: one alias lookup in `system_catalog.csv`. 🔴 **`STD` / `Star S5` return UNRESOLVED, loudly — the refusal IS the feature.** eps Eri's four OBJECT spellings now collapse to one id
- **RYA-957** — GBS params adopted for **tau Ceti + eps Eri** (Heiter+2015 / Jofré+2014); RYA-952's holdings row now REGISTERS. 🔴 **tau Ceti's log g is BRACKETED by Heiter Table 10 — solved, not pinned.** tau Boo verified NOT GBS
- **RYA-952** — CRIRES inventory by header truth: 123 files → **64 distinct**, **all CRIRES+, zero classic on drive**. 🔴 **tau Ceti was hiding under `OBJECT='STD'`**; tau Cet co-add SNR 639→1238. 🔴 **The co-add cross-correlation is TELLURIC-LOCKED across a large dBERV**

## 2026-08-21

- **RYA-944** — Delbouille/Liège Jungfraujoch **disk-center INTENSITY** atlas acquired (BASS2000 visible arm, 3000–10000 Å, 2,380,014 pts), audited and registered as `solar_delbouille_liege`. 🔴 **The ticket's CRITICAL "VACUUM wavelengths, convert vac→air" is REFUTED — it is AIR as delivered**: 9 line centroids sit +3 to +14 mÅ from AIR rest and 1.1–2.4 Å from vacuum, so the prescribed conversion would have shifted the atlas ~−1.8 Å. **The residual is physics, not error** — mean +371 m/s = gravitational redshift minus convective blueshift, and Hα (chromospheric, no blueshift) returns **+626 m/s vs the textbook +633**. ⚠️ **The spec's own smoke test could not have caught it**: `argmin` returns the nearest GRID POINT, within 0.002 Å of any value in any medium — a grid-spacing test wearing a medium test's label. Tellurics MEASURED not assumed, against KP (retains) and IAG-Baker2020 (corrected) with a clean-continuum control: **present and saturated in O2** (A-band 26.3% below 0.5, reaching 0.000) but **H2O 7–8× cleaner than Kitt Peak** — the 3580 m altitude signature, so `line_selection`, NOT `corrected`. 🔴 **Incidental CRITICAL: `iag_fts_solar_atlas` is catalogued `telluric_basis=corrected` but the manifest routes it to the telluric-RETAINING Reiners+2016 file** (46% of the O2 A-band below 0.5) — so `exclusion()` excludes nothing there; `telluric_policy.py`'s own IAG table quotes Baker+2020's numbers. Logged, not fixed — needs its own ticket
- **RYA-945** — Fe I/II primary-lab gf backbone into `canonical_gf`: 1,597 rows rewritten, diagnostic pool 475 → 1,826, cited σ 0.030/0.041 vs the 0.20 blanket. 🔴 **The values barely move — the payoff is the σ.** `RU` is Raassen & Uylings, NOT Ruffoni
- **RYA-933/934 enabler** — **both telluric-corrected Kitt Peak holdings were UNREACHABLE**: registered but in no loader, so the Fe reruns had no way to name them. 🔴 **The RYA-904 guard did not catch it — `audited` holdings are exempt, it only demands `verified` ones.** Wired all three with holding-level `pre_normalised` (Kurucz 2005 is ABSOLUTE irradiance and ships no continuum) and a declared 3000–10000 Å span so the IR cannot be served from it
- **RYA-940** — **the 1984 Kitt Peak atlas is telluric-corrected**, with NO observation metadata to work from. 🔴 **Three O2 bands 1300 Å apart agree on airmass to 4% (1.303/1.306/1.352)** — the referee that replaces 'column must be 1.0'. O2 A median 0.744→**0.991**, %<0.5 32.83→**0.55**, RMS vs Kurucz 2005 0.500→**0.092**. ⚠️ **H2O column is NOT an airmass.** 🔴 **H2O 7160–7340 is an honest NULL** — best fit rms of any band, zero-width LSF from every start. **11120–11560 has no external referee at all**
- **RYA-939** — **molecfit built on Sirius** (4.4.4-5, matched to the Mac; no ESO Ubuntu repo exists); the telluric leg is no longer Mac-only. **Dual-machine agreement to 5e-5.** 🔴 **Product filenames came from `~/.esorex/esorex.rc`** — a fully successful run reported `failed (rc=0)`. Normalised solar spectra now TRACKED + generator-registered; the RYA-686 guard extended to `data/processed/`
- **RYA-938** — RYA-929's Kurucz line table was a **VACUUM-vs-AIR** artifact: measured lag **+1.700 Å** vs predicted air–vacuum offset **+1.697 Å**, IAG control **−0.000 Å**; Al I 6696 recovers **0.008 → 0.248**. 🔴 **`lm0840` was a saved HTTP 500 page and `kp_segments()` reported it as "no segment covers 8420 Å"** — 251 files, 250 inventoried. Re-fetched; integrity guard added. **1984 is the only arm past 10000 Å, and it is uncorrected there**

## 2026-08-19

- **RYA-925** — Near-UV Al retracted: 3057 Å is abundance-insensitive in a 146-transition blend; appendix evidence retained, no product.
- **RYA-925** — Calibration gate corrected: Kitt Peak VIS/red-optical replicate within combined 1σ; IR/model continuation remains open.
- **RYA-925** — Al Kitt Peak matrix banked across EW/synth LTE + Amarsi; non-Fe identity defects fixed; VIS fails literature validation, no tuning.
## 2026-08-20

- **RYA-931** — the Molecfit "MIPAS" failure was an **empty SCIENCE table** (HARPS `ERR` is 100% NaN); the message named the last file touched, not the cause. Correction now runs on all ten exposures and `solar_harps_molecfit_corrected` is registered alongside an untouched `solar_harps`. **O₂ B min 0.0060 → 0.7725, pct_below_0.5 22.06 → 0.00.** 🔴 **The fit is BISTABLE — a zero-width LSF converges too; acceptance is physical (O₂ column ≈ 1), and only the STARTING POINT is retried**
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

## 2026-09-15

- **RYA-587** — Review checkpoint: shared uncertainty gate and 13 recovered Reference xi stamps; full budget migration remains open.

- **RYA-587** — Reference xi integration complete in RYA-1213 0835200f; reuse existing Fe evidence and advance to other elements.
- **RYA-587 / RYA-1213** — Integrated 13 completed VIS xi results; zero Reference NOT_IN_CAMPAIGN entries, no new synthesis.

### 2026-09-15 — RYA-587 merge scope and CNO priority

Canonical covariance API and new/changed-product admission checks land together. Existing legacy products may be retained exactly without certifying completeness; no Fe synthesis was repeated. Per-product migration remains open: CNO with RYA-1220 nitrogen first, followed by Al and Si.
### 2026-09-15 — RYA-1213 provenance-registry integration repair

Reconciled inherited Al/Si result registrations and moved seven audit-only source declarations beside their actual files. Existing generator checks cover all 1185 tracked artifacts; no measurement rerun or value change. Two historical Al checkpoint aggregations explicitly lack a committed harness.

## 2026-09-18

- **RYA-1224** — the `min_paired = 3` xi floor becomes a property of the RULE, not of the artifacts. It was declared by all four band-keyed runs and by no code, so the RYA-1120 campaign route never applied it and **one physical pool carried two honesty standards** — UNMEASURED at Reference, ALIASED at Deep. **ALIASED 6 → 0 in the published feed** (`Fe.json` v1.219 → v1.220), published A moves **0 of 160**. 🔴 **18 artifacts are published at TWO tiers with byte-equal `provenance.sha256`** — the Fe II VIS Deep/Reference rows are ONE CSV whose filename says `DEEPGRADED`, `line_set_resolved` differing only because RYA-1127 derives it from `tier`; 11 of those 18 pairs had been returning two different xi verdicts, now 18 of 18 return one. A cross-tier derivative is served only on proven identity (hash + full physical key + `A`/`n_lines`/`n_excluded`, and `n_paired == n_lines`), never on matching counts. The ticket's stated 4/2 split of the six was really **3 ENGINE-A (n=2) + 3 1D-LTE (n=8)**; the audit artifact had it right and only the prose erred.
- **RYA-1224 / RYA-587** — two defects kept the corrected layer out of the file. `LEGACY_RIDE_ALONG_STAMPS = ("generated_at", "code_commit")`: `enrich` restamps both on every product every run, and the legacy exemption compared whole dicts, so **all 160 rows lost their exemption to a clock tick** and the feed was unwritable by its own emitter (refusals 160 → 16). Then 🔴 **the gate was preserving a KNOWN DEFECT** — a `sigma`/`xi_state` change is exactly what it must refuse without evidence, so the borrowed derivative stayed live because it was already in the file. Ryan's ruling: a narrow `xi_layer_correction_problems()` route, admitting a re-publication only when the changed set is a SUBSET of `XI_LAYER_FIELDS`, `A`/`n_lines`/`n_excluded` also pass a by-name RYA-161 check, an ASSERTION names an artifact that resolves on disk, a WITHDRAWAL leaves nothing readable, and the sigma arithmetic is recomputed because the route skips `validate`. **Evidence is owed to ASSERT a number, not to WITHDRAW one.** The route fires once: a second emit changes only the two clock stamps. Full RYA-587 budget migration for Fe remains open — `validate` makes `stellar.teff`/`stellar.logg` mandatory, which RYA-1112 found the published `sigma_syst` omits, so honest budgets would move `sigma_reported` on all 160 rows.

### 2026-09-18 — RYA-1224 scope and what it does not close

Sirius: identical failure set to main (6 = 6), zero new failures, +57 passes. `data/results/rya1055/fe2_label_audit.json` regenerated because it re-derives feed sigmas under `--check` — six rows move and **zero `within_reported_bar` verdicts flip**, checked because removing a borrowed sigma term NARROWS the bar. Left open and adjudicated on the ticket, not recorded as a landing: RYA-1168 and RYA-1213 measured dA/dxi **2.0–2.8× apart on a byte-identical artifact** at identical `n_paired` (subsequently root-caused by RYA-1225 to RYA-1207's molecular opacity landing between the two runs, `use_molecules` being near-UV only). Also flagged, not fixed: neither a published product feed nor RYA-587's publication gate is a member of `pipeline/state_surfaces.py`, so the register-freshness gate reports "0 state surfaces changed" for a PR that rewrites both.

## 2026-09-19

- **RYA-1227** — the xi campaign CLOSES: **NOT_IN_CAMPAIGN 11 -> 0 feed-wide** (`Fe.json` v1.220 -> v1.221, MEASURED 113 -> 124), published A moves **0 of 160**. 🔴 **`NOT_IN_CAMPAIGN` was never a disposition, it was an unfinished campaign** — those products published a `sigma_reported` with the xi term ABSENT, not zero, on 1D engines where microturbulence plainly applies while the same engines were MEASURED in other bands. **Root cause measured from the artifacts, not assumed:** every affected product artifact is stamped 2026-09-08/09-10 and every band-keyed xi run predates it (NIR/H 2026-08-30, near-UV and red-optical 2026-09-03, VIS older), so RYA-1208's Gerber fan-out landed AFTER the campaigns and its cells never had a xi leg. Measured on Sirius, 11 units x 2 legs, each on the product's OWN pool (`--lines-tier graded` / `--lines-deep-graded`, never `reference`) and on the current synthesis — the two near-UV Fe II cells post-RYA-1207. 11 of 11 MEASURED, zero incomplete, **zero sign disagreements with the aggregate route**. Every 1D engine (models 1-4 **by `model_id`**, because model 7 shares `atlas9` with models 1 and 3 and an atmosphere test would call full 3D a 1D engine) now reads MEASURED or UNMEASURED and nothing else.
- **RYA-1227 / RYA-1224** — the ticket listed 16 owed cells; **it is 11**. Five were already resolved by RYA-1224's same-artifact route (3 Fe II VIS gerber, CRIRES+ H 1D-LTE, H ENGINE-B-NLTE); the list was taken from `Fe.json` v1.219 and v1.220 already read 11. ⚠️ **Three pools agree with RYA-1213's REFERENCE slope to 4 decimals, and that is pool identity, not borrowing** — the H gerber and both near-UV Fe II gerber line sets are byte-for-byte the same lines at GRADED/DEEPGRADED and REFERENCE (compared by wavelength, not by count), so two independent campaigns on one pool agreeing is a free reproducibility check rather than a defect.

### 2026-09-19 — RYA-1227 scope and what it deliberately did not do

The same legs also measured **1D-LTE and ENGINE-A** on these pools, band-keyed and arguably better evidence than what those cells carry today — **not emitted**, because those cells are already MEASURED and a second derivative would move a published sigma on products the ticket never named. They are recorded in the artifact under `base_treatments_measured_but_not_emitted` so the narrowing is auditable and reversible. `data/results/rya1055/fe2_label_audit.json` regenerated (it re-derives feed sigmas under `--check`): two near-UV Fe II rows move, **zero `within_reported_bar` verdicts flip** — checked because removing or adding a sigma term changes a bar's width and can flip an ionisation-balance agreement. Also recorded for the next runner: **a fresh clone cannot run the near-UV legs** until `pipeline.nearuv_linelist.build()` writes the 12 MB untracked `atomic_lines.tsv`; the first pass of those four legs failed on exactly that.

## 2026-09-19 (later)

- **RYA-1223** — Fe site reconcile, Step 1 VERIFY + Step 2 APPLY against feed **v1.221**. 🔴 **No live abundance was wrong.** All **20** stale rows differ ONLY in `sigma_reported` / `sigma_xi` / `xi_state` — rows where `A` differs: **0**. In scope (Ryan's Option B bands VIS/red-optical/NIR/H): **CORRECT 62, STALE 20, NEW 56, ORPHAN 0**; near-UV HELD for RYA-1226 and rendered "systematic under development" with its current post-opacity abundance (no pre-opacity ~7.9xx ships — verified independently on both feed and live site). Three ticket premises did not survive the artifacts: the live snapshot is **v1.163 not v1.158**; the site carried **NO Reference tier at all**, not "Reference in VIS only" (RYA-1213 merged 09-16, site generated 09-11), so all 56 in-scope Reference products are ADDITIONS; and the corrected bars are **not uniformly larger** — 16 widen, **4 NARROW**, and the 4 are the honest ones because RYA-1224 removed a *borrowed* sigma_xi from sub-floor (n=2) holds.
- **RYA-1223 / RYA-587 / RYA-1120** — the Option B **budget-completeness gate does NOT pass, and not differentially**. Measured three ways: `sigma_syst_components` is `(published_syst, sigma_xi)` on **all 160** products with **no named stellar-parameter term in any band**; the sourced stellar-parameter budget still reaches **no product path** (`solar_uncertainty_rya158` is read only by its two stamping scripts plus three auditors, RYA-1112 F1 standing); and `scripts/audit_uncertainty_rya587.py` reports **"full-contract complete set: []"** with **160/160** products failing on `star: provenance is required`. So the 4 bands sit on the byte-identity retention clause exactly as the near-UV rows do — the near-UV difference is *operational* (RYA-1225 has a pending sigma_xi change, which costs the exemption) and not a budget difference. Ryan ruled PROCEED on 2026-09-19: the reconcile is strictly more honest than a live site carrying 16 too-small bars and 7 xi-less products. Teff/logg remains owed to **RYA-1120**, canonical budgets to **RYA-587**.

### 2026-09-19 — RYA-1223 what it deliberately did not do

Minimal diff, everything generated: the site's own `scripts/generate_fe_publication.py` re-run against v1.221, no page hand-edited and no value typed (RYA-914). The per-band best-value headline selection was left exactly as already instructed — it renders the VIS Reference Amarsi 3D-NLTE anchor and **7.466 is not the headline**. ⚠️ Recorded, not fixed: the near-UV `opacity_note` carried in `Fe.json` on main still asserts "the headline, which rests on the VIS anchor 7.466" — stale prose predating RYA-819's retirement and RYA-850's promotion, now contradicted by the rendered anchor. It is historical provenance text rather than a live product value, so it was reported rather than hand-edited in a publish ticket; the fix belongs to the RYA-1178 emitter path.
