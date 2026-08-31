# RYA-1142 — independent QA of the RYA-1136/1131 CNO intake

**3 PASS · 10 FAIL · 5 FLAG** — findings-only; no intake artifact, atomic manifest, molecular constant or gf value was mutated (see `artifact_integrity.csv`).

## Verdict

The CNO intake is **NOT independently verified**. Its census is real and its arithmetic reproduces exactly, but the gate stays closed on 10 findings, of which three are CRITICAL by the ticket's own list: a wavelength-only admission (A6), an ambiguity-tolerant argmin match (A2), and a molecular redistribution labelled primary (A4). A fourth, a missing checksum on the AGSS21 article, is what makes A1 and A7 unclosable.

The blocked verdict is **not overstated** — nothing here reads freeze-ready, and all four of its stated blockers hold under independent recomputation. It is **understated**: three defects this QA found are absent from it, and one of them falsifies its own safety line.

## Per-check results

| Check | Title | Status |
| --- | --- | --- |
| A1 | AGSS21 count reconciliation | **FLAG** |
| A2 | Physical-identity crossmatch, never wavelength-alone | **FAIL** |
| A2b | identity_basis honesty (wavelength's real role) | **FLAG** |
| A2c | RYA-1037/1033 guard still covers the intake | **FAIL** |
| A3 | Molecular identity completeness | **FAIL** |
| A4 | Provenance / grade earned | **FAIL** |
| A4b | Primary-parser transcription validated physically | **PASS** |
| A5 | Molecular constants provenance | **FAIL** |
| A5b | Hand-set C2 energy-origin constant | **FLAG** |
| A6 | Atomic side (C/N/O census, EP-aware joins) | **FAIL** |
| A6b | Atomic species completeness + the [O I] 6300 blend | **FAIL** |
| A7 | Rejected / negative results retained | **FLAG** |
| A8 | Band coverage across both domains (is this really UV-IR?) | **FLAG** |
| A9 | UV molecular transitions held on disk but never read | **FAIL** |
| B1 | Reproduce 408-row inventory + six-bin coverage | **PASS** |
| B2 | No canonical mutation; RYA-1130 separation intact | **FAIL** |
| B3 | BLOCKED verdict honestly derived | **FAIL** |
| INTEGRITY | This QA mutated nothing | **PASS** |

## Reproduced headline claims (B1)

408 used rows — C2 39 / CH 54 / CO 80 / CN 59 / NH 31 / OH 145; VIS 45 / NIR 122 / IR 241; FUV 0 / NUV 0 / RED_OPTICAL 0.

**Scope, stated plainly (A8/A9): this intake is VIS-to-IR, not UV-to-IR.** Zero FUV and zero NUV rows in BOTH domains, against RYA-1136's title 'UV–IR' and RYA-1131's 'across FUV/NUV/IR'. The atomic census spans 5052–10109 Å. And the UV is not simply unavailable — 29,738 ultraviolet molecular transitions are acquired and unread in `data/reference/cno_molecular_primary/`.

## AGSS21 ↔ Amarsi Table 2 reconciliation (A1)

| Cell | AGSS21 (banked, RYA-1131) | Amarsi Table 2 | Δ |
| --- | --- | --- | --- |
| C2 Swan | 39 | 39 | 0 |
| CH X-X dnu=1 | 51 | 48 | -3 |
| CH A-X dnu=0 | 7 | 6 | -1 |
| 12C16O X-X dnu=1 | 28 | 28 | 0 |
| 12C16O X-X dnu=2 | 52 | 52 | 0 |
| NH X-X dnu=0 | 13 | 13 | 0 |
| NH X-X dnu=1 | 15 | 18 | 3 |
| CN A-X dnu=0 | 59 | 59 | 0 |
| OH X-X dnu=0 | 84 | 84 | 0 |
| OH X-X dnu=1 | 50 | 46 | -4 |
| OH X-X dnu=2 | 15 | 15 | 0 |
| TOTAL | 413 | 408 | -5 |

The `(system, |Δν|)` partition is the right decoder: 7 of 11 cells agree exactly, including both CO cells in the published order. The four residual deltas are localised and are not ours — the transcription re-slices byte-exactly to the CDS ReadMe. They cannot be closed because the AGSS21 article was never acquired.

## Independently recomputed match tally (A2)

| Leg | Variant | unique | ambiguous | none |
| --- | --- | --- | --- | --- |
| primary(non-CO) | BASELINE | 278 | 9 | 41 |
| primary(non-CO) | NULL scrambled published_loggf | 17 | 0 | 311 |
| primary(non-CO) | NULL scrambled lower_energy_eV | 3 | 1 | 324 |
| primary(non-CO) | NULL wavelength displaced +20 cm-1 | 0 | 0 | 328 |
| primary(non-CO) | DROP loggf term | 140 | 184 | 4 |
| primary(non-CO) | DROP lower_energy term | 278 | 9 | 41 |
| primary(non-CO) | DROP vibrational band term | 282 | 9 | 37 |
| primary(non-CO) | DROP wavelength window (band+E+gf only) | 58 | 229 | 41 |
| CO | BASELINE | 80 | 0 | 0 |
| CO | NULL scrambled published_loggf | 0 | 0 | 80 |
| CO | NULL scrambled lower_energy_eV | 0 | 0 | 80 |
| CO | NULL wavelength displaced +0.5 A | 0 | 0 | 80 |
| CO | DROP loggf term | 80 | 0 | 0 |
| CO | DROP lower_energy term | 80 | 0 | 0 |
| CO | DROP wavelength window | 76 | 4 | 0 |

Baseline reproduces the artifact exactly on all four quantities. Every tolerance is qualified by a displaced null: scrambling a single published column collapses the match rate by an order of magnitude or more, and displacing wavelength by 20 cm⁻¹ takes it to zero. **Zero wavelength-only admissions exist in either molecular leg** — the wavelength-only defect is on the *atomic* side (A6).

## Findings

### A1 — AGSS21 count reconciliation · FLAG

The (system,|dnu|) partition is the correct decoder for the AGSS21 sub-counts: 7 of 11 cells agree EXACTLY, including both CO cells in the published order, both C2, CN, NH dnu=0, OH dnu=0 and OH dnu=2. Transcription re-slices byte-exactly to the holding's own ReadMe and the columns the crossmatch never uses (Param, |dnu|, System) all fall inside the ReadMe's stated domains, so the 4 residual deltas are NOT ours. They are localised to four cells and total 408 vs 413. They CANNOT be closed here: the AGSS21 article was never acquired (source_bibliography.csv row AGSS21 is asset='article', sha256 EMPTY, status=SOURCE_IDENTIFIED), so the banked side is ticket prose with no checksummed referent to re-read.

Offending rows: `CH X-X dnu=1 banked=51 table2=48; CH A-X dnu=0 banked=7 table2=6; NH X-X dnu=1 banked=15 table2=18; OH X-X dnu=1 banked=50 table2=46`

### A2 — Physical-identity crossmatch, never wavelength-alone · FAIL

REPRODUCED EXACTLY: my independently re-implemented decision logic returns 278 unique / 9 ambiguous / 41 unresolved on the 328 non-CO rows and 80/80 unique on CO -- identical to the artifact on all four quantities. The three-field conjunction is REAL and the tolerances are EARNED: scrambling published_loggf collapses 278 unique to 17, scrambling lower_energy_eV to 3, and displacing wavelength by 20 cm-1 to 0. CO behaves the same way (80 -> 0 on scrambled gf). No wavelength-only admission exists in either molecular leg. BUT the acceptance set ACCEPTED_MOLECULAR_JOINS admits PRIMARY_UNRESOLVED_SUM_MATCH, and 26 of those 32 rows had MORE THAN ONE subset of primary components whose gf sum reproduces the published loggf (up to 16 viable subsets). build_cno_intake_rya1136.py resolves them with min(subsets, key=...) -- an argmin over candidate identities -- and then counts the result as matched coverage. That is the ambiguity-tolerant match the gate forbids: the matcher found A combination, not THE combination, and nothing downstream records that the identity was picked rather than determined.

Offending rows: `1;2;3;5;6;8;9;10;11;13;15;17;19;21;24;29;30;31;32;33;34;36;39;54;65;66`

### A2b — identity_basis honesty (wavelength's real role) · FLAG

ingest_cno_molecular_primary_rya1136.py comments that wavelength 'supplies candidates but can never decide identity'. Measured: removing the wavelength window while keeping band+lower_energy+loggf drops unique matches from 278 to 58 and raises ambiguity from 9 to 229. Wavelength IS what separates 220 of the 278 accepted identities -- expected, since neighbouring J within one band carry near-equal loggf. The join is still a four-field conjunction and the nulls above show it is not wavelength-alone, so this is a documentation defect, not an admission defect. Separately, dropping the lower_energy term changes NOTHING (278 unique either way): E_low excludes zero candidates and acts as a corroborating field, not a discriminating one. Both facts belong in identity_basis; neither is recorded.

### A2c — RYA-1037/1033 guard still covers the intake · FAIL

The guard exits clean and names NEITHER RYA-1136 script (unexpectedly reports one). The reason is structural, not incidental: _enclosing_has_ep() tests the WHOLE enclosing FunctionDef for an EP-like name, so one EP mention anywhere in a function launders every wavelength-only comparison in it (function-scoped test confirmed in the AST). Positive control, run here: an identical lambda-only comparison is FLAGGED with no `ep` in scope (True) and SILENT once an unrelated `ep` is bound earlier in the same function (False). That is exactly the shape of build_cno_intake_rya1136.py:atomic_census(), where the C I/O I loop binds `ep` and the N I loop below it joins on wavelength alone -- see A6. The guard did not pass the intake; it never looked at it.

### A3 — Molecular identity completeness · FAIL

390 rows are treated as matched while 9 of 14 identity fields are absent from every artifact: isotopologue, J_upper, J_lower, branch_or_parity, E_upper, air_or_vacuum_frame, native_intensity_quantity, conversion_provenance, component_vs_band_normalisation. Two are not source limitations but losses in our own code: (1) J'' is parsed for all primary transitions and discarded -- it survives only as a factor inside gf; (2) system and vibrational band exist in primary_molecular_crossmatch.csv and are dropped when the builder merges into molecular_physical_crossmatch.csv. Two 'transition label' columns are not labels at all: every CO row carries an unsplit remainder of the Turbospectrum record, and every CH row carries branch + the o-c residual. Table 2 omitting rotational identity is the honest blocker RYA-1136 names; it does not explain discarding the identity the PRIMARY side does publish.

Offending rows: `isotopologue;J_upper;J_lower;branch_or_parity;E_upper;air_or_vacuum_frame;native_intensity_quantity;conversion_provenance;component_vs_band_normalisation`

### A4 — Provenance / grade earned · FAIL

Two CRITICAL provenance defects. (1) REDISTRIBUTION LABELLED PRIMARY: source_bibliography.csv row Li2015_CO cites 'Li et al. 2015, ApJS 216, 15', role '12C16O wavelengths, energies, transition probabilities', status ACQUIRED -- and points at data/linelists/molecular/turbospectrum/CO/CO_IR_Li2015.dat. That file's own second line reads 'ExoMol Li2015', and the repo's MOLECULAR_MANIFEST.json describes it as an 'RYA-236 conversion of the ExoMol Li2015 CO list to the Turbospectrum babsma .dat format (conversion script external to this repo)'. It is an ExoMol redistribution, twice derived, by a converter that is not in the repo and cannot be re-run -- and it is the sole source behind ALL 80 CO PHYSICAL_TUPLE_MATCH rows, the only clean-match class in the intake. No Li 2015 primary table was ever acquired; data/reference/cno_molecular_primary/ has no CO directory. (2) MISSING CHECKSUM: the AGSS21 row carries asset='article', an EMPTY sha256 and status SOURCE_IDENTIFIED -- the paper whose lineage the whole intake claims was never acquired, which is also what blocks A1 and A7. Everything else verifies: all 9 other assets exist, all 9 recorded checksums recompute EXACTLY, and the CDS ReadMes for CN/CH/Barklem plus the NH reprint's own front matter confirm their bibcodes and DOIs. Also FLAG: 4 acquired C2 supporting archives (ChenEtAl-C2-2015-JChemPhys.zip, ChenEtAl-C2-2016-JChemPhys.zip, RamEtAl-C2-2014-AstroJ.zip, TanabashiEtAl-C2-2007-AstroJSuppl.zip) appear in no bibliography row, with no checksum and no stated role.

Offending rows: `AGSS21;Li2015_CO`

### A5 — Molecular constants provenance · FAIL

molecular_constants_ledger.csv asserts for all six molecules partition_function_source='Barklem & Collet 2016' and verdict='PRIMARY_TABLES_ACQUIRED'. The partition-function table was never acquired. Barklem & Collet publish partition functions in table6.dat and equilibrium constants in table7.dat, per their own ReadMe File Summary; the holding contains ReadMe, table1.dat and list.dat only -- and list.dat is not constants, it is the LIST OF FILENAMES in the table2/ subdirectory, which is likewise absent. So the ledger claims acquisition of two tables that are not on disk and one (table2/*) it never names. What IS genuinely acquired and verifiable is the dissociation energy: table1.dat carries an adopted De for all six molecules (C2 6.371000 eV, CH 3.469600 eV, CN 7.737000 eV, CO 11.117000 eV, NH 3.419000 eV, OH 4.417100 eV) -- and the ledger records none of them. All six ledger rows are byte-identical boilerplate holding no per-molecule value, no table id, no row reference and no checksum, and every row asserts the same isotopic assumption for molecules whose isotopologue the intake never determined. This is an unsourced-constant finding of the kind A5 exists to catch, and the six rows overstate it as PRIMARY_TABLES_ACQUIRED.

Offending rows: `C2;CH;CN;NH;OH;CO`

### A5b — Hand-set C2 energy-origin constant · FLAG

ingest_cno_molecular_primary_rya1136.py adds a hard-coded C2_LOWER_ORIGIN_EV = 0.0753 eV (607.3 cm-1) to every C2 lower energy before the identity test, with no citation and no bibliography row. Swept here: the constant is genuinely data-determined, not arbitrary -- it sits on a plateau of 580.7-645.2 cm-1 (64.5 cm-1 wide, set by the 0.005 eV energy tolerance), the code value lies inside it, and the null is clean: 881 of 1000 sampled offsets yield ZERO matches, with nothing at all outside 565.4-660.6 cm-1. So this is not a free parameter quietly absorbing error. Two things are still wrong with it. It is FITTED rather than sourced -- it was chosen by making the join succeed, which is the shape RYA-161 forbids, and A5 requires typed provenance for exactly this class of constant. And the code's justification, that the shift is 'independently visible across all 39 rows', is overstated: only 15 of the 39 C2 rows ever reach a unique match at the best offset; the other 24 are the sum-matched and strength-mismatched rows, which cannot witness it. Cite it, or derive it from the acquired Brooke/Chen holdings and record the derivation.

### A6 — Atomic side (C/N/O census, EP-aware joins) · FAIL

The C I / O I leg is sound: 17 C I and 26 O I rows from Amarsi 2019 Table 1 all route through nearest_canonical(), which requires wavelength AND excitation potential AND loggf and returns AMBIGUOUS or ABSENT rather than guessing -- 43 rows, EP-aware, correctly refusing. The N I leg is not. The five-line AGSS21 adopted set is selected by `abs(float(r['wavelength_air_A']) - wavelength) <= .05` with NO EP and NO loggf term -- a wavelength-only key, confirmed in the builder's AST (1 wavelength-only comparison(s) vs 0 EP-aware in atomic_census()). 4 of the 5 N I lines are then stamped join_status=PHYSICAL_TUPLE_MATCH, which is exactly the RYA-1034 defect: a lab-tier identity claimed on a wavelength match alone. It is worse than a silent one, because the row then REPORTS lower_EP_eV and published_loggf that were READ OUT OF our own canonical_gf row -- our value round-tripping back as though the primary paper supplied it (RYA-1035's vendor-echo defect), under a column literally named published_loggf. Consequently intake_verdict.json's safety line, 'No abundance derived; no gf tuned; no wavelength-only join admitted', is FALSE on its third clause. The 5th line (10108.90 A) is honestly ABSENT and the N I gap is not silently filled, which is the one thing this leg gets right.

Offending rows: `7442.29A;8216.33A;8629.23A;8683.40A`

### A6b — Atomic species completeness + the [O I] 6300 blend · FAIL

Two sub-clauses of the spec that the wavelength-only finding above must not overshadow. (1) SPECIES: the census is neutrals only -- C I 17, N I 5, O I 26, zero C II / N II / O II. That is honest in origin, not a filter bug: Amarsi 2019 Table 1 contains only CI 17, FeII 142, OI 26, and the FeII rows are rightly excluded as out of scope. The defect is that NOTHING RECORDS IT -- no artifact distinguishes 'the source has no ionised CNO' from 'we did not look', and the ticket asks for a C II / N II / O II census. (2) BLEND: [O I] 6300.300 and 6363.770 are both present and correctly carried, but atomic_source_census.csv has no blend or component column and no Ni I row, so the Ni I blend at 6300 A -- the best-known contaminant of the single most-used solar oxygen diagnostic -- is not retained as a physical component anywhere.

Offending rows: `C II;N II;O II;Ni I @ 6300.300`

### A7 — Rejected / negative results retained · FLAG

All four negative results are RETAINED, not dropped, each with a species, a system, a wavelength region and a stated reason -- the 463 rejected CN A-X red transitions (present) and the three considered-and-rejected UV systems (NH A-X ~340 nm, OH A-X ~320 nm, CN B-X ~390 nm), all three correctly marked count=NOT_PUBLISHED rather than invented. That is the right shape. What cannot be closed: every row's evidence field cites 'Amarsi2021 Sect. 2.1 lines 150-160' / '161-165' -- line numbers into an article body that is not in the repo, so no reader can re-derive the reason or confirm the 463. The 463 does not reconcile against anything held either: Table 2 publishes only the 59 USED CN lines, so the 522 considered total appears in no acquired asset. The ledger is honest but unverifiable, and it should say so.

### A8 — Band coverage across both domains (is this really UV-IR?) · FLAG

The shipped combined_coverage_matrix.csv reproduces exactly in every cell of both domains. The scope claim does not. RYA-1136 is titled 'UV-IR' and RYA-1131 'across FUV/NUV/IR', and the delivered inventory contains 0 UV rows -- ZERO FUV and ZERO NUV, in BOTH domains. Molecular is VIS 45 / NIR 122 / IR 241; atomic is VIS 14 / RED_OPTICAL 33 / NIR 1, spanning only 5052-10109 A. The intake is VIS-to-IR. On the molecular side the UV emptiness is at least CHARACTERISED (Table 2 publishes no indicator below 400 nm, and the three UV systems sit in the rejected ledger) -- though see A9, which shows that characterisation is wrong about availability. On the ATOMIC side it is not characterised at all: rejected_indicator_ledger.csv has four rows and every one is molecular, so nothing anywhere records why a C/N/O census carries no ultraviolet line. C I, N I and O I all have strong solar UV resonance lines; their absence here is a property of the one source table chosen (Amarsi 2019 Table 1), and that is exactly what should be written down rather than left as an empty bin.

### A9 — UV molecular transitions held on disk but never read · FAIL

rejected_indicator_ledger.csv records NH A-X (~340 nm), OH A-X (~320 nm) and CN B-X (~390 nm) as REJECTED with reason 'crowding and continuum/blend limitations; individual list not published'. That conflates two different things. What is unpublished is which subset AMARSI used. The TRANSITIONS are in this repo, acquired and unread: 53231 ultraviolet transitions across those three systems sit in data/reference/cno_molecular_primary/ right now. nh_brooke2014/NH-A-X-linelist.csv and oh_brooke2016/OH-A-X-linelist-final.csv are never opened -- the ingest reads only the X-X members of the sibling archives -- and the CN B-X violet transitions ARE parsed out of table4.dat.gz and then silently dropped, because the index is keyed on (species, system) and no Amarsi row carries a B-X key. Worse for RYA-1148: NH-A-X-linelist.csv publishes J', J", symmetry, branch, v', v", N', N", E_upper, E_lower, f-value AND A -- richer rotational identity than any list the intake does parse, and rotational identity is the intake's own stated blocker. A negative result must say WHICH thing is missing; 'not published' reads as 'not available', and the data is on our disk.

Offending rows: `NH A-X;OH A-X;CN B-X`

### B2 — No canonical mutation; RYA-1130 separation intact · FAIL

MUTATION: clean. Across all five RYA-1136 commits (c314879, 74fee13, ffe67f2, ec4c480, bd9dd08) not one path under data/linelists/, data/audit/rya1129_atomic_intake/ or any canonical_gf file is touched. No transition, gf, grade, manifest row or molecular constant was mutated by the intake, and this QA mutated nothing either (see artifact_integrity.csv). SEPARATION: NOT intact. data/linelists/canonical_gf.csv holds 7800 MOLECULAR rows inside the atomic canonical store -- C2 2654, CH 2265, CN 2766, NH 114, OH 1 -- every one of them seeded linelist(VALD3), loggf_reference VALD3, gf_tier VALD3. RYA-1130 exists to keep molecular transition provenance out of atomic canonical_gf, and it is not being kept out. RYA-1136 did NOT introduce these rows and is not the culprit; but the intake never checked the invariant it was written under, and the consequence is live: the next CNO join that reaches for canonical_gf will find 7,800 VALD3 molecular rows waiting -- exactly the redistribution this intake spent five archives avoiding. Filed separately so it is not lost with this ticket.

Offending rows: `canonical_gf.csv:7800 molecular rows`

### B3 — BLOCKED verdict honestly derived · FAIL

Direction first, because it matters: the verdict does NOT overstate freeze-readiness. frozen_ready_for_measurement is false, summary.json frozen_ready is false, no abundance is derived, and all four blocking_findings HOLD -- I reproduced every count behind them independently. The verdict is nonetheless not honest yet, in three ways. (1) STRING: the ticket asks whether 'BLOCKED_MOLECULAR_DATA' is honestly derived. No artifact contains that string. intake_verdict.json says 'INTAKE_COMPLETE_REVIEW_REQUIRED' and summary.json says 'CROSSMATCH_REVIEW'. (2) UNDERSTATED: the four blockers omit the three defects this QA found -- 26 argmin-resolved identities counted as matched, an ExoMol redistribution standing in for the CO primary, and a wavelength-only N I admission -- and the last of these makes the verdict's safety line, 'no wavelength-only join admitted', FALSE as written. (3) DRIFT: two artifacts in one directory disagree about the same quantity. summary.json reports canonical_matched=0 and crossmatch_review=408; intake_verdict.json reports 390 of 408 in accepted join classes. summary.json is written by the ingest script and intake_verdict.json by the builder, and nothing reconciles them (RYA-1091). A reader who opens summary.json gets a number that is 390 rows stale.

## What each molecule × band still needs to reach FROZEN_READY

| Molecule | Bands used | Blocker | What would close it |
| --- | --- | --- | --- |
| **C₂** | VIS 39 | 23 of 39 resolve only as gf-summed subsets, and ALL 23 were chosen by argmin across up to 16 viable subsets; the 0.0753 eV lower-energy origin shift is fitted, not cited | Amarsi's per-line rotational identity (J″, branch) for the Swan lines, and a citation for the energy-origin offset — measured here to lie on a 580.7–645.2 cm⁻¹ plateau with a clean null, so it is data-determined but unsourced |
| **CH** | NIR/IR 54 | 9 sum-matched (3 of them argmin-chosen), 2 strength-mismatched, 1 ambiguous | J″ from Table 2; the primary side already publishes J, N and parity at bytes 55–95 and we discard them |
| **CN** | NIR 59 | 1 strength mismatch; the 463 rejected red transitions are unverifiable | Acquire the Amarsi 2021 article so the rejection reason and count have a referent |
| **NH** | NIR 31 | 7 ambiguous Λ-doublet pairs the matcher correctly refuses | J″ and parity from Table 2 — nothing else will separate a doublet whose components differ in the 4th decimal of log gf |
| **OH** | NIR/IR 145 | 3 unmatched, 1 energy-mismatched, 1 ambiguous | Reconcile the 4 rows against the acquired Brooke 2016 release; they may be a release-version difference |
| **CO** | IR 80 | **Provenance, not matching.** All 80 join uniquely and survive every null — against an ExoMol→Turbospectrum conversion whose converter is not in the repo | Acquire the Li et al. 2015 ApJS 216, 15 primary tables and re-join against them |
| *all molecular* | red-optical | 0 rows — a genuine published negative (Table 2 lists no molecular indicator between 700 and 1000 nm) | nothing; record it as a negative result |
| **UV — both domains** | FUV / NUV | **0 rows, and NOT the published negative the intake records.** We hold 29,738 unread NUV transitions: NH A-X 6,653 and OH A-X 586 in files never opened, CN B-X 22,499 parsed then dropped on a system key. On the atomic side the UV emptiness is not characterised at all | Read the three held UV lists, and state the negative precisely: Amarsi's UV *selection* is unpublished, the *transitions* are in hand. Record why the atomic census carries no UV line |

## Method and its limits

The decision logic was re-implemented from the published quantities and reproduces the artifact exactly; the vendor parsers were imported rather than re-transcribed, and are instead validated against each holding's own byte-by-byte ReadMe and against `E_up − E_low == 1/λ_vac`, an identity the source never tabulates, with a negative control (A4b). A defect that lives inside a parser *and* inside that physical identity simultaneously would survive this audit; nothing else in the crossmatch would.

This auditor excludes itself by name from every scan it runs (RYA-1116).

