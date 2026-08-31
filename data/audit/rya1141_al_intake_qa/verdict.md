# RYA-1141 - independent QA of the RYA-1132 Al atomic-data intake

**Overall: FAIL.** The measurement gate stays CLOSED. Nothing in the intake was mutated by this audit - every defect below is reported, not fixed (RYA-161).

## Per-check result

| check | title | status |
| --- | --- | --- |
| A1 | Source coverage / truncation (Vujnovic CDS tables 2-5) | **PASS** |
| A1-burheim | Burheim Table 3 = 12 derived log gf; Table 2 does not leak in | **PASS** |
| A1-parse | Fixed-width column extraction, refereed by branching closure | **PASS** |
| A1-flags | CDS limit / note flag columns preserved | **FAIL** |
| A2-control | The identity-comparison test can return True | **PASS** |
| A2 | Crossmatch identity (EP-aware, never wavelength alone) | **FAIL** |
| A2-null | Does the missing identity gate actually mis-assign? | **FAIL** |
| A2-resolution | Would an EP gate alone have been sufficient? | **FLAG** |
| A2-6696 | 6696.015 vs 6696.185 stay distinct; Burheim cannot leak | **PASS** |
| A2-classes | Ambiguity is named, not folded into absence | **FLAG** |
| A3 | HFS component sums independently re-summed and verified | **FAIL** |
| A3-meta | `hfs_n_components` re-verified against the actual component count | **FAIL** |
| A3-rya1001 | The RYA-1001 hfs_n_components defect (3944.006 / 3961.520) | **FAIL** |
| A3-inflation | No HFS/isotope normalisation inflation, no max-component substitution | **PASS** |
| A3-11254 | 11254.9 strong-component vs unresolved-total caveat carried | **PASS** |
| A4 | No air<->vacuum and no cm-1<->Angstrom conflation in the crossmatch | **PASS** |
| A4-medium | Every wavelength carries an explicit medium | **FLAG** |
| A5-lab | No theory graded as primary laboratory | **PASS** |
| A5-bf | No branching-fraction-only row is graded GF-LAB | **PASS** |
| A5-overwrite | No lab/evaluated value silently overwritten by Kurucz/VALD | **PASS** |
| A5-doi-control | The DOI referee accepts correct identifiers | **PASS** |
| A5-doi | Every DOI resolves and names the paper it claims | **FAIL** |
| A5-bibcode | Every ADS bibcode is a resolvable 19-character bibcode | **FLAG** |
| A5-johnson | Johnson 1986 Al II 2669 states the value the intake claims | **PASS** |
| A5-sigma | sigma(log gf) conversions justified with recorded provenance | **FLAG** |
| A6 | Competing gf values retained, never silently dropped | **FAIL** |
| B1 | Inventory reproduced independently from the sources | **PASS** |
| B1-promotions | The 6 Al I promotions + Johnson Al II 2669, and the rejections | **PASS** |
| B1-yield | What the seven new GF-LAB promotions actually unblock | **FLAG** |
| B2 | canonical_gf.csv not mutated by RYA-1132 | **PASS** |
| B3 | Band verdict strings are the ones claimed | **PASS** |
| B3-eligibility | 'Eligible' is derived from the evidence the manifest carries | **FLAG** |
| B3-vocab | `gf_grade` expresses a gf grade | **FLAG** |
| C | Registered holdings reach the coverage module | **FAIL** |
| C-lines | No reachable Al line is reported unreachable | **FAIL** |
| C-bands | RYA-1132's band() covers every band the census does | **FAIL** |
| C-alIII | Nothing is dropped from the census without being recorded | **FLAG** |
| NO-MUTATION | No intake artifact was modified by this QA | **PASS** |

## Detail

### A1 - Source coverage / truncation (Vujnovic CDS tables 2-5): **PASS**

CDS ReadMe declares {2: 29, 3: 22, 4: 24, 5: 31}; files on disk hold {2: 29, 3: 22, 4: 24, 5: 31}; the normalized ledger holds {2: 29, 3: 22, 4: 24, 5: 31}. Total 106 of the claimed 106. The ReadMe's own record counts are an independent referee - it is not the builder's own arithmetic.

### A1-burheim - Burheim Table 3 = 12 derived log gf; Table 2 does not leak in: **PASS**

12 rows, every one carrying both `loggf` and `e_loggf_dex`. No branching-fraction-only row is present, so Table 2 did not leak in as gf.

### A1-parse - Fixed-width column extraction, refereed by branching closure: **PASS**

For every multiplet with more than one finite Aki, A_i/sum(A) reproduces the SEPARATELY PRINTED branching ratio to 0.0031 (n=7), consistent with the source's 2-3 printed figures. The builder never compares these two columns, so this is an identity it cannot have been tuned to.

### A1-flags - CDS limit / note flag columns preserved: **FAIL**

14 flag bytes across 5 documented CDS columns are never read by the parser, so the reference README's claim that 'source limits remain limits' is false for two of them. `l_e_Aki` ('>') turns a LOWER LIMIT on the uncertainty into a determinate sigma, and `n_Lambda` ('*') - which the CDS ReadMe documents as 'the value ... was taken over from Tayal & Hibbert (1984)' - is the only thing distinguishing a theoretical Aki from a Vujnovic measurement, and it is dropped.

### A2-control - The identity-comparison test can return True: **PASS**

The same AST test fires on `nearest`, the builder's own EP-aware matcher (`(frame[epcol] - ep).abs() <= eptol`), so a negative on the ingest path is a real absence and not a broken detector.

### A2 - Crossmatch identity (EP-aware, never wavelength alone): **FAIL**

RYA-1132's `ingest_new_lab_sources` joins every Vujnovic row and the Johnson Al II row to the manifest on `abs(wavelength_air - lambda) <= 0.08` and nothing else. No excitation potential and no level designation appears in any comparison in that function - the level strings it does touch are only WRITTEN into `upper_lower_level_identity`. So no promotion in this ticket was matched on a physical identity. The builder defines an EP-aware matcher, `nearest(..., epcol=..., eptol=0.02)`, and the census loop calls it; the promotion path does not.

### A2-null - Does the missing identity gate actually mis-assign?: **FAIL**

Yes - it is not merely unexercised. 1 manifest row(s) are claimed by two physically DIFFERENT transitions. `alphys_II_3587.0720_0333` is claimed by three source rows: 4f ^3^F^o^_3_ - 3d ^3^D_2_ at 3587.068 A (twice, its true identity) and 4f ^3^F^o^_2_ - 3d ^3^D_3_ at 3587.100 A, a different transition 0.028 A away that the 0.08 A window swallows. That row is MATCHED_NOT_PROMOTED, so no gf is corrupted today - but it is the same code path, with the same window, that wrote all seven GF-LAB values.

### A2-resolution - Would an EP gate alone have been sufficient?: **FLAG**

Not on its own. `lower_EP` is stored at 4 decimal places, and 2 adjacent lower-level pairs in the crossmatch sit closer than 5e-4 eV - the Al II 3d ^3^D_1,2,3_ term spans 0.0003 eV in total, so at the manifest's stored precision EP CANNOT separate the very levels the 3587 collision confuses. The remedy RYA-1034 needs here is the level/J designation plus a uniqueness requirement, not an EP tolerance alone; an EP-only gate would have passed this defect.

### A2-6696 - 6696.015 vs 6696.185 stay distinct; Burheim cannot leak: **PASS**

Two separate manifest rows, dEP = 0.8788 eV. Burheim's laboratory gf sits only on 6696.015 and 6696.185 keeps its fallback source. Note this pair is separated by RYA-1001's census, which IS EP-aware - not by the RYA-1132 ingest path, which would not have distinguished them.

### A2-classes - Ambiguity is named, not folded into absence: **FLAG**

`NO_UNIQUE_MANIFEST_MATCH` covers two different states: 38 rows with NO candidate and 1 with MORE THAN ONE. A reader cannot tell an absent line from an unresolved one. Three classes are needed, not two - the RYA-1072 lesson.

### A3 - HFS component sums independently re-summed and verified: **FAIL**

`HFS_status` is set to the string 'COMPONENT_SUM_VERIFIED' whenever `hfs_n_components > 1` and to 'NO_SPLIT_COMPONENTS_IN_CENSUS' otherwise. No component sum is computed anywhere in the builder, and `component_or_total` is the unconditional constant 'TOTAL_TRANSITION_GF'. The status is ASSERTED, not verified, and RYA-1132's own test (`m[m.HFS_status=='COMPONENT_SUM_VERIFIED'].component_or_total.eq(...)`) compares two constants set three lines apart in the same function.

### A3-meta - `hfs_n_components` re-verified against the actual component count: **FAIL**

The manifest does not carry `hfs_n_components` at all. The metadata RYA-1141 asks to re-verify was dropped at the write, so no reader of the frozen artifact can check it.

### A3-rya1001 - The RYA-1001 hfs_n_components defect (3944.006 / 3961.520): **FAIL**

2 of 2 rows still carry the wrong component count in `canonical_gf.csv` on main today (recorded as 1 while the VALD collapse finds 4 and 6; the summed log gf agrees exactly, so only the metadata is wrong). RYA-1132 read the very census file that records this, stamped both rows 'COMPONENT_SUM_VERIFIED', and PROMOTED both to GF-LAB.

### A3-inflation - No HFS/isotope normalisation inflation, no max-component substitution: **PASS**

Every promoted Vujnovic total sits 0.009-0.036 dex from the census HFS SUM and 0.399-0.515 dex from the strongest component - it is a total-transition gf, as claimed. No log10(n_components) inflation is present.

### A3-11254 - 11254.9 strong-component vs unresolved-total caveat carried: **PASS**

The manifest adopts the blended-feature total (+0.3538) and keeps Burheim's strong component (+0.327) as evidence only; the conflict ledger names the distinction explicitly.

### A4 - No air<->vacuum and no cm-1<->Angstrom conflation in the crossmatch: **PASS**

For all three promoted doublets sharing an upper level, the vacuum wavenumber difference reproduces the Al I 3p ^2^P^o^ ground-term splitting to 0.007 cm-1 against the NIST value 112.061 cm-1. A medium conflation at 2650-3960 A would show as a 5-8 cm-1 offset and a unit conflation as orders of magnitude. This identity is not tabulated anywhere in the intake.

### A4-medium - Every wavelength carries an explicit medium: **FLAG**

There is no `medium` column in the manifest. Below 2000 A the census correctly stores a single vacuum wavelength, but it stores it in the column NAMED `wavelength_air`: for all 32 FUV rows `wavelength_air == wavelength_vac` exactly. The values are right and the label is wrong, so any downstream reader of `wavelength_air` silently receives vacuum wavelengths for those rows. (RYA-835/1001 units trap, RYA-938/944.)

### A5-lab - No theory graded as primary laboratory: **PASS**

All 18 GF-LAB rows carry EXP-BURHEIM23 (11), EXP-VUJNOVIC2002 (6) or EXP-JOHNSON1986 (1). None is an Opacity-Project or other theoretical source: the RYA-1001 `1995JPhB..` class (Mendoza+1995, OP theory) reaches the manifest only through `current_canonical_source`, never through a grade.

### A5-bf - No branching-fraction-only row is graded GF-LAB: **PASS**

Burheim Table 2 (branching fractions, no log gf) is not ingested; only Table 3's 12 derived log gf are. The six Vujnovic promotions each carry a finite 'this work' Aki with a stated uncertainty percentage, normalised by a LABORATORY lifetime (Buurman 1986 / Buurman & Donszelmann 1990 / Davidson 1990, CDS Table 1) - a lab composite, not a bare branching fraction.

### A5-overwrite - No lab/evaluated value silently overwritten by Kurucz/VALD: **PASS**

No GF-LAB row carries a Kurucz or VALD source. The one deliberate exception, 11254.9, is the documented blend-total substitution and is named in the conflict ledger rather than being silent.

### A5-doi-control - The DOI referee accepts correct identifiers: **PASS**

13 of 16 DOIs are confirmed by the same test - it accepts Burheim, Vujnovic, Trabert, Johnson, Kelleher, Papoulia, Roederer, Lind, Jonsson and Chiappino, matching through accented surnames (Vujnovic, Trabert, Jonsson) as well as plain ones. A negative is therefore a finding, not the detector's default.

### A5-doi - Every DOI resolves and names the paper it claims: **FAIL**

3 of 16 DOIs resolve to an UNRELATED paper. The referee is the registered AUTHOR LIST, because a volume comparison is not enough: `10.1086/312738` really is in ApJ 536, the volume the intake claims - it is just a different paper in it. In every case the intake's own prose citation is correct and the DOI beside it is not, so no artifact contradicts itself and RYA-1132's suite (which asserts only that `article_url` starts with 'https://') cannot see it. Resolution: committed cache; re-run with --online.

### A5-bibcode - Every ADS bibcode is a resolvable 19-character bibcode: **FLAG**

2 entries carry something that is not a bibcode: NIST_ASD = 'NIST_ASD'; Chiappino2026 = '2026ApJ...'. A truncated bibcode resolves to nothing and cannot referee its DOI - which matters here, because the bibcode is exactly what would have caught the Griesmann and Nandakumar page/volume mismatches.

### A5-johnson - Johnson 1986 Al II 2669 states the value the intake claims: **PASS**

The paper's abstract reads verbatim: 'The A-value for the intersystem transition is (3.33 +/- 0.23) x 10^3 s^-1 at the 90% confidence level', and 'Because there is only a single decay channel, the transition probability is the inverse of the radiative lifetime'. The intake's g_upper = 3 (3s3p ^3^P^o^_1_) and its log gf = -4.971830 reproduce exactly.

### A5-sigma - sigma(log gf) conversions justified with recorded provenance: **FLAG**

Three different conventions share one `gf_sigma_dex` column with nothing recording which: Burheim's published per-line dex uncertainty; Vujnovic's log10(1 + u) - the ASYMMETRIC upper bound, ~6% below the linear-propagation 1-sigma 0.434*u; and Johnson's 90%-CONFIDENCE bound stored as if it were 1-sigma (conservative by ~1.645x, and deliberately so, but only a code comment says so). The artifact needs a `sigma_basis` column, the lesson RYA-1084 and the sigma_stat/stat_basis finding already paid for.

### A6 - Competing gf values retained, never silently dropped: **FAIL**

6 manifest lines have a Vujnovic 2002 primary-laboratory log gf that THIS TICKET derived, crossmatched and then dropped. Not one reaches the conflict ledger, and `competing_gf_summary` names Vujnovic on 0 of 505 rows - it is hard-coded to 'Burheim=...; canonical=...; NIST=...' and has no slot for a source ingested later. The worst is 13123.416 A, one of the two best-graded Al lines RYA-1003 exists to unblock: Vujnovic +0.1901 +/- 0.0212 against the adopted Burheim +0.2320 +/- 0.0065, a 0.042 dex disagreement between two independent primary-laboratory measurements that the frozen manifest presents as a single unopposed GF-LAB value.

### B1 - Inventory reproduced independently from the sources: **PASS**

505 = 500 RYA-1001 census rows with ion in (I, II) + 1 IGRINS-only candidate + 4 Burheim mid-IR completeness controls. By ion 466/39; by source class 18/19/1/467. Every claimed count reproduces.

### B1-promotions - The 6 Al I promotions + Johnson Al II 2669, and the rejections: **PASS**

Promoted at [2652.475, 2660.386, 3082.153, 3092.71, 3944.006, 3961.52] A (the ticket quotes the Vujnovic source wavelengths 2652.484 / 2660.393; the manifest air wavelengths are 2652.475 / 2660.386, 0.009 and 0.007 A away - the same lines). Each has a finite 'this work' Aki with a stated uncertainty. Correctly NOT promoted: every `<` Aki limit (18 rows), every ratio-only row, and 3092.839, whose Aki carries the CDS `n_Aki` = ')' flag and no independent uncertainty.

### B1-yield - What the seven new GF-LAB promotions actually unblock: **FLAG**

All 7 lines promoted by RYA-1132 have Solar central depth 0.969-0.994, so under RYA-946 every one is DEEP Grade, not Codex Grade, and none enters the 0.05-0.60 measurement window. The intake records this honestly in `measurement_suitability_status`, but the manifest has NO column naming the RYA-946 grade, and it labels these seven 'CANDIDATE_NOT_SELECTED' - a rejection word for lines the contract says should be ROUTED TO SYNTHESIS.

### B2 - canonical_gf.csv not mutated by RYA-1132: **PASS**

PR #478 (merge 04e6afe) touches 25 files and none of them is under `data/linelists/`. `canonical_gf.csv` is byte-identical across the merge. The one data file it does change outside its own audit directory is `data/audit/rya1129_atomic_intake/intake_status_ledger.csv`, one row, as the builder documents.

### B3 - Band verdict strings are the ones claimed: **PASS**

UV/VIS/IR verdicts and `measurement_unblocked: false` match the ticket and `summary.json` exactly.

### B3-eligibility - 'Eligible' is derived from the evidence the manifest carries: **FLAG**

`measurement_suitability_status` is computed from Solar depth and the presence of a canonical id alone. It ignores the census `tier` column the manifest itself was built from, so 8 rows the RYA-1001 census adjudicated CANDIDATE-BLENDED are re-labelled ELIGIBLE_WITH_STATED_GF_TIER - including the GF-LAB line 11253.189 and the RYA-835 lines 7835.309 and 8772.865. The upgrade is silent; no column records that the census disagreed.

### B3-vocab - `gf_grade` expresses a gf grade: **FLAG**

The column mixes three vocabularies: gf provenance tiers (GF-LAB, VALD3, UNRESOLVED), NIST accuracy grades (B, B+, C, C+, D, E) and RYA-1001 SELECTION tiers (GRADEABLE, CANDIDATE-BLENDED, EXCLUDED-SHALLOW/SATURATED/NO-HOME). For 466 of 505 rows the value is a selection state, not a gf grade at all, and the one THEORETICAL row is labelled 'GRADEABLE' - a word that reads as an endorsement. None of RYA-946's four terms (Codex / Deep / Asplund / Consistent) appears anywhere in the manifest.

### C - Registered holdings reach the coverage module: **FAIL**

8 of 13 registered Solar holdings resolve to nothing through `pipeline.coverage.load_registry`, and ALL 5 crires_plus registrations are among them. Each `continue` is individually documented (RYA-776/929/931/945); the aggregate is that the one instrument reaching Al's IR lines is invisible to the module the census reads.

### C-lines - No reachable Al line is reported unreachable: **FAIL**

18 Al manifest lines (5 of them ELIGIBLE) sit inside the wavelength range of a Solar spectrum that is registered in the holdings registry AND present on disk, yet every one carries a BLANK `instrument_reach` in the frozen manifest. The IR verdict `BLOCKED_PIPELINE_COVERAGE` therefore rests partly on registry plumbing rather than on an absence of data.

### C-bands - RYA-1132's band() covers every band the census does: **FAIL**

`band()` leaves three uncovered intervals - 13000-13195.23 A, 17493.69-19510.4 A and >=24857.7 A - and every wavelength in them falls through to 'OUTSIDE_CURRENT_INSTRUMENT_REACH'. That relabels 13 lines the RYA-1001 census calls NIR, including the 2 GF-LAB lines 13123.416 and 13150.753 - the two best-graded Al lines RYA-1003 exists to unblock - and the Nandakumar/Chiappino member 17699.094. All of them carry `instruments_coverage_blind_spot = crires_plus` in the census, i.e. the census says an instrument reaches them. The rows are self-contradictory: band 'OUTSIDE_CURRENT_INSTRUMENT_REACH' beside `measurement_suitability_status = ELIGIBLE_WITH_STATED_GF_TIER`.

### C-alIII - Nothing is dropped from the census without being recorded: **FLAG**

2 Al III census rows are excluded by the builder ('Al III is outside this ticket's atomic scope'), which is a defensible scope call, but the exclusion appears only in a source comment - no artifact records that the 505 is a filtered denominator.

### NO-MUTATION - No intake artifact was modified by this QA: **PASS**

22 audited files hashed before and after every read; 0 changed. This auditor writes only under `data/audit/rya1141_al_intake_qa` and excludes its own source file by name.

## Findings

| severity | check | subject | finding |
| --- | --- | --- | --- |
| CRITICAL | A2 | `scripts/build_al_intake_rya1132.py:ingest_new_lab_sources` | Promotion join is wavelength-only, forbidden by RYA-1034 |
| CRITICAL | A2 | `alphys_II_3587.0720_0333` | One manifest line claimed by two different physical transitions |
| CRITICAL | A3 | `data/audit/rya1132_al_intake/al_line_manifest.csv:HFS_status` | 'COMPONENT_SUM_VERIFIED' is asserted from a count, never from a sum |
| CRITICAL | A3 | `canonical_gf Al I 3944.006` | hfs_n_components is wrong and the intake stamped it verified and promoted it |
| CRITICAL | A3 | `canonical_gf Al I 3961.520` | hfs_n_components is wrong and the intake stamped it verified and promoted it |
| CRITICAL | A5 | `source_bibliography.csv:GriesmannKling2000` | DOI 10.1086/312738 resolves to an unrelated paper |
| CRITICAL | A5 | `source_bibliography.csv:Nandakumar2024` | DOI 10.3847/1538-4357/ad4451 resolves to an unrelated paper |
| CRITICAL | A5 | `web_source_followup.csv:Murphy & Berengut 2014 / Griesmann & Kling 2000` | DOI 10.1093/mnras/stt2120 resolves to an unrelated paper |
| CRITICAL | A6 | `alphys_II_3900.6750_0334 (3900.675 A)` | Competing Vujnovic primary-lab gf derived in this ticket and then dropped |
| CRITICAL | A6 | `alphys_I_13123.4160_0423 (13123.416 A)` | Competing Vujnovic primary-lab gf derived in this ticket and then dropped |
| CRITICAL | C | `alphys_I_13123.4160_0423 (13123.416 A)` | Census NIR line relabelled OUTSIDE_CURRENT_INSTRUMENT_REACH by a band() gap |
| CRITICAL | C | `alphys_I_13150.7530_0425 (13150.753 A)` | Census NIR line relabelled OUTSIDE_CURRENT_INSTRUMENT_REACH by a band() gap |
| HIGH | A1 | `table5.dat row 1` | CDS flag column `n_Lambda` = '*' is dropped by the parser |
| HIGH | A6 | `alphys_I_13150.7530_0425 (13150.753 A)` | Competing Vujnovic primary-lab gf derived in this ticket and then dropped |
| HIGH | A6 | `alphys_I_21093.0290_0480 (21093.029 A)` | Competing Vujnovic primary-lab gf derived in this ticket and then dropped |
| HIGH | A6 | `alphys_I_21163.7550_0481 (21163.755 A)` | Competing Vujnovic primary-lab gf derived in this ticket and then dropped |
| HIGH | A6 | `alphys_I_3092.8390_0331 (3092.839 A)` | Competing Vujnovic primary-lab gf derived in this ticket and then dropped |
| HIGH | C | `holdings:kpno_solar_atlas` | Registered Solar holding silently dropped (SKIPPED_BY_SUFFIX_JSON) |
| HIGH | C | `holdings:crires_plus` | Registered Solar holding silently dropped (SKIPPED_BY_SUFFIX_TXT) |
| HIGH | C | `holdings:crires_plus` | Registered Solar holding silently dropped (SKIPPED_NO_LOADER_COLUMN) |
| HIGH | C | `holdings:crires_plus` | Registered Solar holding silently dropped (SKIPPED_NO_LOADER_COLUMN) |
| HIGH | C | `holdings:harps` | Registered Solar holding silently dropped (SKIPPED_BY_SUFFIX_JSON) |
| HIGH | C | `holdings:kpno_solar_atlas` | Registered Solar holding silently dropped (SKIPPED_BY_SUFFIX_MD) |
| HIGH | C | `holdings:crires_plus` | Registered Solar holding silently dropped (SKIPPED_NO_LOADER_COLUMN) |
| HIGH | C | `holdings:crires_plus` | Registered Solar holding silently dropped (SKIPPED_NO_LOADER_COLUMN) |
| HIGH | C | `alphys_I_13122.3220_0422 (13122.322 A)` | Census NIR line relabelled OUTSIDE_CURRENT_INSTRUMENT_REACH by a band() gap |
| HIGH | C | `alphys_I_13130.1740_0424 (13130.174 A)` | Census NIR line relabelled OUTSIDE_CURRENT_INSTRUMENT_REACH by a band() gap |
| HIGH | C | `alphys_I_17699.0940_0462 (17699.094 A)` | Census NIR line relabelled OUTSIDE_CURRENT_INSTRUMENT_REACH by a band() gap |
| HIGH | C | `alphys_I_17708.0730_0463 (17708.073 A)` | Census NIR line relabelled OUTSIDE_CURRENT_INSTRUMENT_REACH by a band() gap |
| HIGH | C | `alphys_I_18942.3460_0464 (18942.346 A)` | Census NIR line relabelled OUTSIDE_CURRENT_INSTRUMENT_REACH by a band() gap |
| HIGH | C | `alphys_I_18956.9260_0465 (18956.926 A)` | Census NIR line relabelled OUTSIDE_CURRENT_INSTRUMENT_REACH by a band() gap |
| HIGH | C | `alphys_I_19280.1420_0466 (19280.142 A)` | Census NIR line relabelled OUTSIDE_CURRENT_INSTRUMENT_REACH by a band() gap |
| HIGH | C | `alphys_I_19302.1310_0467 (19302.131 A)` | Census NIR line relabelled OUTSIDE_CURRENT_INSTRUMENT_REACH by a band() gap |
| HIGH | C | `alphys_I_19482.2860_0468 (19482.286 A)` | Census NIR line relabelled OUTSIDE_CURRENT_INSTRUMENT_REACH by a band() gap |
| HIGH | C | `alphys_I_19497.7080_0469 (19497.708 A)` | Census NIR line relabelled OUTSIDE_CURRENT_INSTRUMENT_REACH by a band() gap |
| HIGH | C | `alphys_I_24985.9530_0499 (24985.953 A)` | Census NIR line relabelled OUTSIDE_CURRENT_INSTRUMENT_REACH by a band() gap |
| MEDIUM | A1 | `table2.dat row 5` | CDS flag column `l_e_Aki` = '>' is dropped by the parser |
| MEDIUM | A1 | `table2.dat row 6` | CDS flag column `l_e_Aki` = '>' is dropped by the parser |
| MEDIUM | A1 | `table2.dat row 10` | CDS flag column `n_Aki` = ')' is dropped by the parser |
| MEDIUM | A1 | `table5.dat row 10` | CDS flag column `l_e_Aki` = '>' is dropped by the parser |
| MEDIUM | A1 | `table5.dat row 13` | CDS flag column `l_e_Aki` = '>' is dropped by the parser |
| MEDIUM | A4 | `al_line_manifest.csv:wavelength_air` | Column named `wavelength_air` holds VACUUM values below 2000 A |
| MEDIUM | A5 | `source_bibliography.csv:NIST_ASD` | ads_bibcode is a stub, not a resolvable 19-character bibcode |
| MEDIUM | A5 | `source_bibliography.csv:Chiappino2026` | ads_bibcode is a stub, not a resolvable 19-character bibcode |
| MEDIUM | A5 | `al_line_manifest.csv:gf_sigma_dex` | One column carries three different uncertainty conventions, unlabelled |
| MEDIUM | B3 | `alphys_I_6696.7880_0354 (6696.788 A)` | Manifest calls ELIGIBLE a line the RYA-1001 census tiered CANDIDATE-BLENDED |
| MEDIUM | B3 | `alphys_I_6784.2560_0357 (6784.256 A)` | Manifest calls ELIGIBLE a line the RYA-1001 census tiered CANDIDATE-BLENDED |
| MEDIUM | B3 | `alphys_I_6905.6460_0358 (6905.646 A)` | Manifest calls ELIGIBLE a line the RYA-1001 census tiered CANDIDATE-BLENDED |
| MEDIUM | B3 | `alphys_I_6906.2870_0359 (6906.287 A)` | Manifest calls ELIGIBLE a line the RYA-1001 census tiered CANDIDATE-BLENDED |
| MEDIUM | B3 | `alphys_I_7835.3090_0386 (7835.309 A)` | Manifest calls ELIGIBLE a line the RYA-1001 census tiered CANDIDATE-BLENDED |
| MEDIUM | B3 | `alphys_I_8772.8650_0393 (8772.865 A)` | Manifest calls ELIGIBLE a line the RYA-1001 census tiered CANDIDATE-BLENDED |
| MEDIUM | B3 | `alphys_I_11253.1890_0406 (11253.189 A)` | Manifest calls ELIGIBLE a line the RYA-1001 census tiered CANDIDATE-BLENDED |
| MEDIUM | B3 | `alphys_I_15956.6750_0449 (15956.675 A)` | Manifest calls ELIGIBLE a line the RYA-1001 census tiered CANDIDATE-BLENDED |

## Independently reproduced inventory

- rows: **505** (claimed 505)
- by ion: **466 Al I / 39 Al II** (claimed 466 / 39)
- by source class: **18 / 19 / 1 / 467** (claimed 18 / 19 / 1 / 467)

