# Post-release Solar gap reconciliation

Base science commit `c058aa97fa828cfe6ae53a78b839f4bbce635213`. The V I v1.24 and Si I v1.422 releases are live. N I 7468 Å and V II 3760 Å completed follow-up fits and remain diagnostic. Their frozen fit and uncertainty evidence is in the [follow-up commit](https://github.com/damienabraxas/exoplanetcodex/commit/14f707bf).

## Reuse evidence

[candidate_stage_evidence.csv](candidate_stage_evidence.csv) links 1,257 audit candidate IDs to committed nominal stage line fits. Of those IDs, 1,100 have an accepted fit in at least one stage and 157 have only excluded stage evidence. The match uses the canonical physical ID, species, air wavelength within 0.005 Å and excitation within 0.02 eV. It requires a unique canonical identity and an explicit accepted or excluded flag. File counts refer to evidence files, which can include separate estimates for the same spectrum. A previous fit supplies reuse evidence. A different oscillator strength, coverage question or stage rejection can still require adjudication. The [reconciler](reconcile_candidate_stage_evidence.py) runs from the repository root.

## Calcium

The [Aldenius VIS triplet](ca_i_aldenius_triplet.csv) is now joined by both levels and J using the committed NIST Ca I levels and the project air conversion. All three share the 31539.495 cm⁻¹ upper level, so the lifetime term must be common in a future joint uncertainty calculation. The laboratory strengths are 0.057–0.080 dex lower than the current canonical strengths. The canonical rows still carry grade A while the repository's RYA-853 audit reports grade C for the same triplet. An update needs an explicit source and uncertainty decision.

The [new HARPS profile measurements](execution/cai_vis_harps_ew_v1/) returned EWs of 206.09, 228.99 and 266.23 mÅ in ascending wavelength order. All three are above the −4.90 linear-regime ceiling, with REW −4.4715, −4.4271 and −4.3645. The older 120.73 mÅ measurement for the 6122 Å line also falls above the ceiling. No line enters the linear-EW abundance pool. A separate strong-wing analysis would need a validated blend and damping model and a resolved oscillator-strength and covariance decision.

## Oxygen

The 615.8 nm hold is now narrowed using the committed RYA-1353 source snapshot and RYA-1136 Amarsi grid transcription. The three NIST evaluated strengths at 6158.149, 6158.172 and 6158.180 Å are −1.8410, −0.9965 and −0.4090 dex. The Amarsi rows reproduce that sequence to 0.002 dex. Pereira's listed −1.841 is the weak-component strength in both independent tables. It cannot replace the −0.409 strength on the stronger canonical row at 6158.187 Å. The first two components have no separate canonical identity in this intake. The [three-component identity table](o_i_6158_components.csv) also records both energies and J from [NIST ASD](https://physics.nist.gov/cgi-bin/ASD/lines1.pl?A_out=0&J_out=on&allowed_out=1&bibrefs=1&conf_out=on&en_unit=0&enrg_out=on&line_out=0&low_w=615.81&show_calc_wl=1&show_obs_wl=1&spectra=O+I&term_out=on&unit=1&upp_w=615.82). Summing the three linear gf values gives log(gf) = −0.296605, nearly equal to the GES deck's single-row −0.296. This is consistent with a collapsed multiplet. The canonical −0.409 belongs to the strongest component, so applying it to that single deck row loses about 23% of the component-summed strength. The feature remains held for component-complete synthesis, validation of molecular blends on the observed holding, and an applicable departure treatment.

## Run census

The mounted and scratch Sirius inventory (to depth seven) and science main resolve a Solar Ni run in `~/scratch/rya731` in addition to the previously identified committed products. Its latest ORCH checkpoint reports 150 unpublished cells and 240 held cells. It contains 420 Ni product CSV artifacts across routes and vintages. Those are artifact counts. Ni first-publication workflow remains under RYA-731. Its Ni II retry log ends with a preflight failure, and the release staging log was empty when inspected.

The same scan returned no named Solar feed or ORCH latest report for P, S, Co, K, Ba, Y, Cu, Zn, Mn, Sc, Li, Eu, Zr or Sr. This bounded search does not rule out other artifact names or holdings.

## Iron red-optical execution

The four unmatched reference-set candidates are Fe I 8293.513 Å and Fe II 7222.391, 7224.478 and 7515.830 Å. The pinned initial 1D LTE path selected EW inversion and stopped before synthesis for lack of scoped measured EWs. The isolated retry requests synthesis explicitly on the molecfit-corrected Kitt Peak holding and writes under `~/scratch/rya1397-feii-v1/rya1397-fe-redopt-v2`. The run uses the existing canonical gf scale. The 7222 and 7224 Å windows sit inside the failed 7160–7340 Å telluric interval on that holding. An alternate Kurucz 2005 diagnostic covers the same Fe II candidates on a separate registered holding. The literature strengths differ by −0.039 dex for Fe I and +0.100, +0.040 and +0.070 dex for Fe II. Treat the fits as diagnostics pending scale review and uncertainty assembly.

## New laboratory candidates against the installed decks

[Deck coverage](new_lab_deck_coverage.csv) checks 352 lab candidate records against the pinned Sirius deck hashes, matching wavelength, both energies and published J when available. It found 298 outside the scoped decks, 30 without a concordant tuple in the corresponding deck, and 24 with a concordant tuple. These counts are source records and include duplicate literature claims. The deck match clears only the tuple gate.

The Blackwell-Whitehead Ti I line listed at 23176.915 Å is the existing 23176.871 Å feature (`pk_dff5dfdb4f00`). Six committed treatment fits excluded it for nonphysical abundance, an edge-pinned fit or failure to find an interior minimum. The literature strength differs from the adopted NIST strength by 0.053 dex. That difference is too small to account for the LTE fit being 2.82 dex below the solar anchor. Preserve the prior rejection.

Four Wood Cr II near-UV tuples at 3603.828, 3613.124, 3677.898 and 3712.887 Å have the reported strengths in the installed deck and no corresponding committed Cr II stage fit within 0.05 Å and 0.02 eV. They were dispatched as an isolated, non-publishing diagnostic on the corrected Kitt Peak holding. The same preflight found no prior stage rows for two Sobeck Cr I candidates at 3749.027 and 4829.309 Å with matching deck strengths. Their diagnostics were also dispatched. The 4530 and 4535 Å Cr I features have competing measured strengths and remain under source review. Eight other deck-matched Cr I transitions have accepted stage fits.

## Nitrogen follow-up

The 12 Caffau Table 1 source rows crossmatch to six unique red-optical canonical rows by species, air wavelength and excitation, including the existing 8216 and 8683 Å fits, the rejected 7442 Å fit, and the completed 7468 Å diagnostic. The 8223.128 and 8718.837 Å rows use NIST B+ strengths. Neither was in the committed Solar feed, checked stage files or active problem-line registry, so both were dispatched as an isolated synthesis diagnostic on the corrected Kitt Peak holding.

The other six rows match the NIST evaluated snapshot by wavelength and excitation, and were absent from the canonical table and installed deck. The held full-band VALD extract supplies four of them, with excitation, upper energy and both J values. The [Y-band diagnostic runner](run_ni_nir_deck_diagnostic.py) adds those four to a process-local deck and retains its input and output hashes. It uses the project converter's damping defaults where the extract supplies no coefficient and leaves the shared deck unchanged. Fits use the corrected CRIRES+ Y holding, with local telluric residual and width validation still required for any release. The 10507.000 and 10520.580 Å candidates remain held for missing full deck input. The nearby 10108.892 Å deck row is a separate existing transition.


## Additional source and legacy checks

The [legacy Fe I VIS comparison](fe_legacy_vis_sidecar_comparison.csv) verifies four of fourteen summary-only published products against their original source product hash, accepted-line counts and median abundance. Each is a 67-line Kurucz 2005 treatment. Ten still lack an exact per-line sidecar. The alternate RYA-1120 Kitt Peak profile-fit bundle on Sirius has a different product-source hash, so it cannot establish which lines were in the published result.

The held O'Brian 1991 article reports 1,174 transition probabilities from lifetime and branching-fraction measurements and another 640 using interpolated level populations. It prints sample transition results and says the full data were available from the authors on request. Those two evidence classes must stay distinct. The O'Brian near-UV source-complete comparison remains open.

The [Al level join](al_papoulia_level_join.csv) parses the held Papoulia supplement. It binds both published level labels for 587 of 588 transition rows. Observed level energies recover a provisional canonical tuple for 39 Al I and 38 Al II cases. Another 225 fall inside the scoped wavelengths without a unique canonical tuple, 268 lie outside the scoped wavelengths, 17 have no usable observed-level gap and one has an unresolved source level label. The match uses the observed-energy Ritz wavelength and lower excitation. It supplies candidate identity evidence for computed transition strengths and does not supersede laboratory oscillator strengths. The held Belyaev 2013 Al–H paper provides collision rates for projection into an atom, rather than oscillator strengths or ready-made abundance corrections. Existing Solar Al I Engine A uses the Nordlander and Lind model through the GALAH grid.

The held Osorio 2022 Ca II infrared triplet paper's Table 3 grid spans log g 1–4 and metallicities −6 to −2.5 dex. The Solar parameters lie outside those axes. Retain those published corrections as metal-poor evidence. A Solar wing and hydrogen-blend model needs its own validation.

The existing RYA-1136 molecular crossmatch has 358 of 408 rows with a primary or physical tuple match. Six more are matched unresolved CH sums, and 44 carry ambiguity, strength or energy discrepancies or remain unmatched. Those rows need component-level review before a pool change. Existing C and N molecular route products are present in the live feeds. The prior CRIRES+ H OH route produced A(O) = 9.448 and reduced chi-square 130.65 in the RYA-1220 profile audit. Its 14 line-list validation matches are evidence for transition identity, and the abundance fit remains held for known profile residuals.

## Carbon and silicon source levels

The [carbon and silicon level join](c_si_source_level_join.csv) checks 2,798 candidate source rows against the observed-energy level tables. Both labels and the energy gap resolve for 637 C I, 1,169 Si I and 463 Si II source rows. Within the scoped wavelengths, 349 C I, 599 Si I and 45 Si II entries have unique provisional canonical tuples. They still require strength, observability and grade review. The label join deliberately leaves 148 C I and 381 Si I records unresolved rather than guessing on alternative configurations.

The 17 experimental Pehlivan Rhodin Table B3 results have a concrete disposition. For 14, the canonical strength already agrees exactly. Two of those transitions are accepted and published, and 11 already have active saturation exclusions. The 14224.526 Å line is in the J-to-H holding gap. The other three experimental rows lie outside the scoped wavelengths. No repeat fit is warranted for this lab table on the current holdings.


## Den Hartog Fe II literature dispositions

Table 6 contains 131 transitions. The 109 outside-band rows remain out of scope for the existing holdings and nine published rows lie in the paper's rejected 3170–3227 Å interval. Section 4.1 says 11 lines were removed in that interval. The held Table 6 contains nine there and has two more near-UV lines immediately redward at 3255.888 and 3277.349 Å. The text does not establish that the latter two are the missing members of the interval count. Preserve the accepted existing project products for those two.

Of the 13 remaining reference-table identities, 11 are already present in the Fe II feed. The published 4303.170 Å line is independently excluded by the project's active bad-fit registry and by the paper's CH/atomic blend review. The paper also rejects 4351.762 Å, which has accepted fits in project pools. Record the conflict without changing those products. The 4173.451 Å Table 7 feature remains shortward of the installed GES optical synthesis deck.
