# Post-release Solar gap reconciliation

Base science commit `c058aa97fa828cfe6ae53a78b839f4bbce635213`. The V I v1.24 and Si I v1.422 releases are live. N I 7468 Å and V II 3760 Å completed follow-up fits and remain diagnostic. Their frozen fit and uncertainty evidence is in the [follow-up commit](https://github.com/damienabraxas/exoplanetcodex/commit/14f707bf).

## Reuse evidence

[candidate_stage_evidence.csv](candidate_stage_evidence.csv) links 1,257 audit candidate IDs to committed nominal stage line fits. Of those IDs, 1,100 have an accepted fit in at least one stage and 157 have only excluded stage evidence. The match uses the canonical physical ID, species, air wavelength within 0.005 Å and excitation within 0.02 eV. It requires a unique canonical identity and an explicit accepted or excluded flag. File counts refer to evidence files, which can include separate estimates for the same spectrum. A previous fit supplies reuse evidence. A different oscillator strength, coverage question or stage rejection can still require adjudication. The [reconciler](reconcile_candidate_stage_evidence.py) runs from the repository root.

## Calcium

The [Aldenius VIS triplet](ca_i_aldenius_triplet.csv) is now joined by both levels and J using the committed NIST Ca I levels and the project air conversion. All three share the 31539.495 cm⁻¹ upper level, so the lifetime term must be common in a future joint uncertainty calculation. The laboratory strengths are 0.057–0.080 dex lower than the current canonical strengths. The canonical rows still carry grade A while the repository's RYA-853 audit reports grade C for the same triplet. An update needs an explicit source and uncertainty decision.

The committed legacy Ca I 6122 Å equivalent width implies REW −4.705, above the −4.90 linear-regime ceiling. The remaining two have no equivalent width in that same artifact. None of the three is promoted by this identity review. Ca II IR model coverage remains a separate gate.

## Oxygen

The 615.8 nm hold is now narrowed using the committed RYA-1353 source snapshot and RYA-1136 Amarsi grid transcription. The three NIST evaluated strengths at 6158.149, 6158.172 and 6158.180 Å are −1.8410, −0.9965 and −0.4090 dex. The Amarsi rows reproduce that sequence to 0.002 dex. Pereira's listed −1.841 is the weak-component strength in both independent tables. It cannot replace the −0.409 strength on the stronger canonical row at 6158.187 Å. The first two components have no separate canonical identity in this intake. The feature remains held for component-complete synthesis, validation of molecular blends on the observed holding, and an applicable departure treatment.

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

The other six rows match the NIST evaluated snapshot by wavelength and excitation, but have no canonical or synthesis-deck tuple. The nearby 10108.892 Å deck entry is a different transition. All six NIR candidates require a complete atomic deck input with damping and an observed-window telluric check before a fit can count for release. The Y-band holding's nominal wavelength coverage does not settle either gate.
