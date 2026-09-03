# RYA-1187 — (holding x band x engine x grade) applicability matrix

Built from LINE-LEVEL data against `Fe.json` v1.93. Every cell is a live product or a stated reason; a silent empty is the defect this exists to remove.

## The headline question: Deep Grade in red-optical and NIR

**Answered NO LINES, not owed product.** Per-band primary-lab-gf Fe census (the product path's own feature depth, gate 0.60):

| band | species | lab lines | shallow (Codex) | deep (Deep) | deepest feature |
|---|---|---|---|---|---|
| near-UV | Fe I | 59 | 0 | **58** | 0.994 |
| near-UV | Fe II | 12 | 0 | **12** | 0.992 |
| VIS | Fe I | 243 | 66 | **174** | 0.928 |
| VIS | Fe II | 10 | 0 | **10** | 0.869 |
| red-optical | Fe I | 70 | 66 | **0** | 0.583 |
| red-optical | Fe II | 0 | 0 | **0** | None |
| NIR | Fe I | 29 | 29 | **0** | 0.409 |
| NIR | Fe II | 0 | 0 | **0** | None |

Red-optical misses the 0.60 gate by **0.017** (deepest lab feature 0.583) and NIR by **0.191** (0.409). Fe II has no primary-lab gf in either band at all. Emitting a Deep product there would need a gf we do not hold or a grade reassignment — RYA-161 forbids both.

## Cell matrix — holding x band x grade

| holding | band | grade | verdict | n lines | reason |
|---|---|---|---|---|---|
| solar_kpno_kurucz2005_corrected | near-UV | Codex Grade | **ABSENT_NO_LAB_LINE_AT_THIS_DEPTH** | 0 | no primary-lab-gf Fe line with feature depth <=0.6 in 3000-3780 A (deepest lab feature in band = 0.994) |
| solar_kpno_kurucz2005_corrected | near-UV | Deep Grade | **LIVE** | 70 | 70 primary-lab-gf Fe line(s) at this depth grade (>0.6) in 3000-3780 A |
| solar_kpno_kurucz2005_corrected | near-UV | Reference Grade | **ABSENT_NO_REFERENCE_LINE** | 0 | AGSS21 Table A.2 publishes no Fe line inside near-UV ∩ this holding's span (3000-3780 A). The reference set is optical; this is a published SELECTION, not a coverage hole. |
| solar_kpno_kurucz2005_corrected | VIS | Codex Grade | **LIVE** | 66 | 66 primary-lab-gf Fe line(s) at this depth grade (<=0.6) in 3780-6910 A |
| solar_kpno_kurucz2005_corrected | VIS | Deep Grade | **LIVE** | 184 | 184 primary-lab-gf Fe line(s) at this depth grade (>0.6) in 3780-6910 A |
| solar_kpno_kurucz2005_corrected | VIS | Reference Grade | **LIVE** | 36 | 36 AGSS21 Table A.2 line(s) in VIS within this holding's coverage |
| solar_kpno_kurucz2005_corrected | red-optical | Codex Grade | **LIVE** | 66 | 66 primary-lab-gf Fe line(s) at this depth grade (<=0.6) in 6910-9199 A |
| solar_kpno_kurucz2005_corrected | red-optical | Deep Grade | **ABSENT_NO_LAB_LINE_AT_THIS_DEPTH** | 0 | no primary-lab-gf Fe line with feature depth >0.6 in 6910-9199 A (deepest lab feature in band = 0.583) |
| solar_kpno_kurucz2005_corrected | red-optical | Reference Grade | **ABSENT_ENGINE_DOMAIN** | 16 | 16 AGSS21 line(s) fall in red-optical, but ALL of the 12 AGSS21 Fe I line(s) here are OUTSIDE the Amarsi 2022 MLP's training box — every one fails `transition energy Eup-Elo outside training [1.8190, 2.5898] eV`. dE IS the wavelength, so a redder line is below the box by construction. The ENGINE-A-3DNLTE route that produced the four live Reference products cannot reach here without extrapolating the network outside its own domain, which classify_line refuses. MEASURED: amarsi_domain_check.csv. |
| solar_kpno_kurucz2005_corrected | NIR | Codex Grade | **BLOCKED_BAND_TABLE_CONFLICT** | 24 | 24 line(s) at this grade, but this holding's NIR coverage (9199-10000 A) lies entirely inside the slice the two band tables dispute: config/synth_bands.yaml calls 9199-13000 NIR while pipeline.band_policy calls 6910-10000 RED-OPTICAL. A run over this window resolves to red-optical and takes the PROFILEFIT route (measured: it refused with 'no measured EWs ... PROFILEFIT'), so the product cannot be emitted as NIR until RYA-1094's two-table conflict is settled. Not this ticket's to decide. |
| solar_kpno_kurucz2005_corrected | NIR | Deep Grade | **ABSENT_NO_LAB_LINE_AT_THIS_DEPTH** | 0 | no primary-lab-gf Fe line with feature depth >0.6 in 9199-10000 A (deepest lab feature in band = 0.344) |
| solar_kpno_kurucz2005_corrected | NIR | Reference Grade | **ABSENT_ENGINE_DOMAIN** | 1 | 1 AGSS21 line(s) fall in NIR, but ALL of the 1 AGSS21 Fe I line(s) here are OUTSIDE the Amarsi 2022 MLP's training box — every one fails `transition energy Eup-Elo outside training [1.8190, 2.5898] eV`. dE IS the wavelength, so a redder line is below the box by construction. The ENGINE-A-3DNLTE route that produced the four live Reference products cannot reach here without extrapolating the network outside its own domain, which classify_line refuses. MEASURED: amarsi_domain_check.csv. |
| solar_kpno_molecfit_corrected | near-UV | Codex Grade | **ABSENT_NO_LAB_LINE_AT_THIS_DEPTH** | 0 | no primary-lab-gf Fe line with feature depth <=0.6 in 3000-3780 A (deepest lab feature in band = 0.994) |
| solar_kpno_molecfit_corrected | near-UV | Deep Grade | **LIVE** | 70 | 70 primary-lab-gf Fe line(s) at this depth grade (>0.6) in 3000-3780 A |
| solar_kpno_molecfit_corrected | near-UV | Reference Grade | **ABSENT_NO_REFERENCE_LINE** | 0 | AGSS21 Table A.2 publishes no Fe line inside near-UV ∩ this holding's span (3000-3780 A). The reference set is optical; this is a published SELECTION, not a coverage hole. |
| solar_kpno_molecfit_corrected | VIS | Codex Grade | **LIVE** | 66 | 66 primary-lab-gf Fe line(s) at this depth grade (<=0.6) in 3780-6910 A |
| solar_kpno_molecfit_corrected | VIS | Deep Grade | **LIVE** | 184 | 184 primary-lab-gf Fe line(s) at this depth grade (>0.6) in 3780-6910 A |
| solar_kpno_molecfit_corrected | VIS | Reference Grade | **LIVE** | 36 | 36 AGSS21 Table A.2 line(s) in VIS within this holding's coverage |
| solar_kpno_molecfit_corrected | red-optical | Codex Grade | **LIVE** | 66 | 66 primary-lab-gf Fe line(s) at this depth grade (<=0.6) in 6910-9199 A |
| solar_kpno_molecfit_corrected | red-optical | Deep Grade | **ABSENT_NO_LAB_LINE_AT_THIS_DEPTH** | 0 | no primary-lab-gf Fe line with feature depth >0.6 in 6910-9199 A (deepest lab feature in band = 0.583) |
| solar_kpno_molecfit_corrected | red-optical | Reference Grade | **ABSENT_ENGINE_DOMAIN** | 16 | 16 AGSS21 line(s) fall in red-optical, but ALL of the 12 AGSS21 Fe I line(s) here are OUTSIDE the Amarsi 2022 MLP's training box — every one fails `transition energy Eup-Elo outside training [1.8190, 2.5898] eV`. dE IS the wavelength, so a redder line is below the box by construction. The ENGINE-A-3DNLTE route that produced the four live Reference products cannot reach here without extrapolating the network outside its own domain, which classify_line refuses. MEASURED: amarsi_domain_check.csv. |
| solar_kpno_molecfit_corrected | NIR | Codex Grade | **LIVE** | 29 | 29 primary-lab-gf Fe line(s) at this depth grade (<=0.6) in 9199-13000 A |
| solar_kpno_molecfit_corrected | NIR | Deep Grade | **ABSENT_NO_LAB_LINE_AT_THIS_DEPTH** | 0 | no primary-lab-gf Fe line with feature depth >0.6 in 9199-13000 A (deepest lab feature in band = 0.409) |
| solar_kpno_molecfit_corrected | NIR | Reference Grade | **ABSENT_ENGINE_DOMAIN** | 1 | 1 AGSS21 line(s) fall in NIR, but ALL of the 1 AGSS21 Fe I line(s) here are OUTSIDE the Amarsi 2022 MLP's training box — every one fails `transition energy Eup-Elo outside training [1.8190, 2.5898] eV`. dE IS the wavelength, so a redder line is below the box by construction. The ENGINE-A-3DNLTE route that produced the four live Reference products cannot reach here without extrapolating the network outside its own domain, which classify_line refuses. MEASURED: amarsi_domain_check.csv. |
| solar_harps_molecfit_corrected | near-UV | Codex Grade | **ABSENT_NO_COVERAGE** | 0 | this holding covers 3780-6910 A (INSTRUMENT catalog — UPPER BOUND, this holding declares no span of its own) and the near-UV band is 3000-3780 A — they do not overlap. No spectrum exists to measure. |
| solar_harps_molecfit_corrected | near-UV | Deep Grade | **ABSENT_NO_COVERAGE** | 0 | this holding covers 3780-6910 A (INSTRUMENT catalog — UPPER BOUND, this holding declares no span of its own) and the near-UV band is 3000-3780 A — they do not overlap. No spectrum exists to measure. |
| solar_harps_molecfit_corrected | near-UV | Reference Grade | **ABSENT_NO_COVERAGE** | 0 | this holding covers 3780-6910 A (INSTRUMENT catalog — UPPER BOUND, this holding declares no span of its own) and the near-UV band is 3000-3780 A — they do not overlap. No spectrum exists to measure. |
| solar_harps_molecfit_corrected | VIS | Codex Grade | **LIVE** | 66 | 66 primary-lab-gf Fe line(s) at this depth grade (<=0.6) in 3780-6910 A |
| solar_harps_molecfit_corrected | VIS | Deep Grade | **LIVE** | 184 | 184 primary-lab-gf Fe line(s) at this depth grade (>0.6) in 3780-6910 A |
| solar_harps_molecfit_corrected | VIS | Reference Grade | **LIVE** | 36 | 36 AGSS21 Table A.2 line(s) in VIS within this holding's coverage |
| solar_harps_molecfit_corrected | red-optical | Codex Grade | **ABSENT_NO_COVERAGE** | 0 | this holding covers 3780-6910 A (INSTRUMENT catalog — UPPER BOUND, this holding declares no span of its own) and the red-optical band is 6910-9199 A — they do not overlap. No spectrum exists to measure. |
| solar_harps_molecfit_corrected | red-optical | Deep Grade | **ABSENT_NO_COVERAGE** | 0 | this holding covers 3780-6910 A (INSTRUMENT catalog — UPPER BOUND, this holding declares no span of its own) and the red-optical band is 6910-9199 A — they do not overlap. No spectrum exists to measure. |
| solar_harps_molecfit_corrected | red-optical | Reference Grade | **ABSENT_NO_COVERAGE** | 0 | this holding covers 3780-6910 A (INSTRUMENT catalog — UPPER BOUND, this holding declares no span of its own) and the red-optical band is 6910-9199 A — they do not overlap. No spectrum exists to measure. |
| solar_harps_molecfit_corrected | NIR | Codex Grade | **ABSENT_NO_COVERAGE** | 0 | this holding covers 3780-6910 A (INSTRUMENT catalog — UPPER BOUND, this holding declares no span of its own) and the NIR band is 9199-13000 A — they do not overlap. No spectrum exists to measure. |
| solar_harps_molecfit_corrected | NIR | Deep Grade | **ABSENT_NO_COVERAGE** | 0 | this holding covers 3780-6910 A (INSTRUMENT catalog — UPPER BOUND, this holding declares no span of its own) and the NIR band is 9199-13000 A — they do not overlap. No spectrum exists to measure. |
| solar_harps_molecfit_corrected | NIR | Reference Grade | **ABSENT_NO_COVERAGE** | 0 | this holding covers 3780-6910 A (INSTRUMENT catalog — UPPER BOUND, this holding declares no span of its own) and the NIR band is 9199-13000 A — they do not overlap. No spectrum exists to measure. |
| solar_iag | near-UV | Codex Grade | **ABSENT_NO_COVERAGE** | 0 | this holding covers 4047-10650 A (INSTRUMENT catalog — UPPER BOUND, this holding declares no span of its own) and the near-UV band is 3000-3780 A — they do not overlap. No spectrum exists to measure. |
| solar_iag | near-UV | Deep Grade | **ABSENT_NO_COVERAGE** | 0 | this holding covers 4047-10650 A (INSTRUMENT catalog — UPPER BOUND, this holding declares no span of its own) and the near-UV band is 3000-3780 A — they do not overlap. No spectrum exists to measure. |
| solar_iag | near-UV | Reference Grade | **ABSENT_NO_COVERAGE** | 0 | this holding covers 4047-10650 A (INSTRUMENT catalog — UPPER BOUND, this holding declares no span of its own) and the near-UV band is 3000-3780 A — they do not overlap. No spectrum exists to measure. |
| solar_iag | VIS | Codex Grade | **LIVE** | 65 | 65 primary-lab-gf Fe line(s) at this depth grade (<=0.6) in 4047-6910 A |
| solar_iag | VIS | Deep Grade | **LIVE** | 157 | 157 primary-lab-gf Fe line(s) at this depth grade (>0.6) in 4047-6910 A |
| solar_iag | VIS | Reference Grade | **LIVE** | 36 | 36 AGSS21 Table A.2 line(s) in VIS within this holding's coverage |
| solar_iag | red-optical | Codex Grade | **LIVE** | 66 | 66 primary-lab-gf Fe line(s) at this depth grade (<=0.6) in 6910-9199 A |
| solar_iag | red-optical | Deep Grade | **ABSENT_NO_LAB_LINE_AT_THIS_DEPTH** | 0 | no primary-lab-gf Fe line with feature depth >0.6 in 6910-9199 A (deepest lab feature in band = 0.583) |
| solar_iag | red-optical | Reference Grade | **ABSENT_ENGINE_DOMAIN** | 16 | 16 AGSS21 line(s) fall in red-optical, but ALL of the 12 AGSS21 Fe I line(s) here are OUTSIDE the Amarsi 2022 MLP's training box — every one fails `transition energy Eup-Elo outside training [1.8190, 2.5898] eV`. dE IS the wavelength, so a redder line is below the box by construction. The ENGINE-A-3DNLTE route that produced the four live Reference products cannot reach here without extrapolating the network outside its own domain, which classify_line refuses. MEASURED: amarsi_domain_check.csv. |
| solar_iag | NIR | Codex Grade | **LIVE** | 27 | 27 primary-lab-gf Fe line(s) at this depth grade (<=0.6) in 9199-10650 A |
| solar_iag | NIR | Deep Grade | **ABSENT_NO_LAB_LINE_AT_THIS_DEPTH** | 0 | no primary-lab-gf Fe line with feature depth >0.6 in 9199-10650 A (deepest lab feature in band = 0.409) |
| solar_iag | NIR | Reference Grade | **ABSENT_ENGINE_DOMAIN** | 1 | 1 AGSS21 line(s) fall in NIR, but ALL of the 1 AGSS21 Fe I line(s) here are OUTSIDE the Amarsi 2022 MLP's training box — every one fails `transition energy Eup-Elo outside training [1.8190, 2.5898] eV`. dE IS the wavelength, so a redder line is below the box by construction. The ENGINE-A-3DNLTE route that produced the four live Reference products cannot reach here without extrapolating the network outside its own domain, which classify_line refuses. MEASURED: amarsi_domain_check.csv. |
| solar_crires_plus_y_wide_rya1054 | near-UV | Codex Grade | **ABSENT_NO_COVERAGE** | 0 | this holding covers 9479-24855 A (INSTRUMENT (measured from real pixels, RYA-805)) and the near-UV band is 3000-3780 A — they do not overlap. No spectrum exists to measure. |
| solar_crires_plus_y_wide_rya1054 | near-UV | Deep Grade | **ABSENT_NO_COVERAGE** | 0 | this holding covers 9479-24855 A (INSTRUMENT (measured from real pixels, RYA-805)) and the near-UV band is 3000-3780 A — they do not overlap. No spectrum exists to measure. |
| solar_crires_plus_y_wide_rya1054 | near-UV | Reference Grade | **ABSENT_NO_COVERAGE** | 0 | this holding covers 9479-24855 A (INSTRUMENT (measured from real pixels, RYA-805)) and the near-UV band is 3000-3780 A — they do not overlap. No spectrum exists to measure. |
| solar_crires_plus_y_wide_rya1054 | VIS | Codex Grade | **ABSENT_NO_COVERAGE** | 0 | this holding covers 9479-24855 A (INSTRUMENT (measured from real pixels, RYA-805)) and the VIS band is 3780-6910 A — they do not overlap. No spectrum exists to measure. |
| solar_crires_plus_y_wide_rya1054 | VIS | Deep Grade | **ABSENT_NO_COVERAGE** | 0 | this holding covers 9479-24855 A (INSTRUMENT (measured from real pixels, RYA-805)) and the VIS band is 3780-6910 A — they do not overlap. No spectrum exists to measure. |
| solar_crires_plus_y_wide_rya1054 | VIS | Reference Grade | **ABSENT_NO_COVERAGE** | 0 | this holding covers 9479-24855 A (INSTRUMENT (measured from real pixels, RYA-805)) and the VIS band is 3780-6910 A — they do not overlap. No spectrum exists to measure. |
| solar_crires_plus_y_wide_rya1054 | red-optical | Codex Grade | **ABSENT_NO_COVERAGE** | 0 | this holding covers 9479-24855 A (INSTRUMENT (measured from real pixels, RYA-805)) and the red-optical band is 6910-9199 A — they do not overlap. No spectrum exists to measure. |
| solar_crires_plus_y_wide_rya1054 | red-optical | Deep Grade | **ABSENT_NO_COVERAGE** | 0 | this holding covers 9479-24855 A (INSTRUMENT (measured from real pixels, RYA-805)) and the red-optical band is 6910-9199 A — they do not overlap. No spectrum exists to measure. |
| solar_crires_plus_y_wide_rya1054 | red-optical | Reference Grade | **ABSENT_NO_COVERAGE** | 0 | this holding covers 9479-24855 A (INSTRUMENT (measured from real pixels, RYA-805)) and the red-optical band is 6910-9199 A — they do not overlap. No spectrum exists to measure. |
| solar_crires_plus_y_wide_rya1054 | NIR | Codex Grade | **LIVE** | 13 | 13 primary-lab-gf Fe line(s) at this depth grade (<=0.6) in 9479-13000 A |
| solar_crires_plus_y_wide_rya1054 | NIR | Deep Grade | **ABSENT_NO_LAB_LINE_AT_THIS_DEPTH** | 0 | no primary-lab-gf Fe line with feature depth >0.6 in 9479-13000 A (deepest lab feature in band = 0.409) |
| solar_crires_plus_y_wide_rya1054 | NIR | Reference Grade | **ABSENT_ENGINE_DOMAIN** | 1 | 1 AGSS21 line(s) fall in NIR, but ALL of the 1 AGSS21 Fe I line(s) here are OUTSIDE the Amarsi 2022 MLP's training box — every one fails `transition energy Eup-Elo outside training [1.8190, 2.5898] eV`. dE IS the wavelength, so a redder line is below the box by construction. The ENGINE-A-3DNLTE route that produced the four live Reference products cannot reach here without extrapolating the network outside its own domain, which classify_line refuses. MEASURED: amarsi_domain_check.csv. |

## Engine axis — band x treatment

| band | treatment | verdict | reason |
|---|---|---|---|
| near-UV | 1D-LTE | **LIVE** | product exists in the feed |
| near-UV | ENGINE-A | **LIVE** | product exists in the feed |
| near-UV | ENGINE-A-3DNLTE | **ABSENT_ENGINE_BOXED** | the Amarsi MLP is the 3D-NLTE route and every product it has ever emitted is VIS; RYA-1106 ran it on AGSS21's optical set. Whether its grid reaches this band is NOT decidable from this checkout — the network's own domain check is on STELLAR parameters, not wavelength (pipeline/amarsi3d.py). Reported as owed-verification, not as reach. |
| near-UV | ENGINE-B | **GAP** | near-UV has gradeable Fe lines and this treatment is not NLTE-gated here, so a product is applicable and missing |
| near-UV | ENGINE-B-NLTE | **UNDETERMINED_IN_CHECKOUT** | ispec_nearuv_3000_3780/atomic_lines.tsv is deliberately not committed (regenerable on Sirius, 12 MB). NLTE labelling cannot be read here. |
| near-UV | synth-1D-LTE-gerber | **GAP** | near-UV has gradeable Fe lines and this treatment is not NLTE-gated here, so a product is applicable and missing |
| near-UV | synth-mean3D-LTE-gerber-stagger | **GAP** | near-UV has gradeable Fe lines and this treatment is not NLTE-gated here, so a product is applicable and missing |
| near-UV | synth-mean3D-NLTE-gerber-stagger | **UNDETERMINED_IN_CHECKOUT** | ispec_nearuv_3000_3780/atomic_lines.tsv is deliberately not committed (regenerable on Sirius, 12 MB). NLTE labelling cannot be read here. |
| VIS | 1D-LTE | **LIVE** | product exists in the feed |
| VIS | ENGINE-A | **LIVE** | product exists in the feed |
| VIS | ENGINE-A-3DNLTE | **LIVE** | product exists in the feed |
| VIS | ENGINE-B | **LIVE** | product exists in the feed |
| VIS | ENGINE-B-NLTE | **LIVE** | product exists in the feed |
| VIS | synth-1D-LTE-gerber | **LIVE** | product exists in the feed |
| VIS | synth-mean3D-LTE-gerber-stagger | **LIVE** | product exists in the feed |
| VIS | synth-mean3D-NLTE-gerber-stagger | **LIVE** | product exists in the feed |
| red-optical | 1D-LTE | **LIVE** | product exists in the feed |
| red-optical | ENGINE-A | **LIVE** | product exists in the feed |
| red-optical | ENGINE-A-3DNLTE | **ABSENT_ENGINE_BOXED** | the Amarsi MLP is the 3D-NLTE route and every product it has ever emitted is VIS; RYA-1106 ran it on AGSS21's optical set. Whether its grid reaches this band is NOT decidable from this checkout — the network's own domain check is on STELLAR parameters, not wavelength (pipeline/amarsi3d.py). Reported as owed-verification, not as reach. |
| red-optical | ENGINE-B | **GAP** | red-optical has gradeable Fe lines and this treatment is not NLTE-gated here, so a product is applicable and missing |
| red-optical | ENGINE-B-NLTE | **UNDETERMINED_IN_CHECKOUT** | binds `ispec_ges_v6`, an iSpec-bundled list resolved at runtime and not present in the repo. Same list as VIS, where labels ARE present, but the red-optical SLICE cannot be confirmed from this checkout. |
| red-optical | synth-1D-LTE-gerber | **GAP** | red-optical has gradeable Fe lines and this treatment is not NLTE-gated here, so a product is applicable and missing |
| red-optical | synth-mean3D-LTE-gerber-stagger | **GAP** | red-optical has gradeable Fe lines and this treatment is not NLTE-gated here, so a product is applicable and missing |
| red-optical | synth-mean3D-NLTE-gerber-stagger | **UNDETERMINED_IN_CHECKOUT** | binds `ispec_ges_v6`, an iSpec-bundled list resolved at runtime and not present in the repo. Same list as VIS, where labels ARE present, but the red-optical SLICE cannot be confirmed from this checkout. |
| NIR | 1D-LTE | **LIVE** | product exists in the feed |
| NIR | ENGINE-A | **LIVE** | product exists in the feed |
| NIR | ENGINE-A-3DNLTE | **ABSENT_ENGINE_BOXED** | the Amarsi MLP is the 3D-NLTE route and every product it has ever emitted is VIS; RYA-1106 ran it on AGSS21's optical set. Whether its grid reaches this band is NOT decidable from this checkout — the network's own domain check is on STELLAR parameters, not wavelength (pipeline/amarsi3d.py). Reported as owed-verification, not as reach. |
| NIR | ENGINE-B | **GAP** | NIR has gradeable Fe lines and this treatment is not NLTE-gated here, so a product is applicable and missing |
| NIR | ENGINE-B-NLTE | **ABSENT_NO_NLTE_LABELS_IN_BAND_LINELIST** | data/linelists/ispec_ir_9200_13000/atomic_lines.tsv carries 237 Fe I rows and `nlte` reads 'F' on every one. An NLTE engine here would be silently synthesised in LTE (RYA-764), so the NLTE treatments are NOT applicable. |
| NIR | synth-1D-LTE-gerber | **GAP** | NIR has gradeable Fe lines and this treatment is not NLTE-gated here, so a product is applicable and missing |
| NIR | synth-mean3D-LTE-gerber-stagger | **GAP** | NIR has gradeable Fe lines and this treatment is not NLTE-gated here, so a product is applicable and missing |
| NIR | synth-mean3D-NLTE-gerber-stagger | **ABSENT_NO_NLTE_LABELS_IN_BAND_LINELIST** | data/linelists/ispec_ir_9200_13000/atomic_lines.tsv carries 237 Fe I rows and `nlte` reads 'F' on every one. An NLTE engine here would be silently synthesised in LTE (RYA-764), so the NLTE treatments are NOT applicable. |
