#!/usr/bin/env python3
"""
scripts/build_nlte_inventory_beyond28_rya758.py
===============================================
RYA-758 — NLTE grid inventory for the elements OUTSIDE the canonical 28.

Emits ``data/audit/nlte_grid_inventory_beyond28.csv``: one row per naturally
occurring element from Be (Z=4) to U (Z=92) that is NOT in the canonical set
(the 26 symbols of ``config.constants.TARGET_ELEMENTS`` plus Zn, which RYA-757
ratified as the 28th species).  Tc (Z=43) and Pm (Z=61) are excluded — no stable
isotope, no solar photospheric line.

    python scripts/build_nlte_inventory_beyond28_rya758.py

EVERY grid claim in this table carries its source.  The classifications come from
five sources read live on 2026-08-09, never from memory:

 1. Lodders, Bergemann & Palme 2025, Space Sci. Rev. 221, 23 (accepted 2025-02-07)
    -- Table 2 quality index + the per-element text.  Local copy:
    ``Reference documents/AbundancesLodders.pdf``.
 2. Grevesse, Scott, Asplund & Sauval 2015, A&A 573, A27 (Paper III, Cu to Th)
    -- the underlying diagnostic-line list, equivalent widths and the handful of
    published NLTE corrections for heavy elements (arXiv:1405.0288).
 3. MPIA Spectrum Tools (https://nlte.mpia.de) -- the served species were
    ENUMERATED from the live HTML forms, not assumed.
 4. INASAN nLTE2 (https://spectrum.inasan.ru/nLTE2/) -- likewise enumerated.
 5. The Amarsi group's Zenodo "Grid/NLTE" family (concept DOI
    10.5281/zenodo.3888393) -- all 8 versions enumerated through the Zenodo API.

Plus a post-Lodders arXiv sweep (>= 2025-03-01), whose raw hit list is written
alongside as ``data/audit/nlte_grid_inventory_beyond28_arxiv_scan.json`` by
``scripts/scan_arxiv_post_lodders_rya758.py``.

Column semantics that are NOT self-evident:

``grid_tier``      the best PUBLISHED line-formation treatment for the element's
                   solar diagnostic, independent of whether we could obtain the
                   grid file.  ``not_measurable`` means there is no photospheric
                   diagnostic at all (noble gases; sunspot-only species).
``codex_fit_verdict``
                   what it would take for THIS pipeline to derive the element.
                   ``blocked_no_grid`` = the line exists, the grid does not.
                   ``blocked_spectroscopy`` = the grid question is moot because
                   the diagnostic is out of reach (far-UV below the atmospheric
                   cutoff, sunspot-only, or no identified solar line).
``canonical_criterion_gates_passed``
                   the four-gate science rule: alpha / CNO / bio-significant /
                   neutron-capture tracer.  ``bio_significant`` is a codex-
                   internal classification (an established essential role in
                   terrestrial biochemistry, named per row in ``notes``), NOT a
                   quantity lifted from a paper -- see the report.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / 'data' / 'audit' / 'nlte_grid_inventory_beyond28.csv'
TICKET = 'RYA-758'

COLUMNS = [
    'element', 'atomic_number', 'ion_state_diagnostic', 'grid_tier',
    'best_grid_reference', 'grid_reference_doi', 'grid_url_or_source',
    'bio_significant', 'alpha_element', 'cno', 'neutron_capture_tracer',
    'spectroscopic_diagnostic_lines', 'canonical_criterion_gates_passed',
    'codex_fit_verdict', 'notes', 'provenance_ticket',
]
GRID_TIERS = {'3D_NLTE', '1D_NLTE', 'LTE_only', 'not_measurable'}
VERDICTS = {'cheap_to_adopt', 'feasible_high_effort', 'blocked_no_grid',
            'blocked_spectroscopy', 'not_measurable'}

# Recurring source strings, written once so a row can never drift from its twin.
LBP25 = 'Lodders, Bergemann & Palme 2025 (Space Sci. Rev. 221, 23)'
LBP25_DOI = '10.1007/s11214-025-01146-w'
G15 = 'Grevesse, Scott, Asplund & Sauval 2015 (A&A 573, A27)'
G15_DOI = '10.1051/0004-6361/201424111'
NO_GRID = 'no grid published'
NO_SERVICE = (
    'not served by MPIA Spectrum Tools, INASAN nLTE2, or the Amarsi Zenodo '
    'Grid/NLTE family (all three enumerated live 2026-08-09)')
NO_PHOT = (
    f'{LBP25} Table 2 carries no photospheric value for this element (entry '
    '"..."); the solar system abundance is meteoritic (CI-chondrite)')
REE_NO_NLTE = (
    f'{G15} Sect. 4.21: "We are not aware of any other NLTE study of rare Earth '
    'elements in the Sun" (the one exception being Eu, which is canonical).')

# (element, Z, ion, tier, ref, doi, url, bio, alpha, cno, ncap, lines, verdict, notes)
ROWS: list[tuple] = [
    # ---- light elements -----------------------------------------------------
    ('Be', 4, 'Be II', '3D_NLTE',
     'Amarsi, Ogneva, Buldgen, Grevesse, Zhou & Barklem 2024 (A&A 690, A128)',
     '10.1051/0004-6361/202451778',
     'CDS VizieR J/A+A/690/A128 (Table 2); NOT in the Amarsi Zenodo Grid/NLTE '
     'family (concept 10.5281/zenodo.3888393, all 8 versions enumerated '
     '2026-08-09) -- departure grid acquisition is a direct request',
     False, False, False, False,
     '3131.07 [near-UV]', 'feasible_high_effort',
     'THE only outside-28 element Lodders 2025 accepts on a full 3D NLTE '
     'analysis -- and they still grade it C (sigma 0.14, up from Amarsi\'s '
     '0.05). Three named reasons: one single Be II line, a strong unidentified '
     'blend at 3131.02 A, and a UV radiation field that forms in the '
     'chromosphere, which no 1D or 3D photospheric model carries. The line is '
     'inside the Kitt Peak near-UV arm, so the spectroscopy is reachable; the '
     'modelling is what is hard.'),

    ('B', 5, 'B I', 'LTE_only',
     'Cunha & Smith 1999 (ApJ 512, 1006), adopted by ' + LBP25,
     '10.1086/306796',
     NO_GRID,
     True, False, False, False,
     '2496.8 [far-UV]; 2089.6 [far-UV]', 'blocked_spectroscopy',
     'Quality index E in Lodders 2025 (sigma 0.25) -- jointly the worst-graded '
     'element in Table 2. Both diagnostics sit below the atmospheric cutoff, so '
     'this is HST/STIS-only. POST-LODDERS: Spite, Barbuy & Tan 2025 '
     '(arXiv:2510.11594) add a measurable B I 2089.6 A far-UV line in metal-poor '
     'stars -- also STIS. Bio-significant: essential micronutrient in vascular '
     'plants (borate cross-linking of cell-wall rhamnogalacturonan-II).'),

    ('F', 9, 'HF (1-0) R9 vibration-rotation; no F I photospheric line',
     'not_measurable',
     'no NLTE model located; the solar value is a SUNSPOT UMBRA measurement '
     '(Maiorca et al. 2014, ApJ 788, 149)',
     '10.1088/0004-637X/788/2/149',
     'no grid; and no photospheric diagnostic to apply one to',
     False, False, False, False,
     '23358.3 [IR, sunspot umbra only]', 'blocked_spectroscopy',
     'RESOLVES the F line of enquiry. Lodders 2025 Table 2 flags F "sunspot" '
     '(quality D) and defers to Lodders & Fegley 2023 (Geochemistry 83, 125957). '
     'No F line is detectable in the quiet-Sun photospheric spectrum at all: HF '
     'has a dissociation energy of 5.87 eV, so the molecule only survives in cool '
     'umbrae. NLTE is NOT the blocker -- the photospheric diagnostic does not '
     'exist, on the Sun or on any of our FGK targets. F is not a deferred '
     'element, it is an out-of-scope one.'),

    ('Cl', 17, 'Cl I (sunspot)', 'not_measurable',
     f'{LBP25}, deferring to Lodders & Fegley 2023 (Geochemistry 83, 125957)',
     LBP25_DOI,
     'no grid; sunspot-only diagnostic',
     True, False, False, False,
     'sunspot umbral lines; no quiet-Sun photospheric line', 'blocked_spectroscopy',
     'Lodders 2025 Table 2 flags Cl "sunspot" (quality D, sigma 0.20). Same '
     'structural blocker as F. Bio-significant: chloride is a required '
     'electrolyte and the substrate of the photosystem-II oxygen-evolving '
     'complex.'),

    # ---- noble gases: no photospheric route at all --------------------------
    ('Ne', 10, 'no photospheric diagnostic (FIP-biased coronal/SW only)',
     'not_measurable',
     f'{LBP25} Table 3 -- derived from the Ne/O ratio of Young 2018 for the '
     'quiet transition region times the recommended photospheric O',
     LBP25_DOI, 'no grid; derived quantity, not a line measurement',
     False, True, False, False, '', 'not_measurable',
     'Alpha element, but structurally unmeasurable in an FGK photosphere: Ne has '
     'no photospheric lines at solar Teff. Quality index D. Any codex Ne would be '
     'a propagated ratio, not a measurement -- outside what this pipeline claims '
     'to do.'),

    ('Ar', 18, 'no photospheric diagnostic (FIP-biased coronal/SW only)',
     'not_measurable',
     f'{LBP25} Table 3 -- adopted from Lodders 2008 via element ratios',
     LBP25_DOI, 'no grid; derived quantity, not a line measurement',
     False, True, False, False, '', 'not_measurable',
     'Alpha element; same structural blocker as Ne. Quality index D. Lodders '
     'note the adopted 6.50 sits above solar-wind-based values near 6.4, i.e. the '
     'derivation route itself is the dominant error.'),

    ('Kr', 36, 'no photospheric diagnostic', 'not_measurable',
     f'{LBP25} Table 3 -- interpolated from s-process nuclide systematics '
     '(Prantzos et al. 2020)',
     LBP25_DOI, 'no grid; nucleosynthesis interpolation, not a measurement',
     False, False, False, True, '', 'not_measurable',
     'Quality index D. An s-process tracer in principle, but the "measurement" is '
     'an s-process model interpolation -- circular for any codex use that would '
     'test s-process predictions.'),

    ('Xe', 54, 'no photospheric diagnostic', 'not_measurable',
     f'{LBP25} Table 3 -- interpolated from s-process nuclide systematics '
     '(Prantzos et al. 2020)',
     LBP25_DOI, 'no grid; nucleosynthesis interpolation, not a measurement',
     False, False, False, True, '', 'not_measurable',
     'Same as Kr: quality index D, model-interpolated rather than measured.'),

    # ---- trans-Fe, first row ------------------------------------------------
    ('Ga', 31, 'Ga I', 'LTE_only', f'{G15} 3D LTE', G15_DOI,
     NO_GRID + '; ' + NO_SERVICE,
     False, False, False, True, '4172.05 [optical]', 'blocked_no_grid',
     'Quality index C in Lodders 2025 (sigma 0.14), value carried unchanged from '
     'Asplund et al. 2021. One usable line. No Ga model atom located.'),

    ('Ge', 32, 'Ge I', 'LTE_only', f'{G15} 3D LTE', G15_DOI,
     NO_GRID + '; ' + NO_SERVICE,
     False, False, False, True, '3269.49 [near-UV]', 'blocked_no_grid',
     'Quality index C. One near-UV line. No Ge model atom located.'),

    ('As', 33, 'no identified solar photospheric line', 'not_measurable',
     NO_PHOT, LBP25_DOI, 'no grid; no diagnostic',
     False, False, False, True, '', 'blocked_spectroscopy',
     'Lodders 2025 explicitly call out As as a siderophile whose CI-chondrite '
     'analyses need refinement -- the solar side has nothing to refine against.'),

    ('Se', 34, 'no identified solar photospheric line', 'not_measurable',
     NO_PHOT, LBP25_DOI, 'no grid; no diagnostic',
     True, False, False, True, '', 'blocked_spectroscopy',
     'Bio-significant: selenocysteine (the 21st proteinogenic amino acid) and the '
     'glutathione-peroxidase family. Strong bio case, zero spectroscopic route.'),

    ('Br', 35, 'no identified solar photospheric line', 'not_measurable',
     NO_PHOT + '; halogens deferred to Lodders & Fegley 2023', LBP25_DOI,
     'no grid; no diagnostic',
     True, False, False, True, '', 'blocked_spectroscopy',
     'Bio-significant: bromide is required for peroxidasin-catalysed sulfilimine '
     'cross-links in collagen IV basement membranes. Same halogen wall as F/Cl, '
     'one step worse -- not even a sunspot value.'),

    ('Rb', 37, 'Rb I', '1D_NLTE',
     'Korotin 2020 (Astron. Lett. 46, 541) -- 29-level Rb I model atom; solar '
     'NLTE A(Rb) = 2.35 +/- 0.05, NLTE correction -0.12 dex',
     '10.1134/S1063773720080022',
     'no public grid server: ' + NO_SERVICE + '. Corrections are tabulated in the '
     'paper; a departure grid would be a direct request to the author.',
     False, False, False, True,
     '7800.27 [red-optical]; 7947.60 [red-optical]', 'feasible_high_effort',
     'SEED-TABLE REVISION: the first-pass table called Rb cheap_to_adopt. It is '
     'not cheap on two counts. (1) No public grid -- the three public NLTE '
     'services were enumerated and none serves Rb. (2) Both resonance lines sit '
     'past 7000 A, i.e. outside the HARPS arm that carries our entire measured EW '
     'pool; solar Rb needs the Kitt Peak red arm and every other target needs a '
     'red-capable instrument. Genuinely the best-motivated outside-28 n-capture '
     'candidate (s-process branching at 85Kr/86Rb), and the model-atom physics is '
     'done -- the work is acquisition, not modelling.'),

    ('Nb', 41, 'Nb II', 'LTE_only',
     f'{G15} -- 3D correction applied to the HM-model results of Nilsson et al. '
     '2010', G15_DOI, NO_GRID + '; ' + NO_SERVICE,
     False, False, False, True,
     '3194.97 [near-UV]; 3215.59 [near-UV]; 3717.06 [near-UV]; 3740.07 [near-UV]',
     'blocked_no_grid',
     'Quality index C. Four near-UV Nb II lines, all analysed only differentially '
     'against a 1D semi-empirical model.'),

    ('Mo', 42, 'Mo I', 'LTE_only',
     'no Mo model atom exists. Mishenina, Kurtukian-Nieto, Gorbaneva, Amarsi, '
     'Psaltis & Pignatari 2026 (A&A 705, A38; arXiv:2511.21190), the most recent '
     'dedicated Mo study: "we emphasize that defining NLTE corrections for Mo and '
     'Ru abundance measurements is essential, as these remain currently unknown"',
     '10.48550/arXiv.2511.21190',
     NO_GRID + '; ' + NO_SERVICE,
     True, False, False, True,
     '5506.50 [optical]; 5533.03 [optical]', 'blocked_no_grid',
     'THE strongest bio case outside the 28 (nitrogenase FeMo-cofactor and the '
     'molybdopterin enzyme family -- Mo is where biological N2 fixation happens), '
     'and the ONLY outside-28 element that clears two science gates. Blocked '
     'purely by atomic physics. Note the trap: Mishenina et al. quote a +0.15 dex '
     'MPIA correction for the ANALOGOUS Cr I 5208.41/5206.02 lines as a proxy, '
     'NOT for Mo -- MPIA Spectrum Tools was enumerated 2026-08-09 and serves no '
     'Mo. Both Mo I lines are optical and inside the HARPS arm, so the day a Mo '
     'model atom appears, Mo is immediately workable.'),

    ('Ru', 44, 'Ru I', 'LTE_only',
     'no Ru model atom exists (Mishenina et al. 2026, A&A 705, A38)',
     '10.48550/arXiv.2511.21190', NO_GRID + '; ' + NO_SERVICE,
     False, False, False, True,
     '3436.74 [near-UV]; 3498.95 [near-UV]; 3742.28 [near-UV]; 4080.60 [optical]; '
     '4554.52 [optical]; 4584.44 [optical]', 'blocked_no_grid',
     'Quality index C-. Six lines is the richest diagnostic set of any blocked '
     'outside-28 element; Mishenina et al. judge Ru\'s LTE departures likely '
     'insignificant relative to other errors, which makes it the least harmful '
     'LTE-only case here.'),

    ('Rh', 45, 'Rh I', 'LTE_only', f'{G15} 3D LTE', G15_DOI,
     NO_GRID + '; ' + NO_SERVICE,
     False, False, False, True,
     '3434.89 [near-UV]; 3692.36 [near-UV]', 'blocked_no_grid',
     'Quality index C-. Low-excitation lines of a minority species -- exactly the '
     'configuration Grevesse et al. flag as most NLTE-prone, and there is no grid '
     'to check it with.'),

    ('Pd', 46, 'Pd I', '1D_NLTE',
     'Mashonkina & Smogorzhevskii 2025 (A&A 703, A296; arXiv:2510.06968) -- new '
     'comprehensive Pd I model atom; solar NLTE log eps = 1.61 to 1.70',
     '10.48550/arXiv.2510.06968',
     'no public grid located: INASAN nLTE2 was enumerated 2026-08-09 (11 species) '
     'and Pd is not among them -- a request to the authors is owed',
     False, False, False, True,
     '3242.70 [near-UV]; 3404.58 [near-UV]', 'feasible_high_effort',
     'POST-LODDERS (Oct 2025): promotes Pd from the seed table\'s LTE_only. The '
     'authors\' own solar NLTE value spans 1.61-1.70 depending on line, i.e. a '
     '0.09 dex internal spread that already eats our +/-0.10 gate -- adopt only '
     'with that spread propagated.'),

    ('Ag', 47, 'Ag I', '3D_NLTE',
     'Caliskan, Amarsi, Jonsson, Grevesse & Sahoo 2026 (A&A 711, A155; '
     'arXiv:2605.05356) -- first Ag I NLTE model atom, ab-initio MCHF oscillator '
     'strengths; solar 3D NLTE A(Ag) = 1.15 +/- 0.08',
     '10.1051/0004-6361/202659578',
     'https://doi.org/10.5281/zenodo.20037437 -- Amarsi Grid/NLTE v8 (2026-05-05): '
     'nlte_Ag_caliskan_Jan2026_pysme.grd (960 MB), _idlsme.grd, text form, '
     'atmos_Ag.txt, label_Ag.txt; CC-BY-4.0, 4.4 GB total',
     False, False, False, True,
     '3280.68 [near-UV]; 3382.90 [near-UV]', 'cheap_to_adopt',
     'CRITICAL FINDING. A full 3D NLTE outside-28 element that post-dates Lodders '
     '2025 and was NOT in the first-pass table. It is cheap for us specifically: '
     'the departure grid is the SAME PySME .grd family as our already-registered '
     'Cu_Caliskan2024_PySME (Zenodo 15062813) and the Amarsi-2020 grids, from the '
     'same first author, so the derivation path is the one RYA-540 already ran. '
     'Both resonance lines are near-UV, inside the Kitt Peak arm (Grevesse et al. '
     'measure 35 and 22 mA). The new value is +0.19 dex above the previous '
     'reference and cuts the photosphere-meteorite discrepancy from 0.25 to 0.06 '
     'dex. Grid goes on Sirius, never the Mac.'),

    ('Cd', 48, 'Cd I', 'LTE_only', f'{G15} 3D LTE', G15_DOI,
     NO_GRID + '; ' + NO_SERVICE,
     False, False, False, True, '5085.82 [optical]', 'blocked_no_grid',
     'Quality index D. A single line at 7.3 mA equivalent width -- at or below '
     'our own measurable-EW floor even before the missing grid.'),

    ('In', 49, 'In I', 'not_measurable',
     f'{LBP25} Table 2 flags In "Sunspot"; {G15} show the 4511.3 A photospheric '
     'feature is NOT In I (Vitas et al. 2008, MNRAS 384, 370)',
     LBP25_DOI, 'no grid; the photospheric identification itself is refuted',
     False, False, False, True,
     '4511.3 [optical, identification refuted]', 'blocked_spectroscopy',
     'A cautionary row: In carries a tabulated solar abundance (0.80) that rests '
     'on a sunspot measurement, while the photospheric feature historically '
     'attributed to In I is contributed <20% by In I. Quality index D.'),

    ('Sn', 50, 'Sn I', 'LTE_only', f'{G15} 3D LTE', G15_DOI,
     NO_GRID + '; ' + NO_SERVICE,
     False, False, False, True, '3801.02 [near-UV]', 'blocked_no_grid',
     'Quality index C-. One near-UV line at 1.2 mA.'),

    ('Sb', 51, 'no identified solar photospheric line', 'not_measurable',
     NO_PHOT, LBP25_DOI, 'no grid; no diagnostic',
     False, False, False, True, '', 'blocked_spectroscopy',
     'Lodders 2025 raise the CI-chondrite Sb abundance but have no solar value to '
     'compare it against.'),

    ('Te', 52, 'no identified solar photospheric line', 'not_measurable',
     NO_PHOT, LBP25_DOI, 'no grid; no diagnostic',
     False, False, False, True, '', 'blocked_spectroscopy',
     'The post-Lodders arXiv sweep returns Te atomic-structure work '
     '(arXiv:2510.17357, Te I-V photoionization) but it targets kilonova ejecta, '
     'not a solar photospheric diagnostic.'),

    ('I', 53, 'no identified solar photospheric line', 'not_measurable',
     NO_PHOT + '; halogens deferred to Lodders & Fegley 2023', LBP25_DOI,
     'no grid; no diagnostic',
     True, False, False, True, '', 'blocked_spectroscopy',
     'Bio-significant: iodine is the functional atom of the thyroid hormones T3 '
     'and T4. Clears two gates on paper; has no solar line at all.'),

    ('Cs', 55, 'no identified solar photospheric line', 'not_measurable',
     NO_PHOT, LBP25_DOI, 'no grid; no diagnostic',
     False, False, False, True, '', 'blocked_spectroscopy',
     'No solar photospheric value in Lodders 2025 Table 2.'),

    # ---- rare earths --------------------------------------------------------
    ('La', 57, 'La II', 'LTE_only',
     'Lawler, Bonvallet & Sneden 2001a (ApJ 556, 452) 1D LTE on 14 La II lines, '
     f'corrected by {G15} for the 3D-HM difference (-0.03 dex)',
     '10.1086/321549', NO_GRID + '. ' + REE_NO_NLTE,
     False, False, False, True,
     '14 La II near-UV/optical lines (Lawler et al. 2001a); not individually '
     'enumerated in Grevesse et al. 2015', 'blocked_no_grid',
     'Quality index C-. Lodders 2025 state plainly: "No estimates of NLTE '
     'abundances are available for La."'),

    ('Ce', 58, 'Ce II', 'LTE_only',
     'Lawler, Sneden, Cowan, Ivans & Den Hartog 2009 (ApJS 182, 51) 1D LTE, '
     f'3D-HM corrected by {G15}',
     '10.1088/0067-0049/182/1/51', NO_GRID + '. ' + REE_NO_NLTE,
     False, False, False, True,
     'Ce II lines per Lawler et al. 2009; not individually enumerated in '
     'Grevesse et al. 2015', 'blocked_no_grid',
     'Quality index C-, value unchanged since Grevesse et al. 2015.'),

    ('Pr', 59, 'Pr II', '1D_NLTE',
     'Mashonkina, Ryabchikova, Ryabtsev & Kildiyarova 2009 (A&A 495, 297) -- Pr '
     'II/Pr III model atom; a solar NLTE value log(Pr/H) = -11.15 +/- 0.08 is '
     'given as a by-product',
     '10.1051/0004-6361:200810258',
     'no public grid; and the published grid hull is Teff 7250-9500 K (A and Ap '
     'stars), which EXCLUDES the Sun (5772 K) and every codex target',
     False, False, False, True,
     'Pr II lines per Lawler/Sneden et al. 2009; not individually enumerated in '
     'Grevesse et al. 2015', 'blocked_no_grid',
     'SEED-TABLE CORRECTION. The first-pass table cited "1D NLTE, Lodders 2025 '
     'p703" -- there is no such citation in the paper; Lodders assert only that '
     '"the only elements with NLTE estimates of abundances in 1D available are '
     'Sr, Y, Zr, Ba, Pr, Eu, Gd, and Nd", with no reference. The only Pr NLTE '
     'work this audit could locate is an A/Ap-star study whose parameter hull '
     'does not contain any star we analyse. Tier is 1D_NLTE because a solar NLTE '
     'number exists; the verdict is blocked because a grid we could interpolate '
     'does not.'),

    ('Nd', 60, 'Nd II', '1D_NLTE',
     'Dixon, Ezzeddine, Li, Merle, Bautista & Guo 2025 (ApJ 994, 44; '
     'arXiv:2509.22811) -- new Nd I/Nd II model atom and a grid of NLTE '
     'corrections for 122 Nd II lines from the UV to the H band; solar '
     'A(Nd) = 1.44 +/- 0.05 at S_H = 0.1',
     '10.3847/1538-4357/ae0b59',
     'grid existence is stated in the paper; a public download location is NOT '
     'stated in the abstract -- VERIFICATION OWED (check the ApJ data-availability '
     'statement / Zenodo before any adoption)',
     False, False, False, True,
     '122 Nd II lines, UV through H band (Dixon et al. 2025); not enumerated in '
     'the abstract', 'feasible_high_effort',
     'POST-LODDERS (Sep 2025) and the single largest rare-earth advance since '
     'Grevesse et al. 2015. Corrections span -0.3 to +0.3 dex: positive for '
     'optical and UV lines, NEGATIVE for H-band lines -- so an IR-arm Nd would '
     'need the sign handled per band, which is exactly the RYA-708 per-band '
     'discipline. The model atom is calibrated on a Drawin S_H = 0.1 scaling, not '
     'ab-initio H collisions, which is the RYA-546 vintage concern.'),

    ('Sm', 62, 'Sm II', 'LTE_only',
     'Lawler, Den Hartog, Sneden & Cowan 2006 (ApJS 162, 227) + Rehse et al. '
     f'2006 gf-values, 3D-HM corrected by {G15}',
     '10.1086/498213', NO_GRID + '. ' + REE_NO_NLTE,
     False, False, False, True,
     'Sm II lines per Lawler et al. 2006; not individually enumerated in '
     'Grevesse et al. 2015', 'blocked_no_grid', 'Quality index C-.'),

    ('Gd', 64, 'Gd II', 'LTE_only',
     'Den Hartog, Lawler, Sneden & Cowan 2006 (ApJS 167, 292) 1D LTE, 3D-HM '
     f'corrected by {G15}',
     '10.1086/508262', NO_GRID + '. ' + REE_NO_NLTE,
     False, False, False, True,
     'Gd II lines per Den Hartog et al. 2006; not individually enumerated in '
     'Grevesse et al. 2015', 'blocked_no_grid',
     'SEED-TABLE CORRECTION. Listed as 1D NLTE on the strength of the same '
     'uncited Lodders sentence as Pr. This audit could locate NO Gd II model atom '
     'or NLTE study, solar or otherwise, and Grevesse et al. 2015 state the '
     'opposite. Recorded as LTE_only until a citation is produced.'),

    ('Tb', 65, 'Tb II', 'LTE_only',
     f'Lawler, Wickliffe, Cowley & Sneden 2001c (ApJS 137, 341), 3D-HM corrected '
     f'by {G15}', '10.1086/323001', NO_GRID + '. ' + REE_NO_NLTE,
     False, False, False, True, '3659 [near-UV]', 'blocked_no_grid',
     'Quality index C- and the weakest evidence base of the rare earths: '
     'Grevesse et al. retained ONE of the three Lawler lines, dropping the other '
     'two for blends and very wide HFS. sigma 0.10 even in the source paper.'),

    ('Dy', 66, 'Dy II', 'LTE_only',
     f'Sneden, Lawler, Cowan, Ivans & Den Hartog 2009 (ApJS 182, 80), 3D-HM '
     f'corrected by {G15}', '10.1088/0067-0049/182/1/80',
     NO_GRID + '. ' + REE_NO_NLTE,
     False, False, False, True,
     'Dy II lines per Sneden et al. 2009', 'blocked_no_grid', 'Quality index C-.'),

    ('Ho', 67, 'Ho II', 'LTE_only',
     f'Lawler, Sneden & Cowan 2004 (ApJ 604, 850), 3D-HM corrected by {G15}',
     '10.1086/382068', NO_GRID + '. ' + REE_NO_NLTE,
     False, False, False, True,
     'Ho II lines per Lawler et al. 2004', 'blocked_no_grid',
     'Quality index C-; sigma 0.11 in the source analysis. Grevesse et al. note '
     'Ho II needed its own partition functions (Bord et al.).'),

    ('Er', 68, 'Er II', 'LTE_only',
     f'Lawler, Sneden, Cowan et al. 2008 (ApJS 178, 71), 3D-HM corrected by {G15}',
     '10.1086/589834', NO_GRID + '. ' + REE_NO_NLTE,
     False, False, False, True,
     'Er II lines per Lawler et al. 2008', 'blocked_no_grid', 'Quality index C-.'),

    ('Tm', 69, 'Tm II', 'LTE_only',
     f'Sneden et al. 2009 (ApJS 182, 80), 3D-HM corrected by {G15}',
     '10.1088/0067-0049/182/1/80', NO_GRID + '. ' + REE_NO_NLTE,
     False, False, False, True,
     'Tm II lines per Sneden et al. 2009', 'blocked_no_grid', 'Quality index C-.'),

    ('Yb', 70, 'Yb II', 'LTE_only',
     f'Sneden et al. 2009 (ApJS 182, 80), 3D-HM corrected by {G15}',
     '10.1088/0067-0049/182/1/80', NO_GRID + '. ' + REE_NO_NLTE,
     False, False, False, True,
     'Yb II lines per Sneden et al. 2009', 'blocked_no_grid',
     'Quality index C-; sigma 0.11 in the source analysis.'),

    ('Lu', 71, 'Lu II', 'LTE_only',
     f'Sneden et al. 2009 (ApJS 182, 80), 3D-HM corrected by {G15}',
     '10.1088/0067-0049/182/1/80', NO_GRID + '. ' + REE_NO_NLTE,
     False, False, False, True,
     'Lu II lines per Sneden et al. 2009', 'blocked_no_grid',
     'Quality index C-; the lowest solar abundance of any rare earth (A = 0.10).'),

    ('Hf', 72, 'Hf II', 'LTE_only',
     f'Lawler, Den Hartog, Labby et al. 2007 (ApJS 169, 120), 3D-HM corrected by '
     f'{G15}', '10.1086/510368', NO_GRID + '. ' + REE_NO_NLTE,
     False, False, False, True,
     'Hf II lines per Lawler et al. 2007', 'blocked_no_grid',
     'Quality index C+ -- the best-graded of the LTE-only heavies. Grevesse et '
     'al. flag a possible gf problem: their Hf II photospheric value sits 0.14 dex '
     'above the meteoritic one.'),

    # ---- third-row transition metals and beyond -----------------------------
    ('Ta', 73, 'no identified solar photospheric line', 'not_measurable',
     NO_PHOT, LBP25_DOI, 'no grid; no diagnostic',
     False, False, False, True, '', 'blocked_spectroscopy',
     'Lodders 2025 count Ta among the refractory lithophiles whose abundances are '
     '"much improved" -- on the meteoritic side only.'),

    ('W', 74, 'W I', 'LTE_only',
     'Holweger & Werner 1982 HM-model synthesis, gf-updated and 3D-corrected by '
     f'{G15}', G15_DOI, NO_GRID + '; ' + NO_SERVICE,
     True, False, False, True,
     '4008.75 [optical]; 4843.85 [optical]', 'blocked_no_grid',
     'Quality index D, and Lodders name W explicitly as the example of "1D '
     'LTE+3D LTE correction, but heavily blended 1 diagnostic line in the blue". '
     'W II outnumbers W I ~10:1 in the photosphere but has no usable lines, so '
     'the minority species carries the whole measurement. Bio-significant: '
     'tungstoenzymes (formate dehydrogenase, aldehyde ferredoxin oxidoreductase) '
     'in anaerobic archaea and bacteria -- a real but narrow biological role.'),

    ('Re', 75, 'no identified solar photospheric line', 'not_measurable',
     NO_PHOT, LBP25_DOI, 'no grid; no diagnostic',
     False, False, False, True, '', 'blocked_spectroscopy',
     'No solar photospheric value in Lodders 2025 Table 2.'),

    ('Os', 76, 'Os I', 'LTE_only', f'{G15} 3D LTE', G15_DOI,
     NO_GRID + '; ' + NO_SERVICE,
     False, False, False, True, '3301.56 [near-UV]', 'blocked_no_grid',
     'Quality index C. Only very few faint Os I lines are identified at all; '
     'Grevesse et al. retain one.'),

    ('Ir', 77, 'Ir I', 'LTE_only', f'{G15} 3D LTE', G15_DOI,
     NO_GRID + '; ' + NO_SERVICE,
     False, False, False, True, '3220.78 [near-UV]', 'blocked_no_grid',
     'Quality index E (sigma 0.20+), jointly the worst-graded element in Lodders '
     'Table 2 with B. One heavily blended line whose continuum Grevesse et al. '
     'had to re-place by hand, on the Jungfraujoch tracing because no Kitt Peak '
     'spectrum reaches these wavelengths.'),

    ('Pt', 78, 'no identified solar photospheric line', 'not_measurable',
     NO_PHOT, LBP25_DOI, 'no grid; no diagnostic',
     False, False, False, True, '', 'blocked_spectroscopy',
     'No solar photospheric value in Lodders 2025 Table 2.'),

    ('Au', 79, 'Au I', 'LTE_only',
     f'Youssef 1986 HM-model synthesis, re-measured and 3D-corrected by {G15}',
     G15_DOI, NO_GRID + '; ' + NO_SERVICE,
     False, False, False, True, '3122.78 [near-UV]', 'blocked_no_grid',
     'Quality index D. The only useful Au feature in the solar spectrum, and it '
     'is contaminated by a Fe I line at 3122.775 A whose gf is too poor to model '
     'out. Again a Jungfraujoch-only wavelength.'),

    ('Hg', 80, 'no photospheric determination exists', 'not_measurable',
     f'{LBP25}: "photospheric Hg determinations are not available and have their '
     'own challenges"', LBP25_DOI, 'no grid; no diagnostic',
     False, False, False, True, '', 'blocked_spectroscopy',
     'Lodders 2025 devote a section to Hg and conclude the Sun would be the best '
     'source precisely because the meteoritic Hg is unreliable -- and then record '
     'that no photospheric determination exists.'),

    ('Tl', 81, 'Tl I (sunspot)', 'not_measurable',
     f'{LBP25} Table 2 flags Tl "Sunspot"; {G15}: "Tl I is present in the sunspot '
     'spectrum, but heavily blended"', LBP25_DOI,
     'no grid; sunspot-only and blended there',
     False, False, False, True,
     'sunspot Tl I, heavily blended; no quiet-Sun photospheric line',
     'blocked_spectroscopy', 'Quality index D.'),

    ('Pb', 82, 'Pb I', '1D_NLTE',
     'Mashonkina, Ryabtsev & Frebel 2012 (A&A 540, A98), "Non-LTE effects on the '
     'lead and thorium abundance determinations for cool stars" -- NLTE correction '
     '+0.12 dex at S_H = 0.1 on a MARCS model (range +0.15 to +0.07); applied by '
     f'{G15}', '10.1051/0004-6361/201218790',
     'no public grid server: ' + NO_SERVICE + '. The correction is a published '
     'per-line number, not an interpolable grid.',
     False, False, False, True, '3683.48 [near-UV]', 'feasible_high_effort',
     'SEED-TABLE CORRECTION: the first-pass table swept "Hf-Bi" into 1D LTE only. '
     'Pb is the exception -- it carries a real published 1D NLTE correction, and '
     'a large one (+0.12 dex). The blocker is the line, not the physics: Pb I '
     '3683.48 sits in the outer red wing of an extremely strong Co I + Fe I + V I '
     'blend and in the blue wing of another Fe I line, so it needs the in-window '
     'blend-fit treatment built for RYA-581/585, not an EW measurement.'),

    ('Bi', 83, 'no identified solar photospheric line', 'not_measurable',
     NO_PHOT, LBP25_DOI, 'no grid; no diagnostic',
     False, False, False, True, '', 'blocked_spectroscopy',
     'No solar photospheric value in Lodders 2025 Table 2.'),

    ('Th', 90, 'Th II', '1D_NLTE',
     'Mashonkina, Ryabtsev & Frebel 2012 (A&A 540, A98), "Non-LTE effects on the '
     'lead and thorium abundance determinations for cool stars" -- NLTE correction '
     f'+0.01 dex at S_H = 0.1 on a MARCS model; applied by {G15}',
     '10.1051/0004-6361/201218790',
     'no public grid server: ' + NO_SERVICE + '. Published per-line correction '
     'only; it is small enough (+0.01) that a grid is not the binding constraint.',
     False, False, False, True, '4019.13 [optical]', 'feasible_high_effort',
     'The cosmochronology anchor, and the clearest test case for an LTE-tolerant '
     'inclusion clause -- except that Th is not LTE-blocked at all, it is '
     'blend-blocked. The Th II 4019.13 feature is a three-component deblend: of '
     'the 0.56 pm total, Grevesse et al. attribute 0.208 pm to Co I and 0.038 pm '
     'to V I, leaving 0.314 pm for Th II, with ~20% uncertainty on that share. '
     'Lodders 2025 bracket the value ([0.08], from Caffau et al. 2011) and list '
     'it "only for reference". We already have measured solar Co and V, which is '
     'exactly what the deblend needs.'),

    ('U', 92, 'no reliable solar photospheric line', 'not_measurable',
     NO_PHOT, LBP25_DOI, 'no grid; no usable diagnostic',
     False, False, False, True, '', 'blocked_spectroscopy',
     'Lodders 2025 Table 2 gives no photospheric U value at all -- the U side of '
     'the Th/U cosmochronometer is meteoritic only. Any codex Th/U ratio would '
     'therefore be half-photospheric, half-meteoritic; say so or do not report '
     'it.'),
]


def build() -> list[dict]:
    rows = []
    for (el, z, ion, tier, ref, doi, url, bio, alpha, cno, ncap, lines,
         verdict, notes) in ROWS:
        if tier not in GRID_TIERS:
            raise ValueError(f'{el}: bad grid_tier {tier!r}')
        if verdict not in VERDICTS:
            raise ValueError(f'{el}: bad codex_fit_verdict {verdict!r}')
        # LOUD: a cheap_to_adopt row without an acquisition path is the ticket's
        # own CRITICAL condition -- refuse to emit it rather than ship a blank.
        if verdict == 'cheap_to_adopt' and not url.strip():
            raise ValueError(
                f'{el}: cheap_to_adopt with an empty grid_url_or_source -- '
                f'RYA-758 section 3 marks this CRITICAL')
        rows.append({
            'element': el,
            'atomic_number': z,
            'ion_state_diagnostic': ion,
            'grid_tier': tier,
            'best_grid_reference': ref,
            'grid_reference_doi': doi,
            'grid_url_or_source': url,
            'bio_significant': bio,
            'alpha_element': alpha,
            'cno': cno,
            'neutron_capture_tracer': ncap,
            'spectroscopic_diagnostic_lines': lines,
            'canonical_criterion_gates_passed': sum((alpha, cno, bio, ncap)),
            'codex_fit_verdict': verdict,
            'notes': notes,
            'provenance_ticket': TICKET,
        })
    return rows


def main(argv=None) -> int:
    rows = build()

    # The canonical set must not leak into an "outside the canonical set" table.
    sys.path.insert(0, str(REPO))
    from config.constants import TARGET_ELEMENTS
    canonical = set(TARGET_ELEMENTS) | {'Zn'}          # Zn per RYA-757
    overlap = canonical & {r['element'] for r in rows}
    if overlap:
        raise SystemExit(f'[LOUD] canonical elements present in the beyond-28 '
                         f'inventory: {sorted(overlap)}')

    seen = [r['element'] for r in rows]
    if len(seen) != len(set(seen)):
        raise SystemExit('[LOUD] duplicate element rows')
    if len(rows) < 45:
        raise SystemExit(f'[LOUD] only {len(rows)} rows; RYA-758 requires >= 45')

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open('w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(rows)

    tiers, verdicts = {}, {}
    for r in rows:
        tiers[r['grid_tier']] = tiers.get(r['grid_tier'], 0) + 1
        verdicts[r['codex_fit_verdict']] = verdicts.get(r['codex_fit_verdict'], 0) + 1
    print(f'wrote {OUT.relative_to(REPO)}  rows={len(rows)}')
    print('  grid_tier         :', dict(sorted(tiers.items())))
    print('  codex_fit_verdict :', dict(sorted(verdicts.items())))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
