# Line uncertainty verification

[line_fit_terms.csv](line_fit_terms.csv) records statistical curvature terms from the actual LTE profile fits. These are individual fit terms. Systematic terms and profile adequacy require the calculations and dispositions below.

| Holding and product | LTE abundance | Assembled uncertainty | Status |
|---|---:|---:|---|
| Kitt Peak molecfit Fe I 8293.513 Å | 7.452 | 0.0345 dex | Assembled standard budget. Canonical strength uncertainty is unavailable, and the −0.039 dex reference scale shift is outside this budget. |
| Kurucz 2005 Fe II red-optical three-line pool | 7.506 | 0.0544 dex | Assembled standard budget. Individual canonical strength uncertainties are unavailable, and telluric eligibility is pending. |
| Kitt Peak molecfit N I 8223/8718 Å | 8.022 | — | Carbon blend legs running. Fit term 0.0248 dex, transition data 0.0208 dex and atmosphere 0.0343 dex. |

The four Lawler et al. Cr II near-UV line-fit terms range from 0.0474 to 0.0800 dex and the primary table gives transition-strength uncertainties of 0.0460, 0.0405, 0.0262 and 0.0332 dex in increasing wavelength order. The poor profiles hold all four lines independently of those terms.

Cr I 4829.309 and 4829.372 Å are distinct Sobeck et al. laboratory transitions with lower J of 4 and 5, sharing upper J 4. The lower-J line carries 0.0629 dex in the source and the stronger line carries 0.0580 dex. The archived fit's curvature term is 0.0101 dex. The canonical species/EP/wavelength clustering can treat the pair as one physical line, so its printed fit result cannot be assigned to an isolated oscillator strength. A diagnostic restores both lab strengths by the transition's J and energy and perturbs both together to test the shared error. A bound under correlated shifts is needed before a recommendation on this feature.

## Iron strength evidence

The [source comparison](fe_strength_source_comparison.csv) records a direct check of Meléndez and Barbuy Table 1. The three Fe II strengths are the exact values used in AGSS21 and each is marked L. The paper normalizes a multiplet with laboratory evidence and uses calculated relative strengths within the multiplet. The table has no individual sigma column. The 7222 and 7224 Å deck rows share the upper level, so a common normalization must not be treated as two independent terms. The checked [open author preprint](https://arxiv.org/abs/0901.4451) has SHA-256 `2b6eef11ef23cbc6159d69aee05f8424c2c36eda17058036e7316775b3812911`.

Scott et al. Table 1 flags the Fe I 8293.515 Å row for large oscillator-strength uncertainty and downweights it. Treat the initial 0.0345 dex budget as fit plus atmosphere. The [isolated Fe scale comparison](../../run_fe_redopt_literature_scale.py) tests how the published log gf affects these same fits. The first attempt stopped before fitting because the local adapter forwarded a duplicated loader keyword. The v2 run uses the corrected adapter.

## Chromium two-transition check

The active synthesis overlay subtracts 0.050390 dex from both Sobeck values, which verifies the grouping error. The named fit selector is kept at 4829.309 Å. The observed fit window contains both transitions and changes the element abundance jointly. The control reproduces 5.561 dex exactly. Restoring both distinct lab strengths returns 5.511 dex and a statistical fit term of 0.0101 dex. The difference of −0.050 dex is a model-input correction. The separate uncertainty sensitivity still uses the published errors on both transition strengths.
