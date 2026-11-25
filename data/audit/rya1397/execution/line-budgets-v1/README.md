# Line uncertainty verification

[line_fit_terms.csv](line_fit_terms.csv) records statistical curvature terms from the actual LTE profile fits. These are individual fit terms. Systematic terms and profile adequacy require the calculations and dispositions below.

| Holding and product | LTE abundance | Assembled uncertainty | Status |
|---|---:|---:|---|
| Kitt Peak molecfit Fe I 8293.513 Å | 7.452 | 0.0345 dex | Assembled standard budget. Canonical strength uncertainty is unavailable, and the −0.039 dex reference scale shift is outside this budget. |
| Kurucz 2005 Fe II red-optical three-line pool | 7.506 | 0.0544 dex | Assembled standard budget. Individual canonical strength uncertainties are unavailable, and telluric eligibility is pending. |
| Kitt Peak molecfit N I 8223/8718 Å | 8.022 | — | Carbon blend legs running. Fit term 0.0248 dex, transition data 0.0208 dex and atmosphere 0.0343 dex. |

The four Lawler et al. Cr II near-UV line-fit terms range from 0.0474 to 0.0800 dex and the primary table gives transition-strength uncertainties of 0.0460, 0.0405, 0.0262 and 0.0332 dex in increasing wavelength order. The poor profiles hold all four lines independently of those terms.

Cr I 4829.309 and 4829.372 Å are distinct Sobeck et al. laboratory transitions with lower J of 4 and 5, sharing upper J 4. The lower-J line carries 0.0629 dex in the source and the stronger line carries 0.0580 dex. The archived fit's curvature term is 0.0101 dex. The canonical species/EP/wavelength clustering can treat the pair as one physical line, so its printed fit result cannot be assigned to an isolated oscillator strength. A diagnostic restores both lab strengths by the transition's J and energy and perturbs both together to test the shared error. A bound under correlated shifts is needed before a recommendation on this feature.
