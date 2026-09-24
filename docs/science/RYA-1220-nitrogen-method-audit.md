# RYA-1220 — nitrogen method and reproduction audit

**15 September 2026 · Working evidence package · gate OPEN; no abundance admitted**

## Updated dependency — RYA-587 owns the uncertainty contract

The 15 September directive added to [RYA-1220](https://linear.app/ryans-adventure-zone/issue/RYA-1220) makes [RYA-587](https://linear.app/ryans-adventure-zone/issue/RYA-587) the canonical element- and product-independent uncertainty contract. RYA-587 is now merged at `c6529aab`; the nitrogen migration below uses `pipeline.uncertainty_contract` and its `rya587.uncertainty.v1` schema. The existing numerical probes remain diagnostic evidence until every applicable component is populated.

The merged contract supplies the canonical schema/API, component provenance and applicability states, absolute/differential covariance path, and product admission checks. Its audit confirms the completed RYA-1213 shell ξ outputs and ingestion/stamping; no Fe rerun is needed. Experimental mean-3D Gerber is not a completeness prerequisite unless separately ratified as publication-grade.

Nitrogen consumes the base contract unchanged and maps each applicable component to **measured / N/A / HOLD**, with evidence. Its additions are limited to N-specific blend treatment, CN/NH equilibrium and C/O covariance, and justified indicator/model covariance. Differential products require matched Solar/reference transitions or indicators and engine/model routes. The N adapter in [`pipeline/nitrogen_uncertainty.py`](../../pipeline/nitrogen_uncertainty.py) only records these extensions; the shared contract owns totals and validation. Literature, transition curation and diagnostic profile work can proceed independently.

### Reconciliation with the earlier CNO campaign

The earlier `ryandamienschmitt/rya-1214-cno-product-campaign` branch is present in the repository history and its CNO artifacts are already in this working tree. It supplied the multi-band Turbospectrum routes, molecular lists, AGSS21 indicator joins, C/O NLTE grids, CRIRES+ J/K wiring, and near-UV rejection analysis. This audit reuses those results; it does not duplicate them. The new work is the missing uncertainty layer: per-band covariance, residual continuum/profile terms, and explicit separation of accepted 3D anchors from 1D-LTE molecular diagnostics. The prior campaign's `C/O ... NOT_PROPAGATED` rows are therefore treated as known incomplete uncertainty records, not as missing CNO computations.

## Finding

The two supplied articles are the correct primary anchors for the solar atomic and molecular investigations. They do **not** establish that our current flux synthesis should reproduce a single common nitrogen abundance. Their observing geometry, diagnostics, atmospheric treatment, and blend measurements differ from our current products. The present audit establishes several concrete causes of non-comparability and one substantial measured fitting sensitivity; it does not close the nitrogen gate.

The supplied files are `aa37890-20.pdf` (Amarsi et al. 2020, atomic N) and `Amarsi21.pdf` (Amarsi et al. 2021, molecular CNO). Both were read in full. Seven additional primary papers were inspected and copied into the shared reference-documents folder. File hashes, source locations, and the detailed extraction matrix accompany this report. No oscillator strength was adjusted to match the Sun.

### Decisions supported now

1. **Retain atomic N I as a diagnostic route.** Catalog central-depth shares are not integrated absorption fractions and cannot establish that only 5.5% of a feature is nitrogen. A rejection based on that calculation is unsupported.
2. **The exact five-line reproduction is incomplete.** Four reference N I lines match the canonical list. An isolated, non-canonical Amarsi Table 1 overlay now supplies 10108.90 Å together with a Brooke 12C14N list for a diagnostic test; this is not yet an admitted line-list update.
3. **Fit-window sensitivity is measured.** At 8216.336 Å, otherwise identical IAG fits change from A(N)=7.9562 to 8.1722 when the fitting half-width changes from 0.25 to 1.1 Å. This +0.2160 dex effect demonstrates sensitivity to surrounding absorption/model mismatch. It does not identify the correct continuum or abundance by itself.
4. **IAG and the current CRIRES+ J product cannot provide a matched-line CN comparison.** IAG ends below the J product's first rest-frame pixel: zero common transitions are possible. Comparing their aggregate nitrogen values is a diagnostic-family comparison, not line-by-line validation.
5. **The existing uncertainty totals are partial.** They do not demonstrate the complete requested stack or target–solar covariance. The shared contract is now connected through a review-only migration: all 12 historical N records map to `rya587.uncertainty.v1` and remain HOLD; no live feed is changed.

## Evidence and reproduction assets

| Artifact | Purpose |
|---|---|
| [Evidence matrix](../../data/reference/nitrogen_method_rya1220/evidence_matrix.csv) | Ten full-text source records, regime, diagnostics, physics, limitations, and page/table locators |
| [Source manifest](../../data/reference/nitrogen_method_rya1220/source_manifest.json) | Nine acquired PDFs and SHA-256 checksums |
| [Atomic reference rows](../../data/output/rya1220/audit_v1/ni_five_line_reference.csv) | Exact five-line reference, canonical join and source data |
| [Molecular reference rows](../../data/output/rya1220/audit_v1/molecular_reference_lines.csv) | Ninety N-bearing final-paper rows, with identity limitations retained |
| [Prior line audit](../../data/output/rya1220/audit_v1/ni_prior_line_audit.csv) | Twenty-four prior rows: four lines, three holdings, two method labels |
| [Prior budget audit](../../data/output/rya1220/audit_v1/prior_budget_adjudication.csv) | Term-by-term evidence gaps; no inferred zero errors |
| [Observed holding comparison](../../data/output/rya1220/observed_v1/) | Original reader flux, feature integrals, coverage checks, profiles |
| [Measured response table](../../data/reference/nitrogen_method_rya1220/measured_sensitivities.json) | Symmetric numerical probes at 8216 Å; steps are not adopted uncertainties |
| [Draft policy](../../data/reference/nitrogen_method_rya1220/policy.json) | Proposed routes and outstanding admission conditions; not ratified or publisher-enforced |
| [Contract migration](../../data/results/rya1220/contract_migration_v2/summary.csv) | Twelve historical N records mapped to the merged RYA-587 schema; all remain HOLD |
| [CN response runs](../../data/output/rya1220/cn_iag_coupling_v2/fits.json) | Same-pool C/O perturbations and stellar probes; measured responses retained, covariance unresolved |
| [CN RYA-587 candidate](../../data/output/rya1220/cn_contract_candidate_v1/candidate.json) | Measured profile/stellar components plus explicit HOLD components; publication gate refuses it |

Source text and run products are working evidence under `data/output/rya1220`; they are not solar gold. Copyrighted full text is a local research aid, not a proposed public repository payload.

## 1. Atomic solar reference: what must actually be reproduced

Amarsi et al. use **disk-centre intensity**, five selected N I lines, and empirical CN subtraction for four of them. They compare 3D, mean-3D and 1D atmospheres with both LTE and non-LTE. BALDER calculations use a STAGGER 3D solar atmosphere, modern electron collisions and hydrogen-collision treatments; these physics must be recorded independently of an abundance label. Their adopted oscillator strengths are theoretical MCHF values, not solar calibrations. [Amarsi et al. 2020, §§2–3, Tables 1–3](https://doi.org/10.1051/0004-6361/202037890).

| Air wavelength Å | Lower EP eV | log gf | Total EW pm | CN blend pm | N I EW pm | 1D LTE A(N) | 3D NLTE A(N) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 7442.29 | 10.330 | −0.403 | 0.310 | 0.075 | 0.235 | 7.806 | 7.767 |
| 8216.33 | 10.336 | +0.138 | 0.770 | 0 | 0.770 | 7.806 | 7.767 |
| 8629.23 | 10.690 | +0.077 | 0.620 | 0.210 | 0.410 | 7.812 | 7.759 |
| 8683.40 | 10.330 | +0.106 | 0.865 | 0.115 | 0.750 | 7.805 | 7.768 |
| 10108.90 | 11.753 | +0.444 | 0.275 | 0.075 | 0.200 | 7.835 | 7.772 |

These are **published measurements**, not new Codex EWs. One pm equals ten mÅ. The respective five-line means are 7.813 (1D LTE), 7.804 (1D NLTE), and 7.767 (3D NLTE). The tiny line scatter does not eliminate the paper's approximately 0.05 dex systematic uncertainty. The 8216 Å wings still require attention despite its zero tabulated CN subtraction. [Amarsi et al. 2020, Tables 1–3](https://doi.org/10.1051/0004-6361/202037890).

### Exact identity and line-data status

The four canonical matches are unique under the audit's wavelength/EP join (±0.05 Å, ±0.005 eV); these tolerances are screening criteria, not substitutes for quantum-level identity. Their current canonical gf values differ from the paper by at most 0.00201 dex. That difference cannot plausibly explain a roughly 0.2 dex discrepancy on these lines. The 10108.90 Å transition remains absent from the canonical inventory; the isolated overlay is retained separately and does not alter canonical gf.

Li et al. provide newer MCDHF/RCI calculations. Their Table 7 Babushkin/Coulomb log gf pairs for these lines are −0.429/−0.433, +0.116/+0.108, +0.070/+0.068, +0.084/+0.079, and +0.410/+0.406. These are candidates for a source-quality review, not permission to silently replace the reproduction's input scale. Gauge agreement is one quality diagnostic; it is not a complete uncertainty certificate. [Li et al. 2023, Tables 6–7](https://doi.org/10.3847/1538-4365/acb705).

The same primary paper compares laboratory measurements from Musielok et al. (1995) and Bridges & Wiese (2010). Consequently a categorical claim that no experimental N I transition probabilities exist is unjustified. The primary records are now in the bibliography: [Musielok et al. 1995](https://doi.org/10.1103/PhysRevA.51.3588) reports roughly 11–15% expanded uncertainties, while [Bridges & Wiese 2010](https://doi.org/10.1103/PhysRevA.82.024502) reports strong/weak-line agreement differences. Their line tables still need direct acquisition before any experimental values are admitted to canonical data. Theory wavelengths also require level identity checking before a merge. No canonical update is made here.

### Independent solar flux analysis

Mashonkina & Ryabchikova use the same five lines in a solar **flux** atlas, MARCS, DETAIL/SynthV_NLTE and a different model-atom implementation. They report 1D LTE 7.94 and 1D NLTE 7.92, with a hybrid 3D estimate 7.88. Their continuum treatment changes 8216 Å by approximately −0.12 dex and 8683 Å by −0.05 dex. This is direct evidence that geometry, continuum and blend treatment must be reproduced before interpreting an abundance difference as failed N physics. Their hybrid estimate is not a self-consistent 3D NLTE synthesis. [Mashonkina & Ryabchikova 2024, solar analysis and model-atom sections](https://arxiv.org/abs/2406.11367).

## 2. New observed tests on Sirius

Runs used the existing CNO engine in `/home/damienabraxas/scratch/rya1214`, with output under `/mnt/codex-data/outputs/rya1220`. Engine and line-list hashes are in each completed run's `provenance.json`. The remote engine checkout differed from this local audit branch; the source hashes, not an assumption of identical checkouts, establish run identity.

The diagnostic uses 1D LTE Turbospectrum flux, source-normalized IAG observations, fixed broadening, canonical gf, and full available molecular synthesis. Solar inputs were 5772 K, log g=4.438, [Fe/H]=0, ξ=1 km/s; background C=8.46 and O=8.69. The synthetic interval is ±6 Å; only the declared central mask contributes to fitting. The 0.01 flux scale in the objective is a model-adequacy scale, **not measured photon noise**. No formal abundance uncertainty follows from its reduced chi-square.

The source-normalization question is separate from residual-placement uncertainty. IAG is delivered normalized, so the continuum experiment did not renormalize the atlas. It divided the observed flux by `(1+f)` for `f=±0.001` and refit the same CN pool: nominal A(N)=7.716 moved to 7.898 and 7.916, a 0.200 dex span. This is measured diagnostic sensitivity, not an adopted error. The same response remains to be run for HARPS and CRIRES+.

The IAG CN pool was then rerun on Sirius with the accepted 3D C/O anchors (A(C)=8.46, A(O)=8.73). The exact 2,211-pixel pool gives A(N)=7.941, reduced χ²=42.7, with dA(N)/dC=−1.045 and dA(N)/dO=+0.150; all five fits were constrained. This supersedes the earlier raw-1D-LTE-input diagnostic for coupling work. The run remains diagnostic until C/O covariance and the remaining IR terms are attached.

The priority routes are already executed on Sirius with the vendored CN molecular list and Turbospectrum: HARPS CN-red gives A(N)=7.384 (σfit 0.835, reduced χ² 33.5) and CRIRES+ J gives A(N)=8.016 (σfit 0.096, reduced χ² 56.5). Both remain HOLD because fit curvature alone does not close the shared contract. Values and the unresolved list are in [`priority_routes_status.json`](../../data/output/rya1220/priority_routes_status.json).

Near-UV is tracked separately. Its selected OH A-X windows are 3063–3066 and 3122–3125 Å; the existing Kitt Peak result A(O)=9.745 is explicitly a continuum/blend-limited upper bound, not a C/O anchor. The rejected wider OH windows are CH-dominated. Near-UV C is pinned from the optical because the arm has no independent carbon primary. Thus near-UV contributes a diagnostic failure/upper-bound constraint and its own continuum uncertainty, while the accepted C/O anchor remains the 3D-calibrated optical atomic set.

| N I wavelength Å | A(N), ±0.25 Å | A(N), ±1.1 Å | Wide minus narrow dex | Status |
|---:|---:|---:|---:|---|
| 7442.29 | 8.1341 | 8.0934 | −0.0407 | diagnostic |
| 8216.336 | 7.9562 | 8.1722 | +0.2160 | diagnostic; surrounding-profile sensitivity |
| 8629.23 | 8.0744 | 8.0937 | +0.0193 | diagnostic |
| 8683.40 | 8.0655 | 8.1158 | +0.0503 | diagnostic; extended Ca/H wings still owed |
| 10108.90 | 8.4351 | 8.7762 | +0.3411 | isolated Amarsi+Brooke diagnostic; non-canonical overlay |

The ±6 Å synthesis interval does not adequately establish the extended Ca II 8662/H 8665 wing background near 8683 Å. Fixed broadening, unmeasured local continuum and incomplete blend validation are limitations across the experiment. These values therefore cannot be aggregated into an admitted solar abundance or called an exact Amarsi reproduction.

A synthetic smoke test with input A(N)=8.100 recovered 8.0951. This checks numerical inversion using the same forward model. It does not test atomic-data accuracy, observing geometry, atmosphere adequacy or deblending. The driver now also rejects failed, nonfinite and boundary optimizer results; those added checks were not present in the archived first runs.

At 8216 Å, symmetric probes were completed for continuum normalization ±0.001, C and O ±0.1 dex, and ξ ±0.1 km/s. See the response table for both sides, central derivatives and curvature. Responses below the optimizer's 0.005 dex tolerance should be treated as unresolved until step-size/convergence tests are done. None of these probe sizes is an adopted parameter uncertainty.

### Shared-contract CN response run

After RYA-587 merged, the IAG CN A–X (0–0) diagnostic was rerun on its exact 12-window, 2,211-pixel pool with the established CNO engine. The nominal fit is A(N)=7.716, with profile curvature σ=0.016 dex and reduced χ²=43.1. The fit is numerically constrained, but the high reduced χ² means the profile scale and residual correlation require adjudication before the curvature value can be a publication measurement.

| Perturbation | Refitted A(N) | Central response |
|---|---:|---:|
| C −0.100 dex | 8.004 | −0.975 dex per dex |
| C +0.100 dex | 7.809 | |
| O −0.100 dex | 7.896 | +0.150 dex per dex |
| O +0.100 dex | 7.926 | |
| Teff −/+1 K | 7.716 / 7.716 | 0 at 0.001-dex fit precision |
| ξ −/+0.2912 km/s | 7.716 / 7.907 | +0.654 dex per km/s central slope; visibly asymmetric |

The C/O responses are measured on the same physical window pool, but the adopted C/O abundance covariance is not yet available. The molecular-coupling component therefore remains `HOLD`, rather than being converted to a quadrature total. The ξ result uses the shared solar allowance from `uncertainty_stack.params_and_deltas`; the earlier ±0.1 km/s probe was a numerical test and is not substituted for it. Full machine-readable provenance and both sides are in the run directories. The stellar run also recorded Teff and ξ; solar log g and [Fe/H] are definitionally zero in the shared parameter record.

### Item 1: C/O covariance — what is required

The missing object is a 2×2 covariance matrix for `[A(C), A(O)]` in dex², measured from the same star, route and conditioning used by the CN fit. Separate C and O error bars are insufficient because molecular equilibrium couples them. We also need exact-pool CN refits at covariance-supported perturbations and a positive-semidefinite check before attaching `molecular_coupling` to RYA-587. The current VIS inputs are A(C)=8.488, σ=0.034 and A(O)=8.810, σ=0.421; the C/O product explicitly says `NOT_PROPAGATED`, and the O fit has reduced χ²=66.97. Using the measured CN derivatives (−0.975, +0.150), an independent-input diagnostic is 0.0713 dex, with a full correlation range of 0.0300–0.0963 dex. Those numbers are not adopted: they expose the scale of the missing term while the covariance and the poor oxygen constraint remain unresolved. See [`c_o_covariance_assessment.json`](../../data/output/rya1220/c_o_covariance_assessment.json).

This is a measured sensitivity result, not a new solar N product. It also explains why the old CN total of 0.050 dex cannot be carried forward: it omitted the measured C/O coupling, used an incomplete parameter treatment, and did not establish profile-pixel covariance or holding conditioning.

The resulting CN candidate has measured profile, Teff and ξ components, definitionally exact solar log g and [Fe/H], and eleven material HOLD components. Calling the shared validator returns one publication problem: unresolved components. This is the intended RYA-587 behavior; the candidate is review evidence and does not enter the live N feed.

Observed feature integrals were independently computed across holdings with their original normalization. They include blends and are labelled `BLENDED_FEATURE_NOT_ATOMIC_EW`. They must not be compared directly with empirically deblended disk-centre N I EWs. The 10108 Å holding profiles show differing continua/absorption; their comparison remains diagnostic until normalization, telluric and blend terms are separated.

## 3. Molecular solar reference and the CN comparison

The final Amarsi et al. molecular paper has 408 lines in 12 groups and iterative CNO inference. Its nitrogen subset is **59 CN A–X (0–0), 13 NH pure rotational, and 18 NH fundamental rows**. The NH diagnostics retained in this analysis are infrared, not the common NH feature at 3360 Å. [Amarsi et al. 2021, Tables 1–3 and molecular analysis](https://doi.org/10.1051/0004-6361/202141384).

| Group | Rows | Vacuum wavelength nm | Mean 3D A(N) | Mean MARCS A(N) | Mean 3D−MARCS dex |
|---|---:|---|---:|---:|---:|
| CN A–X (0–0) | 59 | 1087.516–1320.778 | 7.8668 | 7.9269 | −0.0601 |
| NH X–X Δv=0 | 13 | 11311.232–15023.045 | 7.8809 | 8.0215 | −0.1405 |
| NH X–X Δv=1 | 18 | 2891.001–3445.428 | 7.9214 | 7.9144 | +0.0071 |

These are descriptive means of the published rows, not a fresh abundance solution. The variation in correction, including its sign, rules out one universal molecular correction. Corrections must be joined to the actual transition, atmosphere, observing geometry and measured line subset. The extracted table lacks full rotational identities: a wavelength match alone does not ratify a molecular transition join.

AGSS21 describes a different nitrogen inventory than the final 2021 molecular paper, including different NH counts and additional CN sets. Keep those inventories separate; the overview's atomic/molecular discrepancy and adopted solar value do not authorize mixing line memberships. [Asplund, Amarsi & Grevesse 2021, §3.2, nitrogen discussion](https://doi.org/10.1051/0004-6361/202140445).

### What the existing holdings can test

The current IAG holding is declared through 11083.46 Å; the sampled edge probe reaches 11082.9878 Å. The current CRIRES+ J rest-frame product begins at 11159.9497 Å. The intervals do not overlap. A matched IAG/J CN comparison therefore has **zero possible lines**, irrespective of abundance agreement. Use a genuinely overlapping holding pair for repeatability, and report IAG/J only as separate subsets after consistent C/O inputs and line-specific corrections.

The previous IAG CN analysis used C=8.488 whereas the inspected J product used C=8.46. Carbon is not an innocuous normalization for CN. CO equilibrium changes the available carbon, and oxygen therefore also affects fitted nitrogen. The present atomic C/O probes do not substitute for molecular C/O measurements. Symmetric molecular refits with a measured C/O covariance remain required.

## 4. Regime-specific proposed routes

These are evidence-supported candidate routes for ratification. Literature sample limits are not validated universal pipeline cutoffs. Each new target still requires coverage, detection, line identity, model-domain and uncertainty checks.

| Regime | Candidate nitrogen route | Required qualifications |
|---|---|---|
| Sun / near solar twins | The five N I lines plus separate CN/NH molecular checks | Match geometry and blends; preserve atomic/molecular results separately; no forced common zero point |
| FGK dwarfs with near-UV coverage | NH around 3345–3375 Å; CN where measurable; N I if demonstrably detected | UV continuum/blends, temperature sensitivity, independent line data, applicable atmosphere and formation treatment |
| Solar twins with high-quality optical spectra | Blue CN features near 4180–4212 Å as a differential route | Joint C/O constraints, exact paired masks; do not import solar-calibrated gf |
| Cool giants | H-band CN alongside CO and OH | Iterative CNO equilibrium, isotope mixture, correlated errors, stellar-evolution interpretation |
| Metal-poor stars | NH 3360 Å and CN 3883/3888 Å where detected | Explicit line-list comparison and 3D sensitivity; nondetections remain limits; no blanket NH-to-CN offset |
| Warmer A/F stars | N I with dedicated non-LTE model atom | Verify radiation field, collision data, rotational blending, parameter coverage; solar sparse corrections cannot be extrapolated |

Suárez-Andrés et al. analyze UV NH in dwarfs (4583–6431 K before their cool-star exclusion), with R≈80,000 and typically S/N≈150. Their solar-calibrated line list cannot be imported into our canonical gf inventory; the study supplies diagnostic and failure-mode evidence. [Suárez-Andrés et al. 2016, analysis and uncertainty sections](https://arxiv.org/abs/1605.03049).

Botelho et al. use HARPS solar twins and CN features at 4180.02, 4192.94, 4193.40, 4195.95 and 4212.25 Å. The exact feature masks are retained in the evidence matrix/source text. Their use of solar-calibrated gf again means that methodological success is not an independent laboratory-data endorsement. [Botelho et al. 2020, spectral analysis](https://arxiv.org/abs/2009.09003).

Smith et al. demonstrate H-band synthesis in five cool giants, with iterative CO/C, OH/O and CN/N determination. Their observed spectra span roughly R=45,000–100,000. This supports coupled molecular inference, not fixed solar C/O in evolved stars. [Smith et al. 2013, CNO analysis](https://arxiv.org/abs/1212.4091).

Spite et al. document NH/CN discrepancies in extremely metal-poor giants. The later analysis shows that changing NH line lists alone can shift inferred abundances by approximately 0.50 dex over 3357–3365 Å, with feature-dependent differences. An empirical universal −0.4 dex correction would conflate transition data, atmosphere and diagnostic effects. [Spite et al. 2005](https://arxiv.org/abs/astro-ph/0409536), [Spite et al. 2022, Table 2 and NH line-list comparison](https://arxiv.org/abs/2209.10219).

The inspected hotter-star study includes an A/F sample spanning approximately 7250–10400 K. This is sample coverage, not proof of a general grid extending across all gravities or metallicities. Its model-atom/collision choices must be assessed for each application. [Mashonkina & Ryabchikova 2024, Table 2](https://arxiv.org/abs/2406.11367).

## 5. Uncertainty and differential contract

The implementation in [`pipeline/nitrogen_uncertainty.py`](../../pipeline/nitrogen_uncertainty.py) propagates measured line responses and validates finite, symmetric, positive-semidefinite covariance matrices. It refuses missing or changing transition membership during perturbations. Required evidence terms are statistical, continuum, line data, stellar parameters, blend abundances, C/O coupling, tellurics, atmosphere, line formation, holding repeatability and diagnostic family.

For line abundances **a** and parameters **p**, record symmetric responses

`J[i,k] = (a_i(p_k+h_k) − a_i(p_k−h_k)) / (2 h_k)`.

Retain the two one-sided slopes and curvature. Refit the same lines at both steps; missing/failed sides remain failures. Select numerical steps through convergence checks, independently of observational parameter errors. The absolute covariance is `Σa = J Σp Jᵀ + Σresidual`. Correlated C/O errors are entries in `Σp`; shared molecular line data and continuum can produce off-diagonal line covariance in the appropriate block. Do not count the same nuisance source in both blocks.

For a paired target-minus-Sun vector, use the joint parameter matrix with both cross-star blocks and `Jpair = [Jtarget, −Jsun]`. Shared gf errors cancel only to the extent supported by the measured responses and justified cross-star covariance. Form each transition's difference **before** aggregating. The helper deliberately uses declared fixed weights: it does not imply that a median, clipping rule or data-dependent weighting has this same uncertainty formula.

The shared RYA-587 suite and N adapter tests pass (37 tests across both suites). This is a contract and mathematical implementation check, not a validated observational uncertainty budget.

### Existing products re-adjudicated

The prior near-UV, IAG CN, KP CN and visible CNO budgets were inspected. For example, their reported N totals include 0.067, 0.050 and 0.059 dex for the first three routes. Those numbers are not complete errors under this gate. The previous one-sided stellar probes and absent covariance do not establish the full requested stack. The term audit distinguishes absent evidence from zero uncertainty; statistical terms are also only reported partial estimates.

The local derived N non-LTE correction table contains 29 rows at three wavelengths (7468.31, 8216.34, 8683.40 Å), over 5100–6200 K, log g 4.0–4.7, and [Fe/H] −0.3 to +0.6. It is neither the full model atom nor complete coverage of the five reference lines. A full Amarsi/GALAH departure binary is held on Sirius; its presence does not prove that the current product route uses it correctly. Distinguish data availability, executable applicability, and successful per-line evaluation. The old RYA-1214 budget text's statement that no primary laboratory N I table exists is superseded by the two primary-source leads above; its historical numeric budget remains unchanged and is not promoted.

## 6. Route uncertainty assembly

The completed Sirius probes are now assembled per instrument and band in [`data/output/rya1220/cno_route_uncertainty_summary.json`](../../data/output/rya1220/cno_route_uncertainty_summary.json). The record keeps the Fe-style local continuum perturbation separate from the other RYA-587 terms and does not turn a failed perturbation into a symmetric error.

The measured continuum responses are: Kitt Peak near-UV NH, 0.310 dex span; IAG CN, 0.200 dex in the earlier residual-continuum probe plus the 3D C/O response; CRIRES+ J CN, 0.029 dex; CRIRES+ H OH, 0.213 dex; and CRIRES+ K CO, 0.107 dex. HARPS VIS CN has a one-sided edge-pinned local-continuum fit (A(N)=6.638 on the negative side), so that term is held rather than symmetrized. The H-band OH fit has red χ²≈130.7 and remains diagnostic only.

CRIRES+ telluric evidence is inherited from the measured RYA-1191/1192 window products and is marked as joined or awaiting mask-level joining in the route record. This avoids rerunning a separate synthesis comparison while preserving the distinction between measured telluric residuals and a clean preflight window. No route is publication-ready: C/O covariance, profile residuals, blend/isotope identity, molecular 3D/model terms, and holding repeatability remain explicit HOLD components where they have not been measured.

The Kitt Peak NIR CN probe is now complete. The existing Brooke CN IR list was linked into the Sirius iSpec molecule directory and the five-fit continuum sweep completed with all fits constrained: nominal A(N)=8.198, local-level span 0.235 dex (8.076–8.311), slope response 0.002 dex, and red χ²=202.1. This is measured route evidence, but the large residual means profile/blend/model terms remain HOLD; it is not a publication value.

That VIS/UV audit is now recorded in [`data/output/rya1220/vis_uv_route_audit.json`](../../data/output/rya1220/vis_uv_route_audit.json). ESPRESSO O I 777/[O I] 6300 agrees at 0.01 dex and remains useful route evidence. ESPRESSO carbon has a 0.112 dex line spread, while UVES N I 8216 versus CN red differs by 0.528 dex; both are flagged disagreements, not averaged. The HST UV loader has a smoke-tested conditioned arm, but production fitting remains gated on the FUV C/O 3D-NLTE grid and pseudo-continuum treatment.

The HST gate is separately captured in [`data/output/rya1220/hst_uv_readiness_gate.json`](../../data/output/rya1220/hst_uv_readiness_gate.json). The whitelist contains 117 target-confirmed science rows spanning E140M, E140H and E230H, and the loader/conditioning path records the correct vacuum-to-air handling and chromospheric masking. The route is data-ready, not abundance-ready: no CNO fit is admitted until the Amarsi C/O grid, FUV pseudo-continuum, and matched target–Sun covariance gates are available.

The complete component-by-route state is now consolidated in [`data/output/rya1220/cno_component_gate_matrix.json`](../../data/output/rya1220/cno_component_gate_matrix.json). Every route has an explicit state for profile measurement, continuum, C/O coupling, tellurics, blends/isotopes, 3D/NLTE/model treatment, holding repeatability, and target–Sun covariance. This closes the accounting step: remaining work is now concrete execution or data admission, rather than an unspecified “uncertainty expansion.”

The blend and molecular identity audit is now captured in [`data/output/rya1220/cno_blend_identity_audit.json`](../../data/output/rya1220/cno_blend_identity_audit.json). It reconciles 328 molecular joins: 278 exact tuple matches, with 50 ambiguous, unresolved, strength-mismatched, or unmatched rows retained as HOLD. It also records the dominant atomic blends that limit C and N lines. These identity results are now part of the route admission state; they are not folded into a numerical error until the transition joins and isotope treatment are resolved.

The model-form audit is recorded in [`data/output/rya1220/cno_model_form_audit.json`](../../data/output/rya1220/cno_model_form_audit.json). Existing same-line CN comparisons support diagnostic 1D-to-3D shifts of −0.081 dex (Kitt Peak) and −0.082 dex (IAG), while the near-UV molecular sizing gives a −0.050 dex missing-opacity lever as a lower-bound diagnostic. These values are not universal corrections. Atomic C/O terms remain admitted where the cited 3D/NLTE grids exist; CN, NH, OH, and CO model terms stay HOLD until transition-specific calculations are available.

Holding repeatability is now measured for the shared atomic pools in [`data/output/rya1220/cno_holding_repeatability_audit.json`](../../data/output/rya1220/cno_holding_repeatability_audit.json): C VIS spans 0.021 dex, N red-optical 0.064 dex, O VIS 0.105 dex, and O red-optical 0.026 dex across the registered solar holdings. Molecular routes and target–Sun repeatability remain open because they do not have matched pools.

The target–Sun pair audit is now recorded in [`data/output/rya1220/cno_target_sun_pair_audit.json`](../../data/output/rya1220/cno_target_sun_pair_audit.json). Procyon and Sun share six HARPS VIS transition keys, so the pair is structurally available. The raw paired differences range from −0.210 to +2.404 dex, but the target/solar parameter Jacobians, cross-star covariance, and profile residual covariance are not present. These values remain diagnostic red flags rather than differential abundances.

For the solar IAG CN route, the C/O covariance propagation is now explicit in [`data/output/rya1220/solar_c_o_covariance_propagation.json`](../../data/output/rya1220/solar_c_o_covariance_propagation.json). Using the accepted atomic-anchor uncertainties (σC=0.038 dex, σO=0.050 dex) with documented independent-source covariance gives a PSD matrix and σN(C/O)=0.0405 dex through the measured CN Jacobian. This is a measured partial term; the same-pool joint C/O-plus-CN fit remains the next solar execution task.

That solar same-holding execution is complete in [`data/output/rya1220/solar_iag_vis_cno_v2`](../../data/output/rya1220/solar_iag_vis_cno_v2). It measured C I 5052=8.518, C I 5380=8.426, C2=8.497, CN red=8.550, and [O I] 6300=8.907. The oxygen fit is poor (red χ²=64.431), so these values provide route context only; the accepted O=8.73 anchor remains the C/O input for the CN product.

The solar profile-quality gate is now consolidated in [`data/output/rya1220/solar_cno_profile_quality_audit.json`](../../data/output/rya1220/solar_cno_profile_quality_audit.json). All completed molecular routes have elevated residuals, with Kitt Peak NIR CN at χ²ᵣ=202.1 and CRIRES+ H OH at χ²ᵣ=130.7 the clearest model failures. These residuals are retained as profile/model HOLDs rather than inflated into arbitrary uncertainties.

The first residual decomposition is complete for IAG CN in [`data/output/rya1220/iag_cn_profile_residual_v2/profile_residuals.json`](../../data/output/rya1220/iag_cn_profile_residual_v2/profile_residuals.json): 2,211 pixels, weighted residual RMS 0.0653 in normalized flux, and mean lag-1 residual correlation 0.473. This is a measured correlated-residual component. It does not close the route because the correlation is evidence of remaining blend/model structure, so the profile term remains measured-partial and the route stays HOLD.

The same decomposition is now complete for Kitt Peak NIR CN in [`data/output/rya1220/kp_cn_profile_residual_v1/profile_residuals.json`](../../data/output/rya1220/kp_cn_profile_residual_v1/profile_residuals.json): 3,215 pixels, weighted residual RMS 0.1422, and mean lag-1 correlation 0.725. The high correlation confirms a structured model/blend failure rather than independent pixel noise; the route remains HOLD.

Blend handling has been adjudicated separately from EW viability. The CNO routes use full Turbospectrum synthesis with the line lists and molecular opacity included, so a weak target line or an EW failure does not by itself block a synthesis route. The remaining gates are transition/isotope provenance and the measured structured residual/model terms; the component matrix now records blends as modeled-by-synthesis rather than treating blend contribution itself as an automatic HOLD.

The line-priority review is in [`data/output/rya1220/solar_blend_priority_review.json`](../../data/output/rya1220/solar_blend_priority_review.json). It identifies the strongest solar C and O atomic anchors and the N lines where synthesis is essential because neighboring species dominate. The review changes no line selection or abundance; it records which lines should carry future profile/model validation.

### Full-synthesis nitrogen blend demonstration

The first direct blend test requested for this continuation is now complete on Sirius. The IAG CN A–X window at 10870.68–10874.00 Å was synthesized with Turbospectrum using the accepted solar C/O anchors (A(C)=8.46, A(O)=8.73) and the current IAG CN nitrogen input A(N)=7.941. The full 112-pixel synthesis has mean line depth 0.02133 and RMS residual 0.09843 against the observed normalized profile.

Two component counterfactuals were run on the identical pixels. Removing molecular opacity changes the mean depth to 0.01803 and increases RMS to 0.09963; the molecular blend contribution is therefore 0.00330 in normalized depth, about 15.5% of the full synthetic depth. This is direct evidence that the CN blend is being carried by the synthesis rather than discarded because an EW decomposition is inconvenient. The machine-readable result is [`nitrogen_blend_synthesis.json`](../../data/output/rya1220/n_blend_synthesis_v2/nitrogen_blend_synthesis.json).

The IAG CN atomic linelist contains zero N I transitions in this molecular window, so the no-N-I counterfactual is identical to the full result. That is expected for this route and does not test the five optical N I lines. The remaining limitation here is the large structured profile residual already measured for IAG CN; this blend demonstration closes the “are blends actually synthesized?” question, while model/profile and isotope provenance terms remain HOLD.

## 7. Gate checklist and next execution

| Requirement | Current state | Evidence still required |
|---|---|---|
| Primary literature and regime matrix | Initial review complete; laboratory originals outstanding | Inspect original experimental N gf sources before canonical selection; expand regimes if target demands it |
| Exact five-line solar reproduction | OPEN | Admit 10108 Å using independently justified identity/gf; geometry-matched intensity data or explicitly separate flux benchmark; exact blend treatment |
| Profile root-cause analysis | Partial measured evidence | Continuum fitting, extended H/Ca wings, species-specific blend contribution and mask/step convergence |
| Molecular CN/NH line audit | Reference inventory extracted; identities incomplete | Quantum-level joins, available holding masks, molecular symmetric C/O and stellar-parameter refits |
| Full uncertainty stack | RYA-587 integrated; IAG residual-continuum response and priority HARPS/CRIRES+ Turbospectrum diagnostics are recorded; budget still HOLD | C/O covariance, priority-route profile/continuum responses, CRIRES+ telluric covariance, blend/isotope admission, molecular 3D/non-LTE treatment, holding repeatability, and matched target–Sun differential run |
| Target–solar differential behavior | Contract and algebra tested | Same-transition, same-method observed target/Sun pilot with joint covariance |
| Regime policy | Draft | Scientific ratification and executable route/admission wiring |
| Prior RYA-1214 promotion | HELD | Complete the above; retain original products and diagnostic status |

The immediate blockers to an exact five-line run are canonical transition admission, complete CN/atomic blend coverage, and geometry/continuum validation for 10108 Å. Repair must follow source-quality and level-identity review rather than solar calibration. The isolated fit establishes useful diagnostic evidence, but the gate remains open. No gold file, canonical gf inventory, production abundance, publication route, or Linear completion status was changed by this audit.

## Reproducibility

Local working branch: `codex/rya-1220-nitrogen-method`, based on `3371e140`. Sirius runs used the existing CNO engine, with per-run source hashes. Original spectra and previous products were retained.

```sh
python -m pytest tests/test_nitrogen_uncertainty_rya1220.py -q
python scripts/rya1220_evidence_audit.py --out NEW_AUDIT_DIRECTORY --prior data/output/rya1220/prior_products
# Sirius only; use a fresh output directory and the configured iSpec environment:
python pipeline/diagnostics/nitrogen_reproduction.py --engine-root /home/damienabraxas/scratch/rya1214 --out NEW_RUN_DIRECTORY --line 8216.336
```

The reproduction driver consumes private engine APIs and requires the registered holdings and engine dependencies. It is a diagnostic harness, not a portable production command. Run directories are exclusive-create to preserve earlier evidence.

### Complete solar nitrogen synthesis sweep status

The Sirius Turbospectrum sweep now covers the four observed atomic N I lines in the Kitt Peak solar atlas: 7442.29, 8216.336, 8629.23 and 8683.40 Å. Narrow versus wing-inclusive profile fits are respectively 8.014/7.896, 8.153/8.286, 8.069/8.087 and 8.085/8.136 dex. These are line-level 1D-LTE diagnostics; the spread is retained as profile/continuum evidence, not averaged into a product.

The molecular nitrogen routes are also full Turbospectrum synthesis: near-UV NH, optical CN on HARPS and IAG, Kitt Peak and IAG NIR CN, and CRIRES+ J CN. The complete machine-readable inventory is [`nitrogen_full_synthesis_sweep.json`](../../data/output/rya1220/nitrogen_full_synthesis_sweep.json).

The Sirius iSpec installation exposes additional radiative-transfer backends (`spectrum`, `moog`, `moog-scat`, `synthe`, and `sme`). Probe runs for the atomic 7442.29 Å window did not return within the configured execution window; they are recorded as `PROBE_TIMEOUT`, rather than being presented as completed results. The `sme` package is absent on Sirius. None of these backends is substituted for Turbospectrum on CN/NH because their molecular opacity path is not wired into this product. Turbospectrum remains the complete blend-capable route; the backend probe timeout is an infrastructure follow-up, not a scientific equivalence claim.
