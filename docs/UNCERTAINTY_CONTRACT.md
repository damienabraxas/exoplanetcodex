# Universal publication uncertainty contract — RYA-587

Status: implemented for review; scientific budget population remains incomplete.
The 2026-09-15 RYA-587 directive governs every element, indicator, holding and model
route. `pipeline.uncertainty_contract` owns schema `rya587.uncertainty.v1`, the
component vocabulary, assembly, validation and covariance propagation. It reuses
`error_budget.Term` semantics. Legacy arithmetic remains diagnostic; historical
subtotals are not proof of publication completeness.

## Machine record and API

Each product supplies `star`, `uncertainty_indicator_ids` (stable physical IDs),
`uncertainty` and `sigma_reported`. The uncertainty document contains:

- Scope: star, canonical ten-axis product key and sorted physical-pool SHA256.
- Components: name, MEASURED/DEFINED/N/A/HOLD, abundance sigma in dex (null for
  N/A/HOLD), source and structured evidence.
- Full component covariance in dex², numeric-component order, source and explicit
  independence/unknown-covariance assumptions.
- Recomputed total only when all required components are resolved. A diagnostic
  `sigma_measured_partial` must never replace `sigma_reported`.

`assemble`/`validate` are the executable schema. Do not duplicate `COMPONENTS` in
an element-specific policy. `ErrorBudget.component_records()` exports existing
terms without inferring applicability for omitted terms. A legacy aggregate
cannot be reverse-engineered into line-specific evidence.

## Required evidence

| Family | Requirement |
| --- | --- |
| Measurement | Vetted independent-line raw scatter, N and SE; small-N flag. Single/correlated lines require an alternative. Profile likelihood records pixel-correlation treatment; never divide by pixel count. |
| Transition data | Per-indicator sigma, authoritative sources, response weights and covariance including shared laboratory scales. No blanket gf replacement. |
| Stellar | Individual Teff/logg/xi/metallicity applicability, sourced adopted uncertainties, exact-pool ±1σ responses, signed sensitivity and asymmetry assessment. |
| Observation | Continuum, EW/profile measurement, pseudo-continuum, telluric residual and holding/instrument assessment. IR telluric cannot be N/A. Record contributions already included in a likelihood to avoid double-counting. |
| Physics | NLTE, atmosphere, HFS/isotopes, blends and molecular coupling: measured/defined, justified N/A or HOLD. No invented term for a material unresolved effect. |
| Covariance | Symmetric positive-semidefinite matrix with component-consistent diagonal. Independence is an explicit sourced assumption. |

`paired_response` preserves both sides and rejects moved indicator pools. Its
linearized response does not ratify nonlinearity; the caller records an assessment.
A full-3D route without xi records N/A with `parameter_exists: false` and route
evidence. Teff/logg remain required. Mean-3D Gerber is experimental and is not a
required completeness route.

## Differential and ratios

`differential_uncertainty` requires absolute budgets recomputed on the matched
physical pool, matching element/ion/route/treatment/line-set axes and an explicit
target/reference cross-covariance block. It retains both absolute budgets. A shared
gf scale can cancel; independent absolute quadrature is never the default.

Target [X/H] records carry `differential_uncertainty`; the publication gate
recomputes it and checks `sigma_differential`. Do not manufacture a Solar
self-differential against a literature abundance.

`ratio_uncertainty` propagates joint abundance covariance to log differences or
number ratios, including ln(10) for C/O. [X/Fe] inputs must already be matched
reference differentials. Arithmetic alone is not a multi-element evidence package:
publication of [X/Fe]/C/O remains held pending that package.

## Wiring and migration

`product_eligibility.evaluate` invokes the same contract for every element, adding
`UNCERTAINTY_INCOMPLETE` alongside existing scientific holds. The common publisher
preserves JSON budgets, physical IDs and differential fields, stamps the star and
checks candidates. The final feed-write boundary and Fe schema stamper also enforce
the contract, including metadata-only updates. Budget-only changes count as updates. Al, Si, atomic C/N/O,
molecular N and future RYA-709/946 children inherit this gate automatically.
A validated profile-likelihood budget need not impersonate scatter SE.

Existing live feeds are not migrated or deployed here. Their old values and Fe
diagnostic row remain byte-identical. The standing live-feed eligibility test
exposes migration debt; it must not be waived or grandfather legacy records.
This branch is not merge-ready while that test is red.

## Reproduce

```bash
python3 scripts/audit_uncertainty_rya587.py \
  --output data/results/rya587/new-audit \
  --reference-repo /path/to/rya1213-checkout \
  --xi-runs /path/to/copied/rya1213_xi
python3 -m pytest tests/test_uncertainty_contract_rya587.py -q
```

The audit refuses an existing output directory, reuses the actual Reference
importer in an isolated process, records source SHA and file hashes and stages
recovered stamps in a separate HOLD review artifact. It never writes a live feed.
See [audit findings](audit/uncertainty_rya587.md).
