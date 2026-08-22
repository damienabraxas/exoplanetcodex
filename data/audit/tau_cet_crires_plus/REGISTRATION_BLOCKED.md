# tau Ceti CRIRES+ — staged and audited, deliberately NOT registered

The 4 frames are on Sirius at `/mnt/codex-data/spectra/tau_cet/CRIRESPlus/` and fully
audited in `tau_cet_crires_plus_manifest.csv` (target confirmed astrometrically, telluric
state measured three ways). They are **not** in `data/catalog/holdings_manifest_registry.csv`,
and that is a refusal, not an oversight.

## Why

`instruments.validate_holdings` resolves a holding's `system_id` against the set of
**non-blank `star_params_key`** values in `data/catalog/system_catalog.csv`. tau Ceti is in
that catalogue (row `tau Ceti`, slug `tau-ceti`) as a **`future_target` with a blank
`star_params_key`**. And `pipeline/system_catalog.py` requires that any non-blank key
**resolve in `STAR_PARAMS`** — so filling it in means adding real Teff / logg / [Fe/H] /
vmic for tau Ceti to `config/stars.yaml`.

That is a science decision needing cited sources, and this ticket is data prep only ("no
abundances"). Inventing parameters to satisfy a foreign key is the shape of defect this
project exists to avoid.

## 🔴 The schema gap this exposes

**The holdings registry cannot express "we hold data for a star we have no parameters
for."** Data acquisition and parameter adoption are independent events, and the registry
assumes the second precedes the first. Every future target we acquire spectra for hits this,
not just tau Ceti — `eps Eri` and `tau Boo` have CRIRES+ data on the drive today and blank
`star_params_key` values too.

## To register it

1. Adopt cited stellar parameters for tau Ceti into `config/stars.yaml`.
2. Set `star_params_key = tau_cet` on the `tau Ceti` row of `system_catalog.csv` and promote
   its role from `future_target`.
3. Append the holding row; its prepared text is in this ticket's Linear comment.

Alternatively, give the registry a way to hold a **parameterless** system, which is the
change that fixes the class rather than this instance.
