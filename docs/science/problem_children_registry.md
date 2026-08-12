# RYA-463 — Master problem-children line registry

_Generated 2026-08-11 by pipeline/problem_children.py. A CATALOG (never mutates a measured value); the cross-star feed for every per-star pre-run audit (RYA-205) and the line-level companion to the RYA-277 acceptance profiles._

**Rows:** 29 (29 curated, 0 auto-aggregated from RYA-458 EW-integrity).

## Schema

`species | lambda_or_scope | problem_class | required_treatment | observed_in | amplifies_with | severity | governing_tickets | status | population_source | notes`

- **problem_class** (curated): ATOMIC_BLEND, BAD_GF, CITED_NOT_MEASURED, CONTINUUM_LIMITED, DATA_GAP, HFS_SUMMING, MOLECULAR_BLEND, MOLECULAR_SYNTH_ONLY, NLTE_OWED, NLTE_VOID, SATURATION_COG, UPPER_LIMIT
- **required_treatment**: 3D, HFS_sum, NLTE_grid, astrophysical_gf_differential, cited_substitution, deblend, exclude, none, per_region_source, synthesis, upper_limit
- **amplifies_with**: the Teff / [Fe/H] axis where it worsens — the prediction key.

## Registry (curated layer)

| species | lambda_or_scope | problem_class | required_treatment | observed_in | amplifies_with | severity | governing_tickets | status | population_source |
|---|---|---|---|---|---|---|---|---|---|
| O I | [O I] 6300 | CONTINUUM_LIMITED | cited_substitution | Sun | [Fe/H]↑↑ | high | 447-455,460 | active | curated |
| C I | 5380.34 | SATURATION_COG | exclude | Sun | Teff↑ | high | 454,458 | active | curated |
| Li I | 6707.84 | MOLECULAR_BLEND | upper_limit | Sun | Teff↓/[Fe/H]↑ | medium | 103,458 | active | curated |
| Eu II | 6645.13 | HFS_SUMMING | HFS_sum | Sun | weak-line | medium | 102,458 | active | curated |
| NH | ~3360 (band head) | CONTINUUM_LIMITED | synthesis | Sun | blue/UV | high | 451,460 | active | curated |
| CN | 3883 (violet B-X) | MOLECULAR_SYNTH_ONLY | synthesis | Sun | Teff↓ | medium | 369,460 | active | curated |
| N I | 7442-8718 (red multiplets) | DATA_GAP | per_region_source | Sun | Teff↑ | high | 369,459,460,491 | owed | curated |
| Cr I | gf pool | BAD_GF | astrophysical_gf_differential | Sun | all | high | 398,399,161 | active | curated |
| Ti I | gf pool | BAD_GF | astrophysical_gf_differential | Sun | all | high | 398,161 | active | curated |
| Ni I | gf pool | BAD_GF | astrophysical_gf_differential | Sun | all | high | 398,161 | active | curated |
| Si I | gf pool | BAD_GF | astrophysical_gf_differential | Sun | all | high | 398,399,161 | active | curated |
| V I | all | NLTE_VOID | none | Sun | all | medium | 404 | active | curated |
| Fe I | Procyon gf outliers | BAD_GF | astrophysical_gf_differential | Procyon(pending) | Teff↑ | medium | 281 | predicted | curated |
| CH/CN/C2/NH/OH/CO | bands | MOLECULAR_SYNTH_ONLY | synthesis | all | Teff↓ | high | methodology | active | curated |
| P I | FUV / near-IR 10581/10596 | DATA_GAP | per_region_source | Sun | — | medium | 119,460 | active | curated |
| K I | 7665/7699 | NLTE_OWED | NLTE_grid | Sun | all | low | 460,462 | resolved | curated |
| Sc II | 4246.82 | HFS_SUMMING | HFS_sum | Sun | blue/UV | medium | 460 | active | curated |
| Co I | 3845.46 (demoted) / red HFS 4813/5212/5352/5483/6771 | CONTINUUM_LIMITED | synthesis | Sun | blue/UV | medium | 460,534,564 | resolved | curated |
| Sr II | 4077/4161/4215/4305 | SATURATION_COG | synthesis | Sun | [Fe/H]↑ | medium | 421,428,429,430,433,520,551 | resolved | curated |
| Mn I | 6013/16/21 (e6S->z6P triplet) | HFS_SUMMING | synthesis | Sun | all | medium | 411,468,473,476 | resolved | curated |
| Cu I | HFS-split lines | HFS_SUMMING | synthesis | Sun | all | medium | 402,466 | active | curated |
| Ba II | 5853 (+HFS, strong) | SATURATION_COG | synthesis | Sun | [Fe/H]↑ | medium | 279,165 | owed | curated |
| Zr II | 4629/5350/5372 (strong) | SATURATION_COG | synthesis | Sun | all | low | 279,458 | owed | curated |
| Y II | dominant ion absent from pool | DATA_GAP | per_region_source | Sun | blue/UV | medium | 458 | active | curated |
| Mg I | b-triplet 5167/72/83 + 5528 (strong) | SATURATION_COG | synthesis | Sun | all | medium | 410,235,534,561,592 | resolved | curated |
| Fe I | 8024.543 (IR) | ATOMIC_BLEND | exclude | Sun | all | critical | 782,760,523 | active | curated |
| Fe I | 7052.715 / 7120.021 (IR) | ATOMIC_BLEND | exclude | Sun | all | high | 782,523 | active | curated |
| Fe I | 7810.814 (IR) | BAD_GF | exclude | Sun | all | high | 782,523 | active | curated |
| Fe I | 7107.459 (IR) | BAD_GF | exclude | Sun | all | medium | 782,712,523 | owed | curated |

## Auto-aggregated layer (RYA-458 EW-integrity)

| species | lambda_or_scope | problem_class | observed_in | notes |
|---|---|---|---|---|

## Prediction layer — sample heads-ups

```
Problem-children heads-up for (5777 K, [Fe/H]=+0.00, Sun):
  - [O I] 6300 OK at this [Fe/H] but WATCH it for metal-rich targets (amplifies [Fe/H]↑↑) -> O I 777 is the safe primary.
  -> 29 entries: 0 amplified, 24 always-on, 0 watch. Top:
     [critical|expected ] Fe I 8024.543 (IR) -> ATOMIC_BLEND (exclude)
     [high  |expected ] O I [O I] 6300 -> CONTINUUM_LIMITED (cited_substitution)
     [high  |expected ] C I 5380.34 -> SATURATION_COG (exclude)
     [high  |expected ] NH ~3360 (band head) -> CONTINUUM_LIMITED (synthesis)
     [high  |expected ] N I 7442-8718 (red multiplets) -> DATA_GAP (per_region_source)
     [high  |expected ] Cr I gf pool -> BAD_GF (astrophysical_gf_differential)
```

```
Problem-children heads-up for (6554 K, [Fe/H]=+0.01, Procyon):
  - F-star (Teff↑): COG saturation knee shifts up to ~114 mA (est); watch strong-line saturation + bad-gf Fe outliers (RYA-281).
  - [O I] 6300 OK at this [Fe/H] but WATCH it for metal-rich targets (amplifies [Fe/H]↑↑) -> O I 777 is the safe primary.
  -> 28 entries: 14 amplified, 4 always-on, 5 watch. Top:
     [critical|amplified] Fe I 8024.543 (IR) -> ATOMIC_BLEND (exclude)
     [high  |amplified] C I 5380.34 -> SATURATION_COG (exclude)
     [high  |amplified] N I 7442-8718 (red multiplets) -> DATA_GAP (per_region_source)
     [high  |amplified] Cr I gf pool -> BAD_GF (astrophysical_gf_differential)
     [high  |amplified] Ti I gf pool -> BAD_GF (astrophysical_gf_differential)
     [high  |amplified] Ni I gf pool -> BAD_GF (astrophysical_gf_differential)
```

```
Problem-children heads-up for (5196 K, [Fe/H]=+0.32, 55 Cnc):
  - cool (Teff↓): molecular bands (CH/CN/C2/NH/OH/CO) strengthen — synthesis-only handling amplified.
  - very metal-rich ([Fe/H]↑↑): [O I] 6300 unreliable -> O I 777 primary; blends + saturation amplified.
  - [O I] 6300 UNRELIABLE at this metallicity -> use O I 777 (primary).
  -> 26 entries: 16 amplified, 3 always-on, 2 watch. Top:
     [critical|amplified] Fe I 8024.543 (IR) -> ATOMIC_BLEND (exclude)
     [high  |amplified] O I [O I] 6300 -> CONTINUUM_LIMITED (cited_substitution)
     [high  |amplified] Cr I gf pool -> BAD_GF (astrophysical_gf_differential)
     [high  |amplified] Ti I gf pool -> BAD_GF (astrophysical_gf_differential)
     [high  |amplified] Ni I gf pool -> BAD_GF (astrophysical_gf_differential)
     [high  |amplified] Si I gf pool -> BAD_GF (astrophysical_gf_differential)
```

_Predictions for Procyon / 55 Cnc come straight from the curated `amplifies_with` axes — they work BEFORE those stars run (walk in knowing the landmines). The auto-layer folds in `procyon_ew_integrity.csv` the moment RYA-281/348/349 land._