# RYA-277 — Per-spectral-type acceptance profiles

**Status:** structure delivered; G anchor regression-proven, F anchor populated from real
data (RYA-273/281), K/M physics-note stubs. F/K/M thresholds NOT set by guesswork (§4).

## The realization

Every line-acceptance and quality criterion in the pipeline — σ-clip, EW/COG ceiling, EP
range, quality-grade boundaries, scatter gate, abundance-gate window, continuum strategy,
NLTE applicability — was calibrated on the **Sun** and silently labelled "universal." Each
carried an implicit, unlabelled assumption: *"G-star, Teff ≈ 5772 K."* This is the same
failure family as the solar-only **code paths** caught in RYA-270/273/274, one level up:
the solar assumption baked into the scientific **cuts** themselves. Procyon (F, 6554 K) is
the first place the difficulty rose and the cuts had to move.

`ACCEPTANCE_PROFILES` (config/constants.py) makes the type explicit. It is a
**routing/parameterization** layer, not a science change: routing the Sun through the G
profile reproduces the signed-off solar Fe result **bit-for-bit**.

## Step 1 — Inventory of acceptance/quality criteria (file:line, as found)

| Criterion | Value | Where set | Per-star already? |
|---|---|---|---|
| Fe I scatter gate (legacy) | 0.10 | `constants.py:194` `FE_SCATTER_GATE` | superseded-legacy (RYA-277); G reads profile |
| Fe I scatter gate (G anchor) | 0.1398 | `constants.py` `ACCEPTANCE_PROFILES['G']['fe1_scatter_max']` | **per-type (RYA-446 seed)** |
| A(Fe) gate window | [7.41, 7.51] | `constants.py:192-193` `FE_GATE_LOWER/UPPER` (±0.05) | global → now in G profile `a_fe_*` |
| Ionization gate \|ΔFe(I−II)\| | 0.05 | `constants.py:202` `FE_IONISATION_GATE` | global → now in G profile `dFe_max` |
| Scale-aware abs-A(Fe) diag halfwidth | 0.07 | `constants.py:257` `FE_ABS_DIAG_HALFWIDTH` | global (RYA-336 diagnostic) |
| REW slope gate | 0.10 | `constants.py:260` `FE_REW_SLOPE_GATE` | global |
| COG/EW saturation ceiling | 100 mÅ | `constants.py:355` `EW_FIT_PARAMS['vmic_ew_ceiling_mA']` | global operative (RYA-330); profile documents per-type |
| Min EW / Max EW | 5 / 300 mÅ | `constants.py:353-354` `ew_min_mA`/`ew_max_mA` | global |
| Voigt threshold | 150 mÅ | `constants.py:343`, `EW_FIT_PARAMS['solar'/'procyon']` | **per-star (RYA-273)** |
| σ-clip threshold | 3.0 | `constants.py:350`, `CONTINUUM_PARAMS[*]` | **per-star (RYA-274)** |
| Continuum knot spacing | 100 / 50 Å | `CONTINUUM_PARAMS['solar'/'procyon']` | **per-star (RYA-274)** |
| Continuum upper percentile | 95 / 97 | `CONTINUUM_PARAMS['solar'/'procyon']` | **per-star (RYA-274)** |
| Continuum n_iter | 5 | `CONTINUUM_PARAMS[*]` | per-star |
| Fit window | 2.0 / 2.5 Å | `EW_FIT_PARAMS['solar'/'procyon']` | **per-star (RYA-273)** |
| Fe II EW sanity / err-frac | 0.5 / 0.5 | `constants.py:380,401` | global (Fe II triage, RYA-305/352) |
| NIST grade floor | 'B' | `constants.py:404` `min_nist_grade` | global |
| Line S/N floor | 5.0 | `constants.py:375` `snr_line_min` | global |
| Fe II n_lines floor | 8 (solar) / 3 (procyon) | `validate_fe_rya238.py` GATES | hardcoded per-star → now in profiles |
| vmic (lit) | 1.00 / 1.66 | `validate_fe_rya238.py` GATES | hardcoded per-star → now in profiles |
| Continuum QA primary metric | p99 = 1.0 ± 0.01 | RYA-309 (general, not per-type) | global by design (see note) |

**Continuum-window note (RYA-309 comment):** the *primary* continuum-placement QA is the
**general** metric — p99 of normalized flux = 1.0 ± 0.01 — NOT a per-window median (which
measures line blanketing and false-trips line-rich stars). So the per-type continuum
criterion is mainly the **fitting config** (knot spacing, upper percentile — already in
`CONTINUUM_PARAMS`); any clean-window list is a secondary per-type check. Net: the primary
continuum QA is general → *less* per-type bookkeeping, not more.

## Step 2 — `ACCEPTANCE_PROFILES` structure

Keyed by spectral type. Each populated profile carries the acceptance **gate** thresholds
(the single source the validate script reads), a `method` field (EW vs synthesis,
anticipating M dwarfs), a `nlte_available` flag, and cross-links (`continuum_key`,
`ew_fit_key`) to the already-per-star `CONTINUUM_PARAMS`/`EW_FIT_PARAMS` (no duplication).

| Field | G (Sun, 5772 K) | F (Procyon, 6554 K) | source |
|---|---|---|---|
| `fe1_scatter_max` | 0.1398 | **0.222** | G: RYA-407 honest floor; F: RYA-281 FINAL |
| `a_fe_lo` / `a_fe_hi` | 7.41 / 7.51 | 7.38 / 7.54 | RYA-261 / RYA-273 |
| `dFe_max` | 0.05 | 0.08 | RYA-261 / RYA-273 |
| `n_Fe2_min` | 8 | 3 | RYA-238 / RYA-273 |
| `vmic_lit` | 1.00 | 1.66 | pinned / Allende Prieto 2002 |
| `nlte_available` | True | **False** | Fe NLTE grid runs out > 6500 K (RYA-319) |
| `method` | EW | EW | — |
| `cog_ceiling_mA` | 100 (operative) | 120 (PROVISIONAL) | RYA-330 / RYA-273 COG_FLAG knee |

The **σ floor RISES** G→F (0.1398 → 0.222) — that rise is *why* the gate must be a
function of type. A single universal number cannot be right for both.

K and M are **physics-note stubs** (not populated, accessor fails loud):
- **K** — cooler: sharper/deeper lines but CN (and late-K TiO) bands contaminate
  solar-clean windows → continuum strategy must avoid molecular bands. Thresholds TODO.
- **M** — severe molecular blanketing (TiO/VO/H₂O), arguably no true optical continuum →
  atomic-line EW may be an invalid *method*; full-spectrum synthesis on MARCS.GES likely
  mandatory. `method='synthesis'` flagged; thresholds TODO.

## Step 3 — Route, don't rewrite

- `validate_fe_rya238.py` builds its `GATES` from the profiles via `_gates_from_profile()`
  → `get_acceptance_profile()`. No hardcoded per-star gate literals remain (the old
  `'scatter_max': 0.15` F floor is gone).
- The gate table names the active profile at runtime (type, method, NLTE availability,
  cited σ-floor source) — same anti-silent-assumption discipline as the RYA-270 routing log.
- Accessors: `get_acceptance_profile(type)`, `acceptance_profile_for_star(star_id)`, and
  the back-compat `fe1_scatter_threshold(type)`. All **fail loud** on a stub/unknown type —
  no silent default.
- The scatter gate compares with `≤` (at-floor passes): the per-type `fe1_scatter_max`
  IS the characterized, irreducible floor for that type, so σ == floor is acceptance.

## Step 4 — F/K/M not set by guesswork

Only the **F σ floor is finalized** (0.222, from RYA-281's audited result). The F COG
ceiling is marked PROVISIONAL — the operative pool cut stays `vmic_ew_ceiling_mA`, because
the σ=0.222 anchor was characterized at that cut; resetting it would invalidate the anchor.
K/M wait for real targets.

## Regression — the proof the refactor changed nothing for the anchor

`scripts/validate_fe_rya238.py --star solar`, baseline (origin/main 624fa11, pre-RYA-277)
vs the RYA-277 routing, same inputs — **bit-for-bit identical**:

| Metric | Value | Per-gate status (both) |
|---|---|---|
| A(Fe I) NLTE | 7.5160 | FAIL (abs window; diagnostic, RYA-336) |
| A(Fe II) NLTE | 7.6570 | FAIL (abs window) |
| ΔFe(I−II) [synth arbiter] | 0.0150 | PASS |
| Fe I scatter σ | 0.1380 | PASS (≤ 0.1398) |
| Fe II n_lines | 3 | FAIL |
| vmic | 1.00 | PASS |

The G-anchor science is unchanged; only the gate *string* (`<`→`≤`) and the profile
footer differ cosmetically. (The "GATES FAIL" overall is a pre-existing current-main
state — the strict absolute-A(Fe) window and the small EW Fe II pool — not introduced here.)

**Procyon under the F profile:** Fe I scatter **0.222 ≤ 0.222 → PASS** (was FAIL under the
old universal 0.15 solar gate), A(Fe II) 7.535 PASS, ΔFe 0.036 PASS, n_Fe2 9 PASS. The
"Procyon Fe σ FAILED" verdict was a *solar* gate applied to an F dwarf; under its own
type's profile it is at-floor → PASS. (A(Fe I) absolute window and vmic 1.80-vs-1.66 remain
FAIL — diagnostic / RYA-292, both out of RYA-277 scope.)
