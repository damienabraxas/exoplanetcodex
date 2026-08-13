# RYA-805 — is telluric correction APPLIED to the 18 CRIRES+ Vesta IDPs?

**Verdict: `telluric_applied = not-applied`, all 18, unanimous.** Confirmed from headers,
corroborated in the flux itself. No corrected variant exists to re-pull. The permanent IR
rule holds in `pipeline/crires_telluric.py` but **not** on the `scripts/measure_band_ew.py`
path — see [The gap](#the-gap-the-rule-is-not-enforced-on-the-path-rya-797-will-use).

Reproduce: `python3 scripts/rya805_telluric_audit.py` (Sirius only — the IDPs live at
`/mnt/codex-data/spectra/vesta/CRIRESPlus`). Per-file evidence:
[`rya805_telluric_status.csv`](rya805_telluric_status.csv).

## Why this ticket existed

RYA-370 asserted "no telluric correction" and RYA-373's spec repeats it, but neither
showed the keyword that proves it. Ryan, 2026-08-12: *"normal CRIRES+ data should carry
telluric correction if memory serves."* An assertion inherited across four tickets is
exactly the kind of thing that is load-bearing and unverified, so this is a confirmation
run — and it confirms the inherited claim rather than overturning it.

## Leg 1 — the recipe chain (header proof)

Identical across all 18 files:

| keyword | value | reading |
|---|---|---|
| `ESO PRO CATG` | `OBS_NODDING_EXTRACTC_IDP` | nodding extraction, not a corrected class |
| `PRODCATG` | `SCIENCE.SPECTRUM` | |
| `ESO PRO REC*` | **exactly one** recipe: `cr2res_obs_nodding` | no `REC2` — the chain ends at extraction |
| `PROCSOFT` | `cr2re/1.6.9` | |
| calibrations | `CAL_DARK_BPM`, `CAL_DETLIN_COEFFS`, `CAL_FLAT_EXTRACT_1D`, `CAL_FLAT_MASTER`, `CAL_WAVE_TW`, `PHOTO_FLUX` | detector + flat + wavelength only; **no telluric standard, no atmospheric model** |
| HDUs | `PRIMARY`, `SPECTRUM` (2 total) | **no transmission / `Recon` HDU** |
| columns | `WAVE FLUX ERR QUAL ORDER DETEC XPOS TRACE` | **no `FLUX_TELL`, no transmission column** |
| `FLUXCAL` / `CONTNORM` | `UNCALIBRATED` / `F` | not flux-calibrated, not normalised |
| `SPECSYS` | `TOPOCENT` | rest frame still owed (RYA-372/796) |

Telluric vocabulary (`molecfit`, `TELLURIC`, `ATM_`, `TRANSMISSION`, `RECON`, `STD_STAR`,
`BEST_FIT_MODEL`, …) across every card of every HDU of all 18 files: **0 genuine hits.**
A naive regex returns 162, and all 162 are false positives — detector board shift
registers (`ESO DET DEV1 BOARD*n* TRANS`), the OB sky-transparency *constraint*
(`ESO OBS AMBI TRANS = 3THN`), and FITS boilerplate. The script filters them explicitly;
counting them as evidence would have inverted the verdict.

> ⚠️ `HIERARCH` (RYA-791): astropy strips the prefix, so looking up
> `HIERARCH ESO PRO CATG` returns empty and **manufactures an absence**. In an audit whose
> finding *is* an absence, that bug is indistinguishable from the answer. All lookups here
> use the bare `ESO ...` form.

**Positive control.** When a telluric correction *is* applied in this repo, the product
says so: `data/audit/crires_co_conditioned/vesta_crires_K_CO_K{2192,2217}_topocent_PROVISIONAL.csv`
(RYA-373/380) carries an **`mtrans`** column — the molecfit transmission model. The 18
IDPs carry no such column. So "not-applied" is not merely "we found no keyword"; it is
"the keyword our own corrected products carry is absent."

## Leg 2 — the absorption is still in the flux

A header proves no recipe *ran*. It cannot prove the absorption is still *there*. Same
target (the Sun), same instrument (CRIRES+), same normalisation, against Elgueta+2026
`sp/Sun_{Y,J,H}_rv.dat`, which **is** telluric-corrected. Percent of good pixels
(`QUAL==0`) below 0.7 of the local continuum:

| window | absorber | control (corrected) | the 18 IDPs | excess |
|---|---|---|---|---|
| **1265–1300 nm** | **O₂ 1.27 µm** | **0.54 %** | **2.02 – 2.81 %** | **4–5×** |
| 1480–1520 nm | H₂O | 5.11 % | 10.13 – 25.03 % | 2–5× |
| 1240–1260 nm | H₂O (dry) | 0.18 % | 0.40 – 0.46 % | ~2× |
| 980–1080 nm | — (dry) | 0.62 % | 0.57 – 0.64 % | none |

The **O₂ 1.27 µm band is decisive**: O₂ is purely terrestrial, so the Sun contributes no
counterpart at all, and the IDPs sit 4–5× deeper than the corrected twin.

> ⚠️ **Window choice is the whole experiment.** A first pass compared *whole settings* and
> found the IDPs ~20× deeper than the control. That was wrong — it compared different
> wavelength ranges, and the "excess" was just whichever band edges each setting happens
> to include. Restricted to the control's own 980–1080 nm the two agree to 0.05 pp,
> because 980–1080 nm is intrinsically dry and discriminates nothing. Only a
> telluric-heavy window that **both** spectra cover is a test. The corrected number to
> quote is 4–5× in the O₂ band, not 20×.

## Leg 3 — the depth tracks the atmosphere (the falsifiable one)

If the absorption is terrestrial, its depth must follow the *weather*. H₂O 1480–1520 nm
across the seven H-band frames, against the header's own precipitable water vapour
(`ESO TEL AMBI IWV START`) and airmass:

| setting | night | IWV | airmass | IWV×AM | % < 0.7 |
|---|---|---|---|---|---|
| H1582 | 2022-11-22 | 2.67 | 1.040 | 2.78 | 10.13 |
| H1559 | 2022-11-23 | 2.68 | 1.051 | 2.82 | 11.26 |
| H1582 | 2022-11-22 | 2.91 | 1.051 | 3.06 | 11.24 |
| H1559 | 2022-11-22 | 2.91 | 1.060 | 3.08 | 12.14 |
| H1567 | 2022-11-24 | 4.32 | 1.090 | 4.71 | 16.40 |
| H1575 | 2022-11-24 | 8.37 | 1.040 | 8.71 | 25.03 |
| H1575 | 2022-11-24 | 8.37 | 1.046 | 8.76 | 24.82 |

**Pearson r = +0.996** (n=7) against IWV×airmass; +0.994 against IWV alone.

The Sun does not know about Paranal's humidity. Absorption depth that is a near-perfect
linear function of the water column above the telescope is atmospheric by construction.
Note **H1575 was observed at IWV = 8.37 mm — 3× the water of the other nights** — and H is
where 53 of the 74 Elgueta-certified Fe I lines live, so RYA-797's worst telluric
contamination sits on its most valuable band.

## Is a corrected variant available? **No.**

ESO TAP (`ivoa.ObsCore`, `dbo.raw`), queried 2026-08-12:

1. **These exposures**: `instrument_name LIKE 'CRIRES%' AND target_name LIKE '%esta%'`
   returns **exactly the 18 `ADP.2025-06-06T16:*` we already hold** — the whole Vesta
   CRIRES+ holding, at any date. No alternative processed product.
2. **The whole collection**: `obs_collection='CRIRESplus'` is **4131 files, one single
   `dataproduct_type`/`obstech` combination** (`spectrum`/`NODDING`). ESO publishes **no
   telluric-corrected CRIRES+ product class at all**, so there is nothing to re-pull for
   any target, not just this one.
3. **No standard-star route either.** Raw calibrations on 2022-11-21…25 are flats (81),
   darks (246), and wavelength references (`LAMP,METROLOGY` 284, `WAVE,UNE` 25,
   `WAVE,FPET` 25, `WAVE,ABSORPTION_N2O` 2, `WAVE,ABSORPTION_SGC` 1). **No telluric
   standard star was observed on any of those nights** — independently reproducing what
   RYA-377 found ("no telluric std").

**Do not re-download.** Nothing better exists. Model-based correction (molecfit) is not
the preferred route; it is the **only** route.

## Does RYA-373 run as-stated, shrink, or become unnecessary?

**Neither as-stated nor unnecessary — it shrinks to a RUN, and the run belongs to RYA-797.**

The telluric requirement is real and confirmed, so RYA-373 does not become unnecessary.
But its *build* is done and on main: `pipeline/crires_telluric.py` (molecfit driver, GDAS
handling, the D1 residual gate) plus `tests/test_crires_telluric_rya373.py`, with banked
output for K2192/K2217. What that banked output covers is **2283.5–2301.2 nm — about
1.2 % of the 947.9–2485.5 nm arm**, K-band CO-overtone only, topocentric, PROVISIONAL.

So the outstanding work is not a molecfit build but a molecfit **run** on the H+J
settings, which is where Fe lives (74 certified Fe I lines: 21 in J, 53 in H, 0 in Y).
Per Ryan's 2026-08-12 reframe that run is folded into **RYA-797**. RYA-373's own residual
blocker is unchanged and unrelated to telluric: RV-insufficiency at sub-floor SNR, needing
a higher-|RV| epoch.

## The gap: the rule is not enforced on the path RYA-797 will use

The ticket asks to confirm the permanent rule — *no IR abundance without confirmed
telluric correction*. It holds in one place and not the other:

- ✅ `pipeline/crires_telluric.py` — `telluric_corrected` defaults `False`;
  `assert_telluric_corrected()` raises `TelluricNotCorrectedError`.
- ✅ `pipeline/telluric_policy.py` — the catalog registers `crires_plus` as
  `telluric_required=yes`, `telluric_basis=correction_required`, and
  `gate('crires_plus', analysis_ready=False)` correctly returns **False**.
- ❌ `scripts/measure_band_ew.py` — the RYA-796 loader RYA-797 will call — **never invokes
  `gate()` or `requires_correction()`.** It imports only `TELLURIC_BANDS` and
  `exclusion()`, and `exclusion()` returns `''` for every CRIRES+ wavelength
  (10280 / 15000 / 22935 Å) because the enumerated band set stops at **11560 Å** — the
  whole J/H/K arm is off the end of the list.

What blocks RYA-797 today is the **rest-frame** gate (`RestFrameNotConditioned`), which
refuses these frames for the RV reason alone. Condition them for RV — exactly what
RYA-372/373 machinery does — and the loader would return telluric-uncorrected H-band flux
with no telluric objection, and `exclusion()` would call every H-band Fe line clean.

**Latent, not live**: no IR abundance has been emitted, so no number is wrong today. But
the guard that should catch this is one RV-conditioning commit away from being bypassed,
and the band at risk (H, IWV 8.37 mm) is the one carrying most of the Fe. This is
**RYA-806**'s scope — `telluric_applied` determined at intake from headers, as the
software switch — and this audit supplies the header basis it needs.

## Sequence note — the RYA-794 holdings rows are not in conflict

`solar_crires_plus_y_rya794` ("our own Vesta IDPs DO NOT EXIST on Sirius… molecfit had
nothing to run on") and `solar_vesta_crires_plus_idp` ("18 IDPs re-pulled to Sirius,
37 MB") are both RYA-794 and read as a live contradiction. They are not: they are
**consecutive states**. The Y product was built from Elgueta's reduced `sp/Sun_Y_rv.dat`
while no IDPs were on Sirius; the 18 IDPs were re-pulled from ESO *afterwards*. Both rows
are now dated so the record reads as a sequence.
