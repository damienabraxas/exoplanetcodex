# RYA-494 — Generalized IR telluric module + α Cen A IR conditioning

## What was built
`pipeline/ir_telluric.py` — the per-star / per-instrument generalization of the RYA-373
Vesta CRIRES+ telluric logic. α Cen A is caller #1; 55 Cnc and every future IR target
inherit it. The strategic point: benchmark work builds the software for all science targets.

Parameterized what RYA-373 hardcoded to Vesta/CRIRES-K:
- **Instrument-agnostic.** CRIRES+ (cr2res IDP: nm, TOPOCENT, no telluric → molecfit,
  reusing the proven `crires_telluric` driver) **and** NIRPS (geneva S1D: already
  telluric-corrected via `FLUX_TELL_*` + `ATM_TRANSM`, already BARYCENT + Å → select the
  column, no molecfit).
- **Star-agnostic.** `IRTarget` descriptor + registry (`vesta`, `alpha_cen_a`).
- **Branched velocity step** (the RYA-373→494 split): `VelocityMode.REFLECTED_SOLAR`
  (asteroid ephemeris, Vesta, RYA-372) vs `VelocityMode.STELLAR` (a direct star's own
  BERV → barycentric, + systemic RV → stellar rest). α Cen A uses STELLAR — **not** asteroid.
- **Permanent IR rules preserved:** telluric in the TOPOCENTRIC frame before any RV shift;
  no number without `telluric_corrected`; nm→Å at the loader; no blind cross-epoch coadd.
- **Generalized telluric-specific gate** (star-agnostic): RYA-373's D1 gate masked solar
  lines; here we detrend (median filter) and compare residual scatter at telluric pixels vs
  clean continuum, **pooled across the frame's segments** (CRIRES tellurics are bimodal per
  segment). The excess is the telluric-specific misfit.

## α Cen A IR — applied + verified
- **CRIRES+: 6 frames** (Y1029/J1232/H1582/H1559/K2192/K2148, 949.6–2485 nm, TOPOCENT).
  molecfit telluric on Y1029 (SNR 302, 1100–1120 nm H₂O), real 3-hourly GDAS:
  **telluric-specific excess residual = 1.58% → PASS (<2%)** (n_telluric=1691, n_clean=2996).
  Stellar-RV branch: BERV +12.8 + systemic −22.3 km/s → α Cen A rest. Conditioned product
  written to `data/audit/acen_a_ir_rya494/acen_a_CRIRES_Y1029_telluric_clean_rya494.fits`.
- **NIRPS: 28 frames**, telluric column `FLUX_TELL_CAL`, SPECSYS=BARYCENT, 9661–19231 Å.
  Telluric already applied by the DRS; stellar branch removes only the systemic RV (no
  double-BERV).

## Findings (honest scope)
1. **🚩 O I 844.6 / 926.6 nm are NOT covered by α Cen A IR.** Both fall *below* the
   CRIRES-Y blue edge (949.6 nm) and the NIRPS blue edge (966.1 nm) — a coverage gap, not a
   telluric gap. The brief's headline atomic lines are not reachable in the in-hand α Cen A
   IR. (They sit in the optical/NIR seam; UVES red, not the IR arm, is the only in-hand path
   — separate question.) Atomic IR lines that *are* covered (e.g. Si I 1.0827 µm, Fe I /
   Mg I in the NIRPS YJH) remain available on the telluric-clean spectra.
2. **🚩 NIRPS attribution INVERTS RYA-479.** The RYA-423 IR RV star-ID (`verdict=A`,
   obs_rv ≈ −26.2 matches α Cen A −25.8, not B −18.4) is authoritative over the *mislabeled*
   `OBJECT='AlphaCenB'` header. So those 28 NIRPS frames are **α Cen A, not B**. RYA-479's
   OBJECT-based optical re-split put them under B; for the IR arm that is wrong. α Cen B has
   **no** NIRPS (RYA-439). The module attributes NIRPS by the star-ID, never the header.
3. **¹³C/CO stays STAGGER-walled.** Molecular CO ¹³C/¹²C needs the Amarsi 3D STAGGER model
   (the RYA-373 collaborator gate). Deferred, flagged — not silently dropped. molecfit
   (telluric) is *not* walled and ran here.
4. **α Cen B IR is a data-acquisition gap** (RYA-439), not a molecfit gap.

## Status
Background task — does **not** block RYA-493. Branch only, not merged. Tests:
`tests/test_ir_telluric_rya494.py` (8) + RYA-373/372 regression (45) green.
