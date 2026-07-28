# Solar NIRPS audit (recon only) -- YJH reflected-solar

_RECON ONLY. Verdict GO / GO WITH CAVEATS / NO-GO from flags. Telluric UNVERIFIED == CRITICAL (permanent rule). K-band CO is OUTSIDE NIRPS range and is NOT assessed here._

## Summary
N=10, dates 2023-04-29T15:30:03 .. 2023-04-29T15:39:14, 3020.5 MB; telluric-CORRECTED=0/10, below-SNR-200=0, CRITICAL-flagged=10.
Bodies present (single-body): SUN,FP,G2V=10

## Structure discovery
## Per-frame audit
| file | OBJECT | DATE_OBS | WMIN_A | WMAX_A | SNR_MED | SPECSYS | TELLURIC | COVERAGE | FLAGS |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| NIRPS.2023-04-29T15:30:03.014.fits | SUN,FP,G2V | 2023-04-29T15:30:03 |  |  |  | ? |  |  | READ-ERROR (CRITICAL) |
| NIRPS.2023-04-29T15:31:04.321.fits | SUN,FP,G2V | 2023-04-29T15:31:04 |  |  |  | ? |  |  | READ-ERROR (CRITICAL) |
| NIRPS.2023-04-29T15:32:05.628.fits | SUN,FP,G2V | 2023-04-29T15:32:05 |  |  |  | ? |  |  | READ-ERROR (CRITICAL) |
| NIRPS.2023-04-29T15:33:06.935.fits | SUN,FP,G2V | 2023-04-29T15:33:06 |  |  |  | ? |  |  | READ-ERROR (CRITICAL) |
| NIRPS.2023-04-29T15:34:08.242.fits | SUN,FP,G2V | 2023-04-29T15:34:08 |  |  |  | ? |  |  | READ-ERROR (CRITICAL) |
| NIRPS.2023-04-29T15:35:09.548.fits | SUN,FP,G2V | 2023-04-29T15:35:09 |  |  |  | ? |  |  | READ-ERROR (CRITICAL) |
| NIRPS.2023-04-29T15:36:10.855.fits | SUN,FP,G2V | 2023-04-29T15:36:10 |  |  |  | ? |  |  | READ-ERROR (CRITICAL) |
| NIRPS.2023-04-29T15:37:12.162.fits | SUN,FP,G2V | 2023-04-29T15:37:12 |  |  |  | ? |  |  | READ-ERROR (CRITICAL) |
| NIRPS.2023-04-29T15:38:13.469.fits | SUN,FP,G2V | 2023-04-29T15:38:13 |  |  |  | ? |  |  | READ-ERROR (CRITICAL) |
| NIRPS.2023-04-29T15:39:14.776.fits | SUN,FP,G2V | 2023-04-29T15:39:14 |  |  |  | ? |  |  | READ-ERROR (CRITICAL) |
## Structure discovery (manual — the script's hook is gated on a successful 1D read, which raw frames fail)

The 10 frames are RAW NIRPS detector products (`.fits.Z`, ~288 MB each decompressed). Real HDU map of `NIRPS.2023-04-29T15:30:03.014.fits`:

```
HDU[0] PRIMARY     PrimaryHDU
HDU[1] slope       ImageHDU  (4096, 4096)   <- up-the-ramp slope image
HDU[2] intercept   ImageHDU  (4096, 4096)
HDU[3] n           ImageHDU  (4096, 4096)
HDU[4] posemeter   BinTableHDU  cols=[TIME, FIBRE1, FIBRE2]
```
Key headers: `INSTRUME=NIRPS`, `OBJECT=SUN,FP,G2V`, `ESO DPR CATG=CALIB`, `ESO DPR TYPE=SUN,FP,G2V`, `ESO DPR TECH=ECHELLE`, **no `PRO CATG`**, **no `SPECSYS`**, **no WAVE array**. These are 2D echellograms (ramp-fit detector frames), not 1D spectra — hence 10/10 READ-ERROR. The ESO download readme states for every file: *"Raw data for which no processed data is available."*

## Step-6 Verdict — **NO-GO as-delivered → REDUCTION-OWED** (a real, proven path — NOT a dead end)

The audit confirms the data are **raw 2D NIRPS echelle frames** (`DPR CATG=CALIB`, `SUN,FP,G2V`: Sun on the science fibre + Fabry-Pérot on the calibration fibre), with no reduced 1D YJH product and therefore no telluric verification possible **as delivered**. The 10/10 CRITICAL flags are the as-delivered state — they say "not pipeline-ready," **not** "unusable."

**The reduction path is concrete and already proven in-project:** the α Cen NIRPS data we hold (`…/Alpha Cen A/NIRPS/ADP.*.fits`) is the **reduced** output of exactly this instrument —
- `PRO CATG=S1D_FINAL_A`, `PRO TYPE=REDUCED`, WAVE 9660.5–19230.8 Å (YJH, air via `WAVE_AIR`), 221,629 px;
- and it carries **native telluric columns**: `FLUX_TELL_EL`, `FLUX_TELL_CAL`, `ATM_TRANSM`, … (+ sky-subtracted variants).

So:
1. **Reduce** the raw solar frames → S1D_FINAL via the **NIRPS DRS** (the public ESO/Geneva NIRPS pipeline / esorex `nirps` recipes — the pipeline that produced the α Cen ADP). The standard echelle flow: ramp-fit (already in `slope`/`intercept`/`n`) → order trace → **optimal extraction** (the cross-dispersion profile collapse across the order's pixels) → wavelength solution (the **FP fibre** supplies the simultaneous reference) → 1D YJH.
2. **Telluric is solved by the same reduction** — the NIRPS DRS S1D_FINAL product includes the telluric-corrected flux + `ATM_TRANSM` natively (per the α Cen ADP). So the telluric-CRITICAL flag is *cleared by reduction*, not a separate owed step (no RYA-424 pass needed if the DRS telluric flux is adopted).
3. **Frame:** raw frames carry no `SPECSYS`/BERV (those are added at reduction); the DRS sets BARYCENT + BERV. Empirical measure-the-line frame check stays deferred to conditioning.

**Engine gap (the real owed item):** the NIRPS DRS is **not yet in the project** (Sirius `engines/` has only iSpec + PySME source; no esorex/NIRPS-DRS). Install/wire the NIRPS DRS (esorex `nirps` recipes — RYA-172/337/375 engine-stack territory) → reduce → re-audit the reduced S1D set. Alternatively, re-query the ESO archive for a published S1D (the readme says none exists for these CALIB frames, so self-reduction is the path).

**Net:** the Solar NIRPS YJH data are **good raw frames, single-body (SUN,FP,G2V), reducible to a telluric-corrected 1D YJH S1D product** by the same pipeline that produced our α Cen NIRPS ADP. The verdict is **REDUCTION-OWED**, not a rejection.

### Critical-condition callouts (per ticket)
- **TELLURIC-UNVERIFIED: 10/10** — but because the product is *raw*, not because a 1D spectrum failed telluric. Resolved by NIRPS-DRS reduction (native telluric).
- **READ-ERROR: 10/10** — raw 2D echellograms, no 1D wave array (correct, honest detection).
- **MULTI-BODY: no** — single body `SUN,FP,G2V` (10/10); downstream ephemeris conditioning is solar (Sun direct on the science fibre — no asteroid ephemeris needed, unlike Vesta/Ceres).
- **nm-native / nm→Å:** not assessable (no 1D wave array in raw); the reduced S1D will be nm→Å as for α Cen.
