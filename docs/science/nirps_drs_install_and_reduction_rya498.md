# NIRPS DRS install + solar reduction runbook (RYA-498)

Stand up the public **NIRPS DRS** on Sirius and reduce the 10 raw HELIOS solar
frames (`Solar Calibration/NIRPS Raw Data/Solar NIRPS/NIRPS.2023-04-29T15:3*.fits.Z`)
to a telluric-corrected 1D YJH `S1D_FINAL`, then re-audit (RYA-497 script).

This is the **build** downstream of the RYA-497 audit (which found the solar
download is RAW, reduction-owed). Reduction path is proven in-project: our α Cen
NIRPS ADP is the reduced `S1D_FINAL_A` output of the same instrument.

Everything lives on **Sirius** (`/mnt/codex-data/...`); the Mac is out of room.

Sources (cite inline, nothing from memory): Mercier et al. 2025 (A&A 700, A8 /
arXiv 2507.21290); Bouchy et al. 2025 (NIRPS on-sky, A&A); Allart et al. 2022
(telluric); ESO NIRPS pipeline kit + manual; Birch & Downs 1994 (vac→air).

---

## 0. Layout on Sirius

```
/mnt/codex-data/engines/nirps_drs/
  nirps-kit-3.3.12-6.tar.gz          # ESO install kit (497M, from ftp.eso.org)
  nirps-kit-3.3.12-6/                # extracted kit + install_pipeline (Perl)
  install/                           # build PREFIX  (bin/esorex, lib, include, nirps-3.3.12/)
  calib/                             # pipeline calibration tree (from nirps-calib-3.3.12)
  nirps_extra_raw_calibs.tar.gz      # 7.4G raw night calibrations (darks/flats/FP/orderdef)
  nirps_static_wave_matrixes_300.tar.gz  # 61M static wavelength matrices
  build.log / dl_*.log
```

ESO FTP source: `ftp://ftp.eso.org/pub/dfs/pipelines/instruments/nirps/`
(kit, `nirps_extra_raw_calibs.tar.gz`, `nirps_static_wave_matrixes_300.tar.gz`,
`nirps-pipeline-manual-3.3.12.pdf`).

---

## 1. Build (EsoRex / CPL kit)

The kit `install_pipeline` is a self-contained Perl builder that compiles, in
order: cfitsio-4.6.2, wcslib-8.4, gsl-2.8, erfa-2.0.1, fftw-3.3.10, cpl-7.4,
esorex-3.13.11, nirps-3.3.12 (recipes) + installs nirps-calib-3.3.12.

Sirius toolchain check (sufficient — bundled libs ship as release tarballs with
pre-generated `configure`, so no autoconf/automake/libtool/m4 bootstrap needed):
perl ✓ make ✓ gcc ✓ g++ ✓ pkg-config ✓ patch ✓ zlib.h ✓ (bzlib.h absent → only
disables optional cfitsio bzip2, not needed). gfortran/cmake/java absent — not
required (all-C build; java only for the deprecated gasgano GUI, skipped).

**Non-interactive invocation gotchas (both real, both hit this session):**

1. The kit refuses to run in a directory where a previous `install_pipeline`
   attempt left artifacts ("Oops, you are trying to rerun the installation").
   Fix: `rm -rf nirps-kit-3.3.12-6 install calib`, re-extract from the tarball,
   then run.
2. `install_pipeline`'s `confirm()` helper hard-requires a tty
   (`-t STDIN || die`). The prefix/calib *path* prompts read fine from a pipe,
   but every `[Y/n]` dies on a non-tty. On the happy path the ONLY `[Y/n]` that
   fires is "directory does not exist, create?" (all others are error-path
   "Should I abort" or gasgano). Fix: **pre-create the prefix + calib dirs** so
   no confirm fires, then feed `prefix\ncalib\n` on stdin and close it.

```bash
BASE=/mnt/codex-data/engines/nirps_drs
mkdir -p "$BASE/install" "$BASE/calib"
cd "$BASE/nirps-kit-3.3.12-6"
printf '%s\n%s\n' "$BASE/install" "$BASE/calib" \
  | perl install_pipeline -makeopts=-j4 > "$BASE/build.log" 2>&1
```

**Gotcha 3 — missing libcurl dev headers.** The kit built cfitsio/wcslib/gsl/
erfa/fftw/cpl/esorex fine, then the **nirps recipe** `configure` died with
`No LIBCURL available` (the recipes need libcurl for the cal-DB/DataLink layer;
esorex itself does not). Sirius had the libcurl runtime but no dev headers. Fix:
`sudo apt-get install -y libcurl4-openssl-dev` (passwordless sudo available),
then resume just the recipe package (no need to rebuild the libs):

```bash
P=$BASE/install
export FFTWDIR=$P ERFADIR=$P GSLDIR=$P CFITSIODIR=$P QFITSDIR=$P CPLDIR=$P TELLURICCORRDIR=$P
export LD_LIBRARY_PATH=$P/lib:$P/lib64:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=$P/lib/pkgconfig:$P/lib64/pkgconfig:$PKG_CONFIG_PATH
export PATH=$P/bin:$PATH
cd "$BASE/nirps-kit-3.3.12-6/nirps-3.3.12"
make distclean; ./configure --prefix=$P && make -j4 && make install
```

### BUILD RESULT — installed + verified ✓

esorex **3.13.11**, CPL **7.4**, nirps recipes **3.3.12**. The full NIRPS/ESPRESSO
recipe set (18 `espdr_*` plugins) registers and runs:

```
espdr_mbias  espdr_single_bias  espdr_mdark  espdr_orderdef  espdr_led_ff
espdr_mflat  espdr_wave_FP  espdr_wave_THAR  espdr_wave_THAR_THAR  espdr_wave_LFC
espdr_wave_LFC_LFC  espdr_wave_TH_drift  espdr_compu_drift  espdr_shifted_extraction
espdr_sci_red  espdr_cal_flux  espdr_cal_eff_ab  espdr_cal_contam
```

```bash
P=$BASE/install
LD_LIBRARY_PATH=$P/lib:$P/lib64 \
  $P/bin/esorex --recipe-dir=$P/lib/esopipes-plugins/nirps-3.3.12 --recipes
```

**Provenance (per RYA-481/461):**
- kit `nirps-kit-3.3.12-6.tar.gz` md5 `ab329a55461ebfc4f1999537df767cec`
- static waves `nirps_static_wave_matrixes_300.tar.gz` md5 `d1b8757bc28e8606c320466c260c3b98`
  (carries `NIRPS_HA_STATIC_DLL_MATRIX_A/B_*` — the DLLDX seed)
- `nirps_extra_raw_calibs.tar.gz` md5 `fe301c105c61a22a1eb65a9fb410ddd0` — **NOT
  our night**: it is the ESO staff (`amodigliops`) demo/reference cal set, 34
  frames dated **2023-01-12/13/15**. A pipeline self-test only, not for our
  2023-04-29 solar reduction.
- all under `/mnt/codex-data/engines/nirps_drs/` (Sirius; nothing big on the Mac).

---

## 2. Findings that change the ticket's assumptions

### Mode: HA, not HE (corrects the comment)

Mercier+2025 used **HE** mode. Our solar frames are **HA**:
`HIERARCH ESO OCS DET1 IMGNAME = NIRPS_HA_SUN`, `INS MODE = HA`,
`OBS PROG ID = 1102.D-0954(K)`. The α Cen ADP is **also HA**
(`INS MODE = HA`). So:
- The night calibrations we associate for 2023-04-29 must be **HA** cals.
- Sun and α Cen are **mode-consistent** (both HA) — good for the differential.

### Version mismatch (the open differential-consistency item — guard #7)

| product            | pipeline       | CPL     |
|--------------------|----------------|---------|
| α Cen archive ADP  | nirps/**3.2.6**| 7.3.2   |
| our install kit    | nirps/**3.3.12**| 7.4    |

Guard #7 (Mercier comment): reduce Sun AND α Cen with the **same DRS version**;
do not mix self-reduced Sun with ESO-archive α Cen. Two ways to satisfy it:

- **(a)** Fetch the matching **3.2.6** kit from the ESO FTP and reduce the Sun
  with it → bit-consistent with the held α Cen ADP. No α Cen re-pull.
- **(b)** Keep 3.3.12 for the Sun and **re-reduce α Cen from raw** with 3.3.12
  (requires pulling the α Cen raw frames + their night cals — a follow-on).

**FINDING: (a) is infeasible.** The ESO FTP
(`ftp://ftp.eso.org/pub/dfs/pipelines/instruments/nirps/`) carries **only**
`nirps-kit-3.3.12-{1..6}` — the 3.2.6 kit that produced the α Cen ADP is no
longer published (ESO keeps only the current 3.3.12 line). So the Sun will be
reduced with **3.3.12**, and guard #7 can only be met via **(b): re-reduce α Cen
from raw with 3.3.12** (pull α Cen raw + HA night cals — a follow-on), or by
**accepting the 3.2.6→3.3.12 delta as a documented caveat** if the same-version
α Cen re-reduction is deferred. The extraction algorithm (Horne optimal) is
unchanged across these point releases, so the differential risk is small but
must be stated, not assumed away.
*Flagged, not silently resolved — the build proceeds on 3.3.12 to prove the
machinery; whether to re-reduce α Cen is the banked-reference decision above.*

### Recipe difference: SUN,FP,G2V (CALIB) vs OBJ_SKY

Solar frames: `DPR CATG=CALIB`, `DPR TYPE=SUN,FP,G2V` (Sun on science fibre A,
**FP** simultaneous reference on fibre B). α Cen ADP: `RAW1 CATG=OBJ_SKY`. The
science recipe + association rules must be driven for the SUN,FP,G2V config
(FP fibre = the simultaneous wavelength reference, per spec §C).

---

## 2b. Night calibrations — acquire via ESO TAP, NOT calSelector (spec §B)

The 2023-04-29 HA night cals were **not** in the RYA-497 download, and the
downloaded `nirps_extra_raw_calibs` is the wrong night (Jan demo set, §1). They
must be fetched. Two findings settle *how*:

- **calSelector auto-association does not work for these frames.** Querying
  `https://archive.eso.org/calselector/v1/associations?dp_id=<solar>&mode=raw2master`
  returns category `CALIB_SUN_FP_G2V` with *"No association could be established:
  this file type is not supported"*. The HELIOS SUN,FP,G2V type is outside the
  standard cal-association tree. (Consistent with RYA-497's "no processed data".)
- **The public ESO TAP service does work** — query the raw archive directly by
  night + `ins_mode='HA'` + `dp_type`. No auth; NIRPS 2023 raw is out of its
  proprietary period, anonymously downloadable from the data portal
  (`https://dataportal.eso.org/dataPortal/file/<dp_id>`, verified — 32 MB
  `.fits.Z` per frame).

Script: `scratch/nirps_drs_rya498/acquire_night_cals_rya498.sh`. Manifest (the
complete self-consistent same-night HA cascade, **31 frames**):
`scratch/nirps_drs_rya498/cal_manifest_20230429_HA.csv` —

| dp_type | n | role |
|---|---|---|
| `DARK` (IMAGE, 61.3s)        | 3  | master dark (NIR has no bias; dark is the zero) |
| `FLAT,DARK,LAMP` (ECHELLE)   | 10 | tungsten flat, fibre A |
| `FLAT,LAMP,DARK` (ECHELLE)   | 10 | tungsten flat, fibre B |
| `FLAT,LAMP,LAMP` (IMAGE)     | 1  | LED flat (blaze) |
| `ORDERDEF,{DARK,LAMP}/{LAMP,DARK}` | 2 | order tracing |
| `CONTAM,DARK,FP` (ECHELLE)   | 1  | contamination check |
| `WAVE,FP,FP` / `WAVE,UN1,FP` / `WAVE,FP,UN1` / `WAVE,UN1,UN1` | 4 | wavelength (FP + **UNe** lamp), taken 14:50–15:12 UT — immediately before the 15:30 solar run |

(The HELIOS day ran 15:30→16:35+ UT, ~60 SUN,FP,G2V frames; we hold 10. Pull more
only if co-added SNR < 200.)

---

## 3. Wavelength / frame conventions → RYA-481 (CRITICAL)

For the NIRPS loader (the convention RYA-481 codifies):

- **Vacuum** wavelengths native. Convert vac→air at the loader boundary
  (Birch & Downs 1994, the RYA-303 convention). The S1D_FINAL carries a
  **`WAVE_AIR`** column — prefer it directly.
- Frame is **STELLAR REST** (BERV + systemic already applied by the DRS).
  **Do NOT re-apply BERV.** This is the **opposite** of the UVES set
  (TOPOCENT, where the loader applies BERV — RYA-272). The loader must branch
  on instrument: NIRPS = already-rest (no-op), UVES = apply-BERV.
- Telluric columns native: `FLUX_TELL_EL`, `FLUX_TELL_CAL`, `ATM_TRANSM`
  (+ sky-subtracted variants).

---

## 4. Reduction recipe (Mercier+2025 §3.1.2)

Steps to copy when building the reference from the DRS S1D:
- Blaze correction from the **Tungsten-lamp** flat.
- Detrend the residual monotonic slope: (a) **DLLDX** product (per-pixel
  wavelength size) for dispersion; (b) a **5778 K blackbody** for the solar IR SED.
- Continuum normalization with **RASSINE** (Cretignier 2020), forced to skip the
  lines of interest (so features are not absorbed into the continuum).
- QA clips: σ-clip S/N and RV at 20σ; drop airmass > 2.
- Co-add the 10 consecutive 1-min frames (15:30–15:39, same target); pull more
  of the 2023-04-29 day only if co-added SNR < 200 floor or telluric
  verification needs multi-airmass leverage.

### 4b. DRS recipe cascade (esorex SOF order) — the remaining build step

Each step is `esorex <recipe> <recipe>.sof`, where the SOF lists the raw inputs
(+ produced master cals from prior steps) with their DO-category tags. Order:

1. `espdr_mdark` — 3× `DARK` → `MASTER_DARK` (NIR: dark serves as the bias zero).
2. `espdr_orderdef` — `ORDERDEF,*` + master dark + static order table → `ORDER_TABLE`.
3. `espdr_mflat` — 20× tungsten `FLAT,*,*` + `FLAT,LAMP,LAMP` (LED) + master dark
   + order table → `MASTER_FLAT` + `BLAZE` (the blaze for step §4 blaze-correction).
4. `espdr_wave_THAR` / `espdr_wave_FP` — `WAVE,UN1,*` (UNe) + `WAVE,FP,*` +
   `NIRPS_HA_STATIC_DLL_MATRIX_*` (the static-wave seed) → `WAVE_MATRIX` +
   `DLL_MATRIX` (the DLLDX product for the §4 dispersion detrend).
5. `espdr_compu_drift` — instrument drift from the simultaneous FP.
6. `espdr_sci_red` — the 10 `SUN,FP,G2V` frames + all masters → **`S2D` → `S1D_FINAL`**,
   telluric-corrected (`FLUX_TELL_*`/`ATM_TRANSM`), BERV/SPECSYS set, WAVE_AIR.

The DRS normally sequences this via the EDPS/reflex workflow (auto-association);
running raw `esorex` requires building each SOF by hand from the manifest above
(the cal DO-category tags are in the pipeline manual,
`nirps-pipeline-manual-3.3.12.pdf`, on Sirius). This is the remaining mechanical
step C; the engine + all inputs are now in place.

---

## 5. Telluric — VERIFY, do not just trust the column (acceptance)

Native correction = **Allart+2022** (HITRAN line-by-line, single-layer,
per-species LM CCF fit) → `FLUX_TELL_*`. Its demonstrated ~2% peak-to-valley
quality is on **visible** F–K spectra; Mercier only *assumed* similar in the
near-IR (the He region is telluric-light). Our S/P/K/CNO targets sit in the
**deep YJH water bands** where that assumption is weakest. So: the DRS clears
"corrected," but quality must be **verified against the Wallace near-IR telluric
atlas** (on disk, RYA-390) in the deep bands before any IR abundance.

ND filter: not a paper concern (the larger blackbody + dispersion slopes are
detrended/normalized, where a smooth ND term also lands). Verify the ND
transmission curve for structure/fringing and eyeball the normalized YJH for
ripple; otherwise benign.

---

## 6. Re-audit (acceptance)

Run `scratch/nirps_solar_audit.py` (RYA-497) on the **reduced** S1D set.
Expected (vs the raw set's 10/10 READ-ERROR): nm→Å fires, telluric-CORRECTED
(named col), SNR assessable, YJH coverage, GO / GO-WITH-CAVEATS. Product is the
solar-NIRPS YJH reference for the IR leg (RYA-162), alongside IAG/Kitt-Peak.

Scope reminders: single body solar (`SUN,FP,G2V`) → **no asteroid ephemeris**
(unlike Vesta/Ceres UVES). 2.3 µm CO / ¹³C is **outside** NIRPS YJH range
(stays CRIRES+/STAGGER, RYA-373). Validate-don't-tune: this is a *reference*,
abundances come later.

---

## 7. Reduction execution status (hand-driven esorex cascade)

Driver: `scratch/nirps_drs_rya498/reduce_solar_cascade_rya498.sh` (step-gated:
`prep|orderdef|mflat|wavefp|wavethar|scired`). Runs on Sirius against
`/mnt/codex-data/engines/nirps_drs/reduce_solar/`. Two install/data gaps fixed
en route:

- **Empty `datastatic`.** `install_pipeline` created the calib dir but never
  populated the static reference set. Extracted `nirps-calib-3.3.12/cal/` (140
  statics: `NIRPS_CCD_geom_config`, `NIRPS_HA_master_inst_config`, hot/bad
  pixels, `master_dark`, `ORDERS_MASK`, `STATIC_WAVE/DLL_MATRIX`,
  `TH_REF_LINE_TABLE`, `NIRPS_G2` = the G2V CCF mask, HITRAN telluric tables) to
  `calib_static/`. Statics valid ≤2023-04-29 are the 2022-11-01 / 2023-01-01 set.
- **One corrupt raw frame.** `NIRPS.2023-04-29T14:50:00.845` (WAVE,FP,FP) was
  the truncated 28 MB partial from the RYA-498 anonymous-DL *test* (the batch
  then skipped it as "exists"). Re-fetched clean (302 MB). Size-scan confirms all
  other cals + science frames intact.

SOF construction: raw→tag via the reflex OCA (`nirps_wkf.oca`) classification;
product→tag via each product's actual `PRO.CATG` (filenames ≠ tags, e.g.
`NIRPS_FLAT_A.fits` carries `FSPECTRUM_A`). Detector-cal steps short-circuited
with the static `master_dark`/`hot_pixels`/`bad_pixels` (the night's 61.3 s darks
don't meet the `mdark` OCA predicate `EXPTIME>830 & NDSAMPLES>=150` — NIRPS long
darks are monthly, not in the science night).

**Cascade result — WORKS through wave_FP; blocked at wave_THAR:**

| recipe | status | products |
|---|---|---|
| `espdr_orderdef` | ✅ 105 s | `ORDER_TABLE_A/B` |
| `espdr_mflat` | ✅ 128 s | `FSPECTRUM_A/B` (master flat), `BLAZE_A/B`, `ORDER_PROFILE_A/B`, backgrounds |
| `espdr_wave_FP` | ✅ 14 s | `FP_SEARCHED_LINE_TABLE_A/B`, `S2D_FP_FP_A/B`, `S2D_BLAZE_FP_FP_A/B` |
| `espdr_wave_THAR` | ❌ blocked | — |
| `espdr_sci_red` | ⛔ gated on wave_THAR | — |

18 valid intermediate products in `reduce_solar/products/`.

**Blocker — `espdr_wave_THAR` edge-order convergence.** With the full, OCA-correct
input set (order tables/profiles, master flat, blaze, FP line tables, S2D FP
blaze, `TH_REF_LINE_TABLE`, `ORDERS_MASK`, `MASTER_DARK`, fibre-B static wave/dll),
the recipe extracts spectra and finds FP lines for most orders but dies:
`No FP line identified for order 1, exiting` →
`espdr_get_all_FP_ll_per_order failed: Input data do not match`
(`espdr_wave_THAR_cal.c:1286`). The absolute UNe/ThAr solution can't bracket the
edge order against the **2022-11-01 static wave first-guess** (the latest the
kit ships — ~6 months stale vs the 2023-04-29 data). `wave_THAR` exposes no
order-range / FP-window relaxation option. This is the wave-solution bootstrap
the operational reflex/EDPS iterates; a single hand-driven pass with a stale
static prior does not converge on the edge order. **No S1D produced yet.**

**Resume options (in preference order):**
1. Install **EDPS**/pyesorex and run the `nirps` workflow — it iterates the wave
   bootstrap and auto-associates (the intended headless path; blocked here only
   because pyesorex/edps aren't installed on Sirius py3.14 — needs a py≤3.12 venv).
2. Seed `wave_THAR` with a **fresher wave prior** — e.g. run `espdr_wave_THAR_THAR`
   on the `WAVE,UN1,UN1` frame first for an FP-independent absolute solution, or
   obtain a same-epoch `WAVE_MATRIX` from a nearby archival reduction.
3. Diagnose the `order 1` gap directly: check whether the wave_FP
   `FP_SEARCHED_LINE_TABLE` populates the edge order, and whether the
   `ORDERS_MASK` order count matches the orderdef `ORDER_TABLE`.

The cascade + all inputs are proven and staged; only the wave-solution
convergence remains between here and `sci_red → S1D_FINAL`. Same blocker will
gate the RYA-500 α Cen re-reduction (same DRS, same wave step).
