# RYA-498 — reflex prep sequence (the 3 DO fixes) — reproducible

Standing up esoreflex is in runbook §8. Beyond install + config, the Data
Organizer needed three fixes to build datasets from the solar HELIOS frames.
All are the standard/documented path, not Kepler-bookkeeping hacking.

Base dir `B=/mnt/codex-data/engines/nirps_drs`. Prereqs: esoreflex 2.11.2 +
openjdk-8-jre + xvfb + python3 numpy/scipy/astropy/matplotlib; `~/.esoreflex/
esoreflex.rc` (sibling); `nirps.xml` placeholders substituted → `reflex_solar`
(ROOT) + `reflex_calib` (CALIB, with `reflex_calib/nirps-3.3.12` → the static
cal set); 41 raw frames (10 science + 31 night cals) in
`reflex_solar/reflex_input/nirps/`.

## Fix 1 — qualifying long darks (so mdark generates HOT_PIXEL_MASK, VIRTUAL=T)
The science OCA select is `HOT_PIXEL_MASK where VIRTUAL==T` (a generated
product), and NIRPS ships no static substitute for it in the association. The
night's 61.3 s darks fail the DARK predicate (`EXPTIME>830 & NDSAMPLES>=150`);
no qualifying long HA darks exist near 2023-04-29. The **5 HA long darks**
(830.42 s, NDSAMPLES=150, `NIRPS_HA_cal_dark`) that made the shipped static
`hot_pixels_2023-01-01` are in the Jan-2023 demo set (`nirps_extra_raw_calibs`):
```
tar xzf $B/nirps_extra_raw_calibs.tar.gz -C $B/jan_democals
for d in T22:38:42.691 T22:52:44.252 T23:06:45.812 T23:20:47.372 T23:34:48.932; do
  ln -sf $B/jan_democals/nirps_extra_raw_calibs/NIRPS.2023-01-12$d.fits \
         $B/reflex_solar/reflex_input/nirps/
done
```
Sound: detector masks are month-stable; these are the exact darks behind the
static.

## Fix 2 — python sci stack
`sudo apt-get install -y python3-numpy python3-scipy python3-astropy python3-matplotlib`
(reflex python actors `import numpy` at load). Run with `MPLBACKEND=Agg`.

## Fix 3 — target the solar frames (1-line OCA toggle, documented deviation)
The `SUN,FP,G2V → OBJ_FP` rule in `nirps_wkf.oca` has `//REFLEX.TARGET = "T"`
commented out (these are `DPR.CATG=CALIB` HELIOS frames the default workflow
reduces only for the He-line RV product, not full YJH). Uncomment it so our
frames seed a dataset (`.orig` backed up):
```
O=$B/install/share/esopipes/nirps-3.3.12/reflex/nirps_wkf.oca
N=$(awk '/RAW.TYPE = "SUN_FP_G2V"/{f=NR} f && /\/\/REFLEX.TARGET/{print NR; exit}' $O)
sed -i "${N}s|//REFLEX.TARGET|REFLEX.TARGET|" $O
```

## Diagnostic tool — DEBUG logging (the key to all three)
Append to `esoreflex-2.11.2-linux/kepler/resources/log4j.properties` (`.orig`
backed up):
```
log4j.logger.org.eso.util.FileSystemDatasetCreator=DEBUG
log4j.logger.org.eso.DataOrganizer=DEBUG
log4j.logger.org.eso.oca=DEBUG
```
This exposed the per-file `classify` output (`RAW.TYPE`/`REFLEX.CATG`/
`REFLEX.TARGET`) and the association `Discarding surplus … / Found 0` lines that
pinned each cause.

## Then
Run `run_reflex_solar.sh`. The DO builds the dataset tree and reflex runs the
cascade `espdr_mdark → orderdef → mflat → wave_FP → wave_THAR (bootstrap) →
sci_red → S1D_FINAL` with reflex-generated SOFs under
`reflex_solar/reflex_book_keeping/nirps/`.
