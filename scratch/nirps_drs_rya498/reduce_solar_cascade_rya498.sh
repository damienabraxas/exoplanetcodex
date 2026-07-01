#!/bin/bash
# RYA-498 — reduce the 10 HELIOS solar SUN,FP,G2V raw frames (2023-04-29, HA) to
# S1D_FINAL via the NIRPS DRS (esorex), hand-driven SOF cascade.
#
# Short-circuits the detector-cal steps with the DRS static reference products
# (master_dark / hot_pixels / bad_pixels / static wave+dll matrices valid
# 2022-11-01..2023-01-01, the latest <= our 2023-04-29 night), and bootstraps the
# night-specific cals (orderdef -> mflat -> wave) from the RYA-498 manifest.
#
# Raw->tag mapping and recipe I/O are from the reflex OCA (nirps_wkf.oca) + the
# esorex --man-page of each recipe. Run on Sirius. Usage:
#   reduce_solar_cascade_rya498.sh <step>   step in {prep,orderdef,mflat,wavefp,wavethar,scired,all}
set -uo pipefail
B=/mnt/codex-data/engines/nirps_drs
P=$B/install
export LD_LIBRARY_PATH=$P/lib:$P/lib64
RD="--recipe-dir=$P/lib/esopipes-plugins/nirps-3.3.12"
ESOREX="$P/bin/esorex $RD --suppress-prefix=TRUE"
CAL=$B/rawcals_20230429                    # decompressed night cals (dp_id.fits)
SCI=$B/rawsci_solar_20230429               # science .fits.Z
STAT=$B/calib_static/nirps-calib-3.3.12/cal
W=$B/reduce_solar                          # workdir
SOF=$W/sof; PROD=$W/products; RAW=$W/raw
mkdir -p "$SOF" "$PROD" "$RAW"

# --- static reference products valid for 2023-04-29 (latest <= date) ---
GEOM=$STAT/NIRPS_CCD_geom_config.fits
ICFG=$STAT/NIRPS_HA_master_inst_config_2022-11-01.fits
HOTP=$STAT/NIRPS_HA_hot_pixels_2023-01-01.fits
BADP=$STAT/NIRPS_HA_bad_pixels_2023-01-01.fits
MDARK=$STAT/NIRPS_HA_master_dark_2023-01-01.fits
SWA=$STAT/NIRPS_HA_STATIC_WAVE_MATRIX_A_2022-11-01.fits
SWB=$STAT/NIRPS_HA_STATIC_WAVE_MATRIX_B_2022-11-01.fits
SDA=$STAT/NIRPS_HA_STATIC_DLL_MATRIX_A_2022-11-01.fits
SDB=$STAT/NIRPS_HA_STATIC_DLL_MATRIX_B_2022-11-01.fits
GMASK=$STAT/NIRPS_G2.fits
EXT=$STAT/NIRPS_EXTINCTION_TABLE.fits
RLA=$STAT/NIRPS_HA_TH_REF_LINE_TABLE_A_2022-11-01.fits
RLB=$STAT/NIRPS_HA_TH_REF_LINE_TABLE_B_2022-11-01.fits
ODM=$STAT/NIRPS_HA_ORDERS_MASK_2022-11-01.fits

man=$CAL/cal_manifest.csv
# files for a given DPR.TYPE. dp_id (field 1) + dp_tech (last) are comma-free;
# dp_type may contain commas and is quoted -> grep the exact delimited form.
files_of() {
  local t="$1" pat
  if [[ "$t" == *,* ]]; then pat="\"$t\""; else pat=",$t,"; fi
  grep -F "$pat" "$man" | cut -d, -f1 | while read -r id; do echo "$CAL/$id.fits"; done
}

# emit "<path> <PRO.CATG>" for a pipeline product (tag downstream inputs by their
# actual PRO.CATG, since product filenames != DO tags)
pt() { local f="$1"; local c; c=$($P/bin/fitshdr "$f" 2>/dev/null | sed -n "s/.*PRO CATG *= *'\([^ ']*\).*/\1/p" | head -1); [ -n "$c" ] && echo "$f $c"; }
# the mflat/orderdef product set that every downstream recipe reuses
prods_common() {
  for pf in NIRPS_ORDER_TABLE_A NIRPS_ORDER_TABLE_B NIRPS_ORDER_PROFILE_A NIRPS_ORDER_PROFILE_B \
            NIRPS_FLAT_A NIRPS_FLAT_B NIRPS_BLAZE_A NIRPS_BLAZE_B; do pt "$PROD/$pf.fits"; done
}

step="${1:-all}"

if [ "$step" = prep ] || [ "$step" = all ]; then
  echo "### PREP: decompress science ###"
  for f in "$SCI"/*.fits.Z; do o="$RAW/$(basename "${f%.Z}")"; [ -s "$o" ] || zcat "$f" > "$o"; done
  echo "science decompressed: $(ls "$RAW"/*.fits 2>/dev/null | wc -l)"
  echo "cals decompressed:    $(ls "$CAL"/*.fits 2>/dev/null | wc -l)"
fi

if [ "$step" = orderdef ] || [ "$step" = all ]; then
  echo "### espdr_orderdef ###"
  s=$SOF/orderdef.sof; : > "$s"
  for f in $(files_of "ORDERDEF,LAMP,DARK"); do echo "$f ORDERDEF_A" >> "$s"; done
  for f in $(files_of "ORDERDEF,DARK,LAMP"); do echo "$f ORDERDEF_B" >> "$s"; done
  { echo "$GEOM CCD_GEOM"; echo "$ICFG INST_CONFIG"; echo "$HOTP HOT_PIXEL_MASK"; echo "$BADP BAD_PIXEL_MASK"; } >> "$s"
  cat "$s"
  ( cd "$PROD" && $ESOREX espdr_orderdef "$s" ) 2>&1 | tail -25
  echo "products:"; ls -t "$PROD"/*ORDER_TABLE* "$PROD"/*ORDER_PROFILE* 2>/dev/null | head
fi

OTA="$PROD/NIRPS_ORDER_TABLE_A.fits"; OTB="$PROD/NIRPS_ORDER_TABLE_B.fits"

if [ "$step" = mflat ] || [ "$step" = all ]; then
  echo "### espdr_mflat ###"
  s=$SOF/mflat.sof; : > "$s"
  for f in $(files_of "FLAT,LAMP,DARK"); do echo "$f FLAT_A" >> "$s"; done
  for f in $(files_of "FLAT,DARK,LAMP"); do echo "$f FLAT_B" >> "$s"; done
  { echo "$GEOM CCD_GEOM"; echo "$ICFG INST_CONFIG"; echo "$HOTP HOT_PIXEL_MASK"; echo "$BADP BAD_PIXEL_MASK";
    echo "$OTA ORDER_TABLE_A"; echo "$OTB ORDER_TABLE_B";
    echo "$SWA STATIC_WAVE_MATRIX_A"; echo "$SWB STATIC_WAVE_MATRIX_B";
    echo "$SDA STATIC_DLL_MATRIX_A"; echo "$SDB STATIC_DLL_MATRIX_B"; } >> "$s"
  ( cd "$PROD" && $ESOREX espdr_mflat "$s" ) 2>&1 | tail -30
  echo "products:"; ls -t "$PROD"/*FSPECTRUM* "$PROD"/*BLAZE* "$PROD"/*ORDER_PROFILE* "$PROD"/*MASTER_FLAT* 2>/dev/null | head
fi

if [ "$step" = wavefp ] || [ "$step" = all ]; then
  echo "### espdr_wave_FP ###"
  s=$SOF/wavefp.sof; : > "$s"
  for f in $(files_of "WAVE,FP,FP"); do echo "$f FP_FP" >> "$s"; done
  { echo "$GEOM CCD_GEOM"; echo "$ICFG INST_CONFIG"; echo "$HOTP HOT_PIXEL_MASK"; echo "$BADP BAD_PIXEL_MASK";
    prods_common;
    echo "$SWA STATIC_WAVE_MATRIX_A"; echo "$SWB STATIC_WAVE_MATRIX_B";
    echo "$SDA STATIC_DLL_MATRIX_A"; echo "$SDB STATIC_DLL_MATRIX_B"; } >> "$s"
  ( cd "$PROD" && $ESOREX espdr_wave_FP "$s" ) 2>&1 | tail -30
  echo "products:"; ls -t "$PROD"/*FP_SEARCHED* "$PROD"/*WAVE_MATRIX* "$PROD"/*DLL_MATRIX* 2>/dev/null | head
fi

if [ "$step" = wavethar ] || [ "$step" = all ]; then
  echo "### espdr_wave_THAR ###"
  s=$SOF/wavethar.sof; : > "$s"
  for f in $(files_of "WAVE,UN1,FP"); do echo "$f THAR_FP" >> "$s"; done
  { echo "$GEOM CCD_GEOM"; echo "$ICFG INST_CONFIG"; echo "$HOTP HOT_PIXEL_MASK"; echo "$BADP BAD_PIXEL_MASK";
    prods_common;
    for w in "$PROD"/*FP_SEARCHED_LINE_TABLE* "$PROD"/*S2D_BLAZE_FP_FP* "$PROD"/*S2D_FP_FP*; do pt "$w"; done
    echo "$RLA REF_LINE_TABLE_A"; echo "$RLB REF_LINE_TABLE_B";
    echo "$ODM ORDERS_MASK"; echo "$MDARK MASTER_DARK";
    echo "$SWB STATIC_WAVE_MATRIX_B"; echo "$SDB STATIC_DLL_MATRIX_B"; } >> "$s"
  ( cd "$PROD" && $ESOREX espdr_wave_THAR "$s" ) 2>&1 | tail -30
  echo "products:"; ls -t "$PROD"/*WAVE_MATRIX* "$PROD"/*DLL_MATRIX* 2>/dev/null | head
fi

if [ "$step" = scired ] || [ "$step" = all ]; then
  echo "### espdr_sci_red (10 SUN,FP,G2V) ###"
  s=$SOF/scired.sof; : > "$s"
  # one science frame per sci_red run is the norm; do all 10 via a loop below.
  # Build the shared cal part once:
  cal=$SOF/scired_cal.part; : > "$cal"
  { echo "$GEOM CCD_GEOM"; echo "$ICFG INST_CONFIG"; echo "$HOTP HOT_PIXEL_MASK"; echo "$BADP BAD_PIXEL_MASK";
    echo "$EXT EXT_TABLE"; echo "$GMASK MASK_TABLE";
    prods_common;
    for w in "$PROD"/*WAVE_MATRIX*A.fits "$PROD"/*WAVE_MATRIX*B.fits \
             "$PROD"/*DLL_MATRIX*A.fits "$PROD"/*DLL_MATRIX*B.fits \
             "$PROD"/*CONTAM* "$PROD"/*REL_EFF* "$PROD"/*ABS_EFF*; do [ -s "$w" ] && pt "$w"; done
  } >> "$cal"
  n=0
  for f in "$RAW"/*.fits; do
    n=$((n+1)); od=$PROD/sci_$(basename "$f" .fits); mkdir -p "$od"
    { echo "$f OBJ_FP"; cat "$cal"; } > "$SOF/scired_$n.sof"
    echo "--- sci_red frame $n: $(basename "$f") ---"
    ( cd "$od" && $ESOREX espdr_sci_red "$SOF/scired_$n.sof" ) 2>&1 | tail -12
    ls -t "$od"/*S1D* 2>/dev/null | head -3
  done
fi
echo "=== step '$step' done ==="
