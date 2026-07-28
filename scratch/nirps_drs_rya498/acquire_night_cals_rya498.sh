#!/bin/bash
# RYA-498 — acquire the 2023-04-29 HA night calibrations for the Solar NIRPS
# reduction, via the PUBLIC ESO TAP service (no auth; NIRPS 2023 raw is out of
# its proprietary period).
#
# WHY TAP and not calSelector: the ESO calSelector auto-association REST API
# returns "this file type is not supported" for the HELIOS SUN,FP,G2V frames
# (category CALIB_SUN_FP_G2V) — it cannot build the cal cascade for them. So we
# query the raw archive directly by night + INS MODE + DPR TYPE, which works.
#
# Run on Sirius (the cruncher; never the Mac). Frames land as *.fits.Z (~32 MB
# each, 31 frames ≈ 1 GB). Solar science frames (10× SUN,FP,G2V, on the Mac) are
# staged separately.
set -euo pipefail
BASE="${1:-/mnt/codex-data/engines/nirps_drs/rawcals_20230429}"
mkdir -p "$BASE"; cd "$BASE"

# 1) Query: all 2023-04-29 HA daytime CALIB frames (exclude the SUN science).
#    Night cal block runs ~13:31–15:12 UT (mjd 60063.55–60063.64). The WAVE
#    (FP + UNe) set at 14:50–15:12 immediately precedes the 15:30 solar run.
curl -s --connect-timeout 30 "https://archive.eso.org/tap_obs/sync" \
  --data-urlencode "REQUEST=doQuery" --data-urlencode "LANG=ADQL" \
  --data-urlencode "FORMAT=csv" --data-urlencode "MAXREC=200" \
  --data-urlencode "QUERY=SELECT dp_id, dp_type, dp_tech FROM dbo.raw \
WHERE instrument='NIRPS' AND dp_cat='CALIB' AND ins_mode='HA' \
AND mjd_obs BETWEEN 60063.55 AND 60063.64 AND dp_type NOT LIKE '%SUN%' \
ORDER BY dp_type, dp_id" > cal_manifest.csv
echo "manifest: $(( $(wc -l < cal_manifest.csv) - 1 )) cal frames"

# 2) Download each frame anonymously from the ESO data portal.
tail -n +2 cal_manifest.csv | sed 's/"//g' | cut -d, -f1 | while read -r id; do
  [ -z "$id" ] && continue
  out="$id.fits.Z"
  [ -s "$out" ] && { echo "skip $id"; continue; }
  if curl -sLf --connect-timeout 30 -o "$out" \
       "https://dataportal.eso.org/dataPortal/file/$id"; then
    echo "ok $id ($(stat -c%s "$out") B)"
  else
    echo "FAIL $id"; rm -f "$out"
  fi
done
echo "=== done: $(ls -1 *.fits.Z 2>/dev/null | wc -l) frames on disk ==="
md5sum *.fits.Z > cal_md5sums.txt
