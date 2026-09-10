#!/bin/bash
# RYA-1208 tail: (a) the H cell with the CANONICAL bounds the existing H product used
# (15007.11-17493.69, no --force-synthesis) -- the queue ran it with narrowed bounds,
# which would mint a different band identity; (b) the two coverage N/A cells, ATTEMPTED
# so the refusal is the harness's own words rather than my inference.
cd "$(dirname "$0")"
export ISPEC_DIR=/mnt/codex-data/engines/ispec_src
PY=/mnt/codex-data/venv312/bin/python
L=/tmp/gerber_tail.log

echo "### h_crires_canonical" >> $L
$PY scripts/derive_band_products.py --element Fe --ion I --lo 15007.11 --hi 17493.69 \
    --instrument crires_plus --holding solar_crires_plus_h_rya1094 --lines-tier graded \
    --engine-b-deck gerber-1d-lte > /tmp/g_h_canon.log 2>&1
echo "### h_crires_canonical exit=$?" >> $L
grep -E "A\(Fe|Traceback|REFUS" /tmp/g_h_canon.log | tail -3 >> $L

echo "### NA_PROBE harps_nearUV" >> $L
$PY scripts/derive_band_products.py --star solar --element Fe --ion I --lo 3000 --hi 3780 \
    --instrument harps --holding solar_harps_molecfit_corrected --lines-deep-graded \
    --engine-b-deck gerber-1d-lte > /tmp/g_na_harps_nuv.log 2>&1
echo "### NA_PROBE harps_nearUV exit=$?" >> $L
tail -4 /tmp/g_na_harps_nuv.log >> $L

echo "### NA_PROBE kur2005_NIR" >> $L
$PY scripts/derive_band_products.py --star solar --element Fe --ion I --lo 9199 --hi 12976 \
    --force-synthesis --instrument kpno_solar_atlas --holding solar_kpno_kurucz2005_corrected \
    --lines-tier graded --engine-b-deck gerber-1d-lte > /tmp/g_na_kur_nir.log 2>&1
echo "### NA_PROBE kur2005_NIR exit=$?" >> $L
tail -4 /tmp/g_na_kur_nir.log >> $L
echo "TAIL DONE" >> $L
