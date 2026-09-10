#!/bin/bash
# RYA-1208 — the four near-UV Gerber cells, DEEPGRADED tier (--lines-deep-graded),
# matching the invocation the existing near-UV products were built with.
cd "$(dirname "$0")"
export ISPEC_DIR=/mnt/codex-data/engines/ispec_src
PY=/mnt/codex-data/venv312/bin/python
run () {  # ion holding tag
  echo "### $3" >> /tmp/gerber_nuv.log
  $PY scripts/derive_band_products.py --star solar --element Fe --ion "$1" \
      --lo 3000 --hi 3780 --instrument kpno_solar_atlas --holding "$2" \
      --lines-deep-graded --engine-b-deck gerber-1d-lte > "/tmp/g_$3.log" 2>&1
  echo "### $3 exit=$?" >> /tmp/gerber_nuv.log
  grep -E "A\(Fe|Traceback|REFUS|error:" "/tmp/g_$3.log" | tail -3 >> /tmp/gerber_nuv.log
}
run I  solar_kpno_kurucz2005_corrected nuv_kur_FeI
run II solar_kpno_kurucz2005_corrected nuv_kur_FeII
run I  solar_kpno_molecfit_corrected   nuv_mol_FeI
run II solar_kpno_molecfit_corrected   nuv_mol_FeII
echo "NUV DONE" >> /tmp/gerber_nuv.log
