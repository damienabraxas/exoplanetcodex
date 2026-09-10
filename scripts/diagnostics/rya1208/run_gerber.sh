#!/bin/bash
# RYA-1208 — every missing synth-1D-LTE-gerber cell, run sequentially.
cd "$(dirname "$0")"
export ISPEC_DIR=/mnt/codex-data/engines/ispec_src
PY=/mnt/codex-data/venv312/bin/python
run () {  # ion lo hi instrument holding tag
  echo "### $6" >> /tmp/gerber_all.log
  $PY scripts/derive_band_products.py --star solar --element Fe --ion "$1" \
      --lo "$2" --hi "$3" --force-synthesis --instrument "$4" --holding "$5" \
      --lines-tier "$7" --engine-b-deck gerber-1d-lte > "/tmp/g_$6.log" 2>&1
  echo "### $6 exit=$?" >> /tmp/gerber_all.log
  grep -E "A\(Fe|Traceback|REFUS|SystemExit" "/tmp/g_$6.log" | tail -3 >> /tmp/gerber_all.log
}
run I  3000 3780  kpno_solar_atlas solar_kpno_kurucz2005_corrected nuv_kur_FeI  deep
run II 3000 3780  kpno_solar_atlas solar_kpno_kurucz2005_corrected nuv_kur_FeII deep
run I  3000 3780  kpno_solar_atlas solar_kpno_molecfit_corrected   nuv_mol_FeI  deep
run II 3000 3780  kpno_solar_atlas solar_kpno_molecfit_corrected   nuv_mol_FeII deep
run I  9199 11083 iag_fts_solar_atlas solar_iag                    nir_iag      graded
run I  9199 12976 kpno_solar_atlas solar_kpno_molecfit_corrected   nir_kpmol    graded
run I  15009 17491 crires_plus solar_crires_plus_h_rya1094         h_crires     graded
echo "ALL DONE" >> /tmp/gerber_all.log
