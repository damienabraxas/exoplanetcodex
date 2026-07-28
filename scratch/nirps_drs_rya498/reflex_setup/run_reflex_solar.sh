#!/bin/bash
# RYA-498 — headless NIRPS reflex reduction of the 10 HELIOS solar frames.
# Prereqs (Sirius): esoreflex 2.11.2 extracted; openjdk-8-jre + xvfb + python3
# numpy/scipy/astropy/matplotlib (apt); ~/.esoreflex/esoreflex.rc (see sibling);
# nirps.xml install-placeholders substituted (ROOT_DATA_PATH_TO_REPLACE ->
# reflex_solar, CALIB_DATA_PATH_TO_REPLACE -> reflex_calib with
# reflex_calib/nirps-3.3.12 -> calib_static/nirps-calib-3.3.12/cal);
# reflex_solar/reflex_input/nirps/ <- 41 raw frames (10 science + 31 cals).
B=/mnt/codex-data/engines/nirps_drs
R=$B/esoreflex-2.11.2-linux
export LD_LIBRARY_PATH=$B/install/lib:$B/install/lib64
export PATH=$B/install/bin:$PATH
export MPLBACKEND=Agg
xvfb-run -a $R/esoreflex/bin/esoreflex -n nirps \
  -SelectDatasetMethod All -RecipeFailureMode Continue -ProductExplorerMode Disabled \
  -GlobalPlotInteractivity false -GlobalCalibPlotInteractivity false -EraseDirs true
echo "=== REFLEX RUN EXIT=$? ==="
