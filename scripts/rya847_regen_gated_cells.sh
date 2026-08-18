#!/usr/bin/env bash
# RYA-847 item 6 — regenerate the cells the NON-MINIMUM check touches, and only those.
#
# Four cells hold a line the gate excludes; the other five synthesis cells are
# deliberately NOT re-run. Regenerating a cell nothing changed in churns its artifact and
# buries the real diff — RYA-845 lost a finding exactly that way.
#
# 🔴 EACH DECK WRITES ITS OWN --out DIRECTORY, and this is not a preference.
# derive_band_products names its products file by band + instrument + route; the DECK is
# NOT in the name. Point both VIS decks at one directory and the second run silently
# overwrites the first products.csv, taking the ENGINE-B row with it — which is the cell
# the gate does the most work in (5333.768, stat bar 0.0323 -> 0.0193). This is the trap
# scripts/rya783_run_matrix.sh already documented in its header ("a shared directory
# would have the second run silently overwrite the first and the matrix would lose a
# column"); the first cut of this script did not follow it and lost exactly that row.
# near-UV and NIR have one route each, so they stay at the top level.
#
# Sequential, never parallel: `_fit_synth_flux` defaults to a SHARED
# /tmp/ispec_codex_synth, so two synthesis jobs read each other's Turbospectrum scratch
# (the RYA-785 stale-workdir class).
#
# Sirius only (grids + iSpec live there). Resumable: a cell with its `.done_` marker is
# skipped, so a killed run picks up where it stopped.
set -u
cd "$(dirname "$0")/.."
export ISPEC_DIR="${ISPEC_DIR:-/mnt/codex-data/engines/ispec_src}"
export PYTHONUNBUFFERED=1
PY="${PY:-/mnt/codex-data/venv312/bin/python}"
OUT="${OUT:-data/results/rya847/gated}"
mkdir -p "$OUT" "$OUT/ts-lte" "$OUT/gerber-nlte"

run () {  # tag outdir [args...]
  local tag=$1 dir=$2; shift 2
  if [ -f "$OUT/.done_$tag" ]; then echo "[skip] $tag"; return; fi
  echo "=== [$(date +%H:%M:%S)] $tag -> $dir ==="
  $PY -u scripts/derive_band_products.py --out "$dir" "$@" >> "regen_${tag}.log" 2>&1 \
    && { touch "$OUT/.done_$tag"; echo "[ok] $tag"; } \
    || { echo "[FAIL] $tag"; tail -5 "regen_${tag}.log"; }
}

run "FeI_vis_ts-lte"      "$OUT/ts-lte"      --element Fe --ion I --lo 3800 --hi 6910 \
    --instrument kpno_solar_atlas --engine-b-deck ts-lte
run "FeI_vis_gerber-nlte" "$OUT/gerber-nlte" --element Fe --ion I --lo 3800 --hi 6910 \
    --instrument kpno_solar_atlas --engine-b-deck gerber-nlte
run "FeI_nearuv"          "$OUT"             --element Fe --ion I --lo 3000 --hi 3780 \
    --instrument kpno_solar_atlas --force-synthesis
run "FeI_nir"             "$OUT"             --element Fe --ion I --lo 10000 --hi 12935 \
    --instrument kpno_solar_atlas --force-synthesis
echo "=== [$(date +%H:%M:%S)] REGEN COMPLETE ==="
