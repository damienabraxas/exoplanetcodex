#!/usr/bin/env bash
# RYA-783 — the Fe product matrix, run SEQUENTIALLY.
#
# Sequential is not conservatism: `_fit_synth_flux` defaults to a SHARED
# /tmp/ispec_codex_synth, so two synthesis jobs at once read each other's Turbospectrum
# scratch files. That is the RYA-785 stale-workdir class of defect.
#
# Each deck writes its OWN --out directory. Both decks emit the same
# {species}_{band}_..._products.csv filename, so a shared directory would have the second
# run silently overwrite the first and the matrix would lose a column.
set -u
cd /mnt/codex-data/codex/rya783
export ISPEC_DIR=/srv/codex/engines/ispec_src
export PYTHONUNBUFFERED=1
PY=/mnt/codex-data/venv312/bin/python

run () {  # ion lo hi deck
  local ion=$1 lo=$2 hi=$3 deck=$4
  local tag="Fe${ion}_${lo}_${hi}_${deck}"
  local out="data/results/band_products/${deck}"
  echo "=== ${tag} $(date +%H:%M:%S) ==="
  timeout 21600 $PY scripts/derive_band_products.py \
      --element Fe --ion "$ion" --lo "$lo" --hi "$hi" \
      --instrument kpno_solar_atlas --engine-b-deck "$deck" \
      --out "$out" > "/tmp/mx_${tag}.log" 2>&1
  local rc=$?
  grep -E "A=|excluded=|not-served" "/tmp/mx_${tag}.log" | sed 's/^/    /'
  echo "    rc=${rc}"
}

# cheapest first, so a structural problem surfaces early rather than after hours
run II 6910 9199 ts-lte
run II 6910 9199 gerber-nlte
run II 3800 6910 ts-lte
run II 3800 6910 gerber-nlte
run I  6910 9199 ts-lte
run I  6910 9199 gerber-nlte
run I  3800 6910 ts-lte
run I  3800 6910 gerber-nlte
echo "=== MATRIX COMPLETE $(date +%H:%M:%S) ==="
