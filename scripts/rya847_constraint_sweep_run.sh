#!/usr/bin/env bash
# RYA-847 item 2 — the nine-cell constraint sweep, REBUILT.
#
# 🔴 WHY THIS FILE EXISTS. The original driver lived only in a Sirius home directory and
# wrote its per-line CSVs to an untracked output directory. A routine tree re-sync after
# the sweep deleted both — 6 h 50 m of compute AND the method that produced it, leaving
# the ticket's central claim ("no metric threshold transfers across bands") supported by
# numbers nobody could regenerate. Committing the driver makes the loss recoverable by
# running a script instead of by reconstructing one from memory.
#
# The companion reporter, scripts/rya847_constraint_sweep.py, is READ-ONLY over the CSVs
# this produces: it globs `*_lines.csv` from ONE FLAT directory and takes each cell's
# identity from the filename. Lines files carry the TREATMENT in their name, so the two
# Engine-B decks do not collide here the way products.csv does — the flat layout is
# deliberate and the reporter depends on it.
#
# Each cell banks its output into the artifact store THE MOMENT IT FINISHES (RYA-461),
# because "preserve it afterwards" is what failed the first time: afterwards is exactly
# when the tree gets cleaned.
#
# Sequential by default. Cells are independent processes and iSpec gives each run its own
# mkdtemp scratch, so PARALLEL=N runs N at a time; keep N <= cores.
set -u
cd "$(dirname "$0")/.."
export ISPEC_DIR="${ISPEC_DIR:-/mnt/codex-data/engines/ispec_src}"
export PYTHONUNBUFFERED=1
PY="${PY:-/mnt/codex-data/venv312/bin/python}"
OUT="${OUT:-data/results/rya847_sweep}"
PARALLEL="${PARALLEL:-1}"
mkdir -p "$OUT"

bank () {
  $PY - "$OUT" <<'PYEOF'
import sys
from pathlib import Path
sys.path.insert(0, str(Path.cwd()))
from pipeline.artifact_store import save_artifact
n = sum(bool(save_artifact(p, "data")) for p in sorted(Path(sys.argv[1]).iterdir())
        if p.is_file() and not p.name.startswith("."))
print(f"    [banked] {n} file(s) into the artifact store")
PYEOF
}

cell () {  # tag ion lo hi deck [extra...]
  local tag=$1 ion=$2 lo=$3 hi=$4 deck=$5; shift 5
  if [ -f "$OUT/.done_$tag" ]; then echo "[skip] $tag"; return; fi
  echo "=== [$(date +%H:%M:%S)] $tag ==="
  $PY -u scripts/derive_band_products.py --out "$OUT" \
      --element Fe --ion "$ion" --lo "$lo" --hi "$hi" \
      --instrument kpno_solar_atlas --engine-b-deck "$deck" "$@" \
      >> "sweep_${tag}.log" 2>&1 \
    && { touch "$OUT/.done_$tag"; echo "[ok] $tag"; bank; } \
    || { echo "[FAIL] $tag"; tail -5 "sweep_${tag}.log"; }
}

throttle () { while [ "$(jobs -rp | wc -l)" -ge "$PARALLEL" ]; do wait -n; done; }

# The nine cells: both ions x VIS and red-optical x both Engine-B decks, plus Fe I NIR,
# which is synthesis-only (the EW route's curve-of-growth linelist stops at 9199.9 A).
run_all () {
  throttle; cell "FeI_vis_ts-lte"       I  3800  6910 ts-lte      &
  throttle; cell "FeII_vis_ts-lte"      II 3800  6910 ts-lte      &
  throttle; cell "FeI_ro_ts-lte"        I  6910  9199 ts-lte      &
  throttle; cell "FeII_ro_ts-lte"       II 6910  9199 ts-lte      &
  throttle; cell "FeI_vis_gerber-nlte"  I  3800  6910 gerber-nlte &
  throttle; cell "FeII_vis_gerber-nlte" II 3800  6910 gerber-nlte &
  throttle; cell "FeI_ro_gerber-nlte"   I  6910  9199 gerber-nlte &
  throttle; cell "FeII_ro_gerber-nlte"  II 6910  9199 gerber-nlte &
  throttle; cell "FeI_nir"              I 10000 12935 ts-lte --force-synthesis &
  wait
}

run_all
echo "=== [$(date +%H:%M:%S)] SWEEP COMPLETE ==="
echo "report with: $PY scripts/rya847_constraint_sweep.py --sweep $OUT"
