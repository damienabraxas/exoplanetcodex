#!/usr/bin/env bash
# RYA-419 — create the canonical folder skeleton on the Sirius data drive.
#
# Run on Sirius AFTER the drive is mounted at /mnt/codex-data (Steps 2-5). It refuses
# to run unless the mount guard passes, so the skeleton can never be created on the
# 256 GB root by mistake. Idempotent (mkdir -p); safe to re-run. Creates NO per-target
# dirs beyond the audited Milestone-1 stars — additional targets (e.g. hd209458) get
# their dir when that target's data is staged.
#
# Usage:  bash scripts/init_codex_data_skeleton.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_ROOT="${CODEX_DATA_ROOT:-/mnt/codex-data}"

# Mount guard is mandatory — never write to /mnt/codex-data unmounted (RYA-419).
bash "$HERE/check_data_mount.sh"

cd "$DATA_ROOT"
mkdir -p spectra/{sol,procyon,alpha_cen_a,alpha_cen_b,tau_boo,55cnc_a}
mkdir -p spectra/_quarantine
mkdir -p linelists/{vald_raw,built}
mkdir -p grids/model_atmospheres/{atlas9_castelli,marcs_ges}
mkdir -p grids/nlte/{amarsi2019_cno,mpia_inspect,pysme,gerber_ts}
mkdir -p outputs

echo "OK: canonical skeleton present under $DATA_ROOT"
find "$DATA_ROOT" -maxdepth 2 -type d | sort
