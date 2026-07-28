#!/usr/bin/env bash
# RYA-419 — Sirius data-drive mount guard (the silent-fallback kill).
#
# The trap this guards against: if the 500 GB data drive ever fails to mount,
# /mnt/codex-data is just an empty directory on the 256 GB root SSD. Reads fail
# loudly (fine), but WRITES silently land on root and look fine — until root fills
# or the data is in the wrong place. This script refuses to proceed unless the
# drive is genuinely mounted AND is the right volume (sentinel present).
#
# Usage:
#   bash scripts/check_data_mount.sh        # exit 0 = mounted+verified, else exit 1
#   CODEX_DATA_ROOT=/some/path bash scripts/check_data_mount.sh   # override root
#
# Wiring this into the Python pipeline preflight is DEFERRED to the transfer /
# CODEX_DATA_ROOT work (RYA-264/332) — this issue stands up the guard only.
set -euo pipefail

DATA_ROOT="${CODEX_DATA_ROOT:-/mnt/codex-data}"
SENTINEL="$DATA_ROOT/.codex_mounted"   # exists ONLY on the real data volume

if ! mountpoint -q "$DATA_ROOT"; then
  echo "CRITICAL: $DATA_ROOT is not a mount point — data drive not mounted. Refusing to read/write." >&2
  exit 1
fi

if [[ ! -f "$SENTINEL" ]]; then
  echo "CRITICAL: sentinel $SENTINEL missing — wrong/empty volume at $DATA_ROOT. Refusing." >&2
  exit 1
fi

echo "OK: $DATA_ROOT mounted and verified."
