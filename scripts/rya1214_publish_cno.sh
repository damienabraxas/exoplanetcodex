#!/usr/bin/env bash
# RYA-1214 — bring the CNO band products back from Sirius and publish them.
#
# `data/results/band_products/` is GITIGNORED and the products were produced on Sirius,
# which is the only box that can run iSpec. So the artifacts are COPIED here and then
# published through `scripts/publish_product.py`, which is the only thing that may write
# `data/products/solar/<El>.json` (RYA-1034: "a product is published HERE or it does not
# exist").
#
# 🔴 --host / --origin-path ARE NOT DECORATION. The sha256 is taken of the bytes read
# locally, but the record must say where the number was PRODUCED. Recording the path of a
# local copy would state that a Sirius-only result was produced on the Mac, which is the
# laundering RYA-772 names.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REMOTE_REPO="${REMOTE_REPO:-\$HOME/scratch/rya1214}"
REMOTE="${REMOTE:-sirius}"
DEST="$REPO/data/results/band_products"
mkdir -p "$DEST"

echo "== copying CNO band products from $REMOTE =="
rsync -a --stats \
  --include='CI_*' --include='NI_*' --include='OI_*' --exclude='*' \
  "$REMOTE:$REMOTE_REPO/data/results/band_products/" "$DEST/"

cd "$REPO"
ORIGIN_DIR="$(ssh "$REMOTE" "cd ${REMOTE_REPO} && pwd")/data/results/band_products"

published=0
for csv in "$DEST"/{CI,NI,OI}_*_SYNTH*_products.csv; do
  [ -e "$csv" ] || continue
  stem="$(basename "$csv" _products.csv)"
  # holding is the stem segment after the instrument id, and it is NEVER inferred from
  # the filename's words: two holdings of one instrument are two different PRODUCTS
  # (RYA-1026), so it is read back out of the stem exactly as the run wrote it.
  case "$stem" in
    *_kpno_solar_atlas_*)        inst=kpno_solar_atlas ;;
    *_harps_*)                   inst=harps ;;
    *_iag_fts_solar_atlas_*)     inst=iag_fts_solar_atlas ;;
    *) echo "SKIP $stem — unrecognised instrument segment"; continue ;;
  esac
  holding="${stem#*_${inst}_}"; holding="${holding%%_SYNTH*}"
  # The SELECTOR is part of the product KEY (RYA-984), and it is what distinguishes an
  # AGSS21-named-set run from the depth-selected one on the same cell. Read back out of
  # the stem the run wrote, never assumed.
  sel=""; case "$stem" in *_SYNTH_*) sel="${stem#*_SYNTH_}";; esac
  # ⚠️ NOT `--line-set asplund`. That name is REGISTERED for the AGSS21 *Fe* table
  # (reference_lineset.SETS['asplund'], element="Fe"), so claiming it here would assert
  # membership in a set this product was not measured on. The selector carries the fact
  # instead; registering an `asplund-cno` set is a separate, stated decision.
  echo "== publishing $stem  (holding=$holding${sel:+, selector=$sel})"
  python3 scripts/publish_product.py \
    --from "$csv" \
    --holding "$holding" \
    --tier ALL \
    ${sel:+--selector "$sel"} \
    --route SYNTH \
    --host Sirius \
    --origin-path "$ORIGIN_DIR/$(basename "$csv")" \
    "$@"
  published=$((published + 1))
done
echo "== $published product artifact(s) offered to the publisher =="
