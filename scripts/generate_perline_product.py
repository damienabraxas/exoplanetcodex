#!/usr/bin/env python3
"""
RYA-870 CLI — emit data/products/<star>/<element>_perline.csv (RYA-489 Section 6).

    python3 scripts/generate_perline_product.py --star solar --element Fe

The generator is a JOIN over committed artifacts and measures nothing; see
pipeline/perline_product.py for the sources and the rules each one carries.

🔴 THERE IS NO --band-products ANY MORE (RYA-1229). It used to take a repeatable list of
roots defaulting to a module constant that named two ticket-scoped snapshots:

    DEFAULT_BAND_PRODUCTS = [ROOT/"data/results/rya847/gated", ROOT/"data/results/rya877"]

The measurements moved to data/results/band_products/ and the constant did not. The
generator went on succeeding while reading 13 of 178 per-line files, so the shipped product
covered ONE of four instruments, carried zero rows for the engine the published headline
rests on, and labelled 932 of its 1039 rows `1D-LTE (ts-lte)` — a vocabulary from the
pre-RYA-906 directory layout that no published product has ever used. Nothing failed,
because nothing compared what was read against what had been published.

⚠️ Repointing that constant at data/results/band_products/ does NOT fix it: that directory
holds every tier, selector and route of every holding, and the emitted key could not tell
them apart. Measured, not assumed — it refuses with "6505 rows share a
(line x instrument x engine) key".

So the product is now a projection OF THE FEED. data/products/<star>/<element>.json is
walked, each product is asked where its own per-line evidence is (provenance.copied_to),
every row carries the full RYA-1127 identity, and a published product whose evidence cannot
be reached is REPORTED BY NAME rather than skipped.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.perline_product import (  # noqa: E402
    PRODUCTS_DIR, PerLineProductError, build_perline_product,
)
from pipeline.perline_sources import PerLineSourceError  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--star", default="solar")
    ap.add_argument("--element", default="Fe")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--damping-source", choices=["synthesis", "linelist"],
                    default="synthesis",
                    help="where the damping constants come from. 'synthesis' (default) is "
                         "the GES list the abundance was actually derived with and needs "
                         "iSpec, so it is Sirius-only; 'linelist' emits from "
                         "linelist_<star>.csv and marks the file NOT replication-grade.")
    a = ap.parse_args()

    try:
        product = build_perline_product(a.star, a.element,
                                        damping_source=a.damping_source)
    except (PerLineProductError, PerLineSourceError) as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        return 1

    out = a.out or (PRODUCTS_DIR / a.star / f"{a.element}_perline.csv")
    product.to_csv(out)

    print(f"=== {a.element} / {a.star} — per-line data product ===")
    for k, v in product.header.items():
        print(f"  {k}: {v}")
    print("\n=== accounting (RYA-844: emitted == measured, never filtered) ===")
    print(json.dumps(product.accounting, indent=2, default=str))

    # 🔴 UNREACHABLE EVIDENCE IS AN OUTCOME, NOT A GAP IN THE PRINTOUT. The defect this
    # generator was rewritten to fix was silent: 33 published products contributed nothing
    # and the run looked clean. Every one of them is named here and in the header.
    un = product.accounting.get("unresolved_products") or []
    if un:
        print(f"\n=== {len(un)} PUBLISHED PRODUCTS HAVE NO REACHABLE PER-LINE EVIDENCE ===")
        for ident, why in un:
            print("  {ion:4} {band:12} {instrument:20} {holding:34} {tier:11} "
                  "{selector:24} {route:11} {treatment}".format(**ident))
            print(f"        {why}")

    print(f"\n[out] {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
