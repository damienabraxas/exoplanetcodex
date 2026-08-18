#!/usr/bin/env python3
"""
RYA-870 CLI — emit data/products/<star>/<element>_perline.csv (RYA-489 Section 6).

    python3 scripts/generate_perline_product.py --star solar --element Fe

The generator is a JOIN over committed artifacts and measures nothing; see
pipeline/perline_product.py for the sources and the rules each one carries.

⚠️ --band-products may be repeated. It defaults to the per-line sets committed under
data/results/, which is what makes the emitted file reproducible from the repo alone.
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

DEFAULT_BAND_PRODUCTS = [
    ROOT / "data" / "results" / "rya847" / "gated",
    ROOT / "data" / "results" / "rya877",
]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--star", default="solar")
    ap.add_argument("--element", default="Fe")
    ap.add_argument("--band-products", action="append", type=Path, default=None,
                    help="root(s) holding *_lines.csv; repeatable")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--damping-source", choices=["synthesis", "linelist"],
                    default="synthesis",
                    help="where the damping constants come from. 'synthesis' (default) is "
                         "the GES list the abundance was actually derived with and needs "
                         "iSpec, so it is Sirius-only; 'linelist' emits from "
                         "linelist_<star>.csv and marks the file NOT replication-grade.")
    a = ap.parse_args()

    roots = a.band_products or DEFAULT_BAND_PRODUCTS
    try:
        product = build_perline_product(a.star, a.element, roots,
                                        damping_source=a.damping_source)
    except PerLineProductError as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        return 1

    out = a.out or (PRODUCTS_DIR / a.star / f"{a.element}_perline.csv")
    product.to_csv(out)

    print(f"=== {a.element} / {a.star} — per-line data product ===")
    for k, v in product.header.items():
        print(f"  {k}: {v}")
    print("\n=== accounting (RYA-844: emitted == measured, never filtered) ===")
    print(json.dumps(product.accounting, indent=2, default=str))
    print(f"\n[out] {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
