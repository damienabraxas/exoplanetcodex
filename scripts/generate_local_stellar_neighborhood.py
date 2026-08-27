#!/usr/bin/env python3
"""Generate the versioned public local-neighbourhood map data product (RYA-1065)."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import sys
from pathlib import Path

import astropy.units as u
from astropy.coordinates import SkyCoord

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from config.constants import STAR_PARAMS
from pipeline.system_catalog import load_catalog

ASTROMETRY = ROOT / "data/reference/local_neighborhood_astrometry.csv"
ASTROMETRY_PROV = ASTROMETRY.with_suffix(".prov.json")
DEFAULT_OUTPUT = ROOT / "data/products/site/local_stellar_neighborhood_v1.json"
SCHEMA_VERSION = "1.0.0"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def finite(value: object) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def generate(site_root: Path | None = None) -> dict:
    catalog = load_catalog()
    # The catalog is the registry. STAR_PARAMS is an optional physics join, never
    # a membership filter: future targets must not disappear because their
    # fundamental parameters have not yet been adopted.
    eligible = catalog
    astrometry = {row["target_id"]: row for row in read_csv(ASTROMETRY)}
    provenance = json.loads(ASTROMETRY_PROV.read_text(encoding="utf-8"))

    ids = [(row["star_params_key"] or row["website_slug"]) for row in eligible]
    slugs = [row["website_slug"] for row in eligible if row["website_slug"]]
    if len(ids) != len(set(ids)):
        raise ValueError("duplicate eligible target IDs")
    if len(slugs) != len(set(slugs)):
        raise ValueError("duplicate eligible target slugs")
    if set(ids) - ({"solar"} | set(astrometry)):
        raise ValueError(f"registered targets missing fixed-position astrometry: {set(ids) - ({'solar'} | set(astrometry))}")

    targets = []
    for row in eligible:
        target_id = row["star_params_key"] or row["website_slug"]
        slug = row["website_slug"] or None
        route = f"/systems/{slug}/" if slug else None
        published = bool(site_root and route and (site_root / "systems" / slug / "index.html").is_file())
        if published and not slug:
            raise ValueError(f"published target {target_id} has no catalog slug")
        base = {
            "id": target_id,
            "name": row["system_name"],
            "slug": slug,
            "role": row["role"],
            "spectral_type": row["spectral_type"] or None,
            "publication_status": "published" if published else "unpublished",
            "url": route if published else None,
        }
        if target_id == "solar":
            base.update({
                "distance_pc": 0.0, "astrometric_source": "IAU heliocentric origin",
                "galactic_l_deg": None, "galactic_b_deg": None,
                "cartesian_pc": {"x_gc": 0.0, "y_rotation": 0.0, "z_ngp": 0.0},
            })
        else:
            astro = astrometry[target_id]
            values = [astro["ra_deg_j2000"], astro["dec_deg_j2000"], astro["parallax_mas"]]
            if not all(finite(value) for value in values) or float(astro["parallax_mas"]) <= 0:
                raise ValueError(f"invalid astrometry for {target_id}")
            ra, dec = float(values[0]), float(values[1])
            if not (0 <= ra < 360 and -90 <= dec <= 90):
                raise ValueError(f"coordinate out of range for {target_id}")
            distance = 1000.0 / float(astro["parallax_mas"])
            gal = SkyCoord(ra=ra * u.deg, dec=dec * u.deg,
                           distance=distance * u.pc, frame="icrs").galactic
            xyz = gal.cartesian.xyz.to_value(u.pc)
            base.update({
                "distance_pc": round(distance, 6),
                "astrometric_source": f"SIMBAD/CDS ({astro['astrometry_bibcode']})",
                "galactic_l_deg": round(gal.l.deg, 6),
                "galactic_b_deg": round(gal.b.deg, 6),
                "cartesian_pc": {
                    "x_gc": round(float(xyz[0]), 6),
                    "y_rotation": round(float(xyz[1]), 6),
                    "z_ngp": round(float(xyz[2]), 6),
                },
            })
        targets.append(base)

    extent = math.ceil(max(t["distance_pc"] for t in targets))
    product = {
        "schema_version": SCHEMA_VERSION,
        "product_id": "codex-local-stellar-neighborhood-v1",
        "coordinate_system": {
            "input_frame": "ICRS", "input_epoch": "J2000",
            "output_frame": "heliocentric Galactic Cartesian",
            "axes": {"x_gc": "toward Galactic center (l=0°, b=0°)",
                     "y_rotation": "direction of Galactic rotation (l=90°, b=0°)",
                     "z_ngp": "toward north Galactic pole (b=+90°)"},
            "projection": "orthographic x_gc-y_rotation; no visual coordinate nudging",
            "extent_pc": extent,
        },
        "provenance": {
            "membership": "all rows in data/catalog/system_catalog.csv; target ID is star_params_key when present, otherwise website_slug",
            "stellar_parameters": "optional join through config/stars.yaml / STAR_PARAMS; never a membership filter",
            "astrometry": provenance,
            "transformation": "Astropy SkyCoord ICRS to Galactic Cartesian",
            "generator": "python3 scripts/generate_local_stellar_neighborhood.py --site-root /path/to/exoplanetcodex-site",
        },
        "targets": targets,
    }
    validate(product, site_root)
    return product


def validate(product: dict, site_root: Path | None = None) -> None:
    if product.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("unsupported schema version")
    if not product.get("provenance"):
        raise ValueError("missing provenance")
    seen_ids, seen_slugs = set(), set()
    for target in product["targets"]:
        if target["id"] in seen_ids or (target["slug"] and target["slug"] in seen_slugs):
            raise ValueError("duplicate target ID/slug")
        seen_ids.add(target["id"])
        if target["slug"]:
            seen_slugs.add(target["slug"])
        if not finite(target["distance_pc"]) or target["distance_pc"] < 0:
            raise ValueError(f"invalid distance for {target['id']}")
        for value in target["cartesian_pc"].values():
            if not finite(value):
                raise ValueError(f"non-finite Cartesian coordinate for {target['id']}")
        if target["publication_status"] == "published":
            expected = f"/systems/{target['slug']}/"
            if target["url"] != expected or not site_root or not (site_root / expected.lstrip("/") / "index.html").is_file():
                raise ValueError(f"published route does not resolve for {target['id']}")
        elif target["url"] is not None:
            raise ValueError(f"unpublished target has fabricated URL: {target['id']}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--site-root", type=Path)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    product = generate(args.site_root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(product, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    if args.site_root:
        destination = args.site_root / "assets/data/local-stellar-neighborhood.v1.json"
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(args.output, destination)
    print(f"wrote {args.output} ({len(product['targets'])} targets)")


if __name__ == "__main__":
    main()
