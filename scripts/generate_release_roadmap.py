#!/usr/bin/env python3
"""Generate public release progress from the canonical product store (RYA-1066)."""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ELEMENTS = ROOT / "data/config/elements_master.json"
PRODUCTS = ROOT / "data/products"
DEFAULT_OUTPUT = ROOT / "data/products/site/release_roadmap_v1.json"
SCHEMA = "codex.release_roadmap/1"

TARGETS = (
    ("solar", "Sun"),
    ("alpha_cen_a", "Alpha Centauri A"),
    ("alpha_cen_b", "Alpha Centauri B"),
)


def canonical_targets() -> list[dict[str, str]]:
    registry = json.loads(ELEMENTS.read_text(encoding="utf-8"))
    targets = []
    for row in registry["elements"]:
        label = row["symbol"]
        targets.append({
            "id": "Fe-II" if label == "Fe II" else label,
            "element": "Fe" if label == "Fe II" else label,
            "ion": "II" if label == "Fe II" else ("I" if label == "Fe" else "any"),
        })
    if len(targets) != registry["n_targets"] or len({t["id"] for t in targets}) != len(targets):
        raise ValueError("canonical element registry count/identity mismatch")
    return targets


def target_progress(star: str, name: str, canonical: list[dict[str, str]]) -> dict:
    completed = []
    for target in canonical:
        path = PRODUCTS / star / f'{target["element"]}.json'
        if not path.is_file():
            continue
        document = json.loads(path.read_text(encoding="utf-8"))
        if document.get("schema") != "codex.element_product/1" or document.get("star") != star:
            raise ValueError(f"invalid product document: {path}")
        current = document.get("products", [])
        # `products` is explicitly the store's CURRENT, publishable channel. Archive,
        # quarantine and gaps are separate arrays and therefore cannot inflate progress.
        if target["ion"] == "any":
            ready = bool(current)
        else:
            ready = any(str(row.get("ion")) == target["ion"] for row in current)
        if ready:
            completed.append(target["id"])
    total = len(canonical)
    return {
        "target_id": star,
        "target_name": name,
        "complete_products": len(completed),
        "total_products": total,
        "progress_fraction": round(len(completed) / total, 6),
        "complete_product_ids": completed,
        "gate": "current_publishable_product",
    }


def generate() -> dict:
    canonical = canonical_targets()
    progress = {star: target_progress(star, name, canonical) for star, name in TARGETS}
    return {
        "schema": SCHEMA,
        "product_id": "codex-public-release-roadmap-v1",
        "sources": {
            "canonical_products": "data/config/elements_master.json",
            "abundance_product_store": "data/products/<target>/<element>.json",
            "publication_gate": "Only rows in a document's current `products` array count; archive, quarantine, gaps, missing files, and superseded products do not.",
        },
        "canonical_product_count": len(canonical),
        "releases": [
            {
                "id": "foundation",
                "name": "Foundation",
                "status": "complete",
                "summary": "Core spectroscopy pipeline, abundance engines, provenance framework, data architecture, validation tooling, and open-science infrastructure established. Ongoing R&D continues as normal scientific development.",
                "progress_fraction": 1.0,
            },
            {
                "id": "alpha",
                "name": "Alpha Release — Full Solar Calibration",
                "status": "in_progress",
                "summary": "Reviewed elemental products for the Sun across the strongest valid wavelength regions and available engines.",
                "notes": ["Fe is nearing full intended band and engine coverage, with remaining infrared due diligence.", "Al is in active completion and audit."],
                "targets": [progress["solar"]],
            },
            {
                "id": "beta",
                "name": "Beta Release — Alpha Centauri A & B",
                "status": "planned",
                "summary": "The solar multi-band, multi-engine method extended independently to both Alpha Centauri components, including Codex-generated 3D model tests where scientifically justified and computationally ready.",
                "targets": [progress["alpha_cen_a"], progress["alpha_cen_b"]],
            },
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--site-root", type=Path)
    args = parser.parse_args()
    product = generate()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(product, indent=2) + "\n", encoding="utf-8")
    if args.site_root:
        destination = args.site_root / "assets/data/release-roadmap.v1.json"
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(args.output, destination)
    print(f'wrote {args.output} ({product["canonical_product_count"]} canonical products)')


if __name__ == "__main__":
    main()
