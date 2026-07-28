#!/usr/bin/env python3
"""Run a deterministic, repository-native O I NLTE correction example."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from config.constants import SOLAR_ASPLUND2021, get_star_params
from pipeline.nlte_cno import CITATION, cno_nlte_delta, resolve_line, select_leg


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path,
        default=ROOT / "results/tables/quickstart_oi777.json",
    )
    args = parser.parse_args()
    star = get_star_params("solar")
    wave_A = 7771.94
    label = resolve_line("OI", wave_A)
    leg = select_leg(float(star["teff"]))
    a_lte = float(SOLAR_ASPLUND2021["O"])
    delta = cno_nlte_delta(
        "OI", label, float(star["teff"]), float(star["logg"]),
        float(star["feh_ref"]), float(star["xi"]), a_lte, leg=leg,
    )
    if not (-1.0 < delta < 0.0):
        raise RuntimeError(f"unexpected O I 777 nm correction: {delta:+.3f} dex")
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except Exception:
        commit = "unknown"
    grid_dir = ROOT / "data/nlte_grids/amarsi2019_cno"
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "repository_commit": commit,
        "example": "solar O I 777.194 nm NLTE correction",
        "scientific_scope": "correction-grid smoke test; not a spectrum synthesis",
        "inputs": {
            "star_registry": "config/stars.yaml:solar",
            "teff_K": float(star["teff"]),
            "logg": float(star["logg"]),
            "feh": float(star["feh_ref"]),
            "vmic_kms": float(star["xi"]),
            "wavelength_air_A": wave_A,
            "grid_label": label,
            "lte_abundance_A_O": a_lte,
            "leg": leg,
        },
        "result": {
            "delta_nlte_minus_lte_dex": float(delta),
            "corrected_abundance_A_O": a_lte + float(delta),
        },
        "assets": {
            path.name: sha256(path)
            for path in sorted(grid_dir.glob("table[13].dat*"))
        },
        "citation": CITATION,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        f"O I {wave_A:.2f} A ({leg} NLTE): "
        f"A(O) {a_lte:.3f} + {delta:+.3f} = {a_lte + delta:.3f}"
    )
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
