#!/usr/bin/env python3
"""Validate repository inputs and, optionally, external synthesis assets."""

from __future__ import annotations

import argparse
import importlib
import json
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def _check(label: str, ok: bool, detail: str, results: list[dict]) -> None:
    results.append({"check": label, "ok": bool(ok), "detail": detail})
    print(f"[{'OK' if ok else 'FAIL'}] {label}: {detail}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--full", action="store_true",
        help="also require iSpec, Turbospectrum, atmosphere grids, and spectra",
    )
    parser.add_argument("--json", type=Path, help="write a machine-readable report")
    args = parser.parse_args()
    results: list[dict] = []

    _check("Python", sys.version_info[:2] in {(3, 9), (3, 10), (3, 11)},
           platform.python_version(), results)
    for module in ("numpy", "scipy", "pandas", "matplotlib", "astropy", "yaml"):
        try:
            imported = importlib.import_module(module)
            version = getattr(imported, "__version__", "unknown")
            _check(f"Python module {module}", True, str(version), results)
        except Exception as exc:
            _check(f"Python module {module}", False, str(exc), results)

    from config.constants import ISPEC_DIR, PATHS, STAR_PARAMS, validate_constants

    try:
        validate_constants()
        _check("canonical constants", True, "config/constants.py", results)
    except Exception as exc:
        _check("canonical constants", False, str(exc), results)

    required_repo = [
        ROOT / "config" / "stars.yaml",
        ROOT / "data" / "method_policy.yaml",
        ROOT / "data" / "catalog" / "system_catalog.csv",
        ROOT / "data" / "catalog" / "instrument_catalog.csv",
        ROOT / "data" / "audit" / "element_status_tracker.csv",
        PATHS["linelist_master"],
        ROOT / "data" / "nlte_grids" / "Fe_Bergemann_MPIA.csv",
        ROOT / "data" / "nlte_grids" / "amarsi2019_cno" / "provenance.json",
    ]
    for path in required_repo:
        _check(f"repository asset {path.relative_to(ROOT)}", path.is_file(),
               "present" if path.is_file() else "missing", results)
    _check("STAR_PARAMS registry", bool(STAR_PARAMS),
           f"{len(STAR_PARAMS)} targets in config/stars.yaml", results)
    from data.catalog.instruments import load_catalog, validate_catalog
    instrument_errors = validate_catalog()
    _check(
        "instrument catalog",
        not instrument_errors,
        f"{len(load_catalog())} records" if not instrument_errors
        else "; ".join(instrument_errors),
        results,
    )

    if args.full:
        ispec_dir = Path(os.environ.get("ISPEC_DIR", ISPEC_DIR)).expanduser()
        _check("ISPEC_DIR", ispec_dir.is_dir(), str(ispec_dir), results)
        full_assets = {
            "ATLAS9.Castelli": ispec_dir / "input/atmospheres/ATLAS9.Castelli",
            "MARCS.GES": ispec_dir / "input/atmospheres/MARCS.GES",
            "Turbospectrum bsyn_lu":
                ispec_dir / "synthesizer/turbospectrum/exec-gf/bsyn_lu",
            "Turbospectrum molecular lists":
                ispec_dir / "input/linelists/turbospectrum/molecules",
        }
        for label, path in full_assets.items():
            _check(label, path.exists(), str(path), results)
        moog = shutil.which("MOOGSILENT") or str(
            ispec_dir / "synthesizer/moog/MOOGSILENT"
        )
        _check("MOOGSILENT (EW comparison)", Path(moog).is_file(), moog, results)
        for label, path in {
            "solar spectra": PATHS["solar_spectra"],
            "Procyon spectra": PATHS["procyon_spectra"],
        }.items():
            count = len(list(path.glob("*.fits"))) if path.is_dir() else 0
            _check(label, count > 0, f"{count} FITS files at {path}", results)

    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip()
    except Exception:
        commit = "unknown"
    report = {
        "repository": str(ROOT),
        "commit": commit,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "mode": "full" if args.full else "core",
        "checks": results,
    }
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(report, indent=2) + "\n")
        print(f"Wrote {args.json}")
    failures = [item for item in results if not item["ok"]]
    print(f"\n{len(results) - len(failures)}/{len(results)} checks passed")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
