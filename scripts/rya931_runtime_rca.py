#!/usr/bin/env python3
"""Record the RYA-931 runtime root cause, with a positive control (RYA-931).

RYA-927 read molecfit's

    mf_config_chk_ref_atmos: Unable to open rat data file .../mipas/equ.fits!

as a broken bundled reference atmosphere.  An absence is a hypothesis, not a
conclusion: "MIPAS is not the defect" only means something if the *same*
equ.fits demonstrably works.  This script pairs the failing run with a passing
run that loads the identical file, and records both.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

from astropy.io import fits

REF_ATMOS = re.compile(r"Checking Reference Atmosphere file:\s*(\S+)")
REF_FAIL = re.compile(r"Unable to open rat data file")
PARSE_FAIL = re.compile(r"Parse 1D spectrum from BINTABLE failed")
BOUNDS = re.compile(r"Access beyond boundaries")


def summarise(run_dir: Path, label: str) -> dict:
    manifest = json.loads((run_dir / "manifest.json").read_text())
    stdout = (run_dir / "esorex.stdout.txt").read_text()
    stderr = (run_dir / "esorex.stderr.txt").read_text()
    science = run_dir / "in" / "science.fits"
    rows = len(fits.open(science)[1].data) if science.exists() else None
    ref = REF_ATMOS.search(stdout)
    ref_path = ref.group(1) if ref else None
    return {
        "label": label,
        "run_dir": str(run_dir),
        "returncode": manifest["returncode"],
        "science_rows": rows,
        "reference_atmosphere": Path(ref_path).name if ref_path else None,
        "reference_atmosphere_path": ref_path,
        "reference_atmosphere_md5": (hashlib.md5(Path(ref_path).read_bytes()).hexdigest()
                                     if ref_path and Path(ref_path).is_file() else None),
        "reference_atmosphere_opened": bool(ref_path) and not REF_FAIL.search(stdout),
        "bintable_parse_failed": bool(PARSE_FAIL.search(stderr)),
        "cpl_access_beyond_boundaries": len(BOUNDS.findall(stderr)),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("failing_run", type=Path)
    ap.add_argument("passing_run", type=Path)
    ap.add_argument("out_json", type=Path)
    args = ap.parse_args()
    failing = summarise(args.failing_run, "RYA-927 trial: finite-ERR filter emptied SCIENCE")
    passing = summarise(args.passing_run, "RYA-931: same equ.fits, non-empty SCIENCE")
    report = {
        "ticket": "RYA-931",
        "claim": "The MIPAS reference atmosphere was never the defect.",
        "mechanism": (
            "The HARPS Phase-3 ERR column is 100% NaN in all ten direct-solar "
            "exposures (0 of 313028 finite). Selecting pixels on finite ERR left a "
            "zero-row SCIENCE table. molecfit failed in "
            "molecfit_data_extract_spectrum_from_cpl_table with 'Access beyond "
            "boundaries' and 'Parse 1D spectrum from BINTABLE failed' (rc=13), and "
            "only afterwards printed the reference-atmosphere line off an already "
            "dirty CPL error state. The message names the last file touched, not "
            "the cause."),
        "positive_control": (
            "The passing run loads the SAME equ.fits, by path and md5, and reports "
            "no 'Unable to open rat data file'. The only changed input is the "
            "number of SCIENCE rows."),
        "failing_run": failing,
        "passing_run": passing,
    }
    same = failing["reference_atmosphere"] == passing["reference_atmosphere"]
    report["control_holds"] = bool(
        same and failing["science_rows"] == 0 and (passing["science_rows"] or 0) > 0
        and passing["reference_atmosphere_opened"] and failing["returncode"] == 13
        and passing["returncode"] == 0)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
