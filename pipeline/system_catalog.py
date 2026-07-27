"""
Codex Systems Catalog -- loader, validator, and pointer-resolver.

The Systems Catalog (data/catalog/system_catalog.csv) is the single per-system
index for the Exoplanet Codex: one row per system.

DESIGN RULE -- INDEX BY POINTER, NEVER DUPLICATE A PARAMETER:
    Stellar physics lives ONLY in STAR_PARAMS (config/stars.yaml, exposed by
    config.constants), the single source of truth. This catalog carries
    star_params_key that resolves INTO STAR_PARAMS. It never carries its own
    copy of a physical parameter. The validator enforces this by rejecting any
    physics-like column name.

Naming convention (RYA-631):
    Catalog   = master enumeration, one row per entity (THIS file)
    Register  = mutable current-state facts ledger (State Register)
    Tracker   = per-item work/progress status (element_status_tracker)
    Reference = frozen validated truth values -- "gold" lives ONLY there
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path
from typing import Mapping

# Single source of truth for stellar physics. RYA-631: the ticket brief named
# `pipeline.constants`; the live repo keeps it in `config.constants` (which
# loads config/stars.yaml, RYA-298). Only this import line was adjusted -- no
# parameter is copied anywhere.
try:
    from config.constants import STAR_PARAMS
except Exception as exc:  # loud fail -- no silent fallback
    raise ImportError(
        "system_catalog requires STAR_PARAMS from config.constants "
        "(single source of truth for stellar parameters). Confirm the "
        f"module path/symbol against the live repo. Import error: {exc!r}"
    ) from exc

CATALOG_PATH = Path(__file__).resolve().parents[1] / "data" / "catalog" / "system_catalog.csv"

REQUIRED_COLUMNS = [
    "system_name", "star_params_key", "hd_id", "hip_id", "simbad_id",
    "spectral_type", "role", "is_gbs", "run_order", "n_known_planets",
    "pipeline_status", "linear_parent", "website_slug", "data_archives",
    "source_refs", "notes",
]

ALLOWED_ROLES = {
    "calibration_anchor", "benchmark", "benchmark_candidate",
    "science_target", "future_target",
}
ALLOWED_STATUS = {"queued", "data_audit", "vald_done", "running", "published"}

# Any column name matching one of these tokens would be a duplicated physical
# parameter -- forbidden. Physics is read THROUGH star_params_key.
FORBIDDEN_COLUMN_TOKENS = (
    "teff", "logg", "log_g", "feh", "fe_h", "xi", "vmic", "microturb",
    "abund", "a_fe", "metallicity", "vsini", "vmac",
)
# Roles exempt from the pointer-must-resolve rule (no params upstream yet).
#
# RYA-631: `benchmark_candidate` is exempt alongside `future_target`. The brief
# listed only future_target, but a candidate is by definition a system that has
# NOT yet been through a parameter solve -- neither tau Boo nor HD 77338 has a
# record in config/stars.yaml, and the permanent rules forbid inventing a key
# or adding a made-up STAR_PARAMS entry to force a pass. Promotion of a
# candidate to `benchmark` is exactly the point at which its pointer must
# resolve, and that promotion re-arms this check automatically.
POINTER_EXEMPT_ROLES = {"future_target", "benchmark_candidate"}


def _forbidden_columns(header):
    bad = []
    for col in header:
        low = col.strip().lower()
        if any(tok in low for tok in FORBIDDEN_COLUMN_TOKENS):
            bad.append(col)
    return bad


def _key_resolves(key: str) -> bool:
    if isinstance(STAR_PARAMS, Mapping):
        return key in STAR_PARAMS
    return hasattr(STAR_PARAMS, key)  # adapt if STAR_PARAMS is not a Mapping


def load_catalog(path: Path = CATALOG_PATH):
    """Load + validate. Returns list[dict]. Raises loudly on any defect."""
    if not path.exists():
        raise FileNotFoundError(f"Systems Catalog not found at {path}")
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        header = reader.fieldnames or []
        rows = [dict(r) for r in reader]
    validate_catalog(rows, header)
    return rows


def validate_catalog(rows, header):
    """Structural + provenance validation. Returns rows on success, else raises."""
    errors = []

    if set(header) != set(REQUIRED_COLUMNS):
        missing = sorted(set(REQUIRED_COLUMNS) - set(header))
        extra = sorted(set(header) - set(REQUIRED_COLUMNS))
        if missing:
            errors.append(f"missing columns: {missing}")
        if extra:
            errors.append(f"unexpected columns: {extra}")

    bad_cols = _forbidden_columns(header)
    if bad_cols:
        errors.append(
            f"forbidden physics columns (params live in STAR_PARAMS only): {bad_cols}"
        )

    seen_names, seen_slugs, seen_run_order = set(), set(), {}
    for i, row in enumerate(rows, start=2):  # +2: header is line 1
        name = (row.get("system_name") or "").strip()
        if not name:
            errors.append(f"line {i}: empty system_name")
            continue
        if name in seen_names:
            errors.append(f"line {i}: duplicate system_name '{name}'")
        seen_names.add(name)

        role = (row.get("role") or "").strip()
        if role not in ALLOWED_ROLES:
            errors.append(f"{name}: role '{role}' not in {sorted(ALLOWED_ROLES)}")

        is_gbs = (row.get("is_gbs") or "").strip().lower()
        if is_gbs not in {"true", "false"}:
            errors.append(f"{name}: is_gbs must be true/false, got '{is_gbs}'")

        status = (row.get("pipeline_status") or "").strip()
        if status and status not in ALLOWED_STATUS:
            errors.append(f"{name}: pipeline_status '{status}' not in {sorted(ALLOWED_STATUS)}")

        run_order = (row.get("run_order") or "").strip()
        if run_order:
            if not run_order.isdigit():
                errors.append(f"{name}: run_order must be a positive integer, got '{run_order}'")
            else:
                if run_order in seen_run_order:
                    errors.append(
                        f"{name}: run_order {run_order} duplicated with '{seen_run_order[run_order]}'"
                    )
                seen_run_order[run_order] = name

        slug = (row.get("website_slug") or "").strip()
        if slug:
            if slug in seen_slugs:
                errors.append(f"{name}: duplicate website_slug '{slug}'")
            seen_slugs.add(slug)
            if not all(c.islower() or c.isdigit() or c == "-" for c in slug):
                errors.append(f"{name}: website_slug '{slug}' must be lower-kebab-case")

        # Pointer must resolve for anything carrying real params.
        key = (row.get("star_params_key") or "").strip()
        if role in POINTER_EXEMPT_ROLES:
            continue
        if not key:
            errors.append(f"{name}: star_params_key required for role '{role}'")
        elif not _key_resolves(key):
            errors.append(
                f"{name}: star_params_key '{key}' does not resolve in STAR_PARAMS "
                "(single source of truth). Fix the key or add params upstream."
            )

    if errors:
        raise ValueError("Systems Catalog validation FAILED:\n  - " + "\n  - ".join(errors))
    return rows


def resolve_system(name: str, path: Path = CATALOG_PATH):
    """Return {'catalog': row, 'params': STAR_PARAMS[key]} for one system.

    Physics is read THROUGH the pointer at call time, never stored in the
    catalog. This is the join the website/report layer should consume.
    """
    rows = load_catalog(path)
    match = next((r for r in rows if (r.get("system_name") or "").strip() == name), None)
    if match is None:
        raise KeyError(f"system '{name}' not in catalog")
    key = (match.get("star_params_key") or "").strip()
    params = None
    if key:
        params = STAR_PARAMS.get(key) if isinstance(STAR_PARAMS, Mapping) else getattr(STAR_PARAMS, key, None)
        if params is None and match.get("role") not in POINTER_EXEMPT_ROLES:
            raise KeyError(f"system '{name}': star_params_key '{key}' present but did not resolve")
    return {"catalog": match, "params": params}


def _cli():
    ap = argparse.ArgumentParser(description="Codex Systems Catalog tool")
    ap.add_argument("--validate", action="store_true", help="validate the catalog and exit")
    args = ap.parse_args()
    if args.validate:
        rows = load_catalog()
        print(f"Systems Catalog OK: {len(rows)} systems, {CATALOG_PATH}")
        for role, n in sorted(Counter((r.get("role") or "").strip() for r in rows).items()):
            print(f"  {role}: {n}")
        return 0
    ap.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(_cli())
