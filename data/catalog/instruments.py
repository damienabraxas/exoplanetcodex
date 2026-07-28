"""Loader and validation contract for the canonical instrument catalog."""

from __future__ import annotations

import csv
from pathlib import Path
from urllib.parse import urlparse

CATALOG = Path(__file__).with_name("instrument_catalog.csv")
STATUSES = {"active", "experimental", "candidate", "rejected", "reference"}
INGEST = {"loader_ready", "loader_partial", "audit_only", "candidate",
          "reference_only", "rejected"}
REQUIRED = {
    "instrument_id", "instrument_name", "archive_name", "archive_url",
    "bands_supported", "wavelength_min_nm", "wavelength_max_nm",
    "resolving_power_min", "resolving_power_max", "pipeline_ingest_status",
    "reduction_requirements", "science_role", "codex_status",
    "source_issue_ids", "source_documentation", "last_verified_date",
}


def load_catalog(path: Path = CATALOG) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def validate_catalog(path: Path = CATALOG) -> list[str]:
    rows = load_catalog(path)
    errors: list[str] = []
    if not rows:
        return ["catalog is empty"]
    missing = REQUIRED - set(rows[0])
    if missing:
        errors.append(f"missing columns: {sorted(missing)}")
    ids: set[str] = set()
    for line, row in enumerate(rows, start=2):
        ident = row.get("instrument_id", "")
        if not ident or ident in ids:
            errors.append(f"line {line}: blank/duplicate instrument_id {ident!r}")
        ids.add(ident)
        if row.get("codex_status") not in STATUSES:
            errors.append(f"{ident}: invalid codex_status")
        if row.get("pipeline_ingest_status") not in INGEST:
            errors.append(f"{ident}: invalid pipeline_ingest_status")
        try:
            lo = float(row["wavelength_min_nm"])
            hi = float(row["wavelength_max_nm"])
            if not 0 < lo < hi:
                errors.append(f"{ident}: invalid wavelength range {lo}..{hi}")
        except (KeyError, TypeError, ValueError):
            errors.append(f"{ident}: non-numeric wavelength range")
        for field in ("resolving_power_min", "resolving_power_max"):
            try:
                if row.get(field) and float(row[field]) <= 0:
                    errors.append(f"{ident}: {field} must be positive")
            except ValueError:
                errors.append(f"{ident}: {field} must be numeric or blank")
        if row.get("codex_status") in {"active", "experimental", "reference"}:
            url = row.get("archive_url", "")
            if urlparse(url).scheme != "https":
                errors.append(f"{ident}: public non-candidate source needs HTTPS archive_url")
            if not row.get("source_documentation"):
                errors.append(f"{ident}: supported source lacks documentation")
        if row.get("telluric_required") == "yes" and not row.get("reduction_requirements"):
            errors.append(f"{ident}: telluric-heavy product lacks reduction guidance")
    return errors


README_SELECTION = (
    "hst_stis", "hst_cos", "harps", "espresso", "uves", "feros",
    "crires_plus", "nirps", "spirou", "ghost", "kpno_solar_atlas",
    "calspec_solar",
)


def render_markdown_table(path: Path = CATALOG) -> str:
    by_id = {row["instrument_id"]: row for row in load_catalog(path)}
    lines = [
        "| Instrument / atlas | Facility · archive | Coverage (nm) | R | Bands | Science role | Key caveats | Codex status |",
        "|---|---|---:|---:|---|---|---|---|",
    ]
    for ident in README_SELECTION:
        row = by_id[ident]
        rlo = row["resolving_power_min"]
        rhi = row["resolving_power_max"]
        resolving = rlo if rlo == rhi else f"{rlo}–{rhi}"
        caveat = row["reduction_requirements"].replace("|", "/")
        lines.append(
            f"| {row['instrument_name']} | {row['facility']} · "
            f"[{row['archive_name']}]({row['archive_url']}) | "
            f"{row['wavelength_min_nm']}–{row['wavelength_max_nm']} | "
            f"{resolving} | {row['bands_supported'].replace('|', '/')} | "
            f"{row['science_role']} | {caveat} | {row['codex_status']} |"
        )
    return "\n".join(lines)


def main() -> int:
    errors = validate_catalog()
    if errors:
        print("\n".join(f"[FAIL] {error}" for error in errors))
        return 1
    rows = load_catalog()
    print(f"[OK] {len(rows)} instrument/data-source records validated")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
