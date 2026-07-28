"""Loader and validation contract for the canonical instrument catalog."""

from __future__ import annotations

import csv
from pathlib import Path
from urllib.parse import urlparse

CATALOG = Path(__file__).with_name("instrument_catalog.csv")
MODES = Path(__file__).with_name("instrument_modes.csv")
HOLDINGS = Path(__file__).with_name("holdings_manifest_registry.csv")
SYSTEMS = Path(__file__).with_name("system_catalog.csv")
STATUSES = {"active", "experimental", "candidate", "rejected", "reference"}
INGEST = {"loader_ready", "loader_partial", "audit_only", "candidate",
          "reference_only", "rejected"}
BANDS = {"FUV", "NUV", "VIS", "red_optical", "NIR", "Y", "J", "H", "K", "L", "M", "IR",
         "FUV_edge"}
MODE_STATUSES = {"active", "experimental", "candidate", "rejected"}
EVIDENCE_STATES = {"verified", "audited", "candidate", "rejected"}
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


def _load_csv(path: Path) -> list[dict[str, str]]:
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
        bands = set(filter(None, row.get("bands_supported", "").split("|")))
        if not bands or not bands <= BANDS:
            errors.append(f"{ident}: invalid bands_supported {sorted(bands - BANDS)}")
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


def validate_modes(path: Path = MODES) -> list[str]:
    rows = _load_csv(path)
    instrument_ids = {row["instrument_id"] for row in load_catalog()}
    errors: list[str] = []
    mode_ids: set[str] = set()
    for line, row in enumerate(rows, start=2):
        ident = row.get("mode_id", "")
        if not ident or ident in mode_ids:
            errors.append(f"mode line {line}: blank/duplicate mode_id {ident!r}")
        mode_ids.add(ident)
        if row.get("instrument_id") not in instrument_ids:
            errors.append(f"{ident}: unknown instrument_id {row.get('instrument_id')!r}")
        if row.get("mode_status") not in MODE_STATUSES:
            errors.append(f"{ident}: invalid mode_status")
        bands = set(filter(None, row.get("bands_supported", "").split("|")))
        if not bands or not bands <= BANDS:
            errors.append(f"{ident}: invalid bands_supported {sorted(bands - BANDS)}")
        try:
            lo, hi = float(row["wavelength_min_nm"]), float(row["wavelength_max_nm"])
            resolution = float(row["resolving_power"])
            if not 0 < lo < hi:
                errors.append(f"{ident}: invalid wavelength range")
            if resolution <= 0:
                errors.append(f"{ident}: resolving_power must be positive")
        except (KeyError, TypeError, ValueError):
            errors.append(f"{ident}: non-numeric wavelength/resolution")
    return errors


def validate_holdings(path: Path = HOLDINGS) -> list[str]:
    rows = _load_csv(path)
    instrument_ids = {row["instrument_id"] for row in load_catalog()}
    system_ids = {row["star_params_key"] for row in _load_csv(SYSTEMS)
                  if row["star_params_key"]}
    root = CATALOG.parents[2]
    errors: list[str] = []
    holding_ids: set[str] = set()
    for line, row in enumerate(rows, start=2):
        ident = row.get("holding_id", "")
        if not ident or ident in holding_ids:
            errors.append(f"holding line {line}: blank/duplicate holding_id {ident!r}")
        holding_ids.add(ident)
        if row.get("instrument_id") not in instrument_ids:
            errors.append(f"{ident}: unknown instrument_id {row.get('instrument_id')!r}")
        if row.get("system_id") not in system_ids:
            errors.append(f"{ident}: unknown system_id {row.get('system_id')!r}")
        manifest = row.get("manifest_path", "")
        if not manifest or Path(manifest).is_absolute() or not (root / manifest).is_file():
            errors.append(f"{ident}: missing/non-portable manifest_path {manifest!r}")
        if row.get("evidence_state") not in EVIDENCE_STATES:
            errors.append(f"{ident}: invalid evidence_state")
    return errors


def validate_all() -> list[str]:
    return validate_catalog() + validate_modes() + validate_holdings()


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
    errors = validate_all()
    if errors:
        print("\n".join(f"[FAIL] {error}" for error in errors))
        return 1
    rows = load_catalog()
    print(
        f"[OK] {len(rows)} sources, {len(_load_csv(MODES))} modes, and "
        f"{len(_load_csv(HOLDINGS))} holdings links validated"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
