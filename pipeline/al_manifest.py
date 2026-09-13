"""RYA-1176 validation and loading for the frozen Solar Al manifest."""
from __future__ import annotations
from pathlib import Path
import pandas as pd
from pipeline.model_registry import LINE_SETS

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data/audit/rya1176_al_manifest/al_line_manifest_v2.csv"
REQUIRED_COLUMNS = ("line_set", "gf_provenance", "selection_state", "telluric_applied", "normalization_state", "observed_conditioning")

def validate_manifest(frame: pd.DataFrame) -> list[str]:
    errors = []
    missing = [c for c in REQUIRED_COLUMNS if c not in frame.columns]
    if missing:
        return ["missing required columns: " + ",".join(missing)]
    bad = sorted(set(frame["line_set"].fillna("").astype(str)) - set(LINE_SETS))
    if bad:
        errors.append(f"line_set values outside model_registry.LINE_SETS: {bad}")
    for column in ("telluric_applied", "normalization_state", "observed_conditioning"):
        if frame[column].fillna("").eq("").any():
            errors.append(f"{column} contains blank values")
    if frame["line_set"].fillna("").eq("").any():
        errors.append("line_set contains unresolved blank values")
    return errors

def load_manifest(path: Path = MANIFEST, *, require_conditioning: bool = False) -> pd.DataFrame:
    frame = pd.read_csv(path).fillna("")
    errors = validate_manifest(frame)
    if errors:
        raise ValueError("RYA-1176 manifest refused: " + " | ".join(errors))
    if require_conditioning and frame[["telluric_applied", "normalization_state", "observed_conditioning"]].eq("unknown").any(axis=1).any():
        raise ValueError("RYA-1176 product refused: holding conditioning is unresolved")
    return frame

def require_product_manifest(frame: pd.DataFrame, *, holding: dict | None = None) -> pd.DataFrame:
    if holding is None:
        raise ValueError("RYA-1176 product refused: no holding supplied")
    out = frame.copy()
    for column in ("telluric_applied", "normalization_state", "observed_conditioning"):
        value = str(holding.get(column, "unknown") or "unknown").strip()
        if value == "unknown":
            raise ValueError(f"RYA-1176 product refused: {column}=unknown")
        out[column] = value
    errors = validate_manifest(out)
    if errors:
        raise ValueError("RYA-1176 product refused: " + " | ".join(errors))
    return out
