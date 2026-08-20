"""Minimal, audited reader for NOAA ARL-packed GDAS1 pressure profiles (RYA-927).

GDAS1 is a 360 x 181 global 1-degree grid in 3-hour slots.  The decoder follows
NOAA ARL's PAKOUT inverse exactly; it is intentionally narrow and reads only the
HGTS/TEMP/RELH pressure-level fields and surface PRSS/T02M/RH2M needed by molecfit.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np

NX, NY = 360, 181
RECORD_LEN = 50 + NX * NY
RECORDS_PER_SLOT = 164
PRESSURE_HPA = (1000, 975, 950, 925, 900, 850, 800, 750, 700, 650, 600,
                550, 500, 450, 400, 350, 300, 250, 200, 150, 100, 50, 20)


@dataclass(frozen=True)
class ArlHeader:
    year: int
    month: int
    day: int
    hour: int
    level: int
    variable: str
    exponent: int
    precision: float
    first_value: float


def parse_header(raw: bytes) -> ArlHeader:
    if len(raw) != 50:
        raise ValueError(f"ARL header must be 50 bytes, got {len(raw)}")
    text = raw.decode("ascii")
    return ArlHeader(
        year=2000 + int(text[0:2]), month=int(text[2:4]), day=int(text[4:6]),
        hour=int(text[6:8]), level=int(text[10:12]), variable=text[14:18].strip(),
        exponent=int(text[18:22]), precision=float(text[22:36]),
        first_value=float(text[36:50]),
    )


def unpack(raw: bytes, header: ArlHeader) -> np.ndarray:
    """Inverse of NOAA PAKOUT, preserving its per-row first-column predictor."""
    if len(raw) != NX * NY:
        raise ValueError(f"ARL field must be {NX * NY} bytes, got {len(raw)}")
    packed = np.frombuffer(raw, dtype=np.uint8).reshape(NY, NX)
    scale = 2.0 ** (7 - header.exponent)
    out = np.empty((NY, NX), dtype=np.float64)
    first_column = header.first_value
    for j in range(NY):
        previous = first_column
        for i in range(NX):
            previous += (int(packed[j, i]) - 127) / scale
            out[j, i] = previous
            if i == 0:
                first_column = previous
    return out


def _nearest_grid(lat: float, lon: float) -> tuple[int, int, float, float]:
    lon360 = lon % 360.0
    i = int(np.floor(lon360 + 0.5)) % NX
    j = int(np.clip(np.floor(lat + 90.0 + 0.5), 0, NY - 1))
    return i, j, float(j - 90), float(i if i <= 180 else i - 360)


def extract_profile(path, when: datetime, lat: float, lon: float) -> dict:
    """Extract one exact 3-hour GDAS slot at the nearest grid point."""
    if when.minute or when.second or when.microsecond or when.hour % 3:
        raise ValueError("GDAS extraction time must be an exact 3-hour UTC slot")
    i, j, grid_lat, grid_lon = _nearest_grid(lat, lon)
    wanted = {"PRSS", "T02M", "RH2M", "HGTS", "TEMP", "RELH"}
    fields: dict[tuple[int, str], float] = {}
    with Path(path).open("rb") as handle:
        size = Path(path).stat().st_size
        if size % RECORD_LEN:
            raise ValueError(f"ARL file size {size} is not divisible by {RECORD_LEN}")
        for recno in range(size // RECORD_LEN):
            header_raw = handle.read(50)
            data = handle.read(NX * NY)
            h = parse_header(header_raw)
            stamp = datetime(h.year, h.month, h.day, h.hour)
            if stamp != when or h.variable not in wanted:
                continue
            fields[(h.level, h.variable)] = float(unpack(data, h)[j, i])

    rows = []
    for level, pressure in enumerate(PRESSURE_HPA, start=1):
        keys = [(level, name) for name in ("HGTS", "TEMP")]
        if not all(key in fields for key in keys):
            raise ValueError(f"GDAS slot lacks pressure-level fields at {pressure} hPa")
        # GDAS1 omits RELH at 50 and 20 hPa.  ESO's own molecfit GDAS ASCII
        # profiles encode those same two absent stratospheric values as 0.0.
        humidity = fields.get((level, "RELH"), 0.0 if pressure <= 50 else None)
        if humidity is None:
            raise ValueError(f"GDAS slot lacks RELH at {pressure} hPa")
        rows.append({"press_hpa": float(pressure), "height_m": fields[(level, "HGTS")],
                     "temp_k": fields[(level, "TEMP")], "relhum_pct": humidity})
    for key in ((0, "PRSS"), (0, "T02M"), (0, "RH2M")):
        if key not in fields:
            raise ValueError(f"GDAS slot lacks surface field {key[1]}")
    return {
        "slot_utc": when.isoformat(), "requested_lat": lat, "requested_lon": lon,
        "grid_lat": grid_lat, "grid_lon": grid_lon,
        # GDAS1 PRSS is hPa in this ARL product.
        "surface_pressure_hpa": fields[(0, "PRSS")],
        "surface_temp_k": fields[(0, "T02M")],
        "surface_relhum_pct": fields[(0, "RH2M")], "levels": rows,
    }


def write_molecfit_ascii(profile: dict, path) -> Path:
    """Write the four-column format accepted by gdas_fetch._gdas_ascii_to_fits."""
    path = Path(path)
    lines = ["# P[hPa] HGT[m] T[K] RELHUM[%]"]
    lines += [f"{r['press_hpa']:7.1f} {r['height_m']:10.3f} {r['temp_k']:8.3f} {r['relhum_pct']:8.3f}"
              for r in profile["levels"]]
    path.write_text("\n".join(lines) + "\n")
    return path


def interpolate_at_height(profile: dict, height_m: float) -> dict:
    """Interpolate GDAS P/T/RH to observatory elevation for telemetry validation."""
    rows = profile["levels"]
    height = np.array([r["height_m"] for r in rows], dtype=float)
    if not height.min() <= height_m <= height.max():
        raise ValueError(f"height {height_m} m lies outside GDAS profile")
    pressure = np.array([r["press_hpa"] for r in rows], dtype=float)
    temp = np.array([r["temp_k"] for r in rows], dtype=float)
    rh = np.array([r["relhum_pct"] for r in rows], dtype=float)
    return {
        "height_m": float(height_m),
        # Pressure is exponential with height; interpolate log(P), not P.
        "pressure_hpa": float(np.exp(np.interp(height_m, height, np.log(pressure)))),
        "temperature_c": float(np.interp(height_m, height, temp) - 273.15),
        "relative_humidity_pct": float(np.interp(height_m, height, rh)),
    }
