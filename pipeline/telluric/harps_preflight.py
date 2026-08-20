"""Exposure-level preflight for HARPS molecfit correction (RYA-927).

This module deliberately does not correct a stacked/normalised CSV.  It inventories
the source Phase-3 FITS products, records the observation and weather telemetry, and
resolves the real per-exposure La Silla GDAS profile.  A correction driver may consume
the resulting records only after every exposure passes this gate.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

from astropy.io import fits

from pipeline.telluric.gdas_fetch import fetch_gdas, nearest_3hourly


@dataclass(frozen=True)
class HarpsExposure:
    path: str
    date_obs: str
    mjd_obs: float
    program_id: str
    object_name: str
    airmass_start: float
    airmass_end: float
    pressure_hpa: float
    temperature_c: float
    relative_humidity_pct: float
    gdas_slot_utc: str
    gdas_profile: str

    def manifest_record(self) -> dict:
        return asdict(self)


def _required(header, key: str):
    value = header.get(key)
    if value is None or str(value).strip() == "":
        raise ValueError(f"HARPS exposure lacks required provenance/telemetry key {key!r}")
    return value


def inspect_exposure(path, *, cache_dir=None) -> HarpsExposure:
    """Read one Phase-3 HARPS exposure and resolve its observation-time GDAS profile."""
    path = Path(path)
    with fits.open(path, mode="readonly", memmap=True) as hdul:
        h = hdul[0].header
        instrument = str(_required(h, "INSTRUME")).strip().upper()
        if instrument != "HARPS":
            raise ValueError(f"{path}: expected INSTRUME=HARPS, got {instrument!r}")
        mjd = float(_required(h, "MJD-OBS"))
        program = str(_required(h, "ESO OBS PROG ID"))
        if program != "1102.D-0954(A)":
            raise ValueError(f"{path}: not the RYA-927 solar sequence ({program!r})")
        date_obs = str(_required(h, "DATE-OBS"))
        obj = str(_required(h, "OBJECT"))
        airmass_start = float(_required(h, "ESO TEL AIRM START"))
        airmass_end = float(_required(h, "ESO TEL AIRM END"))
        pressure = float(_required(h, "ESO TEL AMBI PRES START"))
        temperature = float(_required(h, "ESO TEL AMBI TEMP"))
        humidity = float(_required(h, "ESO TEL AMBI RHUM"))

    profile = fetch_gdas("la_silla", mjd=mjd, cache_dir=cache_dir)
    slot = nearest_3hourly(None, mjd=mjd)
    return HarpsExposure(
        path=str(path), date_obs=date_obs, mjd_obs=mjd, program_id=program,
        object_name=obj, airmass_start=airmass_start, airmass_end=airmass_end,
        pressure_hpa=pressure, temperature_c=temperature,
        relative_humidity_pct=humidity, gdas_slot_utc=slot.isoformat(),
        gdas_profile=str(profile),
    )


def inspect_sequence(source_dir, *, cache_dir=None) -> list[HarpsExposure]:
    """Preflight the canonical ten-exposure solar sequence; loud-fail on drift."""
    files = sorted(Path(source_dir).glob("*.fits"))
    if len(files) != 10:
        raise ValueError(f"expected 10 direct-solar HARPS FITS products, found {len(files)}")
    records = [inspect_exposure(p, cache_dir=cache_dir) for p in files]
    if len({r.gdas_slot_utc for r in records}) != 1:
        raise ValueError("solar sequence spans multiple GDAS slots; correct per slot")
    return records
