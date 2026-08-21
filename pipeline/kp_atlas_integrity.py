"""Structural integrity of the 1984 Kitt Peak Flux Atlas (RYA-938).

The atlas is 251 plain-text segments, and `measure_band_ew.kp_segments()` used to
inventory them behind a bare ``except Exception: continue``.  That is the
[RYA-833] shape: a segment that fails to parse is dropped from the inventory,
and the very next question -- "does anything cover 8420 A?" -- answers "no Kitt
Peak segment covers 8420.000 A".  A corrupt file therefore presents as MISSING
COVERAGE, which is indistinguishable from genuinely not holding the data.

That is not hypothetical.  `lm0840` in both staged copies is a saved HTTP 500
error page (714 bytes), as is the bundled `README`; the download failed and the
failure was archived as if it were data.

The invariant pinned here is what a segment IS, not which file happened to be
broken: three numeric columns, strictly increasing wavelength, and a residual
flux that lives near unity because unity IS the continuum for this product
(NSO Atlas No. 1 README: "the first column contains the wavelength in air, the
second column contains the pseudo-residual flux").
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np

#: NSO Atlas No. 1: each file carries 4 nm plus 0.1 nm of overlap at the red end.
SEGMENT_SPAN_NM = 4.0
SEGMENT_OVERLAP_NM = 0.1
EXPECTED_COLUMNS = 3
#: A pseudo-residual flux may exceed 1 slightly (noise about the continuum) but
#: cannot be negative or wildly super-unity. Bounds are generous on purpose: this
#: separates DATA from NOT-DATA, it does not judge data quality.
FLUX_LO, FLUX_HI = -0.2, 2.0


class KpAtlasCorrupt(RuntimeError):
    """A segment exists but is not readable as atlas data."""


@dataclass(frozen=True)
class SegmentReport:
    path: str
    name: str
    size_bytes: int
    ok: bool
    reason: str = ""
    n_rows: int = 0
    n_columns: int = 0
    lo_A: float = float("nan")
    hi_A: float = float("nan")
    monotonic: bool = False
    flux_min: float = float("nan")
    flux_max: float = float("nan")
    stem_start_nm: float = float("nan")
    stem_matches_data: bool = False

    def as_dict(self) -> dict:
        return asdict(self)


def inspect_segment(path: Path) -> SegmentReport:
    """Read one `lmNNNN` segment and say -- with a reason -- whether it is data."""
    path = Path(path)
    size = path.stat().st_size if path.exists() else 0
    stem_nm = float("nan")
    if path.name.startswith("lm") and path.name[2:].isdigit():
        stem_nm = float(path.name[2:])

    try:
        data = np.loadtxt(path)
    except Exception as exc:                       # noqa: BLE001 - reason is the product
        head = ""
        try:
            head = path.read_text(errors="replace")[:80].replace("\n", " ").strip()
        except OSError:
            pass
        reason = f"not parseable as numeric text ({type(exc).__name__})"
        if "<html" in head.lower() or "<!doctype" in head.lower():
            reason = ("the file is saved HTML, not spectral data -- a failed download "
                      f"was archived as if it were the segment: {head[:60]!r}")
        return SegmentReport(str(path), path.name, size, False, reason,
                             stem_start_nm=stem_nm)

    data = np.atleast_2d(data)
    if data.ndim != 2 or data.shape[1] != EXPECTED_COLUMNS:
        return SegmentReport(str(path), path.name, size, False,
                             f"expected {EXPECTED_COLUMNS} columns, got shape {data.shape}",
                             n_rows=int(data.shape[0]), stem_start_nm=stem_nm)

    wave_nm, flux = data[:, 0], data[:, 1]
    monotonic = bool(np.all(np.diff(wave_nm) > 0))
    lo_a, hi_a = float(wave_nm[0] * 10.0), float(wave_nm[-1] * 10.0)
    stem_ok = bool(np.isnan(stem_nm) or abs(wave_nm[0] - stem_nm) < 0.05)

    reason = ""
    ok = True
    if not monotonic:
        ok, reason = False, "wavelength column is not strictly increasing"
    elif not (FLUX_LO <= float(np.nanmin(flux)) and float(np.nanmax(flux)) <= FLUX_HI):
        ok, reason = False, (f"residual flux outside [{FLUX_LO}, {FLUX_HI}]: "
                             f"{float(np.nanmin(flux)):.4g}..{float(np.nanmax(flux)):.4g}")
    elif not stem_ok:
        ok, reason = False, (f"filename says {stem_nm} nm but data starts at "
                             f"{float(wave_nm[0]):.5f} nm")

    return SegmentReport(str(path), path.name, size, ok, reason,
                         n_rows=int(data.shape[0]), n_columns=int(data.shape[1]),
                         lo_A=lo_a, hi_A=hi_a, monotonic=monotonic,
                         flux_min=float(np.nanmin(flux)), flux_max=float(np.nanmax(flux)),
                         stem_start_nm=stem_nm, stem_matches_data=stem_ok)


def inspect_atlas(directory: Path) -> dict:
    """Inventory every segment, and report coverage gaps between good ones."""
    directory = Path(directory)
    reports = [inspect_segment(p) for p in sorted(directory.glob("lm[0-9]*"))
               if p.is_file()]
    good = sorted([r for r in reports if r.ok], key=lambda r: r.lo_A)
    gaps = []
    for a, b in zip(good, good[1:]):
        if b.lo_A > a.hi_A:
            gaps.append({"after": a.name, "before": b.name,
                         "gap_lo_A": a.hi_A, "gap_hi_A": b.lo_A,
                         "gap_width_A": b.lo_A - a.hi_A})
    corrupt = [r for r in reports if not r.ok]
    return {
        "directory": str(directory),
        "n_files": len(reports),
        "n_ok": len(good),
        "n_corrupt": len(corrupt),
        "corrupt": [r.as_dict() for r in corrupt],
        "coverage_lo_A": good[0].lo_A if good else None,
        "coverage_hi_A": good[-1].hi_A if good else None,
        "n_gaps": len(gaps),
        "gaps": gaps,
        "segments": [r.as_dict() for r in reports],
    }


def require_parseable(reports) -> None:
    """Raise naming every corrupt segment, rather than letting one go quiet.

    Callers that genuinely want degraded operation must ask for it explicitly and
    record what they skipped -- an omission nobody chose is the thing this
    prevents.
    """
    bad = [r for r in reports if not r.ok]
    if not bad:
        return
    lines = [f"{len(bad)} Kitt Peak atlas segment(s) exist but are NOT readable as data.",
             "",
             "This is NOT a coverage gap. Each of these files is present, so a loader that "
             "skips it silently reports the wavelengths it covers as data we do not hold "
             "(RYA-833).",
             ""]
    for r in bad:
        lines.append(f"  {r.name}  ({r.size_bytes} bytes): {r.reason}")
    lines += ["", "Re-fetch the segment (NSO Atlas No. 1, "
                  "https://nispdata.nso.edu/ftp/pub/atlas/fluxatl/) and verify it, or "
                  "register the gap explicitly. Do not leave it silent."]
    raise KpAtlasCorrupt("\n".join(lines))
