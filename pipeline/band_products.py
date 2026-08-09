"""Per-(element × instrument × band × treatment) products — RYA-712/713.

Ryan, 2026-08-09: *"that median is calculating both engines, while cool to mention, is
not the science product we showcase. We showcase Engine A in one plot/product in the IR,
and the same treatment for B, and LTE."*

WHAT THIS MODULE IS FOR
-----------------------
A **product** is one element, measured on one instrument, in one band, under ONE
treatment, carrying its own value, its own sigma and its own line count. 1D-LTE is a
product. Engine A is a product. Engine B is a product. They are reported side by side
and never merged — RYA-712.

A cross-engine difference (B − A) is a **diagnostic**. It may be mentioned. It is not a
product and this module will not emit one as a headline: `Product` has no field for it,
and `combine()` does not exist.

WHY IT IS PARAMETERISED
-----------------------
Ryan: *"no bad copy paste of code lol, Ba kicked our ass last time."* Adapting the Ba
harness by hand for Al produced 13 defects, two of which (`'element': 'Ba'` and
`SOLAR_ASPLUND2021['Ba']`) would have reached an emitted record on an Al result.

So the element is a PARAMETER on every function here — there is no element symbol
anywhere in the logic, and `assert_single_element()` re-checks the emitted rows against
the element that was asked for. That guard is generic: it does not know or care which
element burned us last time, so it keeps working for the next one.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field, asdict
from pathlib import Path

import numpy as np
import pandas as pd

from pipeline._numcompat import trapezoid as _trapezoid  # numpy>=2 removed np.trapz

# Treatments are separate products. This is the whole vocabulary; there is deliberately
# no 'combined' member, because a combined value is not a product (RYA-712).
TREATMENTS = ("1D-LTE", "ENGINE-A", "ENGINE-B")

# Saturation: above this REW the EW->abundance inversion runs along the flat part of the
# curve of growth and is ill-conditioned in BOTH directions. Lines past it are measured
# and reported, never silently dropped -- but they are excluded from the aggregate and
# say why. (RYA-711: quarantined, not culled.)
REW_SATURATION_CEILING = -4.9

_ELEMENT_TOKEN = re.compile(r"\b([A-Z][a-z]?)\s*(?:I{1,3}|IV|VI?|1|2|3)?\b")


@dataclass
class LineMeasurement:
    """One line, one instrument, one treatment. Always carries WHY it was excluded."""
    element: str
    ion: str
    wavelength_air_A: float
    instrument: str
    ew_mA: float
    ew_method: str
    abundance: float | None = None
    rew: float | None = None
    treatment: str = ""
    in_aggregate: bool = True
    excluded_reason: str = ""

    def __post_init__(self) -> None:
        if self.rew is None and self.ew_mA and self.wavelength_air_A:
            self.rew = float(np.log10(self.ew_mA * 1e-3 / self.wavelength_air_A))
        if self.rew is not None and self.rew > REW_SATURATION_CEILING:
            self.in_aggregate = False
            self.excluded_reason = (
                f"REW {self.rew:.3f} above the {REW_SATURATION_CEILING} saturation "
                f"ceiling — the EW->abundance inversion is ill-conditioned here. "
                f"Measured and reported; excluded from the aggregate only.")


@dataclass
class Product:
    """One showcased product. Deliberately has NO field for a cross-engine delta."""
    element: str
    ion: str
    instrument: str
    band: str
    treatment: str
    value: float | None
    sigma: float | None
    n_lines: int
    n_excluded: int
    lines: list[LineMeasurement] = field(default_factory=list)
    provenance: str = ""

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame([asdict(l) for l in self.lines])


def assert_single_element(rows, element: str) -> None:
    """Refuse to emit a record carrying a DIFFERENT element than the one requested.

    Generic on purpose. The Ba->Al incident left `'element': 'Ba'` on an Al result; a
    guard hard-coded to look for 'Ba' would have caught that once and then rotted. This
    compares every emitted element field against the element actually asked for, so it
    catches the next copy-paste too, whichever pair it involves.
    """
    seen = set()
    for r in rows:
        e = r.element if isinstance(r, LineMeasurement) else r.get("element")
        if e:
            seen.add(str(e).strip())
    foreign = seen - {element}
    if foreign:
        raise ValueError(
            f"emitted rows claim element(s) {sorted(foreign)} but this product is for "
            f"{element!r}. This is the copy-paste signature (RYA-701) — a harness "
            f"adapted from another element kept its source's identity. Fix the harness; "
            f"do not filter the rows.")


def assert_no_cross_treatment_mix(lines, treatment: str) -> None:
    """A product is ONE treatment. Mixing is the RYA-712 violation."""
    seen = {l.treatment for l in lines if l.treatment}
    if seen - {treatment}:
        raise ValueError(
            f"product declares treatment {treatment!r} but its lines carry {sorted(seen)}. "
            f"Engines and LTE are separate products and are never combined (RYA-712).")


def local_continuum(w: np.ndarray, f: np.ndarray, centre: float, half_width: float,
                    *, sideband_frac: float = 0.5, pct: float = 95.0) -> float:
    """Continuum from clean side-bands, as the 95th percentile of the flanking regions.

    Side-bands sit OUTSIDE the fitting window so line flux never sets its own continuum.
    A percentile rather than a max, so one hot pixel cannot lift the whole normalisation.
    """
    inner, outer = half_width, half_width * (1.0 + sideband_frac * 2.0)
    d = np.abs(w - centre)
    band = (d > inner) & (d <= outer)
    if band.sum() < 5:
        raise ValueError(f"too few side-band points around {centre:.3f} "
                         f"({int(band.sum())}) — cannot set a local continuum honestly")
    return float(np.percentile(f[band], pct))


def equivalent_width(w: np.ndarray, f: np.ndarray, centre: float, half_width: float,
                     **kw) -> tuple[float, str]:
    """EW in mA over +/-half_width, normalised to a local continuum. Returns the method
    string too, because how the continuum was set IS part of the measurement."""
    cont = local_continuum(w, f, centre, half_width, **kw)
    m = np.abs(w - centre) <= half_width
    if m.sum() < 3:
        raise ValueError(f"too few points inside the window at {centre:.3f}")
    depth = 1.0 - (f[m] / cont)
    ew = float(_trapezoid(depth, w[m]) * 1000.0)
    return ew, (f"integrated over +/-{half_width:.3f} A, local continuum = "
                f"{pct_label(kw)} of flanking side-bands (cont={cont:.5f})")


def pct_label(kw) -> str:
    return f"{kw.get('pct', 95.0):.0f}th pct"


def build_product(element: str, ion: str, instrument: str, band: str, treatment: str,
                  lines: list[LineMeasurement], *, provenance: str = "") -> Product:
    """Aggregate ONE treatment into ONE product. Never touches another treatment."""
    if treatment not in TREATMENTS:
        raise ValueError(f"unknown treatment {treatment!r}; expected one of {TREATMENTS}")
    assert_single_element(lines, element)
    assert_no_cross_treatment_mix(lines, treatment)

    used = [l for l in lines if l.in_aggregate and l.abundance is not None]
    excluded = [l for l in lines if not l.in_aggregate]
    for l in excluded:
        if not l.excluded_reason.strip():
            raise ValueError(
                f"{element} {l.wavelength_air_A:.3f} is excluded with no reason recorded. "
                f"An unexplained exclusion is indistinguishable from a silent drop "
                f"(RYA-429).")

    vals = np.array([l.abundance for l in used], dtype=float)
    value = float(np.median(vals)) if len(vals) else None
    # Sigma is the scatter of the lines that actually entered, not a fitted error bar.
    sigma = float(np.std(vals, ddof=1)) if len(vals) > 1 else None
    return Product(element=element, ion=ion, instrument=instrument, band=band,
                   treatment=treatment, value=value, sigma=sigma,
                   n_lines=len(used), n_excluded=len(excluded),
                   lines=lines, provenance=provenance)


def products_frame(products: list[Product]) -> pd.DataFrame:
    """Side-by-side table of products. One row per product — NOT a merge.

    Any comparison a reader makes between rows is theirs to make and is a diagnostic.
    This function will not compute it for them.
    """
    return pd.DataFrame([
        dict(element=p.element, ion=p.ion, instrument=p.instrument, band=p.band,
             treatment=p.treatment, value=p.value, sigma=p.sigma,
             n_lines=p.n_lines, n_excluded=p.n_excluded, provenance=p.provenance)
        for p in products])
