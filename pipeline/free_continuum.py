"""The IR hyperfine special case: a per-line free-continuum fit, run BESIDE the standard one.

RYA-470 (Ryan, 2026-11-18): "it would fire this method IF the element has hyperfine BS going
on, and we need to resolve it in the IR, this would probably not be applicable for VIS" /
"both methods should be done in cases where it will provide results for both line by line
versus continuum, but continuum should only be done in special cases."

The standard synthesis fit holds the frozen per-band continuum. The special case
(`derive_band_products --free-continuum`) fits a continuum scale jointly with A(X) per line,
as Takeda (1995) and Elgueta et al. (2024) fit a "continuum normalization factor". On V I in
CRIRES+ H/K the standard fit rails at the bracket floor on every line, and the free fit
returns 3.56-4.14 on three lines (Fe I control: 7.476 vs 7.50 standard).

Everything that NAMES the special case lives here, so the descriptor, the cell matcher, the
artifact glob, the budget and coverage cannot spell it differently:
  * the cell carries `continuum="free"`;
  * its products publish under the standard selector + SELECTOR_SUFFIX, so a
    special-case product never takes the standard product's identity;
  * derive_band_products' stem gains STEM_TAG.
"""
from __future__ import annotations

import csv
from functools import lru_cache
from pathlib import Path

FREE = "free"
FLAG = "--free-continuum"
SELECTOR_SUFFIX = "-FREECONT"
STEM_TAG = "_FREECONT"
#: band_policy names that are infrared. VIS / red-optical / near-UV never trigger.
IR_BANDS = frozenset({"NIR", "H", "K"})

_CANONICAL = Path(__file__).resolve().parents[1] / "data" / "linelists" / "canonical_gf.csv"


@lru_cache(maxsize=None)
def _hfs_waves(species: str, path: str = str(_CANONICAL)) -> tuple[float, ...]:
    """Wavelengths of the canonical physical lines of `species` split into >1 HFS component."""
    out = []
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r.get("species") != species:
                continue
            try:
                if float(r.get("hfs_n_components") or 0) > 1:
                    out.append(float(r["wavelength_air_A"]))
            except ValueError:
                continue
    return tuple(sorted(out))


def applies(element: str, ion: str, band: str, lo_A: float, hi_A: float,
            canonical: Path | str = _CANONICAL) -> bool:
    """True when the special case fires: an IR band AND a line of this species in the window
    that the canonical list carries as a hyperfine multiplet (hfs_n_components > 1)."""
    if band not in IR_BANDS:
        return False
    return any(lo_A <= w <= hi_A for w in _hfs_waves(f"{element} {ion}", str(canonical)))


def selector(base: str | None, free: bool) -> str | None:
    """The published selector: the standard one, suffixed for the special case."""
    if not free:
        return base
    return (base or "ALL") + SELECTOR_SUFFIX


def base_selector(sel: str | None) -> str | None:
    """The standard selector a (possibly suffixed) selector was built from."""
    if sel and sel.endswith(SELECTOR_SUFFIX):
        return sel[: -len(SELECTOR_SUFFIX)] or None
    return sel


def is_free_stem(name: str) -> bool:
    return STEM_TAG + "_" in name or name.endswith(STEM_TAG)
