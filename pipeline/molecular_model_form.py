"""
RYA-1220 -- the molecular model-form (1D->3D) term, measured rather than held.

WHY THIS EXISTS. Every molecular CNO route in the gate matrix carried
`3d_nlte_model: HOLD_MOLECULAR`, and the stated reason was that no transition-specific
3D calculation had been run. But Amarsi et al. 2021 Table 2 publishes each of its 408
used lines under FIVE model atmospheres, so the 1D->3D shift is already measured for
every molecule in the reference set. Holding a component that the repo can price is the
RYA-1226 failure mode: "unpriceable" with no cited reason re-buries measured work.

WHAT THE TERM IS. Our molecular syntheses are 1D. The published 3D-minus-MARCS shift is
therefore not a correction we have applied -- it is the size of the systematic we are
carrying by staying 1D. So it enters as a BOUND on `model_atmosphere`, not as an
abundance move. Nothing here changes any A value (RYA-161 firewall).

WHY IT IS PER MOLECULE. The shift runs +0.035 (CH) to -0.140 (CO). A single uniform
"molecular 3D offset" would be wrong by up to 0.18 dex, which is precisely why no one
could write one down and the component stayed held.

WHY BAND IDENTITY GATES IT. The shift is measured on the transitions the reference
analysis actually used. Two of our diagnostics sit outside those ranges entirely --
CN_red at 6125-6200 A against a used CN range of 10875-13208 A, and NH_AX at 3358-3373 A
against a used NH range of 2.9-15 um, a band that analysis explicitly declined to retain.
Borrowing a CN A-X (0-0) NIR model term for the CN red system would be the
neighbouring-pool error: a real measurement of the wrong population. Those routes get
HOLD with that reason attached, not a silent number.
"""
from __future__ import annotations

import functools
import pathlib

import pandas as pd

from pipeline.error_budget import Term

REFERENCE_CSV = (pathlib.Path(__file__).resolve().parents[1] / "data" / "reference"
                 / "molecular_cno_literature_rya1220"
                 / "amarsi2021_molecular_reference_by_species.csv")

SOURCE = ("Amarsi et al. 2021, A&A 656 A113, Table 2 (CDS J/A+A/656/A113), 408 used "
          "lines under 3D / mean-3D / ATMO / MARCS / HM74; ingested at "
          "data/reference/amarsi2021_cno/")

#: Fractional slack when asking whether our window lies inside the used-line range.
BAND_TOLERANCE = 0.02


class ModelFormError(ValueError):
    """The term was requested for something the reference set cannot price."""


@functools.lru_cache(maxsize=1)
def reference_table() -> pd.DataFrame:
    if not REFERENCE_CSV.exists():
        raise ModelFormError(
            f"reference table absent: {REFERENCE_CSV}. Regenerate with "
            f"scripts/rya1220_molecular_literature_matrix.py")
    return pd.read_csv(REFERENCE_CSV).set_index("molecule")


#: Isotopologue spellings map onto the reference table's species names. The reference
#: analysis publishes one CN entry, not one per isotopologue, so 13C14N resolves to the
#: same row -- an isotope term is a SEPARATE component (hfs_isotopes), not this one.
_SPECIES_ALIAS = {
    "CO": "12C16O", "12C16O": "12C16O", "13C16O": "12C16O",
    "CN": "CN", "12C14N": "CN", "13C14N": "CN", "12C15N": "CN",
    "CH": "CH", "12CH": "CH", "13CH": "CH",
    "C2": "C2", "12C12C": "C2", "12C13C": "C2",
    "NH": "NH", "14NH": "NH", "15NH": "NH",
    "OH": "OH", "16OH": "OH", "18OH": "OH",
}


def _row(molecule: str) -> pd.Series:
    table = reference_table()
    key = _SPECIES_ALIAS.get(molecule, molecule)
    if key not in table.index:
        raise ModelFormError(
            f"{molecule}: not in the reference set {sorted(table.index)}. The 1D->3D "
            f"shift is molecule-specific; there is no universal molecular value to fall "
            f"back on.")
    return table.loc[key]


def reference_range_A(molecule: str) -> tuple[float, float]:
    """The air-Angstrom span of the transitions the reference analysis actually used."""
    r = _row(molecule)
    return float(r.lambda_vac_nm_min) * 10.0, float(r.lambda_vac_nm_max) * 10.0


def in_reference_band(molecule: str, lo_A: float, hi_A: float) -> bool:
    lo, hi = reference_range_A(molecule)
    return lo_A >= lo * (1 - BAND_TOLERANCE) and hi_A <= hi * (1 + BAND_TOLERANCE)


def model_form_term(molecule: str, lo_A: float, hi_A: float) -> dict:
    """The `model_atmosphere` component for a 1D molecular route.

    Returns a MEASURED bound when the route sits on the reference transitions, and a
    HOLD carrying its reason when it does not. Never returns a borrowed number.
    """
    r = _row(molecule)
    lo, hi = reference_range_A(molecule)
    shift = float(r.model_form_3D_minus_MARCS)
    common = {
        "component": "model_atmosphere",
        "molecule": molecule,
        "route_window_air_A": [float(lo_A), float(hi_A)],
        "reference_range_air_A": [round(lo, 1), round(hi, 1)],
        "reference_system": str(r.system),
        "n_reference_lines": int(r.n_lines),
        "published_3D_minus_MARCS_dex": shift,
        "published_3D_minus_mean3D_dex": float(r.granulation_3D_minus_mean3D),
        "published_line_to_line_sd_dex": float(r.line_to_line_sd_3D),
        "source": SOURCE,
    }
    if not in_reference_band(molecule, lo_A, hi_A):
        return {**common, "state": "HOLD", "sigma_dex": None,
                "reason": (
                    f"route window {lo_A:.0f}-{hi_A:.0f} A lies outside the reference "
                    f"used-line range {lo:.0f}-{hi:.0f} A for {molecule} "
                    f"{r.system}. The published 1D->3D shift was measured on those "
                    f"transitions; applying it here would be a real measurement of the "
                    f"wrong population. Measure this band or retire the diagnostic.")}
    return {**common, "state": "MEASURED", "sigma_dex": abs(shift),
            "reason": (
                f"our synthesis is 1D (MARCS); the published 1D->3D shift for {molecule} "
                f"{r.system} on these transitions is {shift:+.3f} dex, which bounds the "
                f"model-form systematic we carry by staying 1D. Bound, not a correction: "
                f"no abundance is moved.")}


def as_term(molecule: str, lo_A: float, hi_A: float) -> Term:
    """The MEASURED case as an error-budget Term; HOLD raises rather than returning zero."""
    got = model_form_term(molecule, lo_A, hi_A)
    if got["state"] != "MEASURED":
        raise ModelFormError(f"{molecule} {lo_A:.0f}-{hi_A:.0f} A: {got['reason']}")
    return Term(name="model_atmosphere", sigma=got["sigma_dex"], source=got["source"])
