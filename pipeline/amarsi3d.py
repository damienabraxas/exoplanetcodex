"""Amarsi 2022 3D-NLTE Fe engine — reactivation + LINE-parameter domain check (RYA-817).

WHAT THIS MODULE ADDS
---------------------
`pipeline/nlte_corrections.py` already carries the Amarsi, Liljegren & Nissen 2022
(A&A 668, A68) 3D-NLTE MLP. It was archived in RYA-319 for single-methodology
consistency, not because it is wrong, and Ryan's 2026-08-14 standing rule
(`docs/SCIENCE_STANDARDS.md`) is that a capable model archived for convenience gets
used where it adds value. This module is the reactivation path: it calls the SAME
`_apply_aberr_to_line` the archived leg uses — no second copy of the sign convention,
no second copy of the iteration (RYA-701: one hand-adapted copy produced 13 defects) —
and adds the guard that was missing.

THE GUARD, AND WHY IT IS THE SCIENCE
------------------------------------
The vendored `main_aberr.py` checks only the STELLAR box: Teff 5000-6500 K, log g
4.0-4.5, vmic 0-3 km/s, A(Fe) 4.5-7.5. It does not check the LINE parameters at all. So
it will return a confident float for any (Elo, Eup, log gf) you hand it, including a
combination the network never saw. `_in_grid` in nlte_corrections.py inherits exactly
that blind spot.

RYA-817 runs the network on the near-IR Fe band. The training set (recovered by
`scripts/rya817_recover_amarsi_training_set.py`; the Jofre et al. 2014 'golden' list,
171 Fe I + 12 Fe II) is entirely OPTICAL: 4787.83-6810.26 A. So the run needs a
line-parameter domain check, and this module is it.

FOUR AXES, REPORTED SEPARATELY — never collapsed into one boolean
-----------------------------------------------------------------
  STELLAR   the existing Teff/logg/vmic/A(Fe) box (delegated to nlte_corrections).

            One subtlety, and it bit the first run of this ticket. The grid's fourth
            axis is A(Fe;3N) -- the 3D NON-LTE abundance, ceiling 7.5 -- and the vendor
            README is explicit that the user "should adopt an initial guess for the 3D
            non-LTE abundance and iterate". The solar per-line 1D LTE abundances run to
            ~7.8, so testing the INITIAL GUESS against that ceiling rejects every strong
            line and keeps the weak ones, which biased the first VIS control 0.07 dex
            LOW purely by selection. The axis test therefore applies to the CONVERGED
            A(Fe;3N) = A(1D-LTE) + aberr, which is the quantity the axis actually means.
  FEATURE   Elo, Eup and log gf each inside the min-max the assigned network saw.
  DELTA_E   Eup - Elo inside the training range.

            This axis is the one that matters and it is not obvious. Amarsi derives
            Eup as Elo + hc/lambda_vac -- verified to 4e-6 eV against the vendored
            test_data.csv -- so (Eup - Elo) IS the transition energy, i.e. the
            WAVELENGTH, entering the network as a derived feature. The premise that
            the MLP is wavelength-agnostic because it splits on excitation potential
            is true of its ROUTING and false of its INPUTS. A line can sit comfortably
            inside the Elo and Eup ranges and still be an extrapolation, because the
            PAIR encodes a photon energy the network never saw.

  LEVEL     both Elo and Eup match an energy that actually appears in the training
            set, within `LEVEL_TOL_EV`. RYA-817 asks for this explicitly ("check the
            line's levels are represented, not just the Elo scalar"), and RYA-763's
            lesson applies: an energy is a weak proxy for a level, so this axis is
            reported as evidence, never as the sole basis for admitting a line.

A line is IN-DOMAIN only if all four pass. Anything else is refused with the axis named.
There is no extrapolation switch here on purpose: `main_aberr.py` has one and its own
README says the authors cannot vouch for the result.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from pipeline.nlte_corrections import (  # the archived leg — reused, never re-typed
    _apply_aberr_to_line,
    _in_grid,
    _GRID,
)

_REPO_ROOT = Path(__file__).resolve().parents[1]
TRAINING_CSV = (_REPO_ROOT / "data" / "reference" / "amarsi2022_training"
                / "amarsi2022_training_lines.csv")

CITATION = "Amarsi, Liljegren & Nissen 2022, A&A 668, A68 (3D-NLTE Fe MLP)"
TRAINING_CITATION = ("Jofre et al. 2014, A&A 564, A133, Tables 4/5 "
                     "('golden' Fe I/II line list)")
SOLAR_CONTROL_CITATION = ("Allende Prieto et al. 2002, ApJ 567, 544, Table 2 "
                          "(the solar Fe line list Amarsi+2022 Table 6 was built on)")

#: How close an energy must sit to a training level energy to count as represented.
#: 0.02 eV is ~160 cm^-1 — loose enough to absorb the digit difference between Jofre's
#: printed 4-decimal energies and a VALD extraction, tight enough that a genuinely
#: different level fails. Widen this and the LEVEL axis stops discriminating.
LEVEL_TOL_EV = 0.02

#: Edge tolerance on the min-max boxes below, in eV / dex.
#:
#: The training energies come from Jofre et al.'s printed tables at four decimals, while
#: a measured line's energies come from a VALD extraction at whatever precision VALD
#: carries. The two disagree in the fourth decimal, so the line that DEFINES a box edge
#: can land a few times 1e-5 outside its own box and be refused — which happened, on
#: Fe I 4787.83, the very first line of the training set.
#:
#: 1e-3 is digit noise, not slack, and it cannot rescue the case this module exists for:
#: the near-IR band misses the transition-energy box by 0.028 eV at its closest, 28x this
#: tolerance. Raising it far enough to admit the IR would take a 30-fold increase, which
#: would be visible as exactly what it was.
BOX_TOL = 1e-3

#: The treatment label this engine emits. NOT a variant of ENGINE-A: same EW route,
#: but a 3D-NLTE per-line correction from a different source, so RYA-712 makes it its
#: own product with its own value/sigma/n.
TREATMENT = "ENGINE-A-3DNLTE"


# ── the recovered training domain ─────────────────────────────────────────────

@dataclass(frozen=True)
class NetworkDomain:
    """What one of the three MLPs actually saw during training."""
    network: str
    n_lines: int
    elo: tuple[float, float]
    eup: tuple[float, float]
    loggf: tuple[float, float]
    delta_E: tuple[float, float]
    lambda_air: tuple[float, float]
    level_energies: np.ndarray = field(repr=False)


_domain_cache: dict = {}


def load_training(path: Path | None = None) -> pd.DataFrame:
    p = Path(path) if path is not None else TRAINING_CSV
    if not p.exists():
        raise FileNotFoundError(
            f"Amarsi 2022 training line list not found at {p}. It is not shipped with "
            f"vendor/1L-3NErrors — regenerate it with "
            f"`python3 scripts/rya817_recover_amarsi_training_set.py`. Without it there "
            f"is no line-parameter domain check and this engine must not run (RYA-817).")
    return pd.read_csv(p)


def network_for(ion: str, elo_eV: float) -> str:
    """Which of the three MLPs handles this line. Mirrors nlte_corrections._compute_aberr."""
    if str(ion).strip() in ("II", "2", "Fe2"):
        return "fe2"
    return "lt02" if float(elo_eV) < 2.0 else "gt02"


def domains(path: Path | None = None) -> dict[str, NetworkDomain]:
    key = str(path or TRAINING_CSV)
    if key in _domain_cache:
        return _domain_cache[key]
    df = load_training(path)
    out = {}
    for net, sub in df.groupby("network"):
        species = sub["species"].iloc[0]
        same_species = df[df["species"] == species]
        levels = np.unique(np.concatenate(
            [same_species["elo_eV"].values, same_species["eup_eV"].values]))
        out[str(net)] = NetworkDomain(
            network=str(net),
            n_lines=int(len(sub)),
            elo=(float(sub.elo_eV.min()), float(sub.elo_eV.max())),
            eup=(float(sub.eup_eV.min()), float(sub.eup_eV.max())),
            loggf=(float(sub.loggf.min()), float(sub.loggf.max())),
            delta_E=(float(sub.delta_E_eV.min()), float(sub.delta_E_eV.max())),
            lambda_air=(float(sub.wavelength_air_A.min()),
                        float(sub.wavelength_air_A.max())),
            level_energies=levels,
        )
    _domain_cache[key] = out
    return out


# ── the per-line verdict ──────────────────────────────────────────────────────

@dataclass
class DomainVerdict:
    network: str
    in_domain: bool
    feature_ok: bool
    delta_E_ok: bool
    level_ok: bool
    stellar_ok: bool
    delta_E_eV: float
    reasons: list[str]

    @property
    def reason(self) -> str:
        return "; ".join(self.reasons)


def classify_line(ion: str, elo_eV: float, eup_eV: float, loggf: float, *,
                  teff: float, logg: float, vmic: float, afe: float,
                  path: Path | None = None) -> DomainVerdict:
    """Is this line inside the Amarsi 2022 MLP's training domain? Names every miss.

    `afe` is the value placed on the grid's A(Fe;3N) axis. Pass the CONVERGED 3D-NLTE
    abundance, not the 1D-LTE starting guess — see the module docstring.
    """
    net = network_for(ion, elo_eV)
    d = domains(path)[net]
    reasons: list[str] = []

    def _box(name, val, lo, hi, unit):
        lo, hi = lo - BOX_TOL, hi + BOX_TOL
        if val < lo or val > hi:
            reasons.append(f"{name} {val:.4f} {unit} outside training "
                           f"[{lo:.4f}, {hi:.4f}] ({net}, n={d.n_lines})")
            return False
        return True

    feature_ok = all([
        _box("Elo", float(elo_eV), *d.elo, "eV"),
        _box("Eup", float(eup_eV), *d.eup, "eV"),
        _box("loggf", float(loggf), *d.loggf, "dex"),
    ])

    dE = float(eup_eV) - float(elo_eV)
    delta_E_ok = _box("transition energy Eup-Elo", dE, *d.delta_E, "eV")
    if not delta_E_ok:
        reasons.append(
            f"  -> that is a WAVELENGTH statement: the network saw "
            f"{d.lambda_air[0]:.1f}-{d.lambda_air[1]:.1f} A only, and Eup-Elo is "
            f"hc/lambda_vac by construction")

    lvl = d.level_energies
    lo_ok = bool(np.min(np.abs(lvl - float(elo_eV))) <= LEVEL_TOL_EV)
    up_ok = bool(np.min(np.abs(lvl - float(eup_eV))) <= LEVEL_TOL_EV)
    level_ok = lo_ok and up_ok
    if not lo_ok:
        reasons.append(f"lower level {elo_eV:.4f} eV is not within {LEVEL_TOL_EV} eV of "
                       f"any level in the training set")
    if not up_ok:
        reasons.append(f"upper level {eup_eV:.4f} eV is not within {LEVEL_TOL_EV} eV of "
                       f"any level in the training set")

    stellar_ok = _in_grid(float(teff), float(logg), float(afe), float(vmic))
    if not stellar_ok:
        reasons.append(
            f"stellar parameters (Teff={teff:.0f}, logg={logg:.2f}, vmic={vmic:.2f}, "
            f"A(Fe;3N)={afe:.3f}) outside the published grid box "
            f"Teff {_GRID['teff'][0]:.0f}-{_GRID['teff'][1]:.0f}, "
            f"logg {_GRID['logg'][0]:.1f}-{_GRID['logg'][1]:.1f}, "
            f"vmic {_GRID['vmic'][0]:.0f}-{_GRID['vmic'][1]:.0f}, "
            f"A(Fe) {_GRID['afe'][0]:.1f}-{_GRID['afe'][1]:.1f}")

    return DomainVerdict(
        network=net,
        in_domain=feature_ok and delta_E_ok and level_ok and stellar_ok,
        feature_ok=feature_ok, delta_E_ok=delta_E_ok, level_ok=level_ok,
        stellar_ok=stellar_ok, delta_E_eV=dE, reasons=reasons)


#: How the grid's fourth axis, A(Fe;3N), is filled.
#:   'star' — the STAR's converged 3D non-LTE iron abundance, one value for every line.
#:            This is what the vendor README describes ("adopt an initial guess for the
#:            3D non-LTE abundance and iterate") and what the axis means: it is the
#:            model's iron content, which sets how saturated every line is.
#:   'line' — each line starts from its OWN 1D LTE abundance (this repo's RYA-207
#:            per-line leg). A line whose gf is wrong by 0.2 dex then places the STAR
#:            0.2 dex off on a stellar axis, and on the solar VIS band that pushes 79 of
#:            153 lines past the grid's 7.5 ceiling — a selection, not a measurement.
#: 'star' is the product mode. 'line' is kept and reported as a diagnostic because it is
#: what the archived leg does today and the difference between them is worth seeing.
AFE_AXIS_MODES = ("star", "line")


def aberr_for_line(ion: str, elo_eV: float, eup_eV: float, loggf: float,
                   a_1dlte: float, *, teff: float, logg: float, vmic: float,
                   afe3n_axis: float | None = None,
                   allow_out_of_domain: bool = False,
                   path: Path | None = None) -> tuple[float, DomainVerdict]:
    """A(3D-NLTE) - A(1D-LTE) for one line, or NaN with the reason it was refused.

    `afe3n_axis` pins the grid's A(Fe;3N) axis ('star' mode). Leave it None to keep the
    archived per-line behaviour ('line' mode).

    `allow_out_of_domain=True` exists ONLY so the out-of-domain lines can be quantified
    in a diagnostic table (how far would the network have gone?). Nothing that value
    touches may enter a product: the run script keeps it in a separate column and the
    aggregate is built from in-domain lines only.

    ORDER MATTERS in 'line' mode. The A(Fe;3N) axis test needs the converged abundance,
    and the converged abundance needs a network call — so the correction is computed
    FIRST and the verdict is formed against the result. Every axis that does NOT depend
    on the abundance (Elo/Eup/log gf, transition energy, level representation,
    Teff/logg/vmic) is decided independently of it, so a line is never admitted on a
    number the domain check would have refused.
    """
    axis_probe = float(afe3n_axis) if afe3n_axis is not None \
        else float(np.clip(a_1dlte, *_GRID['afe']))
    probe = classify_line(ion, elo_eV, eup_eV, loggf, teff=teff, logg=logg, vmic=vmic,
                          afe=axis_probe, path=path)
    if not probe.stellar_ok:
        return float("nan"), probe
    line_axes_ok = probe.feature_ok and probe.delta_E_ok and probe.level_ok
    if not line_axes_ok and not allow_out_of_domain:
        return float("nan"), probe

    ab = _apply_aberr_to_line(str(ion), float(elo_eV), float(eup_eV), float(loggf),
                              float(a_1dlte), float(teff), float(logg), float(vmic),
                              afe3n_axis=afe3n_axis)
    if afe3n_axis is not None:
        v = probe                      # the axis was pinned; the probe IS the verdict
    else:
        afe_converged = float(a_1dlte) + float(ab) if np.isfinite(ab) else float(a_1dlte)
        v = classify_line(ion, elo_eV, eup_eV, loggf, teff=teff, logg=logg, vmic=vmic,
                          afe=afe_converged, path=path)
    if not np.isfinite(ab):
        v.reasons.append("the vendored network returned no value for this line")
        v.in_domain = False
        return float("nan"), v
    if not v.in_domain and not allow_out_of_domain:
        return float("nan"), v
    return float(ab), v


def converge_star_abundance(lines: pd.DataFrame, *, teff: float, logg: float,
                            vmic: float, start: float | None = None,
                            max_iter: int = 25, tol: float = 1e-4,
                            path: Path | None = None) -> tuple[float, int, bool]:
    """Iterate the STAR's A(Fe;3N) to self-consistency ('star' axis mode).

    `lines` needs columns ion, elo_eV, eup_eV, loggf, a_1dlte and must already be
    restricted to the lines that will enter the product (in-domain, in-aggregate);
    otherwise the axis is being set by lines the product does not use.

    Returns (A(Fe;3N), iterations, converged). The median is the estimator, matching
    `band_products.build_product` — using a mean here and a median there would make the
    axis inconsistent with the value it is supposed to describe.
    """
    if lines.empty:
        return float("nan"), 0, False
    afe = float(start) if start is not None else float(np.median(lines["a_1dlte"]))
    for it in range(1, max_iter + 1):
        corrected = []
        for _, r in lines.iterrows():
            ab, _ = aberr_for_line(str(r["ion"]), float(r["elo_eV"]), float(r["eup_eV"]),
                                   float(r["loggf"]), float(r["a_1dlte"]),
                                   teff=teff, logg=logg, vmic=vmic,
                                   afe3n_axis=float(np.clip(afe, *_GRID['afe'])),
                                   allow_out_of_domain=True, path=path)
            if np.isfinite(ab):
                corrected.append(float(r["a_1dlte"]) + ab)
        if not corrected:
            return float("nan"), it, False
        new = float(np.median(corrected))
        if abs(new - afe) < tol:
            return new, it, True
        afe = new
    return afe, max_iter, False
