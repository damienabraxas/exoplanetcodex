"""The synthesis window per band — ONE home. RYA-967.

`SYNTH_BANDS` was declared inside `scripts/derive_band_products.py` and
`scripts/rya855_rung_audit.py` imported it FROM THAT SCRIPT. A config constant living in
an executable driver is the second-home defect this project keeps paying for
(RYA-350/353/954): the data cannot be read without importing a driver, and the driver's
import chain pulls in an atlas loader. It lives here now; the values are in
`config/synth_bands.yaml` beside `stars.yaml`, and this module is the accessor.

WHAT DID NOT CHANGE. The three pre-existing entries are carried through byte-for-byte.
Recomputing them from the invariant recovered in the YAML would move published products
(the near-UV 7.487) for a rounding difference, and RYA-832 forbids that without a stated
cause. `tests/test_synth_bands_config_rya967.py` pins them against the values that were
in the script, so the lift cannot have changed a number.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

_ROOT = Path(__file__).resolve().parent.parent
_YAML = Path(__file__).resolve().parent / "synth_bands.yaml"

#: Sentinel for "the iSpec-vendored GES v6 transitions list". It is NOT a repo-relative
#: path — it lives under `ISPEC_DIR`, which is a machine truth resolved at import by
#: `config.constants`. Writing it as a literal here would both break on every other
#: machine and trip `scripts/audit_path_literals.py` (RYA-800/810).
ISPEC_GES_V6 = "ispec_ges_v6"


@dataclass(frozen=True)
class SynthBand:
    """Everything the synthesis route needs that is a property of the BAND — RYA-837.

    The route itself is band-agnostic: `select_lines`, `fit_one` and `build_solar_context`
    take the window as an argument. Hard-coding one band's constants into the route would
    have duplicated 100 lines to change four numbers — the RYA-701 failure mode.
    """
    lo_A: float
    hi_A: float
    #: The wavelength the half-width was SET at — named in each band's original docstring
    #: and not always the band centre. See the YAML header.
    anchor_A: float
    #: The list SPEC, not the resolved path — see `linelist` below.
    linelist_spec: str
    half_width_A: float
    min_sep_A: float
    n_lines: int
    half_width_note: str
    build_hint: str

    @property
    def linelist(self) -> Path:
        """Resolved ON ACCESS, never at import.

        🔴 Eager resolution defeated the entire point of the lift. The VIS entry's list is
        the iSpec-vendored GES v6 file, so resolving it needs `pipeline.abundances_derive`,
        which does `import ispec` at module scope. Reading a band's HALF-WIDTH would then
        have required an engine to be installed — turning a config read into an engine
        dependency, which is worse than the second home this ticket set out to remove.
        Only the synthesis route asks for `.linelist`, and there iSpec is already loaded.
        """
        if self.linelist_spec == ISPEC_GES_V6:
            from pipeline.abundances_derive import _SYNTH_LINELIST_FILE
            return Path(_SYNTH_LINELIST_FILE)
        return _ROOT / self.linelist_spec

    @property
    def doppler_sigma_A(self) -> float:
        """The Doppler width at this band's ANCHOR — the unit the half-width is set in."""
        return self.anchor_A * STELLAR_SIGMA_KMS / 299792.458

    @property
    def half_width_in_doppler_sigma(self) -> float:
        return self.half_width_A / self.doppler_sigma_A


def _load() -> tuple[dict[str, SynthBand], float]:
    if not _YAML.exists():
        raise FileNotFoundError(
            f"synth_bands.yaml not found at {_YAML} — the single source of synthesis "
            f"window definitions (RYA-967). Refusing to fall back to a default: an "
            f"invented half-width is what `synthesis_route` already refuses to do.")
    data = yaml.safe_load(_YAML.read_text(encoding="utf-8"))
    k = float(data["invariant"]["half_width_in_doppler_sigma"])
    bands: dict[str, SynthBand] = {}
    for name, b in data["bands"].items():
        bands[name] = SynthBand(
            lo_A=float(b["lo_A"]), hi_A=float(b["hi_A"]),
            anchor_A=float(b["anchor_A"]),
            linelist_spec=str(b["linelist"]),
            half_width_A=float(b["half_width_A"]),
            min_sep_A=float(b["min_sep_A"]),
            n_lines=int(b["n_lines"]),
            half_width_note=str(b["half_width_note"]).strip(),
            build_hint=str(b["build_hint"]).strip())
    return bands, k


STELLAR_SIGMA_KMS: float = float(
    yaml.safe_load(_YAML.read_text(encoding="utf-8"))["invariant"]["stellar_sigma_kms"])

#: The invariant the three original entries turned out to share: a half-width is this many
#: Doppler sigmas, whatever the wavelength. Used to DERIVE a new band, never to restate an
#: old one.
HALF_WIDTH_IN_DOPPLER_SIGMA: float

SYNTH_BANDS, HALF_WIDTH_IN_DOPPLER_SIGMA = _load()


def derive_half_width_A(centre_A: float) -> float:
    """The half-width the invariant gives at this wavelength, in Angstrom.

    This is how the VIS entry was set, and how the next one should be. It is NOT applied
    at runtime — every band's value is written down in the YAML so a reader sees the
    number that was used, not a formula that might evaluate differently later.
    """
    return HALF_WIDTH_IN_DOPPLER_SIGMA * centre_A * STELLAR_SIGMA_KMS / 299792.458
