"""Handler contract and band→handler routing — RYA-713."""
from __future__ import annotations

import abc
from dataclasses import dataclass
from typing import Any

import numpy as np

from pipeline.band_policy import BandPolicy, resolve, BandPolicyError
from pipeline.band_products import LineMeasurement


class HandlerNotControlled(RuntimeError):
    """A handler tried to run in a band that cannot falsify it, without having passed
    the optical control that can."""


@dataclass(frozen=True)
class ControlResult:
    """What a handler proved in the optical, and therefore what it may claim elsewhere."""
    handler: str
    element: str
    n_lines: int
    median_ratio: float          # measured / reference EW
    mad_ratio: float
    dex_offset: float            # log10(median_ratio) — the harness systematic
    passed: bool
    tolerance_dex: float
    evidence: str

    def summary(self) -> str:
        v = "PASSED" if self.passed else "FAILED"
        return (f"{self.handler}: control {v} on {self.element} — {self.n_lines} lines, "
                f"median ratio {self.median_ratio:.3f}, offset {self.dex_offset:+.4f} dex "
                f"(tolerance ±{self.tolerance_dex})")


# A method may claim a frontier band only if its optical residual is inside this. It is
# NOT a quality target to tune toward -- it is the point past which a harness systematic
# would dominate the science. RYA-561's ratification gate is ±0.10 dex; a harness that
# eats half of that before any physics has not earned a frontier run.
CONTROL_TOLERANCE_DEX = 0.05


class MeasurementHandler(abc.ABC):
    """One measurement METHOD. Bands select it; they do not subclass it.

    Subclasses implement `measure_line`. Everything band-dependent arrives via `policy`.
    """

    #: the `pipeline.band_policy` method name this handler implements
    method: str = ""

    def __init__(self) -> None:
        self._control: ControlResult | None = None

    # ── control ───────────────────────────────────────────────────────────────
    @property
    def control(self) -> ControlResult | None:
        return self._control

    def record_control(self, result: ControlResult) -> None:
        if result.handler != self.name:
            raise ValueError(f"control result is for {result.handler!r}, not {self.name!r}")
        self._control = result

    def assert_controlled(self, policy: BandPolicy) -> None:
        """A handler may run un-controlled ONLY in the band that can falsify it.

        Anywhere else there is no reference value, so an uncontrolled run produces a
        number nobody can check — which is the wild-west failure this exists to prevent.
        """
        from pipeline.band_policy import CONTROL_BAND
        if policy.name == CONTROL_BAND:
            return  # running the control itself, or working where truth exists
        if self._control is None:
            raise HandlerNotControlled(
                f"{self.name} has not been controlled against the known optical answer, so "
                f"it may not run in {policy.name} ({policy.lo_A:.0f}-{policy.hi_A:.0f} A) "
                f"where nothing can falsify it. Run the control in "
                f"{CONTROL_BAND} first.\n"
                f"  Control status does NOT transfer between handlers — each method fails "
                f"in its own way.")
        if not self._control.passed:
            raise HandlerNotControlled(
                f"{self.name} FAILED its optical control and may not run in {policy.name}.\n"
                f"  {self._control.summary()}\n  {self._control.evidence}")

    def systematic_dex(self) -> float:
        """The handler's own measured systematic, for the frontier error budget.

        Not assumed zero: it is the optical residual, measured. A handler with no control
        has no measured systematic and must not pretend to one.
        """
        if self._control is None:
            raise HandlerNotControlled(
                f"{self.name} has no measured systematic because it has not been "
                f"controlled. Assuming zero would understate every frontier error bar.")
        return abs(self._control.dex_offset)

    # ── measurement ───────────────────────────────────────────────────────────
    @property
    def name(self) -> str:
        return type(self).__name__

    @abc.abstractmethod
    def measure_line(self, wav: np.ndarray, flux: np.ndarray, *, element: str, ion: str,
                     wavelength_A: float, instrument: str, policy: BandPolicy,
                     pre_normalised: bool, context: dict[str, Any]) -> LineMeasurement:
        """Measure one line. Must return a LineMeasurement, never raise for a bad line —
        a line that cannot be measured is QUARANTINED with a reason (RYA-711), because a
        raised exception loses the line silently and a silent drop is the defect."""

    def prepare(self, policy: BandPolicy, context: dict[str, Any]) -> None:
        """Optional per-run setup. Default: verify the band permits this method."""
        if self.method not in policy.permitted_methods:
            raise BandPolicyError(
                f"{self.name} implements {self.method!r}, which the {policy.name} band "
                f"does not permit (permitted: {policy.permitted_methods}).\n"
                f"  why: {policy.justification}")


HANDLERS: dict[str, MeasurementHandler] = {}


def register(handler: MeasurementHandler) -> MeasurementHandler:
    HANDLERS[handler.method] = handler
    return handler


def resolve_handler(wavelength_A: float) -> MeasurementHandler:
    """Band → policy → permitted method → handler. The routing Ryan asked for, with the
    seam at the method rather than the band so two bands sharing a method share code."""
    pol = resolve(wavelength_A)
    for m in pol.permitted_methods:
        if m in HANDLERS:
            return HANDLERS[m]
    raise BandPolicyError(
        f"the {pol.name} band permits {pol.permitted_methods} but no handler is "
        f"registered for any of them. Register one rather than relaxing the policy.")
