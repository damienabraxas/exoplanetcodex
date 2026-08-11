"""Handler contract and band→handler routing — RYA-713."""
from __future__ import annotations

import abc
import dataclasses
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
    """What a handler proved in the optical, and therefore what it may claim elsewhere.

    THIS IS METHOD VALIDATION, NOT A SCIENCE CONSTRAINT (RYA-713).

    Ryan, 2026-08-09: *"It should not be forced to require the same exact value as VIS. We
    want to audit what is known in IR... We want to compare to VIS, sure, and hopefully
    within error bars it fits, but also check against what is found from other sources in
    the IR. And we document regardless."*

    The comparison here is **our EW against a banked EW, on the SAME LINES, in the SAME
    band**. It asks one question: does this harness measure what a validated harness
    measured? That is a property of the tool.

    It is emphatically NOT "does the IR abundance equal the optical abundance". Requiring
    that would be circular — it would tune the frontier to the control and erase exactly
    the physics we are trying to observe. Lines in different bands form at different
    depths, carry different NLTE departures, and sit on different gf scales, so a real
    cross-band difference is a RESULT.

    Cross-band abundance comparison lives in `cross_band_comparison()`, which reports and
    never gates."""
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


# ── Cross-band and external comparison: REPORTED, NEVER GATED ────────────────

@dataclass(frozen=True)
class BandComparison:
    """One band's abundance against a comparand. Carries no pass/fail, by construction.

    A cross-band difference is a scientific finding whatever its sign or size. Encoding a
    verdict here would turn an observation into a compliance check, and the first thing it
    would suppress is a real 3D, NLTE or gf-scale difference between formation regimes.
    """
    element: str
    ion: str
    band: str
    value: float
    sigma: float
    n_lines: int
    comparand_label: str        # "VIS (this work)" | a citation | an instrument
    comparand_value: float
    comparand_sigma: float
    comparand_kind: str         # "internal-cross-band" | "external-literature" | "cross-instrument"

    @property
    def difference(self) -> float:
        return self.value - self.comparand_value

    @property
    def combined_sigma(self) -> float:
        return float(np.hypot(self.sigma, self.comparand_sigma))

    @property
    def n_sigma(self) -> float:
        cs = self.combined_sigma
        return float(self.difference / cs) if cs > 0 else float("nan")

    def describe(self) -> str:
        """Report. Deliberately never returns PASS or FAIL."""
        c = "consistent within the combined uncertainty" if abs(self.n_sigma) <= 2.0 \
            else "DISCREPANT beyond the combined uncertainty — a finding, not a failure"
        return (f"{self.element} {self.ion} {self.band}: {self.value:.3f} +/- {self.sigma:.3f} "
                f"(n={self.n_lines})  vs  {self.comparand_label} "
                f"{self.comparand_value:.3f} +/- {self.comparand_sigma:.3f}  "
                f"[{self.comparand_kind}]  d={self.difference:+.3f} "
                f"({self.n_sigma:+.1f} sigma) — {c}")


def cross_band_comparison(*args, **kwargs) -> BandComparison:
    """Build a comparison. Present so the distinction from ControlResult is structural.

    There is intentionally no `passed` field and no threshold argument. If a caller wants
    to act on a discrepancy that is their decision to record, not this layer's to impose.
    """
    return BandComparison(*args, **kwargs)


def assert_not_a_science_gate(obj: Any) -> None:
    """Guard: a cross-band or external comparison must never acquire a verdict field.

    The moment one does, someone will gate a frontier band on reproducing the optical
    value, and the pipeline will start hiding real astrophysics as 'failures'.
    """
    banned = ("passed", "pass", "failed", "verdict", "tolerance", "gate", "required")
    got = {f.name.lower() for f in dataclasses.fields(obj)} if dataclasses.is_dataclass(obj) else set()
    hits = {n for n in got for b in banned if b == n or n.startswith(b)}
    if hits:
        raise RuntimeError(
            f"{obj.__name__ if isinstance(obj, type) else type(obj).__name__} gained "
            f"field(s) {sorted(hits)}. A cross-band or external comparison REPORTS; it "
            f"does not adjudicate. Forcing the IR to reproduce the VIS value would erase "
            f"the formation-depth, NLTE and gf-scale differences we are measuring.")


assert_not_a_science_gate(BandComparison)
