"""Synthesis handler — near-UV and NIR — RYA-713.

NOT YET FUNCTIONAL. This file exists so the routing is honest: the near-UV and NIR bands
resolve to a handler that states precisely what it needs, rather than silently falling
through to a method those bands forbid.

WHY THE UV NEEDS THIS AND CANNOT REUSE THE PROFILE FITTER
---------------------------------------------------------
Median line separation in the near-UV is 0.146 A -- SMALLER than a strong line's own
wings. No interval contains one profile and excludes its neighbours, and there is no
isolated profile to fit. Only a method that models every contributor in the window at
once is valid there.

WHAT IT NEEDS, EXPLICITLY
-------------------------
1. A synthesis engine call. Turbospectrum (`bsyn`) is staged on Sirius and is already
   driven by `scripts/ts_gerber_gate.py`; the per-line abundance-fitting loop is the new
   part, not the engine.
2. A COMPLETE line list for the window -- every contributor, not just the target element.
   Measured completeness (RYA-713): near-UV 77%, VIS 55%, red-optical 37%, NIR 22%. The
   near-UV is the best-catalogued band we have, which is why it is the tractable one.
3. A declared pseudo-continuum treatment. Near-UV median flux runs 0.283-0.805; the true
   continuum is never observed, so the systematic must be carried, not hidden.
4. For the NIR, a telluric model -- `policy.telluric_required` is True and before
   correction the observed flux is not a stellar spectrum at all.

AND IT NEEDS ITS OWN CONTROL
----------------------------
The profile fitter passing at -0.013 dex licenses NOTHING here. Synthesis fails in
different ways: incomplete line lists, wrong blend abundances, broadening, continuum
placement. `assert_controlled()` will refuse to let this run in the near-UV until it has
reproduced the known optical answer in the VIS -- which is the only band that can
falsify it.
"""
from __future__ import annotations

from typing import Any

from pipeline.band_policy import BandPolicy
from pipeline.band_products import LineMeasurement
from pipeline.measure.base import MeasurementHandler, register


class SynthesisHandler(MeasurementHandler):
    method = "synthesis"

    def measure_line(self, wav, flux, *, element, ion, wavelength_A, instrument,
                     policy: BandPolicy, pre_normalised: bool,
                     context: dict[str, Any]) -> LineMeasurement:
        raise NotImplementedError(
            f"SynthesisHandler is not implemented, so {element} {ion} "
            f"{wavelength_A:.3f} A in the {policy.name} band cannot be measured yet.\n"
            f"  Needs: (1) a bsyn per-line abundance fit on Sirius, (2) a complete "
            f"window line list, (3) a declared pseudo-continuum treatment"
            + (", (4) a telluric model" if policy.telluric_required else "") + ".\n"
            f"  And its OWN optical control -- the profile fitter's licence does not "
            f"transfer.")


register(SynthesisHandler())
