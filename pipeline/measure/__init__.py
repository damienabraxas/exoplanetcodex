"""Measurement handlers — one per METHOD, selected by band — RYA-713.

Ryan, 2026-08-09: *"If I was architecting the code, I would have a handler class for each
case, UV, VIS, and IR. Different tools for different work. Am I wrong?"*

Right instinct, one refinement: the handler is keyed on **method**, and the **band selects
the method**. Keying on band directly would give near-UV and NIR each their own copy of
synthesis, and VIS and red-optical each their own copy of profile fitting — four handlers,
two pairs of near-duplicates, which is precisely the drift the Ba→Al copy caused.

    band            method        telluric   continuum
    near-UV         synthesis        —       pseudo-continuum only
    VIS             profile-fit      —       true continuum
    red-optical     profile-fit    required  true continuum
    NIR             synthesis      required  post-correction only

Two handlers cover four bands, and a fifth band routes to an existing one. Everything that
differs BY BAND rather than by method — the telluric requirement, the continuum treatment,
the systematic floor — already lives in `pipeline.band_policy` and is handed to the handler
as parameters.

THE CONTROL REQUIREMENT IS PER-HANDLER
--------------------------------------
Every handler must reproduce the known optical answer before it may run anywhere else.
This does NOT transfer between handlers: the profile fitter passing at −0.013 dex says
nothing about synthesis, which fails in entirely different ways (line-list completeness,
blend modelling, broadening, pseudo-continuum placement).

The optical is the only band where any method can be falsified — A(Fe) is known there and
nowhere else. So `control_status` is a property of the HANDLER, checked before a frontier
run, and a handler that has not passed its control cannot silently produce a number that
nobody can check.
"""
from pipeline.measure.base import (  # noqa: F401
    MeasurementHandler, HandlerNotControlled, ControlResult, resolve_handler, HANDLERS)
from pipeline.measure.profile_fit import ProfileFitHandler  # noqa: F401
from pipeline.measure.synthesis import SynthesisHandler  # noqa: F401
