"""
NumPy compatibility shims — ONE authoritative copy (RYA-313).

Why this exists: `numpy.trapz` was REMOVED in NumPy 2.0 (renamed `numpy.trapezoid`).
`requirements.txt` floors numpy at >=1.26 with no ceiling, and RYA-517 declares the
reference stack as **py3.12 + numpy 2.2** — so on the project's own production runner
every `np.trapz` call site raises `AttributeError`. Twelve call sites carried it
(equivalent-width integration in lines_fit / cno_synthesis / abundances_derive /
pysme_nlte / the O I 6300 audit, plus five Sirius synthesis scripts).

This was invisible until RYA-313 stood up CI: the Mac dev env pins numpy 1.26, where
`trapz` still exists, so the local suite was green while the reference stack was broken.

`trapezoid` is a pure rename — same signature, same result — so this shim resolves the
name once at import and every caller uses it. Do NOT re-inline `np.trapz` anywhere; do
NOT pin numpy<2 to dodge this (that would contradict the ratified RYA-517 stack).
"""

from __future__ import annotations

import numpy as np

# numpy >= 2.0 exposes `trapezoid`; 1.x exposes `trapz`. Resolve once, loudly.
if hasattr(np, "trapezoid"):
    trapezoid = np.trapezoid
elif hasattr(np, "trapz"):
    trapezoid = np.trapz
else:  # pragma: no cover - no numpy release lacks both
    raise ImportError(
        "numpy exposes neither `trapezoid` (>=2.0) nor `trapz` (<2.0); "
        f"installed numpy = {np.__version__}"
    )

__all__ = ["trapezoid"]
