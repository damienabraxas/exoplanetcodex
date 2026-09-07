"""A fit that returned an impossible abundance did not measure anything — RYA-1191.

Ryan, 2026-09-07, on finding A = 4.53 in a live Fe pool: *"Why are we keeping an
abundance line that is 4,53? That seems like a bad line."*

🔴 WHAT WAS ACTUALLY WRONG. `9437.793` returns **A = 4.539 on solar_iag and 10.988 on
solar_kpno_molecfit_corrected** — the same line, the same gf, 6.4 dex apart — with
`red_chi2` of 152 and 246 and `excluded_reason` blank on both. The optimiser railed. It is
not a measurement of anything, and nothing stopped it entering the aggregate: RYA-992
built the synthesis goodness-of-fit cut and it was never wired in (`SYNTH_CONSTRAINT` is
still `None`).

⚠️ THIS IS NOT AN OUTLIER CUT, AND THE DIFFERENCE IS THE WHOLE JUSTIFICATION. Dropping
points because they are far from the others is the RYA-981 error — a quota dressed as a
cut — and RYA-515 says so directly for Fe II: the scatter there is a saturation regime and
"dropping outliers cannot fix it". This bound is on the *physical validity of a fit's
OUTPUT*: A(Fe) = 4.5 is a factor of a thousand below solar iron, A = 11.0 a factor of
thirty above. No line list, atmosphere or telluric residual produces those; a
non-convergent fit does.

TWO GUARDS WERE TRIED FIRST AND BOTH FAILED ON THE DATA, WHICH IS WHY THIS ONE IS ON A:
  * `red_chi2 > 10` flags **1262 of 1366** in-aggregate lines. The chi2 is not calibrated
    — the HARPS S1D `ERR` column is entirely NaN, so molecfit and the fitter run on a
    DEFAULT error and the absolute chi2 means nothing.
  * "no fit-constraint diagnostics" flags 334 lines of which only 3 are bad. Missing
    `sigma_A`/`frac_rise_weaker`/`ew_mA` is NORMAL on the synthesis route, not a fault.

THE WIDTH IS CHOSEN SO IT CANNOT SHAPE A RESULT. At +/-1.5 dex the bound admits everything
from A = 5.96 to A = 8.96 — a factor of a thousand in iron abundance, against a
line-to-line scatter of ~0.2 dex. It removes **15 of 1366** in-aggregate graded lines
(1.1%), and the product medians move by at most 0.021 dex while the reported sigma_stat
collapses (KP NIR ENGINE-A 0.502 -> 0.059). ⚠️ It also REPRODUCES a shipped number that
the current code does not: solar_iag NIR ENGINE-A returns to A = 7.599 / sigma 0.072 / n=6,
exactly the committed product.
"""
from __future__ import annotations

#: Solar A(Fe), the reference the bound is centred on (Asplund-scale gold, RYA-1109).
#: Centred on a FIXED published value, never on the band's own median — a bound that moves
#: with the data it is filtering can be dragged by the very fits it exists to reject.
SOLAR_A_FE = 7.46

#: Half-width, dex. A factor of ~1000 in abundance either way.
VALIDITY_HALF_WIDTH_DEX = 1.5


def fit_is_physical(abundance: float | None, element: str = "Fe") -> bool:
    """False when a fit returned an abundance no stellar atmosphere can produce."""
    if abundance is None:
        return True                       # absence is a different problem (RYA-833)
    try:
        a = float(abundance)
    except (TypeError, ValueError):
        return True
    if a != a:                            # NaN
        return True
    if element != "Fe":
        return True                       # the reference below is Fe's; refuse to guess
    return abs(a - SOLAR_A_FE) <= VALIDITY_HALF_WIDTH_DEX


def rejection_reason(abundance: float, element: str = "Fe") -> str:
    """The `excluded_reason` for a rejected fit, saying what it is and what it is not."""
    return (f"FIT-NOT-PHYSICAL: A({element}) = {float(abundance):.3f} is "
            f"{abs(float(abundance) - SOLAR_A_FE):.2f} dex from solar "
            f"{SOLAR_A_FE}, outside the +/-{VALIDITY_HALF_WIDTH_DEX} dex validity bound. "
            f"This is a NON-CONVERGENT FIT, not an outlier: the bound spans a factor of "
            f"~1000 in abundance and cannot shape a result (RYA-1191).")
