"""A fit that returned an impossible abundance did not measure anything — RYA-1191.

Ryan, 2026-09-07, on finding A = 4.53 in a live Fe pool: *"Why are we keeping an
abundance line that is 4,53? That seems like a bad line."*

🔴 WHAT WAS ACTUALLY WRONG. `9437.793` returns **A = 4.539 on solar_iag and 10.988 on
solar_kpno_molecfit_corrected** — the same line, the same gf, 6.4 dex apart — with
`red_chi2` of 152 and 246 and `excluded_reason` blank on both. The optimiser railed. It is
not a measurement of anything, and nothing stopped it entering the aggregate.

⚠️ AND "RYA-992 SHOULD HAVE CAUGHT THIS" IS NOT QUITE RIGHT — I said it first and it is
worth correcting, because the real reason is more interesting than an omission. RYA-992
DID land: `synth_gof_cut(instrument)` and the per-arm `ARM_SCALE` registry exist and ARE
called, from `pipeline/gf_empirical.py`. What it deliberately did not do is put a
frac_rise THRESHOLD on the band-product route — RYA-847's sweep (9 cells, 581 synthesis
lines) refuted every candidate threshold, and `constraint_gate` documents
`SYNTH_CONSTRAINT` as staying "None PERMANENTLY rather than pending" for that reason. That
was a ratified decision, not a gap.

🔴 THE GAP IS THAT A CONSTRAINT GATE AND A VALIDITY BOUND ARE DIFFERENT QUESTIONS. Both
arms of `constraint_gate.verdict` ask whether the fit was CONSTRAINED; neither asks
whether the answer is POSSIBLE. And 9437.793 slips both: its `frac_rise_weaker` is NaN, so
arm 1's non-minimum check cannot fire (it is guarded by `np.isfinite`), and arm 2 returns
`GateVerdict(True)` by default because no ratified cut exists. A line carrying NO
constraint metrics at all is therefore treated as fine — which is the one case where
"was it constrained?" has no answer and "is 4.53 a possible iron abundance?" still does.

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


def solar_reference(element: str) -> float | None:
    """The published solar A(X) this bound centres on, per element — RYA-1214.

    🔴 THE GUARD WAS INERT FOR EVERY ELEMENT BUT IRON, AND A CNO POOL PAID FOR IT.
    `fit_is_physical` returned True unconditionally for `element != "Fe"` — a deliberate
    refusal to guess a reference it did not have. It has one: `SOLAR_ASPLUND2021` is the
    same fixed published table `SOLAR_A_FE = 7.46` is taken from, for 28 elements. So the
    refusal was protecting against a lookup that already existed.

    What it cost: C I 4890.653 fitted to **A = 5.443** on solar_harps_molecfit_corrected —
    3.0 dex below solar carbon, a factor of a thousand — and entered the VIS C I aggregate
    with a blank `excluded_reason`. The two Kitt Peak holdings excluded the same line as
    `edge_pinned`, so the pool that kept it is the one where no other gate fired either.
    Identical in shape to the A(Fe) = 4.53 that motivated RYA-1191, one element over.

    ⚠️ THE HALF-WIDTH IS NOT RE-DERIVED AND DOES NOT NEED TO BE. RYA-1191's justification
    for +/-1.5 dex is a statement about NON-CONVERGENT FITS, not about iron: it is a factor
    of ~1000 either way against a line-to-line scatter of ~0.2 dex, so it cannot shape a
    result. That argument transfers unchanged; the element-specific part is the CENTRE,
    and the centre comes from a published table rather than from the band's own median for
    exactly the reason the Fe constant does.

    Returns None for an element the table does not carry — which keeps the original
    refusal for that case instead of inventing a centre.
    """
    if str(element).strip() == "Fe":
        return SOLAR_A_FE          # unchanged, and not routed through the table
    from config.constants import SOLAR_ASPLUND2021
    v = SOLAR_ASPLUND2021.get(str(element).strip())
    return float(v) if v is not None else None

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
    ref = solar_reference(element)
    if ref is None:
        return True                       # no published centre for it; refuse to guess
    return abs(a - ref) <= VALIDITY_HALF_WIDTH_DEX


def rejection_reason(abundance: float, element: str = "Fe") -> str:
    """The `excluded_reason` for a rejected fit, saying what it is and what it is not."""
    ref = solar_reference(element)
    if ref is None:                       # unreachable from fit_is_physical; explicit anyway
        raise ValueError(f"no published solar reference for {element!r}; a rejection "
                         f"reason must name the centre it measured against")
    return (f"FIT-NOT-PHYSICAL: A({element}) = {float(abundance):.3f} is "
            f"{abs(float(abundance) - ref):.2f} dex from solar "
            f"{ref}, outside the +/-{VALIDITY_HALF_WIDTH_DEX} dex validity bound. "
            f"This is a NON-CONVERGENT FIT, not an outlier: the bound spans a factor of "
            f"~1000 in abundance and cannot shape a result (RYA-1191).")
