"""RYA-1083 — the NLTE effect is a PER-LINE PAIRED DIFFERENTIAL, and nothing else is.

🔴 THE DEFECT THIS CLOSES. RYA-1040's paired-comparand design says the NLTE effect is
⟨3D⟩-NLTE minus ⟨3D⟩-LTE. The products report MEDIANS, and on the graded solar Fe I run
both members published the same one:

    synth-mean3D-NLTE-gerber-stagger : A = 7.552   n=67
    synth-mean3D-LTE-gerber-stagger  : A = 7.552   n=67

Differencing the two published values gives EXACTLY 0.000. The per-line paired differential
is +0.0320, nonzero on 66 of 67 lines.

⚠️ NOT A QUANTISER -- 62 unique abundance values, minimum gap 0.001. The medians collided
because the 34th of 67 values happened to land on the same number. A COINCIDENCE, which is
worse than a bug: it will not reproduce, so re-running cannot find it.

⚠️ AND NOT CONFINED TO THAT PAIR. Across the RYA-1051 ladder the collision fires on two of
three: the atmosphere leg reads +0.0920 paired against +0.1050 differenced, and nobody had
noticed that one at all.

    median-of-differences  !=  difference-of-medians

Those are different statistics. They coincide only when the pairing is irrelevant, and the
whole point of a paired comparand is that it is not.

🔴 WHY A FUNCTION AND NOT A NOTE. The design INVITES the wrong operation: a reader doing
exactly what the contract says -- take the two products, subtract -- gets a clean,
well-formed, plausible ZERO. Nothing raises. So the effect has to be COMPUTED and EMITTED,
not left to the reader with a warning attached.

⚠️ AND THE FIX IS NOT "USE A MEAN". A mean would stop these particular medians colliding
and would hide the class rather than close it. The pipeline must compute the quantity it
means, whatever the aggregator is.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np
import pandas as pd

#: The two routes round wavelengths differently -- 4372.9817 against 4372.9820 in the same
#: pool -- so an exact join silently drops most of the pairing. Measured: exact merge kept
#: 6 of 55 lines where a tolerance join kept 55.
MATCH_TOL_A = 0.05


class UnpairableProducts(ValueError):
    """The two products share no line, so no differential exists."""


@dataclass(frozen=True)
class PairedDifferential:
    """What a consumer should READ instead of subtracting two published numbers."""
    n_paired: int
    median: float
    mean: float
    sd: float
    min: float
    max: float
    n_nonzero: int
    #: What the naive subtraction WOULD have said. Carried so the trap is visible in the
    #: artifact rather than only in a docstring -- this is the number that read 0.000.
    difference_of_aggregates: float
    aggregator: str

    @property
    def collision(self) -> bool:
        """True when the naive path disagrees with the real one by enough to mislead."""
        return abs(self.median - self.difference_of_aggregates) > 0.005

    def as_dict(self) -> dict:
        d = asdict(self)
        d["collision"] = self.collision
        return d


def paired_differential(hi, lo, *, value_col="abundance",
                        wave_col="wavelength_air_A",
                        aggregator="median", tol_A=MATCH_TOL_A) -> PairedDifferential:
    """(hi - lo) matched line by line, THEN aggregated.

    `hi` and `lo` are per-line frames -- the NLTE member and its comparand. Both are
    filtered to `in_aggregate` when the column is present, because a line excluded from one
    product's aggregate is not part of that product.
    """
    hi, lo = _agg_rows(hi), _agg_rows(lo)
    lw = np.asarray(lo[wave_col], dtype=float)
    lv = np.asarray(lo[value_col], dtype=float)
    deltas = []
    for _, r in hi.iterrows():
        d = np.abs(lw - float(r[wave_col]))
        if len(d) and d.min() <= tol_A:
            deltas.append(float(r[value_col]) - float(lv[int(np.argmin(d))]))
    if not deltas:
        raise UnpairableProducts(
            "no line is in-aggregate in BOTH products, so no differential exists. The "
            "effect is a per-line DIFFERENCE; a line measured in one leg and excluded "
            "from the other cannot contribute to it (RYA-1083).")
    v = np.asarray(deltas, dtype=float)
    agg = {"median": np.median, "mean": np.mean}[aggregator]
    naive = float(agg(np.asarray(hi[value_col], dtype=float))
                  - agg(np.asarray(lo[value_col], dtype=float)))
    return PairedDifferential(
        n_paired=len(v), median=round(float(np.median(v)), 4),
        mean=round(float(np.mean(v)), 4),
        sd=round(float(np.std(v, ddof=1)), 4) if len(v) > 1 else float("nan"),
        min=round(float(v.min()), 4), max=round(float(v.max()), 4),
        n_nonzero=int((np.abs(v) > 1e-6).sum()),
        difference_of_aggregates=round(naive, 4), aggregator=aggregator)


def _agg_rows(d: pd.DataFrame) -> pd.DataFrame:
    if "in_aggregate" in d.columns:
        return d[d.in_aggregate.astype(str).str.lower().isin(("true", "1"))]
    return d
