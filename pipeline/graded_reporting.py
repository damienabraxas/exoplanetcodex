"""
Graded vs ungraded reporting — RYA-850
======================================
Two products per cell, and the GRADED one is what gets reported.

A product is GRADED when every line in it carries a PRIMARY LABORATORY oscillator
strength (RYA-824 in VIS/IR, RYA-836 in the near-UV). It is UNGRADED when the pool
includes Kurucz semi-empirical gf, which is most lines. The distinction is already in
`error_budget.gf_term`:

    graded    0.041 dex   NIST grade B = 10% on the transition probability
    ungraded  0.170 dex   Kurucz semi-empirical, 0.1-0.3 dex (RYA-161)

and it was never wired: `derive_band_products.py` passes `gf_graded=False` unconditionally,
so every cell but one has been reported on the 0.17 term regardless of what its lines
actually are.

WHY THE GRADED VALUE IS THE PRIMARY ONE, AND WHY THAT IS NOT IN TENSION WITH RYA-842
RYA-836 found the near-UV lab-gf pool scatters WORSE than the Kurucz pool (0.652 against
0.413), and RYA-842 ratified that line selection — not gf — drives value and spread. Both
still hold. What the lab gf buys is the SYSTEMATIC, and the reported uncertainty is
    total = sqrt(stat^2 + sys^2)
where the worse scatter enters only through stat = scatter/sqrt(N), which averages down,
while the gf term enters sys, which does not. So a pool can scatter more and still report a
smaller total. That is an arithmetic consequence, not a re-litigation: it is checked per
cell in `pair_products`, and a graded cell that does NOT beat its ungraded twin on total is
flagged rather than silently promoted.

⚠️ UNGRADED IS KEPT, NOT DISCARDED. An ungraded line is usually good data with a
worse-known oscillator strength. The ungraded product stays, labelled, with its wider bar —
it is the BROAD number over every measurable line, and dropping it would destroy the
comparison the two exist for (RYA-712).

NO NEW STATISTICS. stat, sys, total, and nothing else. This deliberately does not
decompose further or invent a weighting; it matches standard solar-abundance practice
(value +/- total, with stat and sys stated).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

#: A treatment is graded when its pool is restricted to primary-lab-gf lines. The suffix is
#: the marker because RYA-712 keys a product on what produced it, and the lab-gf pool is a
#: different pool through the same engine — `1D-LTE` and `1D-LTE-LABGF` are the pair.
GRADED_SUFFIX = "-LABGF"


def is_graded(treatment: str, *, gf: str | None = None) -> bool:
    """Is this product's pool restricted to primary-laboratory gf?

    🔴 RYA-906 — ASK THE AXIS FIRST. `gf == "lab"` is a FIELD test; the `-LABGF` suffix is
    a string test, and a string test against a label whose variant set can grow is the
    RYA-869 bug class exactly. A caller that has the product's `gf` axis passes it and the
    suffix is never consulted.

    The suffix remains the fallback for the 2153 committed rows that predate the axis
    columns — the legacy label is permanent and never rewritten (Ryan: dual-label
    forever), so this must keep reading it. Fallback, not equal partner: when both are
    present the axis wins, because the axis is what the emitter measured and the suffix is
    what someone typed.
    """
    if gf is not None and str(gf).strip():
        return str(gf).strip().lower() == "lab"
    return str(treatment).endswith(GRADED_SUFFIX)


def base_treatment(treatment: str) -> str:
    """The engine a treatment belongs to, with the pool marker removed.

    `1D-LTE-LABGF` and `1D-LTE` share a base, which is what makes them a comparable pair;
    ENGINE-A and ENGINE-B do not pair with anything today because no lab-gf pool has been
    measured through them.

    ⚠️ THE STRIP IS KEYED ON THE SUFFIX ITSELF, NOT ON `is_graded`. They were the same
    test while `1D-LTE-LABGF` was the only lab-gf label; RYA-1106 corrected
    `ENGINE-A-3DNLTE` to gf="lab" and it carries no suffix, so a caller passing the axis
    would have had six characters chopped off a name that never ended in `-LABGF`.
    Removing a suffix is a question about the STRING and is now asked of the string.
    """
    t = str(treatment)
    return t[: -len(GRADED_SUFFIX)] if t.endswith(GRADED_SUFFIX) else t


def total_sigma(stat_dex: float, syst_dex: float) -> float:
    """sqrt(stat^2 + sys^2). The one combination this module performs."""
    if not np.isfinite(stat_dex) or not np.isfinite(syst_dex):
        return float("nan")
    return float(np.hypot(float(stat_dex), float(syst_dex)))


def stat_from_scatter(scatter_dex: float, n_lines: int) -> float:
    """scatter / sqrt(N) — the standard error of the mean, which DOES average down.

    This is the only place the line-to-line scatter enters the reported bar, and it is why
    a pool with a worse spread can still report a smaller total once its systematic drops.
    """
    if n_lines is None or n_lines < 1 or not np.isfinite(scatter_dex):
        return float("nan")
    return float(scatter_dex) / float(np.sqrt(n_lines))


@dataclass
class Pairing:
    """One (element, ion, band, engine) cell, with whichever of its two pools exist."""
    element: str
    ion: str
    band: str
    base: str
    graded: dict | None
    ungraded: dict | None

    @property
    def primary(self) -> dict | None:
        """The graded product when there is one; otherwise the ungraded, clearly labelled.

        A cell with no lab-gf pool is NOT dropped and NOT silently reported as graded — it
        falls back and says so, because most cells have no lab-gf pool measured yet.
        """
        return self.graded if self.graded is not None else self.ungraded

    @property
    def primary_is_graded(self) -> bool:
        return self.graded is not None

    @property
    def graded_beats_ungraded(self) -> bool | None:
        """Does promoting the graded pool actually tighten the reported bar?

        Returns None when the comparison cannot be made (one pool missing). A False here is
        the interesting case: it would mean the graded pool's worse scatter outweighs its
        better gf term, and promoting it would WIDEN the published uncertainty — which the
        ticket's premise assumes does not happen, so it is measured rather than assumed.
        """
        if self.graded is None or self.ungraded is None:
            return None
        g, u = self.graded.get("total_dex"), self.ungraded.get("total_dex")
        if not (np.isfinite(g or np.nan) and np.isfinite(u or np.nan)):
            return None
        return bool(g < u)


def annotate(df: pd.DataFrame) -> pd.DataFrame:
    """Add graded/base/stat/total columns to a product table.

    `stat_dex` and `syst_dex` are taken as given — they are what the producing route
    published — and only the total is derived here, so this cannot silently disagree with
    the cells it is describing.
    """
    d = df.copy()
    d["graded"] = d.treatment.map(is_graded)
    d["base_treatment"] = d.treatment.map(base_treatment)
    d["total_dex"] = [total_sigma(s, y) for s, y in zip(d.stat_dex, d.syst_dex)]
    return d


def pair_products(df: pd.DataFrame) -> list[Pairing]:
    """Group a product table into graded/ungraded pairs per (element, ion, band, engine)."""
    d = annotate(df)
    out: list[Pairing] = []
    keys = ["element", "ion", "band", "base_treatment"]
    for (el, ion, band, base), grp in d.groupby(keys, dropna=False):
        g = grp[grp.graded]
        u = grp[~grp.graded]
        out.append(Pairing(
            element=str(el), ion=str(ion), band=str(band), base=str(base),
            graded=(g.iloc[0].to_dict() if len(g) else None),
            ungraded=(u.iloc[0].to_dict() if len(u) else None),
        ))
    return sorted(out, key=lambda p: (p.element, p.ion, p.band, p.base))


def element_table(df: pd.DataFrame) -> pd.DataFrame:
    """The main Solar-page view: one row per (element, ion, band, engine), graded first.

    ⚠️ This does NOT average across bands or engines. RYA-712 forbids a combined product
    and `pipeline.band_products` deliberately has no `combine()`; the per-band values are
    separate measurements whose spread is the result, not noise to be collapsed. So the
    "headline" is a SELECTED cell, never a computed mean — and choosing which cell is a
    reporting decision, not arithmetic this module should make on its own.
    """
    rows = []
    for p in pair_products(df):
        prim = p.primary
        if prim is None:
            continue
        rows.append({
            "element": p.element, "ion": p.ion, "band": p.band, "engine": p.base,
            "pool": "GRADED (primary lab gf)" if p.primary_is_graded
                    else "UNGRADED (incl. Kurucz)",
            "primary_is_graded": p.primary_is_graded,
            "A": prim.get("A"),
            "n_lines": prim.get("n_lines"),
            "stat_dex": prim.get("stat_dex"),
            "syst_dex": prim.get("syst_dex"),
            "total_dex": prim.get("total_dex"),
            "graded_beats_ungraded": p.graded_beats_ungraded,
            "ungraded_A": (p.ungraded or {}).get("A"),
            "ungraded_total_dex": (p.ungraded or {}).get("total_dex"),
        })
    return pd.DataFrame(rows)


def format_value(a: float, total: float, *, places: int = 3) -> str:
    """`7.577 +/- 0.137` — the +/- is always visible, per the ticket."""
    if a is None or not np.isfinite(a):
        return "n/a"
    if total is None or not np.isfinite(total):
        return f"{a:.{places}f} +/- n/a"
    return f"{a:.{places}f} +/- {total:.{places}f}"
