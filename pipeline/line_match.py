"""
pipeline/line_match.py — RYA-1033
=================================
THE answer to "which atomic-data row is this measured line?", in one place.

🔴 WHY THIS EXISTS: A ROUNDED WAVELENGTH IS NOT A LINE IDENTITY.

Several joins keyed a measured line to its atomic data through a 2-decimal-rounded air
wavelength (`round(wavelength_air_A, 2)`, sometimes stringified). That key is wrong three
different ways, and all three were MEASURED on this repo's own data, not argued:

1. IT SPLITS A MATCHED PAIR. `canonical_gf` and the measured pools store the same line to
   different precision — `sol_ew_results_v1` has Fe I 4787.49462, `canonical_gf` has
   4787.495, 0.38 mA apart. Rounding to 2 dp puts them on 4787.49 and 4787.50, so a line
   present in BOTH tables joins to NOTHING. Exactly 17 Fe I lines do this, every one of
   them within +/-1.2 mA of a canonical row. ZERO are genuinely absent.

2. IT IS NOT EVEN A FUNCTION OF THE VALUE. Half-way cases resolve differently depending on
   WHICH LIBRARY rounds them:

       round(6136.615, 2)            -> 6136.61     (Python: correctly-rounded decimal)
       np.round(6136.615, 2)         -> 6136.62     (numpy/pandas: scale-multiply-round)

   `scripts/promote_solar_ew.py` rounded with pandas and `abundances_derive` with Python.
   Two keys for one wavelength means the join result depends on import choices, so it is
   not reproducible even in principle. This is the sharper half of the defect: a tolerance
   is a declared approximation, but a key that disagrees with itself is a coin flip.

3. IT HIDES THE AMBIGUITY IT CANNOT RESOLVE. See below — wavelength alone does not
   identify an Fe line, and a rounded key silently picks one of the candidates.

WAVELENGTH ALONE IS NOT ENOUGH — EP IS PART OF THE IDENTITY.
`canonical_gf` holds 360 Fe I clusters whose members sit within 5 mA of each other, and in
ALL 360 the members disagree on gf: they are different transitions that happen to coincide.
7 of the 421 measured solar Fe I lines land on such a cluster, e.g.

    6065.48200 -> 6065.4820  EP 2.609  log gf -1.530  NIST-C+     <- the real line
                  6065.4850  EP 4.956  log gf -3.471  KURUCZ      <- 3 mA away, 1.9 dex off

Picking by wavelength alone is a 1.9 dex coin flip. RYA-780/852 found this the hard way and
`perline_product`/`gf_grades` already require EP; this module makes that the shared rule
rather than a habit two call sites happen to share.

WHY NOT JOIN ON `line_id` / `key_z`? Because the measured artifacts do not carry one. Every
EW and per-line product is keyed (element, ion, wavelength_air_A) — adding a stable id to
them is a schema migration across every committed and gitignored product, and it cannot fix
the pools that already exist. The declared-tolerance nearest match is the fix that works on
the data as it is; RYA-1033 names it as the sanctioned alternative.

UNRESOLVED IS LOUD, NEVER NaN. A measured line with no atomic-data row is a single-source
violation (RYA-833): it silently becomes NaN gf_tier / no loggf_reference and passes
straight through the graded/ungraded split. `require_resolved` raises instead.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

#: Wavelength agreement window between two tables describing the SAME line.
#:
#: DERIVED, not chosen: the worst table-to-table disagreement measured across the 17 lines
#: the rounded key split is 1.17 mA, and the largest over the whole measured Fe I pool is
#: under 2 mA. 5 mA clears that by ~4x while staying far below the ~180 mA VIS Fe I line
#: spacing. It is the same constant `anchor_pools._MATCH_TOL_A` and
#: `derive_band_products._EW_MATCH_TOL_A` already use, so "this line" means one thing.
MATCH_TOL_A = 0.005

#: Excitation-potential agreement window. `gf_grades.EP_TOL_EV` uses the same 0.02 eV.
#: The coincident-line clusters above differ by 0.1-2.3 eV, so this separates them with
#: room to spare while absorbing rounding in the EP columns themselves.
EP_TOL_EV = 0.02


class LineMatchError(RuntimeError):
    """A measured line could not be tied to exactly one atomic-data row."""


@dataclass(frozen=True)
class MatchResult:
    """Row indices into the atomic-data frame, one per requested wavelength.

    `index` is -1 wherever nothing matched. `ambiguous` lists the wavelengths that had more
    than one candidate the EP test could not separate — those are reported, never guessed.
    """
    index: np.ndarray
    distance_A: np.ndarray
    unresolved: list
    ambiguous: list

    @property
    def resolved(self) -> np.ndarray:
        return self.index >= 0

    @property
    def n_resolved(self) -> int:
        return int(self.resolved.sum())


def match(want_wl,
          src_wl,
          *,
          want_ep=None,
          src_ep=None,
          tol_A: float = MATCH_TOL_A,
          ep_tol_eV: float = EP_TOL_EV,
          require_ep: bool = False) -> MatchResult:
    """Nearest atomic-data row within `tol_A`, disambiguated by EP where EP is available.

    `src_wl` need not be sorted. When both `want_ep` and `src_ep` are given, candidates
    whose EP disagrees by more than `ep_tol_eV` are discarded BEFORE the nearest-wavelength
    choice — that ordering is the point, since the nearest row in wavelength is routinely
    the wrong transition (see the module docstring).

    Without EP, a wavelength that has more than one candidate is recorded as AMBIGUOUS
    rather than silently resolved to the closest one.

    `require_ep=True` (RYA-1037) additionally refuses to resolve a SINGLE candidate without
    EP. The default stays False so existing callers are unchanged, but the strict mode is
    what new code should ask for, because "only one candidate in the window" is not the same
    as "the right transition": RYA-853 found `crosscheck_nist()` stamping a NIST grade on
    Fe I 6065.490 from a lone in-window row whose EP was 2.608 eV against our line's 4.956.
    One candidate, wrong level, no ambiguity flag to warn anyone.
    """
    if require_ep and (want_ep is None or src_ep is None):
        raise LineMatchError(
            "match(require_ep=True) called without excitation potentials. A wavelength "
            "alone does not identify a transition; emit EP upstream (RYA-871/1036) rather "
            "than relaxing the key (RYA-1037).")
    want_wl = np.asarray(want_wl, dtype=float)
    src_wl = np.asarray(src_wl, dtype=float)
    if src_wl.size == 0:
        return MatchResult(np.full(want_wl.shape, -1, int),
                           np.full(want_wl.shape, np.inf),
                           [float(w) for w in want_wl], [])

    order = np.argsort(src_wl, kind="stable")
    s_wl = src_wl[order]
    s_ep = None
    if src_ep is not None:
        s_ep = np.asarray(src_ep, dtype=float)[order]
    w_ep = None if want_ep is None else np.asarray(want_ep, dtype=float)

    idx = np.full(want_wl.shape, -1, dtype=int)
    dist = np.full(want_wl.shape, np.inf, dtype=float)
    unresolved: list = []
    ambiguous: list = []

    lo_all = np.searchsorted(s_wl, want_wl - tol_A, side="left")
    hi_all = np.searchsorted(s_wl, want_wl + tol_A, side="right")

    for n, (w, lo, hi) in enumerate(zip(want_wl, lo_all, hi_all)):
        cand = np.arange(lo, hi)
        if cand.size == 0:
            # Report the miss distance too: "no row within tolerance" and "no row at all"
            # are different diagnoses and the caller's error message should say which.
            near = int(np.clip(np.searchsorted(s_wl, w), 0, len(s_wl) - 1))
            unresolved.append((float(w), float(abs(s_wl[near] - w))))
            continue

        if cand.size > 1 and s_ep is not None and w_ep is not None and np.isfinite(w_ep[n]):
            keep = cand[np.abs(s_ep[cand] - w_ep[n]) <= ep_tol_eV]
            if keep.size:
                cand = keep

        if cand.size > 1:
            # Same line listed twice (HFS components, duplicate ingests) is NOT ambiguity —
            # ambiguity is candidates that would give DIFFERENT answers. Collapse on the
            # wavelength itself; anything left over is a genuine fork.
            if np.ptp(s_wl[cand]) > 0:
                ambiguous.append((float(w), [float(x) for x in s_wl[cand]]))
                continue

        best = cand[int(np.argmin(np.abs(s_wl[cand] - w)))]
        idx[n] = order[best]
        dist[n] = abs(s_wl[best] - w)

    return MatchResult(idx, dist, unresolved, ambiguous)


def require_resolved(result: MatchResult,
                     *,
                     what: str,
                     species: str = "",
                     source: str = "canonical_gf") -> np.ndarray:
    """Return the row indices, or RAISE naming every line that failed (RYA-833).

    🔴 The whole point of RYA-1033. An unmatched line used to become NaN gf_tier and travel
    on as "ungraded", which is indistinguishable in the output from a real Kurucz-tier line
    — a wrong answer that looks like a normal one. Refusing is the only honest option: the
    caller either fixes the atomic data or declares the exclusion, but nothing publishes a
    provenance it does not have.
    """
    sp = f"{species} " if species else ""
    if result.ambiguous:
        lines = "\n".join(
            f"    {w:.5f} -> {len(c)} candidates within {MATCH_TOL_A} A: "
            + ", ".join(f"{x:.5f}" for x in c)
            for w, c in result.ambiguous[:12])
        raise LineMatchError(
            f"{len(result.ambiguous)} {sp}line(s) in {what} match MORE THAN ONE {source} "
            f"row and no excitation potential was available to separate them. Wavelength "
            f"alone does not identify an Fe line — the candidates below carry different EP "
            f"and different gf, so picking the nearest is a guess, not a match. Supply "
            f"`ep_eV` for these lines or adjudicate them in {source}:\n{lines}")
    if result.unresolved:
        lines = "\n".join(f"    {w:.5f}  (nearest {source} row {1000 * d:.2f} mA away)"
                          for w, d in result.unresolved[:12])
        more = (f"\n    ... and {len(result.unresolved) - 12} more"
                if len(result.unresolved) > 12 else "")
        raise LineMatchError(
            f"{len(result.unresolved)} {sp}line(s) in {what} resolve to NO {source} row "
            f"within {MATCH_TOL_A} A. A measured line with no atomic-data provenance must "
            f"not travel on as NaN gf_tier (RYA-833/RYA-1033) — it would be published as "
            f"'ungraded' and be indistinguishable from a real Kurucz-tier line. Either add "
            f"the line to {source} or exclude it explicitly:\n{lines}{more}")
    return result.index


def match_frames(measured: pd.DataFrame,
                 atomic: pd.DataFrame,
                 *,
                 wl_col: str = "wavelength_air_A",
                 ep_col: str | None = "ep_eV",
                 atomic_wl_col: str = "wavelength_air_A",
                 atomic_ep_col: str = "excitation_potential_eV",
                 tol_A: float = MATCH_TOL_A) -> MatchResult:
    """`match` for two DataFrames, using EP whenever both sides carry it."""
    w_ep = None
    if ep_col and ep_col in measured.columns:
        w_ep = measured[ep_col].to_numpy(float)
    s_ep = (atomic[atomic_ep_col].to_numpy(float)
            if atomic_ep_col in atomic.columns else None)
    return match(measured[wl_col].to_numpy(float),
                 atomic[atomic_wl_col].to_numpy(float),
                 want_ep=w_ep, src_ep=s_ep, tol_A=tol_A)
