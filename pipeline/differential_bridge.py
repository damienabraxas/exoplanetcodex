"""
RYA-969 — the differential bridge: the benchmark ladder as a gf-cancelling chain.

WHY IT IS LEGITIMATE. The gf term cancels because the SAME PHYSICAL LINE is measured in both
stars, so the same unknown `log gf` enters both abundances and subtracts out. That is an
algebraic cancellation of a shared unknown, not a statistical trick and not tuning
(Melendez+2009 ApJ 704 L66; Bedell+2014; Nissen 2015 reach 0.01-0.02 dex this way).

🔴 WHAT IT DOES NOT DO IS GIVE YOU A SCALE. A hop yields `[X/H]_target-ref`. Turning that into an
absolute `A(X)` needs the reference star's own abundance, which comes from RYA-968's absolute
layer or a published benchmark value -- NEVER from fitting the target. The two layers compose in
exactly one direction:

    PRECISION comes from the differential. SCALE comes from the absolute.
    THE SCALE ERROR NEVER SHRINKS, however precise the differential is.

⚠️ THE LADDER IS SHORTER THAN IT LOOKS. Measured from `stars.yaml`, our largest gaps are Procyon
at +782 K and tau Ceti at -0.49 dex, against Jofre+2015's documented failures at ~2000 K
(Sun -> M giant) and ~2.5 dex (Sun -> HD140283). Every current star clears a DIRECT hop to the
Sun, so the two-step bridge is machinery for the M rung (RYA-970), not for the stars we hold.

🔴 NEAREST NEIGHBOUR IS NOT THE BRIDGE REFERENCE. alpha Cen B and 55 Cnc A are each other's
nearest neighbours (scaled distance 1.4) and neither's nearest is alpha Cen A -- a naive
nearest-neighbour chain routes alpha Cen B through 55 Cnc for no gain. Hops accumulate in
quadrature, so the objective is the SHORTEST PATH TO THE SUN. Nearest-neighbour is only the
fallback for choosing an INTERMEDIATE rung when a direct hop fails the gate.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

SUN = "solar"

#: Parameter-space scaling. A hop's difficulty is not one number, so these normalise the three
#: axes onto comparable units: 100 K, 0.1 dex in log g, 0.1 dex in [Fe/H].
SCALE_TEFF_K, SCALE_LOGG, SCALE_FEH = 100.0, 0.1, 0.1

#: Jofre+2015's DOCUMENTED failures, kept as the outer sanity bound and nothing more.
#: ⚠️ These are the points where line sharing was observed to BREAK, not a tolerance we chose --
#: and RYA-968 §3.1 is the standing warning that a borrowed constant is not a control. The real
#: gate is `line_sharing_gate`, which measures rather than predicts.
JOFRE_FAIL_DTEFF_K = 2000.0
JOFRE_FAIL_DFEH_DEX = 2.5


class GateUnevaluable(RuntimeError):
    """The gate could not be applied at all — a missing key, not a hop that is too large."""


class HopTooLarge(RuntimeError):
    """A hop that cannot be bridged. Never silently crossed — RYA-969 requires a loud failure."""


class ThresholdNotDeclared(RuntimeError):
    """RYA-968 F3, applied here too: declared, or nothing."""


@dataclass
class BridgeThresholds:
    """No defaults, for the same reason RYA-968 has none."""
    #: minimum lines shared by both stars before a hop may be bridged
    min_shared_lines: int | None = None
    #: per-star linear-COG window; derived per star, never borrowed (RYA-968 §3.1)
    rew_min: float | None = None
    rew_max: float | None = None
    #: wavelength / EP match tolerances -- BOTH keys, never wavelength alone (RYA-780)
    match_tol_A: float = 0.02
    match_tol_eV: float = 0.05

    def require(self, name):
        v = getattr(self, name)
        if v is None:
            raise ThresholdNotDeclared(
                f"bridge threshold {name!r} is not declared. It must be fixed before any "
                f"abundance is computed (RYA-968 F3), so that the chain cannot be re-pointed "
                f"after seeing which hops it would admit.")
        return v


# ── geometry ─────────────────────────────────────────────────────────────────────────
def hop_distance(a: dict, b: dict) -> float:
    """Scaled parameter-space distance between two stars' fundamental parameters."""
    return math.sqrt(((a["teff"] - b["teff"]) / SCALE_TEFF_K) ** 2
                     + ((a["logg"] - b["logg"]) / SCALE_LOGG) ** 2
                     + ((a["feh_ref"] - b["feh_ref"]) / SCALE_FEH) ** 2)


def within_known_limits(a: dict, b: dict) -> tuple[bool, str]:
    """Is this hop inside the regime where line sharing has been OBSERVED to work?

    A coarse pre-filter only. It can rule a hop OUT; it can never rule one IN — only
    `line_sharing_gate` can do that, because only it looks at the actual lines.
    """
    dT, dF = abs(a["teff"] - b["teff"]), abs(a["feh_ref"] - b["feh_ref"])
    if dT >= JOFRE_FAIL_DTEFF_K:
        return False, (f"dTeff {dT:.0f} K reaches the regime where Jofre+2015 observed line "
                       f"sharing to fail (Sun -> M giant, ~{JOFRE_FAIL_DTEFF_K:.0f} K)")
    if dF >= JOFRE_FAIL_DFEH_DEX:
        return False, (f"d[Fe/H] {dF:.2f} dex reaches the regime where Jofre+2015 observed "
                       f"failure (Sun -> HD140283, ~{JOFRE_FAIL_DFEH_DEX:.1f} dex)")
    return True, f"dTeff {dT:.0f} K, d[Fe/H] {dF:.2f} dex — inside the observed working regime"


# ── the line-sharing gate ────────────────────────────────────────────────────────────
#: 🔴 THE SHARED-LINE MINIMUM IS DERIVED, NOT TYPED (RYA-981).
#: A hop's delta is a MEAN, and a mean is only as trustworthy as the scatter estimate behind
#: it. The relative uncertainty on a sample sd is `1/sqrt(2(n-1))`, so requiring the sd to be
#: known to within `tol` fixes n:
#:
#:     n >= 1 + 1/(2 * tol**2)
#:
#: At the 25% default that is **9 lines**. The number is computed from the criterion rather
#: than chosen, so changing it means changing a stated statistical requirement -- which is a
#: reviewable act -- instead of editing a literal.
SD_RELATIVE_TOLERANCE = 0.25


def derive_min_shared_lines(sd_relative_tolerance: float = SD_RELATIVE_TOLERANCE) -> int:
    """Minimum shared lines for a statistically sound hop mean. See SD_RELATIVE_TOLERANCE."""
    if not (0 < sd_relative_tolerance < 1):
        raise ValueError(f"sd tolerance must be in (0,1), got {sd_relative_tolerance!r}")
    return int(math.ceil(1.0 + 1.0 / (2.0 * sd_relative_tolerance ** 2)))


@dataclass
class SharedSet:
    pairs: list                      # (target_row, reference_row)
    n_shared: int
    n_target: int
    n_reference: int
    reason: str = ""
    #: 🔴 "cannot be evaluated" is a THIRD STATE, distinct from "too few shared lines".
    #: Without it a product that simply does not carry REW reports as a hop that is too large,
    #: which sends a reader looking for an intermediate benchmark that would not help. Found on
    #: RYA-981: the RYA-967 SYNTH product carries no `rew`, `ew_mA` or `observed_depth` at all.
    unevaluable: str = ""

    @property
    def ok(self) -> bool:
        return self.reason == "" and self.unevaluable == ""


def line_sharing_gate(target_lines, reference_lines, th: BridgeThresholds) -> SharedSet:
    """Lines measurable and UNSATURATED in BOTH stars.

    🔴 MATCHED ON WAVELENGTH *AND* EXCITATION POTENTIAL. RYA-780 manufactured a 2.85 dex
    discrepancy from a single 0.06 A coincidence on wavelength alone; a bridge built on
    mis-paired lines would difference two different transitions and look perfectly consistent.

    Below the declared minimum the hop FAILS LOUDLY. The remedy is an intermediate rung, which
    is exactly how the M-dwarf gap (RYA-970) will present itself — never a silent bridge.
    """
    need = th.require("min_shared_lines")
    lo, hi = th.require("rew_min"), th.require("rew_max")

    # Diagnose a MISSING saturation key before diagnosing a short shared set. A synthesis
    # product carries no REW; treating that as "too few shared lines" names the wrong remedy.
    def _has_rew(rows):
        return any(r.get("rew") is not None and np.isfinite(r.get("rew", np.nan))
                   for r in rows)
    t_all, r_all = list(target_lines), list(reference_lines)
    if not _has_rew(t_all) or not _has_rew(r_all):
        s = SharedSet(pairs=[], n_shared=0, n_target=len(t_all), n_reference=len(r_all))
        s.unevaluable = (
            "no usable REW on one or both sides, so the saturation gate cannot be applied. "
            "This is NOT a hop that is too large — it is a product that does not carry the "
            "key. A synthesis product has no equivalent width and therefore no REW; supply "
            "an EW-route product, or declare a saturation proxy for the synthesis route.")
        return s

    def usable(r):
        rew = r.get("rew")
        return (rew is not None and np.isfinite(rew) and lo < rew < hi
                and r.get("abundance") is not None and np.isfinite(r.get("abundance")))

    tgt = [r for r in target_lines if usable(r)]
    ref = [r for r in reference_lines if usable(r)]
    rw = np.array([r["wavelength_air_A"] for r in ref], float) if ref else np.zeros(0)
    re_ = np.array([r.get("ep_eV", np.nan) for r in ref], float) if ref else np.zeros(0)

    pairs = []
    for t in tgt:
        if not len(rw):
            break
        dw = np.abs(rw - t["wavelength_air_A"])
        de = np.abs(re_ - t.get("ep_eV", np.nan))
        m = np.where((dw <= th.match_tol_A) & ((de <= th.match_tol_eV) | ~np.isfinite(de)))[0]
        if m.size == 1:
            pairs.append((t, ref[int(m[0])]))
    s = SharedSet(pairs=pairs, n_shared=len(pairs), n_target=len(tgt), n_reference=len(ref))
    if len(pairs) < need:
        s.reason = (f"only {len(pairs)} lines are measurable and unsaturated in BOTH stars "
                    f"(target {len(tgt)}, reference {len(ref)}), below the declared minimum "
                    f"{need}. This hop is TOO LARGE: insert an intermediate rung. Do not bridge "
                    f"it and do not lower the minimum to make it pass.")
    return s


# ── one hop ──────────────────────────────────────────────────────────────────────────
@dataclass
class Hop:
    target: str
    reference: str
    delta_dex: float
    sigma_dex: float
    n_lines: int
    per_line: list = field(default_factory=list)

    def __repr__(self):
        return (f"Hop({self.target} - {self.reference} = {self.delta_dex:+.4f} "
                f"+/- {self.sigma_dex:.4f}, n={self.n_lines})")


def bridge_hop(target_lines, reference_lines, th: BridgeThresholds, *,
               target: str, reference: str) -> Hop:
    """Line-by-line `[X/H]_target - [X/H]_reference`. The gf cancels HERE, per line.

    🔴 R1 — THE REFERENCE STAR'S ABUNDANCE IS NOT AN ARGUMENT TO THIS FUNCTION. It takes two
    line sets and returns a DIFFERENCE. There is no parameter through which a reference value
    could be adjusted to make a target look right, because no reference value enters at all;
    it is applied later, once, in `chain_to_sun`.
    """
    s = line_sharing_gate(target_lines, reference_lines, th)
    if s.unevaluable:
        raise GateUnevaluable(f"{target} -> {reference}: {s.unevaluable}")
    if not s.ok:
        raise HopTooLarge(f"{target} -> {reference}: {s.reason}")
    d = np.array([t["abundance"] - r["abundance"] for t, r in s.pairs], float)
    return Hop(target=target, reference=reference,
               delta_dex=float(np.mean(d)),
               sigma_dex=float(np.std(d, ddof=1) / math.sqrt(d.size)) if d.size > 1 else float("nan"),
               n_lines=int(d.size),
               per_line=[{"wavelength_air_A": t["wavelength_air_A"],
                          "delta_dex": float(t["abundance"] - r["abundance"])}
                         for t, r in s.pairs])


# ── the chain ────────────────────────────────────────────────────────────────────────
def plan_chain(star: str, params: dict, *, reachable=None) -> list:
    """The chain of hops from `star` down to the Sun.

    🔴 SHORTEST PATH TO THE SUN, NOT NEAREST NEIGHBOUR. Each hop adds its own uncertainty in
    quadrature, so a chain that visits a closer star on the way is WORSE than a direct hop that
    the gate already clears. The rule is therefore:

        take the DIRECT hop to the Sun whenever it is inside the observed working regime;
        only when it is not, insert the nearest intermediate that IS, and recurse.

    On our current star set this returns a single direct hop for every star -- which is the
    honest, measured answer, not a simplification (see the module docstring).
    """
    if star == SUN:
        return []
    if star not in params:
        raise KeyError(f"no fundamental parameters for {star!r}; a bridge cannot be planned "
                       f"for a star whose Teff/logg/[Fe/H] we do not hold (RYA-957)")
    ok, why = within_known_limits(params[star], params[SUN])
    if ok:
        return [(star, SUN, why)]

    # The direct hop is outside the regime: pick the nearest star that is closer to the Sun
    # than `star` is, and that `star` can itself reach. Nearest-neighbour earns its place ONLY
    # here, choosing an intermediate -- never as the primary chain rule.
    pool = [n for n in (reachable or params) if n not in (star, SUN) and n in params]
    d_sun = hop_distance(params[star], params[SUN])
    cands = [(hop_distance(params[star], params[n]), n) for n in pool
             if hop_distance(params[n], params[SUN]) < d_sun
             and within_known_limits(params[star], params[n])[0]]
    if not cands:
        raise HopTooLarge(
            f"{star} cannot reach the Sun directly ({why}) and no intermediate rung exists "
            f"among {sorted(pool)}. This is the M-dwarf-gap shape: the ladder needs a new "
            f"benchmark here (RYA-970), not a longer stride.")
    _, mid = min(cands)
    return [(star, mid, within_known_limits(params[star], params[mid])[1])] + \
        plan_chain(mid, params, reachable=reachable)


@dataclass
class BridgedAbundance:
    """A target's two products. BOTH are always reported (RYA-969 §4)."""
    star: str
    chain: list
    differential_dex: float
    differential_sigma_dex: float
    #: absolute = reference's own absolute + the chained difference
    absolute_dex: float | None
    #: 🔴 from RYA-968's laboratory anchor. Does NOT shrink as the differential gets better.
    zero_point_dex: float | None

    def report(self) -> str:
        path = " -> ".join([self.star] + [h.reference for h in self.chain])
        a = "n/a" if self.absolute_dex is None else f"{self.absolute_dex:.4f}"
        z = "n/a" if self.zero_point_dex is None else f"{self.zero_point_dex:.4f}"
        return (f"{path}\n"
                f"  differential [X/H] = {self.differential_dex:+.4f} "
                f"+/- {self.differential_sigma_dex:.4f} dex  (gf cancelled)\n"
                f"  absolute     A(X)  = {a} dex, zero-point +/- {z} "
                f"(laboratory anchor; does not shrink)")


def chain_to_sun(hops, *, star: str, solar_absolute_dex=None,
                 zero_point_dex=None) -> BridgedAbundance:
    """Compose hops into the two products.

    🔴 R2 — THE CHAIN TRANSPORTS A DIFFERENCE. Each hop is a delta; the absolute value is
    assembled ONCE, at the end, from the chain plus an anchor supplied from outside. Nothing
    here can adjust a reference to improve a target.
    🔴 R3 — the chain is an INPUT. This function cannot re-point it, so a hop cannot be swapped
    after seeing the answer it gives.
    """
    if not hops:
        raise ValueError(f"no hops for {star!r}; a star with no chain has no differential "
                         f"product, and that is a refusal rather than a zero")
    tot = float(sum(h.delta_dex for h in hops))
    sig = float(math.sqrt(sum(
        (h.sigma_dex ** 2) for h in hops if np.isfinite(h.sigma_dex))))
    absolute = None if solar_absolute_dex is None else float(solar_absolute_dex + tot)
    return BridgedAbundance(star=star, chain=list(hops), differential_dex=tot,
                            differential_sigma_dex=sig, absolute_dex=absolute,
                            zero_point_dex=zero_point_dex)


def compare_to_benchmark(bridged: BridgedAbundance, published_dex: float,
                         source: str) -> dict:
    """Validation against a published benchmark abundance.

    🔴 R4 — THIS IS A REPORT, NEVER A TERM. It returns a finding and changes nothing. The moment
    a bridge is adjusted to match a published value, it transports an assumption instead of a
    difference, and the whole layer's justification (gf cancels) is gone.
    """
    if bridged.absolute_dex is None:
        return {"status": "no absolute product to compare", "source": source}
    d = bridged.absolute_dex - published_dex
    return {"star": bridged.star, "ours": bridged.absolute_dex, "published": published_dex,
            "difference_dex": float(d), "source": source,
            "status": "REPORTED — a finding, never fed back into the bridge (RYA-969 R4)"}
