"""Per-band uncertainty budget — RYA-713.

Ryan, 2026-08-09: *"I am ok with IR having bigger Error Bars if the science backs it"* ·
*"same with UV"* · *"ratify the uncertainty if it is legit. Gotta keep track of all that.
We said diff error bars for diff bands, so all good in the hood."*

WHY THIS IS A MODULE AND NOT A NUMBER
-------------------------------------
An uncertainty that is a single float has already lost the information that matters. The
terms behave differently under more lines:

* **RANDOM** terms average down as ``σ/√N``. Line-to-line scatter is the honest one.
* **SYSTEMATIC** terms **do not average down at any N**. An oscillator-strength scale
  offset, a pseudo-continuum that is never observed, a harness residual measured against
  a known answer — collecting a thousand more lines does not shrink any of them.

Reporting only ``σ/√N`` on a frontier band would understate the truth by an order of
magnitude: with 271 lines and σ_gf ≈ 0.2 dex the random part is ≈ 0.012 dex, while the
systematic floor underneath it is several times larger. So the budget carries both, keeps
them apart, and prints them apart.

WHY THE BANDS DIFFER, CONCRETELY
--------------------------------
This is the ratified "different error bars for different bands" made arithmetic. The terms
are not chosen per band — they are *present* per band, because the physics differs:

    term                      VIS    red-optical   near-UV   NIR
    line-to-line scatter       ✓          ✓           ✓        ✓     random
    gf scale (ungraded)        ✓          ✓           ✓        ✓     systematic
    harness residual           ✓          ✓           ✓        ✓     systematic
    pseudo-continuum           –          –           ✓        –     systematic
    telluric residual          –          ✓           –        ✓     systematic

The near-UV carries a pseudo-continuum term because its true continuum is never observed
(median flux 0.283–0.805). The IR bands carry a telluric term because that residual is
epoch- and airmass-dependent and cannot be calibrated once. The VIS carries neither, which
is exactly why it is the control band.

WHAT THIS REFUSES TO DO
-----------------------
``total()`` returns a *pair*. There is deliberately no method that collapses statistical
and systematic into one number, because doing so is how a frontier measurement comes to
look like an optical one.

UNMEASURED IS NOT ZERO — RYA-907
--------------------------------
A term whose value **could not be measured from this product** carries ``dex=None``, and
``None`` never becomes ``0.0`` anywhere below. The distinction is not pedantry: two
published Fe II red-optical cells reported ``stat_dex = 0.0`` at ``n = 1``, because
``build_product`` correctly returned ``sigma=None`` (*one line has no line-to-line
scatter*) and the driver then wrote ``(product.sigma or 0.0)``. One character turned
*"we could not measure this"* into *"we measured this and it is zero"* — different
claims, and only one of them was true. Those cells then advertised the **tightest**
statistical bar in the whole matrix while being its **least** constrained measurement.

This is the same shape RYA-873 fixed for the harness residual, in the other half of the
budget. So the rule is now structural rather than remembered:

* ``Term.contribution()`` **raises** on an unmeasured term. It cannot silently
  contribute zero, because there is no code path where it returns a number.
* ``statistical()`` never returns a value below the RYA-771 quantiser floor, and
  ``stat_basis()`` says which of *measured* / *floored* / *unmeasured* produced it.
* ``assert_stat_publishable()`` refuses to let a zero reach an artifact, naming the cell.
"""
from __future__ import annotations


from decimal import Decimal, ROUND_HALF_UP as _ROUND_HALF_UP

#: Decimal places every published dex bar is rounded to.
DEX_PLACES = 4


def round_dex(x: float, places: int = DEX_PLACES) -> float:
    """Round a dex quantity for publication, deterministically. THE single rounding rule.

    🔴 WHY THIS EXISTS (RYA-1084). `stat_dex=round(stat, 4)` was written at five separate
    call sites, and one product's published bar moved 0.0217 -> 0.0218 between two runs
    that agreed on every other field. The cause is not a Python-vs-numpy disagreement —
    `round()` and `np.round()` return the SAME answer at every point here. It is that

        0.02175 as a float is 0.0217499999999999985, one ULP BELOW the decimal tie

    so `round()` gives 0.0217, while the very next representable float gives 0.0218. A
    difference of one ULP in an upstream RMS — invisible in `stat_basis`, which prints
    five decimals and shows "0.02175" for both — became a visible change in a published
    uncertainty.

    The rule pinned here is HALF-UP on the SHORTEST ROUND-TRIPPING DECIMAL (`repr`), which
    is stable across exactly that pair: both 0.0217499999999999985 and its ULP-neighbour
    repr as values that quantize to 0.0218. Measured against `round()` on 200,000 random
    values it differs on ZERO of them — it changes nothing but the tie behaviour, which is
    the only thing that was ambiguous.

    ⚠️ WHAT IT DOES NOT DO. It does not make a genuinely different number the same. Values
    that straddle a boundary by more than the tie neighbourhood (0.0217499 vs 0.0217501)
    still round apart, correctly. Upstream 1-ULP nondeterminism is a separate exposure and
    a rounding rule cannot stand in for fixing it.
    """
    return float(Decimal(repr(float(x))).quantize(Decimal(1).scaleb(-places),
                                                  rounding=_ROUND_HALF_UP))


import math
from dataclasses import dataclass, field

from pipeline.band_policy import resolve


#: The floor beneath which no published sigma is honest, whatever the line count —
#: RYA-771. iSpec writes the trial abundance with `%.2f` in both `bsyn` and `babsma`, so
#: EW(A) is a STAIRCASE with 0.01 dex treads: two abundances closer together than one
#: tread produce byte-identical synthetic spectra, and no amount of averaging can resolve
#: below the step. Declared HERE and imported, never re-typed at a call site — a number
#: written down twice drifts (RYA-845).
QUANTISER_FLOOR_DEX = 0.01


class UnmeasuredTerm(ValueError):
    """Asked for the numeric contribution of a term that was never measured.

    Deliberately an exception rather than a zero. A zero is a MEASUREMENT — it asserts
    that the quantity was determined and found negligible — and the whole RYA-907 defect
    was a `None` being spent as if it were one.
    """


class UnpublishableStat(ValueError):
    """A statistical uncertainty that must not reach an artifact."""


@dataclass(frozen=True)
class Term:
    """One contribution to the budget, and how it behaves under more lines.

    `dex is None` means NOT MEASURABLE FROM THIS PRODUCT. It is a third state, distinct
    from both a measured value and a measured zero, and `contribution()` refuses to
    collapse it into either.
    """
    name: str
    dex: float | None
    averages_down: bool     # True => random, scales as 1/sqrt(N); False => systematic floor
    source: str             # where the number comes from -- never "assumed"

    @property
    def measured(self) -> bool:
        return self.dex is not None

    def contribution(self, n_lines: int) -> float:
        if self.dex is None:
            raise UnmeasuredTerm(
                f"{self.name}: this term was NOT MEASURED for this product ({self.source}). "
                f"It has no numeric contribution. Treating it as 0.0 would publish "
                f"'measured and negligible' for a quantity nobody determined (RYA-907).")
        if not self.averages_down:
            return abs(self.dex)
        if n_lines < 1:
            raise ValueError(f"{self.name}: cannot scale a random term with n_lines={n_lines}")
        return abs(self.dex) / math.sqrt(n_lines)


@dataclass
class ErrorBudget:
    """The uncertainty on one (element × band × treatment) product."""
    element: str
    band: str
    n_lines: int
    terms: list[Term] = field(default_factory=list)

    def add(self, term: Term) -> "ErrorBudget":
        if not term.source.strip():
            raise ValueError(
                f"term {term.name!r} has no source. Every term must say where its number "
                f"comes from -- an unsourced uncertainty is indistinguishable from a guess, "
                f"and the whole point of this budget is that a large bar is defensible.")
        self.terms.append(term)
        return self

    # ── the two halves, kept apart on purpose ────────────────────────────────
    def unmeasured_terms(self) -> list[Term]:
        """Terms this product could not determine. Never silently worth zero."""
        return [t for t in self.terms if not t.measured]

    def measured_statistical(self) -> float | None:
        """The RMS of the random terms that WERE measured — None if none were.

        Separate from `statistical()` so the floor is applied in exactly one place and
        the raw arithmetic stays inspectable underneath it.
        """
        rs = [t.contribution(self.n_lines) for t in self.terms
              if t.averages_down and t.measured]
        return math.sqrt(sum(r * r for r in rs)) if rs else None

    def statistical(self) -> float:
        """The PUBLISHED statistical bar. Never below the RYA-771 quantiser floor.

        Two ways to land on the floor, and `stat_basis()` distinguishes them:

        * every random term was measured, and their RMS is genuinely below one tread of
          the 0.01 dex staircase — the arithmetic is finer than the engine's resolution,
          so the floor is the honest number;
        * a random term could not be measured at all (n=1 has no line-to-line scatter),
          so there is no arithmetic to be below anything — the floor is a STAND-IN and
          must be reported as one.

        What this can never return is 0.0. That value is now unreachable by construction,
        which is the point: the defect it replaces was not a wrong formula, it was a
        `None` spent as a zero one layer up.
        """
        m = self.measured_statistical()
        return QUANTISER_FLOOR_DEX if m is None else max(m, QUANTISER_FLOOR_DEX)

    def stat_basis(self) -> str:
        """Where the published statistical bar came from — a field, not a footnote.

        Written into the product artifact so a page can never render a number whose
        origin is unstated (RYA-907 §4.3, option (a)).
        """
        unmeasured = [t.name for t in self.unmeasured_terms() if t.averages_down]
        m = self.measured_statistical()
        if m is None:
            return (f"quantiser-floor {QUANTISER_FLOOR_DEX} dex (RYA-771) — "
                    f"UNMEASURED: {', '.join(unmeasured) or 'no random term'} at "
                    f"n_lines={self.n_lines}; the floor STANDS IN for a scatter that was "
                    f"never determined, it is not a measurement of one")
        if m < QUANTISER_FLOOR_DEX:
            return (f"quantiser-floor {QUANTISER_FLOOR_DEX} dex (RYA-771) — measured "
                    f"random RMS {m:.5f} dex sits below one 0.01 dex synthesis tread, "
                    f"so the engine cannot resolve it")
        return f"measured — RMS of the random terms, {m:.5f} dex at n_lines={self.n_lines}"

    def systematic(self) -> float:
        ss = [t.contribution(self.n_lines) for t in self.terms
              if not t.averages_down and t.measured]
        return math.sqrt(sum(s * s for s in ss)) if ss else 0.0

    def total(self) -> tuple[float, float]:
        """(statistical, systematic). Deliberately a PAIR.

        There is no `combined()` here. Collapsing the two is how a frontier number comes
        to look like an optical one: the statistical part shrinks with more lines and the
        systematic part does not, so a single figure invites exactly the wrong inference
        about what another night of observing would buy.
        """
        return self.statistical(), self.systematic()

    def dominant(self) -> Term | None:
        """The term that actually limits this measurement -- what to fix first.

        Ranks the MEASURED terms only. An unmeasured term has no size to compare, and
        calling it dominant-or-not would be the same guess this module refuses to make.
        """
        measured = [t for t in self.terms if t.measured]
        if not measured:
            return None
        return max(measured, key=lambda t: t.contribution(self.n_lines))

    def describe(self) -> str:
        stat, sys_ = self.total()
        lines = [f"{self.element} · {self.band} · n={self.n_lines}",
                 f"  statistical  {stat:.4f} dex   (averages down as 1/sqrt(N))",
                 f"    basis: {self.stat_basis()}",
                 f"  systematic   {sys_:.4f} dex   (does NOT average down)"]
        measured = [t for t in self.terms if t.measured]
        for t in sorted(measured, key=lambda t: -t.contribution(self.n_lines)):
            kind = "random" if t.averages_down else "SYSTEMATIC"
            lines.append(f"    {t.contribution(self.n_lines):.4f}  {t.name:<24s} "
                         f"[{kind}] {t.source}")
        # UNMEASURED TERMS ARE LISTED, IN WORDS, NOT OMITTED. A term dropped from the
        # printout reads as a term the budget does not have; this budget HAS it and
        # could not measure it, which is a different and more useful thing to know.
        for t in self.unmeasured_terms():
            kind = "random" if t.averages_down else "SYSTEMATIC"
            lines.append(f"    UNMEASURED  {t.name:<24s} [{kind}] {t.source}")
        d = self.dominant()
        if d is not None:
            lines.append(f"  dominant: {d.name} -- "
                         + ("collect more lines" if d.averages_down
                            else "more lines will NOT help; fix the source"))
        return "\n".join(lines)


# ── standard terms, each sourced ─────────────────────────────────────────────

# An ungraded Kurucz semi-empirical log gf carries 0.1-0.3 dex (RYA-161). Taken at the
# geometric middle of that range as a SCALE offset: the random part is captured by the
# observed line-to-line scatter, so this term is the part that does not average down.
UNGRADED_GF_SYSTEMATIC_DEX = 0.17

# NIST-graded gf, worst grade we accept into a product (B = 10% => log10(1.10)).
GRADED_GF_SYSTEMATIC_DEX = 0.041


def cited_gf_term(sigma_dex: float, *, n_lines: int, source: str) -> Term:
    """The gf term from the pool's OWN cited laboratory sigmas (RYA-850).

    `GRADED_GF_SYSTEMATIC_DEX` is a BOUND -- the worst grade we accept -- and a bound is
    the right answer only while the actual sigmas are unknown. When the pool's lines carry
    published per-line uncertainties, those are a MEASUREMENT of the same quantity and
    supersede the bound. For the Fe I pools this is larger, not smaller (0.052-0.060 vs
    0.041), so the bound was optimistic; a pool of grade-A lines would legitimately come
    out below it, which is why nothing here clamps to the generic value.

    Combined as an RMS by the caller rather than a median: these enter the budget in
    quadrature, and the median discards exactly the tail a quadrature sum is sensitive to.

    `source` must name the papers, never "assumed" -- and per RYA-853 the Fe I table it
    reads was refereed line-by-line against those papers before this term was wired.
    """
    if not math.isfinite(sigma_dex) or sigma_dex <= 0:
        raise ValueError(f"cited gf sigma must be finite and positive, got {sigma_dex!r}")
    if n_lines < 1:
        raise ValueError(f"cited gf sigma needs at least one line, got n_lines={n_lines}")
    return Term("gf scale (cited lab)", sigma_dex, False,
                f"RMS of the published per-line laboratory sigma over {n_lines} lines "
                f"({source}); measured, so it supersedes the generic graded bound")


def empirical_gf_term(sigma_dex: float, *, n_lines: int, provenance: str) -> Term:
    """The gf term from RYA-968's PER-LINE empirical sigmas.

    `gf_term(graded=...)` is a two-branch switch and cannot express a mixed pool: it returns the
    0.17 blanket the moment one line is ungraded, which is why RYA-855 moved 0 of 36 bars. This
    takes the RMS of the per-line sigma over the lines that actually contribute, so a pool of 7
    laboratory lines among 257 is charged for what its lines are, not for its worst member.

    RMS rather than median, for the reason `cited_gf_term` already gives: these combine in
    quadrature and a median discards the tail a quadrature sum is most sensitive to.

    ⚠️ THIS IS NOT THE WHOLE ABSOLUTE ERROR. It is the gf term only. The absolute zero-point is
    set by the laboratory anchor and does NOT shrink as ungraded lines are admitted -- see
    `zero_point_term`. Reporting this one alone is the manufactured precision RYA-968 §3.2
    measured (an admitted-pool sem of 0.005 dex against a 0.059 dex anchor).
    """
    if not math.isfinite(sigma_dex) or sigma_dex <= 0:
        raise ValueError(f"empirical gf sigma must be finite and positive, got {sigma_dex!r}")
    if n_lines < 1:
        raise ValueError(f"empirical gf sigma needs at least one line, got n_lines={n_lines}")
    if not provenance:
        raise ValueError("empirical gf term must state how each line's sigma was obtained "
                         "(cited / self-reported / inferred / fallback) -- never 'assumed'")
    return Term("gf scale (empirical, per-line)", sigma_dex, False,
                f"RMS of the per-line empirical gf sigma over {n_lines} contributing lines "
                f"({provenance}); RYA-968")


def zero_point_term(zero_point_dex: float, *, n_anchor: int) -> Term:
    """The absolute scale's own floor, from the laboratory anchor alone.

    🔴 IT DOES NOT AVERAGE DOWN, AND IT DOES NOT SHRINK WHEN UNGRADED LINES ARE ADMITTED.
    Admitting lines because they behave like the anchor produces a pool whose scatter is smaller
    than the anchor's own BY CONSTRUCTION -- but those lines carry no independent information
    about the SCALE, which the anchor alone sets. On the RYA-959 VIS cell that is 0.157/sqrt(7)
    = 0.059 dex against an admitted-pool sem of 0.005: a factor of twelve of precision that
    would be manufactured rather than measured.
    """
    if not math.isfinite(zero_point_dex) or zero_point_dex <= 0:
        raise ValueError(f"zero-point must be finite and positive, got {zero_point_dex!r}")
    if n_anchor < 1:
        raise ValueError(f"zero-point needs an anchor, got n_anchor={n_anchor}")
    return Term("absolute zero-point (laboratory anchor)", zero_point_dex, False,
                f"sigma/sqrt(n) over {n_anchor} laboratory-graded anchor lines; the absolute "
                f"scale rests on these alone and admitting ungraded lines cannot improve it "
                f"(RYA-968)")


def gf_term(*, graded: bool) -> Term:
    if graded:
        return Term("gf scale (NIST-graded)", GRADED_GF_SYSTEMATIC_DEX, False,
                    "NIST grade B = 10% on the transition probability = log10(1.10) dex; "
                    "a BOUND -- superseded by cited_gf_term where per-line sigmas exist")
    return Term("gf scale (UNGRADED)", UNGRADED_GF_SYSTEMATIC_DEX, False,
                "ungraded Kurucz semi-empirical loggf, 0.1-0.3 dex (RYA-161); the random "
                "part shows up in the line-to-line scatter, this is the scale offset that "
                "does not")


def harness_term(measured_residual_dex: float, handler: str,
                 provenance: str = "") -> Term:
    """The handler's own systematic, described BY HOW IT WAS OBTAINED — RYA-873.

    🔴 THIS FUNCTION USED TO ASSERT. It printed "{handler} MEASURED against the known
    optical answer, not assumed zero (control/frontier rule)" beside WHATEVER value it
    was handed — including `SynthesisHandler`'s 0.0000, which no control established and
    which that sentence explicitly disclaims. Every synthesis-route budget in the repo
    carried that contradiction: the prose said measured, the arithmetic was an assumption.

    So the sentence is now a function of the provenance rather than a constant. An
    unstated provenance says it is unstated; it does not inherit the claim.
    """
    from pipeline.harness_residual import (
        PROV_MEASURED, PROV_UNCHARGED, PROV_UNSTATED, uncharged_note)
    prov = (provenance or PROV_UNSTATED).strip()
    if prov == PROV_UNCHARGED:
        why = uncharged_note(handler) or "no residual has been established for it"
        src = (f"{handler}: NO RESIDUAL IS CHARGED — {why}. This is an ABSENCE, not a "
               f"measured zero, and it is carried at 0.0 so the bar is not inflated by "
               f"a number nobody has")
    elif prov == PROV_MEASURED:
        src = (f"{handler} MEASURED against the known optical answer, not assumed zero "
               f"(control/frontier rule)")
    else:
        src = (f"{handler}: the provenance of this residual is NOT STATED by the caller, "
               f"so it is not described as measured (RYA-873)")
    return Term("harness residual", measured_residual_dex, False, src)


def scatter_term(observed_scatter_dex: float | None, n_lines: int | None = None) -> Term:
    """The line-to-line scatter — or the honest statement that there is none to measure.

    `None` is what `band_products.build_product` returns when fewer than two lines entered
    the aggregate, and it is CORRECT: the spread of one number is not a small spread, it
    is an undefined quantity. Passing that `None` through as `None` is the entire RYA-907
    fix; the bug was one call site writing `(product.sigma or 0.0)`.
    """
    if observed_scatter_dex is None:
        n = "" if n_lines is None else f" (n_lines={n_lines})"
        return Term("line-to-line scatter", None, True,
                    f"NOT MEASURABLE{n}: the spread of a single accepted line is "
                    f"undefined, not small. `band_products.build_product` returns None "
                    f"here and that None is carried, never spent as a zero (RYA-907)")
    return Term("line-to-line scatter", observed_scatter_dex, True,
                "observed spread of the accepted lines")


def assert_stat_publishable(stat_dex: float, *, cell: str) -> float:
    """Refuse to let a statistical bar reach an artifact unless it can be defended.

    THE GUARD THE DEFECT NEEDED. Two Fe II red-optical ENGINE-A cells shipped
    `stat_dex = 0.0` and nothing objected — the CSV simply carried it, the published
    total became `hypot(0.0, 0.1731)`, and the systematic alone was presented as though
    the statistical half had been measured and found negligible.

    Called at every emit site with the value ACTUALLY ABOUT TO BE WRITTEN, not with the
    budget it came from: a guard that re-derives its own subject cannot catch a driver
    that writes something else. Returns the value so it can wrap an assignment.
    """
    if stat_dex is None or not math.isfinite(stat_dex):
        raise UnpublishableStat(
            f"{cell}: statistical uncertainty is {stat_dex!r}. A product may not be "
            f"emitted without one (RYA-907).")
    if stat_dex <= 0.0:
        raise UnpublishableStat(
            f"{cell}: statistical uncertainty is {stat_dex}. A published 0.0 claims the "
            f"quantity was measured and found negligible; for an n=1 product nobody "
            f"measured it at all. Carry the term as UNMEASURED and publish the "
            f"{QUANTISER_FLOOR_DEX} dex quantiser floor (RYA-771/907).")
    if stat_dex < QUANTISER_FLOOR_DEX:
        raise UnpublishableStat(
            f"{cell}: statistical uncertainty {stat_dex} is below the {QUANTISER_FLOOR_DEX} "
            f"dex quantiser floor. iSpec writes trial abundances with %.2f, so EW(A) is a "
            f"staircase and no bar can honestly be finer than one tread (RYA-771).")
    return stat_dex


def build(element: str, wavelength_A: float, n_lines: int, *,
          scatter_dex: float | None, gf_graded: bool, harness_residual_dex: float,
          handler: str, harness_provenance: str = "",
          cited_gf_sigma_dex: float | None = None,
          cited_gf_source: str = "") -> ErrorBudget:
    """Assemble the budget for a band, adding the terms that band actually has.

    `cited_gf_sigma_dex` supplies the pool's own published per-line gf uncertainty and
    REPLACES the generic graded bound when present (RYA-850). It is only meaningful for a
    graded pool -- an ungraded Kurucz line has no cited sigma to average -- so passing it
    with `gf_graded=False` is a caller error rather than a silent preference.
    """
    pol = resolve(wavelength_A)
    b = ErrorBudget(element=element, band=pol.name, n_lines=n_lines)
    b.add(scatter_term(scatter_dex, n_lines))
    if cited_gf_sigma_dex is not None:
        if not gf_graded:
            raise ValueError(
                "cited_gf_sigma_dex is only defined for a graded pool; an ungraded "
                "Kurucz line has no published per-line sigma (RYA-850)")
        if not cited_gf_source:
            raise ValueError("cited_gf_sigma_dex requires cited_gf_source naming the "
                             "papers -- a budget term is never unsourced")
        b.add(cited_gf_term(cited_gf_sigma_dex, n_lines=n_lines,
                            source=cited_gf_source))
    else:
        b.add(gf_term(graded=gf_graded))
    b.add(harness_term(harness_residual_dex, handler, harness_provenance))

    if "pseudo" in pol.continuum_treatment.lower():
        # The size is set by how far the reachable envelope sits below unity: near-UV
        # median flux 0.283-0.805 means the normalisation itself is uncertain at the
        # ~10% level in the worst windows, and no number of lines fixes that.
        b.add(Term("pseudo-continuum", 0.10, False,
                   f"{pol.name}: the true continuum is never observed (median flux "
                   f"0.283-0.805); the envelope stands in for it and the residual does "
                   f"not average down"))
    if pol.telluric_required:
        b.add(Term("telluric residual", 0.03, False,
                   f"{pol.name}: telluric correction residuals are epoch- and "
                   f"airmass-dependent, so they vary between observations of the SAME "
                   f"star and cannot be calibrated once"))
    return b


# ── THE ZERO-POINT CAP, DECOMPOSED — RYA-987 ─────────────────────────────────
#
# `zero_point_term` above charges sigma/sqrt(n) over the anchor. That is the STATISTICAL
# half only, and while the anchor was 7 lines it was also the whole story: 0.157/sqrt(7) =
# 0.059 dex swamped everything else. RYA-984 grew the anchor to 163, the statistical term
# fell by ~4x, and a second term that had been hiding under it is now the larger one.
#
# 🔴 THE TWO TERMS ANSWER DIFFERENT QUESTIONS AND MUST NOT BE BLENDED INTO ONE NUMBER.
#
#   STATISTICAL  scatter/sqrt(n) — how well the anchor's MEAN is determined. Shrinks with
#                every line added.
#   SYSTEMATIC   the laboratory gf SCALE's own zero point — a shift shared by lines that
#                come from the same measurement. Adding lines from that source cannot
#                reduce it, because they all carry the same shift.
#
# WHY THE SYSTEMATIC IS SOURCE-CORRELATED AND NOT PER-LINE RANDOM. A published per-line
# sigma covers that line's own measurement, but a branching-fraction x lifetime experiment
# also carries a normalisation common to every line it reports. Treating the cited sigmas
# as independent would let them average down as 1/sqrt(n) — which is exactly the
# manufactured precision RYA-968 §3.2 exists to prevent. So they are combined as
# CORRELATED WITHIN A SOURCE and independent ACROSS sources.
#
# ⚠️ THE BLIND SPOT, STATED. An offset COMMON TO ALL THREE LABORATORIES is invisible from
# inside the anchor: every internal comparison differences it away. Detecting it needs an
# external absolute reference, and the RYA-161 firewall forbids using the gold value for
# that — gold carries its own hidden zero point, so "validating" against it would launder
# one unknown into another. The number below is therefore a floor on the systematic, not a
# bound on it, and that is a property of the measurement, not a defect in the code.

@dataclass(frozen=True)
class ZeroPointCap:
    """The absolute-accuracy floor of an abundance scale, split into its two terms."""
    n_anchor: int
    scatter_dex: float
    statistical_dex: float
    systematic_dex: float
    per_source: tuple[tuple[str, int, float], ...]   # (source, n, cited sigma RMS)

    @property
    def combined_dex(self) -> float:
        return math.hypot(self.statistical_dex, self.systematic_dex)

    @property
    def limited_by(self) -> str:
        return "systematic" if self.systematic_dex > self.statistical_dex else "statistical"

    @property
    def variance_share_systematic(self) -> float:
        return self.systematic_dex ** 2 / (self.combined_dex ** 2)

    def lines_to_halve(self) -> float | None:
        """How many anchor lines would halve the cap — `None` when no n can.

        The honest way to say "more lines will not help". Once the systematic exceeds half
        the combined cap, no finite anchor reaches half of it, because the systematic is
        the floor the statistical term is added to in quadrature.
        """
        target = self.combined_dex / 2.0
        if self.systematic_dex >= target:
            return None
        need = math.sqrt(target ** 2 - self.systematic_dex ** 2)
        return self.n_anchor * (self.statistical_dex / need) ** 2

    def describe(self) -> str:
        src = "; ".join(f"{s} n={n} sigma_RMS={g:.4f}" for s, n, g in self.per_source)
        cap = self.combined_dex
        lim = self.limited_by
        need = self.lines_to_halve()
        lever = ("no anchor size halves this cap — the laboratory gf scale is the floor, "
                 "so BETTER gf is the only remaining lever"
                 if need is None else f"~{need:.0f} anchor lines would halve it")
        return (f"zero-point cap {cap:.4f} dex on {self.n_anchor} graded anchor lines "
                f"[{lim.upper()}-LIMITED, systematic carries "
                f"{self.variance_share_systematic:.0%} of the variance]\n"
                f"  statistical {self.statistical_dex:.4f} = scatter {self.scatter_dex:.4f}"
                f" / sqrt({self.n_anchor})\n"
                f"  systematic  {self.systematic_dex:.4f} = laboratory gf zero point, "
                f"correlated within a source ({src})\n"
                f"  {lever}")


def zero_point_cap(abundances, sources, cited_sigmas) -> ZeroPointCap:
    """Decompose an anchor's absolute floor into statistical + laboratory-gf zero point.

    `abundances`, `sources` and `cited_sigmas` are per anchor line and must be parallel.

    🔴 NO REFERENCE VALUE APPEARS HERE, BY CONSTRUCTION. There is no parameter for one.
    The cap is computed from the anchor's OWN scatter and its OWN cited laboratory sigmas,
    never from proximity to a published solar abundance (RYA-161). A function that took a
    reference could be asked "how close are we?", and that question answered against a
    number carrying its own zero point is how a scale gets talked into a precision it does
    not have.
    """
    import numpy as _np
    a = _np.asarray(abundances, dtype=float)
    ok = _np.isfinite(a)
    a = a[ok]
    src = _np.asarray(list(sources), dtype=object)[ok]
    sig = _np.asarray(cited_sigmas, dtype=float)[ok]
    n = a.size
    if n < 2:
        raise ValueError(f"a zero-point cap needs at least two anchor lines, got {n}")
    if not _np.isfinite(sig).all():
        raise ValueError("every anchor line must carry a cited laboratory sigma — a "
                         "missing one is not a zero, and silently treating it as one "
                         "would shrink the systematic (RYA-873)")

    scatter = float(a.std(ddof=1))
    statistical = scatter / math.sqrt(n)

    per_source, var = [], 0.0
    for s in sorted({str(x) for x in src}):
        m = _np.array([str(x) == s for x in src])
        rms = float(_np.sqrt((sig[m] ** 2).mean()))
        frac = float(m.sum()) / n
        per_source.append((s, int(m.sum()), rms))
        # Correlated within the source: its sigma enters weighted by its share of the
        # anchor and does NOT get a 1/sqrt(n_s). Independent across sources: quadrature.
        var += (frac * rms) ** 2
    return ZeroPointCap(n_anchor=n, scatter_dex=scatter, statistical_dex=statistical,
                        systematic_dex=math.sqrt(var), per_source=tuple(per_source))


def _main(argv=None) -> int:
    """`python3 -m pipeline.error_budget --anchor <name> --zero-point-cap --decompose`.

    Re-runnable BY DESIGN: the anchor is still growing (the HARPS deep leg is owed, RYA-977
    adds ~65 lines, the UV/IR bands are not folded in), so the cap is a snapshot keyed to a
    named pool rather than a constant anyone can quote out of date.
    """
    import argparse
    ap = argparse.ArgumentParser(description="RYA-987 zero-point cap")
    ap.add_argument("--anchor", required=True)
    ap.add_argument("--zero-point-cap", action="store_true")
    ap.add_argument("--decompose", action="store_true",
                    help="report the two terms separately. RYA-968 §3.2: a single blended "
                         "number cannot say whether more lines would help.")
    a = ap.parse_args(argv)
    if not a.zero_point_cap:
        ap.error("nothing to do; pass --zero-point-cap")

    from pipeline.anchor_pools import ANCHORS, load
    pool = ANCHORS[a.anchor] if a.anchor in ANCHORS else None
    if pool is None:
        raise SystemExit(f"unknown anchor {a.anchor!r}; known: {sorted(ANCHORS)}")
    df = load(a.anchor)
    cap = zero_point_cap(df.abundance, df.lab_source, df.cited_sigma_dex)

    print(f"\nANCHOR {a.anchor}  ({pool.species}, n={cap.n_anchor})")
    print(f"  {pool.note}\n")
    print(cap.describe())
    if a.decompose:
        print("\n  per-source (correlated within, independent across):")
        for s, n, g in cap.per_source:
            print(f"    {s:8s} n={n:4d}  cited sigma RMS {g:.4f}  "
                  f"share {n / cap.n_anchor:5.1%}  contributes "
                  f"{(n / cap.n_anchor) * g:.4f} dex")
        print(f"\n  statistical {cap.statistical_dex:.4f}   "
              f"systematic {cap.systematic_dex:.4f}   "
              f"CAP {cap.combined_dex:.4f} dex")
    print("\n  🔒 no reference abundance enters this computation (RYA-161): the cap is "
          "the anchor's own scatter and its own cited laboratory sigmas.")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
