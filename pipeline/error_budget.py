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
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

from pipeline.band_policy import resolve


@dataclass(frozen=True)
class Term:
    """One contribution to the budget, and how it behaves under more lines."""
    name: str
    dex: float
    averages_down: bool     # True => random, scales as 1/sqrt(N); False => systematic floor
    source: str             # where the number comes from -- never "assumed"

    def contribution(self, n_lines: int) -> float:
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
    def statistical(self) -> float:
        rs = [t.contribution(self.n_lines) for t in self.terms if t.averages_down]
        return math.sqrt(sum(r * r for r in rs)) if rs else 0.0

    def systematic(self) -> float:
        ss = [t.contribution(self.n_lines) for t in self.terms if not t.averages_down]
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
        """The term that actually limits this measurement -- what to fix first."""
        if not self.terms:
            return None
        return max(self.terms, key=lambda t: t.contribution(self.n_lines))

    def describe(self) -> str:
        stat, sys_ = self.total()
        lines = [f"{self.element} · {self.band} · n={self.n_lines}",
                 f"  statistical  {stat:.4f} dex   (averages down as 1/sqrt(N))",
                 f"  systematic   {sys_:.4f} dex   (does NOT average down)"]
        for t in sorted(self.terms, key=lambda t: -t.contribution(self.n_lines)):
            kind = "random" if t.averages_down else "SYSTEMATIC"
            lines.append(f"    {t.contribution(self.n_lines):.4f}  {t.name:<24s} "
                         f"[{kind}] {t.source}")
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


def scatter_term(observed_scatter_dex: float) -> Term:
    return Term("line-to-line scatter", observed_scatter_dex, True,
                "observed spread of the accepted lines")


def build(element: str, wavelength_A: float, n_lines: int, *,
          scatter_dex: float, gf_graded: bool, harness_residual_dex: float,
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
    b.add(scatter_term(scatter_dex))
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
