"""
RYA-968 — the empirical gf-grade framework: per-LINE tier, reason-code and sigma.

WHAT THIS REPLACES. `gf_rung.decide()` grades a POOL all-or-nothing: rung 1 (the 0.17 blanket)
the moment one line in it is ungraded. That rule is CORRECT -- a pool must not inherit a
laboratory pedigree from a subset -- and it is also why RYA-855 moved 0 of 36 bars, because
every real Fe cell is mixed (VIS: 7 GF-LAB among 257). The fix is not to relax it. **This module
grades LINES**, so a mixed pool becomes a quadrature over per-line sigma instead of a category
that collapses to its worst member.

THE METHOD IS THE FIELD'S, NOT OURS (design spec §1). The GBS third version (Soubiran+2024) and
OCCASO (Carbajo-Hijarrubia+2024) use the laboratory-graded lines as a REFERENCE DISTRIBUTION and
ADMIT ungraded lines that behave like it. The stage structure below is Elgueta+2026's physical
decision tree -- Depth, Saturation, Purity, goodness-of-fit -- carrying Gaia-ESO/Heiter+2021's
two axes (`gf_flag` and `synflag` are separate questions and one tier cannot answer both).

🔴 THE FIREWALL IS STRUCTURAL, NOT ASPIRATIONAL (RYA-161).
  F1  No external reference is an input. Nothing here imports a gold/reference abundance, and
      `test_gf_empirical_rya968` greps this module to keep it that way. The anchor is built from
      OUR OWN laboratory-graded lines, in the same spectrum.
  F2  Stages 1-5 are INVARIANT UNDER A CONSTANT OFFSET: add delta to every abundance in a cell
      and every verdict and every sigma comes out bit-identical. Scatter is offset-invariant,
      agreement is not, and this is the property that makes "we graded on precision" CHECKABLE.
      ⚠️ STAGE 6 IS DELIBERATELY *NOT* OFFSET-INVARIANT -- that is what it is for -- so it owes a
      per-line audit list, which `stage6_anchor_consistency` returns rather than hides.
  F3  Every threshold is DECLARED, and undeclared means LOUD-FAIL, never a guessed default. See
      `Thresholds`.
  F4  The Cr canary (+0.402) is a blocking regression test, asserted in the test module.

⚠️ WHY THERE ARE NO DEFAULT THRESHOLD VALUES. RYA-968 §3.1 measured what a BORROWED constant does
here: the GBS window -6.7 < REW < -4.5 admits 249 of our 249 VIS lines -- it controls nothing --
and the GBS +/-0.05 dex admission tolerance rejects 6 of our own 7 anchor lines. A borrowed number
is not a control. `Thresholds` therefore has no defaults and raises until each value is set.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Iterable, Sequence

import numpy as np

# ── vocabulary ───────────────────────────────────────────────────────────────────────
#: Tiers. GOLDEN never moves off laboratory gf; the bottom two are DOCUMENTED, never deleted
#: (RYA-711) and render through RYA-707/224/851.
TIER_GOLDEN = "GOLDEN"
TIER_CONSISTENT = "UNGRADED-CONSISTENT"
TIER_SCATTERED = "UNGRADED-SCATTERED"
TIER_INVALID = "INVALID"

#: 🔴 REASON CODES ARE RYA-809's, NOT A PARALLEL VOCABULARY. These are the values already used in
#: `data/registry/problem_children.csv`; inventing a second set here would split the taxonomy the
#: appendix renderers already know how to display.
RC_BAD_GF = "BAD_GF"
RC_ATOMIC_BLEND = "ATOMIC_BLEND"
RC_MOLECULAR_BLEND = "MOLECULAR_BLEND"
RC_TELLURIC = "TELLURIC_ADJACENT"
RC_SATURATION = "SATURATION_COG"
RC_CONTINUUM = "CONTINUUM_LIMITED"
RC_NON_MINIMUM = "NON_MINIMUM"
RC_OUTLIER = "ABUNDANCE_OUTLIER"
RC_BAD_FIT = "BAD_FIT"
RC_OK = ""

#: Stage verdicts, Elgueta+2026's alphabet.
YES, NO, UNKNOWN = "Y", "N", "U"

#: gf tiers as `canonical_gf.gf_tier` writes them (RYA-945).
LAB_TIERS = frozenset({"LAB"})
GRADED_TIERS = frozenset({"LAB", "NIST-C+"})


class ThresholdNotDeclared(RuntimeError):
    """A threshold was read before it was set. RYA-968 F3: declared, or nothing."""


@dataclass
class Thresholds:
    """Every number this framework needs, with NO DEFAULTS.

    🔴 A BORROWED CONSTANT IS NOT A CONTROL (design spec §3.1). Each field starts as None and
    `require()` raises with the ticket reference rather than substituting a published value that
    was tuned on somebody else's distribution.
    """
    #: weak-end noise floor: below this a line is not measurable, which is a different
    #: question from saturation. DERIVED from our own EW distribution -- not GBS's -6.7,
    #: which admits 100% of our pool and therefore controls nothing.
    rew_min: float | None = None
    #: 🔴 `rew_max` IS RETIRED (RYA-1043). Saturation is a PER-LINE property: the ensemble
    #: COG for solar Fe I has no resolvable knee, so no scalar can express it. Stage 2 now
    #: gates on the line's OWN measured dREW/dA against the flattening threshold below.
    #: Keeping a scalar here would invite someone to set it to Gray's -4.9 again.
    min_d_rew_dA: float | None = None
    #: minimum core depth for stage 1
    min_depth: float | None = None
    #: reduced chi2 above which stage 4 fails — 🔴 EW ROUTE ONLY (RYA-1041). The synth
    #: route gates on RYA-992's per-arm `frac_rise_weaker` instead, because the graded
    #: anchor's own synth red_chi2 runs to 3404 and any cap here would delete it.
    max_red_chi2: float | None = None
    #: |A - anchor_mean| tolerance for stage 6 admission.
    #:
    #: 🔴 SET AT 0.0119 dex — the Fe anchor's own DIFFERENTIAL line-to-line scatter,
    #: measured 2026-08-25 (RYA-1041 D2, artifact
    #: `data/audit/rya1041/rya1041_threshold_diagnostics.json`): same-instrument pair
    #: KP1984-molecfit vs KP2005, 67 lines paired per line, robust sigma 1.4826*MAD =
    #: 0.0119, median offset +0.0080.
    #:
    #: DIFFERENTIAL, not solo, and the distinction is an order of magnitude: solo A(Fe)
    #: scatter on the same pool is 0.1290 dex. Every published value is
    #: [X/H] = A(X)star - A(X)sun and works because the systematics CANCEL on a shared
    #: chain (RYA-898); sizing this from solo scatter would admit lines on a spread the
    #: differential never pays. GBS's 0.05 admits only 28.3% of our own laboratory anchor.
    #:
    #: SAME-INSTRUMENT is the right comparison because the product is per
    #: (instrument x band) (RYA-1026), so a line is graded against its OWN cell's anchor
    #: and no cross-instrument ambiguity arises. Cross-instrument pairs run 0.0163
    #: (IAG-KP) to 0.0667 (HARPS-KP) — the HARPS figure being the arm that does not
    #: cancel, which is a finding about HARPS, not a tolerance.
    #:
    #: SET AT IT, NEVER TIGHTER (RYA-1041/161): tighter selects for agreement with our own
    #: median, which manufactures precision the anchor cannot support.
    admission_tol_dex: float | None = 0.0119
    #: below this the anchor cannot define a distribution and stage 6 must not run at all.
    min_anchor_lines: int | None = None
    #: |sigma| above which an admitted line is SCATTERED rather than CONSISTENT.
    consistency_bound_dex: float | None = None
    #: last-resort per-line sigma. NOT Kurucz's 0.17 -- our own pooled measurement.
    fallback_sigma_dex: float | None = None

    def require(self, name: str):
        v = getattr(self, name)
        if v is None:
            raise ThresholdNotDeclared(
                f"threshold {name!r} has not been declared. RYA-968 F3 requires every threshold "
                f"to be fixed in config BEFORE the data is seen, and §3.1 measured what happens "
                f"when a published constant is borrowed instead: the GBS REW window admits 100% "
                f"of our pool and the GBS admission tolerance rejects 6 of our own 7 anchor "
                f"lines. Set it deliberately; there is no default.")
        return v


@dataclass
class LineVerdict:
    """One line's verdict. `stages` is the Elgueta Y/N/U tree, in order."""
    wavelength_air_A: float
    stages: dict = field(default_factory=dict)
    tier: str = ""
    reason_code: str = RC_OK
    sigma_dex: float | None = None
    sigma_source: str = ""
    note: str = ""
    #: stage-6 audit (F2): the distance that decided admission, or None if stage 6 did not run.
    anchor_distance_dex: float | None = None

    @property
    def in_value(self) -> bool:
        """Does this line contribute to the reported abundance?"""
        return self.tier in (TIER_GOLDEN, TIER_CONSISTENT)


# ── stages 1-4: PHYSICAL, and they never see an abundance ────────────────────────────
# 🔴 The signatures below take no abundance at all. That is F1/F2 enforced by construction
# rather than by discipline: a function that is not given the number cannot grade toward it.

def stage1_depth(observed_depth, th: Thresholds) -> tuple[str, str]:
    """Measurable against the noise?"""
    lo = th.require("min_depth")
    if observed_depth is None or not np.isfinite(observed_depth):
        return UNKNOWN, RC_OK
    return (YES, RC_OK) if observed_depth >= lo else (NO, RC_CONTINUUM)


def stage2_saturation(rew, th: Thresholds, *,
                      d_rew_dA: float | None = None) -> tuple[str, str]:
    """Is THIS line still on the linear part of its own curve of growth? — RYA-1043.

    🔴 THE SATURATION CEILING IS NOT A SCALAR, and that is measured rather than argued.
    An ensemble COG for solar Fe I has no resolvable knee: the clean-and-untruncated
    1805-line intersection IS a curve of growth (corr +0.362, medians monotone) and its
    local dREW/dX still carries three sign changes, because damping, blending and core
    depth vary line to line and smear the transition. A single-line COG is sharp; an
    ensemble COG is broad BY CONSTRUCTION (RYA-1041/1043).

    So every line saturates at its OWN reduced EW, and the gate is not "is REW below the
    one knee" but "has THIS line's own response flattened". `d_rew_dA` is that response,
    measured per line by `scripts/rya1041_perline_cog.py`: 1 on the linear part where EW
    is proportional to the number of absorbers, collapsing toward 0 as the core saturates.
    Measured over 617 VIS lines it runs 1.044 down to -0.129 and falls monotonically with
    line strength (corr -0.807).

    The scalar `rew_max` is therefore GONE. `rew_min` survives as the weak-end noise floor,
    which is a different question -- can we measure this line at all -- and is not what a
    curve of growth answers.

    A line with no measured derivative is UNKNOWN, never a pass: the per-line COG is
    expensive and incomplete, and treating "not yet computed" as "linear" would admit
    exactly the saturated lines this stage exists to catch (RYA-833).
    """
    lo = th.require("rew_min")
    if rew is not None and np.isfinite(rew) and rew <= lo:
        return NO, RC_SATURATION            # too weak to measure — the noise floor
    if d_rew_dA is None or not np.isfinite(d_rew_dA):
        return UNKNOWN, RC_OK
    return ((YES, RC_OK) if float(d_rew_dA) >= th.require("min_d_rew_dA")
            else (NO, RC_SATURATION))


def stage3_purity(problem_class, excluded_reason="") -> tuple[str, str]:
    """Unblended -- Gaia-ESO's `synflag` axis.

    Reads the standing per-line registry (RYA-463/809) rather than re-deriving blending, so one
    taxonomy governs and a line already adjudicated as blended cannot be silently re-admitted.
    """
    pc = str(problem_class or "").strip().upper()
    blend = {RC_ATOMIC_BLEND, RC_MOLECULAR_BLEND, RC_TELLURIC}
    if pc in blend:
        return NO, pc
    txt = str(excluded_reason or "").upper()
    for rc in (RC_TELLURIC, RC_ATOMIC_BLEND, RC_MOLECULAR_BLEND):
        if rc in txt:
            return NO, rc
    return (YES, RC_OK) if pc == "" else (UNKNOWN, RC_OK)


class RouteNotDeclared(RuntimeError):
    """Stage 4 was handed a line that does not say which route measured it — RYA-1041.

    Refused rather than defaulted. The two routes disagree about `red_chi2` by five and a
    half orders of magnitude (RYA-981; measured again at band scale in RYA-1041: EW median
    2.2e-4 against synth 67), so a line that falls into the wrong branch is not slightly
    mis-graded -- a synth line entering the EW branch is deleted by a cap it was never
    meant to face.
    """


def stage4_fit(red_chi2, ew_mA, observed_depth, wavelength_A, th: Thresholds,
               fit_window_A: float | None = None, *,
               route: str | None = None, instrument: str | None = None,
               frac_rise_weaker: float | None = None) -> tuple[str, str]:
    """Does the profile actually describe the feature? ROUTE-AWARE — RYA-1041.

    ⚠️ DELEGATES to `pipeline.line_width` rather than re-deriving the ceilings. RYA-959 built
    those checks in one module precisely so both fitters -- and now this -- resolve to the same
    functions; a second copy here is how RYA-906/911's floor reached only one of two sites.

    🔴 `red_chi2` IS DISQUALIFIED AS A SYNTH GATE, and this is measured rather than argued.
    The GRADED laboratory anchor -- lines whose gf is KNOWN-good -- has synth `red_chi2`
    median 67 and max 3404 (RYA-1041 D1). Any conventional cap deletes the reference
    itself. RYA-981 measured why: corr(red_chi2, |A - meanA|) = 0.025, i.e. it does not
    predict abundance error on this route at all. So the synth branch gates on RYA-992's
    `frac_rise_weaker` at its PER-ARM cut instead -- a statement about the objective (chi2
    doubling as A is driven weaker), not a quota -- and never touches the cap.

    `max_red_chi2` is therefore an EW-ROUTE-ONLY threshold. It is not consulted on synth.
    """
    from pipeline import line_width
    rt = str(route or "").strip().upper()
    if not rt:
        raise RouteNotDeclared(
            "stage4_fit needs the line's route ('ew' or 'synth'). It is not defaulted: "
            "a synth line falling into the EW branch meets a red_chi2 cap that would "
            "delete the laboratory anchor itself (RYA-1041/981).")
    if ew_mA is not None and observed_depth is not None \
            and np.isfinite(ew_mA) and np.isfinite(observed_depth) and observed_depth > 0:
        try:
            if line_width.implied_width_exceeds_ceiling(
                    float(ew_mA), float(observed_depth), float(wavelength_A),
                    fit_window_A=fit_window_A):
                return NO, RC_BAD_FIT
        except TypeError:
            if line_width.implied_width_exceeds_ceiling(
                    float(ew_mA), float(observed_depth), float(wavelength_A)):
                return NO, RC_BAD_FIT
    if rt.startswith("SYNTH") or rt == "EW-3D":
        # RYA-992's statistic, at the cut derived for THIS arm. `synth_gof_cut` takes the
        # instrument from the caller by design (RYA-913/922: no module-level instrument
        # constant), so a missing arm is a refusal, not a fallback to Kitt Peak's 1.0.
        from pipeline.synth_gof import synth_gof_ok, synth_gof_reason
        if not instrument:
            raise RouteNotDeclared(
                "the synth branch needs the instrument: the frac_rise cut is PER-ARM "
                "(RYA-992), and defaulting it would silently apply Kitt Peak's cut to "
                "another arm.")
        if frac_rise_weaker is None or not np.isfinite(frac_rise_weaker):
            return UNKNOWN, RC_OK
        return ((YES, RC_OK) if synth_gof_ok(float(frac_rise_weaker), instrument)
                else (NO, RC_BAD_FIT))
    # EW route — unchanged. This is the only surviving use of max_red_chi2.
    cap = th.require("max_red_chi2")
    if red_chi2 is None or not np.isfinite(red_chi2):
        return UNKNOWN, RC_OK
    return (YES, RC_OK) if red_chi2 <= cap else (NO, RC_BAD_FIT)


def stage5_gf(gf_tier) -> tuple[str, str]:
    """Gaia-ESO's `gf_flag` axis, read from RYA-945's provenance ingest."""
    t = str(gf_tier or "").strip()
    if t in LAB_TIERS:
        return YES, RC_OK
    if t in GRADED_TIERS:
        return YES, RC_OK
    return (NO, RC_BAD_GF) if t else (UNKNOWN, RC_BAD_GF)


# ── the anchor ───────────────────────────────────────────────────────────────────────
@dataclass
class Anchor:
    """The laboratory-graded reference distribution, built from OUR OWN spectrum.

    🔴 THIS IS THE ONLY THING STAGE 6 MAY COMPARE AGAINST (F1). It is built from lines whose gf
    an apparatus measured, in the same cell, so "does this line behave like lines whose gf we
    know?" is a gf-QUALITY question. Comparing against a literature A(X) would be the circular
    version and is what RYA-161 forbids.
    """
    n: int
    mean: float
    sd: float
    #: 🔴 THE ZERO-POINT THE WHOLE ABSOLUTE PRODUCT RESTS ON, and it does NOT shrink when
    #: ungraded lines are admitted -- they carry no independent information about the SCALE.
    #: Measured on our VIS cell: 0.157/sqrt(7) = 0.059 dex, against an admitted-pool sem of
    #: 0.005. Quoting the latter as the absolute uncertainty is the circularity in a new costume
    #: (design spec §3.2).
    zero_point_dex: float
    cited_sigma_rms: float | None = None

    @property
    def usable(self) -> bool:
        return self.n >= 2 and np.isfinite(self.sd) and self.sd > 0


def build_anchor(abundances, th: Thresholds, cited_sigmas=None) -> Anchor:
    """The anchor, or a loud refusal if there are too few laboratory lines to define one."""
    a = np.asarray([x for x in abundances if x is not None and np.isfinite(x)], float)
    need = th.require("min_anchor_lines")
    if a.size < need:
        raise ValueError(
            f"anchor has {a.size} laboratory-graded lines, below the declared minimum {need}. "
            f"An anchor this thin cannot define an admission distribution, and the quadrature "
            f"the zero-point rests on is least stable exactly when the anchor is poorly known "
            f"(RYA-968 §3.2/§7). Grow the measured laboratory pool -- do NOT lower the minimum "
            f"to make this pass.")
    sd = float(np.std(a, ddof=1))
    cs = None
    if cited_sigmas is not None:
        c = np.asarray([x for x in cited_sigmas if x is not None and np.isfinite(x)], float)
        if c.size:
            cs = float(np.sqrt(np.mean(c ** 2)))
    return Anchor(n=int(a.size), mean=float(np.mean(a)), sd=sd,
                  zero_point_dex=float(sd / math.sqrt(a.size)), cited_sigma_rms=cs)


# ── stage 6: the ONLY stage that compares abundances ─────────────────────────────────
def stage6_anchor_consistency(abundance, anchor: Anchor, th: Thresholds) -> tuple[str, float]:
    """Does the line behave like the laboratory-anchored distribution?

    Returns (verdict, distance). ⚠️ **THE DISTANCE IS RETURNED, NOT DISCARDED.** Stage 6 is the
    one stage that is deliberately NOT offset-invariant, so F2 cannot cover it; the compensating
    control is that every line's distance is auditable. A caller that drops it has removed the
    only check on the only stage that can see an abundance.
    """
    tol = th.require("admission_tol_dex")
    if abundance is None or not np.isfinite(abundance) or not anchor.usable:
        return UNKNOWN, float("nan")
    d = float(abundance) - anchor.mean
    return (YES if abs(d) <= tol else NO), d


# ── per-line sigma, in strict precedence ─────────────────────────────────────────────
def resolve_sigma(*, cited_sigma_dex=None, cross_measurement_abundances=None,
                  inferred_sigma_dex=None, th: Thresholds) -> tuple[float, str]:
    """cited -> self-reported -> inferred -> fallback. First available wins.

    The ordering is not a preference, it is a hierarchy of evidence:
      1. CITED -- someone measured this line's gf with an apparatus and published a sigma.
      2. SELF-REPORTED -- the line's own spread across independent measurements. Needs NO anchor
         and NO confound model, because the line is compared only to itself. That is why it
         outranks anything inferred, and why it is the route that still works for elements with
         no laboratory anchor at all (design spec §8).
      3. INFERRED -- confound-controlled excess over the floor, in matched REW/EP bins.
      4. FALLBACK -- our own pooled measurement. NOT Kurucz's 0.17.
    """
    if cited_sigma_dex is not None and np.isfinite(cited_sigma_dex) and cited_sigma_dex > 0:
        return float(cited_sigma_dex), "cited"
    if cross_measurement_abundances is not None:
        a = np.asarray([x for x in cross_measurement_abundances
                        if x is not None and np.isfinite(x)], float)
        if a.size >= 3:
            return float(np.std(a, ddof=1)), f"self-reported (n={a.size})"
    if inferred_sigma_dex is not None and np.isfinite(inferred_sigma_dex) \
            and inferred_sigma_dex > 0:
        return float(inferred_sigma_dex), "inferred"
    return float(th.require("fallback_sigma_dex")), "fallback"


# ── the emitter ──────────────────────────────────────────────────────────────────────
def grade_line(row, th: Thresholds, anchor: Anchor | None = None) -> LineVerdict:
    """One line -> tier + reason-code + sigma. `row` is a mapping of the per-line columns.

    🔴 `row` MUST NOT carry a reference abundance and this function never asks for one. The only
    abundance that enters is the line's OWN, and only at stage 6.
    """
    g = (lambda k, d=None: row.get(k, d))
    v = LineVerdict(wavelength_air_A=float(g("wavelength_air_A", float("nan"))))

    s1, r1 = stage1_depth(g("observed_depth"), th)
    s2, r2 = stage2_saturation(g("rew"), th, d_rew_dA=g("d_rew_dA"))
    s3, r3 = stage3_purity(g("problem_class"), g("excluded_reason", ""))
    s4, r4 = stage4_fit(g("red_chi2"), g("ew_mA"), g("observed_depth"),
                        g("wavelength_air_A"), th,
                        route=g("route"), instrument=g("instrument"),
                        frac_rise_weaker=g("frac_rise_weaker"))
    s5, r5 = stage5_gf(g("gf_tier"))
    v.stages = {"depth": s1, "saturation": s2, "purity": s3, "fit": s4, "gf": s5}

    # A physical failure is INVALID -- an artifact, not the line -- and it needs a NAMED cause.
    # 🔴 "Scattered" is never one of them: that boundary is what keeps a documentation tier from
    # becoming a delete button (design spec §3).
    for verdict, rc in ((s1, r1), (s2, r2), (s3, r3), (s4, r4)):
        if verdict == NO:
            v.tier, v.reason_code = TIER_INVALID, rc
            v.note = "failed a physical stage; documented and excluded from the value, never deleted"
            v.sigma_dex, v.sigma_source = None, "not applicable — line excluded"
            return v

    cited = g("gf_sigma_dex")
    sigma, src = resolve_sigma(cited_sigma_dex=cited,
                               cross_measurement_abundances=g("cross_measurements"),
                               inferred_sigma_dex=g("inferred_sigma_dex"), th=th)
    v.sigma_dex, v.sigma_source = sigma, src

    if str(g("gf_tier") or "").strip() in LAB_TIERS:
        v.tier, v.note = TIER_GOLDEN, "primary laboratory gf; the showcase never moves off it"
        return v

    if anchor is None:
        v.tier = TIER_SCATTERED
        v.reason_code = RC_BAD_GF
        v.note = ("no laboratory anchor available for this cell, so consistency cannot be "
                  "established; documented, excluded from the value")
        return v

    s6, d = stage6_anchor_consistency(g("abundance"), anchor, th)
    v.stages["anchor"] = s6
    v.anchor_distance_dex = None if not np.isfinite(d) else d

    if s6 != YES:
        v.tier = TIER_SCATTERED
        v.reason_code = RC_BAD_GF
        v.note = ("does not behave like the laboratory-anchored distribution; documented and "
                  "excluded from the value, NOT deleted (RYA-711)")
        return v

    # 🔴 THE CONSISTENCY BOUND GATES A *MEASURED* SIGMA, NEVER THE FALLBACK.
    # The fallback is a placeholder for IGNORANCE, not a measurement of width. Applying the
    # bound to it would tier every line we know nothing about as "measurably wide" -- a claim
    # the data does not support -- and would make the admitted tier unreachable for exactly the
    # lines this framework exists to rescue. A line admitted at stage 6 with no measured sigma
    # enters the value carrying the fallback, which is honest and large; that is the design's
    # last-resort route (§4.2), not an exclusion.
    measured = v.sigma_source != "fallback"
    if measured and sigma > th.require("consistency_bound_dex"):
        v.tier = TIER_SCATTERED
        v.reason_code = RC_OUTLIER
        v.note = (f"admitted by the anchor but its own measured sigma ({sigma:.3f} dex) exceeds "
                  f"the consistency bound; documented, excluded from the value, sized to its "
                  f"own scatter")
        return v

    v.tier = TIER_CONSISTENT
    v.note = ("behaves like the laboratory-anchored distribution"
              + ("" if measured else "; carries the fallback sigma — no measured width"))
    return v


# ── pool roll-up ─────────────────────────────────────────────────────────────────────
@dataclass
class PoolGrade:
    """What a cell reports. TWO uncertainties, always."""
    verdicts: list
    anchor: Anchor | None
    n_in_value: int
    #: RMS of the per-line gf sigma over the lines that DO contribute. RMS not median: these
    #: enter the budget in quadrature and the median discards exactly the tail a quadrature sum
    #: is sensitive to (RYA-850's stated reason, reused).
    gf_sigma_dex: float | None
    #: 🔴 does NOT shrink as ungraded lines are admitted.
    zero_point_dex: float | None
    tier_counts: dict
    reason_counts: dict

    def summary(self) -> str:
        z = "n/a" if self.zero_point_dex is None else f"{self.zero_point_dex:.4f}"
        s = "n/a" if self.gf_sigma_dex is None else f"{self.gf_sigma_dex:.4f}"
        return (f"{self.n_in_value} lines in value | gf sigma {s} dex | "
                f"zero-point {z} dex | tiers {self.tier_counts}")


def grade_pool(rows: Iterable, th: Thresholds, anchor: Anchor | None = None) -> PoolGrade:
    """Grade every line in a cell and roll up, WITHOUT collapsing to the worst member.

    🔴 THIS IS THE WHOLE POINT. `gf_rung.decide()` returns rung 1 -- the 0.17 blanket -- the
    instant one line in the pool is ungraded, which is why RYA-855 moved 0 of 36 bars on pools
    that are 7 laboratory lines among 257. Here a mixed pool is a quadrature over per-line sigma
    and the graded lines keep their pedigree instead of surrendering it to their neighbours.
    """
    vs = [grade_line(r, th, anchor) for r in rows]
    keep = [v for v in vs if v.in_value and v.sigma_dex is not None]
    sig = (float(np.sqrt(np.mean(np.asarray([v.sigma_dex for v in keep], float) ** 2)))
           if keep else None)
    tc, rc = {}, {}
    for v in vs:
        tc[v.tier] = tc.get(v.tier, 0) + 1
        if v.reason_code:
            rc[v.reason_code] = rc.get(v.reason_code, 0) + 1
    return PoolGrade(verdicts=vs, anchor=anchor, n_in_value=len(keep), gf_sigma_dex=sig,
                     zero_point_dex=(anchor.zero_point_dex if anchor else None),
                     tier_counts=tc, reason_counts=rc)


def audit_table(pool: PoolGrade) -> list:
    """Per-line rows for the RYA-707/224/851 renderers, including the stage-6 audit.

    ⚠️ Two fields per line is the whole contract with the appendix layer (design spec §6): a
    TIER and a REASON-CODE, in RYA-809's vocabulary. Everything else here is evidence for them.
    """
    out = []
    for v in pool.verdicts:
        out.append({
            "wavelength_air_A": v.wavelength_air_A,
            "tier": v.tier,
            "reason_code": v.reason_code,
            "in_value": v.in_value,
            "gf_sigma_dex": v.sigma_dex,
            "sigma_source": v.sigma_source,
            "anchor_distance_dex": v.anchor_distance_dex,
            "note": v.note,
            **{f"stage_{k}": s for k, s in v.stages.items()},
        })
    return out
