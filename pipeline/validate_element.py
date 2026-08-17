#!/usr/bin/env python3
"""
pipeline/validate_element.py — RYA-813 (child of RYA-812)
=========================================================
THE PER-ELEMENT VALIDATE-AND-REPORT STAGE. Gold is OUTPUT ONLY: this stage reads
the OUTSIDE WORLD (litscan literature + Asplund) and the element's own per-band x
per-engine products, and NEVER reads the frozen gold reference.

WHY THAT ONE RULE IS THE WHOLE POINT
------------------------------------
Today `scripts/phase_c_verdict_rya371.py` reads the frozen gold to classify — both
the VALUE (`_abs_value(mrow)`) and the SCALE LABEL (`apply_reported_scale_correction`
on the gold row). Gold gates the thing that produces gold:

    gold vN --> phase_c verdict --> gold candidate --> freeze --> gold vN+1

A stale label in ONE cell therefore vetoes ALL 28 elements, which is exactly what
blocked Fe. Validating against literature instead means the circle cannot form —
not "is discouraged from forming". `assert_no_gold_read()` makes that structural
rather than aspirational, and `tests/test_validate_element_rya813.py` asserts it.

THE TWO TIERS (RYA-812 principles 1, 2, 7)
------------------------------------------
    VIS         VALIDATE  -- must land inside the element's litscan band, else
                             pass-with-exception (documented) or fail-with-reason.
    near-UV     REPORT    -- honest bars; compare to literature ONLY where a
    red-optical            reference exists; NEVER gated, NEVER required to hit.
    NIR

A frontier band is never failed for missing a literature number. That is not
leniency — it is refusing to punish being first. The discipline is being honest
about WHICH tier a band is in, and the appendix carries the burden of proof.

An element with no usable literature anywhere is UN-ANCHORABLE: report-only across
every band, flagged loudly. An honest gap, not a failure.

NOTHING HERE IS SILENT
----------------------
Every exception, every un-anchorable, every scale mismatch and every missing
product is a named field on the returned verdict object and appears in the
appendix. A skipped check that leaves no trace is the RYA-786 defect class.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path as _Path
from typing import Any, Optional

from pipeline.litscan import LiteratureRange, literature_range

# ── the two tiers ────────────────────────────────────────────────────────────
VALIDATING_BANDS = ("VIS",)
REPORTING_BANDS = ("near-UV", "red-optical", "NIR")
ALL_BANDS = VALIDATING_BANDS + REPORTING_BANDS

# verdict vocabulary — deliberately distinct strings, never booleans
PASS = "pass"
PASS_WITH_EXCEPTION = "pass-with-exception"
FAIL = "fail-with-reason"
REPORT = "report"
DEVIATE = "deviate"
UN_ANCHORABLE = "un-anchorable"
NO_PRODUCT = "no-product"

#: The curated exclusion registry (RYA-807). Named with a leading underscore and read
#: lazily inside the appendix so this module's global namespace stays as small as
#: `assert_no_gold_read` expects -- it inspects the bound globals, and a module that
#: accumulates imports is one where that check gets harder to trust. The registry is not
#: gold: it records which LINES are excluded, never an abundance.
_PROBLEM_CHILDREN = (_Path(__file__).resolve().parents[1]
                     / "data" / "registry" / "problem_children.csv")

#: Modules/paths this stage must never touch. Gold is output-only (principle 5).
FORBIDDEN_GOLD_SURFACES = (
    "read_solar_reference",
    "differential_denominator",
    "solar_abundances_v",
    "data/reference/solar/CURRENT",
)


class GoldReadInValidationError(RuntimeError):
    """Raised if the validation path touches the frozen gold reference."""


def assert_no_gold_read() -> None:
    """
    Structural guarantee that gold is not reachable from this stage.

    Checks that `pipeline.data_namespace` — the ONLY module that opens a frozen
    gold reference — has not been imported into this module's namespace, and that
    no gold accessor is bound here. This is a cheap, honest check: it cannot stop
    a determined caller from importing gold elsewhere, and it does not pretend to.
    What it does stop is this file quietly growing a gold read later, which is how
    the original loop formed in the first place.
    """
    bound = set(globals())
    for name in FORBIDDEN_GOLD_SURFACES:
        if name in bound:
            raise GoldReadInValidationError(
                f"validate_element bound {name!r} — gold is OUTPUT ONLY (RYA-812 "
                f"principle 5). Validation reads literature, never the frozen "
                f"reference it will eventually be assembled into.")
    mod = sys.modules.get(__name__)
    src_globals = getattr(mod, "__dict__", {})
    if "data_namespace" in src_globals or "ns" in src_globals:
        raise GoldReadInValidationError(
            "validate_element imported pipeline.data_namespace. That module reads "
            "the frozen gold reference; importing it here re-opens the "
            "gold -> verdict -> gold loop RYA-812 exists to kill.")


# ── measurement inputs ───────────────────────────────────────────────────────
@dataclass(frozen=True)
class BandProduct:
    """
    One (band x engine) measurement product. RYA-712: products are SEPARATE data
    products and are NEVER combined — this stage reads them side by side and
    reports each on its own terms. There is deliberately no method that averages
    two engines into one number.
    """
    band: str
    engine: str
    value: float
    sigma_statistical: Optional[float] = None
    sigma_systematic: Optional[float] = None
    n_lines: Optional[int] = None
    scale: Optional[str] = None
    note: str = ""

    @property
    def bars(self) -> str:
        """Honest bars: statistical and systematic stay APART (pipeline.error_budget)."""
        st = "?" if self.sigma_statistical is None else f"{self.sigma_statistical:.3f}"
        sy = "?" if self.sigma_systematic is None else f"{self.sigma_systematic:.3f}"
        return f"stat {st} / sys {sy}"


@dataclass
class BandVerdict:
    band: str
    tier: str                      # 'validate' | 'report'
    engine: str = ""
    value: Optional[float] = None
    verdict: str = NO_PRODUCT
    reason: str = ""
    exception_reason: str = ""
    literature_central: Optional[float] = None
    deviate: bool = False
    best_external: str = ""
    literature_range: Optional[tuple[float, float]] = None
    offset_vs_literature: Optional[float] = None
    excess_beyond_band: Optional[float] = None
    bars: str = ""
    scale: Optional[str] = None
    scale_mismatch: str = ""
    references_used: tuple[str, ...] = ()

    @property
    def is_exception(self) -> bool:
        return self.verdict == PASS_WITH_EXCEPTION

    @property
    def is_loud(self) -> bool:
        """Anything a reader must not scroll past."""
        return self.verdict in (PASS_WITH_EXCEPTION, FAIL, UN_ANCHORABLE) or bool(
            self.scale_mismatch) or self.deviate


@dataclass
class ElementVerdict:
    element: str
    asplund2021: Optional[float]
    un_anchorable: bool
    bands: list[BandVerdict] = field(default_factory=list)
    litscan_present: bool = False
    litscan_derivation: str = ""
    litscan_is_statistical: bool = False
    verification_owed: tuple[str, ...] = ()
    references_used: tuple[str, ...] = ()

    @property
    def exceptions(self) -> list[BandVerdict]:
        return [b for b in self.bands if b.is_exception]

    @property
    def failures(self) -> list[BandVerdict]:
        return [b for b in self.bands if b.verdict == FAIL]

    @property
    def loud(self) -> list[BandVerdict]:
        return [b for b in self.bands if b.is_loud]

    @property
    def validating(self) -> list[BandVerdict]:
        return [b for b in self.bands if b.tier == "validate"]

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["exception_count"] = len(self.exceptions)
        d["failure_count"] = len(self.failures)
        return d


# ── the stage ────────────────────────────────────────────────────────────────
def validate_element(
    element: str,
    products: list[BandProduct],
    *,
    asplund: Optional[float] = None,
    exceptions: Optional[dict[str, str]] = None,
) -> ElementVerdict:
    """
    Validate (VIS) and report (frontier) one element against the LITERATURE.

    `products`   the element's per-band x per-engine products (RYA-712, unmerged).
    `asplund`    the Asplund-2021 reference, passed in rather than imported so this
                 stage has no opinion about where constants live.
    `exceptions` band -> documented reason. A VIS miss WITH a reason is
                 `pass-with-exception`; a VIS miss WITHOUT one is
                 `fail-with-reason`. An exception is never silent: it is recorded on
                 the verdict and printed in the appendix.

    NEVER reads gold. See `assert_no_gold_read`.
    """
    assert_no_gold_read()
    exceptions = exceptions or {}

    lit = literature_range(element)
    un_anchorable = lit is None

    ev = ElementVerdict(
        element=element,
        asplund2021=asplund,
        un_anchorable=un_anchorable,
        litscan_present=lit is not None,
        litscan_derivation=(lit.derivation if lit else ""),
        litscan_is_statistical=(lit.is_statistical_spread if lit else False),
        verification_owed=(lit.verification_owed if lit else ()),
        references_used=(lit.citations if lit else ()),
    )

    by_band: dict[str, list[BandProduct]] = {}
    for p in products:
        by_band.setdefault(p.band, []).append(p)

    for band in ALL_BANDS:
        tier = "validate" if band in VALIDATING_BANDS else "report"
        prods = by_band.get(band, [])
        if not prods:
            continue                      # a band with no product is simply absent
        for p in prods:
            ev.bands.append(_judge(band, tier, p, lit, exceptions, un_anchorable))

    # a product in a band we do not recognise is a loud finding, never dropped
    for band, prods in by_band.items():
        if band in ALL_BANDS:
            continue
        for p in prods:
            ev.bands.append(BandVerdict(
                band=band, tier="report", engine=p.engine, value=p.value,
                verdict=REPORT, bars=p.bars,
                reason=(f"band {band!r} is not in the declared vocabulary "
                        f"{ALL_BANDS} — reported, never validated, and surfaced so "
                        f"the vocabulary gap is visible rather than silent.")))
    return ev


def _judge(band: str, tier: str, p: BandProduct, lit: Optional[LiteratureRange],
           exceptions: dict[str, str], un_anchorable: bool) -> BandVerdict:
    bv = BandVerdict(band=band, tier=tier, engine=p.engine, value=p.value,
                     bars=p.bars, scale=p.scale,
                     references_used=(lit.citations if lit else ()))

    # scale mismatch is recorded on EVERY band, and it never silently invalidates a
    # comparison -- it says so instead (a 1D value vs a 3D band is a scale error
    # dressed as a disagreement).
    if lit is not None and p.scale and lit.scale and p.scale != lit.scale:
        bv.scale_mismatch = (
            f"product is on '{p.scale}' but the literature range is quoted on "
            f"'{lit.scale}'. The comparison below is NOT scale-consistent and must "
            f"be corrected onto '{lit.scale}' before it means anything "
            f"(pipeline.solar_scale_provenance owns that transform).")

    if un_anchorable:
        bv.verdict = UN_ANCHORABLE
        bv.reason = (
            f"{p.engine}: no usable literature range for this element, so there is "
            f"nothing to validate against in ANY band. Reported at {p.value:.3f} "
            f"({p.bars}). This is an honest gap, NOT a failure (RYA-812 principle 2).")
        return bv

    bv.literature_central = lit.central
    bv.literature_range = (lit.min, lit.max)
    bv.best_external = lit.best_external
    bv.deviate = lit.is_deviate(p.value)
    bv.offset_vs_literature = round(lit.offset(p.value), 4)
    bv.excess_beyond_band = round(lit.excess(p.value), 4)

    if tier == "report":
        inside = lit.contains(p.value)
        bv.verdict = REPORT
        bv.reason = (
            f"{p.engine}: frontier band — REPORTED, not gated (RYA-777). "
            f"{p.value:.3f} ({p.bars}); literature {lit.central:.3f} "
            f"[{lit.min:.3f}, {lit.max:.3f}] — "
            + ("agrees with the reference range."
               if inside else
               (f"sits {bv.excess_beyond_band:+.3f} dex outside it"
                + (" and beyond 2-sigma (RYA-714 DEVIATE)" if bv.deviate else "")
                + ", NOTED as a gap and NOT a failure."))
            + (f" {p.note}" if p.note else ""))
        return bv

    # ── VIS: the one band that validates ─────────────────────────────────────
    if lit.contains(p.value):
        bv.verdict = PASS
        bv.reason = (
            f"{p.engine}: {p.value:.3f} is inside the literature range "
            f"[{lit.min:.3f}, {lit.max:.3f}] (central {lit.central:.3f}, offset "
            f"{bv.offset_vs_literature:+.3f}).")
        return bv

    reason = exceptions.get(band) or exceptions.get(f"{band}:{p.engine}")
    if reason:
        bv.verdict = PASS_WITH_EXCEPTION
        bv.exception_reason = reason
        bv.reason = (
            f"{p.engine}: {p.value:.3f} is {bv.excess_beyond_band:+.3f} dex OUTSIDE "
            f"the literature range [{lit.min:.3f}, {lit.max:.3f}]. Passed by "
            f"DOCUMENTED EXCEPTION — see exception_reason; recorded in the appendix.")
        return bv

    bv.verdict = FAIL
    bv.reason = (
        f"{p.engine}: {p.value:.3f} is {bv.excess_beyond_band:+.3f} dex outside the "
        f"literature range [{lit.min:.3f}, {lit.max:.3f}]"
        + (f" and beyond 2*sigma_ext — RYA-714 §4 codes this DEVIATE, a stronger "
           f"claim than a near miss" if bv.deviate else "")
        + f" (best external: {lit.best_external or 'unnamed'}). NO exception was "
        f"documented. A VIS miss is a fail-with-reason, never a silent miss "
        f"(RYA-812 principle 1).")
    return bv


# ── appendix (the burden of proof) ───────────────────────────────────────────
def excluded_lines_section(element: str) -> list[str]:
    """Name every line this element currently excludes, with its stated cause.

    RYA-844 puts the burden of proof on the appendix, and RYA-847 items 5 and 7 require
    each excluded line to be NAMED here with its physical cause. Until now the appendix
    reported band verdicts and never said which lines had been removed to reach them — so
    a reader could see A(Fe) and its bars without seeing that a line had been dropped, and
    an exclusion and a tuning look identical in the output. Only the stated reason
    separates them, which is exactly why the reason has to be printed where the number is.

    SOURCED FROM THE REGISTRY, not from a list kept here. `data/registry/
    problem_children.csv` is the single declaration (RYA-807); a second copy beside the
    report is how the two come to disagree (RYA-845).

    The discriminator is `status`, NOT `required_treatment` (RYA-807): only `exclude` +
    `active` is actually excluded. A row that is `owed` is KEPT AND FLAGGED, because its
    cause is not established and removing it on a hypothesis would be tuning (RYA-161).
    Both are listed, under headings that say which is which -- a flagged line the reader
    cannot see is the same defect one step removed.
    """
    reg = _PROBLEM_CHILDREN
    if not reg.exists():
        return ["## Excluded and flagged lines", "",
                f"> ⚠️ `{reg.name}` is absent, so this section can make NO claim about "
                f"exclusions. That is a missing input, not an empty exclusion set.", ""]
    import pandas as _pd
    d = _pd.read_csv(reg, dtype=str).fillna("")
    sp = d.species.astype(str).str.strip()
    mine = d[(sp == element) | sp.str.startswith(f"{element} ")]
    if mine.empty:
        return ["## Excluded and flagged lines", "",
                f"None registered for {element}.", ""]

    _NO_REASON = ("⚠️ NO REASON RECORDED — an exclusion without a stated cause cannot "
                  "be distinguished from tuning (RYA-844).")

    def _rows(sub, kind):
        out = []
        for _, r in sub.iterrows():
            out += [f"- **{r['species']} {r['lambda_or_scope']}** — "
                    f"`{r['problem_class']}` / {kind}, severity {r['severity'] or 'n/a'} "
                    f"[{r['governing_tickets']}]",
                    f"  - {r['notes'] or _NO_REASON}"]
        return out

    excluded = mine[(mine.required_treatment.str.strip() == "exclude")
                    & (mine.status.str.strip() == "active")]
    flagged = mine[mine.status.str.strip() == "owed"]
    L = ["## Excluded and flagged lines", "",
         "Sourced from `data/registry/problem_children.csv`. The discriminator is "
         "`status`, not `required_treatment` (RYA-807).", ""]
    if len(excluded):
        L += [f"### Excluded from the aggregate ({len(excluded)})", ""] + \
             _rows(excluded, "EXCLUDED") + [""]
    else:
        L += ["### Excluded from the aggregate (0)", "", "None.", ""]
    if len(flagged):
        L += [f"### Kept and flagged — cause not established ({len(flagged)})", "",
              "Retained in the aggregate deliberately: removing a line on an "
              "undiagnosed cause is tuning (RYA-161/844).", ""] + \
             _rows(flagged, "KEPT+FLAGGED") + [""]
    return L


def appendix(ev: ElementVerdict) -> str:
    """Per element, per band: measured vs literature, the call, and every reason."""
    L = [f"# Validation appendix — {ev.element}",
         "",
         f"Asplund 2021 reference: "
         f"{'(none supplied)' if ev.asplund2021 is None else f'{ev.asplund2021:.3f}'}"]
    if ev.un_anchorable:
        L += ["", "**UN-ANCHORABLE** — no usable literature range. Report-only across "
                  "all bands. This is an honest gap, not a failure."]
    else:
        L += ["", f"Litscan band derivation: `{ev.litscan_derivation}` "
                  f"(statistical spread across determinations: "
                  f"{ev.litscan_is_statistical})"]
        # The caveat must match the DERIVATION. An earlier version printed
        # "ratified POLICY window" for every non-statistical band, which became
        # wrong the moment the band came from the literature's own sigma -- the
        # note would have understated a perfectly good band.
        if ev.litscan_derivation == "best-external-plus-sigma":
            L += ["", "> The band is the NAMED best external determination ± its "
                      "published σ, which is the litscan dossier's own ratified rule "
                      "(PASS if |Δ| ≤ σ_ext). It is literature-derived — but it is NOT "
                      "the envelope of every determination, which is wider. Both are "
                      "recorded in the litscan so the choice is visible."]
        elif not ev.litscan_is_statistical:
            L += ["", "> ⚠️ The pass band is NOT a literature spread and NOT a published "
                      "σ — check the litscan's `derivation` before relying on it. A band "
                      "of unclear provenance must not silently act as a pass criterion."]
    L += ["", "## Bands", ""]
    for b in ev.bands:
        L += [f"### {b.band} — {b.tier.upper()} — `{b.verdict}`",
              "",
              f"- engine: `{b.engine}`",
              f"- measured: {'n/a' if b.value is None else f'{b.value:.3f}'}  ({b.bars})"]
        if b.literature_central is not None:
            L += [f"- literature: {b.literature_central:.3f} "
                  f"[{b.literature_range[0]:.3f}, {b.literature_range[1]:.3f}]",
                  f"- offset: {b.offset_vs_literature:+.3f} dex"
                  + (f"  (excess beyond band {b.excess_beyond_band:+.3f})"
                     if b.excess_beyond_band else "")]
        L += [f"- call: {b.reason}"]
        if b.exception_reason:
            L += ["", f"  **DOCUMENTED EXCEPTION.** {b.exception_reason}"]
        if b.scale_mismatch:
            L += ["", f"  **⚠️ SCALE MISMATCH.** {b.scale_mismatch}"]
        L += [""]
    L += excluded_lines_section(ev.element)
    if ev.references_used:
        L += ["## References used", ""] + [f"- {r}" for r in ev.references_used] + [""]
    if ev.verification_owed:
        L += ["## Verification owed (litscan)", ""] + [
            f"- {v}" for v in ev.verification_owed] + [""]
    L += ["## Provenance", "",
          "- Gold reference: **NOT READ.** This stage validates against literature "
          "only (RYA-812 principle 5); `assert_no_gold_read()` enforces it and "
          "`tests/test_validate_element_rya813.py` asserts it.",
          "- Products are read per (band × engine) and are NEVER merged (RYA-712)."]
    return "\n".join(L)
