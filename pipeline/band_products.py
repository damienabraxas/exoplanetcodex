"""Per-(element × instrument × band × treatment) products — RYA-712/713.

Ryan, 2026-08-09: *"that median is calculating both engines, while cool to mention, is
not the science product we showcase. We showcase Engine A in one plot/product in the IR,
and the same treatment for B, and LTE."*

WHAT THIS MODULE IS FOR
-----------------------
A **product** is one element, measured on one instrument, in one band, under ONE
treatment, carrying its own value, its own sigma and its own line count. 1D-LTE is a
product. Engine A is a product. Engine B is a product. They are reported side by side
and never merged — RYA-712.

A cross-engine difference (B − A) is a **diagnostic**. It may be mentioned. It is not a
product and this module will not emit one as a headline: `Product` has no field for it,
and `combine()` does not exist.

WHY IT IS PARAMETERISED
-----------------------
Ryan: *"no bad copy paste of code lol, Ba kicked our ass last time."* Adapting the Ba
harness by hand for Al produced 13 defects, two of which (`'element': 'Ba'` and
`SOLAR_ASPLUND2021['Ba']`) would have reached an emitted record on an Al result.

So the element is a PARAMETER on every function here — there is no element symbol
anywhere in the logic, and `assert_single_element()` re-checks the emitted rows against
the element that was asked for. That guard is generic: it does not know or care which
element burned us last time, so it keeps working for the next one.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field, asdict
from pathlib import Path

import numpy as np
import pandas as pd

from pipeline._numcompat import trapezoid as _trapezoid  # numpy>=2 removed np.trapz

# Treatments are separate products. This is the whole vocabulary; there is deliberately
# no 'combined' member, because a combined value is not a product (RYA-712).
#
# ENGINE-B-NLTE is its own member and NOT a variant of ENGINE-B: it is a different
# measurement (TS-native Gerber departures, and a MARCS atmosphere rather than ATLAS9),
# so folding it into ENGINE-B would combine two engines under one label — the exact thing
# this tuple exists to prevent. RYA-785 validated the deck, RYA-798 wired it, and the
# cross-engine spread against ENGINE-A stays an RYA-525 DIAGNOSTIC, never an error bar.
#
# ENGINE-A-3DNLTE (RYA-817) is likewise its own member, not an ENGINE-A variant. It runs
# the same EW route, but the per-line departure comes from the Amarsi+2022 3D-NLTE MLP
# instead of the MPIA/Bergemann 1D-NLTE grid — a different DIMENSIONALITY (3D vs 1D), so
# the two cannot share a label or an aggregate. It is admitted here rather than being
# invented at a call site because RYA-798 shipped ENGINE-B-NLTE without adding it to this
# tuple and the product died at `build_product` after the synthesis had already run.
# 1D-LTE-LABGF (RYA-836) is the same engine and the same route as 1D-LTE, differing
# in ONE input: the oscillator strength comes from a primary laboratory measurement
# instead of the Kurucz semi-empirical value. It is a separate member rather than a
# variant because RYA-712 keys a product on what produced it, and because the two
# must be reportable side by side — the Kurucz pool is the BROAD number over every
# measurable line, the lab-gf pool is the TIGHT number over the few with an
# independent gf. Averaging them would destroy exactly the comparison they exist for.
TREATMENTS = ("1D-LTE", "ENGINE-A", "ENGINE-B", "ENGINE-B-NLTE", "ENGINE-A-3DNLTE",
              "1D-LTE-LABGF")

# Saturation: above this REW the EW->abundance inversion runs along the flat part of the
# curve of growth and is ill-conditioned in BOTH directions. Lines past it are measured
# and reported, never silently dropped -- but they are excluded from the aggregate and
# say why. (RYA-711: quarantined, not culled.)
REW_SATURATION_CEILING = -4.9

_ELEMENT_TOKEN = re.compile(r"\b([A-Z][a-z]?)\s*(?:I{1,3}|IV|VI?|1|2|3)?\b")


@dataclass
class LineMeasurement:
    """One line, one instrument, one treatment. Always carries WHY it was excluded."""
    element: str
    ion: str
    wavelength_air_A: float
    instrument: str
    ew_mA: float
    ew_method: str
    abundance: float | None = None
    rew: float | None = None
    treatment: str = ""
    in_aggregate: bool = True
    excluded_reason: str = ""
    #: RYA-807 — the curated registry's verdict on this line, carried so the per-line
    #: output shows WHY, not merely that it was excluded (RYA-711 quarantine-not-cull,
    #: RYA-429 reported-never-dropped). `problem_action` is 'exclude' | 'flag' | ''.
    #: A FLAGGED line stays in the aggregate: its cause is not established, and removing
    #: it on a hypothesis would be tuning (RYA-161). These default empty, so every
    #: existing caller and artifact is unchanged until the registry actually matches.
    problem_class: str = ""
    problem_status: str = ""
    problem_tickets: str = ""
    problem_action: str = ""
    #: RYA-871 — THE EXCITATION POTENTIAL OF THE TRANSITION THAT WAS MEASURED.
    #:
    #: A wavelength does not identify a line. `gf_rung.resolve_lines` matched a measured
    #: line back to the loaded line list on wavelength ALONE because this column did not
    #: exist, and 16 of 152 VIS Fe I lines did not resolve — 14 with no row inside the
    #: 0.005 A window at all and 2 with two rows inside it. Both halves are the same
    #: missing key: the measured wavelength comes from
    #: `data/audit/line_accounting/per_line.csv`, whose rows are FEATURES rather than
    #: lines (`line_accounting_rya709.features()` groups list rows within 0.05 A and
    #: reports the group MEAN), so a blended feature sits BETWEEN its components by
    #: construction. Widening the window to reach the component then buys a CHOICE rather
    #: than an identification — measured: at 0.020 A with no EP key, 7 of the 136 lines
    #: that already resolved change which row they resolve to. The EP is what makes the
    #: widening legitimate.
    #:
    #: The emitter always had it: the accounting row carries `ep_eV` and
    #: `measure_band_profilefit` copied the wavelength off that same row and dropped this.
    #:
    #: ⚠️ ON A CLUSTERED FEATURE THIS IS THE MINIMUM EP OVER THE CLUSTER, which is what
    #: the accounting table reports. That is a real transition's EP, not an average — and
    #: it is the key that resolved 13 of the 16 with ZERO re-identifications, measured in
    #: `data/results/rya871/`.
    #:
    #: None means "this route does not carry it", never "0 eV" — `gf_rung` falls back to
    #: the narrow wavelength-only rule for such a line rather than widening blind.
    ep_eV: float | None = None
    #: 🔴 RYA-911 — THE CONTINUUM THE HARNESS ACTUALLY PLACED. Read, never reconstructed.
    #:
    #: RYA-897's RCA had to REFIT the window edges to approximate what the profile-fit
    #: harness does, in order to ask whether the HARPS continuum was sitting low. It said
    #: so, and it was right to: an inference about what code does is not a measurement of
    #: what it did (the RYA-869/875 lesson, and the decorative-flag lesson of RYA-904 —
    #: three flags in one file that nothing downstream read).
    #:
    #: The harness has always HAD this number: `_local_renorm` returns `(flux/cont, cont)`
    #: and the driver bound `cont` and threw it away. So the -0.34 dex HARPS Fe II deficit
    #: was diagnosable from the artifact all along, and was not diagnosed, because the one
    #: quantity that decides it never reached a file.
    #:
    #: UNITS ARE THE SOURCE SPECTRUM'S OWN. Kitt Peak residual flux gives 1.0 by
    #: construction; HARPS raw counts give ~1e5. A ratio against `continuum_ref` is the
    #: comparable quantity, and it is deliberately NOT stored — it is two stored readings
    #: divided, and storing a derived number beside its inputs is how RYA-845's
    #: double-count survived.
    continuum_level: float | None = None
    #: HOW that continuum was arrived at, in the harness's own words. `pre_normalised`
    #: and `local-linear-refit` are not two settings of one method: one asserts the
    #: product's continuum IS unity, the other fits a new one and divides.
    continuum_method: str = ""
    #: The SOURCE PRODUCT's own continuum at this wavelength, where it ships one (HARPS
    #: does; the atlases do not). None means "no second opinion exists here", never
    #: "they agree".
    continuum_ref: float | None = None
    #: RYA-880 — THE NLTE CORRECTION, RECORDED WHERE IT IS APPLIED.
    #:
    #: RYA-489 §6.4 requires the correction be SHOWN, "not folded silently", and until now
    #: nothing carried it. RYA-870's reproduce-or-fail guard found the hole by failing: it
    #: re-inverted an EW (an LTE operation) and compared the result against an ENGINE-A
    #: number, and the 0.011 dex gap on 6094.372 WAS the departure — invisible in every
    #: artifact.
    #:
    #: 🔴 RECORDED, NEVER RECONSTRUCTED. Differencing ENGINE-A against 1D-LTE downstream
    #: would assume NLTE is the only thing separating those two rows, and the pools are not
    #: even the same lines (ts-lte deck: ENGINE-A 148 vs 1D-LTE 159), so a line served by
    #: one and not the other would have no delta at all. The deriver knows the number at
    #: the moment it adds it; that is the only honest place to write it down.
    #:
    #: ⚠️ `nlte_delta_dex is None` DOES NOT MEAN ZERO. It means "no additive per-line
    #: correction exists on this route" — which is the truth for ENGINE-B-NLTE, where the
    #: departures enter the radiative transfer itself and the product is a SEPARATE fit,
    #: never a corrected LTE value (RYA-712). An LTE row carries 0.0 and says so.
    nlte_delta_dex: float | None = None
    #: The deck/grid that supplied the departure, or an explicit LTE marker. NEVER blank:
    #: a blank cannot distinguish "LTE" from "nobody recorded it" (RYA-833).
    nlte_source: str = ""
    #: Did this line's abundance come from an EW -> abundance INVERSION?
    #:
    #: RYA-770/342. The REW saturation ceiling below is a property of that inversion —
    #: on the flat part of the curve of growth EW barely responds to abundance, so
    #: inverting it is ill-conditioned. It is NOT a property of the line, and it does
    #: not apply to a flux-space synthesis fit, which measures the profile directly and
    #: exists precisely to recover the strong/blended lines EW saturation kills
    #: (`_run_synthesis_v2_mode` says so in as many words: "the flux-space path does
    #: NOT gate on the EW saturation ceiling — that is an EW-path concept").
    #:
    #: Applying it regardless is how RYA-342 came back. Measured: on the banked
    #: synth-v2 lines the synthesis handler fitted 5 of 6 at red_chi2 1.7-8.8 with a
    #: median A(Fe I) = 7.529 against a banked 7.520 — and this gate quarantined all
    #: but one of them as "saturated", so the control reported "1 line, +0.181 dex".
    #:
    #: Defaults True so every EW-domain caller is unchanged; only the synthesis
    #: handler opts out, and it says why at the call site.
    ew_inversion: bool = True
    #: RYA-847 — WAS THE ABUNDANCE ACTUALLY CONSTRAINED BY THE DATA?
    #:
    #: These are carried because the alternative was measured and it fails silently.
    #: `_fit_synth_flux` has always returned `red_chi2`, and the band-product synthesis
    #: route dropped it on the floor — there was no field here to put it in — which is
    #: how two NIR lines whose chi2 moves 2.2% and 1.4% across EIGHT DEX of iron reached
    #: a published aggregate at 7.833 and 7.979 (RYA-843). A quantity with nowhere to
    #: live is a quantity nobody can gate on.
    #:
    #: None on every EW-route line, and that is correct rather than missing: there is no
    #: chi2 surface behind an EW inversion, so the question does not apply. A consumer
    #: must treat None as "not applicable", never as "unconstrained".
    sigma_A: float | None = None
    frac_rise_weaker: float | None = None
    edge_distance_dex: float | None = None
    red_chi2: float | None = None

    def __post_init__(self) -> None:
        if self.rew is None and self.ew_mA and self.wavelength_air_A:
            self.rew = float(np.log10(self.ew_mA * 1e-3 / self.wavelength_air_A))
        if (self.ew_inversion and self.rew is not None
                and self.rew > REW_SATURATION_CEILING):
            self.in_aggregate = False
            self.excluded_reason = (
                f"REW {self.rew:.3f} above the {REW_SATURATION_CEILING} saturation "
                f"ceiling — the EW->abundance inversion is ill-conditioned here. "
                f"Measured and reported; excluded from the aggregate only.")


@dataclass
class Product:
    """One showcased product. Deliberately has NO field for a cross-engine delta."""
    element: str
    ion: str
    instrument: str
    band: str
    treatment: str
    value: float | None
    sigma: float | None
    n_lines: int
    n_excluded: int
    lines: list[LineMeasurement] = field(default_factory=list)
    provenance: str = ""
    #: RYA-869 — WHICH MEASUREMENT HANDLER PRODUCED THIS NUMBER (`MeasurementHandler.name`,
    #: i.e. "ProfileFitHandler" or "SynthesisHandler").
    #:
    #: The harness residual in the error budget is that handler's own measured optical
    #: systematic, so the budget has to know which handler ran. It used to be inferred
    #: from `treatment`, which is not a function of it: `ENGINE-B-NLTE` is the same flux
    #: fit as `ENGINE-B` and an equality test missed it (four published bars charged the
    #: profile fitter's residual), and in the other direction the near-UV `1D-LTE`
    #: product is a flux fit while the VIS `1D-LTE` product of the same treatment name is
    #: an EW inversion. No mapping from the label exists, so the producing route declares
    #: it here and `pipeline.harness_residual` reads it.
    #:
    #: Defaults empty so a Product can still be constructed directly in a test that is
    #: not about the budget; `build_product` REQUIRES it, and `harness_residual.
    #: for_product` refuses an empty one rather than choosing a default.
    handler: str = ""

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame([asdict(l) for l in self.lines])


def assert_single_element(rows, element: str) -> None:
    """Refuse to emit a record carrying a DIFFERENT element than the one requested.

    Generic on purpose. The Ba->Al incident left `'element': 'Ba'` on an Al result; a
    guard hard-coded to look for 'Ba' would have caught that once and then rotted. This
    compares every emitted element field against the element actually asked for, so it
    catches the next copy-paste too, whichever pair it involves.
    """
    seen = set()
    for r in rows:
        e = r.element if isinstance(r, LineMeasurement) else r.get("element")
        if e:
            seen.add(str(e).strip())
    foreign = seen - {element}
    if foreign:
        raise ValueError(
            f"emitted rows claim element(s) {sorted(foreign)} but this product is for "
            f"{element!r}. This is the copy-paste signature (RYA-701) — a harness "
            f"adapted from another element kept its source's identity. Fix the harness; "
            f"do not filter the rows.")


def assert_no_cross_treatment_mix(lines, treatment: str) -> None:
    """A product is ONE treatment. Mixing is the RYA-712 violation."""
    seen = {l.treatment for l in lines if l.treatment}
    if seen - {treatment}:
        raise ValueError(
            f"product declares treatment {treatment!r} but its lines carry {sorted(seen)}. "
            f"Engines and LTE are separate products and are never combined (RYA-712).")


def local_continuum(w: np.ndarray, f: np.ndarray, centre: float, half_width: float,
                    *, sideband_frac: float = 0.5, pct: float = 95.0) -> float:
    """Continuum from clean side-bands, as the 95th percentile of the flanking regions.

    Side-bands sit OUTSIDE the fitting window so line flux never sets its own continuum.
    A percentile rather than a max, so one hot pixel cannot lift the whole normalisation.
    """
    inner, outer = half_width, half_width * (1.0 + sideband_frac * 2.0)
    d = np.abs(w - centre)
    band = (d > inner) & (d <= outer)
    if band.sum() < 5:
        raise ValueError(f"too few side-band points around {centre:.3f} "
                         f"({int(band.sum())}) — cannot set a local continuum honestly")
    return float(np.percentile(f[band], pct))


# A side-band pseudo-continuum this far below the atlas continuum is itself absorbed.
# Dividing by it removes REAL line flux and biases EW low. Measured on Fe I 6910-9199:
# lines whose side-bands sat at 0.90/0.94 lost 71% and 60% of their EW to re-normalisation.
SIDEBAND_CLEAN_MIN = 0.97   # 0.99 fired on ordinary IR continuum (0.98) and diluted
                            # itself; the real offenders measured 0.902/0.936/0.949


def equivalent_width(w: np.ndarray, f: np.ndarray, centre: float, half_width: float,
                     *, pre_normalised: bool = False, **kw) -> tuple[float, str, str]:
    """EW in mA over +/-half_width. Returns (ew, method, concern).

    `pre_normalised` matters and is not cosmetic. The Kitt Peak atlas ships column 1 as
    RESIDUAL FLUX -- Kurucz already divided by his continuum. Applying our own local
    continuum on top is a SECOND normalisation, and across 4600-9000 A the atlas
    continuum is already excellent (95th pct 0.986-0.997), so that second pass is not
    correcting an error, it is introducing one.

    HARPS is different: it arrives un-normalised and our pipeline sets its continuum.
    Two instruments, two normalisation histories -- which is exactly why a cross-
    instrument abundance difference can be methodological rather than physical.

    So: on pre-normalised data we trust the atlas continuum, and we MEASURE whether the
    side-bands support that. If the side-bands are themselves absorbed, we say so rather
    than silently re-normalising into a pseudo-continuum.
    """
    m = np.abs(w - centre) <= half_width
    if m.sum() < 3:
        raise ValueError(f"too few points inside the window at {centre:.3f}")

    sideband = local_continuum(w, f, centre, half_width, **kw)
    concern = ""
    if pre_normalised:
        cont = 1.0
        how = "atlas continuum trusted (data ship pre-normalised as residual flux)"
        if sideband < SIDEBAND_CLEAN_MIN:
            concern = (f"side-band {pct_label(kw)} is {sideband:.4f}, below "
                       f"{SIDEBAND_CLEAN_MIN} — the local pseudo-continuum is itself "
                       f"absorbed. The atlas continuum is used, but this window sits in "
                       f"crowded spectrum and the EW may include unresolved neighbours.")
    else:
        cont = sideband
        how = (f"local continuum = {pct_label(kw)} of flanking side-bands "
               f"(cont={cont:.5f})")

    depth = 1.0 - (f[m] / cont)
    ew = float(_trapezoid(depth, w[m]) * 1000.0)
    return ew, f"integrated over +/-{half_width:.3f} A, {how}", concern


def carried_ep(row, *, wavelength_A: float, element: str, ion: str) -> float:
    """The excitation potential of the line-accounting row a candidate came from.

    RYA-871 — the EW emitters always HAD this and dropped it, so a measured line could
    only be identified downstream by its wavelength, and 16 of 152 VIS Fe I lines could
    not be identified at all. It is read off the SAME row the wavelength is read off, so
    the two describe one transition by construction.

    🔴 NO SILENT FALLBACK. A candidate with no stateable EP RAISES with the line named.
    Emitting a null EP would put the line back on the wavelength-only rule while looking
    like it carried a key, and once the column exists but is empty a consumer cannot tell
    "this route carries no EP" from "the EP is missing for this line".

    Lives here rather than in either emitter because BOTH drivers
    (`measure_band_profilefit`, `measure_band_ew`) need it and a rule written at two call
    sites drifts between them — the RYA-845/855/869 shape, three tickets deep now.
    """
    ep = getattr(row, "ep_eV", None)
    if ep is None or not np.isfinite(float(ep)):
        raise ValueError(
            f"{element} {ion} {wavelength_A:.4f} A: the line-accounting row carries no "
            f"ep_eV, so the transition cannot be identified downstream on anything but "
            f"its wavelength (RYA-871). Regenerate "
            f"data/audit/line_accounting/per_line.csv "
            f"(scripts/line_accounting_rya709.py) rather than emitting a null EP.")
    return float(ep)


def pct_label(kw) -> str:
    return f"{kw.get('pct', 95.0):.0f}th pct"


def build_product(element: str, ion: str, instrument: str, band: str, treatment: str,
                  lines: list[LineMeasurement], *, handler: str,
                  provenance: str = "") -> Product:
    """Aggregate ONE treatment into ONE product. Never touches another treatment.

    `handler` names the measurement handler that produced these lines and is REQUIRED
    (RYA-869). It is not derivable from `treatment` -- see `Product.handler` -- and the
    route that ran the handler is the only caller that knows it, so it is asked for here
    rather than guessed downstream in the error budget.
    """
    if treatment not in TREATMENTS:
        raise ValueError(f"unknown treatment {treatment!r}; expected one of {TREATMENTS}")
    # Checked against the residual registry rather than against a local list: a handler
    # this product could be built under but whose harness systematic nobody has declared
    # is a product whose error bar cannot be assembled, and finding that out here names
    # the route instead of failing three calls later inside `error_budget.build`.
    from pipeline.harness_residual import for_handler
    for_handler(handler)
    assert_single_element(lines, element)
    assert_no_cross_treatment_mix(lines, treatment)

    used = [l for l in lines if l.in_aggregate and l.abundance is not None]
    excluded = [l for l in lines if not l.in_aggregate]
    for l in excluded:
        if not l.excluded_reason.strip():
            raise ValueError(
                f"{element} {l.wavelength_air_A:.3f} is excluded with no reason recorded. "
                f"An unexplained exclusion is indistinguishable from a silent drop "
                f"(RYA-429).")

    vals = np.array([l.abundance for l in used], dtype=float)
    value = float(np.median(vals)) if len(vals) else None
    # Sigma is the scatter of the lines that actually entered, not a fitted error bar.
    sigma = float(np.std(vals, ddof=1)) if len(vals) > 1 else None
    return Product(element=element, ion=ion, instrument=instrument, band=band,
                   treatment=treatment, value=value, sigma=sigma,
                   n_lines=len(used), n_excluded=len(excluded),
                   lines=lines, provenance=provenance, handler=handler)


def products_frame(products: list[Product]) -> pd.DataFrame:
    """Side-by-side table of products. One row per product — NOT a merge.

    Any comparison a reader makes between rows is theirs to make and is a diagnostic.
    This function will not compute it for them.
    """
    return pd.DataFrame([
        dict(element=p.element, ion=p.ion, instrument=p.instrument, band=p.band,
             treatment=p.treatment, value=p.value, sigma=p.sigma,
             n_lines=p.n_lines, n_excluded=p.n_excluded, handler=p.handler,
             provenance=p.provenance)
        for p in products])
