"""RYA-1079: can this line be MEASURED here? -- the per-line telluric observability gate.

THE POLICY, ratified by Ryan 2026-08-27
---------------------------------------
**Correction is the standard, not exclusion.** Now that molecfit is mature, skipping a
line because it sits in a telluric zone throws away recoverable data. The question is not
"is this line inside an enumerated band" (RYA-460/786) but "is this line's region
telluric-corrected?" -- yes, use it; no, correct it, then use it. Exclusion is the rare
exception.

**The exception is PHYSICS, not policy.** Where O2 absorbs to ~0 there are no stellar
photons left, and dividing signal by a number that went to zero recovers nothing. That
line is genuinely lost. Measured here, on our own atlases: the O2 B core at 6870 A drives
the Kitt Peak/HARPS ratio to a transmission of -0.002 / 0.006. Nothing is recoverable
from that.

So there are THREE per-line states, and they are DECIDED FROM THE MEASURED TELLURIC DEPTH
over the line's own window -- never asserted from `instrument_catalog.telluric_basis`:

    CLEAN        no telluric detectable above the spectrum's own noise -> measure as-is
    RECOVERABLE  a real telluric, but correction leaves usable S/N     -> CORRECT, then
                 measure. 🔴 THIS IS THE CLASS BEING RESCUED: today it is silently skipped.
    SATURATED    the core takes transmission below what correction can survive -> exclude,
                 to problem_children with reason `telluric_saturated_core` and the
                 measured depth attached.

TWO AXES, KEPT ORTHOGONAL, AND OBSERVABILITY RUNS FIRST
-------------------------------------------------------
GRADE is about the *gf*: was the transition probability lab-measured. TELLURIC STATE is
about the *observation*: can an EW be extracted here at all. A line can carry a pristine
Ruffoni lab gf and still be unmeasurable because the sky ate the photons.

So `partition_pool` runs BEFORE tier assignment, and a SATURATED line never reaches the
grader. It is recorded **unmeasured, with a reason** -- NOT ungraded. Recording it ungraded
would say "we looked and the gf was weak" when the truth is "we could not look": the same
label-vs-reality defect as RYA-833's un-normalised/unmeasured collapse and RYA-1072's
ambiguity-reads-as-no.

WHERE THE DEPTH COMES FROM -- MEASURED, NOT MODELLED
-----------------------------------------------------
Transmission over a window is `uncorrected_flux / corrected_flux` for a SIBLING PAIR: two
holdings of the same Sun on the same instrument, one telluric-corrected and one not. That
ratio is the atmosphere and nothing else -- the solar spectrum divides out, because it is
the same Sun through the same instrument.

The pairs are DERIVED from `holdings_manifest_registry.csv` (same instrument_id, one
`telluric_applied=applied` and one `not-applied`), not listed here. Measured controls, on
the pairs that read on this workstation:

    window            KP raw/corrected        HARPS raw/corrected
    6700 A (clean)    T = 1.0000 exactly      T = 1.0002
    6870 A (O2 B)     T_min = -0.0022         T_min = 0.0060

The clean window returning EXACTLY 1.0 is the negative control that matters: it says the
ratio is measuring the atmosphere and not a normalisation difference between the two
products.

THE THRESHOLD IS DERIVED, NOT BORROWED (RYA-161)
-------------------------------------------------
🔴 NO MAGIC NUMBER. Both edges fall out of one photon-noise argument and one constant the
repo already declares.

Dividing an observed spectrum by transmission T recovers the signal and amplifies the
noise. Photon-noise limited, observed counts are N0*T, so the observed S/N is S0*sqrt(T)
with S0 the continuum S/N; correction scales signal and noise together and the post-
correction S/N is the same S0*sqrt(T). Requiring it to clear the repo's own science floor:

    S0 * sqrt(T) >= PIPELINE['snr_min_science']   =>   T >= (SNR_min / S0)^2

so `saturated_min_depth = 1 - (SNR_min/S0)^2`. Below that transmission, correction is
ill-conditioned: the number exists, the S/N behind it does not.

The CLEAN edge is the other end of the same statistic. A telluric shallower than the
spectrum's own per-pixel noise is not detectable in that spectrum and cannot bias an EW
by more than the noise already does, so `clean_max_depth = 1 / S0`. No multiplier is
applied -- a "3 sigma" would be a borrowed convention, and this needs none.

BOTH EDGES ARE PER-HOLDING, because S0 is measured per holding. That is not a detail: on
our own two arms the derived saturation edge differs by a factor of 30 (Kitt Peak
S0 ~ 2510 -> SATURATED below T = 0.006; HARPS S0 ~ 458 -> SATURATED below T = 0.19). A
single global cut would be wrong for both.

WHAT THIS RETIRES
-----------------
`telluric_basis = line_selection` stops being a DECISION INPUT. It stays in the catalog as
description; the per-line verdict decides. That is where RYA-928's red-optical lines get
CORRECTED rather than kept-uncorrected-or-excluded, and it supersedes the RYA-460/786
enumerated-band exclusion method as the matured-molecfit successor.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
HOLDINGS = ROOT / "data" / "catalog" / "holdings_manifest_registry.csv"

#: The three per-line observability states. Ratified RYA-1079.
CLEAN = "CLEAN"
RECOVERABLE = "RECOVERABLE"
SATURATED = "SATURATED"

#: What a caller must DO with the line. Distinct from the state, because the same state
#: means different work depending on whether this holding is already corrected.
MEASURE = "MEASURE"
NEEDS_CORRECTION = "NEEDS_CORRECTION"
EXCLUDE = "EXCLUDE"

#: The problem_children reason for the one genuinely lost class. Named once, here.
SATURATED_CORE_REASON = "telluric_saturated_core"


class TransmissionUnavailable(RuntimeError):
    """No sibling pair -- the telluric depth here has not been measured.

    🔴 NEVER DEFAULTED. An unmeasured depth is not CLEAN. Reading it as clean would let a
    line inside an O2 core through as if nothing were there, which is the exact failure
    this module exists to make impossible; reading it as SATURATED would throw away the
    recoverable class this module exists to rescue. It is refused, and the caller records
    that the observation has not been characterised.
    """


def snr_floor() -> float:
    """The science S/N floor -- READ from the repo's constants, never typed here."""
    from config.constants import PIPELINE
    return float(PIPELINE["snr_min_science"])


def continuum_snr(flux) -> float:
    """Per-pixel continuum S/N, from the SUCCESSIVE DIFFERENCES of the flux.

    Differencing removes the solar spectrum -- real line structure is smooth on the pixel
    scale while noise is not -- so this measures the noise without needing a continuum fit
    or a line mask. Robust (MAD) so a cosmic ray or a saturated core cannot set the scale.
    """
    f = np.asarray(flux, float)
    f = f[np.isfinite(f)]
    if f.size < 3:
        raise TransmissionUnavailable(
            f"fewer than 3 finite pixels ({f.size}); the noise cannot be measured and "
            f"both observability edges are derived from it.")
    # 1.4826 converts MAD to sigma for a Gaussian and sqrt(2) undoes the variance
    # doubling of a first difference -- both are DEFINITIONS of this estimator, not
    # thresholds, and neither is tunable.
    #
    # ⚠️ KNOWN BIAS, STATED: where real structure varies on the pixel scale the difference
    # picks it up as noise, so S0 is UNDER-estimated. That lowers `saturated_min_depth`
    # and calls MORE lines saturated -- the conservative direction, costing coverage
    # rather than admitting an unmeasurable line. Measured on a synthetic control it
    # recovers S/N 200 to within 4% and S/N 1000 to within 25% once structure approaches
    # the noise scale.
    d = np.diff(f)
    sigma = 1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2.0)
    if not np.isfinite(sigma) or sigma <= 0:
        raise TransmissionUnavailable(
            "the flux has zero measurable pixel-to-pixel scatter, so its S/N is not "
            "determined; refusing rather than assuming a noiseless spectrum.")
    return float(1.0 / sigma)


@dataclass(frozen=True)
class Thresholds:
    """The two edges, DERIVED from one measured S/N and one declared constant."""
    snr_continuum: float
    snr_min_science: float

    @property
    def clean_max_depth(self) -> float:
        """A telluric below the per-pixel noise is not detectable in this spectrum."""
        return 1.0 / self.snr_continuum

    @property
    def saturated_min_depth(self) -> float:
        """Below this transmission, post-correction S/N falls under the science floor."""
        return 1.0 - (self.snr_min_science / self.snr_continuum) ** 2

    def cite(self) -> str:
        return (f"S0={self.snr_continuum:.1f} measured; SNR_min="
                f"{self.snr_min_science:.0f} (PIPELINE['snr_min_science']); "
                f"CLEAN at depth<={self.clean_max_depth:.5f} (=1/S0); SATURATED at "
                f"depth>={self.saturated_min_depth:.4f} (=1-(SNR_min/S0)^2)")


def thresholds(snr_continuum: float) -> Thresholds:
    return Thresholds(snr_continuum=float(snr_continuum), snr_min_science=snr_floor())


@dataclass
class Observability:
    """The per-line verdict, with the measurement it was read off."""
    holding_id: str
    wavelength_A: float
    verdict: str
    disposition: str
    depth: float
    transmission_min: float
    thresholds: Thresholds
    corrected: bool
    reason: str = ""
    evidence: str = ""

    @property
    def measurable(self) -> bool:
        """Can an EW be extracted here AT ALL? Not 'is it ready' -- see `disposition`.

        A RECOVERABLE line on an uncorrected holding is measurable; it is simply not
        measurable YET, and the difference between those two is the whole ticket.
        """
        return self.verdict != SATURATED


# ── the sibling pair, DERIVED from the registry ──────────────────────────────

def _registry() -> list[dict]:
    return list(csv.DictReader(HOLDINGS.open(encoding="utf-8")))


def corrected_siblings(holding_id: str) -> list[str]:
    """EVERY telluric-corrected holding registered for the same instrument, in order.

    Derived from the registry, never tabulated here: a hand-written pair table would be a
    second home for a fact the registry already carries, and it would go stale the first
    time a corrected product is registered (RYA-350/353/967).

    A LIST, not one id, because "registered" and "readable on this machine" are different
    facts. `kpno_solar_atlas` registers two corrected siblings and only one of them is
    staged on the workstation; returning the first and letting the read explode would
    report a MISSING FILE as if it were a telluric finding -- an environment fault wearing
    a science verdict's clothes (RYA-1064 friction #3, same shape).
    """
    rows = _registry()
    by_id = {r["holding_id"]: r for r in rows}
    if holding_id not in by_id:
        raise KeyError(f"holding {holding_id!r} is not registered; its telluric depth "
                       f"cannot be measured (RYA-806).")
    inst = by_id[holding_id]["instrument_id"]
    return [r["holding_id"] for r in rows
            if r["instrument_id"] == inst and r["holding_id"] != holding_id
            and r.get("telluric_applied", "").strip() == "applied"]


def is_corrected(holding_id: str) -> bool:
    from pipeline.telluric_policy import applied_state
    return applied_state(holding_id) == "applied"


def measure_transmission(harness, instrument: str, holding_id: str,
                         centre_A: float, pad_A: float):
    """(transmission array, evidence) over this window. MEASURED from a sibling pair.

    `harness` is the measurement harness module (`scripts/measure_band_ew.py`), passed in
    rather than imported so this module never inherits its import-time atlas resolution
    (RYA-1064 friction #3) -- a policy module must be importable with no data staged.
    """
    siblings = corrected_siblings(holding_id)
    if not siblings:
        raise TransmissionUnavailable(
            f"{holding_id}: no telluric-corrected sibling is registered for "
            f"{instrument}, so the telluric depth over this window has never been "
            f"measured. It is NOT assumed clean and NOT assumed saturated -- correct this "
            f"holding (RYA-424) or register a corrected sibling, then re-ask.")
    if is_corrected(holding_id):
        raise TransmissionUnavailable(
            f"{holding_id} is itself the corrected product; the depth must be measured "
            f"against its UNCORRECTED sibling, which is the holding that carries the "
            f"atmosphere. Ask about the raw holding.")
    raw = holding_id

    a = harness.load_window_ex(instrument, centre_A, pad_A, holding=raw,
                               allow_uncorrected=True)
    b, sibling, tried = None, None, []
    for cand in siblings:
        try:
            b = harness.load_window_ex(instrument, centre_A, pad_A, holding=cand,
                                       allow_uncorrected=True)
            sibling = cand
            break
        except Exception as exc:                       # not staged here, or no coverage
            tried.append(f"{cand}: {type(exc).__name__}")
    if b is None:
        raise TransmissionUnavailable(
            f"{holding_id} @ {centre_A:.3f} A: {len(siblings)} corrected sibling(s) are "
            f"registered and NONE could be read here -- {'; '.join(tried)}. That is a "
            f"staging fact about this machine, not a telluric verdict, so the depth stays "
            f"unmeasured rather than being defaulted either way.")

    grid = np.linspace(centre_A - pad_A, centre_A + pad_A, 800)
    fa = np.interp(grid, a.wave, a.flux)
    fb = np.interp(grid, b.wave, b.flux)
    # A corrected flux at or below zero cannot divide -- that is a quarantined core, not a
    # transmission. Those pixels are dropped, not floored: a floored ratio is a number we
    # made up (RYA-844's rule, on pixels).
    #
    # ⚠️ THE TEST IS `> 0`, NOT A FLOOR. An earlier draft used `fb > 0.05`, which is a
    # threshold with no derivation behind it -- exactly what RYA-161 forbids, smuggled in
    # as a housekeeping guard. It is not needed: the statistic taken from T is its
    # MINIMUM, and a vanishing denominator drives T UP, away from the saturated edge. A
    # near-zero corrected flux can therefore never manufacture a SATURATED verdict, so
    # dropping only the undividable pixels is both sufficient and derivation-free.
    ok = np.isfinite(fa) & np.isfinite(fb) & (fb > 0.0)
    if ok.sum() < 3:
        raise TransmissionUnavailable(
            f"{holding_id} @ {centre_A:.3f} A: fewer than 3 pixels where the corrected "
            f"sibling is positive, so no ratio can be formed.")
    t = np.clip(fa[ok] / fb[ok], 0.0, None)
    return t, (f"T = {raw} / {sibling} over {centre_A - pad_A:.2f}-{centre_A + pad_A:.2f} A, "
               f"{int(ok.sum())} px; min T={t.min():.4f}, median T={np.median(t):.4f}")


_SNR_CACHE: dict = {}


class ContinuumBelowScienceFloor(RuntimeError):
    """This spectrum cannot reach the science S/N floor even through a clear sky.

    🔴 A SEPARATE FAILURE FROM `SATURATED`, and keeping them apart is the whole discipline
    of this ticket. `saturated_min_depth = 1 - (SNR_min/S0)^2` goes NEGATIVE when
    S0 < SNR_min, which would classify every line in the window -- including ones under a
    completely clear sky -- as `telluric_saturated_core`. That is a false accusation
    against the atmosphere for what is a property of the DATA, and it is the same
    label-vs-reality defect as recording an unobservable line "ungraded".
    """


def band_continuum_snr(harness, instrument: str, holding_id: str, lo_A: float,
                       hi_A: float, pad_A: float, n_probes: int = 16) -> tuple[float, str]:
    """(S0, evidence) for this holding IN THIS BAND -- measured where the sky is quiet.

    🔴 SCOPE AND SITE BOTH MATTER, and the first two attempts got one each wrong.

    PER-WINDOW WAS WRONG. Taking the noise from each line's own window self-defeats inside
    a telluric forest: H2O lines are narrow, so successive differences count them as noise,
    S0 collapses, and the saturation edge collapses with it. Fe I 8198.921 -- a 0.53-deep
    line in the H2O 8100-8400 complex, squarely the RECOVERABLE class this ticket exists to
    rescue -- came back SATURATED, its threshold set by the very contamination it was meant
    to adjudicate.

    PER-HOLDING WAS ALSO WRONG. S0 is strongly wavelength-dependent -- measured on our own
    Kitt Peak atlas it runs 85 at 3200 A, 2510 at 6700 A, 316 at 11000 A -- so one number
    for a 2960-13000 A holding is not a property, it is an average of incomparable things.
    A whole-holding scan picked a near-UV probe and set the red-optical threshold from
    blanketed blue noise.

    So the scope is the BAND (the unit the census already runs in, over which S0 varies
    slowly) and the SITE is the cleanest window found IN that band -- candidates spread
    across it, transmission measured for each, the highest minimum transmission wins. No
    wavelength is written down: "measure the noise where the sky is quiet" is the rule, and
    the measured transmission is what says where that is.
    """
    key = (instrument, holding_id, round(lo_A, 1), round(hi_A, 1))
    if key in _SNR_CACHE:
        return _SNR_CACHE[key]

    lo, hi = lo_A + pad_A, hi_A - pad_A
    best = None
    for i in range(n_probes):
        centre = lo + (hi - lo) * (i + 0.5) / n_probes
        try:
            t, _ev = measure_transmission(harness, instrument, holding_id, centre, pad_A)
            win = harness.load_window_ex(instrument, centre, pad_A, holding=holding_id,
                                         allow_uncorrected=True)
            snr = continuum_snr(win.flux)
        except Exception:
            continue
        if best is None or float(t.min()) > best[0]:
            best = (float(t.min()), snr, centre)
    if best is None:
        raise TransmissionUnavailable(
            f"{holding_id} in {lo_A:.0f}-{hi_A:.0f} A: no probe window could be both read "
            f"and ratioed against a corrected sibling, so the continuum S/N is "
            f"undetermined and neither observability edge can be derived. Refusing rather "
            f"than assuming an S/N.")
    cleanliness, snr, centre = best
    floor = snr_floor()
    if snr < floor:
        raise ContinuumBelowScienceFloor(
            f"{holding_id} in {lo_A:.0f}-{hi_A:.0f} A: continuum S/N is {snr:.0f} at the "
            f"cleanest window found ({centre:.1f} A, min T={cleanliness:.4f}), below the "
            f"declared floor of {floor:.0f} (PIPELINE['snr_min_science']). No telluric "
            f"verdict is issued here: even a perfectly clear sky would not reach the "
            f"science floor, so calling these lines telluric-saturated would blame the "
            f"atmosphere for a property of the data.")
    ev = (f"S0={snr:.1f} measured at {centre:.1f} A, the cleanest of {n_probes} probes "
          f"across {lo_A:.0f}-{hi_A:.0f} A (min T there = {cleanliness:.4f})")
    _SNR_CACHE[key] = (snr, ev)
    return snr, ev


def observe(harness, instrument: str, holding_id: str, wavelength_A: float,
            pad_A: float, band_A: tuple[float, float]) -> Observability:
    """The per-line verdict. THE entry point -- measured depth in, three states out."""
    corrected = is_corrected(holding_id)
    probe = holding_id
    if corrected:
        # Measure the atmosphere on the holding that still carries it, then report the
        # verdict for the CORRECTED product -- same sky, same window.
        for r in _registry():
            if (r["instrument_id"] == instrument
                    and r.get("telluric_applied", "").strip() == "not-applied"):
                probe = r["holding_id"]
                break
    t, evidence = measure_transmission(harness, instrument, probe, wavelength_A, pad_A)

    snr, snr_ev = band_continuum_snr(harness, instrument, probe,
                                     band_A[0], band_A[1], pad_A)
    th = thresholds(snr)

    tmin = float(t.min())
    depth = 1.0 - tmin

    if depth <= th.clean_max_depth:
        verdict, reason = CLEAN, ""
    elif depth >= th.saturated_min_depth:
        verdict = SATURATED
        reason = SATURATED_CORE_REASON
    else:
        verdict, reason = RECOVERABLE, ""

    if verdict == SATURATED:
        disposition = EXCLUDE
    elif verdict == RECOVERABLE and not corrected:
        disposition = NEEDS_CORRECTION
    else:
        disposition = MEASURE

    return Observability(
        holding_id=holding_id, wavelength_A=float(wavelength_A), verdict=verdict,
        disposition=disposition, depth=depth, transmission_min=tmin, thresholds=th,
        corrected=corrected, reason=reason,
        evidence=f"{evidence}; {snr_ev}; {th.cite()}")


# ── OBSERVABILITY BEFORE TIER (spec 3) ───────────────────────────────────────

@dataclass
class Partition:
    """A pool split by OBSERVABILITY, before any grade is looked at.

    `unmeasured` is deliberately not called `ungraded`. Nothing here has been graded yet;
    that is the point of the ordering.
    """
    measurable: list = field(default_factory=list)
    needs_correction: list = field(default_factory=list)
    unmeasured: list = field(default_factory=list)
    uncharacterised: list = field(default_factory=list)

    def summary(self) -> str:
        return (f"CLEAN/RECOVERABLE-corrected={len(self.measurable)} "
                f"RECOVERABLE-uncorrected={len(self.needs_correction)} "
                f"SATURATED(unmeasured)={len(self.unmeasured)} "
                f"uncharacterised={len(self.uncharacterised)}")


def partition_pool(harness, instrument: str, holding_id: str, wavelengths,
                   pad_A: float, band_A: tuple[float, float]) -> Partition:
    """Split a candidate pool by observability. RUN THIS BEFORE TIER ASSIGNMENT.

    🔴 THE ORDERING IS THE POINT. A SATURATED line must never reach the grader, because
    whatever the grader writes about it will be a statement about its gf -- and the truth
    is that we could not look. It leaves here as `unmeasured` carrying
    `telluric_saturated_core` and its measured depth, and a downstream tier pool that
    contains it is a bug (see the tests).
    """
    p = Partition()
    for w in wavelengths:
        try:
            o = observe(harness, instrument, holding_id, float(w), pad_A,
                        band_A)
        except (TransmissionUnavailable, ContinuumBelowScienceFloor,
                LookupError) as exc:
            p.uncharacterised.append((float(w), f"{type(exc).__name__}: {exc}"))
            continue
        if o.verdict == SATURATED:
            p.unmeasured.append(o)
        elif o.disposition == NEEDS_CORRECTION:
            p.needs_correction.append(o)
        else:
            p.measurable.append(o)
    return p


def problem_child_row(o: Observability, species: str, governing: str = "RYA-1079") -> dict:
    """The `problem_children.csv` row for a SATURATED line -- carrying the measurement.

    Schema is `pipeline.problem_children.SCHEMA_COLUMNS`; the measured depth travels in
    `notes` so the exclusion can be re-checked rather than believed (RYA-711/807).
    """
    return {
        "species": species,
        "lambda_or_scope": f"{o.wavelength_A:.3f}",
        "problem_class": "DATA_GAP",
        "required_treatment": "exclude",
        "observed_in": o.holding_id,
        "amplifies_with": "",
        "severity": "high",
        "governing_tickets": governing,
        "status": "active",
        "population_source": "measured",
        "notes": (f"{SATURATED_CORE_REASON}: measured telluric transmission over the line "
                  f"window reaches T={o.transmission_min:.4f} (depth {o.depth:.4f}), at or "
                  f"past the point where correction is ill-conditioned for this holding. "
                  f"{o.evidence}. Excluded because the photons are gone, NOT because the "
                  f"gf is weak -- this line was never graded (RYA-1079)."),
    }
