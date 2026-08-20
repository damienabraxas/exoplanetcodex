#!/usr/bin/env python3
"""Derive per-(band × treatment) abundance products from measured EWs — RYA-712/713.

    python3 scripts/derive_band_products.py --element Fe --ion I --lo 6910 --hi 9199

Ryan: *"run the IR products now, and the UV"* · *"the two engines should never be added
together and presented, they are separate data products"* · *"diff error bars for diff
bands."*

WHAT THIS PRODUCES
------------------
One row per **treatment**, never merged:

  * **1D-LTE**   — the measured EW inverted through a curve of growth.
  * **ENGINE-A** — 1D-LTE plus the Bergemann/MPIA per-line departure correction.
  * **ENGINE-B** — the spectrum fitted directly by `pipeline.measure.SynthesisHandler`.
                   NO equivalent width is inverted anywhere on this path: 1D-LTE and
                   ENGINE-A invert a measured EW, ENGINE-B fits flux. The EW table supplies
                   only the LINE SET and a per-line strength hint that sets each fitting
                   window; it never supplies the answer.

                   WHICH synthesis is selected by `--engine-b-deck` (RYA-784):
                     `ts-lte`       Turbospectrum LTE flux-fit. What the register calls
                                    Engine B in production for all 27 species. Default.
                     `gerber-nlte`  the TS-native Gerber departure deck (RYA-533/785), wired by
                                    RYA-798. Runs on MARCS.GES and emits a SEPARATE
                                    ENGINE-B-NLTE product (RYA-712).
                                    it today fails LOUDLY rather than silently returning the
                                    LTE number under an NLTE label, which is the one
                                    confusion this flag exists to prevent.

Each carries its own value, line count, and a full `ErrorBudget` with statistical and
systematic kept apart. There is no combined product: `pipeline.band_products` has no
`combine()` and `ErrorBudget.total()` returns a pair.

WHERE THE ABUNDANCES COME FROM
------------------------------
The EWs are whatever the band's controlled handler produced (`*_PROFILEFIT_ew.csv`). The
EW→abundance inversion is the project's own `_bisect_synth_abundance` — Turbospectrum
curve-of-growth, the same engine the EW path uses — not a formula written here.

WHAT IT REFUSES
---------------
* A band whose handler has not passed its optical control (control/frontier rule).
* A line that the measurement quarantined: those are carried into the output with their
  reason and excluded from the aggregate, never dropped (RYA-711).
* An Engine-A product for lines MPIA does not serve: reported as reduced coverage, with
  the count, rather than silently averaging fewer lines.
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.band_policy import resolve as resolve_band  # noqa: E402
from pipeline.treatment_axes import axes_for  # RYA-906
from pipeline.band_products import build_product, LineMeasurement, products_frame  # noqa: E402
from pipeline.error_budget import build as build_budget  # noqa: E402
from pipeline.error_budget import assert_stat_publishable  # noqa: E402  RYA-907
from pipeline import gf_rung  # noqa: E402  RYA-855
from pipeline import harness_residual  # noqa: E402  RYA-869
from pipeline.fit_constraint import as_float_or_none as _f  # noqa: E402  RYA-847
from pipeline.constraint_gate import verdict as constraint_verdict  # noqa: E402
from pipeline.constraint_gate import describe as constraint_describe  # noqa: E402

EW_DIR = ROOT / "data" / "measured" / "band_ew"
OUT = ROOT / "data" / "results" / "band_products"

#: RYA-869 — re-exported, NOT restated. The measured optical residual of the profile
#: fitter against the banked HARPS pool (RYA-713) was written out here AND at
#: `rya850_graded_products.py:73` — two homes for one number, the RYA-845 shape — while
#: `scripts/rya855_rung_audit.py` imports it from this module. It is now declared once in
#: `pipeline.harness_residual`, keyed by the HANDLER that earned it rather than by the
#: script that happened to need it first; this name is kept so that import still resolves.
PROFILE_FIT_RESIDUAL_DEX = harness_residual.HANDLER_RESIDUAL_DEX["ProfileFitHandler"]


def mpia_delta(element: str, ion: str, waves: np.ndarray,
               cache: str | None = None) -> dict[float, float]:
    """Per-line Engine-A departure corrections. Missing = not served.

    The MPIA scraper needs `bs4`, which is absent from the Sirius interpreter that has
    iSpec. Rather than installing into a pinned environment mid-run, the query can be run
    where bs4 lives and handed over as a cache -- the numbers are identical, and the file
    records which lines MPIA actually served.
    """
    if cache:
        d = json.loads(Path(cache).read_text())
        return {float(k): float(v) for k, v in d.items()}
    sys.path.insert(0, str(ROOT))
    from scripts.build_nlte_grids_mpia import _submit_batch
    from config.constants import STAR_PARAMS
    p = STAR_PARAMS["solar"]
    node = [dict(name="sun", teff_K=int(round(float(p["teff"]))),
                 logg=float(p["logg"]), feh=0.0)]
    code = {"I": "26.01", "II": "26.02"}.get(ion) if element == "Fe" else None
    if code is None:
        return {}
    out: dict[float, float] = {}
    for i in range(0, len(waves), 40):          # MPIA dislikes very long batches
        chunk = [float(w) for w in waves[i:i + 40]]
        try:
            r = _submit_batch(chunk, code, node).get("sun", {})
        except Exception as e:
            print(f"    MPIA batch {i//40} failed: {type(e).__name__}: {str(e)[:80]}")
            continue
        for w in chunk:
            if not r:
                continue
            k = min(r, key=lambda x: abs(x - w))
            v = r[k]
            if abs(k - w) < 0.6 and v is not None and np.isfinite(v) and v != 0.0:
                out[w] = float(v)
    return out


# ── the synthesis route (RYA-832) ─────────────────────────────────────────────
#
# The rest of this driver is EW-keyed: it reads measured equivalent widths, inverts them
# through a curve of growth, and keys every synthesis window on the measured EW. The
# near-UV cannot be driven that way and `pipeline.band_policy` says so at intake — it
# FORBIDS profile-fit and interval-integration in 3000-3800 A, because at a median line
# gap of 0.146 A there is no interval that contains one profile and excludes its
# neighbours. RYA-759 tested profile fitting there anyway and falsified it: 901 Fe I
# candidates, 0 measurable.
#
# So the near-UV product came out on its own artifact and never entered the frontier
# matrix — present, but not first-class. This route closes that.
#
# IT CALLS RYA-759's VALIDATED ROUTE RATHER THAN REIMPLEMENTING IT. `select_lines` and
# `fit_one` are imported from that harness, not copied: a hand-adapted copy keeps its
# source's identity in places nobody looks (RYA-701, where one Ba->Al copy produced 13
# defects), and here it would also let the wired value drift from the published 7.487
# silently. Reusing the functions makes drift impossible by construction.

#: RYA-759's own defaults, named here so a reader can see they are NOT re-chosen.
#: Changing any of these would move the product away from the published 7.487 and the
#: ticket forbids that without a stated cause — so they are imported-in-spirit, and the
#: functions below are imported literally.
NEARUV_HALF_WIDTH_A = 0.40      # no EW exists to key the production wing-wide rule
NEARUV_MIN_SEP_A = 4.0          # keeps the set spread instead of piling into one complex
NEARUV_N_LINES = 40
#: 🔴 THE PSEUDO-CONTINUUM SYSTEMATIC IS DELIBERATELY NOT DECLARED HERE (RYA-845).
#: It lives in exactly one place, `pipeline/error_budget.py`, which adds it for any band
#: whose policy says "pseudo-continuum". A second declaration next to a route is what let
#: the term be applied twice: the route could not see what the budget already held, and
#: the number agreed with itself in both places while the product was wrong.


@dataclass(frozen=True)
class SynthBand:
    """Everything the synthesis route needs that is a property of the BAND — RYA-837.

    The route itself is band-agnostic: `select_lines`, `fit_one` and
    `build_solar_context` were already fully parameterised by RYA-759. Only the
    CONSTANTS were near-UV, sitting as module globals. Adding the IR band by copying
    the route would have duplicated 100 lines to change four numbers — the RYA-701
    failure mode (one Ba->Al copy, 13 defects). They are a lookup instead.

    ⚠️ THERE IS DELIBERATELY NO `pseudo_continuum_dex` FIELD (RYA-845, merge of main).
    An earlier cut of this dataclass carried one and set it to 0.100 for all three
    bands. Nothing ever read it — the route used the module global directly — so it
    stated an intent the code did not honour, which is the same second-declaration
    shape that let the term be added twice. `pipeline/error_budget.build()` owns it and
    adds it for any band whose policy says "pseudo-continuum".

    That is the near-UV only. red-optical and NIR get NO continuum systematic, and
    RYA-843 measured the evidence that they should: their fit windows sit at median
    flux 0.73-0.95 against a synthesis normalised to unity, and the fitter spends
    A(Fe) closing that gap. Wiring the field would have hidden the gap behind a
    number nobody derived; it is left OWED and visible instead.
    """
    linelist: Path
    half_width_A: float
    min_sep_A: float
    n_lines: int
    half_width_note: str
    build_hint: str


#: A fixed synthesis half-width has to contain the PROFILE, so it scales with wavelength:
#: the Doppler width is lambda*v/c, and 0.40 A at 3400 A is ~17 Doppler widths while the
#: same 0.40 A at 12000 A is only ~5. Holding the ANGSTROM value fixed across a 3.5x
#: wavelength lever would quietly clip the IR wings and bias A(Fe) low. The IR values
#: below are 0.40 A scaled by lambda, rounded — and swept at runtime (`--sweep-half-width`)
#: rather than asserted, exactly as RYA-759 swept 0.25/0.40/0.60 in the near-UV.
#:
#: The band's own line density says the wider window is affordable: median gap 0.146 A in
#: the near-UV, but 1.872 A (red-optical) and 3.989 A (NIR) per `pipeline.band_policy`.
SYNTH_BANDS: dict[str, SynthBand] = {
    "near-UV": SynthBand(
        linelist=ROOT / "data" / "linelists" / "ispec_nearuv_3000_3780" / "atomic_lines.tsv",
        half_width_A=NEARUV_HALF_WIDTH_A,
        min_sep_A=NEARUV_MIN_SEP_A,
        n_lines=NEARUV_N_LINES,
        half_width_note="RYA-759's own value; swept 0.25/0.40/0.60 there (spread 0.070 dex)",
        build_hint="build it with `pipeline.nearuv_linelist.build()`",
    ),
    "red-optical": SynthBand(
        linelist=ROOT / "data" / "linelists" / "ispec_ir_9200_13000" / "atomic_lines.tsv",
        half_width_A=1.10,          # 0.40 * (9600/3400), rounded
        min_sep_A=4.0,
        n_lines=40,
        half_width_note="0.40 A scaled by lambda from the near-UV anchor (9600/3400)",
        build_hint="RYA-762's VALD extract; see data/linelists/ispec_ir_9200_13000/",
    ),
    "NIR": SynthBand(
        linelist=ROOT / "data" / "linelists" / "ispec_ir_9200_13000" / "atomic_lines.tsv",
        half_width_A=1.40,          # 0.40 * (12000/3400), rounded
        min_sep_A=4.0,
        n_lines=40,
        half_width_note="0.40 A scaled by lambda from the near-UV anchor (12000/3400)",
        build_hint="RYA-762's VALD extract; see data/linelists/ispec_ir_9200_13000/",
    ),
}


def synthesis_route(a, pol) -> None:
    """Derive a band product by SYNTHESIS, for a band where EW measurement is undefined.

    This is a THIN WRAPPER over RYA-759's harness, not a reimplementation of it. Every
    scientific choice — the line list, the gf-provenance decision, the solar context, the
    atlas segments, the selection rule, the fit — is the function 759 validated, called
    with 759's own defaults. The only thing added here is the wrapping of its per-line
    results into `LineMeasurement`/`build_product` so the near-UV becomes a first-class
    matrix cell (RYA-712) instead of a side artifact.

    That structure is deliberate. A hand-adapted copy would keep its source's identity in
    places nobody looks (RYA-701: one Ba->Al copy produced 13 defects), and here it would
    additionally let the wired value drift from the published 7.487 without anyone
    noticing. Calling the same functions makes that drift impossible by construction.
    """
    if "synthesis" not in pol.permitted_methods:
        raise SystemExit(
            f"band {pol.name!r} does not permit synthesis (permitted: "
            f"{pol.permitted_methods}) — refusing to derive a product by a method the "
            f"band policy forbids")

    cfg = SYNTH_BANDS.get(pol.name)
    if cfg is None:
        raise SystemExit(
            f"no synthesis configuration for band {pol.name!r} — refusing to invent a "
            f"half-width and line-selection rule for a regime nobody has characterised. "
            f"Add an entry to SYNTH_BANDS with its reasoning.")

    sys.path.insert(0, str(ROOT / "scripts"))
    from pipeline.nearuv_synth import build_solar_context, gf_provenance
    from rya759_nearuv_fe_product import select_lines, fit_one, species_token
    # RYA-913: `_kp_segments` is a PRIVATE alias of the same Kitt-Peak-specific
    # helper. A guard that bans only the public names would have passed this by,
    # which is why the Part B guard matches behaviour, not spelling.
    from rya759_nearuv_synth import _kp_segments
    from measure_band_ew import load_window_ex, select_holding
    from pipeline.telluric_policy import applied_state as _applied_state

    if not cfg.linelist.exists():
        # Loud-fail, never a silent default (RYA-832). A missing list would otherwise
        # surface downstream as "0 lines measured", which reads like a physics result.
        raise SystemExit(
            f"{pol.name} line list missing at {cfg.linelist}\n"
            f"{cfg.build_hint}. Refusing to emit an empty {pol.name} product that "
            f"would read as a measurement of nothing.")

    cat = pd.read_csv(ROOT / "data" / "catalog" / "instrument_catalog.csv")
    row = cat[cat.iloc[:, 0].astype(str) == a.instrument]
    if row.empty:
        raise SystemExit(f"{a.instrument!r} absent from data/catalog/instrument_catalog.csv")
    R = float(row.iloc[0]["resolving_power_max"])

    # Ask about the lines that will actually be SYNTHESISED, not about the nominal band
    # bound (RYA-837). `--hi 12935` is a request; this band's reddest real line is
    # 12934.67. Passing the nominal 12935.0 put it 0.33 A past canonical_gf's last entry
    # and switched single-sourcing OFF for the WHOLE band over a sliver containing no
    # lines — discarding RYA-834's 28 lab-adjudicated gf in the process. The blueward
    # side has the same shape and absorbs it with a 0.01 A tolerance justified against a
    # ~0.18 A line spacing; that constant cannot be reused here, where the NIR spacing is
    # 3.99 A. Using the real extent needs no tolerance at all.
    _ll = pd.read_csv(cfg.linelist, sep="\t", usecols=["wave_A"], low_memory=False)
    _w = _ll.wave_A[(_ll.wave_A >= a.lo) & (_ll.wave_A <= a.hi)]
    if _w.empty:
        raise SystemExit(
            f"no line in {cfg.linelist.name} lies within {a.lo}-{a.hi} A — refusing to "
            f"state gf provenance for a band this list does not cover.")
    prov_gf = gf_provenance(float(_w.min()), float(_w.max()))
    print(f"\n[synthesis route — {pol.name}]  R={R:.0f}")
    print(f"  [gf] {prov_gf['detail']}")
    ctx = build_solar_context(a.element, R, linelist_file=str(cfg.linelist),
                              apply_canonical_gf=prov_gf["apply_canonical_gf"])
    segs = _kp_segments()
    hw = float(getattr(a, "half_width_A", None) or cfg.half_width_A)

    # ── 🔴 WHICH SPECTRUM DOES `--instrument` ACTUALLY MEAN? — RYA-904 ────────
    # It meant NOTHING until now. This route tagged the product with `a.instrument` and
    # then called `fit_one`, which read Kitt Peak unconditionally. `--instrument
    # crires_plus` would have emitted a CRIRES+-labelled product measured on the KPNO
    # atlas — and converged, because KPNO covers these windows perfectly well. Same
    # defect class as the loader dispatch this ticket fixes: a product naming a source it
    # was not measured from.
    #
    # The window now comes from the HOLDING dispatch, which selects per line, refuses
    # rather than falling back, and names what it served. Two properties are asserted
    # here rather than hoped for:
    #
    #  1. PRE-NORMALISED. `_fit_synth_flux` compares observed flux against a synthesis
    #     normalised to unity, so an un-normalised holding (the raw CRIRES+ Vesta IDPs,
    #     adu to ~1.9e5) is not merely inaccurate on this route, it is undefined. Refused
    #     by name rather than fitted — the RYA-713 direction that matters here.
    #  2. ONE HOLDING FOR THE WHOLE BAND. Per-line selection could in principle serve
    #     line 1 from one holding and line 2 from another; the product would then be an
    #     average over two normalisation histories with one label. Recorded and refused.
    _served: dict[str, int] = {}
    _served_specs: dict[str, object] = {}
    #: median observed flux per fit window. MEASURED, because the alternative was to
    #: quote RYA-843's 0.73-0.95 — which it measured on the telluric-UNCORRECTED Kitt
    #: Peak NIR windows. Transplanting that number onto a corrected holding is exactly
    #: the kind of borrowed claim this ticket exists to stop.
    _window_median: list[float] = []

    def _observed(centre: float, pad: float):
        win = load_window_ex(a.instrument, centre, pad)
        h = win.holding
        if not h.pre_normalised:
            raise LookupError(
                f"holding {h.holding_id} is NOT continuum-normalised, and the synthesis "
                f"route fits observed flux against a synthesis normalised to unity. "
                f"Fitting it would spend A({a.element}) closing a continuum offset "
                f"(RYA-713/843). Normalise the product first.")
        _served[h.holding_id] = _served.get(h.holding_id, 0) + 1
        _served_specs[h.holding_id] = h
        _window_median.append(float(np.median(win.flux)))
        return win.wave, win.flux, win.provenance

    _probe = select_holding(a.instrument, 0.5 * (a.lo + a.hi), hw + 0.4)
    # Asked at band centre so the refusal costs nothing and arrives BEFORE any synthesis.
    # `_observed` checks it too, per line, because a band can in principle select a
    # different holding line by line -- but a per-line failure there surfaces as
    # `status: no_atlas` on each line, which describes coverage rather than the real
    # fault. The loud version comes first (RYA-832).
    if not _probe.pre_normalised:
        raise SystemExit(
            f"{a.instrument} -> {_probe.holding_id} is NOT continuum-normalised, and "
            f"this route fits observed flux against a synthesis normalised to unity. "
            f"The fit would spend A({a.element}) closing a continuum offset rather than "
            f"measuring an abundance (RYA-713/843). Normalise the product first, or "
            f"derive this band through a route that sets its own continuum.")
    print(f"  [holding] {a.instrument} -> {_probe.holding_id} "
          f"(pre_normalised={_probe.pre_normalised})")
    if _probe.note:
        print(f"            {_probe.note}")
    # RYA-904 — the SPECIES, not just the element. `--ion` was decorative on this
    # route: `select_lines` was pinned to Fe I, so `--ion II` returned Fe I lines and
    # every one was labelled Fe II. `species_token` spells the ion the way the linelist
    # does ('Fe 2', not 'Fe II') and the selection refuses loudly when the band holds
    # none of that species — which is the honest answer where it holds none.
    cand = select_lines(ctx["linelist"], lo_A=a.lo, hi_A=a.hi, n=cfg.n_lines,
                        teff=float(ctx["teff"]), min_sep_A=cfg.min_sep_A,
                        species=species_token(a.element, a.ion))
    # 🔴 THE TELLURIC MASK BELONGS AT SELECTION, AND THIS ROUTE NEVER APPLIED IT
    # (RYA-843). The EW route calls `telluric_reason` and skipped 29 NIR lines; the
    # synthesis route inherited `select_lines`, which knows nothing about tellurics. So
    # the first NIR run fitted lines sitting inside the H2O 11120-11560 A band, where the
    # atlas flux is not stellar at all -- measured observed depths of 0.99, 1.00, and one
    # at 1.001, i.e. FLUX AT OR BELOW ZERO. There is no Fe profile to fit there, so the
    # optimizer ran the abundance to its bound and returned ~12.49.
    #
    # That is a LINE-SELECTION failure, not a fitter one (RYA-842: line selection
    # dominates). Fixing it here makes the rail-detector tolerance moot for these lines
    # instead of tuning a threshold to hide them. RYA-762's own atlas survey already said
    # the two H2O bands are severe and masked; the mask simply never reached this path.
    from measure_band_ew import telluric_reason
    _tell = [(float(r.wave_A), telluric_reason(float(r.wave_A), a.instrument))
             for r in cand.itertuples()]
    _bad = [(w, why) for w, why in _tell if why]
    if _bad:
        print(f"  {len(_bad)} candidate(s) QUARANTINED-TELLURIC before fitting "
              f"(flux there is not stellar; excluded, never corrected):")
        for w, why in _bad:
            print(f"      {w:10.3f}  {why.split(' and ')[0]}")
        keep = {w for w, why in _tell if not why}
        cand = cand[cand.wave_A.astype(float).isin(keep)].reset_index(drop=True)
    print(f"  {len(cand)} {a.element} {a.ion} candidates by theoretical depth "
          f"(half-width +/-{hw} A, min separation {cfg.min_sep_A} A)")
    print(f"  [half-width] {cfg.half_width_note}")

    tmp = f"/tmp/ispec_synth_{pol.name.replace(chr(47), chr(45))}"
    Path(tmp).mkdir(parents=True, exist_ok=True)
    lines: list[LineMeasurement] = []
    for r in cand.itertuples():
        w = float(r.wave_A)
        res = fit_one(ctx, segs, w, hw, tmp, load=_observed)
        a_x = float(res.get("a_synth", float("nan")))
        lm = LineMeasurement(
            element=a.element, ion=a.ion, wavelength_air_A=w,
            instrument=a.instrument, ew_mA=float("nan"),
            ew_method=(f"synthesis flux-fit, FIXED half-width +/-{hw} A "
                       f"(RYA-759 route; no EW exists in this band to key the "
                       f"production wing-wide rule); observed from "
                       # RYA-904 — the per-line record names the HOLDING it was fitted
                       # against, not just the instrument. `crires_plus` alone does not
                       # distinguish the corrected Y arm from the raw Vesta IDPs.
                       f"{res.get('obs_source', 'UNRECORDED')}"),
            abundance=(a_x if np.isfinite(a_x) else None),
            treatment="1D-LTE",
            # RYA-880 — an LTE row states that it is LTE. A blank cannot distinguish
            # "no departure applied" from "nobody recorded one" (RYA-833).
            nlte_delta_dex=0.0, nlte_source=LTE_NLTE_SOURCE,
            # RYA-871 — `select_lines` already returns the EP of the row it picked, and
            # this route picks AT the list's own wavelength, so the identity is exact
            # rather than recovered. Carried anyway, for two reasons: the near-UV pool
            # still had 2 lines with two same-species rows inside 0.005 A, which only an
            # EP key separates; and a per-line artifact whose identity column is
            # populated on one route and blank on another is a schema a consumer has to
            # special-case.
            ep_eV=float(r.ep_eV),
            # The REW saturation ceiling is an EW-INVERSION concept and there is no EW
            # here at all. Inheriting it would quarantine lines on a quantity that does
            # not exist (RYA-770/342).
            ew_inversion=False,
            # RYA-847 — carried on every line, and now GATED on one of them. The sweep
            # across nine cells found no transferable threshold for any of these four
            # (per-band cuts span x24-x489), so no metric carries a cut; what does gate
            # is the sign of `frac_rise_weaker` (see below). They stay on the line
            # because the sweep is only reproducible from per-line metrics.
            sigma_A=_f(res.get("sigma_A")),
            frac_rise_weaker=_f(res.get("frac_rise_weaker")),
            edge_distance_dex=_f(res.get("edge_distance_dex")),
            red_chi2=_f(res.get("red_chi2")))
        # 🔴 A NON-'ok' FIT IS NOT A MEASUREMENT — RYA-837.
        # RYA-759's own harness aggregates `status == 'ok'` only, and this route dropped
        # that filter when it lifted the functions. The near-UV happened not to notice;
        # the IR did, loudly. `_fit_synth_flux` bounds the search at a_hi = A_solar + 5
        # and returns THE EDGE VALUE for a line that runs to it, flagged `edge_pinned`
        # — "carry the edge value and are marked, never substituted", per its docstring.
        # Aggregated anyway, those pile up at ~12.49 and wreck the dispersion: the first
        # NIR run reported 7.588 +/- 1.977 dex, where the 1.977 was SEVEN railed lines
        # and the honest MAD over the real ones was 0.416.
        #
        # A railed fit means the window could not constrain the abundance from BELOW —
        # it is an absence of information, and averaging it in states the opposite.
        status = str(res.get("status", "failed"))
        if lm.abundance is None or status != "ok":
            lm.in_aggregate = False
            lm.excluded_reason = (f"SYNTHESIS: {status} "
                                  f"{str(res.get('reason', ''))[:60]}").strip()
        else:
            # RYA-847 — THE BYPASS, CLOSED. This route's entire accept/reject was the
            # `status != "ok"` above: it never consulted `red_chi2` and never asked
            # whether the fit constrained anything, which is how two lines whose chi2
            # moves 2.2% and 1.4% across EIGHT DEX of iron entered the published
            # aggregate at 7.833 and 7.979 (RYA-843).
            #
            # It now calls the SAME decider the Engine-B handler calls. ⚠️ THIS IS NO
            # LONGER NUMERICALLY INERT and the comment here used to say it was: with the
            # non-minimum check adopted, this route drops 3617.318 and the near-UV
            # product moves 7.488 -> 7.498 (n 40 -> 39). `SYNTH_CONSTRAINT` is still None
            # and always will be — no threshold survived the sweep — but `frac_rise <= 0`
            # is a correctness test, not a threshold, and it fires here.
            _cv = constraint_verdict(res)
            if not _cv.ok:
                lm.in_aggregate = False
                lm.excluded_reason = _cv.reason
        lines.append(lm)
    lines.sort(key=lambda l: (l.wavelength_air_A, l.element, l.ion))

    # 🔴 THE gf RUNG IS DECIDED FROM THE LINES, NOT HARDCODED — RYA-855.
    # This route passed `gf_graded=False` to `build_budget` unconditionally, so the
    # product charged the ungraded Kurucz 0.17 whatever its lines actually were. The
    # decision now comes from `pipeline.gf_rung`, which the EW route below calls too:
    # one implementation, so the two routes cannot disagree (RYA-845/847 defect shape).
    # It is computed HERE, before the provenance, because the provenance used to STATE
    # the rung in prose — "THE BAND IS UNGRADED" — as a second, hand-maintained
    # declaration of the same fact the budget was charging. It now quotes the decision.
    rung = gf_rung.for_lines(a.element, a.ion, lines, linelist=ctx["linelist"])
    print(f"  [gf rung] {rung.describe()}")

    # 🔴 THE PROSE SAID "near-UV" FOR EVERY BAND — RYA-904.
    # This string is `Product.provenance`, and it is what a reader gets. Written for
    # RYA-759 and never re-read when RYA-837 made the route band-agnostic, it told a NIR
    # product that its median line gap was 0.146 A (it is 3.989), that RYA-759 had
    # falsified profile fitting "in band" (a different band), that its Fe I was UV
    # over-ionised, and — worst, because it is arithmetic and not colour — that it
    # carried a 0.100 dex PSEUDO-CONTINUUM SYSTEMATIC. `error_budget.build()` charges
    # that term only where the BandPolicy says "pseudo-continuum", which is the near-UV
    # alone; the NIR budget does not carry it. So the provenance asserted a term the
    # product was not charged, which is the RYA-845 shape (two declarations of one fact)
    # pointing the other way. Band-specific claims are now sourced from the band.
    _band_specific = {
        "near-UV": (
            "SYNTHESIS-ONLY: band_policy forbids profile-fit and interval-integration "
            "here (median line gap 0.146 A leaves no interval that contains one profile "
            "and excludes its neighbours), and RYA-759 falsified profile fitting in band "
            "— 901 candidates, 0 measurable. 1D-LTE ONLY: UV Fe I is heavily "
            "over-ionised, so the missing NLTE correction is large and POSITIVE, and a "
            "low value here is expected physics rather than a defect. "
            "PSEUDO-CONTINUUM SYSTEMATIC 0.100 dex, which does NOT average down and is "
            "NOT in the scatter reported here. "),
    }.get(pol.name,
          f"SYNTHESIS-ONLY: band_policy forbids {'/'.join(pol.forbidden_methods)} in "
          f"{pol.name} — {pol.justification} "
          # RYA-904 — READ THE POLICY AGAINST THE HOLDING. The NIR justification is
          # written about UNCORRECTED data ("until a telluric correction is applied the
          # observed flux is not a stellar spectrum"), which is precisely the condition
          # a telluric_applied=applied holding satisfies. Quoting it beside a corrected
          # product without saying so reads as if the product were uncorrected. The
          # policy is NOT relaxed here — its route choice (synthesis) still stands, and
          # changing a band policy on the strength of one holding is a separate decision
          # nobody has taken.
          + ("The clause about correction is SATISFIED by this product rather than "
             "waived: every holding that served it is registered "
             "telluric_applied=applied. The band's ROUTE choice is unchanged — "
             "synthesis — because a band policy is not re-decided by one holding. "
             if pol.telluric_required and _served and all(
                 _applied_state(k) == "applied" for k in _served)
             else ("Note that the holding(s) serving this product are NOT registered "
                   "telluric-corrected, so the clause above applies to them literally. "
                   if pol.telluric_required else ""))
          + f"NO pseudo-continuum systematic is charged "
          f"in this band and none is claimed here. RYA-843 measured the reason one is "
          f"OWED — fit windows sitting at median flux 0.73-0.95 against a synthesis "
          f"normalised to unity — but it measured that on the telluric-UNCORRECTED Kitt "
          f"Peak NIR, so it is NOT quoted for another holding. THIS pool's own windows "
          f"have median observed flux "
          + (f"{np.median(_window_median):.4f} (range {min(_window_median):.4f}-"
             f"{max(_window_median):.4f}, n={len(_window_median)}), measured here. "
             if _window_median else "UNMEASURED (no window loaded). ")
          + "The term stays owed and visible rather than invented. ")
    prov = (
        f"{a.element} {a.ion} {pol.name} 1D-LTE by flux-fit synthesis (Turbospectrum via "
        "iSpec, ATLAS9.Castelli) — the RYA-759 validated route called directly, same "
        "functions and same defaults, so the value cannot drift from that ticket's. "
        # RYA-904 — WHICH HOLDING. `instrument` does not identify the spectrum: CRIRES+
        # has a telluric-corrected Y arm and raw uncorrected Vesta IDPs, and a product
        # measured from either would carry the same instrument tag. The holding is
        # reported from what the loader actually served, per line, not from an intention.
        + "OBSERVED FROM: "
        + ", ".join(f"{k} (n={v})" for k, v in sorted(_served.items()))
        + (" [holding names its telluric state in "
           "data/catalog/holdings_manifest_registry.csv; the gate consumes it and this "
           "route never asserts one of its own, RYA-786/904]. ")
        + _band_specific
        # RYA-904 spec 6 — the holding's own caveat, carried into the product's reason
        # field so it reaches the page. Taken from the holding that ACTUALLY SERVED the
        # lines, not from the probe at band centre: a caveat attached to a holding the
        # product did not use would be worse than none.
        + "".join(getattr(_served_specs[k], "caveat", "") + " "
                  for k in sorted(_served) if getattr(_served_specs[k], "caveat", ""))
        + "gf: " + str(prov_gf["detail"]) + ". "
        + "Half-width is FIXED and must be swept. "
        + constraint_describe() + " " +
        # RYA-855 — the rung is QUOTED from the decider, never restated. The sentence
        # that stood here asserted "THE BAND IS UNGRADED" in prose while the budget
        # asserted it in arithmetic; two declarations of one fact is how RYA-845's
        # double-count survived, and a hand-written one would go stale the first time a
        # pool's lines changed. Context on WHY Fe I in this band is largely ungraded:
        # RYA-822 grades only 6 of the 4,274 in-band Fe I lines as primary-lab, and its
        # GF-NIST class (604 lines) is a COMPILATION grade 822 deliberately keeps
        # outside `is_graded` because FMW *is* NIST and VALD copies it (RYA-760).
        "gf TERM: " + rung.describe() + ". Never coadded with another band (RYA-712).")
    # RYA-904 — ONE PRODUCT, ONE HOLDING. Per-line dispatch could legitimately serve
    # two holdings inside one band (a fixed-span product covering part of it, another
    # covering the rest). The aggregate would then average two normalisation and telluric
    # histories under a single instrument label, and no reader could unpick which line
    # came from which. RYA-712's rule is that different products are never combined; two
    # holdings are two products. Refused rather than reported, because the number would
    # already have been computed by the time anyone read the note.
    if len(_served) > 1:
        raise SystemExit(
            f"this band was served by {len(_served)} holdings — "
            + ", ".join(f"{k} (n={v})" for k, v in sorted(_served.items()))
            + f". One product may not average two holdings (RYA-712/904). Narrow --lo/--hi "
              f"to a range one holding covers, or derive them as separate products.")
    if not _served:
        raise SystemExit(
            "no line in this band was served by any holding, so nothing was measured "
            "against an observed spectrum. Refusing to emit a product (RYA-832).")

    # RYA-869 — the handler is DECLARED by the route that ran it. This route fits flux
    # (`fit_one` -> `_fit_synth_flux`); no profile fitter is on it. Note that its
    # treatment is `1D-LTE`, the same label the VIS EW route wears — that pair is the
    # proof that the handler is not a function of the treatment name, and no wider set
    # of treatment names could have got both right.
    product = build_product(a.element, a.ion, a.instrument, pol.name, "1D-LTE",
                            lines, handler="SynthesisHandler", provenance=prov)

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    stem = f"{a.element}{a.ion}_{int(a.lo)}_{int(a.hi)}_{a.instrument}_SYNTH"
    pd.DataFrame([asdict_line(l) for l in lines]).to_csv(
        out / f"{stem}_1D-LTE_lines.csv", index=False)

    # The products.csv SCHEMA is the matrix's input contract, not this function's choice.
    # `products_frame` emits `value`; the matrix reads `A`, `stat_dex`, `syst_dex`. The
    # first version of this route wrote the products_frame schema and the near-UV cell
    # duly appeared in the matrix with n=40 and a NaN abundance — present, and empty.
    # A row that arrives but carries no number is worse than one that does not arrive,
    # because it looks like a measurement that failed rather than a wiring mismatch.
    b = build_budget(a.element, 0.5 * (a.lo + a.hi), product.n_lines,
                     # 🔴 RYA-907 — this read `(product.sigma or 0.0)`. `build_product`
                     # returns None for n=1 on purpose (the spread of one line is
                     # undefined, not small) and `or 0.0` converted that refusal into a
                     # measured zero. Two published Fe II red-optical cells then carried
                     # stat_dex=0.0 — the tightest statistical bar in the matrix on its
                     # least constrained products. The None is now carried as a None and
                     # the budget records the term UNMEASURED.
                     scatter_dex=product.sigma,
                     # RYA-855 — the gf rung, decided above from this product's own
                     # lines. `budget_kwargs()` carries the graded flag and the cited
                     # sigma TOGETHER so a caller cannot pass one and forget the other.
                     **rung.budget_kwargs(),
                     # RYA-869 — the harness term, decided from THIS PRODUCT'S OWN
                     # handler rather than restated here. Charging the profile fitter's
                     # measured residual to a route it never touched would attribute
                     # someone else's systematic to it; value and LABEL travel together
                     # so the budget file's prose cannot name a different handler than
                     # its arithmetic does.
                     **harness_residual.for_product(product).budget_kwargs())
    stat, syst = b.total()
    assert_stat_publishable(
        stat, cell=f"{a.element} {a.ion} {pol.name} 1D-LTE "
                   f"({a.instrument}, n={product.n_lines})")
    # 🔴 RYA-845: the pseudo-continuum term is NOT added here, and adding it was a bug.
    # `error_budget.build()` already carries it for this band — the near-UV BandPolicy's
    # `continuum_treatment` says "pseudo-continuum only", which fires the branch in
    # error_budget.py. Adding it again here counted it TWICE and inflated the published
    # near-UV systematic 0.1972 -> 0.2211.
    #
    # The comment that used to sit here said the term "is NOT in the line scatter" —
    # which is true, and is not the question. It was already in the BUDGET. The budget
    # owns its terms; a route that hand-adds one cannot know what the budget already has.
    pd.DataFrame([dict(
        element=a.element, ion=a.ion, band=pol.name, instrument=a.instrument,
        treatment="1D-LTE", handler=product.handler,
        A=round(product.value, 3) if product.value is not None else None,
        n_lines=product.n_lines, n_excluded=product.n_excluded,
        stat_dex=round(stat, 4), syst_dex=round(syst, 4),
        # RYA-907 — the bar NAMES ITS OWN ORIGIN. A reader of the matrix could not
        # previously tell a measured scatter from the quantiser floor from a stand-in
        # for a scatter nobody measured, and those three are different claims.
        stat_basis=b.stat_basis(),
        dominant=(b.dominant().name if b.dominant() else ""),
        # 🔴 RYA-906 CANARY CELL. This route emits treatment="1D-LTE" and it is a
        # SYNTHESIS flux fit — the same label the VIS EW cells use. The axes are resolved
        # from `product.handler` (SynthesisHandler), so this row renders `Synth · 1D-LTE`.
        # If it ever renders `EW`, something started reading route off the label and the
        # migration is wrong.
        **axes_for("1D-LTE", handler=product.handler or None).as_columns(),
    )]).to_csv(out / f"{stem}_products.csv", index=False)
    # `describe()` already lists every term the budget holds, pseudo-continuum included.
    # The epilogue that used to be appended here announced the term as "added in
    # quadrature above" — which is exactly the double-add RYA-845 removed, and it made
    # the artifact assert the bug in prose as well as in arithmetic.
    # RYA-855 — the budget names the TERM; this names the DECISION that chose it. A
    # reader must be able to see not only that the cell was charged 0.17 but why this
    # pool was not entitled to anything better, without re-running the deriver.
    (out / f"{stem}_budgets.txt").write_text(
        b.describe() + f"\n  gf rung: {rung.describe()}\n")
    # RYA-847 — WRITE THE PROVENANCE. `build_product` has always accepted it and
    # `products_frame` emits it, but this route writes the MATRIX schema instead
    # (element/ion/band/.../dominant, RYA-832) and that schema has no provenance column —
    # so every word of reasoning assembled above, including whether a constraint cut was
    # applied at all, reached no artifact. Found by grepping the control run's output for
    # "CONSTRAINT GATE" and getting nothing, after a commit message had already claimed
    # otherwise. Same shape as the defect this ticket exists to fix: a quantity with
    # nowhere to live is a quantity nobody can check.
    (out / f"{stem}_provenance.txt").write_text(product.provenance + "\n")
    v = f"{product.value:.3f}" if product.value is not None else "n/a"
    s = f"{product.sigma:.3f}" if product.sigma is not None else "n/a"
    print(f"\n  A({a.element} {a.ion}; {pol.name}, 1D-LTE) = {v} +/- {s}  "
          f"(n={product.n_lines}, excluded {product.n_excluded})")
    print(f"  NOTE sigma here is the line SCATTER (std, ddof=1) that build_product "
          f"reports for\n  every band; RYA-759 published a MAD. Same lines, different "
          f"dispersion statistic —\n  the VALUE is the number to compare.")
    print(f"  wrote {out}/{stem}_products.csv")


#: Fields of `LineMeasurement` that `asdict_line` deliberately does NOT write, with the
#: reason. Anything not listed here and not emitted makes `test_asdict_line_emits_every_
#: field_or_declares_why` fail — see RYA-847.
#:
#: This list exists because a hand-written projection is a second declaration of the
#: schema, and the first one already cost us: `red_chi2` was returned by the fitter for
#: its whole life and never reached the per-line CSV, because there was no field for it
#: and no check that would notice (RYA-843). Adding a field to the dataclass must not
#: silently fail to reach the artifact.
ASDICT_LINE_OMITTED = {
    # RYA-807's registry verdict. Carried on the object for the aggregation decision;
    # the per-line CSV records the OUTCOME (`in_aggregate` + `excluded_reason`), and
    # duplicating the registry's own columns here would give them a second home that can
    # disagree with the registry.
    "problem_class", "problem_status", "problem_tickets", "problem_action",
}


def _intake_ep(row, wavelength_A: float, element: str, ion: str) -> float:
    """The EP the EW artifact carries for this line — RYA-871.

    🔴 LOUD ON AN OLD ARTIFACT. An EW file written before RYA-871 has no `ep_eV` column,
    and silently emitting None would put every line of that band back on the
    wavelength-only rule while the products.csv showed an `ep_eV` column full of blanks.
    A consumer cannot tell "this route does not carry an EP" from "the EP is missing here"
    once the column exists but is empty — so the missing column is named and refused, and
    the fix is to re-measure the band, not to widen blind.
    """
    ep = getattr(row, "ep_eV", None)
    if ep is None:
        raise SystemExit(
            f"the EW artifact for {element} {ion} carries no `ep_eV` column, so its lines "
            f"cannot be identified by anything but their wavelength (RYA-871). It predates "
            f"RYA-871; re-run scripts/measure_band_profilefit.py for this band.")
    if not np.isfinite(float(ep)):
        raise SystemExit(
            f"{element} {ion} {wavelength_A:.4f} A: the EW artifact carries a null ep_eV. "
            f"A blank identity key is not a key — re-measure rather than emit it.")
    return float(ep)


#: RYA-880 — the NLTE provenance strings, declared ONCE (RYA-845). Each says what was
#: applied AND how, because RYA-207 left "per-line vs mean-field" as a live architectural
#: question and a reader cannot tell which they are looking at from a number alone.
MPIA_NLTE_SOURCE = ("Bergemann MPIA per-line delta_nlte (live query, solar node); "
                    "PER-LINE additive correction")
LTE_NLTE_SOURCE = "none — LTE, no departure applied"
#: ⚠️ ENGINE-B-NLTE has NO additive per-line delta to record. The Gerber departures enter
#: the radiative transfer, and RYA-712 makes this a SEPARATE PRODUCT rather than a
#: correction applied to Engine-B-LTE. Reporting a number here would require differencing
#: the two products, which is precisely what RYA-880 forbids and what RYA-712 forbids
#: conceptually. So the delta is None and the source says why.
GERBER_NLTE_SOURCE_FMT = ("gerber:{atom} (TS-native departure deck, RYA-798); departures "
                          "applied INSIDE the synthesis — no additive per-line delta "
                          "exists, this is a separate product, not a corrected LTE value")


def asdict_line(l: LineMeasurement) -> dict:
    return dict(element=l.element, ion=l.ion, wavelength_air_A=l.wavelength_air_A,
                instrument=l.instrument, ew_mA=l.ew_mA, ew_method=l.ew_method,
                abundance=l.abundance, rew=l.rew, treatment=l.treatment,
                in_aggregate=l.in_aggregate, excluded_reason=l.excluded_reason,
                # RYA-871 — the excitation potential of the transition measured. Without
                # it a per-line row can only be identified by its wavelength, which does
                # not identify a line: 16 of 152 VIS Fe I lines could not be resolved
                # back to the list at all.
                ep_eV=l.ep_eV,
                ew_inversion=l.ew_inversion,
                # RYA-847: the constraint metrics. None on EW-route lines — no chi2
                # surface exists behind an EW inversion, so the question does not apply,
                # and a consumer must read None as "not applicable", not "unconstrained".
                sigma_A=l.sigma_A, frac_rise_weaker=l.frac_rise_weaker,
                edge_distance_dex=l.edge_distance_dex, red_chi2=l.red_chi2,
                # RYA-880 — and it must be added HERE too: this hand-written
                # projection is the lossy one of the object's two serialisers, and
                # RYA-847 already lost `red_chi2` to exactly this omission.
                nlte_delta_dex=l.nlte_delta_dex, nlte_source=l.nlte_source,
                # RYA-911 — the continuum the harness PLACED, carried to the product.
                # The EW artifact records it; a product built from that EW has to carry
                # it or the decisive quantity stops at the intermediate file and the next
                # RCA reconstructs it again.
                continuum_level=l.continuum_level, continuum_method=l.continuum_method,
                continuum_ref=l.continuum_ref,
                # RYA-911 — the fit's width and the floor it was judged against, same
                # reasoning as the continuum directly above: the quantity that DECIDED
                # the width verdict has to reach the product, or the selection-bias
                # question can only ever be asked of the intermediate file.
                profile_sigma_A=l.profile_sigma_A,
                profile_sigma_floor_A=l.profile_sigma_floor_A,
                # RYA-906 — emitted BESIDE sigma, never without it. On its own
                # `profile_sigma_A` cannot be read: the two are degenerate in a Voigt
                # fit, and judging a line by sigma alone is the defect this replaced.
                profile_gamma_A=l.profile_gamma_A)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--element", required=True)
    ap.add_argument("--ion", default="I")
    ap.add_argument("--lo", type=float, required=True)
    ap.add_argument("--hi", type=float, required=True)
    ap.add_argument("--instrument", default="kpno_solar_atlas")
    ap.add_argument("--force-synthesis", action="store_true",
                    help="drive a band through the SYNTHESIS route even where its policy "
                         "also permits profile-fit (RYA-837). Needed for red-optical "
                         "9199.9-10000 A: the EW route's curve-of-growth inversion reads "
                         "the 4200-9199.99 A GES synth linelist and returns "
                         "NOT-IN-SYNTH-LINELIST for every line above 9199.9. Explicit, so "
                         "the method a product used is never a silent consequence of "
                         "which input happened to be missing.")
    ap.add_argument("--half-width-A", type=float, default=None,
                    help="override the band's synthesis half-width (RYA-837). Used to "
                         "SWEEP it and report sensitivity rather than assert one value, "
                         "as RYA-759 swept 0.25/0.40/0.60 in the near-UV.")
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--mpia-cache", help="pre-computed per-line delta_nlte JSON")
    ap.add_argument("--engine-b-deck", choices=["ts-lte", "gerber-nlte"], default="ts-lte",
                    help="which synthesis is Engine B. 'ts-lte' is the production "
                         "Turbospectrum LTE flux-fit; 'gerber-nlte' is the TS-native "
                         "Gerber departure deck (validated RYA-785, wired RYA-798). "
                         "It runs on MARCS.GES to match the deck and emits a SEPARATE "
                         "ENGINE-B-NLTE product -- never a correction to Engine-B-LTE.")
    ap.add_argument("--skip-engine-b", action="store_true",
                    help="derive 1D-LTE and ENGINE-A only. Engine B refits the spectrum "
                         "per line and is much slower than the EW inversion.")
    a = ap.parse_args()

    pol_early = resolve_band(0.5 * (a.lo + a.hi))
    if a.force_synthesis:
        if "synthesis" not in pol_early.permitted_methods:
            raise SystemExit(
                f"--force-synthesis refused: band {pol_early.name!r} does not permit "
                f"synthesis ({pol_early.permitted_methods}). The flag chooses BETWEEN "
                f"permitted methods; it does not overrule the policy.")
        return synthesis_route(a, pol_early)
    if "profile-fit" in getattr(pol_early, "forbidden_methods", ()):
        # The band policy has already ruled that profile fitting is undefined here, so
        # there is no EW file to read and there never will be. Take the synthesis route
        # rather than failing on a missing input that is missing BY DESIGN (RYA-832).
        return synthesis_route(a, pol_early)

    stem = f"{a.element}{a.ion}_{int(a.lo)}_{int(a.hi)}_{a.instrument}_PROFILEFIT"
    src = EW_DIR / f"{stem}_ew.csv"
    if not src.exists():
        raise SystemExit(f"no measured EWs at {src}. Run measure_band_profilefit.py first.")
    ew = pd.read_csv(src)
    pol = resolve_band(0.5 * (a.lo + a.hi))
    ok = ew[ew.in_aggregate].copy()
    print(f"{a.element} {a.ion}  {a.lo:.0f}-{a.hi:.0f} A  band={pol.name}  "
          f"{len(ok)} in-aggregate of {len(ew)} measured")

    # ── RYA-807: consume the curated problem-children registry ───────────────
    # RYA-463 built the registry, RYA-774 landed this harness, and nothing connected
    # them: Fe I 8024.543 — registered ATOMIC_BLEND / exclude / critical / active,
    # "MISIDENTIFIED, the absorber is another species" — entered the merged RYA-783
    # product at A = 9.678, ~150x solar. The docstring at the top of this file already
    # PROMISED that quarantined lines are "carried into the output with their reason and
    # excluded from the aggregate, never dropped (RYA-711)"; until now that was true only
    # of the in-harness quarantine (REW ceiling, tellurics, RYA-342 fit quality), never of
    # the curated registry.
    #
    # The exclude-vs-flag decision lives in `problem_children.aggregate_action`, NOT here,
    # so every consumer gets the same answer. Its discriminator is `status`: a row that is
    # `exclude` but `owed` has an undiagnosed cause, and dropping it would be tuning
    # (RYA-161) rather than provenance.
    from pipeline import problem_children as _pc
    _pc_table = _pc.line_dispositions()
    _pc_hits = {}
    for _w in ok.wavelength_air_A.astype(float):
        _d = _pc.disposition_for_line(a.element, a.ion, _w, table=_pc_table)
        if _d is not None:
            _pc_hits[round(_w, 4)] = (_d, _pc.aggregate_action(_d))
    _n_excl = sum(1 for _, act in _pc_hits.values() if act == "exclude")
    _n_flag = sum(1 for _, act in _pc_hits.values() if act == "flag")
    print(f"    registry: {len(_pc_hits)} of {len(ok)} pool lines are registered "
          f"({_n_excl} excluded from the aggregate, {_n_flag} flagged and KEPT)")
    for _w, (_d, _act) in sorted(_pc_hits.items()):
        if _act == "exclude":
            print(f"      EXCLUDE {_w:9.3f}  {_d['problem_class']}/{_d['status']} "
                  f"[{_d['governing_tickets']}]")

    def _stamp(lm):
        """Carry the registry verdict onto a measurement, and honour it.

        🔴 RYA-880 also makes this the choke point for NLTE provenance. Every row in this
        driver passes through here, so a route that forgets to record what it applied
        fails LOUDLY at the point of omission rather than shipping a blank column that
        nobody can tell from "LTE".
        """
        if not lm.nlte_source:
            raise SystemExit(
                f"RYA-880: {lm.element} {lm.ion} {lm.wavelength_air_A} "
                f"({lm.treatment or 'no treatment'}) has no `nlte_source`. Every row must "
                f"state what departure was applied, and an LTE row must say so explicitly "
                f"— a blank cannot be told apart from an unrecorded correction (RYA-833).")
        hit = _pc_hits.get(round(float(lm.wavelength_air_A), 4))
        if hit is None:
            return lm
        d, act = hit
        lm.problem_class = str(d.get("problem_class", ""))
        lm.problem_status = str(d.get("status", ""))
        lm.problem_tickets = str(d.get("governing_tickets", ""))
        lm.problem_action = act
        if act == "exclude":
            lm.in_aggregate = False
            # Prepend, so an in-harness reason already set is not lost.
            why = (f"REGISTRY-{d.get('problem_class', '')}: {d.get('required_treatment', '')}"
                   f"/{d.get('status', '')} per data/registry/problem_children.csv "
                   f"[{d.get('governing_tickets', '')}] — carried, not dropped (RYA-711)")
            lm.excluded_reason = why if not lm.excluded_reason else f"{why} | {lm.excluded_reason}"
        return lm
    if not len(ok):
        raise SystemExit("no lines survived measurement -- nothing to derive")

    # ── 1D-LTE ────────────────────────────────────────────────────────────────
    print("\n[1] 1D-LTE via the project's Turbospectrum curve-of-growth...")
    from scripts.control_synthesis_handler import build_context
    ctx = build_context(a.element, a.ion, 500000.0)
    from pipeline.abundances_derive import _bisect_synth_abundance
    acc = pd.read_csv(ROOT / "data" / "audit" / "line_accounting" / "per_line.csv")

    rows: list[LineMeasurement] = []
    for _, r in ok.iterrows():
        c = float(r.wavelength_air_A)
        w_nm = np.linspace(c - 0.6, c + 0.6, 300) / 10.0
        try:
            A, conv, _ = _bisect_synth_abundance(
                w_nm, float(r.ew_mA), ctx["atmosphere"], ctx["teff"], ctx["logg"],
                ctx["feh"], ctx["vturb"], ctx["linelist"], ctx["isotopes"],
                ctx["solar_abund"], a.element, ctx["atom_code"])
        except Exception as e:
            A, conv = float("nan"), False
        lm = LineMeasurement(element=a.element, ion=a.ion, wavelength_air_A=c,
                             instrument=a.instrument, ew_mA=float(r.ew_mA),
                             ew_method=str(r.ew_method), abundance=(A if conv else None),
                             # RYA-871 — carried through from the EW artifact, which now
                             # carries it from the accounting row its wavelength came
                             # from. `_intake_ep` is loud about an EW file written before
                             # RYA-871 rather than emitting a null and widening blind.
                             ep_eV=_intake_ep(r, c, a.element, a.ion),
                             nlte_delta_dex=0.0, nlte_source=LTE_NLTE_SOURCE,
                             treatment="1D-LTE",
                             # RYA-911 — carried through from the EW artifact that
                             # measured it. `_f` yields None for an EW file written
                             # before RYA-911; None means "this artifact predates the
                             # instrumentation", never "the continuum was unity".
                             continuum_level=_f(getattr(r, "continuum_level", None)),
                             continuum_method=str(getattr(r, "continuum_method", "") or ""),
                             continuum_ref=_f(getattr(r, "continuum_ref", None)),
                             # RYA-911 — same None-means-predates-instrumentation rule.
                             profile_sigma_A=_f(getattr(r, "profile_sigma_A", None)),
                             profile_sigma_floor_A=_f(
                                 getattr(r, "profile_sigma_floor_A", None)),
                             red_chi2=_f(getattr(r, "red_chi2", None)))
        if not conv or not np.isfinite(A):
            lm.in_aggregate = False
            lm.excluded_reason = ("COG-INVERSION: bisection did not converge, so this EW "
                                  "does not map to an abundance here")
        rows.append(_stamp(lm))       # RYA-807

    # RYA-869 — these EWs came out of the profile fitter (`*_PROFILEFIT_ew.csv`, written
    # by `measure_band_ew`). The COG inversion below turns them into abundances but does
    # not re-measure the spectrum, so the systematic this pool carries is still the
    # profile fitter's.
    p_lte = build_product(a.element, a.ion, a.instrument, pol.name, "1D-LTE", rows,
                          handler="ProfileFitHandler",
                          provenance=f"EWs from {src.name}; COG via _bisect_synth_abundance")
    print(f"    1D-LTE: A={p_lte.value}  n={p_lte.n_lines}  excluded={p_lte.n_excluded}")

    # ── ENGINE-A ──────────────────────────────────────────────────────────────
    print("\n[2] ENGINE-A — live MPIA per-line departure corrections...")
    used = [l for l in rows if l.in_aggregate and l.abundance is not None]
    deltas = mpia_delta(a.element, a.ion,
                        np.array([l.wavelength_air_A for l in used]),
                        cache=a.mpia_cache)
    rows_a: list[LineMeasurement] = []
    for l in used:
        # Tolerance match -- the cache keys are rounded, and two catalogues quoting the
        # same transition differ by up to ~0.03 A (RYA-704).
        d = None
        if deltas:
            k = min(deltas, key=lambda x: abs(x - l.wavelength_air_A))
            if abs(k - l.wavelength_air_A) < 0.06:
                d = deltas[k]
        la = LineMeasurement(element=l.element, ion=l.ion,
                             wavelength_air_A=l.wavelength_air_A,
                             instrument=l.instrument, ew_mA=l.ew_mA,
                             ew_method=l.ew_method, treatment="ENGINE-A",
                             # Same transition, same identity (RYA-871): ENGINE-A adds a
                             # departure term to THIS line's inversion, it does not
                             # re-measure anything.
                             ep_eV=l.ep_eV,
                             # RYA-880 — `d` IS the correction, recorded in the same
                             # expression that applies it, so the two can never drift.
                             nlte_delta_dex=(float(d) if d is not None else None),
                             nlte_source=(MPIA_NLTE_SOURCE if d is not None
                                          else "MPIA per-line delta_nlte (NOT SERVED)"),
                             # RYA-911 — ENGINE-A rides THIS line's EW, so it rides THIS
                             # line's continuum. Carrying it makes that inheritance
                             # visible in the artifact instead of something a reader has
                             # to know about the route.
                             continuum_level=l.continuum_level,
                             continuum_method=l.continuum_method,
                             continuum_ref=l.continuum_ref,
                             # RYA-911 — ENGINE-A rides THIS line's EW, so it rides the
                             # fit that produced it too.
                             profile_sigma_A=l.profile_sigma_A,
                             profile_sigma_floor_A=l.profile_sigma_floor_A,
                             red_chi2=l.red_chi2,
                             abundance=(l.abundance + d) if d is not None else None)
        if d is None:
            la.in_aggregate = False
            la.excluded_reason = ("ENGINE-A-NOT-SERVED: MPIA returns no usable delta_nlte "
                                 "for this line (absent, nan, or a placeholder zero). "
                                 "Reduced coverage, not a failed correction.")
        rows_a.append(_stamp(la))     # RYA-807
    # A departure term added to the SAME inversion — the same profile-fit EWs, so the
    # same handler and the same measured residual (RYA-869).
    p_a = build_product(a.element, a.ion, a.instrument, pol.name, "ENGINE-A", rows_a,
                        handler="ProfileFitHandler",
                        provenance="Bergemann MPIA per-line delta_nlte, live query, solar node")
    print(f"    ENGINE-A: A={p_a.value}  n={p_a.n_lines}  not-served={p_a.n_excluded}")

    # ── ENGINE-B ──────────────────────────────────────────────────────────────
    #
    # RYA-784. This block did not exist: the driver reserved Engine B and produced no
    # value, so RYA-783's Engine-B column could not resolve. The flux-fit engine itself
    # has been Done since RYA-287/338 and its adapter was fixed by RYA-770 (24/30 lines,
    # -0.026 dex against the banked optical answer); nothing had connected the two.
    #
    # It is a DIFFERENT MEASUREMENT, not a correction applied to the previous one:
    # 1D-LTE inverts an EW and ENGINE-A adds a departure term to that same inversion,
    # while this fits the observed flux. So it re-derives from the spectrum and shares
    # only the line set — which is exactly what makes the per-engine spread meaningful
    # rather than tautological (RYA-712).
    #
    # RYA-342: the gate here is FIT QUALITY, not the EW saturation ceiling. That ceiling
    # is an EW-path concept and the handler must never inherit it — flux-space synthesis
    # exists precisely to measure the lines EW filtering kills.
    # RYA-288: broadening comes from `build_context` (per-star, refuses a vmac='fit' rail),
    # so no solar default leaks in.
    p_b = None
    if a.skip_engine_b:
        print("\n[3] ENGINE-B — skipped (--skip-engine-b)")
    elif a.engine_b_deck == "gerber-nlte":
        # RYA-798 — WIRED. This branch used to refuse, and its stated reason changed twice
        # as the truth did: first "not provisioned" (wrong once RYA-785 registered the
        # deck), then "validated but not wired" (right until now). What was actually
        # missing was never the physics:
        #
        #   * iSpec's Turbospectrum wrapper has always accepted NLTE — `generate_spectrum(
        #     ..., nlte_departure_coefficients=...)` writes the departure + nlteinfo files;
        #   * but `_synth_flux_at_abund`, the ONE generator the v1 EW path and the v2 flux
        #     fit both call, had no NLTE parameter, so there was nowhere to hand them in;
        #   * and iSpec's own interpolator reads `input/dep-grid/{El}_nlte_grid_data.h5`,
        #     which is EMPTY on Sirius, while our deck is TS-native (atom.* + auxData +
        #     .bin). `pipeline/gerber_nlte.py` is the adapter between the two.
        #
        # ⚠️ iSpec fails SOFT: an element whose linelist rows carry no NLTE label lands in
        # `nlte_ignored` and the synthesis proceeds in LTE WITHOUT RAISING — the RYA-764
        # trap (2,644 in-band Fe I lines at nlte_label_up='none'). iSpec never reports
        # `nlte_available` back, so the check cannot be made after the call: it is made
        # before, by `assert_linelist_supports_nlte`, which raises instead.
        #
        # ⚠️ The departures do NOT track the fitted abundance. Measured, not assumed: the
        # interpolator run at A = 7.36 / 7.46 / 7.56 returns byte-identical files apart
        # from the abundance stamp on line 8, and the Gerber aux table has no abundance
        # axis at all (A(X) is exactly 7.50 + [Fe/H] across all 15229 nodes). So the
        # coefficients are computed ONCE per node and held fixed across the chi2 loop while
        # only the stamp follows each trial value — which is what the deck does anyway. It
        # remains an approximation of the physics, second-order over the ~0.3 dex a fit
        # explores, and it is stated on the product rather than buried in a comment.
        #
        # This is a SEPARATE PRODUCT from Engine-B-LTE, never a correction applied to it
        # (RYA-712). The Gerber-vs-MPIA spread stays an RYA-525 diagnostic.
        pass    # RYA-798 wired it; the run is below, sharing the LTE block.
    if not a.skip_engine_b:
        nlte = a.engine_b_deck == "gerber-nlte"
        treatment = "ENGINE-B-NLTE" if nlte else "ENGINE-B"
        # RYA-880. ⚠️ The NLTE deck has NO additive per-line delta: the departures enter
        # the radiative transfer and RYA-712 makes this a separate product, not a
        # corrected LTE value. None means "no additive correction exists on this route",
        # which is a different statement from 0.0, and `nlte_source` says which.
        _b_delta = None if nlte else 0.0
        _b_source = None if nlte else LTE_NLTE_SOURCE
        print(f"\n[3] {treatment} — Turbospectrum flux-fit (deck={a.engine_b_deck})...")
        from pipeline.measure import resolve_handler
        # RYA-913: the DISPATCH, never an instrument-specific loader. This block
        # imported kp_segments/load_kp_window by name and called them
        # unconditionally, so `--instrument harps` measured Kitt Peak and
        # labelled it HARPS -- the retracted 7.486 (RYA-911). An engine is a
        # function of (lines, instrument); it may not know one instrument's
        # loader exists.
        from scripts.measure_band_ew import load_window_ex, select_holding

        # ── RYA-798: the Gerber deck is MARCS, and the production default is not ────
        # NLTEgrid4TS_Fe_MARCS was computed on MARCS atmospheres. `_load_atmosphere`
        # defaults to ATLAS9.Castelli: 72 layers and only 10 columns, so index 7 — the
        # column iSpec writes into the departure file as tau — is not a tau at all and
        # reads as zeros. iSpec OVERWRITES the departure tau with that column, so ATLAS9
        # would pair MARCS departures with a degenerate depth scale and still return a
        # spectrum. MARCS.GES matches the deck exactly: 56 layers, tau -4.9154..1.7709
        # against the departure grid's -4.9177..1.7744.
        #
        # ⚠️ THIS ALSO MEANS THE LTE COMPARAND MUST BE MARCS. Comparing a MARCS Engine-B
        # NLTE against an ATLAS9 Engine-B LTE would fold an ATMOSPHERE difference into the
        # reported NLTE effect — the exact confound RYA-542 had to disentangle for Ti
        # (MARCS +0.203 reproducing deck +0.221 => ATMOSPHERE). So the atmosphere is
        # swapped for the whole Engine-B leg whenever the NLTE deck is selected, and the
        # product records which atmosphere it ran on.
        ctx_b = dict(ctx)
        if nlte:
            from pipeline.abundances_derive import _load_atmosphere
            from pipeline import gerber_nlte as gnlte
            ctx_b["atmosphere"] = _load_atmosphere(
                float(ctx["teff"]), float(ctx["logg"]), float(ctx["feh"]),
                float(ctx["vturb"]), model_grid="MARCS.GES")
            ctx_b["nlte_deck"] = "gerber"
            dep = gnlte.for_node("Fe" if a.element == "Fe" else a.element,
                                 float(ctx["teff"]), float(ctx["logg"]),
                                 float(ctx["feh"]))
            gnlte.assert_depth_match(dep, ctx_b["atmosphere"])
            n_lab = gnlte.assert_linelist_supports_nlte(
                ctx["linelist"], int(ctx["atom_code"]), a.element)
            _b_source = GERBER_NLTE_SOURCE_FMT.format(
                atom=dep['atom_path'].split('/')[-1])
            print(f"    deck atom={dep['atom_path'].split('/')[-1]} "
                  f"ndep={dep['ndep']} nk={dep['nk']} A_deck={dep['deck_abundance']}")
            print(f"    atmosphere MARCS.GES ({len(ctx_b['atmosphere'])} layers), "
                  f"{n_lab} NLTE-labelled {a.element} lines in the list")

        # TELLURIC — RYA-786. An earlier version of this block DECLARED
        # `telluric_corrected: True` from the instrument catalog to get past the handler's
        # gate. That was wrong in kind and is removed: `telluric_required = no` for a
        # reference atlas means the tellurics are handled by per-line clean-line SELECTION,
        # NOT that the atlas has been telluric-divided. The KPNO atlas HAS tellurics in it,
        # so declaring it corrected asserted something that never happened.
        #
        # The determination now lives in `pipeline.telluric_policy` and the handler
        # consumes it: the instrument decides whether a CORRECTION STAGE is required, and
        # the enumerated O2/H2O bands quarantine individual lines inside them. Passing the
        # instrument through the context is all this driver has to do.
        handler = resolve_handler(3400.0)      # the synthesis handler
        handler.prepare(pol, {**ctx_b, "instrument": a.instrument})
        # Ask once, at band centre, so a refusal costs nothing and arrives before
        # any synthesis runs (the synthesis_route pattern, RYA-904/832).
        _b_probe = select_holding(a.instrument, 0.5 * (a.lo + a.hi), 1.4)
        if not _b_probe.pre_normalised:
            raise SystemExit(
                f"{a.instrument} -> {_b_probe.holding_id} is NOT "
                f"continuum-normalised, and this route fits observed flux "
                f"against a synthesis normalised to unity. The fit would spend "
                f"A({a.element}) closing a continuum offset rather than "
                f"measuring an abundance (RYA-713/843). Normalise the product "
                f"first, or derive this band through a route that sets its own "
                f"continuum. REFUSING rather than falling back to another "
                f"atlas -- that is how the RYA-911 mislabel happened.")
        print(f"  [holding] {a.instrument} -> {_b_probe.holding_id} "
              f"(pre_normalised={_b_probe.pre_normalised})")
        _b_holdings: set[str] = set()
        segs = None
        rows_b: list[LineMeasurement] = []
        for _, r in ok.iterrows():
            c = float(r.wavelength_air_A)
            try:
                _win = load_window_ex(a.instrument, c, pad=1.4, segs=segs)
                w_obs, f_obs = _win.wave, _win.flux
                # RYA-913: record WHICH holding served each line, so a consumer
                # can verify the product was measured on what its label claims.
                _b_holdings.add(_win.holding.holding_id)
            except Exception as e:
                lb = LineMeasurement(element=a.element, ion=a.ion, wavelength_air_A=c,
                                     instrument=a.instrument, ew_mA=float("nan"),
                                     ew_method="synthesis flux-fit",
                                     ep_eV=_intake_ep(r, c, a.element, a.ion),
                                     nlte_delta_dex=_b_delta, nlte_source=_b_source,
                                     treatment=treatment)
                lb.in_aggregate = False
                lb.excluded_reason = f"WINDOW-LOAD: {type(e).__name__}: {str(e)[:60]}"
                rows_b.append(_stamp(lb))

                continue
            lb = handler.measure_line(
                w_obs, f_obs, element=a.element, ion=a.ion, wavelength_A=c,
                instrument=a.instrument, policy=pol,
                # Kitt Peak is atlas residual flux; HARPS arrives normalised by our own
                # pipeline. Both are already on a continuum and neither wants a second.
                pre_normalised=True,
                context={**ctx_b, "ew_hint_mA": float(r.ew_mA)})
            lb.treatment = treatment
            # RYA-871 — the flux fit re-measures the same TRANSITION from the same line
            # set, so it carries the same identity. The handler cannot know it: it is
            # handed a wavelength and a window, not an accounting row.
            lb.ep_eV = _intake_ep(r, c, a.element, a.ion)
            # RYA-880 — the handler builds this row and cannot know the deck, so the
            # driver that CHOSE the deck records it.
            lb.nlte_delta_dex, lb.nlte_source = _b_delta, _b_source
            rows_b.append(_stamp(lb))   # RYA-807

        # RYA-768: deterministic row order, so the artifact byte-diffs clean.
        rows_b.sort(key=lambda l: (l.wavelength_air_A, l.element, l.ion))
        # The product's declared treatment MUST match the treatment on its lines, or
        # assert_no_cross_treatment_mix refuses it (RYA-712). Declaring "ENGINE-B" while
        # the lines carry ENGINE-B-NLTE is not a cosmetic mismatch: it is a mislabelled
        # engine, which is the exact failure that guard exists to catch.
        prov_deck = (
            "the TS-native Gerber departure deck (atom.fe607a, validated RYA-785, wired "
            "RYA-798), on MARCS.GES to match the atmosphere family the departures were "
            "computed on; departures are node-fixed because the deck has no abundance "
            "axis, so only the bsyn abundance stamp follows each trial value"
            if nlte else
            "Turbospectrum LTE, NOT the Gerber TS-native NLTE deck")
        # 🔴 RYA-913 — PROVENANCE MUST MATCH INVOCATION. The defect this ticket exists
        # for produced a product tagged `harps` whose every line was read off the Kitt
        # Peak atlas, and nothing objected: the tag came from the CLI flag and the data
        # came from a hardcoded loader, so the two could disagree indefinitely. Asserting
        # the tag against the holdings that actually served the lines is the only check
        # that would have failed. It is cheap and it runs on every product.
        _b_expected = {h.holding_id for h in
                       [select_holding(a.instrument, 0.5 * (a.lo + a.hi), 1.4)]}
        _b_unexpected = _b_holdings - _b_expected
        if _b_unexpected:
            raise SystemExit(
                f"PROVENANCE MISMATCH: this product is tagged instrument "
                f"{a.instrument!r} but its lines were measured on holding(s) "
                f"{sorted(_b_unexpected)}, which {a.instrument!r} does not select. "
                f"A number labelled with an instrument it was not measured on is the "
                f"RYA-911/913 defect. Refusing to emit it.")
        if _b_holdings:
            print(f"  [provenance] {len(_b_holdings)} holding(s) served this product: "
                  f"{sorted(_b_holdings)} — matches --instrument {a.instrument}")
        p_b = build_product(
            a.element, a.ion, a.instrument, pol.name, treatment, rows_b,
            # RYA-869 — read off the handler OBJECT that actually ran, not asserted.
            # `treatment` here is ENGINE-B or ENGINE-B-NLTE depending only on the deck,
            # and BOTH are this one handler; the budget loop below used to tell them
            # apart by name and charge the NLTE one the profile fitter's residual.
            handler=handler.name,
            provenance=(f"{handler.name} flux-fit against the {a.instrument} spectrum; "
                        f"line set + window hint from {src.name}; deck={a.engine_b_deck} "
                        f"({prov_deck})"))
        print(f"    {treatment}: A={p_b.value}  n={p_b.n_lines}  excluded={p_b.n_excluded}")
        if p_b.n_excluded:
            why = pd.Series([l.excluded_reason.split(":")[0]
                             for l in rows_b if not l.in_aggregate]).value_counts()
            for k, v in why.items():
                print(f"      {k:<30} {v}")

    # ── budgets, per product ──────────────────────────────────────────────────
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    summary, budgets = [], {}
    for prod in (p_lte, p_a, p_b):
        if prod is None or prod.value is None:
            continue
        # 🔴 RYA-869 — the harness residual is a property of the HANDLER that produced
        # the number, and it is READ FROM THE PRODUCT now instead of inferred here.
        # What stood here was `is_b = prod.treatment == "ENGINE-B"`, under a comment
        # stating the handler rule correctly. The comment was right; the code was an
        # equality against a treatment NAME, and that name gained the variant
        # `ENGINE-B-NLTE` in RYA-798 — the same SynthesisHandler flux fit on a different
        # departure deck. Every NLTE product was therefore charged the profile fitter's
        # 0.0129 dex AND labelled `ProfileFitHandler` in its own budget file, beside the
        # LTE product of the same handler charged 0.0000 and labelled correctly. Four
        # published Fe matrix bars carried it (RYA-783/807/832).
        # A WIDER SET WOULD NOT HAVE FIXED IT: the near-UV route above wears the
        # treatment `1D-LTE` and is a flux fit, so no mapping from the label exists.
        harness = harness_residual.for_product(prod)
        print(f"    [harness] {prod.treatment}: {harness.describe()}")
        # 🔴 RYA-855 — PER PRODUCT, not per band and not per element. This call passed
        # `gf_graded=False` unconditionally, the mirror of the synthesis route above, so
        # every EW-route cell charged the ungraded 0.17 no matter what its lines were.
        # It matters that this is per PRODUCT: 1D-LTE, ENGINE-A and ENGINE-B share a
        # band and an element but NOT a line set — Engine A drops every line MPIA does
        # not serve — so the three can legitimately land on different rungs, and a
        # band-level or element-level answer would be wrong for at least one of them.
        rung = gf_rung.for_product(prod, linelist=ctx["linelist"])
        print(f"    [gf rung] {prod.treatment}: {rung.describe()}")
        # 🔴 RYA-907 — the SECOND home of the identical idiom. Fixing only the synthesis
        # route above would have left the bug live on the EW route, which is exactly
        # where both offending Fe II cells were produced. Same shape as RYA-855 and
        # RYA-845: one wrong expression, two call sites, and only one of them read.
        b = build_budget(a.element, 0.5 * (a.lo + a.hi), prod.n_lines,
                         scatter_dex=prod.sigma, **rung.budget_kwargs(),
                         **harness.budget_kwargs())
        stat, syst = b.total()
        assert_stat_publishable(
            stat, cell=f"{a.element} {a.ion} {pol.name} {prod.treatment} "
                       f"({a.instrument}, n={prod.n_lines})")
        summary.append(dict(element=a.element, ion=a.ion, band=pol.name,
                            instrument=a.instrument, treatment=prod.treatment,
                            # RYA-869 — the artifact RECORDS which handler produced the
                            # cell. The budget file named it in prose and the products
                            # table did not carry it at all, so a downstream consumer
                            # (the RYA-783 matrix, RYA-850's recharge) had to guess it
                            # back from the treatment — which is the defect, one layer out.
                            handler=prod.handler,
                            A=round(prod.value, 3), n_lines=prod.n_lines,
                            n_excluded=prod.n_excluded,
                            stat_dex=round(stat, 4), syst_dex=round(syst, 4),
                            stat_basis=b.stat_basis(),   # RYA-907
                            dominant=b.dominant().name if b.dominant() else "",
                            # RYA-906 — the axes, ADDED ALONGSIDE `treatment`, which is never
                            # rewritten. `route` comes from the HANDLER, never the label:
                            # `1D-LTE` names both an EW inversion and a synthesis flux
                            # fit, so on that label the legacy string is FALSE, not just
                            # lossy. This is the same hand-written projection that lost
                            # `red_chi2` in RYA-847 and had to be patched again in RYA-880:
                            # it is the lossy one of the object's two serialisers, so a new
                            # column must be added HERE as well as in `products_frame`.
                            **axes_for(prod.treatment,
                                       handler=prod.handler or None).as_columns(),
                            ))
        budgets[prod.treatment] = b.describe() + f"\n  gf rung: {rung.describe()}"
        prod.to_frame().to_csv(out / f"{stem}_{prod.treatment}_lines.csv", index=False)

    df = pd.DataFrame(summary)
    df.to_csv(out / f"{stem}_products.csv", index=False)
    (out / f"{stem}_budgets.txt").write_text("\n\n".join(budgets.values()))

    print("\n=== PRODUCTS (separate, never combined — RYA-712) ===")
    print(df.to_string(index=False))
    print()
    for t, d in budgets.items():
        print(d); print()
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()
