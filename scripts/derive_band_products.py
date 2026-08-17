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
from pipeline.band_products import build_product, LineMeasurement, products_frame  # noqa: E402
from pipeline.error_budget import build as build_budget  # noqa: E402
from pipeline.fit_constraint import as_float_or_none as _f  # noqa: E402  RYA-847
from pipeline.constraint_gate import verdict as constraint_verdict  # noqa: E402
from pipeline.constraint_gate import describe as constraint_describe  # noqa: E402

EW_DIR = ROOT / "data" / "measured" / "band_ew"
OUT = ROOT / "data" / "results" / "band_products"

# Measured optical residual of the profile fitter against the banked HARPS pool (RYA-713).
# Carried as the harness systematic in every frontier band rather than assumed zero.
PROFILE_FIT_RESIDUAL_DEX = 0.0129


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
    from rya759_nearuv_fe_product import select_lines, fit_one
    from rya759_nearuv_synth import _kp_segments

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
    cand = select_lines(ctx["linelist"], lo_A=a.lo, hi_A=a.hi, n=cfg.n_lines,
                        teff=float(ctx["teff"]), min_sep_A=cfg.min_sep_A)
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
    print(f"  {len(cand)} Fe I candidates by theoretical depth "
          f"(half-width +/-{hw} A, min separation {cfg.min_sep_A} A)")
    print(f"  [half-width] {cfg.half_width_note}")

    tmp = f"/tmp/ispec_synth_{pol.name.replace(chr(47), chr(45))}"
    Path(tmp).mkdir(parents=True, exist_ok=True)
    lines: list[LineMeasurement] = []
    for r in cand.itertuples():
        w = float(r.wave_A)
        res = fit_one(ctx, segs, w, hw, tmp)
        a_x = float(res.get("a_synth", float("nan")))
        lm = LineMeasurement(
            element=a.element, ion=a.ion, wavelength_air_A=w,
            instrument=a.instrument, ew_mA=float("nan"),
            ew_method=(f"synthesis flux-fit, FIXED half-width +/-{hw} A "
                       f"(RYA-759 route; no EW exists in this band to key the "
                       f"production wing-wide rule)"),
            abundance=(a_x if np.isfinite(a_x) else None),
            treatment="1D-LTE",
            # The REW saturation ceiling is an EW-INVERSION concept and there is no EW
            # here at all. Inheriting it would quarantine lines on a quantity that does
            # not exist (RYA-770/342).
            ew_inversion=False,
            # RYA-847 — carried, not yet gated on. The sweep chooses the metric and the
            # threshold across every synthesis band; this route must not pick either
            # from its own product (RYA-161).
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
            # It now calls the SAME decider the Engine-B handler calls. While
            # `SYNTH_CONSTRAINT` is None this changes nothing numerically — deliberately,
            # and provably: the near-UV product reproduces 7.488 with the call in place.
            # The cut comes from a cross-band sweep, not from this band (RYA-161).
            _cv = constraint_verdict(res)
            if not _cv.ok:
                lm.in_aggregate = False
                lm.excluded_reason = _cv.reason
        lines.append(lm)
    lines.sort(key=lambda l: (l.wavelength_air_A, l.element, l.ion))

    prov = (
        "Fe I near-UV 1D-LTE by flux-fit synthesis (Turbospectrum via iSpec, "
        "ATLAS9.Castelli) — the RYA-759 validated route called directly, same functions "
        "and same defaults, so the value cannot drift from that ticket's. "
        "SYNTHESIS-ONLY: band_policy forbids profile-fit and interval-integration here "
        "(median line gap 0.146 A leaves no interval that contains one profile and "
        "excludes its neighbours), and RYA-759 falsified profile fitting in band — 901 "
        "candidates, 0 measurable. 1D-LTE ONLY: UV Fe I is heavily over-ionised, so the "
        "missing NLTE correction is large and POSITIVE, and a low value here is expected "
        "physics rather than a defect. gf: " + str(prov_gf["detail"]) + ". "
        "PSEUDO-CONTINUUM SYSTEMATIC 0.100 dex, which does NOT average down and is NOT "
        "in the scatter reported here. Half-width is FIXED and must be swept. "
        + constraint_describe() + " " +
        "gf REMAINS THE DOMINANT SYSTEMATIC AND THE BAND IS UNGRADED: RYA-822 grades "
        "only 6 of the 4,274 in-band Fe I lines as primary-lab, and its GF-NIST class "
        "(604 lines) is a COMPILATION grade that 822 deliberately keeps outside "
        "`is_graded` because FMW *is* NIST and VALD copies it (RYA-760). Of the lines "
        "actually fitted here, 17 of 40 carry a citable NIST accuracy class and most of "
        "those are poor (C/C+/D/D+/E). Never coadded with another band (RYA-712).")
    product = build_product(a.element, a.ion, a.instrument, pol.name, "1D-LTE",
                            lines, provenance=prov)

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
                     scatter_dex=(product.sigma or 0.0), gf_graded=False,
                     # There is no profile fitter anywhere on this route, so charging it
                     # the profile fitter's measured residual would attribute someone
                     # else's systematic to it.
                     harness_residual_dex=0.0, handler="SynthesisHandler")
    stat, syst = b.total()
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
        treatment="1D-LTE", A=round(product.value, 3) if product.value is not None else None,
        n_lines=product.n_lines, n_excluded=product.n_excluded,
        stat_dex=round(stat, 4), syst_dex=round(syst, 4),
        dominant=(b.dominant().name if b.dominant() else ""),
    )]).to_csv(out / f"{stem}_products.csv", index=False)
    # `describe()` already lists every term the budget holds, pseudo-continuum included.
    # The epilogue that used to be appended here announced the term as "added in
    # quadrature above" — which is exactly the double-add RYA-845 removed, and it made
    # the artifact assert the bug in prose as well as in arithmetic.
    (out / f"{stem}_budgets.txt").write_text(b.describe() + "\n")
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


def asdict_line(l: LineMeasurement) -> dict:
    return dict(element=l.element, ion=l.ion, wavelength_air_A=l.wavelength_air_A,
                instrument=l.instrument, ew_mA=l.ew_mA, ew_method=l.ew_method,
                abundance=l.abundance, rew=l.rew, treatment=l.treatment,
                in_aggregate=l.in_aggregate, excluded_reason=l.excluded_reason,
                ew_inversion=l.ew_inversion,
                # RYA-847: the constraint metrics. None on EW-route lines — no chi2
                # surface exists behind an EW inversion, so the question does not apply,
                # and a consumer must read None as "not applicable", not "unconstrained".
                sigma_A=l.sigma_A, frac_rise_weaker=l.frac_rise_weaker,
                edge_distance_dex=l.edge_distance_dex, red_chi2=l.red_chi2)


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
        """Carry the registry verdict onto a measurement, and honour it."""
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
                             treatment="1D-LTE")
        if not conv or not np.isfinite(A):
            lm.in_aggregate = False
            lm.excluded_reason = ("COG-INVERSION: bisection did not converge, so this EW "
                                  "does not map to an abundance here")
        rows.append(_stamp(lm))       # RYA-807

    p_lte = build_product(a.element, a.ion, a.instrument, pol.name, "1D-LTE", rows,
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
                             abundance=(l.abundance + d) if d is not None else None)
        if d is None:
            la.in_aggregate = False
            la.excluded_reason = ("ENGINE-A-NOT-SERVED: MPIA returns no usable delta_nlte "
                                 "for this line (absent, nan, or a placeholder zero). "
                                 "Reduced coverage, not a failed correction.")
        rows_a.append(_stamp(la))     # RYA-807
    p_a = build_product(a.element, a.ion, a.instrument, pol.name, "ENGINE-A", rows_a,
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
        print(f"\n[3] {treatment} — Turbospectrum flux-fit (deck={a.engine_b_deck})...")
        from pipeline.measure import resolve_handler
        from scripts.measure_band_ew import kp_segments, load_kp_window

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
        segs = kp_segments()
        rows_b: list[LineMeasurement] = []
        for _, r in ok.iterrows():
            c = float(r.wavelength_air_A)
            try:
                w_obs, f_obs, _ = load_kp_window(segs, c, pad=1.4)
            except Exception as e:
                lb = LineMeasurement(element=a.element, ion=a.ion, wavelength_air_A=c,
                                     instrument=a.instrument, ew_mA=float("nan"),
                                     ew_method="synthesis flux-fit",
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
        p_b = build_product(
            a.element, a.ion, a.instrument, pol.name, treatment, rows_b,
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
        # The harness residual is a property of the HANDLER that produced the number.
        # Engine B never touches the profile fitter, so charging it the profile fitter's
        # measured residual would be attributing someone else's systematic to it.
        is_b = prod.treatment == "ENGINE-B"
        b = build_budget(a.element, 0.5 * (a.lo + a.hi), prod.n_lines,
                         scatter_dex=(prod.sigma or 0.0), gf_graded=False,
                         harness_residual_dex=(0.0 if is_b else PROFILE_FIT_RESIDUAL_DEX),
                         handler=("SynthesisHandler" if is_b else "ProfileFitHandler"))
        stat, syst = b.total()
        summary.append(dict(element=a.element, ion=a.ion, band=pol.name,
                            instrument=a.instrument, treatment=prod.treatment,
                            A=round(prod.value, 3), n_lines=prod.n_lines,
                            n_excluded=prod.n_excluded,
                            stat_dex=round(stat, 4), syst_dex=round(syst, 4),
                            dominant=b.dominant().name if b.dominant() else ""))
        budgets[prod.treatment] = b.describe()
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
