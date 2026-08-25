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
from dataclasses import replace
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.band_policy import resolve as resolve_band  # noqa: E402
from pipeline.treatment_axes import axes_for  # RYA-906
from pipeline import treatment_axes as taxes  # RYA-1040: <3D> tokens, never typed
from pipeline.band_products import build_product, LineMeasurement, products_frame  # noqa: E402
from pipeline.error_budget import build as build_budget  # noqa: E402
from pipeline.error_budget import assert_stat_publishable  # noqa: E402  RYA-907
from pipeline import gf_rung  # noqa: E402  RYA-855
from pipeline import harness_residual  # noqa: E402  RYA-869
from pipeline.fit_constraint import as_float_or_none as _f  # noqa: E402  RYA-847
from pipeline.constraint_gate import verdict as constraint_verdict  # noqa: E402
from pipeline.constraint_gate import describe as constraint_describe  # noqa: E402
from pipeline.prenormalised_guard import assert_not_renormalising  # noqa: E402  RYA-1026

EW_DIR = ROOT / "data" / "measured" / "band_ew"
OUT = ROOT / "data" / "results" / "band_products"

#: RYA-869 — re-exported, NOT restated. The measured optical residual of the profile
#: fitter against the banked HARPS pool (RYA-713) was written out here AND at
#: `rya850_graded_products.py:73` — two homes for one number, the RYA-845 shape — while
#: `scripts/rya855_rung_audit.py` imports it from this module. It is now declared once in
#: `pipeline.harness_residual`, keyed by the HANDLER that earned it rather than by the
#: script that happened to need it first; this name is kept so that import still resolves.
PROFILE_FIT_RESIDUAL_DEX = harness_residual.HANDLER_RESIDUAL_DEX["ProfileFitHandler"]


def engine_a_delta(element: str, ion: str, waves: np.ndarray,
                   cache: str | None = None, star: str = "solar") -> dict[float, float]:
    """Per-line Engine-A departure corrections. Missing = not served.

    Fe is served by the live MPIA query.  Every other registered element is served by
    its vendored per-line grid through the same resolver used by the production
    abundance path.  The old implementation returned ``{}`` for every non-Fe element,
    so an Al run labelled its registered Amarsi grid "not served" without consulting it.

    The MPIA scraper needs `bs4`, which is absent from the Sirius interpreter that has
    iSpec. Rather than installing into a pinned environment mid-run, the Fe query can be
    run where bs4 lives and handed over as a cache -- the numbers are identical, and the
    file records which lines MPIA actually served.
    """
    if cache:
        d = json.loads(Path(cache).read_text())
        return {float(k): float(v) for k, v in d.items()}
    sys.path.insert(0, str(ROOT))
    from scripts.build_nlte_grids_mpia import _submit_batch
    from config.constants import get_star_params
    p = get_star_params(star)          # RYA-985: loud-fails on an unknown star
    is_sun = star == "solar"
    # 🔴 THE GRID NODE MUST BE THIS STAR'S NODE. `feh` was pinned to 0.0 here, so every
    # request fetched the SOLAR departure coefficients whatever star was being measured.
    # For tau Ceti at [Fe/H] = -0.49 that is not a rounding difference — it is the NLTE
    # correction of a different star, arriving as a plausible number. Mirrors the is_sun
    # pattern already in `pipeline/nlte_cno.py`.
    feh = 0.0 if is_sun else float(p["feh_ref"])
    node = [dict(name=("sun" if is_sun else star),
                 teff_K=int(round(float(p["teff"]))),
                 logg=float(p["logg"]), feh=feh)]
    if element != "Fe":
        from pipeline.nlte_corrections import _mpia_element_delta
        from pipeline.species import parse_ion
        from config.constants import NLTE_CORRECTION_ELEMENTS

        spec = NLTE_CORRECTION_ELEMENTS.get(element)
        if spec is None or parse_ion(ion) != int(spec["ion"]):
            return {}
        out = {}
        for wave in waves:
            delta = _mpia_element_delta(
                element, float(wave), float(p["teff"]), float(p["logg"]), feh,
                tol=0.06)
            if delta is not None and np.isfinite(delta):
                out[float(wave)] = float(delta)
        return out

    code = {"I": "26.01", "II": "26.02"}.get(ion)
    if code is None:
        return {}
    out: dict[float, float] = {}
    for i in range(0, len(waves), 40):          # MPIA dislikes very long batches
        chunk = [float(w) for w in waves[i:i + 40]]
        try:
            r = _submit_batch(chunk, code, node).get(node[0]["name"], {})
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

from config.synth_bands import SynthBand, SYNTH_BANDS  # noqa: E402,F401

#: How far a synthesis list may fall short of the requested window before the run is
#: refused. NOT a tuned number: a line list's edge is where its last line sits, so the
#: first/last line of a band can legitimately be an Angstrom inside the nominal boundary.
#: It is small enough that the 420 A GES-v6 shortfall at the VIS blue edge cannot pass.
_LIST_COVERAGE_TOL_A = 5.0

#: RYA-759's own defaults, named here so a reader can see they are NOT re-chosen.
#: Changing any of these would move the product away from the published 7.487 and the
#: ticket forbids that without a stated cause — so they are imported-in-spirit, and the
#: functions below are imported literally.
#: 🔴 RYA-967 — DERIVED FROM THE CONFIG, NOT RE-TYPED. These were literals here AND
#: values in the near-UV `SynthBand`; once the bands moved to `config/synth_bands.yaml`
#: the literals became a second declaration of the same three numbers, free to drift from
#: the entry they describe. `test_nearuv_first_class_rya832` still pins them to
#: 0.40 / 4.0 / 40, so if the YAML is edited that test fails — which is the correct
#: alarm, and the one a duplicated literal would have silenced.
_NEARUV = SYNTH_BANDS["near-UV"]
NEARUV_HALF_WIDTH_A = _NEARUV.half_width_A   # no EW exists to key the wing-wide rule
NEARUV_MIN_SEP_A = _NEARUV.min_sep_A         # keeps the set spread, not piled in one complex
NEARUV_N_LINES = _NEARUV.n_lines
#: 🔴 THE PSEUDO-CONTINUUM SYSTEMATIC IS DELIBERATELY NOT DECLARED HERE (RYA-845).
#: It lives in exactly one place, `pipeline/error_budget.py`, which adds it for any band
#: whose policy says "pseudo-continuum". A second declaration next to a route is what let
#: the term be applied twice: the route could not see what the budget already held, and
#: the number agreed with itself in both places while the product was wrong.


#: RYA-967 — `SynthBand` and `SYNTH_BANDS` MOVED TO `config/synth_bands.py`.
#:
#: They were declared here, in an executable driver, and `scripts/rya855_rung_audit.py`
#: imported them FROM THIS SCRIPT. That is a config constant whose only home is a program
#: — it cannot be read without importing this module, and this module's import chain loads
#: the Kitt Peak atlas. Same second-home shape as RYA-350/353/954.
#:
#: The four entries and every word of their reasoning live in `config/synth_bands.yaml`,
#: which also records the invariant the three original values turned out to share and
#: which the new VIS entry was derived from: a half-width is 20.51 Doppler sigma at the
#: band's anchor wavelength, held to +/-0.30 across a 3.5x wavelength lever.


#: Wavelength key for matching an EW artifact's lines into the synthesis list. RYA-945
#: MEASURED the safe window at 5 mA by binning graded NIST matches by residual and asking
#: how often the sources then disagree; below that the pairing is real, above it chance
#: pairings pile up. Both files describe the SAME transitions here, so the match should be
#: exact — the tolerance exists to absorb rounding between two writers, not to search.
_EW_MATCH_TOL_A = 0.005


def _match_into_list(want: np.ndarray, w_sorted: np.ndarray,
                     idx_sorted: np.ndarray) -> tuple[list, list]:
    """Wavelengths -> row indices in the synthesis list. Shared by every selector.

    Nearest-within-tolerance, never a rounded key (RYA-703/704). Factored out rather than
    copied into the second selector: RYA-701 measured that one Ba->Al copy produced 13
    defects, and a matcher that drifts between two selectors would silently change which
    lines a comparison is run on.
    """
    keep, missing = [], []
    for w in want:
        j = int(np.searchsorted(w_sorted, w))
        best, bi = np.inf, -1
        for k in (j - 1, j, j + 1):
            if 0 <= k < len(w_sorted) and abs(w_sorted[k] - w) < best:
                best, bi = abs(w_sorted[k] - w), k
        if bi >= 0 and best <= _EW_MATCH_TOL_A:
            keep.append(idx_sorted[bi])
        else:
            missing.append((float(w), float(best)))
    return keep, missing


def _feature_depth(waves: np.ndarray) -> np.ndarray:
    """Each line's FEATURE depth, grouped exactly as `line_accounting_rya709.features()`.

    🔴 THE GROUPING IS THE POINT, NOT AN IMPLEMENTATION DETAIL. That script groups list
    rows within `GROUP_A` and takes the MAX central_depth over the group, so "this line's
    depth" is the depth of the blended FEATURE it belongs to. Recomputing it any other way
    here would select a different population from the one the EW route gated on, and the
    whole claim of this run — that these are the lines EW could never reach — rests on the
    two agreeing. Thresholds are imported, never re-typed.
    """
    from line_accounting_rya709 import GROUP_A, DEPTH_HI            # noqa: F401
    ls = pd.read_csv(ROOT / "data" / "linelists" / "linelist_solar.csv", low_memory=False)
    a = ls[ls.element == "Fe"].sort_values("wavelength_air_A").copy()
    a["_k"] = (a.wavelength_air_A.diff().fillna(9e9) > GROUP_A).cumsum()
    f = a.groupby("_k").agg(w=("wavelength_air_A", "mean"),
                            d=("central_depth", "max")).reset_index(drop=True)
    fw, fd = f.w.values, f.d.values
    return np.array([fd[int(np.argmin(np.abs(fw - w)))] for w in waves])


def _cand_deep_graded(linelist, *, lo_A: float, hi_A: float, species: str) -> pd.DataFrame:
    """The graded lines EW can NEVER reach: laboratory gf AND too deep to measure — RYA-984.

    🔴 WHY A THIRD SELECTOR. `--lines-from-ew` reads the EW artifact, which is written
    AFTER `line_accounting_rya709`'s [0.05, 0.60] depth triage, so the deep population is
    invisible to it BY CONSTRUCTION (measured on RYA-967: the 55 it returned top out at
    depth 0.595). `select_lines` would pick its own strongest-N and vary selection as well
    as method (RYA-842). Neither can express "every graded line above the depth gate".

    These are the anchor-grower. EW dies on them because the curve of growth goes flat;
    synthesis fits flux-space chi2 and inverts no equivalent width, so the REW saturation
    ceiling never applies to any of them (`pipeline/measure/synthesis.py:277`).
    """
    from line_accounting_rya709 import DEPTH_HI
    cg = pd.read_csv(ROOT / "data" / "linelists" / "canonical_gf.csv", low_memory=False)
    lab = cg[(cg.species == species.replace(" 1", " I").replace(" 2", " II"))
             & cg.gf_tier.astype(str).str.contains("LAB", na=False)
             & cg.wavelength_air_A.between(lo_A, hi_A)]
    if lab.empty:
        raise SystemExit(
            f"no LAB-tier {species} lines in {lo_A}-{hi_A} A of canonical_gf — refusing "
            f"to run a 'graded' product on a pool that is not graded.")
    depth = _feature_depth(lab.wavelength_air_A.values.astype(float))
    deep = lab[depth > DEPTH_HI]
    print(f"  [deep-graded] {len(lab)} LAB-tier {species} lines in band; "
          f"{len(deep)} above the {DEPTH_HI} depth gate that EW could never attempt "
          f"({len(lab) - len(deep)} are in or below the EW window and belong to RYA-967)")
    if deep.empty:
        raise SystemExit("no graded line in this band sits above the EW depth gate")

    names = linelist.dtype.names
    w_A = np.asarray(linelist["wave_A"] if "wave_A" in names
                     else linelist["wave_nm"] * 10.0, dtype=float)
    el = np.asarray([str(x).strip() for x in linelist["element"]])
    idx_all = np.flatnonzero(el == species)
    if not idx_all.size:
        raise SystemExit(f"no {species!r} rows in the synthesis list")
    order = np.argsort(w_A[idx_all])
    idx_sorted, w_sorted = idx_all[order], w_A[idx_all][order]

    keep, missing = _match_into_list(
        np.sort(deep.wavelength_air_A.values.astype(float)), w_sorted, idx_sorted)
    if missing:
        # LOUD (RYA-711/833). A graded line the synthesis list does not carry is not
        # "not recovered" — it was never a candidate, and the distinction is the whole
        # reason RYA-977 exists.
        print(f"  [deep-graded] {len(missing)} of {len(deep)} are NOT in the synthesis "
              f"list (>{_EW_MATCH_TOL_A} A from any row) — never candidates, reported "
              f"rather than counted as failures (RYA-977 territory)")
    if not keep:
        raise SystemExit("no deep graded line matched the synthesis list")
    keep = np.array(sorted(set(keep)))
    return pd.DataFrame({
        "wave_A": w_A[keep],
        "loggf": np.asarray(linelist["loggf"], dtype=float)[keep],
        "ep_eV": np.asarray(linelist["lower_state_eV"], dtype=float)[keep],
        "theo_depth": np.asarray(linelist["theoretical_depth"], dtype=float)[keep],
    }).sort_values("wave_A").reset_index(drop=True)


def _cand_graded(linelist, *, lo_A: float, hi_A: float, species: str,
                 deep: bool = False) -> pd.DataFrame:
    """Every LAB-tier line in the band, split on the EW depth gate — RYA-933.

    🔴 WHY THIS EXISTS. `--lines-tier graded` was consulted ONLY inside the
    `--lines-from-ew` branch. Passing it alone fell through to `select_lines`, which
    picks the strongest N by theoretical depth and applies NO gf filter — so the run
    asked for graded, measured an UNGRADED mixed pool, and the stem did not even record
    that graded had been requested. Measured on Kurucz 2005 VIS: 2 of 37 lines GF-LAB,
    rung 1, syst 0.170, reported as a graded run. That is the RYA-833 shape (a silent
    no-op reading as a result) applied to the one axis RYA-161 says dominates.

    `deep=False` returns the graded lines AT OR BELOW the gate — graded without the deep
    population, which is a separate and much slower product (RYA-984).
    """
    from line_accounting_rya709 import DEPTH_HI
    cg = pd.read_csv(ROOT / "data" / "linelists" / "canonical_gf.csv", low_memory=False)
    lab = cg[(cg.species == species.replace(" 1", " I").replace(" 2", " II"))
             & cg.gf_tier.astype(str).str.contains("LAB", na=False)
             & cg.wavelength_air_A.between(lo_A, hi_A)]
    if lab.empty:
        raise SystemExit(
            f"no LAB-tier {species} lines in {lo_A}-{hi_A} A of canonical_gf — refusing "
            f"to run a 'graded' product on a pool that is not graded.")
    depth = _feature_depth(lab.wavelength_air_A.values.astype(float))
    sel = lab[depth > DEPTH_HI] if deep else lab[depth <= DEPTH_HI]
    print(f"  [graded] {len(lab)} LAB-tier {species} lines in band; using the "
          f"{len(sel)} {'ABOVE' if deep else 'AT OR BELOW'} the {DEPTH_HI} depth gate")
    if len(sel) < 2:
        # 🔴 RYA-1031 — A ONE-LINE POOL IS NOT A MEASUREMENT, and it reached the tracker.
        # Near-UV Fe I: 58 of the 59 LAB-tier lines sit ABOVE the 0.6 depth gate, so the
        # shallow selector returned exactly ONE line (3026.056 A) and the band's "graded"
        # value was that single fit -- A = 8.529, about 1.08 dex above every VIS arm,
        # in_aggregate with no exclusion reason and no EW or REW recorded. Nothing caught
        # it because nothing was wrong per-line; the POOL was the defect.
        #
        # The floor is 2 and it is not a threshold anyone chose: line-to-line scatter has
        # n-1 degrees of freedom, so at n=1 the statistical term cannot be computed from
        # the data at all. The 0.01 dex it reported was the RYA-771 quantiser floor
        # standing in for a number that does not exist -- precision manufactured by
        # arithmetic rather than measured (RYA-968).
        #
        # The DEPTH GATE IS THE REAL STORY and the message says so: 0.6 is calibrated on
        # VIS, where the split is meaningful (55 shallow / 108 deep). In near-UV the
        # population inverts and almost everything is saturated, so "graded" and
        # "deep-graded" stop being two views of one band and the deep pool IS the band's
        # graded pool.
        n_above = int((depth > DEPTH_HI).sum())
        n_below = int((depth <= DEPTH_HI).sum())
        raise SystemExit(
            f"--lines-tier graded selected {len(sel)} {species} line(s) in "
            f"{lo_A:.0f}-{hi_A:.0f} A -- refusing: a pool of fewer than 2 lines cannot "
            f"carry a line-to-line scatter, so its statistical term would be invented "
            f"rather than measured.\n"
            f"  LAB-tier lines in band : {len(lab)}\n"
            f"  at/below the {DEPTH_HI} gate : {n_below}   (what 'graded' selects)\n"
            f"  above the {DEPTH_HI} gate    : {n_above}   (what '--lines-deep-graded' selects)\n"
            f"  The depth gate is calibrated on VIS. Where the population is mostly "
            f"saturated the deep pool IS this band's graded pool -- use "
            f"--lines-deep-graded, and say which selection the product used.")
    names = linelist.dtype.names
    w_A = np.asarray(linelist["wave_A"] if "wave_A" in names
                     else linelist["wave_nm"] * 10.0, dtype=float)
    el = np.asarray([str(x).strip() for x in linelist["element"]])
    idx_all = np.flatnonzero(el == species)
    if not idx_all.size:
        raise SystemExit(f"no {species!r} rows in the synthesis list")
    order = np.argsort(w_A[idx_all])
    idx_sorted, w_sorted = idx_all[order], w_A[idx_all][order]
    keep, missing = _match_into_list(
        np.sort(sel.wavelength_air_A.values.astype(float)), w_sorted, idx_sorted)
    if missing:
        print(f"  [graded] {len(missing)} of {len(sel)} are NOT in the synthesis list "
              f"(>{_EW_MATCH_TOL_A} A from any row) — never candidates, reported rather "
              f"than counted as failures (RYA-977 territory)")
    if not keep:
        raise SystemExit("no graded line matched the synthesis list")
    keep = np.array(sorted(set(keep)))
    return pd.DataFrame({
        "wave_A": w_A[keep],
        "loggf": np.asarray(linelist["loggf"], dtype=float)[keep],
        "ep_eV": np.asarray(linelist["lower_state_eV"], dtype=float)[keep],
        "theo_depth": np.asarray(linelist["theoretical_depth"], dtype=float)[keep],
    }).sort_values("wave_A").reset_index(drop=True)


def _cand_from_ew_artifact(linelist, ew_csv: Path, *, lo_A: float, hi_A: float,
                           species: str, tier: str) -> pd.DataFrame:
    """The lines an EW artifact ATTEMPTED, as synthesis candidates — RYA-967.

    🔴 WHY THIS EXISTS. `select_lines` returns the strongest N lines the synthesis list
    holds. Running that against RYA-959's EW pool would compare two different LINE SETS as
    well as two methods, and RYA-842 measured that line selection dominates gf — so the
    comparison would answer the wrong question. The controlled experiment needs the same
    lines, and only the method varying.

    Every row the EW file attempted is a candidate, INCLUDING the ones it quarantined:
    those are the whole point. "Did synthesis recover the lines EW dropped, and which gate
    had killed them?" cannot be answered from the survivors.
    """
    if not ew_csv.exists():
        raise SystemExit(f"--lines-from-ew: no EW artifact at {ew_csv}")
    ew = pd.read_csv(ew_csv)
    want = ew.wavelength_air_A.astype(float)
    want = want[(want >= lo_A) & (want <= hi_A)]
    if want.empty:
        raise SystemExit(
            f"--lines-from-ew: {ew_csv.name} has no line in {lo_A}-{hi_A} A. Refusing to "
            f"synthesise a set selected from somewhere else.")

    names = linelist.dtype.names
    w_A = np.asarray(linelist["wave_A"] if "wave_A" in names
                     else linelist["wave_nm"] * 10.0, dtype=float)
    el = np.asarray([str(x).strip() for x in linelist["element"]])
    m = (el == species)
    if not m.any():
        raise SystemExit(f"no {species!r} rows in the synthesis list")
    idx_all = np.flatnonzero(m)
    order = np.argsort(w_A[idx_all])
    idx_sorted = idx_all[order]
    w_sorted = w_A[idx_sorted]

    keep, missing = _match_into_list(want.values, w_sorted, idx_sorted)
    # LOUD, never silent (RYA-711). A line the EW leg measured that the synthesis list
    # does not contain cannot be compared, and the count is itself a result: it is the
    # NOT-IN-SYNTH-LINELIST population RYA-959's Engine-B already reported.
    if missing:
        print(f"  [lines-from-ew] {len(missing)} of {len(want)} EW lines are NOT in the "
              f"synthesis list (>{_EW_MATCH_TOL_A} A from any row) — they cannot be "
              f"compared and are reported, not dropped quietly. First few: "
              f"{[f'{w:.3f}(+{d:.3f})' for w, d in missing[:4]]}")
    if not keep:
        raise SystemExit("--lines-from-ew: no EW line matched the synthesis list")

    keep = np.array(sorted(set(keep)))
    df = pd.DataFrame({
        "wave_A": w_A[keep],
        "loggf": np.asarray(linelist["loggf"], dtype=float)[keep],
        "ep_eV": np.asarray(linelist["lower_state_eV"], dtype=float)[keep],
        "theo_depth": np.asarray(linelist["theoretical_depth"], dtype=float)[keep],
    }).sort_values("wave_A").reset_index(drop=True)

    if tier != "all":
        graded = _graded_mask(df.wave_A.values)
        df = df[graded if tier == "graded" else ~graded].reset_index(drop=True)
        print(f"  [tier] {tier}: {len(df)} lines (RYA-946 — graded and ungraded are "
              f"SEPARATE products, never merged)")
    return df


def _graded_mask(waves: np.ndarray) -> np.ndarray:
    """True where the line carries a PRIMARY LABORATORY gf in `canonical_gf`.

    Keyed on `gf_tier == LAB`, the column RYA-945 wrote, so "graded" means the same thing
    here as everywhere else rather than being re-derived from `loggf_reference` prose.
    """
    cg = pd.read_csv(ROOT / "data" / "linelists" / "canonical_gf.csv", low_memory=False)
    lab = cg[cg.gf_tier.astype(str).str.contains("LAB", na=False)]
    lw = np.sort(lab.wavelength_air_A.astype(float).values)
    i = np.searchsorted(lw, waves)
    best = np.full(waves.shape, np.inf)
    for off in (-1, 0, 1):
        j = np.clip(i + off, 0, len(lw) - 1)
        best = np.minimum(best, np.abs(lw[j] - waves))
    return best <= _EW_MATCH_TOL_A


def _selector_tag(a) -> str:
    """How the lines were chosen, as part of the artifact name — RYA-984.

    A product is identified by what it MEASURED, and the line set is part of that. Two
    runs differing only in selector are two different products and must not collide.
    The default selector is unlabelled so every pre-RYA-984 artifact keeps its name.
    """
    if getattr(a, "lines_deep_graded", False):
        return "_DEEPGRADED"
    if getattr(a, "lines_from_ew", None):
        tier = getattr(a, "lines_tier", "all")
        return "_FROMEW" + ("" if tier == "all" else f"-{tier.upper()}")
    if getattr(a, "lines_tier", "all") == "graded":
        return "_GRADED"
    return ""


def conditioning_tag(a) -> str:
    """How the OBSERVED spectrum was conditioned before the fit — RYA-1006.

    🔴 THE AXIS THAT WAS MISSING, AND IT DESTROYED AN ANCHOR BEFORE IT WAS ADDED.

    `_selector_tag` closed the SELECTION axis (RYA-984) after RYA-933/934 closed the
    HOLDING axis. Two axes opened after them and neither reached the stem: `--local-renorm`
    (RYA-1000) divides every fit window by its own 95th percentile, and `--degrade-to-R`
    (RYA-995) convolves the observed spectrum down. Both change WHAT WAS MEASURED while
    leaving element, band, instrument, holding, route and selection identical.

    On 2026-08-23 the consequence landed: two `--local-renorm` runs wrote
    `FeI_4200_6910_*_SYNTH_DEEPGRADED_*` — the exact filenames
    `pipeline.anchor_pools.ANCHORS['rya984_graded_163']` names — and Kitt Peak's anchor
    value moved **7.417 -> 7.337** under the anchor's own name. ⚠️ AND THE PROVENANCE FILE
    WAS BYTE-IDENTICAL: the renorm left no trace anywhere except the number, so nothing
    downstream could have detected it. That is why this function exists AND why
    `_conditioning_note` now writes the treatment into the provenance and
    `observed_conditioning` into every per-line row.

    The default is the empty string, so every product made without a conditioning flag
    keeps the name it has always had (RYA-984's rule, and the reason the anchor's own name
    is still correct once it is regenerated).
    """
    bits = []
    if getattr(a, "local_renorm", False):
        bits.append("LOCALRENORM")
    dR = getattr(a, "degrade_to_R", None)
    if dR:
        bits.append(f"R{int(round(float(dR)))}")
    return "".join(f"_{b}" for b in bits)


def conditioning_id(a) -> str:
    """The same fact as a machine-readable value for the per-line rows and the guard.

    `native` is the default, spelled out rather than left blank: a consumer must be able
    to tell "this spectrum was not conditioned" from "this product predates the column"
    (RYA-833 — an absence is a hypothesis, never a conclusion). `pipeline.anchor_pools`
    depends on exactly that distinction.
    """
    tag = conditioning_tag(a)
    return tag.lstrip("_").replace("_", "+").lower() if tag else "native"


def _conditioning_note(a) -> str:
    """The provenance sentence. Empty when nothing was done to the observed spectrum."""
    if not conditioning_tag(a):
        return ""
    bits = []
    if getattr(a, "local_renorm", False):
        bits.append("each fit window divided by its OWN 95th-percentile local continuum "
                    "(RYA-1000 --local-renorm)")
    dR = getattr(a, "degrade_to_R", None)
    if dR:
        bits.append(f"observed spectrum convolved down to R={float(dR):.0f} before "
                    f"fitting (RYA-995 --degrade-to-R)")
    return (" ⚠️ OBSERVED SPECTRUM CONDITIONED — DIAGNOSTIC PRODUCT, NOT A SCIENCE "
            "PRODUCT: " + "; ".join(bits) + ". This is not the same measurement as the "
            "unconditioned product of the same element, band, instrument and selection, "
            "and it carries a distinct artifact stem for that reason (RYA-1006).")


def synthesis_stem(a) -> str:
    """The synthesis product's filename stem — EVERY independently-varying axis, in one place.

    Hoisted out of `synthesis_route` by RYA-1006 so the guard test drives THIS expression
    rather than a copy of it. A reconstructed stem agrees with itself while the route
    writes something else (RYA-845), which is precisely the class of defect that let two
    runs share a filename in the first place.

    The axes, and the ticket that had to learn each one the hard way:
      element / ion / band   — always
      instrument             — always
      holding                — RYA-933/934 (telluric-corrected vs raw, same instrument)
      route (`_SYNTH`)       — always
      selection              — RYA-984 (shallow-graded vs deep-graded, same everything)
      observed conditioning  — RYA-1006 (local-renorm / degraded-R, same everything)
    """
    stem = f"{a.element}{a.ion}_{int(a.lo)}_{int(a.hi)}_{a.instrument}"
    if a.holding:
        stem += f"_{a.holding}"
    stem += "_SYNTH"
    stem += _selector_tag(a)
    stem += conditioning_tag(a)
    return stem


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
    # 🔴 RYA-995 — THE CONTROLLED RESOLUTION TEST. Degrading means TWO things together,
    # and doing only one of them would be a different experiment: the observed spectrum is
    # convolved down (below, in `_observed`) AND the fit is told the new R, so the
    # synthesis is broadened to match. Overriding R alone would broaden the model against
    # an un-degraded observation and manufacture a mismatch.
    _R_native = R
    _degrade_to = getattr(a, "degrade_to_R", None)
    if _degrade_to:
        if _degrade_to >= R:
            raise SystemExit(
                f"--degrade-to-R {_degrade_to:.0f} is not below {a.instrument}'s own "
                f"R={R:.0f}. Degrading is a one-way operation: resolution cannot be added "
                f"back, and a 'degrade' upward would silently be a no-op.")
        print(f"  [RYA-995 DEGRADE] observed spectrum convolved {R:.0f} -> "
              f"{_degrade_to:.0f}; the fit is told the degraded R. DIAGNOSTIC ONLY.")
        R = float(_degrade_to)

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
    # 🔴 RYA-967 — PARTIAL COVERAGE IS A MISLABEL, NOT A SMALL PROBLEM.
    #
    # The check above refuses only when the list has NO line in the window. That was
    # sufficient while every synth band shipped a list built for it; the VIS entry does
    # not. It reads the iSpec-vendored GES v6 list — deliberately, so Part B varies method
    # alone against the EW route — and that list starts at 4200 A. A VIS run over
    # 3780-6910 would therefore have synthesised 4200-6910, converged, and emitted a
    # product labelled 3780-6910: a number naming a range it was not measured over, which
    # is the RYA-911/913 defect exactly.
    #
    # Refused rather than trimmed. Silently narrowing the band would hide the gap in a
    # product nobody re-reads; the caller is told to ask for the range the list covers.
    _lo_cov, _hi_cov = float(_ll.wave_A.min()), float(_ll.wave_A.max())
    _gap_lo = max(0.0, _lo_cov - a.lo)
    _gap_hi = max(0.0, a.hi - _hi_cov)
    if _gap_lo > _LIST_COVERAGE_TOL_A or _gap_hi > _LIST_COVERAGE_TOL_A:
        raise SystemExit(
            f"{pol.name} synthesis list {cfg.linelist.name} covers "
            f"{_lo_cov:.1f}-{_hi_cov:.1f} A but the run asks for {a.lo:.1f}-{a.hi:.1f} A "
            f"— uncovered: {_gap_lo:.1f} A at the blue edge, {_gap_hi:.1f} A at the red. "
            f"Emitting a product labelled {a.lo:.0f}-{a.hi:.0f} A that was synthesised "
            f"over a narrower range is a number naming a range it was not measured over "
            f"(RYA-911/913). Re-run over the covered range, or extend the list. "
            f"{cfg.build_hint}")
    prov_gf = gf_provenance(float(_w.min()), float(_w.max()))
    print(f"\n[synthesis route — {pol.name}]  R={R:.0f}")
    print(f"  [gf] {prov_gf['detail']}")
    ctx = build_solar_context(a.element, R, linelist_file=str(cfg.linelist),
                              apply_canonical_gf=prov_gf["apply_canonical_gf"],
                              star=a.star)
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

    #: RYA-933/938 — band-wide continuum for a holding that ships NONE. Built lazily on
    #: first use and reused for every window, because the continuum is a property of the
    #: BAND, not of a 1.2 A fit window. Fitting it per window would place the "continuum"
    #: inside the line being measured.
    _placed_cont: dict = {}

    def _band_continuum(h):
        """CONTINUUM_PARAMS['solar'] spline over the whole band, fitted once.

        🔴 THE KNOT SPACING IS THE WHOLE CORRECTNESS. 100 A knots with a 95th-percentile
        upper envelope and 3-sigma rejection is the project's solar convention; a local
        estimator (the 95th percentile of a single fit window, as `--local-renorm` uses)
        is a DIFFERENT operation and is wrong here -- at ~1.2 A it rides down the wings of
        the line it is about to measure and returns a continuum biased low, which the
        fitter then closes by lowering A(X).
        """
        # RYA-1026: the LAST gate before a continuum of ours touches the product. The
        # callers already check `h.pre_normalised` -- but that flag is exactly what
        # RYA-929 got wrong, and a flag nothing cross-checks re-normalised KP2005 for
        # months in silence. This checks the ratified registry AGAINST the flag, so a
        # future flip raises DRIFT here rather than quietly re-enabling the placement.
        assert_not_renormalising(
            h.holding_id, pre_normalised=h.pre_normalised, fitting_continuum=True,
            where="derive_band_products._band_continuum (--place-continuum)")
        if "f" in _placed_cont:
            return _placed_cont["f"]
        from config.constants import CONTINUUM_PARAMS
        from scipy.interpolate import CubicSpline
        cp = CONTINUUM_PARAMS["solar"]
        bw = load_window_ex(a.instrument, 0.5 * (a.lo + a.hi),
                            0.5 * (a.hi - a.lo) + 1.0, holding=a.holding)
        W = np.asarray(bw.wave, float); F = np.asarray(bw.flux, float)
        keep = np.ones(W.size, bool)
        knot = float(cp["knot_spacing_A"]); pct = float(cp["upper_percentile"])
        cont = None
        for _ in range(int(cp["n_iter"])):
            edges = np.arange(W[0], W[-1] + knot, knot)
            kx, ky = [], []
            for lo_, hi_ in zip(edges[:-1], edges[1:]):
                m = keep & (W >= lo_) & (W < hi_)
                if m.sum() < 20:
                    continue
                kx.append(0.5 * (lo_ + hi_)); ky.append(float(np.percentile(F[m], pct)))
            if len(kx) < 4:
                raise SystemExit(
                    f"--place-continuum: only {len(kx)} usable {knot:.0f} A knots across "
                    f"{a.lo:.0f}-{a.hi:.0f} A -- too few to constrain a spline. Refusing "
                    f"to place a continuum the band cannot support.")
            cont = CubicSpline(kx, ky)(W)
            resid = F - cont
            keep = resid > -float(cp["sigma_clip"]) * float(np.std(resid[keep]))
        _placed_cont["f"] = (W, cont)
        print(f"  [continuum] PLACED on {h.holding_id} (ships none): "
              f"{len(kx)} knots x {knot:.0f} A, {pct:.0f}th pct, "
              f"{cp['n_iter']} iters, {cp['sigma_clip']:.1f} sigma "
              f"(CONTINUUM_PARAMS['solar'])")
        return _placed_cont["f"]

    def _observed(centre: float, pad: float):
        win = load_window_ex(a.instrument, centre, pad, holding=a.holding)
        h = win.holding
        if not h.pre_normalised and not getattr(a, "place_continuum", False):
            raise LookupError(
                f"holding {h.holding_id} is NOT continuum-normalised, and the synthesis "
                f"route fits observed flux against a synthesis normalised to unity. "
                f"Fitting it would spend A({a.element}) closing a continuum offset "
                f"(RYA-713/843). Normalise the product first, or pass --place-continuum "
                f"if this holding genuinely ships no continuum to double (RYA-911).")
        _served[h.holding_id] = _served.get(h.holding_id, 0) + 1
        _served_specs[h.holding_id] = h
        _window_median.append(float(np.median(win.flux)))
        _w, _f, _prov = win.wave, win.flux, win.provenance
        if not h.pre_normalised and getattr(a, "place_continuum", False):
            _CW, _CF = _band_continuum(h)
            _c = np.interp(np.asarray(_w, float), _CW, _CF)
            _bad = ~np.isfinite(_c) | (_c <= 0)
            if _bad.any():
                raise SystemExit(
                    f"--place-continuum: the band continuum is non-positive at "
                    f"{int(_bad.sum())} of {_bad.size} pixels near {centre:.3f} A. "
                    f"Refusing to divide by it.")
            _f = np.asarray(_f, float) / _c
            _prov = (f"{_prov} | RYA-933 CONTINUUM PLACED "
                     f"(CONTINUUM_PARAMS['solar'], band-wide spline)")
        if getattr(a, "local_renorm", False):
            # 🔴 RYA-1000 — ONE CONVENTION FOR BOTH ARMS, applied IDENTICALLY.
            #
            # Kitt Peak ships residual flux where unity IS the continuum by construction;
            # HARPS ships its own fitted continuum. Measured over the 108 common deep
            # lines the local flux sits at 0.9281 (KP) against 0.8913 (HARPS) — and the
            # synthesis route fits observed flux against a synthesis normalised to UNITY,
            # so a depressed local continuum reads as extra absorption and the fitter
            # spends A(Fe) closing it (RYA-843's mechanism, RYA-911's -0.34 dex).
            #
            # The 95th percentile is `verify_feature`'s own continuum estimator, reused
            # rather than invented so the two places cannot disagree about what
            # "continuum" means in a window.
            #
            # ⚠️ SYMMETRIC BY DESIGN. Applying this to HARPS alone would be tuning one arm
            # toward the other. Both arms get the same operation; if the gap is a
            # convention artifact it collapses, and if it is real it survives.
            _c95 = float(np.nanpercentile(_f, 95))
            if np.isfinite(_c95) and _c95 > 0:
                _f = np.asarray(_f, float) / _c95
                _prov = f"{_prov} | RYA-1000 LOCAL-RENORM /p95={_c95:.4f}"
        if _degrade_to:
            # Gaussian kernel that takes R_native -> R_target:
            #   FWHM_deg^2 = (lam/R_target)^2 - (lam/R_native)^2
            # Convolving by the QUADRATURE DIFFERENCE, not by lam/R_target, is the whole
            # correctness of this: the spectrum already carries its native instrumental
            # profile, and convolving by the full target width would over-broaden it.
            from scipy.ndimage import gaussian_filter1d
            _lam = float(np.median(_w))
            _fw = _lam * np.sqrt(1.0 / _degrade_to ** 2 - 1.0 / _R_native ** 2)
            _px = float(np.median(np.diff(_w)))
            _sig_px = (_fw / 2.35482) / _px
            _f = gaussian_filter1d(np.asarray(_f, float), _sig_px, mode="nearest")
            _prov = f"{_prov} | RYA-995 DEGRADED to R={_degrade_to:.0f} (sigma {_sig_px:.2f} px)"
        return _w, _f, _prov

    _probe = select_holding(a.instrument, 0.5 * (a.lo + a.hi), hw + 0.4,
                            holding=a.holding)
    # Asked at band centre so the refusal costs nothing and arrives BEFORE any synthesis.
    # `_observed` checks it too, per line, because a band can in principle select a
    # different holding line by line -- but a per-line failure there surfaces as
    # `status: no_atlas` on each line, which describes coverage rather than the real
    # fault. The loud version comes first (RYA-832).
    if not _probe.pre_normalised and not getattr(a, "place_continuum", False):
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
    if getattr(a, "lines_deep_graded", False):
        cand = _cand_deep_graded(ctx["linelist"], lo_A=a.lo, hi_A=a.hi,
                                 species=species_token(a.element, a.ion))
    elif getattr(a, "lines_from_ew", None):
        cand = _cand_from_ew_artifact(ctx["linelist"], Path(a.lines_from_ew),
                                      lo_A=a.lo, hi_A=a.hi,
                                      species=species_token(a.element, a.ion),
                                      tier=a.lines_tier)
    elif getattr(a, "lines_tier", "all") == "graded":
        cand = _cand_graded(ctx["linelist"], lo_A=a.lo, hi_A=a.hi,
                            species=species_token(a.element, a.ion))
    elif getattr(a, "lines_tier", "all") == "ungraded":
        raise SystemExit(
            "--lines-tier ungraded is only meaningful with --lines-from-ew (it splits "
            "an EW artifact's pool). Standalone, use the default 'all'.")
    else:
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
    # 🔴 RYA-1024 — THE QUARANTINE NOW ASKS THE HOLDING, NOT JUST THE INSTRUMENT.
    # `telluric_reason` keys on the INSTRUMENT's telluric_basis, which for
    # kpno_solar_atlas is `line_selection`. That is right for the raw atlas and WRONG for
    # a corrected sibling: the corrected bands and the quarantined bands are the SAME
    # bands, so the only lines telluric correction unlocks were exactly the lines the
    # selector refused, and `solar_kpno_molecfit_corrected` could never measure anything.
    #
    # The lift is granted on EVIDENCE, never on the holding's `applied` label -- five
    # holdings declare `applied` and solar_iag is one of them while reading the RAW
    # Reiners file. `serves_corrected_flux` answers the narrow, checkable question:
    # is there a RYA-940 corrected PRODUCT covering this wavelength? H2O 7160-7340 got no
    # admissible fit, so it has no file, so the quarantine still stands there.
    from measure_band_ew import telluric_reason, serves_corrected_flux
    _tell = [(float(r.wave_A), telluric_reason(float(r.wave_A), a.instrument))
             for r in cand.itertuples()]
    _lifted = [(w, serves_corrected_flux(a.holding, w)) for w, why in _tell if why]
    _lifted = [(w, prov) for w, prov in _lifted if prov]
    _lift_at = {w for w, _ in _lifted}
    if _lifted:
        print(f"  {len(_lifted)} candidate(s) RESTORED by a corrected holding "
              f"(RYA-1024 -- the band is corrected, so the flux IS stellar):")
        for w, prov in _lifted:
            print(f"      {w:10.3f}  served corrected by {prov}")
    _bad = [(w, why) for w, why in _tell if why and w not in _lift_at]
    if _bad:
        print(f"  {len(_bad)} candidate(s) QUARANTINED-TELLURIC before fitting "
              f"(flux there is not stellar; excluded, never corrected):")
        for w, why in _bad:
            print(f"      {w:10.3f}  {why.split(' and ')[0]}")
    if _bad:
        _drop = {w for w, _ in _bad}
        cand = cand[~cand.wave_A.astype(float).isin(_drop)].reset_index(drop=True)
    print(f"  {len(cand)} {a.element} {a.ion} candidates by theoretical depth "
          f"(half-width +/-{hw} A, min separation {cfg.min_sep_A} A)")
    print(f"  [half-width] {cfg.half_width_note}")

    tmp = f"/tmp/ispec_synth_{pol.name.replace(chr(47), chr(45))}"
    Path(tmp).mkdir(parents=True, exist_ok=True)
    # 🔴 RYA-1044 — ONE FIT LOOP, CALLED TWICE. The Engine-B leg below re-fits THESE SAME
    # lines with the departures applied, and it must do so through THIS code rather than a
    # copy of it. RYA-701 is the standing reason (one Ba->Al copy produced thirteen
    # defects); here there is a sharper one -- the two legs are DIFFERENCED against each
    # other, so any drift between them lands directly in the reported NLTE effect rather
    # than in one product's value. A copy could not be wrong in a way that cancels.
    def _fit_lines(treatment: str, **fit_kw) -> list[LineMeasurement]:
        lines: list[LineMeasurement] = []
        for r in cand.itertuples():
            w = float(r.wave_A)
            res = fit_one(ctx, segs, w, hw, tmp, load=_observed, **fit_kw)
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
                treatment=treatment,   # RYA-1044: the leg decides, not the loop
                # RYA-1006 — the row says what was done to the spectrum it was fitted against.
                observed_conditioning=conditioning_id(a),
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
        return lines

    # The 1D-LTE leg -- unchanged. This is the call RYA-759 published against, and the
    # only difference from before RYA-1044 is that its body now lives in a function the
    # Engine-B leg can also call.
    lines = _fit_lines("1D-LTE")

    # ── RYA-1044: the Engine-B leg's INPUTS ─────────────────────────────────────────
    # ⚠️ OPT-IN, AND THE DEFAULT-OFF IS THE BLAST-RADIUS GUARD, NOT TIMIDITY. Every
    # `--force-synthesis` band in the repo reaches this function, and making the leg
    # automatic would (a) add a product row to bands whose GENERATORS entries do not
    # declare one -- which is a HARD CI failure, `check_result_generators` fails in both
    # directions -- and (b) roughly double every such run. So the leg is REACHABLE (the
    # ticket's word) rather than unconditional, and turning it on by default is a
    # follow-up that has to land its GENERATORS entries in the same change.
    eb_lines = eb_treatment = eb_prov = None
    if getattr(a, "synth_engine_b", False) and not a.skip_engine_b:
        _mean3d = a.engine_b_deck in ("gerber-mean3d", "gerber-mean3d-lte")
        _nlte = a.engine_b_deck in ("gerber-nlte", "gerber-mean3d")
        _fit_kw: dict = {}
        if _mean3d:
            from pipeline import gerber_nlte as gnlte
            from pipeline import mean3d_atmosphere as m3d
            _deck_key = f"{a.element}@mean3D"
            if _deck_key not in gnlte.DECKS:
                raise SystemExit(
                    f"no <3D> deck registered for {a.element} (looked for "
                    f"{_deck_key!r}). Registered: {sorted(gnlte.DECKS)}. Staging is "
                    f"RYA-710; which aux file to register against is a committed row in "
                    f"data/results/rya1035/mean3d_aux_defect_sweep.csv -- Fe and Mn "
                    f"REQUIRE the plain aux.")
            _layers, _model3d = m3d.load(gnlte.deck_atmosphere(_deck_key))
            _dep = gnlte.for_node(_deck_key, float(ctx["teff"]), float(ctx["logg"]),
                                  float(ctx["feh"]))
            gnlte.assert_depth_match(_dep, _layers)
            _dtau = m3d.assert_tau_consistent(_dep, _model3d)
            # The <3D> atmosphere REPLACES ctx's for this leg, and the model goes in as a
            # FILE: a five-column mul23 model cannot be written in the MARCS form iSpec's
            # writer produces (see pipeline.mean3d_atmosphere).
            _fit_kw = dict(atmosphere=_layers,
                           atmosphere_layers_file=_model3d["path"])
            if _nlte:
                _fit_kw.update(nlte_deck="gerber", nlte_deck_key=_deck_key)
            eb_treatment = (taxes.MEAN3D_NLTE_STAGGER if _nlte
                            else taxes.MEAN3D_LTE_STAGGER).token
            eb_prov = MEAN3D_SOURCE_FMT.format(
                atom=_dep["atom_path"].split("/")[-1],
                model=Path(_model3d["path"]).name,
                mode=("NLTE, departures applied" if _nlte
                      else "LTE, departures WITHHELD -- the paired comparand"))
            print(f"\n[Engine-B] <3D> deck {_deck_key} on "
                  f"{Path(_model3d['path']).name} (ndep={_dep['ndep']}, "
                  f"nlev={_dep['nk']}, max|dlogtau|={_dtau:.2e}, "
                  f"departures={'ON' if _nlte else 'OFF'})")
        elif _nlte:
            from pipeline import gerber_nlte as gnlte
            _dep = gnlte.for_node(a.element, float(ctx["teff"]), float(ctx["logg"]),
                                  float(ctx["feh"]))
            gnlte.assert_depth_match(_dep, ctx["atmosphere"])
            _fit_kw = dict(nlte_deck="gerber", nlte_deck_key=a.element)
            eb_treatment = "ENGINE-B-NLTE"
            eb_prov = GERBER_NLTE_SOURCE_FMT.format(
                atom=_dep["atom_path"].split("/")[-1])
        else:
            eb_treatment = "ENGINE-B"
            eb_prov = ("Turbospectrum LTE flux fit on this route's own atmosphere "
                       "(RYA-1044): the same lines, fitter and window as the 1D-LTE leg "
                       "above, emitted as a SEPARATE product rather than a correction "
                       "to it (RYA-712).")
        if _nlte:
            # iSpec fails SOFT -- an unlabelled element lands in `nlte_ignored` and
            # synthesises in LTE WITHOUT raising (RYA-764). Checked before the fit.
            from pipeline import gerber_nlte as _gn
            _gn.assert_linelist_supports_nlte(
                ctx["linelist"], int(ctx["atom_code"]), a.element)
        print(f"[Engine-B] fitting {len(cand)} lines as {eb_treatment} ...")
        eb_lines = _fit_lines(eb_treatment, **_fit_kw)

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
            "— 901 Fe I candidates, 0 measurable. 1D-LTE ONLY: this route has no "
            f"in-synthesis NLTE deck registered for {a.element} {a.ion}; the LTE value "
            "must therefore be reported as its own product and not interpreted as an "
            "NLTE expectation. "
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
        "iSpec, ATLAS9.Castelli) — the RYA-759 validated route called directly, retaining "
        "its measured window and constraint choices while fitting the requested element. "
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
        "gf TERM: " + rung.describe() + ". Never coadded with another band (RYA-712)."
        # RYA-1006 — LAST, so it is the final thing a reader sees. Empty on an
        # unconditioned product, so existing provenance text is byte-unchanged.
        + _conditioning_note(a))
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
    # RYA-933/934 -- the HOLDING belongs in the stem. Two holdings of one instrument
    # differing by whether tellurics were removed would otherwise write the same
    # filename, and the second would silently overwrite the first.
    # 🔴 RYA-984 — THE SELECTOR BELONGS IN THE STEM, and it did not used to be.
    #
    # This route now has three line-selection modes, and two of them produce a DIFFERENT
    # SCIENTIFIC PRODUCT over the same element, band, instrument and holding: RYA-967's
    # shallow graded set (55 lines, the EW-comparable ones) and RYA-984's deep graded set
    # (109 lines, the ones EW can never attempt). Both wrote
    # `FeI_4200_6910_kpno_solar_atlas_SYNTH_*`, so the second silently overwrote the first
    # — the RYA-933/934 defect exactly, which that ticket fixed for the HOLDING axis and
    # left open for the SELECTION axis.
    #
    # Caught 30 seconds into the first deep run, before it clobbered a committed product.
    stem = synthesis_stem(a)
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

    # ── RYA-1044: the Engine-B leg, REACHABLE at last ────────────────────────────────
    # 🔴 THE DEFECT THIS CLOSES. `--force-synthesis` returns from this function, and this
    # function had NO Engine-B code -- so the leg RYA-784 ("wire Engine-B into
    # derive_band_products") and RYA-798 ("wire the Gerber decks into the production
    # flux-fit path") both closed as Done was unreachable from any synthesis-only band.
    # Measured across every committed *_products.csv before this change: 1D-LTE x12,
    # ENGINE-A x7, and NOT ONE Engine-B product, for any element, band or deck. Two
    # tickets, both Done, zero products -- the wiring landed in the MAIN route while the
    # bands that can only synthesise take this one.
    #
    # 🔴 AND THERE IS NO MEASURED-EW PRECONDITION HERE, DELIBERATELY. The main route
    # refuses without a profile-fit EW file; a flux fit does not consume an EW, so that
    # is a precondition of the EW ROUTE and requiring it to emit a synthesis product was
    # false coupling. Removing it is the point of this ticket, not a side effect.
    if eb_lines is not None:
        eb_product = build_product(a.element, a.ion, a.instrument, pol.name,
                                   eb_treatment, eb_lines,
                                   handler="SynthesisHandler", provenance=eb_prov)
        eb_rung = gf_rung.for_lines(a.element, a.ion, eb_lines, linelist=ctx["linelist"])
        eb_b = build_budget(a.element, 0.5 * (a.lo + a.hi), eb_product.n_lines,
                            scatter_dex=eb_product.sigma,
                            **eb_rung.budget_kwargs(),
                            **harness_residual.for_product(eb_product).budget_kwargs())
        eb_stat, eb_syst = eb_b.total()
        assert_stat_publishable(
            eb_stat, cell=f"{a.element} {a.ion} {pol.name} {eb_treatment} "
                          f"({a.instrument}, n={eb_product.n_lines})")
        # 🔴 ITS OWN FILES, exactly as the ENGINE-A leg above does. Appending a row to
        # `{stem}_products.csv` would change the ROW COUNT of a file every existing band
        # writes with one row -- which is a shift in an existing output, and the ticket's
        # blast-radius guard forbids exactly that. A separate file is purely additive:
        # a band that does not opt in is byte-identical, and one that does gains files
        # rather than altering them.
        pd.DataFrame([asdict_line(l) for l in eb_lines]).to_csv(
            out / f"{stem}_{eb_treatment}_lines.csv", index=False)
        (out / f"{stem}_{eb_treatment}_budgets.txt").write_text(
            eb_b.describe() + f"\n  gf rung: {eb_rung.describe()}\n")
        (out / f"{stem}_{eb_treatment}_provenance.txt").write_text(
            eb_product.provenance + "\n")
        pd.DataFrame([dict(
            element=a.element, ion=a.ion, band=pol.name, instrument=a.instrument,
            treatment=eb_treatment, handler=eb_product.handler,
            A=round(eb_product.value, 3) if eb_product.value is not None else None,
            n_lines=eb_product.n_lines, n_excluded=eb_product.n_excluded,
            stat_dex=round(eb_stat, 4), syst_dex=round(eb_syst, 4),
            stat_basis=eb_b.stat_basis(),
            dominant=(eb_b.dominant().name if eb_b.dominant() else ""),
            # The axes come from the SAME registry the token does, so a product cannot
            # carry a label the vocabulary has never heard of (the RYA-798 failure, where
            # a treatment was emitted that `TREATMENTS` did not contain and the product
            # died at `build_product` after the synthesis had already run).
            **axes_for(eb_treatment, handler=eb_product.handler or None).as_columns(),
        )]).to_csv(out / f"{stem}_{eb_treatment}_products.csv", index=False)
        print(f"\n  A({a.element} {a.ion}; {pol.name}, {eb_treatment}) = "
              f"{eb_product.value if eb_product.value is not None else float('nan'):.3f}"
              f"  (n={eb_product.n_lines}, excluded {eb_product.n_excluded})")
        if product.value is not None and eb_product.value is not None:
            # The differential this pair EXISTS for: both legs on ONE atmosphere, the
            # same lines, the same fitter -- so the difference is the departures and
            # nothing else (RYA-542).
            print(f"    delta vs 1D-LTE on this route = "
                  f"{eb_product.value - product.value:+.4f} dex")

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
    # ── ENGINE-B-NLTE on the SYNTHESIS route ─────────────────────────────────
    #
    # 🔴 RYA-1040/1043 — `--engine-b-deck` WAS SILENTLY INERT ON THIS ROUTE.
    # `main()` returns into `synthesis_route` before the EW route's Engine-B block, and
    # `fit_one` never forwarded `nlte_deck` to `_fit_synth_flux` (which has accepted it
    # all along). So the flag parsed, was accepted, and did nothing. That is why the
    # Gerber column is empty on EVERY graded and deep-graded product we hold: every one
    # of them is a synthesis-route product.
    #
    # This is a SEPARATE PRODUCT, never a correction applied to the LTE value (RYA-712) —
    # the departures enter the fit itself, so the number is a different measurement of the
    # same line rather than a delta on top of one.
    rows_bn: list[LineMeasurement] = []
    if getattr(a, "engine_b_deck", "ts-lte") != "ts-lte":
        deck = a.engine_b_deck
        print(f"\n[3] ENGINE-B-NLTE — re-fitting in-synthesis with deck={deck}")
        for l in [x for x in lines if x.in_aggregate and x.abundance is not None]:
            try:
                rb = fit_one(ctx, segs, l.wavelength_air_A, hw, tmp,
                             load=_observed, nlte_deck=deck)
            except Exception as exc:
                # LOUD per line. A deck that cannot serve a line is a stated absence,
                # never a silent fallback to the LTE number (RYA-833).
                rows_bn.append(replace(l, treatment="ENGINE-B-NLTE", abundance=None,
                                       in_aggregate=False,
                                       excluded_reason=f"ENGINE-B-NLTE: {str(exc)[:160]}",
                                       nlte_source=f"{deck} (deck error)"))
                continue
            ab = float(rb.get("a_synth", float("nan")))
            rows_bn.append(replace(
                l, treatment="ENGINE-B-NLTE",
                abundance=(ab if np.isfinite(ab) else None),
                in_aggregate=bool(np.isfinite(ab)),
                # ⚠️ NOT a delta. RYA-798: the departures are applied INSIDE the fit, so
                # there is no additive per-line correction to record, and writing one
                # would invite reading it as an Engine-B-LTE offset.
                nlte_delta_dex=None,
                nlte_source=f"in-synthesis departures, deck={deck} (RYA-798)"))
        n_ok = sum(1 for x in rows_bn if x.in_aggregate)
        print(f"    ENGINE-B-NLTE: {n_ok} fitted, {len(rows_bn) - n_ok} refused")
        lines = lines + rows_bn

    # ── ENGINE-A on the SYNTHESIS route ──────────────────────────────────────
    #
    # 🔴 RYA-1002 — THIS BLOCK DID NOT EXIST, AND ITS ABSENCE WAS SILENT.
    #
    # `[2] ENGINE-A` below is on the EW route only, so a band measured with
    # `--force-synthesis` emitted 1D-LTE and stopped — no NLTE product, no message, no
    # disposition. RYA-525's two-engine floor requires the departure product, and a
    # reader of the matrix could not tell "this element has no registered grid" from
    # "the route never asked". Found on Al: the Amarsi-2020 grid serves all six Al I
    # lines AT SOLAR PARAMETERS (RYA-773 landed), and the red-optical run still returned
    # 1D-LTE alone.
    #
    # A departure term added to the SAME inversion. The handler is this route's own
    # (`SynthesisHandler`) and NOT the EW route's ProfileFitHandler — ENGINE-A rides
    # whatever measured the line, and charging the profile fitter's residual to a flux
    # fit would attribute someone else's systematic to it (RYA-869).
    used_a = [l for l in lines if l.in_aggregate and l.abundance is not None]
    deltas = engine_a_delta(a.element, a.ion,
                            np.array([l.wavelength_air_A for l in used_a]),
                            cache=a.mpia_cache, star=a.star)
    rows_a: list[LineMeasurement] = []
    for l in used_a:
        # Tolerance match — two catalogues quoting one transition differ by up to
        # ~0.03 A (RYA-704), and the grid keys are rounded.
        d = None
        if deltas:
            k = min(deltas, key=lambda x: abs(x - l.wavelength_air_A))
            if abs(k - l.wavelength_air_A) < 0.06:
                d = deltas[k]
        la = LineMeasurement(
            element=l.element, ion=l.ion, wavelength_air_A=l.wavelength_air_A,
            instrument=l.instrument, ew_mA=l.ew_mA, ew_method=l.ew_method,
            treatment="ENGINE-A", ep_eV=l.ep_eV,
            nlte_delta_dex=(float(d) if d is not None else None),
            nlte_source=(engine_a_source(a.element) if d is not None
                         else "registered per-line delta_nlte (NOT SERVED)"),
            continuum_level=l.continuum_level, continuum_method=l.continuum_method,
            continuum_ref=l.continuum_ref,
            observed_depth=l.observed_depth, implied_width_A=l.implied_width_A,
            red_chi2=l.red_chi2,
            abundance=(l.abundance + d) if d is not None else None)
        if d is None:
            la.in_aggregate = False
            la.excluded_reason = (
                "ENGINE-A-NOT-SERVED: the registered source returns no usable "
                "delta_nlte for this line (absent, nan, or a placeholder zero). "
                "Reduced coverage, not a failed correction.")
        # `_stamp` is a closure over the EW route's `main()` locals and is not reachable
        # here. Both invariants it enforces are still honoured, explicitly:
        #   RYA-880 — every row must STATE what departure was applied, because a blank
        #             cannot be told apart from an unrecorded correction (RYA-833).
        #             Satisfied by construction above, and asserted rather than assumed.
        #   RYA-711 — the problem-children registry disposition is carried, not dropped.
        #             ENGINE-A is the SAME transition as its parent line (RYA-871), so it
        #             inherits that line's verdict; re-deriving it would be a second
        #             source for one fact.
        if not la.nlte_source:
            raise SystemExit(
                f"RYA-880: {la.element} {la.ion} {la.wavelength_air_A} (ENGINE-A) has no "
                f"`nlte_source`. Every row must state what departure was applied.")
        la.problem_class = l.problem_class
        la.problem_status = l.problem_status
        la.problem_tickets = l.problem_tickets
        la.problem_action = l.problem_action
        if not l.in_aggregate:
            la.in_aggregate = False
            la.excluded_reason = (l.excluded_reason if not la.excluded_reason
                                  else f"{l.excluded_reason} | {la.excluded_reason}")
        rows_a.append(la)
    p_a = build_product(a.element, a.ion, a.instrument, pol.name, "ENGINE-A", rows_a,
                        handler="SynthesisHandler",
                        provenance=engine_a_source(a.element))
    if p_a.value is not None:
        pd.DataFrame([asdict_line(l) for l in rows_a]).to_csv(
            out / f"{stem}_ENGINE-A_lines.csv", index=False)
        b_a = build_budget(a.element, 0.5 * (a.lo + a.hi), p_a.n_lines,
                           scatter_dex=p_a.sigma, **rung.budget_kwargs(),
                           **harness_residual.for_product(p_a).budget_kwargs())
        stat_a, syst_a = b_a.total()
        pd.DataFrame([dict(
            element=a.element, ion=a.ion, band=pol.name, instrument=a.instrument,
            treatment="ENGINE-A", handler=p_a.handler,
            A=round(p_a.value, 3), n_lines=p_a.n_lines, n_excluded=p_a.n_excluded,
            stat_dex=round(stat_a, 4), syst_dex=round(syst_a, 4),
            stat_basis=b_a.stat_basis(),
            dominant=(b_a.dominant().name if b_a.dominant() else ""),
            **axes_for("ENGINE-A", handler=p_a.handler or None).as_columns(),
        )]).to_csv(out / f"{stem}_ENGINE-A_products.csv", index=False)
        (out / f"{stem}_ENGINE-A_budgets.txt").write_text(
            b_a.describe() + f"\n  gf rung: {rung.describe()}\n")
        (out / f"{stem}_ENGINE-A_provenance.txt").write_text(p_a.provenance + "\n")
        _d = [float(r.nlte_delta_dex) for r in rows_a
              if r.nlte_delta_dex is not None]
        print(f"\n  A({a.element} {a.ion}; {pol.name}, ENGINE-A) = {p_a.value:.3f}  "
              f"(n={p_a.n_lines}, not-served {p_a.n_excluded})")
        print(f"    delta_nlte applied: mean {np.mean(_d):+.4f} dex over "
              f"{len(_d)} line(s) — {engine_a_source(a.element)}")
    else:
        # An UNSERVED element is a stated disposition, not a missing row. The old
        # behaviour — printing nothing — is what made this gap invisible.
        print(f"\n  ENGINE-A: NOT PRODUCED — {engine_a_source(a.element)}; "
              f"{p_a.n_excluded} of {len(rows_a)} line(s) unserved. The 1D-LTE product "
              f"above stands alone and is NOT an NLTE value.")

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
#: RYA-1040 — the <3D> pair's provenance string. It names the ATMOSPHERE as well as the
#: atom, because that is the axis distinguishing this product from its 1D sibling, and it
#: states in words whether the departures were applied: the LTE member runs the SAME deck
#: setup with them withheld, and a reader must not have to infer that from the label.
MEAN3D_SOURCE_FMT = (
    "gerber:{atom} on <3D> {model} (STAGGERmean3D deck read directly, RYA-821/1035; "
    "no vendor interpolator); {mode}. The NLTE effect is (<3D>-NLTE minus <3D>-LTE) on "
    "this one atmosphere -- never against a 1D-LTE product, which would report the "
    "1D->mean-3D atmosphere shift as non-LTE physics (RYA-542).")

GERBER_NLTE_SOURCE_FMT = ("gerber:{atom} (TS-native departure deck, RYA-798); departures "
                          "applied INSIDE the synthesis — no additive per-line delta "
                          "exists, this is a separate product, not a corrected LTE value")


def engine_a_source(element: str) -> str:
    """Name the registered source actually consulted for an Engine-A product."""
    if element == "Fe":
        return MPIA_NLTE_SOURCE
    from config.constants import NLTE_CORRECTION_ELEMENTS
    spec = NLTE_CORRECTION_ELEMENTS.get(element)
    if spec is None:
        return "no registered Engine-A source"
    return f"{spec['grid']} ({spec['ref']}); PER-LINE additive correction"


def engine_a_model(element: str) -> str:
    """Physics-model axis for the registered Engine-A source."""
    if element == "Fe":
        return "bergemann"
    from config.constants import NLTE_CORRECTION_ELEMENTS
    spec = NLTE_CORRECTION_ELEMENTS.get(element, {})
    source = f"{spec.get('grid', '')} {spec.get('ref', '')}".lower()
    if "amarsi" in source or "nordlander" in source:
        return "amarsi"
    if "bergemann" in source or "mpia" in source:
        return "bergemann"
    raise ValueError(f"no RYA-906 model-axis mapping for {element}: {source!r}")


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
                # RYA-1006 — and it must be added HERE too, for the reason the comment
                # above already gives: this projection is the lossy serialiser, and it is
                # the ONE the anchor guard reads.
                observed_conditioning=l.observed_conditioning,
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
                profile_gamma_A=l.profile_gamma_A,
                # RYA-959 — the observed core and the width the integrated EW implies
                # against it, carried for the same reason as the two fields above and
                # with a sharper one: RYA-958 diagnosed a contaminated pool by
                # RECONSTRUCTING this quantity by hand, from a MODEL depth out of
                # linelist_solar.csv, because no artifact in the chain recorded the
                # observed one. Stopping it at the EW file would guarantee the next RCA
                # reconstructs it again — the RYA-843 `red_chi2` omission exactly.
                observed_depth=l.observed_depth, implied_width_A=l.implied_width_A)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--star", default="solar",
                    help="star id from config/stars.yaml (RYA-298 single source). Default "
                         "'solar' keeps every existing invocation bit-identical; any other "
                         "value drives BOTH the atmosphere and the NLTE grid node (RYA-985).")
    ap.add_argument("--element", required=True)
    ap.add_argument("--ion", default="I")
    ap.add_argument("--lo", type=float, required=True)
    ap.add_argument("--hi", type=float, required=True)
    ap.add_argument("--instrument", default="kpno_solar_atlas")
    ap.add_argument("--holding", default=None,
                    help="Name ONE holding explicitly instead of taking the instrument's "
                         "first covering candidate (RYA-933/934). Without it, an "
                         "instrument that serves several products silently answers with "
                         "whichever is listed first -- which is how two telluric-corrected "
                         "Kitt Peak holdings sat unmeasurable while `--instrument "
                         "kpno_solar_atlas` looked like it covered them.")
    ap.add_argument("--lines-from-ew", default=None, metavar="EW_CSV",
                    help="RYA-967: drive the SYNTHESIS route over the lines a named EW "
                         "artifact attempted, instead of `select_lines`' strongest-N. "
                         "This is what makes an EW-vs-synth comparison CONTROLLED: the "
                         "two legs then differ by method and by nothing else. Without it "
                         "the synth leg picks its own (stronger) lines and any difference "
                         "in A(X) conflates method with selection — RYA-842, line "
                         "selection dominates.")
    ap.add_argument("--degrade-to-R", type=float, default=None, metavar="R",
                    help="RYA-995: convolve the OBSERVED spectrum down to this resolving "
                         "power before fitting, and fit at it. The controlled "
                         "instrument test — same spectrum, same lines, same pipeline, "
                         "only the resolution changes, so nothing else can explain a "
                         "difference. NOT for producing science products.")
    ap.add_argument("--local-renorm", action="store_true",
                    help="RYA-1000: divide each fit window by its OWN local continuum "
                         "(95th percentile) before fitting, so both arms share one "
                         "continuum convention. DIAGNOSTIC ONLY — the controlled test of "
                         "whether the KP-vs-HARPS arm gap is a normalisation artifact.")
    ap.add_argument("--place-continuum", action="store_true",
                    help="⚠️ CHECK FOR A SHIPPED CONTINUUM FIRST. Kurucz 2005 was served "
                         "through this flag until RYA-933 found the distribution also "
                         "ships `irradrelwl.dat`, a residual atlas with Kurucz's own "
                         "continuum divided out -- the readme does not list it, so the "
                         "intake missed it. The placed spline tilted the band 4% "
                         "blue-to-red and cost 0.0218 dex. A product that ships a "
                         "continuum must use it (RYA-911/938). "
                         "RYA-933/938: place a continuum on a holding that ships NONE "
                         "(pre_normalised=False), instead of refusing it. Fitted ONCE "
                         "across the whole band with CONTINUUM_PARAMS['solar'] (100 A "
                         "knots, 95th pct, 5 iters, 3 sigma) -- never per fit window, "
                         "which at ~1.2 A would ride down the wings of the very lines "
                         "being measured and bias A(X) low. Only legitimate where the "
                         "product genuinely ships no continuum to double (RYA-911): "
                         "Kurucz 2005 is absolute irradiance in W/m2/nm.")
    ap.add_argument("--lines-deep-graded", action="store_true",
                    help="RYA-984: select the laboratory-graded lines that sit ABOVE the "
                         "EW depth gate — the ones EW can never attempt because the curve "
                         "of growth is flat there. Synthesis inverts no EW, so the "
                         "saturation ceiling does not apply to them.")
    ap.add_argument("--lines-tier", choices=["all", "graded", "ungraded"], default="all",
                    help="RYA-946 two-tier: emit the graded (primary-laboratory gf) lines "
                         "as their own product, never merged with the ungraded ones.")
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
    ap.add_argument("--engine-b-deck",
                    choices=["ts-lte", "gerber-nlte", "gerber-mean3d", "gerber-mean3d-lte"],
                    default="ts-lte",
                    help="which synthesis is Engine B. 'ts-lte' is the production "
                         "Turbospectrum LTE flux-fit; 'gerber-nlte' is the TS-native "
                         "Gerber departure deck (validated RYA-785, wired RYA-798). "
                         "It runs on MARCS.GES to match the deck and emits a SEPARATE "
                         "ENGINE-B-NLTE product -- never a correction to Engine-B-LTE. "
                         "RYA-1040: 'gerber-mean3d' runs the <3D> STAGGERmean3D deck on "
                         "the <3D> STAGGER atmosphere those departures were solved on, "
                         "and 'gerber-mean3d-lte' runs THE SAME deck's setup with the "
                         "departures OFF. Those two are a MANDATORY PAIR: the NLTE effect "
                         "is (<3D>-NLTE minus <3D>-LTE) on ONE atmosphere, because "
                         "differencing against 1D-LTE would report the 1D->mean-3D "
                         "ATMOSPHERE shift as non-LTE physics (RYA-542).")
    ap.add_argument("--synth-engine-b", action="store_true",
                    help="RYA-1044: run the Engine-B leg on the SYNTHESIS route too "
                         "(--force-synthesis bands). That leg was unreachable -- "
                         "`synthesis_route` had no Engine-B code, so RYA-784 and RYA-798 "
                         "both closed Done while the leg had never emitted a product for "
                         "any element, band or deck. Unlike the EW route it needs NO "
                         "measured-EW file: a flux fit does not consume an EW. ⚠️ OFF by "
                         "default so the blast radius on existing --force-synthesis "
                         "bands is exactly zero -- switching it on adds a product row, "
                         "and a row with no GENERATORS entry is a hard CI failure. "
                         "Which product it emits is chosen by --engine-b-deck.")
    ap.add_argument("--skip-engine-b", action="store_true",
                    help="derive 1D-LTE and ENGINE-A only. Engine B refits the spectrum "
                         "per line and is much slower than the EW inversion.")
    a = ap.parse_args()

    pol_early = resolve_band(0.5 * (a.lo + a.hi))
    # 🔴 RYA-1006 — the conditioning flags are honoured ONLY by `synthesis_route`. On the
    # EW route they were parsed and silently ignored, so `--local-renorm` there produced
    # an UNCONDITIONED product from a command line that says it is conditioned. The stem
    # would carry no conditioning tag either, which is correct for what ran and wrong for
    # what was asked — the same silent-wrong class this ticket exists to close. Refuse.
    if conditioning_tag(a) and not (
            a.force_synthesis
            or "profile-fit" in getattr(pol_early, "forbidden_methods", ())):
        raise SystemExit(
            f"--local-renorm / --degrade-to-R condition the OBSERVED spectrum before a "
            f"synthesis fit, and only `synthesis_route` applies them. This invocation "
            f"takes the EW/profile-fit route, which would ignore them silently and emit "
            f"an UNCONDITIONED product. Add --force-synthesis, or drop the flag "
            f"(RYA-1006).")
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

    # RYA-933/934 -- the HOLDING belongs in the stem. Two holdings of one instrument
    # differing by whether tellurics were removed would otherwise write the same
    # filename, and the second would silently overwrite the first.
    stem = f"{a.element}{a.ion}_{int(a.lo)}_{int(a.hi)}_{a.instrument}"
    if a.holding:
        stem += f"_{a.holding}"
    stem += "_PROFILEFIT"
    # ⚠️ ORDER MATTERS. `src` is derived from the stem, and the EW ARTIFACT carries no
    # tier -- stage 1 measures every line it can and knows nothing about grading. The
    # input path resolves FIRST; the tier is appended afterwards, for the OUTPUT only.
    src = EW_DIR / f"{stem}_ew.csv"
    if not src.exists():
        raise SystemExit(f"no measured EWs at {src}. Run measure_band_profilefit.py first.")
    # The tier is part of what was MEASURED, so it belongs in the product name (RYA-984):
    # without it a graded and an ungraded EW product write the SAME filename and the
    # second silently overwrites the first -- the RYA-1006 collision, on this route.
    # ⚠️ UNDERSCORE, not hyphen. The tracker's stem parser reads the handler and then
    # ONE OR MORE selector segments each introduced by `_`
    # (`(?:_[A-Z][A-Z0-9]*(?:-[A-Z]+)?)*`). `PROFILEFIT-GRADED` leaves `-GRADED` with no
    # leading underscore, so parse_stem returns None and the product is INVISIBLE to the
    # page -- caught by test_every_committed_band_product_is_visible_to_the_tracker
    # before merge. The synth route already emits `_SYNTH_GRADED`; this matches it.
    if getattr(a, "lines_tier", "all") != "all":
        stem += f"_{a.lines_tier.upper()}"
    ew = pd.read_csv(src)
    pol = resolve_band(0.5 * (a.lo + a.hi))
    ok = ew[ew.in_aggregate].copy()
    print(f"{a.element} {a.ion}  {a.lo:.0f}-{a.hi:.0f} A  band={pol.name}  "
          f"{len(ok)} in-aggregate of {len(ew)} measured")

    # 🔴 RYA-1031: `--lines-tier` WAS A SILENT NO-OP ON THIS ROUTE. It was consulted only
    # inside `_cand_from_ew_artifact`, which drives the SYNTHESIS route over an EW
    # artifact's lines -- the EW/profile-fit route never looked at it and took every
    # in-aggregate line. Measured on HARPS: `--lines-tier graded` produced a 247-line
    # MIXED pool, 7 of them GF-LAB, rung 1, carrying the 0.17 dex ungraded placeholder --
    # and reported it as a graded run. Same shape RYA-933 fixed on the synthesis side,
    # still open here, and the reason the EW arm could not be compared with the graded
    # synth arms: they were not the same pool.
    _tier = getattr(a, "lines_tier", "all")
    if _tier != "all":
        _g = _graded_mask(ok.wavelength_air_A.astype(float).values)
        _keep = _g if _tier == "graded" else ~_g
        if not _keep.any():
            raise SystemExit(
                f"--lines-tier {_tier}: none of the {len(ok)} in-aggregate lines in "
                f"{src.name} carry a primary laboratory gf. Refusing to emit a "
                f"'{_tier}' product with no {_tier} line in it -- the pool IS the claim, "
                f"and an empty one is a different measurement, not a smaller one.")
        ok = ok[_keep].reset_index(drop=True)
        print(f"    [tier] {_tier}: {len(ok)} lines kept (RYA-946 — graded and ungraded "
              f"are separate products, never merged)")

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
    ctx = build_context(a.element, a.ion, 500000.0, star=a.star)
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
    print("\n[2] ENGINE-A — registered per-line departure corrections...")
    used = [l for l in rows if l.in_aggregate and l.abundance is not None]
    deltas = engine_a_delta(a.element, a.ion,
                            np.array([l.wavelength_air_A for l in used]),
                            cache=a.mpia_cache,
                            star=a.star)
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
                             nlte_source=(engine_a_source(a.element) if d is not None
                                          else "registered per-line delta_nlte (NOT SERVED)"),
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
            la.excluded_reason = ("ENGINE-A-NOT-SERVED: the registered source returns no usable delta_nlte "
                                 "for this line (absent, nan, or a placeholder zero). "
                                 "Reduced coverage, not a failed correction.")
        rows_a.append(_stamp(la))     # RYA-807
    # A departure term added to the SAME inversion — the same profile-fit EWs, so the
    # same handler and the same measured residual (RYA-869).
    p_a = build_product(a.element, a.ion, a.instrument, pol.name, "ENGINE-A", rows_a,
                        handler="ProfileFitHandler",
                        provenance=engine_a_source(a.element))
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
        # RYA-1040 — the <3D> route is a THIRD shape here, not a variant of the other two.
        # `mean3d` says which ATMOSPHERE and which DECK; `nlte` says whether the departures
        # are applied. The LTE member of the pair is mean3d=True, nlte=False -- the same
        # deck's setup with the departures off -- which is the whole reason it is a valid
        # comparand.
        mean3d = a.engine_b_deck in ("gerber-mean3d", "gerber-mean3d-lte")
        nlte = a.engine_b_deck in ("gerber-nlte", "gerber-mean3d")
        if mean3d:
            treatment = (taxes.MEAN3D_NLTE_STAGGER if nlte
                         else taxes.MEAN3D_LTE_STAGGER).token
        else:
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
        # ── RYA-1040: the <3D> leg runs on the <3D> atmosphere, not on MARCS ─────────
        # 🔴 AND THE DECK KEY IS SELECTED, NOT GUESSED FROM THE ELEMENT. The line below
        # used to read `gnlte.for_node("Fe" if a.element == "Fe" else a.element, ...)` --
        # an identity expression that resolves the deck from the ELEMENT alone, so it
        # picks `Fe` and can never reach `Fe@mean3D` no matter what deck was asked for.
        # The deck-provenance axis Ryan ratified is what fixes it: the CLI names a deck,
        # the key is built from it, and nothing infers a deck from an element again.
        #
        # ⚠️ The LTE MEMBER TAKES THIS BRANCH TOO (mean3d and not nlte). It must run on
        # the SAME atmosphere as its NLTE partner or the pair stops being a comparand --
        # which is the entire reason the pair exists (RYA-542).
        if mean3d:
            from pipeline import gerber_nlte as gnlte
            from pipeline import mean3d_atmosphere as m3d
            deck_key = f"{a.element}@mean3D"
            if deck_key not in gnlte.DECKS:
                raise SystemExit(
                    f"no <3D> deck registered for {a.element} (looked for {deck_key!r}). "
                    f"Registered: {sorted(gnlte.DECKS)}. Staging is RYA-710; the aux file "
                    f"to register against is in data/results/rya1035/"
                    f"mean3d_aux_defect_sweep.csv -- Fe and Mn REQUIRE the plain aux.")
            layers, model3d = m3d.load(gnlte.deck_atmosphere(deck_key))
            ctx_b["atmosphere"] = layers
            # iSpec must READ the real structure from the mul23 file: the array above
            # carries only the three fields iSpec reads off it (see mean3d_atmosphere).
            ctx_b["atmosphere_layers_file"] = model3d["path"]
            ctx_b["nlte_deck"] = "gerber" if nlte else None
            ctx_b["nlte_deck_key"] = deck_key if nlte else None
            dep = gnlte.for_node(deck_key, float(ctx["teff"]), float(ctx["logg"]),
                                 float(ctx["feh"]))
            # BOTH gates: the depths must pair index-for-index, and the two tau scales
            # must actually be the same scale -- iSpec overwrites the departure tau with
            # the atmosphere's, so a disagreement would be applied silently.
            gnlte.assert_depth_match(dep, layers)
            dtau = m3d.assert_tau_consistent(dep, model3d)
            print(f"    <3D> deck {deck_key} on {Path(model3d['path']).name} "
                  f"(ndep={dep['ndep']}, nlev={dep['nk']}, "
                  f"max|dlogtau|={dtau:.2e}, departures={'ON' if nlte else 'OFF'})")
            _b_source = MEAN3D_SOURCE_FMT.format(
                atom=dep['atom_path'].split('/')[-1],
                model=Path(model3d['path']).name,
                mode=("NLTE, departures applied" if nlte
                      else "LTE, departures WITHHELD -- the paired comparand"))
            if nlte:
                # ⚠️ iSpec fails SOFT: an element whose linelist rows carry no NLTE label
                # lands in `nlte_ignored` and the synthesis proceeds in LTE WITHOUT
                # raising (RYA-764). That check belongs ONLY to the NLTE member -- running
                # it on the LTE member would demand NLTE labels of a run that deliberately
                # applies none, and would refuse the comparand for lacking what it does
                # not use.
                n_lab = gnlte.assert_linelist_supports_nlte(
                    ctx["linelist"], int(ctx["atom_code"]), a.element)
                print(f"    {n_lab} NLTE-labelled {a.element} lines in the list")
            else:
                # The comparand: same deck, same atmosphere, departures withheld. Passing
                # no departures is what makes it the LTE limit OF THIS SETUP rather than
                # an unrelated LTE run on a different atmosphere.
                dep = None
        elif nlte:
            from pipeline.abundances_derive import _load_atmosphere
            from pipeline import gerber_nlte as gnlte
            ctx_b["atmosphere"] = _load_atmosphere(
                float(ctx["teff"]), float(ctx["logg"]), float(ctx["feh"]),
                float(ctx["vturb"]), model_grid="MARCS.GES")
            ctx_b["nlte_deck"] = "gerber"
            dep = gnlte.for_node(a.element,
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
        # 🔴 RYA-959 — `holding=a.holding` WAS MISSING HERE, AND ON THE TWO SITES BELOW.
        # RYA-933/934 added `--holding` so a caller can name ONE product instead of
        # taking the instrument's first covering candidate; the synthesis route honours
        # it (see `synthesis_route`) and this leg did not. Measured: a run launched with
        # `--holding solar_harps_molecfit_corrected` printed
        # `[holding] harps -> solar_harps` and derived Engine-B on the UNCORRECTED
        # sibling — the one pair this repo must never collapse (RYA-911/927).
        _b_probe = select_holding(a.instrument, 0.5 * (a.lo + a.hi), 1.4,
                                  holding=a.holding)
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
                _win = load_window_ex(a.instrument, c, pad=1.4, segs=segs,
                                      holding=a.holding)
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
                       [select_holding(a.instrument, 0.5 * (a.lo + a.hi), 1.4,
                                       holding=a.holding)]}
        _b_unexpected = _b_holdings - _b_expected
        if _b_unexpected:
            raise SystemExit(
                f"PROVENANCE MISMATCH: this product is tagged instrument "
                f"{a.instrument!r} but its lines were measured on holding(s) "
                f"{sorted(_b_unexpected)}, which {a.instrument!r} does not select. "
                f"A number labelled with an instrument it was not measured on is the "
                f"RYA-911/913 defect. Refusing to emit it.")
        # 🔴 RYA-959 — AND THE GUARD COULD NOT HAVE CAUGHT IT. `_b_expected` above is
        # resolved by the same call the lines were loaded with, so the check proved the
        # leg was self-consistent, never that it honoured what the CALLER asked for. When
        # `--holding` is given, the answer is not "whatever the instrument selects" — it
        # is that holding, and nothing else may serve a single line.
        if a.holding and _b_holdings and _b_holdings != {a.holding}:
            raise SystemExit(
                f"HOLDING IGNORED: --holding {a.holding!r} was requested but Engine-B's "
                f"lines were measured on {sorted(_b_holdings)}. A run asked for one "
                f"product and got another; for the HARPS pair that is the difference "
                f"between the telluric-corrected sibling and the raw one (RYA-931/927). "
                f"Refusing to emit it.")
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
                                       handler=prod.handler or None,
                                       model=(engine_a_model(a.element)
                                              if prod.treatment == "ENGINE-A" else None),
                                       atmos=("marcs-ges"
                                              if prod.treatment == "ENGINE-B-NLTE"
                                              else "atlas9")).as_columns(),
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
