#!/usr/bin/env python3
"""RYA-1230 -- assemble the RYA-587 budget for every re-run C/N/O band product.

    python3 scripts/rya1230_cno_budget.py --legs DIR --nominal N_DIR CO_DIR --stage DIR

Ryan, 2026-09-26: a re-run owes its uncertainty, not just its value. Template: RYA-1226's
`rya1226_migrate_nearuv_budget.build`, the one path that has carried products through the
contract. Every term below is MEASURED on the product's OWN pool by a paired leg
(`scripts/rya1230_cno_budget_legs.py`), or is the sourced value RYA-1226 already carried,
with its evidence and its flag. Anything the evidence cannot price is left out and the
CONTRACT holds it -- `validate()` issues the verdict, not this script.

  measurement         std(ddof=1)/sqrt(N) of the accepted pool (line_scatter)
  transition_data     canonical_gf per-line sigma (NIST class for C/N/O), joined on lambda
                      AND EP (RYA-1037); lines sharing one source are fully correlated
  stellar.xi          central paired dA/dxi at 0.90/1.10 km/s, stamped legs (RYA-1178 A),
                      x RYA-1089's delta_xi = 0.2912 km/s
  stellar.teff/logg/  DEFINED zero, the RYA-1226 solar convention (pinned; Teff 1 K)
  metallicity
  continuum           |paired median(A[q97] - A[q80])| / 2 -- the model-guided rule's
                      pixel-selection quantile either side of the nominal q90
  profile_ew          |paired median(A[+/-0.25 A core] - A[nominal window])|
  model_atmosphere    |paired median(A[MARCS.GES] - A[ATLAS9.Castelli])|
  telluric            per line: measured sky depth (raw / independent-correction pair) ->
                      measured residual of THIS holding against an independent correction
                      on the telluric pixels -> dex through the line's own dA/df (the
                      +0.1% continuum leg). 0 where the sky has no pixel above its clean edge.
  blends              N I only: the CN inside the profile follows A(C); central paired
                      response to A(C) +/- 0.10 scaled to AGSS21's sigma(C) = 0.04
  nlte                N/A on a 1D-LTE leg; on a departure leg ASPLUND+2021's rule: half the
                      pool's own non-LTE correction, floor 0.03 (RYA-1232; was RYA-1032's
                      cross-element 0.043)
  model_atmosphere    (+) on a 3D leg, half the pool's own 3D effect (Asplund+2021's
                      inhomogeneity term), in quadrature
  pseudo_continuum    N/A: the continuum is observed in these bands and placed per line
  hfs_isotopes        N/A where every pool line has hfs_n_components == 1
  molecular_coupling  N/A: atomic selector
  holding_instrument  MEASURED 0: SynthesisHandler harness residual (RYA-869)
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import re
import shlex
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from pipeline.uncertainty_contract import (  # noqa: E402
    UncertaintyError, assemble, pool_digest, product_scope, transition_data, validate)
from pipeline.paired_differential import paired_differential  # noqa: E402

GF = ROOT / "data/linelists/canonical_gf.csv"
DELTA_XI = 0.2912
XI_STEP = 0.10
C_STEP, SIGMA_C = 0.10, 0.04
SIGMA_C_SOURCE = "Asplund, Amarsi & Grevesse 2021 (A&A 653, A141) Table 2: A(C) = 8.46 +/- 0.04"
CSCALE = 0.001
SPECIES = {"C": "C I", "N": "N I", "O": "O I", "Fe": "Fe I"}
#: RYA-1232: Fe II units carry --ion II; the canonical_gf species is ion-specific.
SPECIES_ION = {("Fe", "II"): "Fe II"}


def _selector(args: list[str]) -> str | None:
    if "--lines-from-set" in args:
        return "SET-" + args[args.index("--lines-from-set") + 1].split("=")[0]
    return None


def _arg(args, flag, default=None):
    return args[args.index(flag) + 1] if flag in args else default


def _acc(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["in_aggregate"] == True]  # noqa: E712


def _paired(leg: pd.DataFrame, nom: pd.DataFrame, n_product: int) -> dict:
    r = paired_differential(leg, nom)
    return {"median": float(r.median), "n_paired": int(r.n_paired),
            "moved": int(r.n_paired) != n_product}


def _leg_lines(legdir: Path, stem: str) -> pd.DataFrame | None:
    """The leg's per-line file for this product. The stem carries the EFFECTIVE range
    (`lo_hi`), which a leg can clip differently (the +/-0.25 A core leg reads 4200_6909
    where the nominal reads 4200_6908), so the match ignores those two numbers and must be
    unique."""
    f = legdir / stem
    if f.exists():
        return pd.read_csv(f)
    pat = re.sub(r"^([A-Z][a-z]?I+)_\d+_\d+_", r"\1_*_*_", stem)
    hits = [h for h in glob.glob(str(legdir / pat))
            if re.sub(r"^([A-Z][a-z]?I+)_\d+_\d+_", "", Path(h).name)
            == re.sub(r"^([A-Z][a-z]?I+)_\d+_\d+_", "", stem)]
    return pd.read_csv(hits[0]) if len(hits) == 1 else None


# ── telluric ────────────────────────────────────────────────────────────────────
_SKY = {"kpno_solar_atlas": ("kpno_solar_atlas", "solar_kpno", "solar_kpno_kurucz2005_corrected"),
        "iag_fts_solar_atlas": ("kpno_solar_atlas", "solar_kpno", "solar_kpno_kurucz2005_corrected"),
        #: HARPS raw cannot set a clean edge (continuum S/N ~31 < the 200 science floor,
        #: telluric_observability refuses a verdict), so its sky MAP is Kitt Peak's too --
        #: line positions are the atmosphere's; the residual is still HARPS's own.
        "harps": ("kpno_solar_atlas", "solar_kpno", "solar_kpno_kurucz2005_corrected"),
        "crires_plus": ("kpno_solar_atlas", "solar_kpno", "solar_kpno_kurucz2005_corrected")}
#: INDEPENDENT corrections to difference each holding against on the telluric pixels, in
#: order. A candidate byte-identical to raw in THIS window is not a correction there (KP
#: molecfit leaves the 9000-9300 A H2O band untouched, RYA-1190) and is skipped -- using it
#: measured the whole sky as kurucz2005's "residual" on C I 9061 (0.33 flux, 5 dex).
_INDEPENDENT = {"solar_kpno_molecfit_corrected": [("kpno_solar_atlas", "solar_kpno_kurucz2005_corrected"),
                                                  ("iag_fts_solar_atlas", "solar_iag"),
                                                  ("crires_plus", "solar_crires_plus_j_rya1219")],
                "solar_crires_plus_j_rya1219": [("kpno_solar_atlas", "solar_kpno_molecfit_corrected"),
                                                ("iag_fts_solar_atlas", "solar_iag")],
                "solar_kpno_kurucz2005_corrected": [("kpno_solar_atlas", "solar_kpno_molecfit_corrected"),
                                                    ("iag_fts_solar_atlas", "solar_iag")],
                "solar_iag": [("kpno_solar_atlas", "solar_kpno_kurucz2005_corrected"),
                              ("kpno_solar_atlas", "solar_kpno_molecfit_corrected"),
                              ("crires_plus", "solar_crires_plus_j_rya1219")],
                "solar_harps_molecfit_corrected": [("kpno_solar_atlas", "solar_kpno_kurucz2005_corrected"),
                                                   ("iag_fts_solar_atlas", "solar_iag")],
                #: RYA-1232: the CRIRES+ Y arms (9802-10794 A) were missing -> KeyError held 4 Fe
                #: NIR products. Independent corrections that cover Y: KP molecfit (to 13000),
                #: IAG (to 11086), Kurucz 2005 (to 10008).
                "solar_crires_plus_y_wide_rya1054": [("kpno_solar_atlas", "solar_kpno_molecfit_corrected"),
                                                     ("iag_fts_solar_atlas", "solar_iag"),
                                                     ("kpno_solar_atlas", "solar_kpno_kurucz2005_corrected")],
                "solar_crires_plus_y_rya794": [("kpno_solar_atlas", "solar_kpno_molecfit_corrected"),
                                               ("iag_fts_solar_atlas", "solar_iag"),
                                               ("kpno_solar_atlas", "solar_kpno_kurucz2005_corrected")],
                #: CRIRES+ H has NO independent correction on the Mac (KP/IAG stop short of H);
                #: its residual needs the raw/corrected CRIRES+ pair on Sirius (RYA-1192).
                "solar_crires_plus_h_rya1094": [],
                #: RYA-1232: our full-arm corrected H. Independent correction = Elgueta+2026's
                #: own reduction of the same night (a different telluric removal).
                "solar_crires_plus_h_rya1232": [("crires_plus", "solar_crires_plus_h_rya1094")]}
#: RYA-1232 -- holdings that carry their OWN molecfit transmission per pixel (`min_mtrans`
#: in the rest-frame CSV). KP's raw/corrected sky pair stops at 13000 A, so no H-band product
#: ever had its telluric component resolved; for these the sky is the model the correction
#: divided by. No reduction-difference null applies (it is a model, not a ratio of two
#: reductions); a pixel is telluric where the model absorbs more than OWN_SKY_EDGE.
_OWN_SKY = {"solar_crires_plus_h_rya1232":
            "data/results/rya1214_crires_jk/solar_crires_plus_h_rya1232_rest.csv"}
OWN_SKY_EDGE = 0.01
_EDGE: dict = {}
_SPAN: dict = {}


def _full_span(H, inst, hold):
    if (inst, hold) not in _SPAN:
        #: a holding can be registered more than once (solar_crires_plus_j_rya1219 has a
        #: RYA-1219 molecfit reader and an Elgueta reader); use the first span that loads
        last = None
        for spec in (sp for specs in H._INSTRUMENT_HOLDINGS.values() for sp in specs
                     if sp.holding_id == hold and sp.span_A):
            lo, hi = spec.span_A
            try:
                _SPAN[(inst, hold)] = H.load_window_ex(inst, 0.5 * (lo + hi), 0.5 * (hi - lo),
                                                       holding=hold, allow_uncorrected=True)
                break
            except LookupError as exc:
                last = exc
        else:
            raise last
    return _SPAN[(inst, hold)]


def _telluric_own_sky(instrument, holding, w0, hw) -> dict:
    import measure_band_ew as H
    if holding not in _OWN_SKY_CACHE:
        _OWN_SKY_CACHE[holding] = pd.read_csv(ROOT / _OWN_SKY[holding])
    d = _OWN_SKY_CACHE[holding]
    g = np.linspace(w0 - hw, w0 + hw, 400)
    t = np.interp(g, d["wavelength_air_A"], d["min_mtrans"])
    P = (1.0 - t) > OWN_SKY_EDGE
    ev = {"sky_pair": f"{holding} min_mtrans (molecfit model)", "clean_edge": OWN_SKY_EDGE,
          "clean_window_null": 0.0, "sky_baseline": 1.0, "n_px": int(g.size),
          "n_telluric_px": int(P.sum()), "max_depth": float(max(0.0, 1.0 - t.min()))}
    if not P.any():
        ev["residual_flux"] = 0.0
        ev["sky_absorption"] = 0.0
        return ev
    ev["sky_absorption"] = float(np.mean(np.clip(1.0 - t[P], 0, None)) * P.mean())
    own = np.interp(g, d["wavelength_air_A"], d["flux_normalized"])
    for ind_inst, ind in _INDEPENDENT[holding]:
        try:
            x = H.load_window_ex(ind_inst, w0, hw + 0.5, holding=ind, allow_uncorrected=True)
        except LookupError:
            continue
        ref = np.interp(g, x.wave, x.flux)
        ratio = own / ref
        rbase = float(np.median(ratio[~P])) if (~P).sum() >= 3 else float(np.median(ratio))
        ev["residual_basis"] = f"{holding} vs independent correction {ind}, baseline-removed"
        ev["residual_flux"] = float(np.mean(np.abs(ratio[P] - rbase)) * P.mean())
        return ev
    ev["residual_flux"] = None
    ev["residual_basis"] = f"NO independent correction covers {w0:.3f} A for {holding}"
    return ev


_OWN_SKY_CACHE: dict = {}


def telluric_line(instrument, holding, w0, hw, band_lo, band_hi) -> dict:
    import measure_band_ew as H
    from pipeline import telluric_observability as T
    if holding in _OWN_SKY:
        return _telluric_own_sky(instrument, holding, w0, hw)
    sky_inst, raw, cor = _SKY.get(instrument, _SKY["kpno_solar_atlas"])
    if w0 > 10000.0:
        #: Kurucz 2005 stops at 10000 A; beyond it the sky is raw KP over the RYA-1230
        #: full-coverage molecfit holding (0 A raw, 505bdf0f).
        cor = "solar_kpno_molecfit_corrected"
    key = (raw, cor, band_lo, band_hi)
    if key not in _EDGE:
        snr, _ = T.band_continuum_snr(H, sky_inst, raw, band_lo, band_hi, max(hw, 2.0))
        _EDGE[key] = T.thresholds(snr).clean_max_depth
    if ("null", raw, cor) not in _EDGE:
        #: 🔴 THE CLEAN-WINDOW NULL. raw / corrected differs by REDUCTION everywhere, not
        #: only where the sky absorbs (C I 4269 read "telluric" pixels in the blue). Measure
        #: that pair's own dip statistic where no molecular telluric line exists --
        #: 4300-4900 A, blueward of every O2/H2O band -- and require a telluric pixel to
        #: clear it (RYA-1192: "a template needs a clean-window null").
        dips = []
        for c0 in (4300.0, 4450.0, 4600.0, 4750.0, 4900.0):  # null measured on THIS pair
            gg = np.linspace(c0 - 1.0, c0 + 1.0, 400)
            a_ = H.load_window_ex(sky_inst, c0, 1.5, holding=raw, allow_uncorrected=True)
            b_ = H.load_window_ex(sky_inst, c0, 1.5, holding=cor, allow_uncorrected=True)
            tt = np.interp(gg, a_.wave, a_.flux) / np.interp(gg, b_.wave, b_.flux)
            dips.append(float(np.percentile(np.median(tt) - tt, 99.5)))
        _EDGE[("null", raw, cor)] = max(dips)
    edge = max(_EDGE[key], _EDGE[("null", raw, cor)])
    g = np.linspace(w0 - hw, w0 + hw, 400)

    def flux(inst, hold):
        try:
            x = H.load_window_ex(inst, w0, hw + 0.5, holding=hold, allow_uncorrected=True)
        except LookupError:
            #: CRIRES+ readers validate frame-control lines that only exist across the full
            #: arm, so a few-A window cannot pass; read the holding's full span once (the
            #: SAME flux its own fit read) and slice it
            x = _full_span(H, inst, hold)
            if not (x.wave.min() <= g.min() and g.max() <= x.wave.max()):
                raise
        return np.interp(g, x.wave, x.flux)

    t = flux(sky_inst, raw) / flux(sky_inst, cor)
    base = float(np.median(t))
    P = (base - t) > edge
    ev = {"sky_pair": f"{raw}/{cor}", "clean_edge": edge,
          "clean_window_null": _EDGE[("null", raw, cor)], "sky_baseline": base,
          "n_px": int(g.size), "n_telluric_px": int(P.sum()),
          "max_depth": float(max(0.0, base - t.min()))}
    if not P.any():
        ev["residual_flux"] = 0.0
        ev["sky_absorption"] = 0.0
        return ev
    #: the window's mean sky absorption (fraction-weighted), so a residual-to-absorption
    #: ratio measured on covered windows can bound an uncovered one (RYA-1230 CN)
    ev["sky_absorption"] = float(np.mean(np.clip(base - t[P], 0, None)) * P.mean())
    try:
        own = flux(instrument, holding)
    except LookupError as exc:
        #: the sky is measured (KP pair) but this holding's window loader has a gap here
        ev["residual_flux"] = None
        ev["residual_basis"] = f"holding loader cannot serve this window: {str(exc)[:100]}"
        return ev
    if np.array_equal(own, flux(sky_inst, raw)) if instrument == sky_inst else False:
        #: this holding is the raw flux here: nothing was corrected, the residual IS the sky
        ev["residual_basis"] = "holding byte-identical to raw in this window: full depth"
        ev["residual_flux"] = float(np.mean(np.clip(base - t[P], 0, None)) * P.mean())
        return ev
    raw_f = flux(sky_inst, raw)
    ref = None
    for ind_inst, ind in _INDEPENDENT[holding]:
        try:
            cand = flux(ind_inst, ind)
        except Exception:                                    # noqa: BLE001 -- not staged / no coverage
            continue
        if ind_inst == sky_inst and np.array_equal(cand, raw_f):
            continue
        ref = cand
        break
    if ref is None:
        ev["residual_flux"] = None
        ev["residual_basis"] = f"NO independent correction covers {w0:.3f} A for {holding}"
        return ev
    ratio = own / ref
    rbase = float(np.median(ratio[~P])) if (~P).sum() >= 3 else float(np.median(ratio))
    ev["residual_basis"] = f"{holding} vs independent correction {ind}, baseline-removed"
    ev["residual_flux"] = float(np.mean(np.abs(ratio[P] - rbase)) * P.mean())
    return ev


# ── one product ────────────────────────────────────────────────────────────────
def build(unit_args, stem, nominal_dir: Path, unit_dir: Path, gf: pd.DataFrame):
    from config.synth_bands import SYNTH_BANDS
    from pipeline.band_policy import resolve as band_of
    from publish_product import normalise

    lines_stem = stem
    prod_stem = stem.replace("_1D-LTE_lines.csv", "_products.csv").replace(
        "_lines.csv", "_products.csv")
    nom = pd.read_csv(nominal_dir / lines_stem)
    prod = pd.read_csv(nominal_dir / prod_stem)
    acc = _acc(nom)
    n = len(acc)
    # 🔴 RYA-1232 -- LEG-UNSTABLE LINES. A line the nominal fit accepts but a budget leg's
    # fit REJECTS (NON-MINIMUM at frac_rise ~1e-5, edge_pinned, FIT-NOT-PHYSICAL) is not
    # robustly measured in the nominal either: its acceptance hinges on a trivial
    # perturbation. One such line out of 176 held 16 Fe products ("pool moved"). It is
    # excluded from the pool as LEG-UNSTABLE -- recorded with the leg and the reason -- and
    # the product re-aggregated the way the product is formed (median of the accepted
    # lines' `abundance`, stat = std/sqrt(n)). Bounded: at most 5% of the pool and >= 2
    # lines left, otherwise the moved pool still HOLDS (it is then not a marginal line).
    unstable = {}
    for _lg in ("xi_minus", "xi_plus", "core", "q80", "q97", "contref", "marcs", "cscale",
                "c_minus", "c_plus"):
        _L = _leg_lines(unit_dir / _lg, lines_stem)
        if _L is None:
            continue
        _ok = set(_acc(_L)["wavelength_air_A"].round(3))
        _all = {round(float(w), 3): r for w, r in zip(_L["wavelength_air_A"],
                                                      _L.get("excluded_reason", pd.Series([""] * len(_L))))}
        for _w in set(acc["wavelength_air_A"].round(3)) - _ok:
            unstable.setdefault(float(_w), []).append(f"{_lg}: {str(_all.get(_w, ''))[:90]}")
    if unstable and len(unstable) <= max(1, int(0.05 * n)) and n - len(unstable) >= 2:
        acc = acc[~acc["wavelength_air_A"].round(3).isin(list(unstable))]
        n = len(acc)
        _ab = acc["abundance"].astype(float)
        restat = {"A": round(float(_ab.median()), 3),
                  "stat_dex": round(float(_ab.std(ddof=1) / math.sqrt(n)), 4),
                  "n_lines": n, "leg_unstable": {f"{k:.3f}": v for k, v in unstable.items()}}
    else:
        restat = None
    _pool_w = set(acc["wavelength_air_A"].round(3))

    def _accp(df):
        a_ = _acc(df)
        return a_[a_["wavelength_air_A"].round(3).isin(_pool_w)]
    element = str(prod["element"].iloc[0])
    holding = _arg(unit_args, "--holding")
    instrument = _arg(unit_args, "--instrument")
    selector = _selector(unit_args)
    _rows = normalise(prod, holding=holding, tier="ALL", route="SYNTH", selector=selector)
    if not _rows:
        #: an EMPTY product (no accepted line, no value) -- recorded, never a crash
        return {"row": {"A": None, "treatment": str(prod.get("treatment", pd.Series([None])).iloc[0]),
                        "holding": holding, "selector": selector, "element": str(prod["element"].iloc[0]),
                        "band": None},
                "skip": "empty product: no accepted line, no value to budget"}
    row = _rows[0]
    row["star"] = "solar"
    if restat:
        row["A"], row["n_lines"] = restat["A"], restat["n_lines"]
        if "sigma_stat" in row:
            row["sigma_stat"] = restat["stat_dex"]
    band = row["band"]
    hw = float(_arg(unit_args, "--half-width-A", SYNTH_BANDS[band].half_width_A))
    lo, hi = float(_arg(unit_args, "--lo")), float(_arg(unit_args, "--hi"))
    notes, comps = [], []
    # RYA-1232 (Ryan: "we do what other scientists do") -- a SINGLE-line pool is priced the
    # way Asplund+2021 Sect. 2.1 prices one or two lines: from the goodness of the profile
    # fit. That is the line's own chi2-curvature sigma (fit_constraint.curvature_sigma,
    # rescaled to red_chi2 = 1, worse side), recorded on the 1D-LTE fit; a departure leg
    # adds a per-line delta to that SAME fit, so it carries the same sigma.
    single_sigma = None
    if n == 1:
        sa = pd.to_numeric(acc["sigma_A"], errors="coerce").iloc[0] if "sigma_A" in acc else float("nan")
        if not np.isfinite(sa):
            lte = nominal_dir / re.sub(r"_(ENGINE-A-3DNLTE|ENGINE-A|ENGINE-B[^_]*)_lines\.csv$",
                                       "_1D-LTE_lines.csv", lines_stem)
            if lte.exists():
                l1 = pd.read_csv(lte)
                w0 = float(acc["wavelength_air_A"].iloc[0])
                m = l1[(l1.wavelength_air_A - w0).abs() < 1e-3]
                if len(m):
                    sa = pd.to_numeric(m["sigma_A"], errors="coerce").iloc[0]
        if not np.isfinite(sa):
            return {"row": row, "skip": "n_lines=1 and the line's fit records no curvature "
                                        "sigma -- nothing to price it from; HOLD"}
        single_sigma = float(sa)
    elif n < 1:
        return {"row": row, "skip": "no accepted line"}

    # transition data: lambda AND EP join
    g = gf[gf["species"] == SPECIES_ION.get((element, str(_arg(unit_args, "--ion", "I"))), SPECIES[element])]
    ids, sig, src = [], [], []
    for _, l in acc.iterrows():
        m = g[((g.wavelength_air_A - l.wavelength_air_A).abs() < 0.01)
              & ((g.excitation_potential_eV - l.ep_eV).abs() < 0.005)]
        if len(m) != 1:
            return {"row": row, "skip": f"{l.wavelength_air_A}: {len(m)} canonical_gf rows "
                                        f"on lambda+EP -- refusing a lambda-only join (RYA-1037)"}
        m = m.iloc[0]
        ids.append(str(m.physical_id))
        s = float(m.gf_sigma_dex)
        ref = str(m.gf_source_doi if pd.notna(m.gf_source_doi) and str(m.gf_source_doi).strip()
                  else m.loggf_reference)
        if math.isnan(s):
            #: No stored sigma: the row's NIST accuracy class through the SSOT bridge
            #: (gf_grades.nist_sigma_dex). Where `nist_grade` and `gf_tier` name DIFFERENT
            #: classes (C I 5380.325: B vs NIST-C+), charge the WORSE one and say so.
            from pipeline.gf_grades import nist_sigma_dex
            classes = {str(m.nist_grade).strip()} | (
                {str(m.gf_tier).replace("NIST-", "").strip()}
                if str(m.gf_tier).startswith("NIST-") else set())
            classes = {c for c in classes if not math.isnan(nist_sigma_dex(c))}
            if classes:
                worst = max(classes, key=nist_sigma_dex)
                s = nist_sigma_dex(worst)
                ref = (f"{ref}; sigma from NIST class {worst} via gf_grades.nist_sigma_dex"
                       + (f" (row names classes {sorted(classes)}; the worse is charged)"
                          if len(classes) > 1 else ""))
        sig.append(s)
        src.append(ref)
    if any(math.isnan(s) for s in sig):
        return {"row": row, "skip": "a pool line carries no published gf sigma"}
    digest = pool_digest(ids)
    if single_sigma is not None:
        sys.path.insert(0, str(ROOT / "scripts"))
        from measure_band_ew import load_window_ex
        w0 = float(acc["wavelength_air_A"].iloc[0])
        npx = int(np.isfinite(np.asarray(load_window_ex(instrument, w0, hw, holding=holding).flux)).sum())
        from pipeline.uncertainty_contract import fit_measurement
        comps.append(fit_measurement(
            single_sigma,
            source=("Asplund+2021 Sect. 2.1: one line -> statistical error from the goodness of "
                    "the profile fit; chi2-curvature sigma of this line, rescaled to red_chi2 = 1"),
            likelihood="synthesis chi2 over the fit window, curvature at the best abundance",
            correlation_treatment=("red_chi2 rescaling absorbs pixel correlation and model "
                                   "inadequacy; the worse-constrained side is taken"),
            n_pixels=npx))
    else:
        raw = float(acc["abundance"].astype(float).std(ddof=1))
        comps.append(dict(name="measurement", sigma_dex=raw / math.sqrt(n), state="MEASURED",
                          source="per-line scatter of the accepted pool (RYA-1230 re-run artifact)",
                          evidence={"method": "line_scatter", "independent": True, "n_lines": n,
                                    "raw_sigma": raw, "pool_sha256": digest}))
    w = [1.0 / n] * n
    cov = [[sig[a] * sig[b] if src[a] == src[b] else 0.0 for b in range(n)] for a in range(n)]
    comps.append(transition_data(ids, sig, w, covariance=cov, sources=src,
                                 covariance_source=("canonical_gf per-line sigma; lines sharing "
                                                    "one source are fully correlated")))
    for name, why in (("stellar.logg", "solar log g is pinned; delta_logg is definitionally zero (RYA-1089)"),
                      ("stellar.metallicity", "solar [Fe/H] is pinned; delta_feh is definitionally zero (RYA-1089)")):
        comps.append(dict(name=name, sigma_dex=0.0, state="DEFINED", source=why,
                          evidence={"delta": 0.0, "parameter_exists": True}))
    comps.append(dict(name="stellar.teff", sigma_dex=0.0, state="DEFINED",
                      source=("sourced DEFINED zero (RYA-1226 convention): the solar Teff "
                              "allowance is 1 K"),
                      evidence={"parameter_exists": True, "delta_K": 1.0,
                                "owed_measurement": "a Teff perturbation on this pool"}))

    # xi
    from rya1120_xi_campaign import dA_dxi
    lo_d, hi_d = unit_dir / "xi_minus", unit_dir / "xi_plus"
    L, Hh = _leg_lines(lo_d, lines_stem), _leg_lines(hi_d, lines_stem)

    def _restrict(df):
        #: the xi legs are paired leg-vs-leg, so they must see the SAME stable pool as
        #: every other term (a leg-unstable line both xi legs accept still moved it)
        if df is None:
            return None
        df = df.copy()
        df.loc[~df["wavelength_air_A"].round(3).isin(_pool_w), "in_aggregate"] = False
        return df
    L, Hh = _restrict(L), _restrict(Hh)
    if L is not None and Hh is not None:
        d = dA_dxi(Hh, L, step_kms=XI_STEP, minus_dir=lo_d, plus_dir=hi_d, xi_nominal=1.0)
        pl = paired_differential(Hh, L)
        if int(pl.n_paired) == n:
            slope = float(d["dA_dxi_paired"])
            signed = slope * DELTA_XI
            comps.append(dict(name="stellar.xi", sigma_dex=abs(signed), state="MEASURED",
                              source=("per-line paired differential at xi 0.90/1.10 on this "
                                      "product's OWN pool, stamped legs (RYA-1178 A)"),
                              evidence={"pool_sha256": digest,
                                        "parameter_source": "RYA-1089 sourced solar delta_xi = 0.2912 km/s",
                                        "response_assessment": (
                                            f"central paired median over {pl.n_paired} lines; "
                                            f"difference of aggregates {d['dA_dxi_from_aggregates']:.4f}/km/s"),
                                        "signed_response_dex": signed, "delta_plus_dex": signed,
                                        "delta_minus_dex": -signed, "delta_parameter": DELTA_XI,
                                        "dA_dxi": slope, "parameter_exists": True}))
        else:
            notes.append(f"xi pool moved: {pl.n_paired} paired vs n={n}")

    def lever(name, legname, source, extra=None):
        leg = _leg_lines(unit_dir / legname, lines_stem)
        if leg is None:
            notes.append(f"{name}: leg {legname} missing")
            return
        p = _paired(_accp(leg), acc, n)
        if p["moved"]:
            notes.append(f"{name}: pool moved ({p['n_paired']} vs {n})")
            return
        comps.append(dict(name=name, sigma_dex=abs(p["median"]), state="MEASURED",
                          source=source, evidence={"pool_sha256": digest,
                                                   "paired_median_dex": p["median"],
                                                   "n_paired": p["n_paired"], **(extra or {})}))

    # continuum: central difference of the envelope estimator, p99 vs p90 around the
    # nominal p95 -- bins and window unchanged, so the envelope stays constrained
    lo_l, hi_l = _leg_lines(unit_dir / "q80", lines_stem), _leg_lines(unit_dir / "q97", lines_stem)
    if lo_l is None or hi_l is None:
        notes.append("continuum: q80/q97 legs missing")
    else:
        pc = _paired(_accp(hi_l), _accp(lo_l), n)
        if pc["moved"]:
            notes.append(f"continuum: pool moved ({pc['n_paired']} vs {n})")
        else:
            # RYA-1232: + the REFERENCE spread (IAG atlas vs synthesis), in quadrature --
            # Amarsi+2021's two-atlas spread. Outside 5001-11086 A the leg equals nominal.
            ref_l = _leg_lines(unit_dir / "contref", lines_stem)
            pr = _paired(_accp(ref_l), acc, n) if ref_l is not None else None
            if ref_l is None:
                notes.append("continuum: contref leg missing")
            elif pr["moved"]:
                notes.append(f"continuum: contref pool moved ({pr['n_paired']} vs {n})")
            else:
                comps.append(dict(
                    name="continuum", sigma_dex=float(np.hypot(pc["median"] / 2.0, pr["median"])),
                    state="MEASURED",
                    source=("standing continuum rule (synthesis reference, 3-MAD clipped): "
                            "placement = q97 vs q80 central half-difference, REFERENCE = "
                            "IAG-atlas leg minus nominal, in quadrature, paired on this pool"),
                    evidence={"pool_sha256": digest, "paired_median_q97_minus_q80": pc["median"],
                              "paired_median_iag_minus_synthesis": pr["median"],
                              "n_paired": pc["n_paired"]}))
    lever("profile_ew", "core", "fit window +/-0.25 A core vs the band's fixed window, paired on this pool",
          {"varied": "fit half-width", "nominal_A": hw, "alternate_A": 0.25})
    lever("model_atmosphere", "marcs", "MARCS.GES vs ATLAS9.Castelli, one axis varied, paired on this pool",
          {"varied": "atmosphere grid"})

    # telluric: per line, through the line's own dA/df
    cs = _leg_lines(unit_dir / "cscale", lines_stem)
    if cs is None:
        notes.append("telluric: cscale leg missing")
    else:
        m = acc.merge(cs[["wavelength_air_A", "abundance"]], on="wavelength_air_A",
                      suffixes=("", "_cs"))
        if len(m) != n or m["abundance_cs"].isna().any():
            notes.append("telluric: cscale pool moved")
        else:
            per, tot = [], 0.0
            for _, l in m.iterrows():
                dadf = (float(l.abundance_cs) - float(l.abundance)) / (-CSCALE)
                try:
                    ev = telluric_line(instrument, holding, float(l.wavelength_air_A), hw, lo, hi)
                    if ev.get("residual_flux") is None:
                        raise RuntimeError(ev["residual_basis"])
                except Exception as exc:                                  # noqa: BLE001
                    notes.append(f"telluric {l.wavelength_air_A}: {type(exc).__name__}: {str(exc)[:160]}")
                    per = None
                    break
                s = abs(dadf) * ev["residual_flux"]
                per.append({"wavelength_air_A": float(l.wavelength_air_A), "dA_df": dadf,
                            "sigma_line_dex": s, **ev})
                tot += (s / n) ** 2
            if per is not None:
              comps.append(dict(name="telluric", sigma_dex=math.sqrt(tot), state="MEASURED",
                              source=("per-line measured sky depth -> measured residual of this "
                                      "holding vs an independent correction on the telluric "
                                      "pixels -> dex via the line's own dA/df (+0.1% leg); "
                                      "lines independent, weight 1/n"),
                              evidence={"pool_sha256": digest, "per_line": per}))

    # blends (N I: CN follows A(C))
    if element == "N":
        cm, cp = _leg_lines(unit_dir / "c_minus", lines_stem), _leg_lines(unit_dir / "c_plus", lines_stem)
        if cm is None or cp is None:
            notes.append("blends: A(C) legs missing")
        else:
            pp = _paired(_accp(cp), _accp(cm), n)
            if pp["moved"]:
                notes.append(f"blends: pool moved ({pp['n_paired']} vs {n})")
            else:
                resp = pp["median"] / (2 * C_STEP)          # dA(N)/dA(C)
                comps.append(dict(name="blends", sigma_dex=abs(resp) * SIGMA_C, state="MEASURED",
                                  source=("CN blended into the N I profiles, full molecular synthesis: "
                                          "paired response to A(C) +/- 0.10 on this pool x sigma(C)"),
                                  evidence={"pool_sha256": digest, "dAN_dAC": resp,
                                            "sigma_C": SIGMA_C, "sigma_C_source": SIGMA_C_SOURCE}))
    else:
        comps.append(dict(name="blends", sigma_dex=None, state="N/A",
                          source=("full atomic + molecular synthesis in-window (RYA-1230 turned "
                                  "molecular opacity on in every band); no catalogued-blend deficit"),
                          evidence={"use_molecules": True, "owed_measurement":
                                    "a paired molecule-strength leg on this C/O pool"}))

    # RYA-1232 -- NLTE and 3D sized by ASPLUND+2021's rules (Sect. 2.1 "Uncertainties"),
    # measured on THIS pool instead of RYA-1032's cross-element 0.043:
    #   nlte = 1/2 x the pool's own non-LTE correction, floor 0.03 dex;
    #   3D   = 1/2 x the pool's own 3D effect (3D-NLTE minus 1D-NLTE), added to
    #          `model_atmosphere` in quadrature (their "atmospheric inhomogeneities").
    nd = set(nom.loc[nom.in_aggregate == True, "nlte_delta_dex"].dropna().round(6))  # noqa: E712
    three_d = None
    if nd <= {0.0}:
        comps.append(dict(name="nlte", sigma_dex=None, state="N/A",
                          source="this leg applies no departures (nlte_delta_dex = 0 on every line)",
                          evidence={"nlte_delta_dex": 0.0}))
    else:
        def _deltas(df):
            a = df[df.in_aggregate == True]                                # noqa: E712
            return dict(zip(a.wavelength_air_A.round(3), a.nlte_delta_dex.astype(float)))
        d_here = _deltas(nom)
        d_1d = d_here
        if "ENGINE-A-3DNLTE" in lines_stem:
            sib = nominal_dir / lines_stem.replace("ENGINE-A-3DNLTE", "ENGINE-A")
            d_1d = _deltas(pd.read_csv(sib)) if sib.exists() else None
        if d_1d is None:
            notes.append("nlte: 1D-NLTE sibling missing, cannot split NLTE from 3D")
        else:
            common = sorted(set(d_here) & set(d_1d))
            corr_nlte = float(np.median([d_1d[w] for w in common])) if common else None
            if corr_nlte is None:
                notes.append("nlte: no line common to the 1D-NLTE sibling")
            else:
                comps.append(dict(
                    name="nlte", sigma_dex=max(abs(corr_nlte) / 2.0, 0.03), state="MEASURED",
                    source=("Asplund+2021 Sect. 2.1: half the non-LTE abundance correction, minimum "
                            "0.03 dex -- the correction measured on this pool (median 1D-NLTE "
                            "delta of the accepted lines)"),
                    evidence={"pool_sha256": digest, "median_nlte_correction_dex": corr_nlte,
                              "n_lines": len(common), "floor_dex": 0.03}))
                if "ENGINE-A-3DNLTE" in lines_stem:
                    three_d = float(np.median([d_here[w] - d_1d[w] for w in common]))
    comps.append(dict(name="pseudo_continuum", sigma_dex=None, state="N/A",
                      source=("outside the near-UV the continuum is observed; it is placed per line "
                              "by the standing rule and priced on `continuum`"),
                      evidence={"band": band, "continuum_levels": sorted(
                          round(float(x), 5) for x in acc["continuum_level"].dropna())}))
    hfs = sorted({int(x) for x in g[g.physical_id.isin(ids)]["hfs_n_components"].dropna()})
    if hfs == [1]:
        comps.append(dict(name="hfs_isotopes", sigma_dex=None, state="N/A",
                          source="every pool line has hfs_n_components == 1 (canonical_gf)",
                          evidence={"hfs_n_components": hfs}))
    comps.append(dict(name="molecular_coupling", sigma_dex=None, state="N/A",
                      source="atomic product; the selector is not a MOL- indicator set",
                      evidence={"selector": selector}))
    comps.append(dict(name="holding_instrument", sigma_dex=0.0, state="MEASURED",
                      source="SynthesisHandler harness residual MEASURED against the known optical answer (RYA-869)",
                      evidence={"harness_residual_dex": 0.0}))

    if three_d is not None:
        ma = next((c for c in comps if c["name"] == "model_atmosphere"), None)
        if ma is None:
            notes.append("model_atmosphere: marcs leg missing, 3D term not attached")
        else:
            one_d = float(ma["sigma_dex"])
            ma["sigma_dex"] = float(np.hypot(one_d, three_d / 2.0))
            ma["source"] = (ma["source"] + "; (+) Asplund+2021 Sect. 2.1 atmospheric-inhomogeneity "
                            "term = half the pool's own 3D effect (3D-NLTE minus 1D-NLTE deltas)")
            ma["evidence"] = {**ma["evidence"], "one_d_grid_dex": one_d,
                              "median_3d_effect_dex": three_d}

    scope = product_scope(row, star="solar", indicator_ids=ids)
    numeric = [c for c in comps if c["state"] in {"MEASURED", "DEFINED"}]
    doc = assemble(scope, comps, covariance=np.diag([c["sigma_dex"] ** 2 for c in numeric]).tolist(),
                   covariance_source=("components independent of one another; the shared-source "
                                      "gf correlation lives inside transition_data"),
                   assumptions="solar log g and [Fe/H] pinned; Teff allowance 1 K")
    verdict = None
    try:
        validate(doc, scope=scope)
    except UncertaintyError as exc:
        verdict = str(exc)
    if restat:
        notes.append("leg-unstable excluded: " + ", ".join(restat["leg_unstable"]))
    return {"row": row, "budget": doc, "ids": ids, "notes": notes, "verdict": verdict,
            "restat": restat,
            "prod_stem": prod_stem, "nominal_dir": str(nominal_dir)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--legs", type=Path, required=True)
    ap.add_argument("--nominal", type=Path, nargs="*", default=[])
    ap.add_argument("--stage", type=Path, required=True)
    ap.add_argument("--report", type=Path, required=True)
    a = ap.parse_args()
    gf = pd.read_csv(GF, low_memory=False)
    units = [shlex.split(l) for l in (a.legs / "units.txt").read_text().splitlines() if l.strip()]
    a.stage.mkdir(parents=True, exist_ok=True)
    out = []
    for i, args in enumerate(units):
        udir = a.legs / f"unit{i:02d}"
        stems = sorted(Path(f).name for f in glob.glob(str(udir / "nominal" / "*_lines.csv")))
        for stem in stems:
            nd = (udir / "nominal") if (udir / "nominal" / stem).exists() else \
                next((d for d in a.nominal if (d / stem).exists()), None)
            if nd is None:
                out.append({"unit": i, "stem": stem, "skip": "no nominal artifact"})
                continue
            r = build(args, stem, nd, udir, gf)
            rec = {"unit": i, "stem": stem, "A": r["row"]["A"], "key_treatment": r["row"]["treatment"],
                   "holding": r["row"]["holding"], "selector": r["row"]["selector"],
                   "element": r["row"]["element"], "band": r["row"]["band"]}
            if "skip" in r:
                rec["skip"] = r["skip"]
                out.append(rec)
                continue
            b = r["budget"]
            rec.update({"holds": b["holds"], "verdict": r["verdict"], "notes": r["notes"],
                        "prod_stem": r["prod_stem"],
                        "sigma_reported": b["sigma_reported"],
                        "components": {c["name"]: (c["state"], c["sigma_dex"]) for c in b["components"]},
                        "telluric_per_line": [
                            {k: v for k, v in pl.items() if k in ("wavelength_air_A", "dA_df", "sigma_line_dex",
                                                                  "n_telluric_px", "max_depth", "residual_flux",
                                                                  "residual_basis")}
                            for c in b["components"] if c["name"] == "telluric"
                            for pl in c["evidence"].get("per_line", [])]})
            out.append(rec)
            if r["verdict"] is None:
                df = pd.read_csv(Path(r["nominal_dir"]) / r["prod_stem"])
                if r.get("restat"):
                    rs = r["restat"]
                    df["n_excluded"] = df["n_excluded"] + (df["n_lines"] - rs["n_lines"])
                    df["A"], df["stat_dex"], df["n_lines"] = rs["A"], rs["stat_dex"], rs["n_lines"]
                    df["leg_unstable_excluded"] = json.dumps(rs["leg_unstable"])
                df["uncertainty"] = json.dumps(b)
                df["uncertainty_indicator_ids"] = json.dumps(r["ids"])
                df["sigma_reported"] = b["sigma_reported"]
                df.to_csv(a.stage / r["prod_stem"], index=False)
    a.report.parent.mkdir(parents=True, exist_ok=True)
    a.report.write_text(json.dumps({"ticket": "RYA-1230", "products": out}, indent=1, default=str) + "\n")
    ok = [o for o in out if o.get("verdict") is None and "skip" not in o]
    print(f"{len(out)} products: {len(ok)} complete, "
          f"{sum('skip' in o for o in out)} skipped, "
          f"{sum(bool(o.get('verdict')) for o in out)} held")
    for o in out:
        tag = o.get("skip") or o.get("verdict") or f"sigma_reported={o['sigma_reported']:.4f}"
        print(f"  u{o['unit']:02d} {o.get('element','?')} {o.get('band','?'):11} "
              f"{str(o.get('holding'))[:24]:24} {str(o.get('selector')):11} "
              f"{str(o.get('key_treatment')):16} A={o.get('A')}  {str(tag)[:90]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
