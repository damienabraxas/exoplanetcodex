#!/usr/bin/env python3
"""RYA-1230 -- RYA-587 budgets for the molecular CN A-X (0-0) products.

    python3 scripts/rya1230_cn_budget.py --legs DIR --stage DIR --report R.json

One indicator per product (the band fit). Components, each from the region's OWN legs
(`scripts/rya1230_cn_budget_legs.py`) or a cited measurement:

  measurement         fit curvature sigma (pipeline.fit_constraint), profile likelihood
  transition_data     CN A-X band-strength scale: Brooke et al. 2014 (ApJS 210, 23) Table 5
                      computes tau(A2Pi, v=0) = 11.08 us against the measured 8.50 +/- 0.05 us
                      (Taherian & Slanger 1984); f ~ 1/tau, so log10(11.08/8.50) = 0.115 dex,
                      fully correlated across the band
  stellar.xi          paired xi 0.90/1.10, scaled to RYA-1089 delta_xi = 0.2912
  stellar.teff/logg/  DEFINED zero (solar convention, RYA-1226)
  metallicity
  continuum           |A(q97) - A(q80)| / 2, model-guided continuum
  profile_ew          |A(windows x0.6) - A(nominal)|
  model_atmosphere    max(|A(MARCS) - A(ATLAS9)|, 0.061 published 1D->3D CN A-X shift,
                      Amarsi 2021 via RYA-1220) -- the larger, flagged when it is the bound
  telluric            per fit window: measured residual vs an independent correction ->
                      dex via the band's dA/df (+0.1% leg); windows independent
  molecular_coupling  paired A(C), A(O) +/- 0.10 Jacobians x our published C/O budgets
                      (C 8.460 +/- 0.068, O 8.643 +/- 0.076), C-O correlation unmeasured ->
                      evaluated at the WORST rho (RYA-1220 policy)
  holding_instrument  std of the same diagnostic across the three holdings (IAG, KP, CRIRES+ J)
  nlte                N/A: molecular band synthesised in LTE (cno_synthesis lte_molecular_band)
  blends / pseudo_continuum / hfs_isotopes   N/A with evidence
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from pipeline import uncertainty_contract as uc  # noqa: E402
from pipeline import nitrogen_uncertainty as nu  # noqa: E402

REGIONS = {
    "nir_cn_iag": ("iag_fts_solar_atlas", "solar_iag"),
    "nir_cn_kp": ("kpno_solar_atlas", "solar_kpno_molecfit_corrected"),
    "j_cn_crires": ("crires_plus", "solar_crires_plus_j_rya1219"),
    #: RYA-1232: CRIRES+ K CO first overtone -- the only CRIRES+ K carbon indicator
    "k_co_crires": ("crires_plus", "solar_crires_plus_k_rya1219"),
}
DELTA_XI, XI_STEP, CO_STEP, CSCALE = 0.2912, 0.10, 0.10, 0.001
SIGMA_C, SIGMA_O = 0.068, 0.076
CO_SOURCE = ("RYA-1230 published solar C (HARPS VIS AGSS21 3D-NLTE 8.460 +/- 0.068) and O "
             "(IAG red-optical AGSS21 3D-NLTE 8.643 +/- 0.076), full RYA-587 budgets")
TD_SIGMA = math.log10(11.08 / 8.50)
TD_SOURCE = ("Brooke et al. 2014, ApJS 210, 23, Table 5: computed tau(A2Pi1/2, v=0) = 11.08 us vs "
             "measured 8.50 +/- 0.05 us (Taherian & Slanger 1984); f ~ 1/tau -> log10(11.08/8.50)")
MODEL_FORM_BOUND = 0.061
#: RYA-1232 -- per-region indicator spec. CN keeps exactly the values above; CO (measured
#: element C, the less abundant one, as Amarsi+2021 fit it) carries its own sources.
CO_TD_SIGMA = math.log10(1.005)
CO_TD_SOURCE = ("Li et al. 2015, ApJS 216, 15, Sect. 3: the 2-0 band intensities reproduce Malathy "
                "Devi et al. 2012b 'within approximately 0.1%' P(15)-R(25), up to 0.5% at J=30; the "
                "worst stated residual (0.5%) is adopted -> log10(1.005)")
_CN = dict(element="N", ion="I", diag="CN_AX_IR", molecule="12C14N", selector="MOL-CN_AX_IR",
           gf="brooke2014", td=(TD_SIGMA, TD_SOURCE), stem="NI_MOL_CN_AX_IR",
           bound=MODEL_FORM_BOUND,
           bound_source="published 1D->3D shift for 12C14N A-X (0-0), Amarsi et al. 2021 via RYA-1220",
           partners=("C", "O"), family="CN A-X (0-0)")
RATIO_FROM = {"k_co_crires": "j_cn_crires"}
MEASURED_RATIO: dict = {}
SPEC = {"nir_cn_iag": _CN, "nir_cn_kp": _CN, "j_cn_crires": _CN,
        "k_co_crires": dict(element="C", ion="I", diag="CO_K", molecule="12C16O", selector="MOL-CO_K",
                            gf="li2015", td=(CO_TD_SIGMA, CO_TD_SOURCE), stem="CI_MOL_CO_K",
                            bsyn_species="16O12C",
                            bound=0.140,
                            bound_source=("published 1D->3D shift for 12C16O X-X, Amarsi et al. 2021 "
                                          "Table 2 via data/reference/molecular_cno_literature_rya1220 "
                                          "(3D 8.47 vs MARCS 8.61)"),
                            partners=("O",), family="12C16O X-X dv=2")}


def leg(out: Path, region: str, name: str) -> dict | None:
    d = out / region / name
    f = d / f"solar_{region}_cno_per_band.csv"
    if not (d / "DONE").exists() or not f.exists():
        return None
    r = pd.read_csv(f)
    r = r[r.element == SPEC[region]["element"]].iloc[0]
    return {"A": float(r.A_X), "sigma_fit": float(r.sigma_fit), "n_pix": int(r.n_pix),
            "red_chi2": float(r.red_chi2), "constrained": bool(r.constrained), "dir": d}


def windows_of(region: str):
    from pipeline.cno_synthesis import REGION_DIAGNOSTICS
    return [w for dg in REGION_DIAGNOSTICS[region] if dg.element == SPEC[region]["element"]
            for w in dg.windows_A]


def n_cn_lines(region: str) -> int:
    from rya1230_n_product_hygiene import _bsyn_species_in
    #: the .bsyn lists name species their own way (CO is '16O12C'); count by that key
    sp = SPEC[region]
    return int(_bsyn_species_in(windows_of(region)).get(sp.get("bsyn_species", sp["molecule"]), 0))


def build(region: str, legs: Path, others: dict) -> dict:
    import rya1230_cno_budget as B
    instrument, holding = REGIONS[region]
    sp = SPEC[region]
    L = {k: leg(legs, region, k) for k in ("nominal", "xi_minus", "xi_plus", "q80", "q97", "contref", "marcs",
                                           "cscale", "win60", "c_minus", "c_plus", "o_minus", "o_plus")}
    nom = L["nominal"]
    win = windows_of(region)
    ind = (f"{sp['diag']}:{sp['molecule']}:{sp['family']}:" + region + ":windows_sha256:"
           + hashlib.sha256(json.dumps(win).encode()).hexdigest()[:24])
    ids = [ind]
    digest = uc.pool_digest(ids)
    notes, comps = [], []
    comps.append(uc.fit_measurement(
        nom["sigma_fit"], source=f"{nom['dir'].name} fit, pipeline.fit_constraint curvature",
        likelihood="chi2 curvature of the single-abundance band fit, rescaled to red_chi2 = 1",
        correlation_treatment=("pixels treated as independent after the red_chi2 rescale; "
                               "residual pixel covariance not modelled"),
        n_pixels=int(nom["n_pix"])))
    _tds, _tdsrc = sp["td"]
    comps.append(uc.transition_data(ids, [_tds], [1.0], covariance=[[_tds ** 2]],
                                    sources=[_tdsrc], covariance_source=_tdsrc))
    for name, why in (("stellar.logg", "solar log g is pinned (RYA-1089)"),
                      ("stellar.metallicity", "solar [Fe/H] is pinned (RYA-1089)")):
        comps.append(dict(name=name, sigma_dex=0.0, state="DEFINED", source=why,
                          evidence={"delta": 0.0, "parameter_exists": True}))
    comps.append(dict(name="stellar.teff", sigma_dex=0.0, state="DEFINED",
                      source="sourced DEFINED zero (RYA-1226 convention): solar Teff allowance 1 K",
                      evidence={"parameter_exists": True, "delta_K": 1.0}))
    if L["xi_minus"] and L["xi_plus"]:
        s = DELTA_XI / XI_STEP
        dp, dm = (L["xi_plus"]["A"] - nom["A"]) * s, (L["xi_minus"]["A"] - nom["A"]) * s
        signed = (dp - dm) / 2
        comps.append(dict(name="stellar.xi", sigma_dex=abs(signed), state="MEASURED",
                          source="paired xi 0.90/1.10 refits of this band, stamped legs",
                          evidence={"pool_sha256": digest, "parameter_source": "RYA-1089 delta_xi 0.2912 km/s",
                                    "response_assessment": f"A(0.9)={L['xi_minus']['A']}, A(1.0)={nom['A']}, A(1.1)={L['xi_plus']['A']}",
                                    "signed_response_dex": signed, "delta_plus_dex": dp,
                                    "delta_minus_dex": dm, "delta_parameter": DELTA_XI,
                                    "parameter_exists": True}))
    else:
        notes.append("xi legs missing")

    def lever(name, key, src, bound=None):
        if not L[key]:
            notes.append(f"{name}: leg {key} missing"); return
        v = abs(L[key]["A"] - nom["A"])
        ev = {"pool_sha256": digest, "nominal_A": nom["A"], "leg_A": L[key]["A"]}
        if bound is not None and bound > v:
            ev.update(bound=True, measured_leg_dex=v, bound_dex=bound, bound_source=sp["bound_source"])
            v = bound
        comps.append(dict(name=name, sigma_dex=v, state="DEFINED" if ev.get("bound") else "MEASURED",
                          source=src, evidence=ev))

    if L["q80"] and L["q97"] and L["contref"]:
        _pl = (L["q97"]["A"] - L["q80"]["A"]) / 2
        _rf = L["contref"]["A"] - L["nominal"]["A"]
        comps.append(dict(name="continuum", sigma_dex=float((_pl ** 2 + _rf ** 2) ** 0.5),
                          state="MEASURED",
                          source=("continuum (synthesis reference, 3-MAD clipped): placement q97 vs q80 "
                                  "half-difference + REFERENCE (IAG-atlas leg minus nominal), in quadrature"),
                          evidence={"pool_sha256": digest, "A_q80": L["q80"]["A"], "A_q97": L["q97"]["A"],
                                    "A_contref_iag": L["contref"]["A"]}))
    else:
        notes.append("continuum legs missing")
    lever("profile_ew", "win60", "fit sub-windows scaled x0.6 about their centres (continuum on unscaled windows)")
    lever("model_atmosphere", "marcs", "MARCS.GES vs ATLAS9.Castelli; bounded by the published 1D->3D shift",
          bound=sp["bound"])
    if L["cscale"]:
        dadf = (L["cscale"]["A"] - nom["A"]) / (-CSCALE)
        res, lost = [], []
        for lo, hi in win:
            try:
                ev = B.telluric_line(instrument, holding, 0.5 * (lo + hi), 0.5 * (hi - lo),
                                     min(w[0] for w in win), max(w[1] for w in win))
                res.append({"window_A": [lo, hi], **{k: ev[k] for k in ev if k in (
                    "n_telluric_px", "max_depth", "residual_flux", "residual_basis",
                    "sky_pair", "sky_absorption")}})
            except Exception as exc:                                        # noqa: BLE001
                lost.append({"window_A": [lo, hi], "why": f"{type(exc).__name__}: {str(exc)[:140]}"})
        covered = [r for r in res if r["residual_flux"] is not None]
        uncovered = [r for r in res if r["residual_flux"] is None]
        ratio_num = sum(r["residual_flux"] for r in covered if r["sky_absorption"] > 0)
        ratio_den = sum(r["sky_absorption"] for r in covered if r["sky_absorption"] > 0)
        transferred = None
        if uncovered and ratio_den <= 0 and region in RATIO_FROM and RATIO_FROM[region] in MEASURED_RATIO:
            #: RYA-1232: no K-band correction exists to difference against. Transfer the
            #: residual/absorption ratio MEASURED on the J arm -- same instrument, same night,
            #: same RYA-1219 molecfit pipeline -- and bound K by its own molecfit sky depth.
            transferred = RATIO_FROM[region]
            ratio_num, ratio_den = MEASURED_RATIO[transferred], 1.0
        if lost or (uncovered and ratio_den <= 0):
            notes.append(f"telluric: {len(lost)} window(s) unreadable, {len(uncovered)} without an "
                         f"independent reference and no covered window to scale from: {lost[:2]}")
        else:
            ratio = ratio_num / ratio_den if ratio_den > 0 else 0.0
            if covered:
                MEASURED_RATIO[region] = ratio
            for r in uncovered:
                r["residual_flux"] = ratio * r["sky_absorption"]
                r["residual_basis"] += (
                    f"; BOUNDED by {transferred}'s measured residual/absorption ratio {ratio:.4f} "
                    "(same instrument, night and molecfit pipeline)" if transferred else
                    f"; BOUNDED by this holding's measured residual/absorption "
                    f"ratio {ratio:.4f} on its {len(covered)} covered window(s)")
            n = len(res)
            sig = math.sqrt(sum((abs(dadf) * r["residual_flux"] / n) ** 2 for r in res))
            comps.append(dict(name="telluric", sigma_dex=sig,
                              state="DEFINED" if uncovered else "MEASURED",
                              source=("per window: measured sky depth -> residual vs an independent "
                                      "correction -> dex via the band's dA/df (+0.1% leg)"
                                      + ("; windows with no independent reference bounded by the "
                                         "holding's own measured residual/absorption ratio" if uncovered else "")),
                              evidence={"pool_sha256": digest, "dA_df": dadf, "per_window": res,
                                        "bound": bool(uncovered),
                                        "residual_to_absorption_ratio": ratio,
                                        "n_windows_bounded": len(uncovered)}))
    else:
        notes.append("cscale leg missing")
    if sp["partners"] == ("O",) and L["o_minus"] and L["o_plus"]:
        #: CO: C is the fitted element; its only partner is O (pinned)
        jo = (L["o_plus"]["A"] - L["o_minus"]["A"]) / (2 * CO_STEP)
        comps.append(uc.component(
            "molecular_coupling", abs(jo) * SIGMA_O, state="MEASURED",
            source="paired A(O) +/- 0.10 refits of this band (pin) x published sigma(O)",
            evidence={"responses": {"O": uc.paired_response(
                {ind: nom["A"]}, {ind: L["o_minus"]["A"]}, {ind: L["o_plus"]["A"]}, delta=CO_STEP,
                source="paired A(O) pins", parameter_source="Numerical A(O) probe only; not adopted sigma")},
                "abundance_order": ["O"], "abundance_covariance": [[SIGMA_O ** 2]],
                "covariance_source": CO_SOURCE,
                "response_assessment": f"dA(C)/dA(O) = {jo:+.3f} per dex",
                "pool_sha256": digest}))
    elif sp["partners"] == ("C", "O") and all(L[k] for k in ("c_minus", "c_plus", "o_minus", "o_plus")):
        jc = (L["c_plus"]["A"] - L["c_minus"]["A"]) / (2 * CO_STEP)
        jo = (L["o_plus"]["A"] - L["o_minus"]["A"]) / (2 * CO_STEP)
        rho = 1.0 if jc * jo >= 0 else -1.0                     # the worst-case correlation
        cov = [[SIGMA_C ** 2, rho * SIGMA_C * SIGMA_O], [rho * SIGMA_C * SIGMA_O, SIGMA_O ** 2]]
        comps.append(nu.molecular_coupling(
            {ind: nom["A"]}, {ind: L["c_minus"]["A"]}, {ind: L["c_plus"]["A"]},
            {ind: L["o_minus"]["A"]}, {ind: L["o_plus"]["A"]},
            carbon_step=CO_STEP, oxygen_step=CO_STEP,
            source="paired A(C), A(O) +/- 0.10 refits of this band (pins)",
            abundance_covariance=cov,
            covariance_source=(CO_SOURCE + f"; C-O correlation unmeasured, evaluated at the WORST "
                               f"rho = {rho:+.0f} for these Jacobians (RYA-1220 policy)"),
            response_assessment=(f"dA(N)/dA(C) = {jc:+.3f}, dA(N)/dA(O) = {jo:+.3f} per dex; "
                                 f"both sides retained in the responses")))
    else:
        notes.append("C/O legs missing")
    vals = [v for v in others.values() if v is not None]
    if len(vals) < 2 and region == "k_co_crires":
        #: one holding observes 12C16O dv=2 (KP and IAG stop short of K): the band route's
        #: definition of the term applies (RYA-869 harness residual)
        comps.append(dict(name="holding_instrument", sigma_dex=0.0, state="MEASURED",
                          source=("SynthesisHandler harness residual MEASURED against the known optical "
                                  "answer (RYA-869); no second holding covers 12C16O dv=2"),
                          evidence={"harness_residual_dex": 0.0, "per_holding_A": others}))
    elif len(vals) >= 2:
        comps.append(dict(name="holding_instrument", sigma_dex=float(np.std(vals, ddof=1)), state="MEASURED",
                          source=f"std of {sp['family']} across holdings, same star, same diagnostic family",
                          evidence={"per_holding_A": others}))
    comps.append(dict(name="nlte", sigma_dex=None, state="N/A",
                      source=("molecular band synthesised in LTE (cno_synthesis: lte_molecular_band, no NLTE grid); "
                              "Amarsi+2021 carry no non-LTE term for molecules either"),
                      evidence={"nlte_flag": "lte_molecular_band"}))
    comps.append(dict(name="blends", sigma_dex=None, state="N/A",
                      source="full atomic + molecular synthesis in every fit window",
                      evidence={"use_molecules": True}))
    comps.append(dict(name="pseudo_continuum", sigma_dex=None, state="N/A",
                      source="continuum placed per window by the model-guided rule; priced on `continuum`",
                      evidence={}))
    comps.append(dict(name="hfs_isotopes", sigma_dex=None, state="N/A",
                      source=(f"{sp['molecule']} windows; minor isotopologues enter the synthesis at "
                              "the solar isotope ratios and are not the fitted lines"),
                      evidence={"owed_measurement": "a paired 12C/13C leg on this band"}))
    from pipeline.band_policy import resolve as _band_of
    _w = windows_of(region)
    _band = _band_of(0.5 * (min(w[0] for w in _w) + max(w[1] for w in _w))).name \
        if region == "k_co_crires" else "NIR"
    row = dict(element=sp["element"], ion=sp["ion"], band=_band, instrument=instrument, treatment="1D-LTE",
               handler="CNOSynthesis", A=nom["A"], n_lines=n_cn_lines(region), n_excluded=0,
               stat_dex=nom["sigma_fit"], syst_dex=np.nan,
               stat_basis=("measured -- 1 sigma from the chi2 curvature of THIS band's fit "
                           "(pipeline.fit_constraint), rescaled to red_chi2 = 1"),
               dominant="", route="synth", scale="1D-LTE", model="none", atmos="atlas9",
               gf=sp["gf"], route_basis="handler", deck="none")
    from publish_product import normalise
    pub = normalise(pd.DataFrame([row]), holding=holding, tier="ALL", route="SYNTH",
                    selector=sp["selector"])[0]
    pub["star"] = "solar"
    scope = uc.product_scope(pub, star="solar", indicator_ids=ids)
    numeric = [c for c in comps if c["state"] in {"MEASURED", "DEFINED"}]
    doc = uc.assemble(scope, comps, covariance=np.diag([c["sigma_dex"] ** 2 for c in numeric]).tolist(),
                      covariance_source="components independent of one another",
                      assumptions="solar log g and [Fe/H] pinned; Teff allowance 1 K")
    verdict = None
    try:
        uc.validate(doc, scope=scope)
    except uc.UncertaintyError as exc:
        verdict = str(exc)
    return {"region": region, "holding": holding, "row": row, "budget": doc, "ids": ids,
            "verdict": verdict, "notes": notes}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--legs", type=Path, required=True)
    ap.add_argument("--stage", type=Path, required=True)
    ap.add_argument("--report", type=Path, required=True)
    a = ap.parse_args()
    nominal = {r: (leg(a.legs, r, "nominal") or {}).get("A") for r in REGIONS}
    #: holding_instrument compares a diagnostic FAMILY across holdings: CN with CN only
    fam = {r: SPEC[r]["family"] for r in REGIONS}
    a.stage.mkdir(parents=True, exist_ok=True)
    out = []
    for region in REGIONS:
        if nominal[region] is None:
            out.append({"region": region, "skip": "no nominal leg"}); continue
        r = build(region, a.legs, {k: v for k, v in nominal.items() if fam[k] == fam[region]})
        b = r["budget"]
        rec = {"region": region, "holding": r["holding"], "A": r["row"]["A"], "verdict": r["verdict"],
               "element": SPEC[region]["element"], "selector": SPEC[region]["selector"],
               "key_treatment": r["row"]["treatment"],
               "notes": r["notes"], "sigma_reported": b["sigma_reported"], "holds": b["holds"],
               "components": {c["name"]: (c["state"], c["sigma_dex"]) for c in b["components"]}}
        if r["verdict"] is None:
            _sp = SPEC[region]
            stem = f"{_sp['stem']}_{r['row']['instrument']}_{r['holding']}_SYNTH_{_sp['selector']}_products.csv"
            row = dict(r["row"])
            #: the published sigma_syst is the systematic part of the canonical RYA-587 total
            row["syst_dex"] = round(math.sqrt(max(b["sigma_reported"] ** 2 - row["stat_dex"] ** 2, 0.0)), 4)
            df = pd.DataFrame([row])
            df["uncertainty"] = json.dumps(b)
            df["uncertainty_indicator_ids"] = json.dumps(r["ids"])
            df["sigma_reported"] = b["sigma_reported"]
            df.to_csv(a.stage / stem, index=False)
            rec["prod_stem"] = stem
        out.append(rec)
        print(f"{region:12} A={rec['A']}  " + (f"sigma_reported={rec['sigma_reported']:.4f}"
              if rec["verdict"] is None else f"HELD: {rec['verdict']} {rec['notes']}"))
        for k, (st, v) in rec["components"].items():
            print(f"      {k:20} {st:9} {v}")
    a.report.parent.mkdir(parents=True, exist_ok=True)
    a.report.write_text(json.dumps({"ticket": "RYA-1230", "products": out}, indent=1, default=str) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
