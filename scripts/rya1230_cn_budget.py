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
}
DELTA_XI, XI_STEP, CO_STEP, CSCALE = 0.2912, 0.10, 0.10, 0.001
SIGMA_C, SIGMA_O = 0.068, 0.076
CO_SOURCE = ("RYA-1230 published solar C (HARPS VIS AGSS21 3D-NLTE 8.460 +/- 0.068) and O "
             "(IAG red-optical AGSS21 3D-NLTE 8.643 +/- 0.076), full RYA-587 budgets")
TD_SIGMA = math.log10(11.08 / 8.50)
TD_SOURCE = ("Brooke et al. 2014, ApJS 210, 23, Table 5: computed tau(A2Pi1/2, v=0) = 11.08 us vs "
             "measured 8.50 +/- 0.05 us (Taherian & Slanger 1984); f ~ 1/tau -> log10(11.08/8.50)")
MODEL_FORM_BOUND = 0.061


def leg(out: Path, region: str, name: str) -> dict | None:
    d = out / region / name
    f = d / f"solar_{region}_cno_per_band.csv"
    if not (d / "DONE").exists() or not f.exists():
        return None
    r = pd.read_csv(f)
    r = r[r.element == "N"].iloc[0]
    return {"A": float(r.A_X), "sigma_fit": float(r.sigma_fit), "n_pix": int(r.n_pix),
            "red_chi2": float(r.red_chi2), "constrained": bool(r.constrained), "dir": d}


def windows_of(region: str):
    from pipeline.cno_synthesis import REGION_DIAGNOSTICS
    return [w for dg in REGION_DIAGNOSTICS[region] if dg.element == "N" for w in dg.windows_A]


def n_cn_lines(region: str) -> int:
    from rya1230_n_product_hygiene import _bsyn_species_in
    return int(_bsyn_species_in(windows_of(region)).get("12C14N", 0))


def build(region: str, legs: Path, others: dict) -> dict:
    import rya1230_cno_budget as B
    instrument, holding = REGIONS[region]
    L = {k: leg(legs, region, k) for k in ("nominal", "xi_minus", "xi_plus", "q80", "q97", "contref", "marcs",
                                           "cscale", "win60", "c_minus", "c_plus", "o_minus", "o_plus")}
    nom = L["nominal"]
    win = windows_of(region)
    ind = ("CN_AX_IR:12C14N:A-X(0-0):" + region + ":windows_sha256:"
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
    comps.append(uc.transition_data(ids, [TD_SIGMA], [1.0], covariance=[[TD_SIGMA ** 2]],
                                    sources=[TD_SOURCE], covariance_source=TD_SOURCE))
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
            ev.update(bound=True, measured_leg_dex=v, bound_dex=bound,
                      bound_source="published 1D->3D shift for 12C14N A-X (0-0), Amarsi et al. 2021 via RYA-1220")
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
          bound=MODEL_FORM_BOUND)
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
        if lost or (uncovered and ratio_den <= 0):
            notes.append(f"telluric: {len(lost)} window(s) unreadable, {len(uncovered)} without an "
                         f"independent reference and no covered window to scale from: {lost[:2]}")
        else:
            ratio = ratio_num / ratio_den if ratio_den > 0 else 0.0
            for r in uncovered:
                r["residual_flux"] = ratio * r["sky_absorption"]
                r["residual_basis"] += (f"; BOUNDED by this holding's measured residual/absorption "
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
    if all(L[k] for k in ("c_minus", "c_plus", "o_minus", "o_plus")):
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
    if len(vals) >= 2:
        comps.append(dict(name="holding_instrument", sigma_dex=float(np.std(vals, ddof=1)), state="MEASURED",
                          source="std of CN A-X (0-0) across holdings, same star, same diagnostic family",
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
                      source=("12C14N windows (AGSS21 CN 0-0 positions); 13C14N and 12C15N enter the "
                              "synthesis at the solar isotope ratios and are not the fitted lines"),
                      evidence={"owed_measurement": "a paired 12C/13C leg on this band"}))
    row = dict(element="N", ion="I", band="NIR", instrument=instrument, treatment="1D-LTE",
               handler="CNOSynthesis", A=nom["A"], n_lines=n_cn_lines(region), n_excluded=0,
               stat_dex=nom["sigma_fit"], syst_dex=np.nan,
               stat_basis=("measured -- 1 sigma from the chi2 curvature of THIS band's fit "
                           "(pipeline.fit_constraint), rescaled to red_chi2 = 1"),
               dominant="", route="synth", scale="1D-LTE", model="none", atmos="atlas9",
               gf="brooke2014", route_basis="handler", deck="none")
    from publish_product import normalise
    pub = normalise(pd.DataFrame([row]), holding=holding, tier="ALL", route="SYNTH",
                    selector="MOL-CN_AX_IR")[0]
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
    a.stage.mkdir(parents=True, exist_ok=True)
    out = []
    for region in REGIONS:
        if nominal[region] is None:
            out.append({"region": region, "skip": "no nominal leg"}); continue
        r = build(region, a.legs, nominal)
        b = r["budget"]
        rec = {"region": region, "holding": r["holding"], "A": r["row"]["A"], "verdict": r["verdict"],
               "notes": r["notes"], "sigma_reported": b["sigma_reported"], "holds": b["holds"],
               "components": {c["name"]: (c["state"], c["sigma_dex"]) for c in b["components"]}}
        if r["verdict"] is None:
            stem = f"NI_MOL_CN_AX_IR_{r['row']['instrument']}_{r['holding']}_SYNTH_MOL-CN_AX_IR_products.csv"
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
