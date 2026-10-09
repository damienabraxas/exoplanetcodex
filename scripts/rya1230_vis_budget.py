#!/usr/bin/env python3
"""RYA-1230 -- RYA-587 budgets for the HARPS `vis` region diagnostics (cno_synthesis joint
C/N/O fit): CH G-band, C2 Swan, C I 5052, C I 5380, CN red, [O I] 6300.

    python3 scripts/rya1230_vis_budget.py --legs DIR --stage DIR --report R.json

Each diagnostic is one indicator. Components from the `vis` region's OWN legs
(`scripts/rya1230_cn_budget_legs.py --regions vis`) or a cited/measured source:

  measurement         fit curvature sigma, profile likelihood
  transition_data     C I / [O I]: canonical_gf per-line sigma (NIST class, worse class on
                      conflict). CN red: Brooke 2014 A-state lifetime scale (0.115 dex).
                      CH, C2: CONSISTENCY BOUND |A(mol) - A(C I)| (+) sigma(C I) against our
                      HARPS AGSS21 C I 1D-LTE product on the SAME spectrum and model class;
                      the laboratory lifetime papers (Luque & Crosley 1996; Brooke 2013) are
                      not obtainable, so the band-strength scale is bounded by the Sun itself
  stellar.xi / continuum / profile_ew / model_atmosphere / telluric   paired legs, as CN
                      (molecules: model_atmosphere bounded by Amarsi 2021's 3D-MARCS shift)
  molecular_coupling  MOL- selectors: paired pins of the OTHER elements (CH/C2 -> N, O;
                      CN red -> C, O) x our published budgets, worst rho
  blends              [O I] 6300: paired Ni pins x AGSS21 sigma(Ni); others N/A (full synthesis)
  holding_instrument  std of the diagnostic across vis / vis_kp_k05 / vis_kp_mf / vis_iag
"""
from __future__ import annotations

import argparse
import csv
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

KEYS = {"CH_Gband": ("C", "MOL-CH_Gband", "CH"), "C2_Swan": ("C", "MOL-C2_Swan", "C2"),
        "CI_5052": ("C", "ATOM-CI_5052", None), "CI_5380": ("C", "ATOM-CI_5380", None),
        "CN_red": ("N", "MOL-CN_red", "CN"), "OI_6300": ("O", "FORB-OI_6300", None)}
TWINS = ("vis", "vis_kp_k05", "vis_kp_mf", "vis_iag")
INSTRUMENT, HOLDING = "harps", "solar_harps_molecfit_corrected"
DELTA_XI, XI_STEP, STEP, CSCALE = 0.2912, 0.10, 0.10, 0.001
#: our published budgets for the coupling partners (RYA-1230 headlines)
SIGMA = {"C": 0.068, "N": 0.081, "O": 0.076}
SIGMA_NI = 0.04
SIGMA_NI_SOURCE = "Asplund, Amarsi & Grevesse 2021 (A&A 653, A141) Table 2: A(Ni) = 6.20 +/- 0.04"
CN_TD = math.log10(11.08 / 8.50)
CN_TD_SOURCE = ("Brooke et al. 2014, ApJS 210, 23, Table 5: tau(A2Pi, v=0) computed 11.08 us vs "
                "measured 8.50 +/- 0.05 us (Taherian & Slanger 1984); f ~ 1/tau")
MODEL_FORM = ROOT / "data/reference/molecular_cno_literature_rya1220/amarsi2021_molecular_reference_by_species.csv"


def per_band(legs: Path, region: str, leg: str) -> pd.DataFrame | None:
    d = legs / region / leg
    f = d / f"solar_{region}_cno_per_band.csv"
    return pd.read_csv(f) if (d / "DONE").exists() and f.exists() else None


def val(df, key):
    if df is None:
        return None
    r = df[df.key == key]
    return None if r.empty or not np.isfinite(r.A_X.iloc[0]) else float(r.A_X.iloc[0])


#: the band each vis molecular diagnostic fits -- Amarsi 2021's 3D-MARCS shift applies only
#: where it was measured on the SAME system (C2 Swan). Its CH entry is X-X (IR
#: vibration-rotation), not the A-X G-band, and its CN entry is A-X (0-0) in the IR, not the
#: red-system bands here: borrowing those is the wrong population (RYA-1220's CN-red ruling).
VIS_SYSTEM = {"C2": "Swan"}


def model_form(mol: str) -> float | None:
    system = VIS_SYSTEM.get(mol)
    if system is None:
        return None
    with MODEL_FORM.open() as fh:
        for r in csv.DictReader(fh):
            if r["molecule"] == mol and r["system"] == system:
                return abs(float(r["model_form_3D_minus_MARCS"]))
    return None


def gf_sigma(species: str, lo: float, hi: float):
    from pipeline.gf_grades import nist_sigma_dex
    g = pd.read_csv(ROOT / "data/linelists/canonical_gf.csv", low_memory=False)
    g = g[(g.species == species) & g.wavelength_air_A.between(lo, hi)]
    if g.empty:
        return None, None
    c = 0.5 * (lo + hi)
    r = g.iloc[(g.wavelength_air_A - c).abs().argsort()].iloc[0]
    s, ref = float(r.gf_sigma_dex), str(r.loggf_reference)
    if math.isnan(s):
        classes = {str(r.nist_grade).strip()} | ({str(r.gf_tier).replace("NIST-", "")}
                                                 if str(r.gf_tier).startswith("NIST-") else set())
        classes = {k for k in classes if not math.isnan(nist_sigma_dex(k))}
        if not classes:
            return None, None
        worst = max(classes, key=nist_sigma_dex)
        s, ref = nist_sigma_dex(worst), f"{ref}; NIST class {worst} via gf_grades.nist_sigma_dex"
    return s, f"{species} {float(r.wavelength_air_A):.3f} A ({r.physical_id}): {ref}"


def ci_reference():
    """Our HARPS AGSS21 C I 1D-LTE product: same spectrum, same 1D-LTE model class."""
    d = json.loads((ROOT / "data/products/solar/C.json").read_text())
    for p in d["products"]:
        if (p["band"] == "VIS" and p["holding"] == HOLDING and p.get("selector") == "SET-AGSS21"
                and p["treatment"] == "1D-LTE" and "uncertainty" in p):
            return float(p["A"]), float(p["sigma_reported"])
    return None, None


def build(key, legs, twins_A):
    import rya1230_cno_budget as B
    from pipeline.cno_synthesis import REGION_DIAGNOSTICS
    element, selector, mol = KEYS[key]
    L = {k: per_band(legs, "vis", k) for k in ("nominal", "xi_minus", "xi_plus", "q80", "q97", "contref", "marcs",
                                               "cscale", "win60", "c_minus", "c_plus", "n_minus",
                                               "n_plus", "o_minus", "o_plus", "ni_minus", "ni_plus")}
    nomdf = L["nominal"]
    r0 = nomdf[nomdf.key == key].iloc[0]
    A0 = float(r0.A_X)
    diag = next(d for d in REGION_DIAGNOSTICS["vis"] if d.key == key)
    win = list(diag.windows_A)
    ind = f"{key}:vis:{HOLDING}:windows_sha256:" + hashlib.sha256(json.dumps(win).encode()).hexdigest()[:24]
    ids = [ind]
    digest = uc.pool_digest(ids)
    comps, notes = [], []
    comps.append(uc.fit_measurement(
        float(r0.sigma_fit), source="vis nominal joint fit, pipeline.fit_constraint curvature",
        likelihood="chi2 curvature of the single-abundance fit, rescaled to red_chi2 = 1",
        correlation_treatment="pixels independent after the red_chi2 rescale; residual covariance not modelled",
        n_pixels=int(r0.n_pix)))
    # transition data
    if key in ("CI_5052", "CI_5380", "OI_6300"):
        s, src = gf_sigma("C I" if element == "C" else "O I", win[0][0], win[-1][1])
        if s is None:
            notes.append("transition_data: no canonical_gf sigma")
        else:
            comps.append(uc.transition_data(ids, [s], [1.0], covariance=[[s * s]], sources=[src],
                                            covariance_source=src))
    elif key == "CN_red":
        comps.append(uc.transition_data(ids, [CN_TD], [1.0], covariance=[[CN_TD ** 2]],
                                        sources=[CN_TD_SOURCE], covariance_source=CN_TD_SOURCE))
    else:
        a_ci, s_ci = ci_reference()
        if a_ci is None:
            notes.append("transition_data: no budgeted C I reference product")
        else:
            bound = math.hypot(A0 - a_ci, s_ci)
            src = (f"CONSISTENCY BOUND on the {mol} band-strength scale: |A({mol}) - A(C I)| (+) "
                   f"sigma(C I), A(C I) = {a_ci} +/- {s_ci:.4f} (HARPS AGSS21 C I 1D-LTE, same "
                   f"spectrum and model class, NIST-graded gf). Laboratory lifetime sources "
                   f"(Luque & Crosley 1996 JCP 104 2146; Brooke et al. 2013 JQSRT 124 11) not "
                   f"obtainable; a bound, never a lab value")
            comps.append(uc.transition_data(ids, [bound], [1.0], covariance=[[bound * bound]],
                                            sources=[src], covariance_source=src))
    for name, why in (("stellar.logg", "solar log g is pinned (RYA-1089)"),
                      ("stellar.metallicity", "solar [Fe/H] is pinned (RYA-1089)")):
        comps.append(dict(name=name, sigma_dex=0.0, state="DEFINED", source=why,
                          evidence={"delta": 0.0, "parameter_exists": True}))
    comps.append(dict(name="stellar.teff", sigma_dex=0.0, state="DEFINED",
                      source="sourced DEFINED zero (RYA-1226 convention): solar Teff allowance 1 K",
                      evidence={"parameter_exists": True, "delta_K": 1.0}))
    am, ap = val(L["xi_minus"], key), val(L["xi_plus"], key)
    if am is not None and ap is not None:
        s = DELTA_XI / XI_STEP
        dp, dm = (ap - A0) * s, (am - A0) * s
        signed = (dp - dm) / 2
        comps.append(dict(name="stellar.xi", sigma_dex=abs(signed), state="MEASURED",
                          source="paired xi 0.90/1.10 refits, stamped legs",
                          evidence={"pool_sha256": digest, "parameter_source": "RYA-1089 delta_xi 0.2912 km/s",
                                    "response_assessment": f"A(0.9)={am}, A(1.0)={A0}, A(1.1)={ap}",
                                    "signed_response_dex": signed, "delta_plus_dex": dp,
                                    "delta_minus_dex": dm, "delta_parameter": DELTA_XI, "parameter_exists": True}))
    else:
        notes.append("xi legs missing")
    q80, q97 = val(L["q80"], key), val(L["q97"], key)
    cref, nom = val(L["contref"], key), val(L["nominal"], key)
    if q80 is not None and q97 is not None and cref is not None and nom is not None:
        comps.append(dict(name="continuum", sigma_dex=float((((q97 - q80) / 2) ** 2 + (cref - nom) ** 2) ** 0.5),
                          state="MEASURED",
                          source=("continuum (synthesis reference, 3-MAD clipped): placement q97 vs q80 "
                                  "half-difference + REFERENCE (IAG-atlas leg minus nominal), in quadrature"),
                          evidence={"pool_sha256": digest, "A_q80": q80, "A_q97": q97, "A_contref_iag": cref}))
    else:
        notes.append("continuum legs missing")
    w60 = val(L["win60"], key)
    if w60 is not None:
        comps.append(dict(name="profile_ew", sigma_dex=abs(w60 - A0), state="MEASURED",
                          source="fit windows x0.6 about their centres (continuum on unscaled windows)",
                          evidence={"pool_sha256": digest, "A_win60": w60}))
    else:
        notes.append("win60 leg missing")
    am_ = val(L["marcs"], key)
    if am_ is not None:
        v, ev, st = abs(am_ - A0), {"pool_sha256": digest, "A_marcs": am_}, "MEASURED"
        mf = model_form(mol) if mol else None
        if mf is not None and mf > v:
            ev.update(bound=True, measured_leg_dex=v, bound_dex=mf,
                      bound_source=f"Amarsi et al. 2021 3D - MARCS shift for {mol} (RYA-1220 table)")
            v, st = mf, "DEFINED"
        comps.append(dict(name="model_atmosphere", sigma_dex=v, state=st,
                          source="MARCS.GES vs ATLAS9.Castelli" + ("; bounded by the published 1D->3D shift" if st == "DEFINED" else ""),
                          evidence=ev))
    else:
        notes.append("marcs leg missing")
    acs = val(L["cscale"], key)
    if acs is not None:
        dadf = (acs - A0) / (-CSCALE)
        res = []
        for lo, hi in win:
            ev = B.telluric_line(INSTRUMENT, HOLDING, 0.5 * (lo + hi), 0.5 * (hi - lo), 3780.0, 6910.0)
            res.append({"window_A": [lo, hi], **{k: ev[k] for k in ev if k in (
                "n_telluric_px", "max_depth", "residual_flux", "residual_basis", "sky_absorption")}})
        if any(r["residual_flux"] is None for r in res):
            notes.append(f"telluric: window without an independent reference: {res}")
        else:
            n = len(res)
            comps.append(dict(name="telluric", sigma_dex=math.sqrt(sum((abs(dadf) * r["residual_flux"] / n) ** 2 for r in res)),
                              state="MEASURED", source="per window measured residual vs an independent correction, via dA/df",
                              evidence={"pool_sha256": digest, "dA_df": dadf, "per_window": res}))
    else:
        notes.append("cscale leg missing")
    # coupling to the OTHER elements (molecules) / blends
    partners = {"CH_Gband": ("N", "O"), "C2_Swan": ("N", "O"), "CN_red": ("C", "O")}.get(key)
    if partners:
        J, resp = [], {}
        for el in partners:
            lo_, hi_ = val(L[f"{el.lower()}_minus"], key), val(L[f"{el.lower()}_plus"], key)
            if lo_ is None or hi_ is None:
                notes.append(f"coupling: {el} pin legs missing"); J = None; break
            J.append((hi_ - lo_) / (2 * STEP)); resp[el] = {"minus": lo_, "plus": hi_}
        if J is not None:
            s0, s1 = SIGMA[partners[0]], SIGMA[partners[1]]
            rho = 1.0 if J[0] * J[1] >= 0 else -1.0
            cov = np.array([[s0 * s0, rho * s0 * s1], [rho * s0 * s1, s1 * s1]])
            sig = float(math.sqrt(max(np.array(J) @ cov @ np.array(J), 0.0)))
            comps.append(dict(name="molecular_coupling", sigma_dex=sig, state="MEASURED",
                              source=(f"paired pins of A({partners[0]}), A({partners[1]}) +/- 0.10 on this joint fit x "
                                      f"our published budgets, worst rho = {rho:+.0f} (RYA-1220 policy)"),
                              evidence={"pool_sha256": digest, "jacobian": dict(zip(partners, J)),
                                        "sigmas": {partners[0]: s0, partners[1]: s1}, "responses": resp}))
    else:
        comps.append(dict(name="molecular_coupling", sigma_dex=None, state="N/A",
                          source="atomic/forbidden product; the selector is not a MOL- indicator set",
                          evidence={"selector": selector}))
    if key == "OI_6300":
        lo_, hi_ = val(L["ni_minus"], key), val(L["ni_plus"], key)
        if lo_ is None or hi_ is None:
            notes.append("blends: Ni pin legs missing")
        else:
            jn = (hi_ - lo_) / (2 * STEP)
            comps.append(dict(name="blends", sigma_dex=abs(jn) * SIGMA_NI, state="MEASURED",
                              source="Ni I 6300.34 blend in [O I] 6300: paired A(Ni) +/- 0.10 x sigma(Ni)",
                              evidence={"pool_sha256": digest, "dAO_dANi": jn, "sigma_Ni": SIGMA_NI,
                                        "sigma_Ni_source": SIGMA_NI_SOURCE}))
    else:
        comps.append(dict(name="blends", sigma_dex=None, state="N/A",
                          source="full atomic + molecular synthesis in every fit window", evidence={}))
    vals = {t: twins_A[t].get(key) for t in TWINS if twins_A.get(t) and twins_A[t].get(key) is not None}
    if len(vals) >= 2:
        comps.append(dict(name="holding_instrument", sigma_dex=float(np.std(list(vals.values()), ddof=1)),
                          state="MEASURED", source="std of this diagnostic across the vis twin regions (HARPS, KP x2, IAG)",
                          evidence={"per_holding_A": vals}))
    else:
        notes.append(f"holding_instrument: only {len(vals)} holding(s)")
    for name, why in (("nlte", "vis region is LTE by design (lte_by_design); the product's treatment is 1D-LTE"),
                      ("pseudo_continuum", "continuum placed per window by the model-guided rule; priced on `continuum`"),
                      ("hfs_isotopes", "no HFS in C I/O I/CH/C2/CN lines here; isotopologues at solar ratios")):
        comps.append(dict(name=name, sigma_dex=None, state="N/A", source=why, evidence={}))
    # n_lines: a molecular band counts its OWN isotopologue's lines inside the fit windows
    # (the RYA-1214/1230 recount: CH 184, C2 167, CN red 194); an atomic/forbidden
    # diagnostic is one line. RYA-1232 fix: this was hard-coded to 1 for every diagnostic.
    if mol:
        from rya1230_n_product_hygiene import _bsyn_species_in
        iso = {"CH": "12CH", "C2": "12C12C", "CN": "12C14N"}[mol]
        n_lines = int(_bsyn_species_in(win).get(iso, 0))
    else:
        n_lines = 1
    row = dict(element=element, ion="I", band="VIS", instrument=INSTRUMENT, treatment="1D-LTE",
               handler="CNOSynthesis", A=A0, n_lines=n_lines, n_excluded=0, stat_dex=float(r0.sigma_fit),
               syst_dex=np.nan, stat_basis=("measured -- 1 sigma from the chi2 curvature of THIS fit "
                                            "(pipeline.fit_constraint), rescaled to red_chi2 = 1"),
               dominant="", route="synth", scale="1D-LTE", model="none", atmos="atlas9", gf="", route_basis="handler", deck="none")
    from publish_product import normalise
    pub = normalise(pd.DataFrame([row]), holding=HOLDING, tier="ALL", route="SYNTH", selector=selector)[0]
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
    return row, doc, ids, verdict, notes, selector


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--legs", type=Path, required=True)
    ap.add_argument("--stage", type=Path, required=True)
    ap.add_argument("--report", type=Path, required=True)
    a = ap.parse_args()
    twins_A = {}
    for t in TWINS:
        df = per_band(a.legs, t, "nominal")
        twins_A[t] = {k: val(df, k) for k in KEYS} if df is not None else None
    a.stage.mkdir(parents=True, exist_ok=True)
    out = []
    for key in KEYS:
        row, b, ids, verdict, notes, selector = build(key, a.legs, twins_A)
        rec = {"key": key, "selector": selector, "element": row["element"], "A": row["A"],
               "verdict": verdict, "notes": notes, "sigma_reported": b["sigma_reported"],
               "components": {c["name"]: (c["state"], c["sigma_dex"]) for c in b["components"]}}
        if verdict is None:
            r = dict(row)
            r["syst_dex"] = round(math.sqrt(max(b["sigma_reported"] ** 2 - r["stat_dex"] ** 2, 0.0)), 4)
            stem = f"{row['element']}I_{selector.replace('-', '_', 1)}_{INSTRUMENT}_{HOLDING}_SYNTH_{selector}_products.csv"
            df = pd.DataFrame([r])
            df["uncertainty"] = json.dumps(b)
            df["uncertainty_indicator_ids"] = json.dumps(ids)
            df["sigma_reported"] = b["sigma_reported"]
            df.to_csv(a.stage / stem, index=False)
            rec["prod_stem"] = stem
        rec.setdefault("holding", HOLDING)
        rec.setdefault("key_treatment", "1D-LTE")
        out.append(rec)
        print(f"{key:9} {selector:13} A={row['A']}  " + (f"sigma_reported={b['sigma_reported']:.4f}"
              if verdict is None else f"HELD: {verdict} {notes}"))
    a.report.parent.mkdir(parents=True, exist_ok=True)
    a.report.write_text(json.dumps({"ticket": "RYA-1230", "twins_A": twins_A, "products": out},
                                   indent=1, default=str) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
