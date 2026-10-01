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
  nlte                N/A on a 1D-LTE leg; on a departure leg RYA-1032's measured
                      model-family spread (0.043), FLAGGED cross-element (RYA-1226 precedent)
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
SPECIES = {"C": "C I", "N": "N I", "O": "O I"}


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
        "harps": ("kpno_solar_atlas", "solar_kpno", "solar_kpno_kurucz2005_corrected")}
#: INDEPENDENT corrections to difference each holding against on the telluric pixels, in
#: order. A candidate byte-identical to raw in THIS window is not a correction there (KP
#: molecfit leaves the 9000-9300 A H2O band untouched, RYA-1190) and is skipped -- using it
#: measured the whole sky as kurucz2005's "residual" on C I 9061 (0.33 flux, 5 dex).
_INDEPENDENT = {"solar_kpno_molecfit_corrected": [("kpno_solar_atlas", "solar_kpno_kurucz2005_corrected"),
                                                  ("iag_fts_solar_atlas", "solar_iag")],
                "solar_kpno_kurucz2005_corrected": [("kpno_solar_atlas", "solar_kpno_molecfit_corrected"),
                                                    ("iag_fts_solar_atlas", "solar_iag")],
                "solar_iag": [("kpno_solar_atlas", "solar_kpno_kurucz2005_corrected"),
                              ("kpno_solar_atlas", "solar_kpno_molecfit_corrected")],
                "solar_harps_molecfit_corrected": [("kpno_solar_atlas", "solar_kpno_kurucz2005_corrected"),
                                                   ("iag_fts_solar_atlas", "solar_iag")]}
_EDGE: dict = {}


def telluric_line(instrument, holding, w0, hw, band_lo, band_hi) -> dict:
    import measure_band_ew as H
    from pipeline import telluric_observability as T
    sky_inst, raw, cor = _SKY[instrument]
    key = (raw, band_lo, band_hi)
    if key not in _EDGE:
        snr, _ = T.band_continuum_snr(H, sky_inst, raw, band_lo, band_hi, max(hw, 2.0))
        _EDGE[key] = T.thresholds(snr).clean_max_depth
    edge = _EDGE[key]
    g = np.linspace(w0 - hw, w0 + hw, 400)

    def flux(inst, hold):
        x = H.load_window_ex(inst, w0, hw + 0.5, holding=hold, allow_uncorrected=True)
        return np.interp(g, x.wave, x.flux)

    t = flux(sky_inst, raw) / flux(sky_inst, cor)
    base = float(np.median(t))
    P = (base - t) > edge
    ev = {"sky_pair": f"{raw}/{cor}", "clean_edge": edge, "sky_baseline": base,
          "n_px": int(g.size), "n_telluric_px": int(P.sum()),
          "max_depth": float(max(0.0, base - t.min()))}
    if not P.any():
        ev["residual_flux"] = 0.0
        return ev
    own = flux(instrument, holding)
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
        raise RuntimeError(f"no independent correction covers {w0:.3f} A for {holding}")
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
    element = str(prod["element"].iloc[0])
    holding = _arg(unit_args, "--holding")
    instrument = _arg(unit_args, "--instrument")
    selector = _selector(unit_args)
    row = normalise(prod, holding=holding, tier="ALL", route="SYNTH", selector=selector)[0]
    row["star"] = "solar"
    band = row["band"]
    hw = float(_arg(unit_args, "--half-width-A", SYNTH_BANDS[band].half_width_A))
    lo, hi = float(_arg(unit_args, "--lo")), float(_arg(unit_args, "--hi"))
    notes, comps = [], []
    if n < 2:
        return {"row": row, "skip": f"n_lines={n}: a single-line pool has no line-scatter "
                                    f"measurement; RYA-587 keeps it HOLD"}

    # transition data: lambda AND EP join
    g = gf[gf["species"] == SPECIES[element]]
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
        p = _paired(_acc(leg), acc, n)
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
        pc = _paired(_acc(hi_l), _acc(lo_l), n)
        if pc["moved"]:
            notes.append(f"continuum: pool moved ({pc['n_paired']} vs {n})")
        else:
            comps.append(dict(name="continuum", sigma_dex=abs(pc["median"]) / 2.0, state="MEASURED",
                              source=("standing model-guided continuum rule, pixel-selection "
                                      "quantile q97 vs q80 around the nominal q90, central "
                                      "half-difference, paired on this pool"),
                              evidence={"pool_sha256": digest, "paired_median_q97_minus_q80": pc["median"],
                                        "n_paired": pc["n_paired"],
                                        "withdrawn_leg": ("+/-1.5 A envelope: < 5 bins, went "
                                                          "CONTINUUM_UNCONSTRAINED and measured "
                                                          "the whole correction, not placement")}))
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
            pp = _paired(_acc(cp), _acc(cm), n)
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

    nd = set(nom.loc[nom.in_aggregate == True, "nlte_delta_dex"].dropna().round(6))  # noqa: E712
    if nd <= {0.0}:
        comps.append(dict(name="nlte", sigma_dex=None, state="N/A",
                          source="this leg applies no departures (nlte_delta_dex = 0 on every line)",
                          evidence={"nlte_delta_dex": 0.0}))
    else:
        comps.append(dict(name="nlte", sigma_dex=0.043, state="DEFINED",
                          source="RYA-1032 model_family_spread (Gerber - Bergemann, 1D-NLTE)",
                          evidence={"bound": True, "status": "measured on the VIS Fe I pool, carried "
                                    "cross-element (RYA-1226 precedent)",
                                    "owed_measurement": "two independent departure treatments on this pool"}))
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
    return {"row": row, "budget": doc, "ids": ids, "notes": notes, "verdict": verdict,
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
