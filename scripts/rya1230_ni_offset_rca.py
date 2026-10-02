"""RYA-1230 -- decompose the solar atomic N I offset, line by line. DIAGNOSTIC ONLY.

Our published 1D-LTE atomic N I (N.json v1.14: IAG 8.351, KP 8.293 / 8.287) sits ~+0.5
dex above Asplund, Amarsi & Grevesse 2021's own 1D-LTE MARCS N I (7.813) on the same
lines. 3D and NLTE are ~-0.05 and ~-0.01, so the offset is in the INPUTS. This script
localises it by re-running the PRODUCTION fit (the RYA-759 synthesis route:
`pipeline.abundances_derive._fit_synth_flux` on the `build_solar_context` context that
`derive_band_products` builds) and then changing ONE input at a time.

THE LADDER (cumulative; each rung differs from the one before it by one input):

  S0 PROD   production: ATLAS9, GES v6 atomic list + canonical gf, NO molecules,
            +/-1.1 A window, observed flux taken as already normalised (continuum = 1)
  S1 CORE   window +/-1.1 A -> +/-0.25 A (the N I core; the lines are 2-8 mA)
  S2 MOL    + CN/CH/C2 molecular opacity (the RYA-236/RYA-360 vendored Turbospectrum
            lists, byte-identical to the iSpec bundle the synthesis reads)
  S3 CONT   + local continuum: observed / p95(observed within +/-CONT_HW_A), the
            RYA-1000 `--local-renorm` / verify_feature estimator, reused not invented.
            The lines are 0.3-1% deep, so the continuum is not a detail here: a second
            estimator (median obs/model over the model's own continuum pixels) and the
            slope dA/dcontinuum are recorded beside it, and the continuum term is reported
            as the BRACKET the two estimators span, never as one number.
  S4 MARCS  ATLAS9.Castelli -> MARCS.GES (Asplund's 1D column is MARCS)

Every lever is ALSO measured alone from S0 ("alone") so an order-dependent attribution is
visible rather than hidden by the ladder order. gf is not refitted: our log gf and the
published Amarsi 2020 / Tachiev & Froese Fischer 2002 value are compared directly and the
difference is the gf term (a log gf offset moves A by minus the same amount).

red_chi2 IS NOT A QUALITY MEASURE HERE: `_fit_synth_flux` is called directly, without the
SynthesisHandler's edge trim, so the one or two observed pixels that land on iSpec's zeroed
synthesis edge add a CONSTANT to chi2 (it does not move with A(N), so the fitted A is
unaffected -- S0 reproduces production to 0.005 dex). Measured on the IAG S4 core fits with
the edge excluded: RMS residual 1.1-1.4% on 7442/8629/8683, 9.9% on 8216.

RYA-161 VALIDATE-DON'T-TUNE: nothing here changes a gf, a blend strength, a line list, or
a line selection. 7.813 is a CHECK column in the output, never an input.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

#: Asplund, Amarsi & Grevesse 2021 (A&A 653, A141) Table 3, N I, 1D LTE MARCS. CHECK only.
A21_N_I_1D_LTE_MARCS = 7.813
#: AGSS21 five-line set in air A with their published log gf (Amarsi et al. 2020, A&A 636,
#: A120, Tables 1-3, as transcribed by RYA-1220 into audit_v1/ni_five_line_reference.csv).
REFERENCE_CSV = REPO / "data/output/rya1220/audit_v1/ni_five_line_reference.csv"
LINES = (7442.298, 8216.336, 8629.235, 8683.403)
HOLDINGS = (("iag_fts_solar_atlas", "solar_iag"),
            ("kpno_solar_atlas", "solar_kpno_kurucz2005_corrected"),
            ("kpno_solar_atlas", "solar_kpno_molecfit_corrected"))
PROD_HW_A = 1.1          # config/synth_bands.yaml red-optical half_width_A
CORE_HW_A = 0.25
CONT_HW_A = 1.5
CONT_MODEL_MIN = 0.998
A_LO, A_HI = 6.3, 9.5
R_PROD = 700000.0        # what derive_band_products printed for this route (repro log)


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=REPO / "data/output/rya1230")
    ap.add_argument("--holdings", default="all",
                    help="comma list of holding ids, or 'all'")
    ap.add_argument("--lines", default="all")
    ap.add_argument("--blend-only", type=Path, nargs="+", default=None,
                    help="recompute only the blend/EW fields of existing ladder JSON(s)")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    from pipeline.abundances_derive import (_fit_synth_flux, _synth_flux_at_abund,
                                            _load_atmosphere, _MARCS, _ATLAS9)
    from pipeline.nearuv_synth import build_solar_context, gf_provenance
    from config.synth_bands import SYNTH_BANDS
    from measure_band_ew import load_window_ex

    cfg = SYNTH_BANDS["red-optical"]
    gf = gf_provenance(6910.0, 9199.0)
    ctx = build_solar_context("N", R_PROD, linelist_file=str(cfg.linelist),
                              apply_canonical_gf=gf["apply_canonical_gf"])
    atm_marcs = _load_atmosphere(ctx["teff"], ctx["logg"], ctx["feh"], ctx["vturb"],
                                 model_grid=_MARCS)
    ll = ctx["linelist"]
    import tempfile
    tmp = tempfile.mkdtemp(prefix="rya1230_ispec_")   # never inside the output tree

    ref = {round(float(r["wavelength_air_A"]), 1): r
           for r in csv.DictReader(REFERENCE_CSV.open())}

    holdings = HOLDINGS if args.holdings == "all" else tuple(
        h for h in HOLDINGS if h[1] in args.holdings.split(","))
    lines = LINES if args.lines == "all" else tuple(float(x) for x in args.lines.split(","))

    def kw(atm, mol):
        return dict(atmosphere=atm, teff=ctx["teff"], logg=ctx["logg"], feh=ctx["feh"],
                    vturb=ctx["vturb"], linelist=ll, isotopes=ctx["isotopes"],
                    solar_abund=ctx["solar_abund"], element="N",
                    atom_code=ctx["atom_code"], R=ctx["resolving_power"],
                    macroturbulence=ctx["macroturbulence"], vsini=ctx["vsini"],
                    use_molecules=mol, tmp_dir=tmp)

    def synth(w_nm, A, atm, mol):
        return _synth_flux_at_abund(w_nm, trial_A=A, **kw(atm, mol))

    def fit(ow_nm, of, c, hw, atm, mol):
        r = _fit_synth_flux(ow_nm, of, atm, ctx["teff"], ctx["logg"], ctx["feh"],
                            ctx["vturb"], ll, ctx["isotopes"], ctx["solar_abund"], "N",
                            ctx["atom_code"], (c - hw) / 10.0, (c + hw) / 10.0,
                            A_LO, A_HI, ctx["resolving_power"], ctx["macroturbulence"],
                            ctx["vsini"], use_molecules=mol, tmp_dir=tmp)
        return dict(A=float(r["A_X"]), red_chi2=float(r["red_chi2"]),
                    n_pix=int(r["n_pix"]), status=r["status"])

    def p95_continuum(ow_nm, of, c):
        m = np.abs(ow_nm * 10.0 - c) <= CONT_HW_A
        return float(np.nanpercentile(of[m], 95))

    def local_continuum(ow_nm, of, c, A, atm, mol):
        """median(obs/model) over this window's model-continuum pixels."""
        sw = np.arange((c - CONT_HW_A) / 10.0, (c + CONT_HW_A) / 10.0, 0.0002)
        sf = synth(sw, A, atm, mol)
        m = np.abs(ow_nm * 10.0 - c) <= CONT_HW_A
        mod = np.interp(ow_nm[m], sw, sf)
        sel = mod >= CONT_MODEL_MIN
        if sel.sum() < 10:
            return float("nan"), int(sel.sum())
        return float(np.median(of[m][sel] / mod[sel])), int(sel.sum())

    #: Every EW below is integrated over the line CORE of a WIDE synthesis. A narrow
    #: synthesis is unusable: iSpec zeroes the edge pixels of every synthesis, and at
    #: +/-0.1 A those edges ARE the window.
    SYN_HW_A, CORE_EW_HW_A = 3.0, 0.25

    def _ew(c, depth, sw):
        m = np.abs(sw * 10.0 - c) <= CORE_EW_HW_A
        return float(np.trapz(depth[m], sw[m] * 10.0) * 1000.0)

    def ni_ew_mA(c, A, atm):
        """N-I-only synthetic flux EW. Molecules OFF on both sides: with them on, A(N)=3
        also deletes every CN line, and the 'N I' EW silently becomes N I + CN."""
        sw = np.arange((c - SYN_HW_A) / 10.0, (c + SYN_HW_A) / 10.0, 0.0002)
        return _ew(c, synth(sw, 3.0, atm, False) - synth(sw, A, atm, False), sw)

    def blend_fraction(c, A, atm):
        """Composition of the model feature in the +/-0.25 A core AT FIXED A(N).

        molecular = (molecules on) - (molecules off) at the SAME A(N) -- the CN/CH/C2
        absorption the production synthesis omits; N I = (N at A) - (N removed), molecules
        off; other atomic = what is left with N removed and molecules off."""
        sw = np.arange((c - SYN_HW_A) / 10.0, (c + SYN_HW_A) / 10.0, 0.0002)
        full = 1.0 - synth(sw, A, atm, True)
        nomol = 1.0 - synth(sw, A, atm, False)
        noN_nomol = 1.0 - synth(sw, 3.0, atm, False)
        f, mol, ni = _ew(c, full, sw), _ew(c, full - nomol, sw), _ew(c, nomol - noN_nomol, sw)
        return dict(feature_ew_mA=f, molecular_ew_mA=mol, ni_ew_mA=ni,
                    atomic_other_ew_mA=_ew(c, noN_nomol, sw),
                    molecular_fraction=mol / f, ni_fraction=ni / f,
                    core_half_width_A=CORE_EW_HW_A)

    if args.blend_only:
        # Recompute ONLY the blend/EW fields on an existing ladder (no refit).
        for f in args.blend_only:
            doc = json.loads(f.read_text())
            for rec in doc["rows"]:
                c, a4 = rec["wavelength_air_A"], rec["S4_MARCS"]["A"]
                rec["blend_at_fit"] = blend_fraction(c, a4, atm_marcs)
                rec["ni_ew_model_at_A21_mA"] = ni_ew_mA(c, A21_N_I_1D_LTE_MARCS, atm_marcs)
                rec["ni_ew_model_at_fit_mA"] = ni_ew_mA(c, a4, atm_marcs)
                print(rec["holding"], c, json.dumps(rec["blend_at_fit"], default=float), flush=True)
            doc["blend_fields_recomputed"] = "RYA-1230 --blend-only (fixed-A(N) molecular term)"
            f.write_text(json.dumps(doc, indent=1, default=float) + "\n")
        return 0

    rows = []
    t0 = time.time()
    for inst, hold in holdings:
        for c in lines:
            win = load_window_ex(inst, c, CONT_HW_A + 0.5, holding=hold)
            ow = np.asarray(win.wave, float) / 10.0
            of = np.asarray(win.flux, float)
            ok = np.isfinite(ow) & np.isfinite(of)
            ow, of = ow[ok], of[ok]
            A_ = ctx["atmosphere"]
            rec = dict(holding=hold, instrument=inst, wavelength_air_A=c,
                       provenance=str(win.provenance))
            # ladder
            s0 = fit(ow, of, c, PROD_HW_A, A_, False)
            s1 = fit(ow, of, c, CORE_HW_A, A_, False)
            s2 = fit(ow, of, c, CORE_HW_A, A_, True)
            k3 = p95_continuum(ow, of, c)
            s3 = fit(ow, of / k3, c, CORE_HW_A, A_, True)
            s4 = fit(ow, of / k3, c, CORE_HW_A, atm_marcs, True)
            # the continuum bracket, on the S2 state
            k_ratio, n_ratio = local_continuum(ow, of, c, s2["A"], A_, True)
            s3_ratio = fit(ow, of / k_ratio, c, CORE_HW_A, A_, True)
            up = fit(ow, of / 1.001, c, CORE_HW_A, A_, True)
            dn = fit(ow, of / 0.999, c, CORE_HW_A, A_, True)
            # levers alone from S0
            a_core = s1
            a_mol = fit(ow, of, c, PROD_HW_A, A_, True)
            a_cont = fit(ow, of / k3, c, PROD_HW_A, A_, False)
            a_marcs = fit(ow, of, c, PROD_HW_A, atm_marcs, False)
            # gf
            r = ref[round(c, 1)]
            ours = ll[(np.abs(ll["wave_A"] - c) < 0.01) & (ll["element"] == "N 1")]
            our_loggf = float(ours["loggf"][0]) if len(ours) else float("nan")
            pub_loggf = float(r["published_loggf"])
            # EW + blend, at the final state and at the reference abundance
            bf = blend_fraction(c, s4["A"], atm_marcs)
            rec.update(
                S0_PROD=s0, S1_CORE=s1, S2_MOL=s2, S3_CONT=s3, S4_MARCS=s4,
                continuum_p95=k3, continuum_model_ratio=k_ratio,
                continuum_model_ratio_npix=n_ratio, S3_CONT_model_ratio=s3_ratio,
                continuum_bracket_dex=sorted([s3["A"] - s2["A"], s3_ratio["A"] - s2["A"]]),
                dA_per_plus0p1pct_continuum=(up["A"] - dn["A"]) / 2.0,
                alone=dict(window=a_core["A"] - s0["A"], molecules=a_mol["A"] - s0["A"],
                           continuum=a_cont["A"] - s0["A"], marcs=a_marcs["A"] - s0["A"]),
                ladder=dict(window=s1["A"] - s0["A"], molecules=s2["A"] - s1["A"],
                            continuum=s3["A"] - s2["A"], marcs=s4["A"] - s3["A"]),
                our_loggf=our_loggf, published_loggf_amarsi2020=pub_loggf,
                # a higher gf needs less N: our A minus theirs = -(ours - theirs)
                gf_term=-(our_loggf - pub_loggf),
                ni_ew_model_at_A21_mA=ni_ew_mA(c, A21_N_I_1D_LTE_MARCS, atm_marcs),
                ni_ew_model_at_fit_mA=ni_ew_mA(c, s4["A"], atm_marcs),
                amarsi2020_ni_ew_intensity_mA=float(r["ni_ew_mA"]),
                amarsi2020_feature_ew_intensity_mA=float(r["feature_ew_pm"]) * 10.0,
                amarsi2020_1d_lte_intensity_A=float(r["A_1d_lte"]),
                blend_at_fit=bf,
                residual_vs_A21=s4["A"] - A21_N_I_1D_LTE_MARCS,
                total_vs_A21=s0["A"] - A21_N_I_1D_LTE_MARCS)
            rows.append(rec)
            print(json.dumps({k: rec[k] for k in ("holding", "wavelength_air_A", "ladder",
                                                  "alone", "residual_vs_A21")},
                             default=float), flush=True)
            print(f"   S0 {s0['A']:.3f} S1 {s1['A']:.3f} S2 {s2['A']:.3f} "
                  f"S3 {s3['A']:.3f} (p95={k3:.5f}; ratio {k_ratio:.5f}->{s3_ratio['A']:.3f}) "
                  f"S4 {s4['A']:.3f} dA/dc(+0.1%)={(up['A'] - dn['A']) / 2:+.3f}  "
                  f"[{time.time() - t0:.0f}s]", flush=True)
            (args.out / "ni_offset_ladder.json").write_text(
                json.dumps(dict(ticket="RYA-1230", rows=rows), indent=1, default=float) + "\n")

    prov = dict(
        ticket="RYA-1230", status="DIAGNOSTIC -- no product, gf, or line list changed",
        route="RYA-759 synthesis route (_fit_synth_flux) on build_solar_context('N')",
        R=R_PROD, prod_half_width_A=PROD_HW_A, core_half_width_A=CORE_HW_A,
        continuum_rule=f"S3: p95(observed) within +/-{CONT_HW_A} A (RYA-1000 convention); "
                       f"bracket partner: median(obs/model) over pixels with model>="
                       f"{CONT_MODEL_MIN} within +/-{CONT_HW_A} A",
        broadening=dict(vmac=ctx["macroturbulence"], vsini=ctx["vsini"], xi=ctx["vturb"]),
        atmospheres=[_ATLAS9, _MARCS], gf=gf,
        check_value=dict(A21_N_I_1D_LTE_MARCS=A21_N_I_1D_LTE_MARCS,
                         source="Asplund, Amarsi & Grevesse 2021, A&A 653, A141, Table 3"),
        sources={str(p.relative_to(REPO)): _sha(p) for p in (
            REPO / "pipeline/abundances_derive.py", REPO / "pipeline/nearuv_synth.py",
            REFERENCE_CSV, Path(__file__).resolve())},
        linelist="<ISPEC_DIR>/" + str(cfg.linelist).split("ispec_src/", 1)[-1])
    (args.out / "ni_offset_provenance.json").write_text(json.dumps(prov, indent=1, default=str) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
