#!/usr/bin/env python3
"""RYA-1190 — frontier-band uncatalogued opacity: what is resolvable, and what is not.

    python3 scripts/rya1190_frontier_opacity.py

DIAGNOSTIC ONLY. Identifies, classifies, counts and estimates. It writes nothing outside
`data/results/rya1190/`, adds no line to any production list, and changes no published
value (RYA-161).

WHAT IS MEASURED, AND WHY IT IS NOT A SYNTHESIS DIFFERENTIAL
--------------------------------------------------------------
The ticket asks for "the residual between the observed side-band absorption and the
current synthetic built from linelist_solar". A true synthesis differential needs the
production engine, and **Turbospectrum's `bsyn_lu` is absent from this Mac while MOOG is a
Linux ELF** — `spectrum` is the only arm64 binary present, and it is NOT the production
engine (RYA-289 is precisely about not substituting it silently). So instead of a
different engine's synthetic, this computes a CATALOGUED-OPACITY ACCOUNTING from the same
line list the synthesis reads:

    tau0_i = -ln(1 - central_depth_i)        so one line reproduces its own depth exactly
    predicted flux = exp( -SUM_i tau0_i * Gauss(lambda - lambda_i, sigma) )

The exponential of the SUM (rather than a sum of depths) is what makes overlap and
saturation behave: a hundred weak lines on top of each other stop adding linearly, which
is exactly the regime the near-UV is in.

⚠️ IT IS AN ACCOUNTING, NOT RADIATIVE TRANSFER. No atmosphere, no ionisation balance, no
proper source function. It answers "how much absorption does our line list ACCOUNT FOR
here", which is the question the classification needs, and it must not be quoted as a
synthetic spectrum. `sigma` is set physically (thermal + microturbulent Doppler, scaling
with lambda) and swept x0.6/x2.0, because it CANNOT be calibrated on the VIS control --
the control's side-bands are line-free by selection, so every sigma predicts 1.0000 there
and the fit is degenerate. A number that survives a x3.3 sigma sweep is not a sigma
artefact; that is the whole reason the sweep is reported beside it.

THE CONTROL CARRIES A FLOOR, AND IT IS SUBTRACTED
--------------------------------------------------
VIS is line-free in its side-bands and still sits ~1.1% below the accounting. That floor
is either the method's or a real small KP continuum offset, and this does not pretend to
know which -- it is subtracted, so every band is quoted as EXCESS over the control rather
than as an absolute deficit.

FIREWALL (RYA-161): opacity is identified from vendored catalogues with cited gf. Nothing
is fitted, tuned, or adjusted to move Fe toward 7.466.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

RYA1189 = ROOT / "data/results/rya1189/rya1189_per_line.csv"
LINELIST = ROOT / "data/linelists/linelist_solar.csv"
MOLDIR = ROOT / "data/linelists/molecular/turbospectrum"

#: The vendored VALD3 "extract all" pulls, per band. These are the lists the ticket's
#: "pull a full-species line list" step asks for -- already on disk, so the completeness
#: question is answerable offline and reproducibly.
VALD = {
    "near-UV":     ("data/linelists/vald_solar_nearuv_2000_3780_hfson_raw.txt", 3000.0, 3780.0),
    "VIS":         ("data/linelists/vald_solar_raw.txt",                        4200.0, 6910.0),
    "red-optical": ("data/linelists/vald_solar_redopt_6910_9500_hfson_raw.txt", 6910.0, 9199.0),
    "NIR-H":       ("data/linelists/vald_solar_ir_9500_17000_hfson_raw.txt",   15007.0, 17000.0),
}
HOLDING = {
    "near-UV":     ("kpno_solar_atlas", "solar_kpno_molecfit_corrected", "solar_kpno"),
    "VIS":         ("kpno_solar_atlas", "solar_kpno_molecfit_corrected", "solar_kpno"),
    "red-optical": ("kpno_solar_atlas", "solar_kpno_molecfit_corrected", "solar_kpno"),
    "NIR-H":       ("crires_plus",      "solar_crires_plus_h_rya1094",   None),
}
CONTROL = "VIS"
HALF, SB_OUT = 0.18, 0.36          # the RYA-1189 window, so the two tickets compare
DOPPLER_KMS = 1.6                  # thermal(Fe, 5772K) ~1.28 + microturbulence 1.0, quadrature
C_KMS = 2.998e5

#: 🔴 THE WIDTH MUST CARRY THE INSTRUMENT, AND OMITTING IT NEARLY COST THE NIR VERDICT.
#: A first cut used the Doppler term alone for every band. That is fine at Kitt Peak --
#: R 300-500k puts sigma_inst at 0.004 A against a 0.018 A Doppler width, a 2% correction
#: -- and badly wrong for CRIRES+, where R 50-100k gives sigma_inst 0.091 A against a
#: 0.085 A Doppler width. Under-broadening a prediction makes the catalogue look like it
#: absorbs LESS than it does, which manufactures an "uncatalogued deficit" out of nothing.
#: NIR-H's apparent deficit moved from -0.071 to -0.013 across the sigma sweep, which is
#: what made the omission visible; it is now computed rather than swept over.
#:
#: Resolving powers are read from data/catalog/instrument_catalog.csv (geometric mean of
#: the declared min/max), not typed in here.
def _sigma_A(centre: float, instrument: str, resolving: dict) -> float:
    sig_dop = centre * DOPPLER_KMS / C_KMS
    R = resolving.get(instrument)
    if not R:
        return sig_dop
    fwhm_inst = centre / R
    sig_inst = fwhm_inst / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return float(np.hypot(sig_dop, sig_inst))


def _resolving_powers() -> dict:
    d = pd.read_csv(ROOT / "data/catalog/instrument_catalog.csv")
    out = {}
    for _, r in d.iterrows():
        lo, hi = r.get("resolving_power_min"), r.get("resolving_power_max")
        if pd.notna(lo) and pd.notna(hi) and lo > 0 and hi > 0:
            out[str(r.instrument_id)] = float(np.sqrt(float(lo) * float(hi)))
    return out
OPTICAL_ANCHOR = 7.466


def _catalogue():
    d = pd.read_csv(LINELIST, low_memory=False)[["wavelength_air_A", "central_depth"]].dropna()
    w = d.wavelength_air_A.to_numpy()
    tau0 = -np.log(1.0 - np.clip(d.central_depth.to_numpy(), 0.0, 0.999))
    return w, tau0


def predict(wgrid, W0, TAU0, sigma, pad=1.0):
    m = (W0 > wgrid[0] - pad) & (W0 < wgrid[-1] + pad)
    tau = np.zeros_like(wgrid)
    for a, t in zip(W0[m], TAU0[m]):
        z = (wgrid - a) / sigma
        k = np.abs(z) < 6
        if k.any():
            tau[k] += t * np.exp(-0.5 * z[k] ** 2)
    return np.exp(-tau)


def vald_completeness(band: str) -> dict:
    """Is our production list already the VALD extract for this band?

    🔴 THE QUESTION THE TICKET'S PART A2 IS REALLY ASKING. "Pull a full-species VALD3
    extract and see what could account for the residual" presumes the production list is
    a SUBSET of what VALD offers. Measured rather than presumed.
    """
    from build_linelist import parse_vald
    rel, lo, hi = VALD[band]
    p = ROOT / rel
    if not p.exists():
        return {"available": False, "path": rel}
    v = pd.DataFrame(parse_vald(p))
    v = v[v.wavelength_air_A.between(lo, hi)]
    s = pd.read_csv(LINELIST, low_memory=False)
    s = s[s.wavelength_air_A.between(lo, hi)]

    def key(d):
        return set(zip(d.element.astype(str), d.ion.astype(str),
                       d.wavelength_air_A.round(4), d.log_gf.round(3)))
    kv, ks = key(v), key(s)
    depth = v.central_depth
    return {
        "available": True, "path": rel, "window_A": [lo, hi],
        "vald_rows": int(len(v)), "production_rows": int(len(s)),
        "in_vald_only": int(len(kv - ks)), "in_production_only": int(len(ks - kv)),
        "identical_keys": int(len(kv & ks)),
        "extraction_depth_threshold": (round(float(depth.min()), 6) if len(depth) else None),
        "verdict": ("production list IS the VALD extract — nothing to add at this "
                    "threshold" if not (kv - ks) else
                    f"{len(kv - ks)} VALD lines are absent from the production list"),
    }


def molecular_coverage() -> dict:
    """What the vendored Turbospectrum molecular lists actually span — read from the FILES.

    ⚠️ The filenames encode NANOMETRES (`12C14N_400-450.bsyn`), and the manifest declares
    a HARPS range starting at 3800 A. Neither is trusted here: the first column of each
    lowest-region file is read, because a unit misread would silently answer the near-UV
    question the wrong way in either direction.
    """
    out = {}
    for sub in sorted(p.name for p in MOLDIR.iterdir() if p.is_dir()):
        lo = hi = None
        n = 0
        for f in sorted((MOLDIR / sub).glob("*.bsyn")):
            for line in f.read_text(errors="replace").splitlines():
                t = line.strip()
                if not t or t.startswith(("'", "*")):
                    continue
                try:
                    w = float(t.split()[0])
                except (ValueError, IndexError):
                    continue
                n += 1
                lo = w if lo is None else min(lo, w)
                hi = w if hi is None else max(hi, w)
        if n:
            out[sub] = {"n_lines": n, "min_A": round(lo, 2), "max_A": round(hi, 2)}
    return out


def measure(band: str, W0, TAU0, resolving: dict) -> list[dict]:
    from measure_band_ew import load_window_ex
    inst, hold, raw = HOLDING[band]
    per = pd.read_csv(RYA1189)
    per = per[(per.band == band) & (per.status == "ok")]
    rows = []
    for c in per.wavelength_air_A:
        sig = _sigma_A(c, inst, resolving)
        try:
            win = load_window_ex(inst, c, SB_OUT + 0.18, holding=hold)
        except Exception as e:
            rows.append(dict(band=band, wavelength_air_A=c, status=str(e)[:80]))
            continue
        w = np.asarray(win.wave, float); f = np.asarray(win.flux, float)
        k = np.isfinite(w) & np.isfinite(f); w, f = w[k], f[k]
        d = np.abs(w - c); sb = (d > HALF) & (d <= SB_OUT)
        if sb.sum() < 5:
            rows.append(dict(band=band, wavelength_air_A=c, status="side-band too thin"))
            continue
        obs = float(np.percentile(f[sb], 95))
        pr = {lab: float(np.percentile(predict(w, W0, TAU0, sig * m)[sb], 95))
              for lab, m in (("lo", 0.6), ("phys", 1.0), ("hi", 2.0))}

        # 🔴 DID THE TELLURIC CORRECTION ACT IN THIS WINDOW AT ALL? A byte comparison
        # against the uncorrected sibling. "Telluric residual after molecfit" presumes a
        # correction happened; where the corrected product IS the raw one, that
        # explanation is unavailable and must not be offered.
        touched = None
        if raw:
            try:
                a = load_window_ex(inst, c, SB_OUT + 0.18, holding=raw, allow_uncorrected=True)
                fa = np.asarray(a.flux, float)
                n = min(len(fa), len(f))
                touched = bool(n and np.nanmax(np.abs(fa[:n] - f[:n])) > 1e-9)
            except Exception:
                touched = None
        rows.append(dict(
            band=band, wavelength_air_A=round(c, 4), sigma_A=round(sig, 5),
            resolving_power=round(resolving.get(inst, float("nan")), 0),
            observed_p95=round(obs, 5),
            catalogued_p95=round(pr["phys"], 5),
            catalogued_p95_sigma_lo=round(pr["lo"], 5),
            catalogued_p95_sigma_hi=round(pr["hi"], 5),
            uncatalogued_residual=round(obs - pr["phys"], 5),
            residual_sigma_hi=round(obs - pr["hi"], 5),
            telluric_correction_applied=touched,
            status="ok"))
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default="data/results/rya1190")
    a = ap.parse_args(argv)
    out_dir = ROOT / a.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    W0, TAU0 = _catalogue()
    resolving = _resolving_powers()
    rows = []
    for band in ("near-UV", "VIS", "red-optical", "NIR-H"):
        got = measure(band, W0, TAU0, resolving)
        rows.extend(got)
        print(f"[{band:<12}] {sum(1 for g in got if g.get('status')=='ok'):>2} lines")
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "rya1190_per_line.csv", index=False)

    ok = df[df.status == "ok"]
    floor = float(ok[ok.band == CONTROL].uncatalogued_residual.mean())

    per_band = []
    for band in ("near-UV", "VIS", "red-optical", "NIR-H"):
        s = ok[ok.band == band]
        if s.empty:
            per_band.append({"band": band, "n": 0}); continue
        tc = s.telluric_correction_applied
        per_band.append({
            "band": band, "n": int(len(s)),
            "observed_p95": round(float(s.observed_p95.mean()), 5),
            "catalogued_p95": round(float(s.catalogued_p95.mean()), 5),
            "uncatalogued_residual": round(float(s.uncatalogued_residual.mean()), 5),
            "excess_over_control": round(float(s.uncatalogued_residual.mean()) - floor, 5),
            "residual_at_sigma_x2": round(float(s.residual_sigma_hi.mean()), 5),
            # ⚠️ counted on the NON-NULL subset only. `tc` is object dtype carrying None
            # where no uncorrected sibling is wired (CRIRES has none), and a naive
            # `(~tc.fillna(True)).sum()` on that returned NEGATIVE counts -- -10 windows
            # "uncorrected" out of 10.
            "n_windows_telluric_corrected": (int(tc.dropna().astype(bool).sum())
                                             if tc.notna().any() else None),
            "n_windows_uncorrected": (int((~tc.dropna().astype(bool)).sum())
                                      if tc.notna().any() else None),
            # 🔴 THE INTERNAL CONTROL FOR "IS IT TELLURIC RESIDUAL". If the depression is
            # residual telluric, windows molecfit ACTED on should differ from windows it
            # left byte-identical to the raw atlas. Same band, same night, same reduction
            # -- the only difference is whether a correction happened.
            "residual_where_corrected": (
                round(float(s[tc.fillna(False).astype(bool)].uncatalogued_residual.mean()), 5)
                if tc.notna().any() and tc.fillna(False).astype(bool).any() else None),
            "residual_where_untouched": (
                round(float(s[~tc.fillna(True).astype(bool)].uncatalogued_residual.mean()), 5)
                if tc.notna().any() and (~tc.fillna(True).astype(bool)).any() else None),
        })

    # ── the ONE remaining completeness lever, bounded ──────────────────────────────
    #
    # 🔴 BOTH "COMPLETE THE LIST" ROUTES REDUCE TO THE SAME POPULATION. The VALD extract
    # IS our production list, so VALD adds nothing AT ITS THRESHOLD; and the vendored
    # molecular surplus (CH MoLLIST, 4,975 lines in band that our list lacks) is material
    # VALD itself evaluated -- its own reference footer credits "Masseron 2014 obs: CH" --
    # and dropped below that same threshold. So the only opacity left to catalogue is the
    # SUB-THRESHOLD WEAK-LINE TAIL, and it can be bounded without fetching anything.
    #
    # VALD processed 2,386,002 lines over 2000-3780 A and selected 161,526 at depth
    # >= 0.001. The discarded ~2.22 M each have tau0 < 0.001. Spread over 1780 A that is
    # ~1250 lines/A, and a Gaussian of width sigma contributes sqrt(2*pi)*sigma of
    # equivalent width each:
    #
    #     tau_haze  <=  (lines/A) * tau0_max * sqrt(2*pi) * sigma
    #
    # ⚠️ UPPER BOUND, and deliberately generous: it puts EVERY discarded line at exactly
    # the threshold when almost all sit far below it, and assumes the 2000-3780 discards
    # are spread uniformly into 3000-3780. If the bound does not reach the measured
    # deficit, the deficit is not catalogueable from this source -- which is the only
    # direction this estimate is used in.
    hdr = (ROOT / VALD["near-UV"][0]).read_text(errors="replace").split("\n", 1)[0]
    parts = [x.strip() for x in hdr.split(",")]
    lo_hdr, hi_hdr = float(parts[0]), float(parts[1])
    n_sel, n_proc = int(parts[2]), int(parts[3])
    n_disc = n_proc - n_sel
    dens = n_disc / (hi_hdr - lo_hdr)
    sig_uv = 3400.0 * DOPPLER_KMS / C_KMS
    tau_haze = dens * 0.001 * np.sqrt(2 * np.pi) * sig_uv
    haze = {
        "vald_header": {"window_A": [lo_hdr, hi_hdr], "selected": n_sel, "processed": n_proc,
                        "discarded_below_threshold": n_disc},
        "discarded_line_density_per_A": round(dens, 1),
        "sigma_A_at_3400": round(sig_uv, 5),
        "max_tau_from_the_haze": round(float(tau_haze), 5),
        "max_depression_it_could_explain": round(float(1.0 - np.exp(-tau_haze)), 5),
        "basis": ("UPPER BOUND: every discarded line placed at exactly the 0.001 depth "
                  "threshold, and the whole 2000-3780 discard attributed to 3000-3780."),
    }

    vald = {b: vald_completeness(b) for b in VALD}
    mol = molecular_coverage()
    mol_min = min((v["min_A"] for v in mol.values()), default=None)

    # ── verdicts, computed ─────────────────────────────────────────────────────────
    ROBUST = 0.5     # a deficit that loses more than half its size across the sigma sweep,
                     # or flips sign, is not separable from the width assumption
    verdicts, payoff = {}, {}
    for row in per_band:
        b = row["band"]
        if not row.get("n") or b == CONTROL:
            verdicts[b] = "CONTROL — side-bands line-free by selection; sets the floor" \
                if b == CONTROL else "NO DATA"
            continue
        exc, hi = row["excess_over_control"], row["residual_at_sigma_x2"]
        robust = (hi < 0) and (abs(hi) >= ROBUST * abs(row["uncatalogued_residual"]))
        if not robust:
            verdicts[b] = (
                f"UNDETERMINED — the deficit is not separable from the line-width "
                f"assumption: {row['uncatalogued_residual']:+.4f} at the physical sigma "
                f"but {hi:+.4f} at sigma x2"
                + (", where it CHANGES SIGN and the catalogue over-predicts absorption"
                   if hi > 0 else "")
                + ". Classifying it either way would be reading the sigma, not the data.")
            continue
        corr, unt = row.get("residual_where_corrected"), row.get("residual_where_untouched")
        if corr is not None and unt is not None and unt < corr:
            verdicts[b] = (
                f"OPACITY — and specifically UNCORRECTED telluric, not residual-after-"
                f"molecfit. The {row['n_windows_uncorrected']} windows the correction "
                f"never touched are MORE depressed ({unt:+.4f}) than the "
                f"{row['n_windows_telluric_corrected']} it acted on ({corr:+.4f}). A "
                f"per-band continuum redo here would divide by real absorption.")
        elif abs(exc) > 0.10:
            verdicts[b] = (
                f"OPACITY-DOMINATED — {exc:+.4f} excess over the control, robust across "
                f"the sigma sweep, in a band whose catalogued opacity is already complete "
                f"to the VALD threshold. The deficit is real and is NOT catalogueable "
                f"from what we hold.")
        else:
            verdicts[b] = (f"SMALL RESIDUAL DEFICIT {exc:+.4f}, robust to sigma; no "
                           f"telluric split available to attribute it.")

    # payoff: the catalogueable FRACTION of the near-UV deficit, times RYA-1189's dA
    uv = next((r for r in per_band if r["band"] == "near-UV"), None)
    if uv and uv.get("n"):
        frac = min(1.0, haze["max_depression_it_could_explain"] / abs(uv["excess_over_control"]))
        payoff = {
            "near_uv_deficit_excess_over_control": uv["excess_over_control"],
            "max_catalogueable_fraction": round(frac, 3),
            "rya1189_dA_if_all_of_it_were_deblended": -0.086,
            "estimated_payoff_dex": round(-0.086 * frac, 4),
            "published_gap_to_anchor_dex": [0.13, 0.19],
            "reading": ("UPPER BOUND on both factors. Even crediting the weak-line haze "
                        "with its maximum and crediting a full deblend with RYA-1189's "
                        "whole -0.086 dex, a completeness build closes ~"
                        f"{abs(round(-0.086 * frac, 3)):.3f} dex of a 0.13-0.19 dex gap. "
                        "⚠️ The -0.086 is itself a lower bound in magnitude (the near-UV "
                        "lines are saturated), so this estimate is not tight in either "
                        "direction and is offered as an ORDER, not a number to plan on."),
        }

    doc = {
        "ticket": "RYA-1190",
        "verdict_per_band": verdicts,
        "part_A_payoff_estimate": payoff,
        "kind": "DIAGNOSTIC — nothing added to any list, no published value changed",
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "method": {
            "what": "catalogued-opacity accounting, NOT a synthesis differential",
            "why": ("Turbospectrum's bsyn_lu is absent from this Mac and MOOG is a Linux "
                    "ELF; `spectrum` is arm64 but is not the production engine (RYA-289). "
                    "A true synthesis differential is owed as a Sirius run."),
            "sigma_basis": (f"quadrature of a {DOPPLER_KMS} km/s Doppler width and the "
                            f"INSTRUMENTAL width from instrument_catalog.csv's declared "
                            f"resolving power (geometric mean of min/max); swept x0.6 and x2.0"),
            "resolving_powers_used": {k: round(v) for k, v in resolving.items()
                                      if k in {h[0] for h in HOLDING.values()}},
            "control_floor_subtracted": round(floor, 5),
        },
        "part_A_vald_completeness": vald,
        "part_A_weak_line_haze_bound": haze,
        "part_A_molecular_coverage": {
            "by_species": mol,
            "lowest_wavelength_A": mol_min,
            "covers_near_uv_3000_3780": bool(mol_min is not None and mol_min < 3780.0),
        },
        "per_band": per_band,
        "control_band": CONTROL,
        "optical_anchor_for_the_closing_check_only": OPTICAL_ANCHOR,
    }
    (out_dir / "rya1190_frontier_opacity.json").write_text(json.dumps(doc, indent=2) + "\n")
    pd.DataFrame(per_band).to_csv(out_dir / "rya1190_per_band.csv", index=False)

    print("\n=== per-band ===")
    print(pd.DataFrame(per_band).to_string(index=False))
    print(f"\ncontrol floor subtracted: {floor:+.5f}  (VIS, side-bands line-free by selection)")
    print("\n=== Part A: is the production list already the VALD extract? ===")
    for b, v in vald.items():
        if v.get("available"):
            print(f"  {b:<12} VALD {v['vald_rows']:>6} | ours {v['production_rows']:>6} | "
                  f"VALD-only {v['in_vald_only']:>5} | ours-only {v['in_production_only']:>5} "
                  f"| threshold {v['extraction_depth_threshold']}")
        else:
            print(f"  {b:<12} no vendored VALD extract at {v['path']}")
    print(f"\n=== molecular lists span from {mol_min} A "
          f"(covers near-UV: {doc['part_A_molecular_coverage']['covers_near_uv_3000_3780']}) ===")
    for k, v in mol.items():
        print(f"  {k:<18} n={v['n_lines']:>7}  {v['min_A']}..{v['max_A']} A")
    print("\n=== verdicts ===")
    for b, v in verdicts.items():
        print(f"  {b:<12} {v}")
    if payoff:
        print(f"\n=== Part A payoff (upper bound) ===")
        print(f"  catalogueable fraction of the deficit : {payoff['max_catalogueable_fraction']}")
        print(f"  estimated dA if built                 : {payoff['estimated_payoff_dex']} dex")
        print(f"  against a published gap of            : {payoff['published_gap_to_anchor_dex']} dex")
    print(f"\nwrote {a.out_dir}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
