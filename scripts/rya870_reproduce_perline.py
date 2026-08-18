#!/usr/bin/env python3
"""
RYA-870 THE GUARD — can a published row reproduce its own A(X) from its own constants?

    # Sirius, venv312 (iSpec present)
    python3 scripts/rya870_reproduce_perline.py --star solar --element Fe --sample 8

This is what makes the product replication-grade rather than merely wide. RYA-489 §6.1
states the test the product must pass: *someone with the same stack can reproduce any
per-line A(X) from that row plus the file header, with no access to our pipeline.*

🔴 SO THIS READS THE PRODUCT, NOT THE PIPELINE. Every physical input to the re-inversion —
log gf, excitation potential, the three damping constants, the equivalent width, Teff,
log g, [Fe/H], microturbulence — is taken from the emitted CSV and its header. If a row
cannot reproduce its number from what it publishes, the row is not replication-grade, and
that is a FAILURE of the product, not of this script.

⚠️ WHAT IS BORROWED, AND WHY IT IS NOT CHEATING. iSpec's linelist array has a fixed dtype,
so the ARRAY STRUCTURE is loaded from the project's linelist and the row's own values are
written into it. Structure is plumbing; values are the claim. A reproducer that also took
the VALUES from the pipeline would be testing the pipeline against itself and would pass
even if the published columns were wrong — which is the failure this guard exists to catch.

⚠️ EW ROUTE ONLY, STATED PLAINLY. A `synthesis_fit` row is reproduced by re-fitting an
observed spectrum, which needs the spectrum as well as the row, and the row alone cannot
carry it. Those rows are reported as NOT-COVERED rather than quietly counted as passing —
an untested row must never look like a tested one.

⚠️ THE FLOOR IS MEASURED PER LINE, NOT ASSUMED. Two bounds apply and the larger wins:

  * the RYA-771 QUANTISER floor, 0.01 dex — iSpec rounds trial abundances in bsyn, so a
    tighter threshold would be measuring the quantiser rather than the product; and
  * the INVERSION's own convergence. `_bisect_synth_abundance` stops when the synthetic
    EW is within `tol` (relative) of its target — it converges on EQUIVALENT WIDTH, not on
    abundance. Near saturation dEW/dA flattens, so the same 1% in EW buys a much larger
    interval in A, and demanding 0.01 dex there asks the reproduction to be tighter than
    the number it is reproducing ever was.

🔴 That second bound is MEASURED for each row — synthesise either side of the reproduced
abundance, read the local dEW/dA off the result, and convert the EW tolerance into dex:

        dA_tol = tol * max(ew_target, 1) / |dEW/dA|

which is why a weak line is held to 0.01 dex and a saturated one is held to what its own
slope permits. The alternative — one global threshold loose enough for the worst line —
would hide a real defect on every other line, and picking it by eye would be tuning.
`tol` and the bracket floor are read from the deriver's own signature so the number is
declared once (RYA-845).
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

QUANTISER_FLOOR_DEX = 0.01   # RYA-771
OUT = ROOT / "data" / "results" / "rya870"


def _read_product(path: Path):
    header = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if ": " in line:
                k, v = line[1:].split(": ", 1)
                header[k.strip()] = v.strip()
    rows = pd.read_csv(path, comment="#")
    return header, rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--star", default="solar")
    ap.add_argument("--element", default="Fe")
    ap.add_argument("--product", type=Path, default=None)
    ap.add_argument("--sample", type=int, default=8)
    ap.add_argument("--seed", type=int, default=847)
    a = ap.parse_args()

    product = a.product or (ROOT / "data" / "products" / a.star /
                            f"{a.element}_perline.csv")
    header, rows = _read_product(product)

    # ── everything physical comes from the PRODUCT ────────────────────────────────
    teff = float(header["Teff"])
    logg = float(header["log_g"])
    feh = float(header["Fe_H"])
    xi = float(header["xi_microturb"])

    # 🔴 COMPARE LIKE WITH LIKE. This guard re-inverts an equivalent width, which is an
    # LTE operation. Comparing that against an ENGINE-A number — the same inversion PLUS
    # the MPIA departure — charges the row for a correction the reproduction never
    # applied, and it is exactly how RYA-880 was found: a stubborn 0.011 dex on 6094.372
    # that turned out to BE the departure. Two ways to be honest about it, in order of
    # preference:
    #   1. the product records the correction (RYA-880) -> ADD it back, and the comparison
    #      is exact on every engine;
    #   2. the product predates that column -> restrict to rows that are LTE by
    #      construction, and SAY the coverage is narrower rather than quietly widening a
    #      threshold to swallow the difference.
    ew_rows = rows[(rows["status"] == "in_aggregate")
                   & (rows["method"] == "ew_integration")
                   & rows["ew_mA"].notna() & rows["A_X_line"].notna()
                   & rows["log_gf"].notna() & rows["excitation_potential_eV"].notna()]
    has_delta = ("nlte_delta_dex" in rows.columns
                 and rows["nlte_delta_dex"].notna().any())
    if not has_delta:
        n_before = len(ew_rows)
        ew_rows = ew_rows[ew_rows["engine"].astype(str).str.startswith("1D-LTE")]
        print(f"  [coverage] the product carries no recorded NLTE correction, so this run "
              f"is restricted to the {len(ew_rows)} rows that are LTE by construction "
              f"(of {n_before} EW-route rows). RYA-880 removes this restriction.")
    not_covered = int(((rows["status"] == "in_aggregate")
                       & (rows["method"] == "synthesis_fit")).sum())

    if ew_rows.empty:
        print(json.dumps({"status": "NO-EW-ROWS", "n_not_covered": not_covered}, indent=2))
        return 2

    sample = ew_rows.sample(n=min(a.sample, len(ew_rows)), random_state=a.seed)

    from pipeline import abundances_derive as ad          # noqa: E402  (iSpec import)
    import ispec                                          # noqa: E402

    linelist, isotopes, chem = ad._load_synth_resources()
    solar_abund = ispec.read_solar_abundances(ad._ISPEC_SOLAR_ABUND_FILE)
    atmosphere = ad._load_atmosphere(teff, logg, feh, xi)

    # the deriver's own convergence parameters — read, never re-typed (RYA-845)
    import inspect
    _sig = inspect.signature(ad._bisect_synth_abundance).parameters
    BISECT_TOL = float(_sig["tol"].default)
    BISECT_A_LO = float(_sig["a_lo"].default)

    results = []
    for _, r in sample.iterrows():
        wl_A = float(r["wavelength_air_A"])
        # 🔴 CALL THE DERIVER'S CONFIGURATION, DO NOT RECONSTRUCT ONE (RYA-832).
        # Two earlier cuts of this guard invented their own: a single isolated line (three
        # bisections ran out of bracket, one landed +0.49 dex HIGH — the target charged
        # with its blends' EW) and then a +/-1.5 A linelist subset (all five converged LOW
        # by 0.07-1.6 dex — the opposite error). Neither measured the product; both
        # measured my reconstruction. derive_band_products inverts with the FULL linelist
        # over linspace(c-0.6, c+0.6, 300), so this does exactly that, and the ONLY
        # difference from the production call is that the target line's constants are
        # taken from the PUBLISHED ROW instead of from the pipeline's linelist. That
        # difference is the entire claim under test.
        # 🔴 WAVELENGTH *AND* EP. This guard documented the trap in the generator and then
        # walked into it here: a 0.05 A window alone selected a high-excitation neighbour
        # (5852.218 came back at EP 8.548 / log gf -5.965 against the published 4.549 /
        # -1.23), so the override was applied to the wrong transition and the row was
        # blamed for it. RYA-780/852.
        near = np.abs(linelist["wave_A"] - wl_A) < 0.05
        ep_pub = float(r["excitation_potential_eV"])
        cand = np.where(near & (np.abs(linelist["lower_state_eV"] - ep_pub) < 0.05))[0]
        if cand.size == 0:
            results.append({"wavelength_air_A": wl_A, "published_ep_eV": ep_pub,
                            "outcome": "LINE-NOT-IN-SYNTH-LINELIST-AT-THIS-EP"})
            continue
        one = linelist.copy()
        idx = int(cand[np.argmin(np.abs(linelist["wave_A"][cand] - wl_A))])
        one["loggf"][idx] = float(r["log_gf"])
        one["lower_state_eV"][idx] = float(r["excitation_potential_eV"])
        for col, key in (("rad", "damping_rad"), ("stark", "damping_stark"),
                         ("waals", "damping_vdW")):
            if col in one.dtype.names and pd.notna(r.get(key)):
                one[col][idx] = float(r[key])

        wave_nm = np.linspace(wl_A - 0.6, wl_A + 0.6, 300) / 10.0
        a_repro, converged, _ = ad._bisect_synth_abundance(
            wave_nm, float(r["ew_mA"]), atmosphere, teff, logg, feh, xi,
            one, isotopes, solar_abund, str(r["element"]), 26)

        published = float(r["A_X_line"])
        # RYA-880: put the recorded departure back, so an NLTE-corrected published value
        # is compared against an LTE re-inversion on equal terms.
        applied = r.get("nlte_delta_dex")
        if applied is not None and not pd.isna(applied) and float(applied) != 0.0:
            published = published - float(applied)
        delta = (float("nan") if not np.isfinite(a_repro)
                 else abs(a_repro - published))

        # ── the per-row bound, MEASURED ──────────────────────────────────────────
        kw = dict(atmosphere=atmosphere, teff=teff, logg=logg, feh=feh, vturb=xi,
                  linelist=one, isotopes=isotopes, solar_abund=solar_abund,
                  element=str(r["element"]), atom_code=26)
        dA_tol = float("nan")
        if np.isfinite(a_repro):
            try:
                h = 0.02
                ew_floor = ad._synth_ew_at_abund(wave_nm, trial_A=BISECT_A_LO, **kw)
                ew_target = ew_floor + float(r["ew_mA"])
                ew_hi = ad._synth_ew_at_abund(wave_nm, trial_A=a_repro + h, **kw)
                ew_lo_ = ad._synth_ew_at_abund(wave_nm, trial_A=a_repro - h, **kw)
                slope = abs(ew_hi - ew_lo_) / (2 * h)      # mA per dex
                if slope > 0:
                    dA_tol = BISECT_TOL * max(ew_target, 1.0) / slope
            except Exception:
                dA_tol = float("nan")
        bound = (QUANTISER_FLOOR_DEX if not np.isfinite(dA_tol)
                 else max(QUANTISER_FLOOR_DEX, dA_tol))
        results.append({
            "wavelength_air_A": float(r["wavelength_air_A"]),
            "engine": str(r["engine"]),
            "published_A": published,
            "reproduced_A": (None if not np.isfinite(a_repro) else float(a_repro)),
            "delta_dex": (None if not np.isfinite(delta) else float(delta)),
            "converged": bool(converged),
            "reduced_ew": (None if pd.isna(r.get("reduced_ew"))
                           else float(r["reduced_ew"])),
            "ew_tolerance_bound_dex": (None if not np.isfinite(dA_tol) else float(dA_tol)),
            "bound_applied_dex": float(bound),
            "bound_set_by": ("inversion EW tolerance" if np.isfinite(dA_tol)
                             and dA_tol > QUANTISER_FLOOR_DEX else "RYA-771 quantiser"),
            "outcome": ("PASS" if np.isfinite(delta) and delta <= bound else "FAIL"),
        })
        print(f"  {r['wavelength_air_A']:9.3f}  published {published:.4f}  "
              f"reproduced {a_repro:.4f}  delta {delta:.4f}  "
              f"bound {bound:.4f} ({results[-1]['bound_set_by']})  "
              f"{results[-1]['outcome']}")

    tested = [x for x in results if x.get("outcome") in ("PASS", "FAIL")]
    failed = [x for x in tested if x["outcome"] == "FAIL"]
    report = {
        "ticket": "RYA-870",
        "product": str(product.relative_to(ROOT)),
        "product_commit_sha": header.get("commit_sha"),
        "quantiser_floor_dex": QUANTISER_FLOOR_DEX,
        "bisect_rel_tol": BISECT_TOL,
        "bound_rule": ("max(RYA-771 quantiser floor, tol * max(ew_target,1) / |dEW/dA|) "
                       "— the reproduction cannot be held tighter than the convergence "
                       "criterion that produced the number"),
        "n_rows_in_product": int(len(rows)),
        "n_ew_route_in_aggregate": int(len(ew_rows)),
        "n_sampled": int(len(sample)),
        "n_tested": len(tested),
        "n_passed": len(tested) - len(failed),
        "n_failed": len(failed),
        "n_synthesis_rows_NOT_COVERED": not_covered,
        "coverage_note": (
            "synthesis_fit rows are not reproducible from the row alone — they need the "
            "observed spectrum. Reported as not-covered, never counted as passing."),
        "results": results,
        "verdict": "PASS" if tested and not failed else ("FAIL" if failed else "NO-COVERAGE"),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "rya870_reproducibility.json").write_text(json.dumps(report, indent=2) + "\n")
    bounds = [x["bound_applied_dex"] for x in tested if x.get("bound_applied_dex")]
    print(f"\n{report['verdict']}: {report['n_passed']}/{report['n_tested']} reproduced "
          f"within their own measured bound "
          f"({min(bounds):.4f}-{max(bounds):.4f} dex); "
          f"{not_covered} synthesis rows not covered")
    print(f"[out] {OUT}/rya870_reproducibility.json")
    return 0 if report["verdict"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
