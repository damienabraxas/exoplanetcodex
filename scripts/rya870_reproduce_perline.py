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

⚠️ THE FLOOR IS 0.01 dex, and it is not a tolerance anyone chose: iSpec quantises trial
abundances to 0.01 dex in bsyn (RYA-771), so a smaller threshold would be measuring the
quantiser rather than the product.
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

    ew_rows = rows[(rows["status"] == "in_aggregate")
                   & (rows["method"] == "ew_integration")
                   & rows["ew_mA"].notna() & rows["A_X_line"].notna()
                   & rows["log_gf"].notna() & rows["excitation_potential_eV"].notna()]
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
        delta = (float("nan") if not np.isfinite(a_repro)
                 else abs(a_repro - published))
        results.append({
            "wavelength_air_A": float(r["wavelength_air_A"]),
            "engine": str(r["engine"]),
            "published_A": published,
            "reproduced_A": (None if not np.isfinite(a_repro) else float(a_repro)),
            "delta_dex": (None if not np.isfinite(delta) else float(delta)),
            "converged": bool(converged),
            "outcome": ("PASS" if np.isfinite(delta) and delta <= QUANTISER_FLOOR_DEX
                        else "FAIL"),
        })
        print(f"  {r['wavelength_air_A']:9.3f}  published {published:.4f}  "
              f"reproduced {a_repro:.4f}  delta {delta:.4f}  "
              f"{results[-1]['outcome']}")

    tested = [x for x in results if x.get("outcome") in ("PASS", "FAIL")]
    failed = [x for x in tested if x["outcome"] == "FAIL"]
    report = {
        "ticket": "RYA-870",
        "product": str(product.relative_to(ROOT)),
        "product_commit_sha": header.get("commit_sha"),
        "quantiser_floor_dex": QUANTISER_FLOOR_DEX,
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
    print(f"\n{report['verdict']}: {report['n_passed']}/{report['n_tested']} reproduced "
          f"within {QUANTISER_FLOOR_DEX} dex; {not_covered} synthesis rows not covered")
    print(f"[out] {OUT}/rya870_reproducibility.json")
    return 0 if report["verdict"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
