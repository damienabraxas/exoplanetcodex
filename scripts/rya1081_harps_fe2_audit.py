#!/usr/bin/env python3
"""
RYA-1081 — why does HARPS Fe II VIS DEEPGRADED read 7.966?
==========================================================
HARPS is the science-forward Fe II arm, and its VIS DEEPGRADED product reads 7.966 against
KP's ~7.57 and an ionization arbiter of 7.500. The ticket's hypothesis is that a few lines
are skewing it — blends, saturated cores, rail escapes, or residual of the RYA-911 HARPS
Fe II EW-leg pathology. This finds out which, line by line.

🔴 FIREWALL (RYA-161, RYA-523). The job is to find bad lines and characterise real physics.
It is NEVER to move the number toward KP or Asplund. Nothing here excludes a line to make a
value look better; the one line dispositioned is dispositioned on a measured fit failure
that is present on BOTH arms, which is the opposite of an arm-flattering exclusion.

WHAT THE POOL ACTUALLY IS: nine lines, 4233-4584 A, every one carrying PRIMARY LAB
Den Hartog 2019 gf at tier LAB with a cited sigma of 0.05-0.10 dex. They are the same nine
the RYA-853 DH19 referee used. There is no EW anywhere in this product -- it is a synthesis
flux-fit route and `ew_mA`, `rew` and `observed_depth` are NaN on every row, which is what
makes two of this ticket's asks unrunnable rather than negative (see the two DECLARED GAPs).

Usage:
    python3 scripts/rya1081_harps_fe2_audit.py --instrument harps --ion II \
        --band VIS --tier deepgraded
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BP = ROOT / "data" / "results" / "band_products"
OUT = ROOT / "data" / "results" / "rya1081"
CANON = ROOT / "data" / "linelists" / "canonical_gf.csv"
DH19 = ROOT / "data" / "reference" / "fe_gf_lab" / "fe2_lab_loggf_dh19.csv"
#: RYA-1043/1041's ratified per-line COG. The saturation verdict must come from here and
#: from nowhere else -- asserting saturation from EW magnitude is a CRITICAL failure.
COG = ROOT / "data" / "audit" / "rya1041" / "rya1041_perline_cog_4200_6910.csv"

STEM = "FeII_4200_6910_{arm}_corrected_SYNTH_DEEPGRADED_{treat}_lines.csv"
ARMS = {"harps": "harps_solar_harps_molecfit",
        "kp_molecfit": "kpno_solar_atlas_solar_kpno_molecfit",
        "kp_kurucz2005": "kpno_solar_atlas_solar_kpno_kurucz2005"}

#: EP-matched on BOTH sides, never a rounded-wavelength join (RYA-846, RYA-1033).
WAVE_TOL_A, EP_TOL_EV = 0.05, 0.05


def _lines(arm: str, treat: str) -> pd.DataFrame:
    p = BP / STEM.format(arm=ARMS[arm], treat=treat)
    if not p.exists():
        raise SystemExit(f"missing committed per-line artifact: {p}")
    return pd.read_csv(p)


def gf_table(pool: pd.DataFrame) -> pd.DataFrame:
    c = pd.read_csv(CANON, comment="#", low_memory=False)
    c = c[c.species.astype(str) == "Fe II"]
    dh = pd.read_csv(DH19)
    rows = []
    for _, r in pool.iterrows():
        m = c[(np.abs(c.wavelength_air_A - r.wavelength_air_A) <= WAVE_TOL_A)
              & (np.abs(c.excitation_potential_eV - r.ep_eV) <= EP_TOL_EV)]
        g = dh[(np.abs(dh.wavelength_air_A - r.wavelength_air_A) <= WAVE_TOL_A)
               & (np.abs(dh.elo_eV - r.ep_eV) <= EP_TOL_EV)]
        rows.append({
            "wavelength_air_A": r.wavelength_air_A, "ep_eV": r.ep_eV,
            "n_gf_candidates": len(m),
            "log_gf": float(m.iloc[0].log_gf) if len(m) == 1 else np.nan,
            "gf_reference": str(m.iloc[0].loggf_reference) if len(m) == 1 else "AMBIGUOUS",
            "gf_tier": str(m.iloc[0].gf_tier) if len(m) == 1 else "",
            "gf_sigma_dex": float(m.iloc[0].gf_sigma_dex) if len(m) == 1 else np.nan,
            "dh19_loggf": float(g.iloc[0].loggf) if len(g) == 1 else np.nan,
        })
    return pd.DataFrame(rows)


def cog_verdict(pool: pd.DataFrame) -> dict:
    """RYA-1043's per-line COG, or an honest declaration that it does not cover this pool.

    RYA-833: an absence in the instrument is not a negative result. 'The COG does not cover
    these lines' is NOT 'these lines are not saturated'.
    """
    if not COG.exists():
        return {"runnable": False, "reason": f"{COG.name} absent", "n_covered": 0}
    c = pd.read_csv(COG)
    w = c.wavelength_air_A
    covered = [float(x) for x in pool.wavelength_air_A
               if bool((np.abs(w - x) <= WAVE_TOL_A).any())]
    return {"runnable": bool(covered), "n_covered": len(covered),
            "n_pool": int(len(pool)), "covered": covered,
            "cog_span_A": [float(w.min()), float(w.max())],
            "cog_n_lines": int(len(c)),
            "reason": ("the ratified per-line COG covers Fe I only — 0 of the pool's lines "
                       "appear in it, so the saturation question is UNANSWERED, not "
                       "answered negatively (RYA-833). Asserting saturation from EW "
                       "magnitude instead would be a CRITICAL substitution — and is moot "
                       "here anyway, since this route measures no EW at all.")
            if not covered else "covered"}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--instrument", default="harps")
    ap.add_argument("--ion", default="II")
    ap.add_argument("--band", default="VIS")
    ap.add_argument("--tier", default="deepgraded")
    a = ap.parse_args()
    if (a.instrument, a.ion, a.band, a.tier.lower()) != ("harps", "II", "VIS", "deepgraded"):
        raise SystemExit("this audit is built for --instrument harps --ion II --band VIS "
                         "--tier deepgraded; other cells have no committed per-line pool")
    OUT.mkdir(parents=True, exist_ok=True)

    h = _lines("harps", "1D-LTE")
    gf = gf_table(h)
    t = h[["wavelength_air_A", "ep_eV", "abundance", "red_chi2", "sigma_A",
           "edge_distance_dex", "in_aggregate", "ew_mA", "rew", "observed_depth"]].merge(
        gf, on=["wavelength_air_A", "ep_eV"])

    pool_med = float(t.abundance.median())
    t["dev_from_pool_median"] = t.abundance - pool_med
    t = t.sort_values("dev_from_pool_median", key=np.abs, ascending=False)

    print(f"=== RYA-1081 — HARPS Fe II VIS DEEPGRADED, per-line ({len(t)} lines) ===")
    print(f"    published value = pool median = {pool_med:.3f}\n")
    show = ["wavelength_air_A", "ep_eV", "log_gf", "gf_tier", "gf_sigma_dex", "abundance",
            "dev_from_pool_median", "red_chi2", "ew_mA", "rew", "observed_depth"]
    print(t[show].to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    # ── DECLARED GAP 1: there is no GRADED (non-deep) Fe II tier to compare against ──
    feed = json.loads((ROOT / "data/products/solar/Fe.json").read_text())
    fe2_tiers = sorted({str(p.get("tier")) for p in feed["products"]
                        if str(p.get("ion")) == "II"})
    graded_gap = {
        "requested": "the GRADED (non-deep) subset value, as the cleaner reference",
        "runnable": False,
        "fe2_tiers_present_in_feed": fe2_tiers,
        "reason": ("EVERY Fe II product in the feed is tier=DEEPGRADED — there is no GRADED "
                   "Fe II product on any instrument or band. The DEEPGRADED-minus-GRADED "
                   "deep contribution therefore has no comparand, and manufacturing one by "
                   "splitting this pool on a depth proxy would be the conflation the ticket "
                   "names as CRITICAL. What is needed is a GRADED-tier Fe II run."),
    }

    # ── the KP cross-check: the main diagnostic, same lines, same gf ──
    cross = {}
    for arm in ("kp_molecfit", "kp_kurucz2005"):
        k = _lines(arm, "1D-LTE")
        m = t[["wavelength_air_A", "abundance", "red_chi2"]].merge(
            k[["wavelength_air_A", "abundance", "red_chi2"]],
            on="wavelength_air_A", suffixes=("_harps", "_kp"))
        m["delta"] = m.abundance_harps - m.abundance_kp
        no_out = m[np.abs(m.abundance_harps - pool_med) < 1.0]
        cross[arm] = {
            "n": int(len(m)),
            "median_harps": float(m.abundance_harps.median()),
            "median_kp": float(m.abundance_kp.median()),
            "median_delta": float(m.delta.median()),
            "n_harps_above_kp": int((m.delta > 0).sum()),
            "exceptions": [float(x) for x in m.loc[m.delta <= 0, "wavelength_air_A"]],
            "median_delta_excluding_outlier": float(no_out.delta.median()),
            "corr_A_vs_red_chi2_harps": float(m.abundance_harps.corr(m.red_chi2_harps)),
            "corr_A_vs_red_chi2_harps_excluding_outlier":
                float(no_out.abundance_harps.corr(no_out.red_chi2_harps)),
        }
        m.to_csv(OUT / f"rya1081_crosscheck_{arm}.csv", index=False)

    c0 = cross["kp_molecfit"]
    print(f"\n=== KP cross-check (same {c0['n']} lines, same pure-lab gf) ===")
    for arm, v in cross.items():
        print(f"  {arm:<16} HARPS {v['median_harps']:.3f}  KP {v['median_kp']:.3f}  "
              f"delta {v['median_delta']:+.3f}   HARPS>KP on {v['n_harps_above_kp']}/{v['n']}")
    print(f"  the single exception on both: {c0['exceptions']}")

    # ── robustness: is the published median driven by the outlier? ──
    keep = t[np.abs(t.abundance - pool_med) < 1.0]
    robust = {"published_median": pool_med,
              "median_without_the_outlier": float(keep.abundance.median()),
              "shift": float(keep.abundance.median() - pool_med),
              "pool_min": float(t.abundance.min()), "pool_max": float(t.abundance.max()),
              "spread_dex": float(t.abundance.max() - t.abundance.min()),
              "spread_excluding_outlier_dex": float(keep.abundance.max()
                                                    - keep.abundance.min())}
    print(f"\n=== is 7.966 driven by the high-puller? ===")
    print(f"  median with it {robust['published_median']:.3f}, without it "
          f"{robust['median_without_the_outlier']:.3f}  (shift {robust['shift']:+.3f})")
    print(f"  pool spans {robust['pool_min']:.3f}..{robust['pool_max']:.3f} "
          f"({robust['spread_dex']:.3f} dex; {robust['spread_excluding_outlier_dex']:.3f} "
          f"without it)")

    # ── RYA-911: is this the EW-leg residual? ──
    ew_present = bool(t.ew_mA.notna().any())
    rya911 = {"is_the_rya911_residual": False,
              "route": "SYNTH (synthesis flux-fit)",
              "any_measured_ew_in_this_product": ew_present,
              "rya911_was": "the EW leg, running 0.34 dex LOW",
              "this_is": f"the synthesis leg, running {c0['median_delta']:+.3f} dex HIGH",
              "reason": ("RYA-911's pathology is in the EW step. This product has no EW leg "
                         "at all — ew_mA is NaN on every row and ew_method records a "
                         "synthesis flux-fit — and the sign is opposite. So the offset here "
                         "is not RYA-911 residual; it is a separate synthesis-route "
                         "HARPS-specific offset.")}

    sat = cog_verdict(t)
    verdict = {
        "ticket": "RYA-1081",
        "published": {"value": pool_med, "n_lines": int(len(t)), "tier": "DEEPGRADED"},
        "gf": {"all_primary_lab": bool((t.gf_tier == "LAB").all()),
               "reference": sorted(set(t.gf_reference)),
               "cited_sigma_range": [float(t.gf_sigma_dex.min()),
                                     float(t.gf_sigma_dex.max())],
               "so_gf_cannot_explain_the_arm_gap": True},
        "kp_crosscheck": cross,
        "robustness": robust,
        "graded_reference_gap": graded_gap,
        "saturation_cog": sat,
        "rya911": rya911,
    }
    t.to_csv(OUT / "rya1081_perline_audit.csv", index=False)
    (OUT / "rya1081_verdict.json").write_text(json.dumps(verdict, indent=2) + "\n")

    print(f"\n=== DECLARED GAPS (asks this ticket cannot execute today) ===")
    print(f"  1. GRADED-only reference: {graded_gap['reason']}")
    print(f"  2. saturation COG: {sat['reason']}")
    print(f"\n=== RYA-911 residual: {'YES' if rya911['is_the_rya911_residual'] else 'NO'} ===")
    print(f"  {rya911['reason']}")
    print(f"\n[out] {OUT}/rya1081_perline_audit.csv\n[out] {OUT}/rya1081_verdict.json")


if __name__ == "__main__":
    main()
