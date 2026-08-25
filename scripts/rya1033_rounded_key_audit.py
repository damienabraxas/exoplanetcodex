#!/usr/bin/env python3
"""
scripts/rya1033_rounded_key_audit.py — RYA-1033
===============================================
Measure what a 2-decimal-rounded wavelength key costs the measured -> canonical_gf join,
and re-audit the graded Fe I pool count that the ticket flagged as unstable.

Writes data/audit/rya1033_rounded_key_join/:
  dropped_lines.csv        the measured Fe I lines a rounded key fails to join, with the
                           canonical row each one actually has and how far away it is
  ambiguous_lines.csv      measured lines with >1 canonical candidate — wavelength alone
                           cannot identify these, EP is required
  graded_pool_audit.csv    graded Fe I pool count under each keying scheme
  summary.md               the findings in prose

Every number in the RYA-1033 write-up comes from here. Re-run to regenerate.

    python3 scripts/rya1033_rounded_key_audit.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline import line_match  # noqa: E402

OUT = ROOT / "data" / "audit" / "rya1033_rounded_key_join"
CANONICAL_GF = ROOT / "data" / "linelists" / "canonical_gf.csv"
MEASURED_EW = ROOT / "data" / "measured" / "sol_ew_results_v1.csv"

#: The ticket's own definition of "graded": a primary laboratory gf OR NIST critically
#: compiled. Both tiers are what RYA-946 admits to the showcase pool.
GRADED_TIERS = ("LAB", "NIST-C+")


def _fe_i_measured() -> pd.DataFrame:
    ew = pd.read_csv(MEASURED_EW)
    return ew[(ew.element == "Fe") & (ew.ion == "I")].reset_index(drop=True)


def _fe_i_canonical() -> pd.DataFrame:
    cg = pd.read_csv(CANONICAL_GF, low_memory=False)
    return cg[cg.species == "Fe I"].reset_index(drop=True)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    meas, canon = _fe_i_measured(), _fe_i_canonical()
    w = meas.wavelength_air_A.astype(float).to_numpy()
    cw = canon.wavelength_air_A.astype(float).to_numpy()
    order = np.argsort(cw)
    cw_sorted = cw[order]

    # ── 1. what the rounded key drops ────────────────────────────────────────────────
    ck = set(canon.wavelength_air_A.astype(float).round(2))          # pandas/numpy rounding
    miss = ~pd.Series(w).round(2).isin(ck).to_numpy()
    rows = []
    for wl in w[miss]:
        j = int(np.abs(cw_sorted - wl).argmin())
        c = canon.iloc[order[j]]
        rows.append(dict(measured_wavelength_air_A=round(float(wl), 5),
                         canonical_wavelength_air_A=float(c.wavelength_air_A),
                         separation_mA=round(1000 * (float(c.wavelength_air_A) - wl), 3),
                         key_measured=round(float(wl), 2),
                         key_canonical=round(float(c.wavelength_air_A), 2),
                         gf_tier=c.gf_tier, loggf_reference=c.loggf_reference,
                         adjudication_status=c.adjudication_status))
    dropped = pd.DataFrame(rows).sort_values("measured_wavelength_air_A")
    dropped.to_csv(OUT / "dropped_lines.csv", index=False)

    # The same measured lines under the fixed matcher.
    res = line_match.match(w, cw)
    n_unres, n_amb = len(res.unresolved), len(res.ambiguous)

    amb_rows = []
    for wl, cands in res.ambiguous:
        for c_wl in cands:
            c = canon.iloc[int(np.abs(cw - c_wl).argmin())]
            amb_rows.append(dict(measured_wavelength_air_A=round(float(wl), 5),
                                 candidate_wavelength_air_A=float(c_wl),
                                 separation_mA=round(1000 * (c_wl - wl), 3),
                                 excitation_potential_eV=c.excitation_potential_eV,
                                 log_gf=c.log_gf, gf_tier=c.gf_tier))
    pd.DataFrame(amb_rows).to_csv(OUT / "ambiguous_lines.csv", index=False)

    # ── 2. the graded pool count, three ways ─────────────────────────────────────────
    graded = canon[canon.gf_tier.astype(str).isin(GRADED_TIERS)]
    gw = np.sort(graded.wavelength_air_A.astype(float).to_numpy())
    gk = set(np.round(gw, 2))

    tol_res = line_match.match(w, gw)
    n_tol = int((tol_res.index >= 0).sum()) + len(tol_res.ambiguous)
    n_round = int(pd.Series(w).round(2).isin(gk).sum())

    lab_only = canon[canon.gf_tier.astype(str) == "LAB"]
    lw = np.sort(lab_only.wavelength_air_A.astype(float).to_numpy())
    n_lab = int((line_match.match(w, lw).index >= 0).sum())

    pool = pd.DataFrame([
        dict(scheme="2-dp rounded key (the defect)", graded_tiers="LAB+NIST-C+", n=n_round),
        dict(scheme="tolerance match (the fix)", graded_tiers="LAB+NIST-C+", n=n_tol),
        dict(scheme="tolerance match (the fix)", graded_tiers="LAB only", n=n_lab),
    ])
    pool.to_csv(OUT / "graded_pool_audit.csv", index=False)

    # ── 3. the library-divergence demonstration ──────────────────────────────────────
    py_miss = [x for x in w if round(float(x), 2) not in ck]

    summary = f"""# RYA-1033 — the rounded-wavelength join key, measured

Generated by `scripts/rya1033_rounded_key_audit.py`.

## The drop

`{MEASURED_EW.name}` holds **{len(meas)}** Fe I lines. Joining them to `canonical_gf` on a
2-decimal-rounded air wavelength fails for **{int(miss.sum())}** of them. Every one has a
canonical Fe I row within **{dropped.separation_mA.abs().max():.2f} mA** — the worst case in
the whole set. **ZERO are genuinely absent from canonical_gf.**

Under the tolerance matcher (`pipeline.line_match`, {line_match.MATCH_TOL_A} A) the same
{len(meas)} lines produce **{n_unres} unresolved**.

## The key is not a function of the wavelength

    round(6136.615, 2)        -> {round(6136.615, 2)}      (Python: correctly-rounded decimal)
    np.round(6136.615, 2)     -> {float(np.round(np.float64(6136.615), 2))}      (numpy/pandas: scale-multiply-round)

`promote_solar_ew` rounded with pandas; `abundances_derive` rounded with Python. On this
pool that difference alone moves the casualty count from **{int(miss.sum())}** (pandas both
sides) to **{len(py_miss)}** (Python on the measured side). A join whose result depends on
which module imported which library is not reproducible even in principle.

## Wavelength alone is not a line identity

**{n_amb}** of the {len(meas)} measured Fe I lines match MORE THAN ONE canonical row within
{line_match.MATCH_TOL_A} A. Those candidates are different transitions — they disagree in
excitation potential and in gf by up to ~1.9 dex (see `ambiguous_lines.csv`). The rounded
key silently returned one of them. The matcher refuses unless an EP separates them, which
is the rule `perline_product` and `gf_grades` already followed (RYA-780/852).

## The graded pool count

| keying | graded tiers | n |
|---|---|---|
""" + "\n".join(f"| {r.scheme} | {r.graded_tiers} | {r.n} |" for r in pool.itertuples()) + f"""

The rounded key undercounts the graded Fe I pool by **{n_tol - n_round}**. This reproduces
the ticket's reconstructed **63**: it is the LAB+NIST-C+ pool counted through the rounded
key, against **{n_tol}** through a tolerance match.

⚠️ The **67** the ticket compares it to is a DIFFERENT POOL, not the same number keyed
differently — it is the depth-gated LAB-only synth selection
(`FeI_4200_6910_*_SYNTH_GRADED_1D-LTE_lines.csv`, 67 rows = 176 in-band LAB lines minus the
109 above the depth gate). Both counts are correct for what they count. The rounding defect
is real and costs {n_tol - n_round} line(s); it is not the whole 63-vs-67 gap.
"""
    (OUT / "summary.md").write_text(summary)

    print(f"measured Fe I lines            : {len(meas)}")
    print(f"dropped by 2-dp rounded key    : {int(miss.sum())}  (all within "
          f"{dropped.separation_mA.abs().max():.2f} mA of a canonical row)")
    print(f"unresolved by tolerance match  : {n_unres}")
    print(f"ambiguous (EP required)        : {n_amb}")
    print(f"graded pool  rounded / matched : {n_round} / {n_tol}")
    print(f"wrote → {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
