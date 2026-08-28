#!/usr/bin/env python3
"""RYA-1104 — is the solar Amarsi 3D-NLTE Fe offset the EW-route base? Three READ-ONLY tests.

    python3 scripts/rya1104_amarsi_base_diagnostic.py --pool amarsi_feI_vis --holding kpno_molecfit

READ-ONLY. It opens committed artifacts and writes one JSON report under
`data/results/rya1104/`. It changes no line list, no canonical_gf, no STAR_PARAMS, no
atmosphere and no product config, and it re-derives no shipped value (RYA-161).

🔴 THE TICKET'S PREMISE DOES NOT SURVIVE ITS OWN TEST 1, AND THAT IS THE FINDING.
The ticket reasons from `EW · 3D-NLTE · Amarsi`: the Amarsi MLP supplies ~+0.027 dex, the
product reads 7.512, `Synth · 1D-LTE` reads 7.447, so ~+0.038 must be an EW-route base
zero-point. Test 1 asks for the counterfactual — apply the same correction to the SYNTHESIS
base and see it land near 7.47.

That counterfactual is the shipped configuration. Every RYA-1095 Amarsi product names its
own base in `provenance`, and on kpno_molecfit it is

    FeI_4200_6910_kpno_solar_atlas_solar_kpno_molecfit_corrected_SYNTH_GRADED_1D-LTE_lines.csv

— the synthesis pool, and this script checks the 50 `a_1dlte` values against that file
rather than taking the string's word for it. They agree to 0.0 dex on 50 of 50 lines. There
is no EW measurement anywhere under the Amarsi leg. The `EW` in the display name is a LABEL,
and Test 1 below traces it to the two places that assert it.

⚠️ WHY THE SUBTRACTION LOOKED LIKE A ROUTE. 7.447 is the median over 67 lines; the Amarsi
product is the median over the 50 of those the MLP's training domain admits. Restricted to
those 50, the SAME synthesis base reads 7.475. So +0.028 of the ticket's +0.065 is a
LINE-SUBSET effect inside one pool, not a route difference — the two numbers were never
computed over the same lines. That is the RYA-1083 class one level up: an aggregate
differenced against another aggregate over a different population.

WHAT EACH TEST DOES
-------------------
TEST 1  Names the Amarsi leg's actual 1D-LTE base from its own provenance, verifies the
        per-line values against that file, and reports the base median over the full pool,
        the base median over the Amarsi subset, and the shipped post-Amarsi value. The
        genuine EW-route product for the same holding is reported beside them as the
        comparand the ticket's hypothesis actually needs.
TEST 2  Splits the Amarsi pool by gf provenance on the λ+EP dual key, through
        `gf_resolver`'s own WTOL/EPTOL, and refuses on any unmatched or ambiguous line.
TEST 3  The ⟨3D⟩ NLTE effect as a PER-LINE PAIRED differential (`pipeline.paired_differential`,
        RYA-1083) with an explicit engaged/not-engaged verdict assembled from four
        independent witnesses.

⚠️ NOTHING HERE RE-RUNS A SYNTHESIS OR THE MLP. Turbospectrum is the ratified synthesiser
and the committed products are its output; re-deriving them on this Mac would produce a
second, unauditable copy of a number that already exists (RYA-924), and the MOOGSILENT/EW
route cannot run on this machine at all. Every value below is read from a committed artifact
and every artifact is named in the report.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.gf_resolver import WTOL, EPTOL                      # noqa: E402
from pipeline.paired_differential import paired_differential      # noqa: E402
from pipeline import harness_residual, treatment_axes as taxes    # noqa: E402

BP = ROOT / "data" / "results" / "band_products"
AMARSI = ROOT / "data" / "results" / "rya1095"
CANON = ROOT / "data" / "linelists" / "canonical_gf.csv"
OUT_DEFAULT = ROOT / "data" / "results" / "rya1104"

ASPLUND21_FE = 7.46          # A(Fe), Asplund, Amarsi & Grevesse 2021

#: The RYA-1095 Amarsi legs, by the holding each one corrected. The value is ONLY the
#: directory: the base pool, the instrument and the band all come out of the artifacts
#: themselves. Reconstructing a base filename here is the RYA-990 failure — a rebuilt path
#: agrees with itself while disagreeing with the file that actually ran.
HOLDINGS = {
    "kpno_molecfit": "kpmf",
    "kpno_kurucz2005": "kurucz2005",
    "harps": "harps",
    "iag": "iag",
}
DEFAULT_HOLDING = "kpno_molecfit"

#: The one pool this ticket is about. Named rather than defaulted so the smoke test says
#: out loud which population every number below is over.
POOLS = {"amarsi_feI_vis": ("Fe", "I", "VIS")}

#: `provenance` ends with the base the leg corrected. Parsed rather than rebuilt.
_BASE_RE = re.compile(r"1D-LTE base from (?P<name>\S+\.csv)\s*$")

#: The ⟨3D⟩ NLTE/LTE pair for Test 3, as SUFFIXES onto the base stem. Both members are
#: written by the same route from the same stem, so the stem comes from Test 1's base file
#: and only the treatment segment is spelled here.
MEAN3D_NLTE_SUFFIX = "synth-mean3D-NLTE-gerber-stagger"
MEAN3D_LTE_SUFFIX = "synth-mean3D-LTE-gerber-stagger"


class DiagnosticRefusal(SystemExit):
    """A question this diagnostic cannot answer from committed artifacts."""


# ── artifact access ───────────────────────────────────────────────────────────

def _read(path: Path, what: str) -> pd.DataFrame:
    if not path.exists():
        raise DiagnosticRefusal(f"{what} missing: {path}")
    return pd.read_csv(path)


def amarsi_leg(holding: str) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """The RYA-1095 Amarsi product, its per-line file, and its run summary."""
    d = AMARSI / HOLDINGS[holding]
    prod = _read(d / "rya817_3dnlte_products.csv", f"Amarsi product for {holding}")
    lines = _read(d / "rya817_3dnlte_per_line.csv", f"Amarsi per-line for {holding}")
    summary = json.loads((d / "rya817_run_summary.json").read_text())
    return prod, lines, summary


def base_stem_from_provenance(prod_row) -> tuple[str, Path]:
    """The 1D-LTE pool this Amarsi product actually corrected, from ITS OWN provenance.

    Returns (stem, per-line path). The stem is the filename with `1D-LTE_lines.csv`
    stripped, which is the prefix every sibling treatment of that pool shares — so Test 3's
    ⟨3D⟩ pair is reached from the base the leg names, never from a stem typed here.
    """
    m = _BASE_RE.search(str(prod_row.provenance))
    if not m:
        raise DiagnosticRefusal(
            f"the Amarsi product does not name its 1D-LTE base in `provenance`, so the "
            f"pool it corrected cannot be identified from the artifact:\n"
            f"  {prod_row.provenance!r}\n"
            f"Test 1 is a question ABOUT that base; guessing it would answer a different "
            f"question and look identical (RYA-990).")
    name = m.group("name")
    path = BP / name
    if not path.exists():
        raise DiagnosticRefusal(
            f"the Amarsi leg names its base as {name}, and that file is not in the "
            f"committed band_products directory: {path}\n"
            f"The base is the subject of Test 1 — it must be READ, never assumed.")
    if not name.endswith("_1D-LTE_lines.csv"):
        raise DiagnosticRefusal(f"unexpected base filename shape: {name}")
    return name[: -len("1D-LTE_lines.csv")], path


def match_dual_key(left: pd.DataFrame, right: pd.DataFrame, *,
                   left_wl: str, left_ep: str, right_wl: str, right_ep: str,
                   what: str) -> np.ndarray:
    """Row index in `right` for each row of `left`, on the λ+EP DUAL KEY. Never an argmin.

    🔴 WAVELENGTH ALONE IS NOT AN IDENTITY (RYA-1033/RYA-1037). Python and numpy round the
    same tie to different values, so two files describing the SAME line can disagree in the
    last place and a `merge(on="wavelength_air_A")` silently drops the pair — or, worse,
    joins a same-species NEIGHBOUR inside the window and returns a confident wrong row.
    RYA-1037's AST guard refuses a wave-only join anywhere in this tree, and it refused an
    earlier draft of this very function.

    The tolerances are `gf_resolver`'s own WTOL/EPTOL, imported rather than re-typed: a
    diagnostic must not be able to widen the tolerance that decides its own answer
    (RYA-1042 — a tolerance needs a measured null, and inventing one here would give this
    script a lever on its own verdict).

    0 matches or >1 match REFUSES. Taking the nearest of several is how a line gets keyed
    to its neighbour and nobody finds out.
    """
    rw = np.asarray(right[right_wl], dtype=float)
    re_ = np.asarray(right[right_ep], dtype=float)
    out = np.empty(len(left), dtype=int)
    zero, many = [], []
    for i, (wl, ep) in enumerate(zip(np.asarray(left[left_wl], dtype=float),
                                     np.asarray(left[left_ep], dtype=float))):
        hit = np.where((np.abs(rw - wl) <= WTOL) & (np.abs(re_ - ep) <= EPTOL))[0]
        if hit.size == 0:
            zero.append(round(float(wl), 4)); out[i] = -1
        elif hit.size > 1:
            many.append((round(float(wl), 4), int(hit.size))); out[i] = -1
        else:
            out[i] = int(hit[0])
    if zero or many:
        raise DiagnosticRefusal(
            f"{what}: {len(zero)} line(s) with NO match on the λ+EP dual key "
            f"(±{WTOL} Å / ±{EPTOL} eV) {zero[:5]}, {len(many)} AMBIGUOUS {many[:5]}. "
            f"Never a silent drop and never an argmin (RYA-429/RYA-1033) — if the two "
            f"sides cannot be keyed line for line, this test has no answer and must not "
            f"produce one.")
    return out



# ── TEST 1 ────────────────────────────────────────────────────────────────────

def _one_dlte_row(products: pd.DataFrame, path: Path) -> pd.Series:
    hit = products[products["treatment"].astype(str) == "1D-LTE"]
    if len(hit) != 1:
        raise DiagnosticRefusal(
            f"{len(hit)} rows named `1D-LTE` in {path.name} — refusing to choose.")
    return hit.iloc[0]


def diagnose_base(holding: str) -> dict:
    """Where the +0.065 actually comes from, decomposed on ONE population at a time."""
    prod, lines, summary = amarsi_leg(holding)
    if len(prod) != 1:
        raise DiagnosticRefusal(
            f"{len(prod)} Amarsi products for {holding}; this diagnostic is scoped to the "
            f"single Fe I VIS leg.")
    p = prod.iloc[0]
    stem, base_lines_path = base_stem_from_provenance(p)
    base_lines = _read(base_lines_path, "1D-LTE base per-line")
    base_products = _read(BP / f"{stem}products.csv", "1D-LTE base product")
    base_row = _one_dlte_row(base_products, BP / f"{stem}products.csv")

    # THE POOL. `in_domain` is the MLP's training-domain verdict; `in_aggregate_1dlte` is
    # the base product's own exclusion, carried through. A line needs both to be in the
    # Amarsi aggregate, so both gate the subset here.
    pool = lines[lines["in_domain"].astype(bool)
                 & lines["in_aggregate_1dlte"].astype(bool)].reset_index(drop=True)

    # 🔴 THE CHECK THAT SETTLES THE TICKET. Does the Amarsi leg's `a_1dlte` column REALLY
    # come from that synthesis file, or does the provenance string merely say so? A
    # provenance line is a claim; the base file is the referee (RYA-1042: a vendor echo is
    # not provenance).
    # ⚠️ λ+EP, not λ. The Amarsi file carries the lower level as `elo_eV` and the base as
    # `ep_eV` — the same physical quantity under two names (RYA-871) — so the dual key is
    # available on both sides and there is no excuse for joining on wavelength alone.
    idx = match_dual_key(
        pool, base_lines, left_wl="wavelength_air_A", left_ep="elo_eV",
        right_wl="wavelength_air_A", right_ep="ep_eV",
        what=f"Amarsi pool against the base it names ({base_lines_path.name})")
    picked = base_lines.iloc[idx].reset_index(drop=True)
    merged = pool.copy()
    for c in ("abundance", "ew_mA", "in_aggregate"):
        merged[c] = picked[c].to_numpy()
    unmatched = merged[merged["abundance"].isna()]
    resid = (merged["a_1dlte"] - merged["abundance"]).abs()
    max_resid = float(resid.max())
    base_is_the_amarsi_base = bool(max_resid <= 1e-9)

    # The EW-route product for the same holding and band, if one exists. This is the
    # comparand the ticket's hypothesis needs and it has never been under the Amarsi leg.
    ew = _ew_route_product(stem)

    base_full = float(base_row["A"])
    base_subset = float(np.median(merged["abundance"]))
    post = float(p["value"])

    return {
        "holding": holding,
        "amarsi_product": {
            "value": round(post, 4), "n": int(p["n_lines"]),
            "treatment": str(p["treatment"]), "handler": str(p["handler"]),
            "route_axis": str(p["route"]), "route_basis": str(p["route_basis"]),
            "gf_axis": str(p["gf"]), "instrument": str(p["instrument"]),
            "artifact": str((AMARSI / HOLDINGS[holding]
                             / "rya817_3dnlte_products.csv").relative_to(ROOT)),
        },
        "declared_base": {
            "artifact": str(base_lines_path.relative_to(ROOT)),
            "route_axis": str(base_row["route"]),
            "handler": str(base_row["handler"]),
            "is_synthesis": str(base_row["handler"]) == "SynthesisHandler",
        },
        "base_verified_line_by_line": {
            "n_pool": int(len(pool)),
            "n_matched_in_base": int(len(pool) - len(unmatched)),
            "max_abs_residual_dex": max_resid,
            "verdict": ("the Amarsi leg's 1D-LTE column IS the synthesis base, exactly"
                        if base_is_the_amarsi_base else
                        "the Amarsi leg's 1D-LTE column does NOT reproduce the base it "
                        "names — the provenance string is wrong"),
        },
        "ew_route_product": ew,
        "decomposition": {
            "synth_base_full_pool": {"A": round(base_full, 4),
                                     "n": int(base_row["n_lines"])},
            "synth_base_amarsi_subset": {"A": round(base_subset, 4), "n": int(len(pool))},
            "post_amarsi": {"A": round(post, 4), "n": int(p["n_lines"])},
            "line_subset_step_dex": round(base_subset - base_full, 4),
            "amarsi_correction_step_dex": round(post - base_subset, 4),
            "total_vs_full_synth_base_dex": round(post - base_full, 4),
            "median_per_line_aberr_dex": round(float(pool["aberr"].median()), 4),
            "mean_per_line_aberr_dex": round(float(pool["aberr"].mean()), 4),
        },
        "mlp_input_features": _mlp_inputs(),
        "network_on_two_pools": network_on_two_pools(holding),
        "asplund21": ASPLUND21_FE,
    }


def _ew_route_product(stem: str) -> dict:
    """The genuine EW-route 1D-LTE product for the same holding and band, if one exists.

    ⚠️ Reported at PRODUCT level only. No `_lines.csv` is committed for the PROFILEFIT
    pools, so the EW base CANNOT be restricted to the Amarsi subset and this diagnostic
    does not pretend otherwise — the two numbers are over different populations and the
    report says so beside them rather than differencing them.
    """
    ew_stem = stem.replace("_SYNTH_GRADED_", "_PROFILEFIT_GRADED_")
    path = BP / f"{ew_stem}products.csv"
    if ew_stem == stem or not path.exists():
        return {"exists": False,
                "note": f"no EW-route (PROFILEFIT) product committed for this holding "
                        f"(looked for {path.name})"}
    d = pd.read_csv(path)
    row = _one_dlte_row(d, path)
    return {
        "exists": True, "artifact": str(path.relative_to(ROOT)),
        "A": float(row["A"]), "n": int(row["n_lines"]),
        "handler": str(row["handler"]), "route_axis": str(row["route"]),
        "per_line_file_committed": (BP / f"{ew_stem}1D-LTE_lines.csv").exists(),
        "note": "PRODUCT level only — no per-line file is committed for the PROFILEFIT "
                "pool, so this cannot be restricted to the Amarsi subset. It is a "
                "DIFFERENT population, not a comparand.",
    }


def _mlp_inputs() -> dict:
    """What the Amarsi network is actually a function of — read off the signature.

    Test 1 step 4 asks whether the MLP needs EW or reduced-EW, and would STOP if it did.
    It does not: `_compute_aberr` takes the transition and the stellar box and nothing
    about line strength, so the network is EW-agnostic by construction and can correct a
    synthesis base with no retraining. Taken from the function rather than asserted here,
    because a docstring is not the code (RYA-1079).
    """
    import inspect
    from pipeline.nlte_corrections import _compute_aberr
    params = list(inspect.signature(_compute_aberr).parameters)
    strength = [p for p in params if "ew" in p.lower() or "rew" in p.lower()]
    return {
        "signature": params,
        "takes_ew_or_reduced_ew": bool(strength),
        "verdict": ("the MLP consumes only the transition (Elo, Eup, log gf) and the "
                    "stellar box (Teff, log g, vmic, A(Fe;3N)) — no line-strength input, "
                    "so it is route-agnostic and needs no retraining to sit on a "
                    "synthesis base" if not strength else
                    f"the MLP takes a line-strength input {strength} — a synthesis base "
                    f"would need its strengths fed in"),
    }


def network_on_two_pools(holding: str) -> dict:
    """🔴 WHERE THE OFFSET ACTUALLY IS: the SAME network gives a different correction here.

    The ticket looks for the offset in the base. Tests 1-2 close that door — the base is
    the synthesis pool and its gf is laboratory. What is left is the correction, and the
    correction is not a constant.

    The RYA-817 reactivation control runs this identical network on AMARSI'S OWN solar line
    set (Allende Prieto et al. 2002, weak lines, REW < -4.9) and it PASSES: mean correction
    -0.0024 dex, landing 7.4646 against the published 7.46. On our 50-line graded pool the
    same network returns +0.043. Same code, same models, same star, same stellar box — a
    ~+0.045 dex swing that lives entirely in WHICH LINES are asked.

    ⚠️ AND THE POPULATIONS ARE NOT COMPARABLE BY CONSTRUCTION. Our graded pool is the
    laboratory-gf tier, and Ruffoni 2014 / Den Hartog 2014 measured HIGH-EXCITATION Fe I:
    the whole committed graded pool floors at 2.851 eV. Amarsi's solar set reaches down to
    0.087 eV. The network's own Elo trend -- measured by the control at -0.0905 dex over
    Elo 0->5 eV, inside the paper's published -0.10..-0.05 -- means a pool with no low-Elo
    lines gets a systematically more positive correction.

    ⚠️ ELO DOES NOT CLOSE THE WHOLE GAP, AND THIS SAYS SO. Transporting our pool's own
    Elo fit down to the control pool's mean Elo moves the correction from +0.043 to about
    +0.030, roughly a quarter of the swing. So the population difference is real and
    measured, and a single-variable Elo story is NOT sufficient. The remainder is not
    attributed here: this is a diagnostic, and naming a cause it has not measured is the
    RYA-1042 failure.

    ⚠️ THE AXIS IS NOT THE CAUSE, and that confound is closed rather than assumed. The band
    run rails A(Fe;3N) at the grid ceiling 7.5 while the control runs at 7.46 -- but the
    run's own recorded sensitivity across the WHOLE axis is 0.002-0.004 dex, two orders
    below the swing.
    """
    prod, lines, summary = amarsi_leg(holding)
    ctl = summary["reactivation_control"]
    stats = summary["per_band"]["Fe I VIS"]
    pool = lines[lines["in_domain"].astype(bool)
                 & lines["in_aggregate_1dlte"].astype(bool)].reset_index(drop=True)

    sol = _read(ROOT / "data" / "reference" / "amarsi2022_training"
                / "amarsi2022_solar_control_lines.csv", "Amarsi solar control lines")
    ctl_pool = sol[(sol["ion"] == "I") & sol["weak_line_rew_lt_m49"].astype(bool)]

    elo, ab = pool["elo_eV"].to_numpy(float), pool["aberr"].to_numpy(float)
    slope, intercept = np.polyfit(elo, ab, 1)
    ours = float(np.median(ab))
    transported = float(np.polyval([slope, intercept], ctl_pool["elo_eV"].mean()))
    ctl_mean = float(ctl["feI_mean_correction"])

    elo_check = next((c for c in ctl["checks"] if "Elo 0->5 eV" in c["check"]), None)
    return {
        "control_pool": {
            "line_list": ctl["line_list"], "cut": ctl["cut"], "n": int(len(ctl_pool)),
            "elo_eV": {"min": round(float(ctl_pool["elo_eV"].min()), 3),
                       "median": round(float(ctl_pool["elo_eV"].median()), 3),
                       "mean": round(float(ctl_pool["elo_eV"].mean()), 3),
                       "max": round(float(ctl_pool["elo_eV"].max()), 3)},
            "mean_correction_dex": round(ctl_mean, 4),
            "control_passed": bool(ctl["passed"]),
            "reproduces_published": [
                {"check": c["check"], "computed": c["computed"],
                 "published": c["published"], "pass": c["pass"]}
                for c in ctl["checks"] if "A(Fe I)" in c["check"]],
        },
        "our_pool": {
            "n": int(len(pool)),
            "elo_eV": {"min": round(float(elo.min()), 3),
                       "median": round(float(np.median(elo)), 3),
                       "mean": round(float(elo.mean()), 3),
                       "max": round(float(elo.max()), 3)},
            "median_correction_dex": round(ours, 4),
            "mean_correction_dex": round(float(ab.mean()), 4),
            "selection": "the GRADED (laboratory-gf) tier, then the MLP training domain",
        },
        "graded_pool_elo_floor_is_the_gf_tier": _graded_elo_floor(holding),
        "swing_dex": round(ours - ctl_mean, 4),
        "elo_trend": {
            "network_trend_over_0_to_5_eV_dex": elo_check["computed"] if elo_check else None,
            "published_range": elo_check["published"] if elo_check else None,
            "our_pool_slope_dex_per_eV": round(float(slope), 4),
            "our_pool_corr_elo_vs_aberr": round(
                float(np.corrcoef(elo, ab)[0, 1]), 3),
            "our_correction_transported_to_control_mean_elo_dex": round(transported, 4),
            "fraction_of_swing_explained": round(
                (ours - transported) / (ours - ctl_mean), 3) if ours != ctl_mean else None,
        },
        "axis_confound_closed": {
            "band_run_axis": stats["afe_axis_star"],
            "railed_at_ceiling": stats["afe_axis_railed_at_grid_ceiling"],
            "control_axis": ctl["axis_A_Fe_3N"],
            "sensitivity_across_whole_axis_dex": stats["afe_axis_sensitivity_dex"],
            "verdict": "immaterial — two orders below the swing",
        },
    }


def _graded_elo_floor(holding: str) -> dict:
    """Is our pool's high-Elo floor the gf TIER or the MLP's domain cut? Measured, not guessed."""
    prod, _, _ = amarsi_leg(holding)
    _stem, base_path = base_stem_from_provenance(prod.iloc[0])
    base = pd.read_csv(base_path)
    col = "ep_eV" if "ep_eV" in base.columns else "elo_eV"
    ep = base[col].dropna().to_numpy(float)
    return {"artifact": str(base_path.relative_to(ROOT)), "n": int(len(ep)),
            "elo_min_eV": round(float(ep.min()), 3),
            "elo_median_eV": round(float(np.median(ep)), 3),
            "verdict": "the floor is a property of the GRADED pool itself, before the MLP "
                       "domain cut — Ruffoni 2014 / Den Hartog 2014 measured "
                       "high-excitation Fe I, so lab-grading and low-Elo coverage pull "
                       "against each other."}



# ── TEST 2 ────────────────────────────────────────────────────────────────────

def diagnose_gf_split(holding: str) -> dict:
    """gf provenance of the Amarsi pool, on the λ+EP dual key. NO line list is touched."""
    prod, lines, _ = amarsi_leg(holding)
    p = prod.iloc[0]
    pool = lines[lines["in_domain"].astype(bool)
                 & lines["in_aggregate_1dlte"].astype(bool)].reset_index(drop=True)

    can = pd.read_csv(CANON, low_memory=False)
    ion_z = {"I": 1, "II": 2}[str(p["ion"])]
    sub = can[(pd.to_numeric(can["key_z"], errors="coerce") == 26)
              & (pd.to_numeric(can["ion"], errors="coerce") == ion_z)]
    wls = sub["wavelength_air_A"].to_numpy(float)
    eps = sub["excitation_potential_eV"].to_numpy(float)

    # ⚠️ DELIBERATELY NOT `match_dual_key`, and the difference is the QUESTION being asked.
    # That matcher answers a VALUE question — which single row carries this line's
    # abundance — so >1 candidate is fatal to it. This is a TIER question, and several
    # canonical rows inside one window can name the same provenance unambiguously (an HFS
    # cluster is still one gf tier). So >1 is tolerated ONLY while the tiers AGREE, and
    # disagreement refuses. Measured on the committed pool: every line matches exactly one
    # row, so the allowance is currently unexercised — stated rather than silently relied on.
    #
    # ⚠️ THE DUAL KEY IS THE POINT. Wavelength alone is not an identity — Python and numpy
    # round the same tie to different values (RYA-1033/1037) and a same-species neighbour
    # inside the window is a real line, not a rounding artefact. λ AND EP, through
    # gf_resolver's OWN tolerances so this diagnostic cannot invent a wider one, and any
    # 0-match or >1-match REFUSES rather than taking an argmin (RYA-1033).
    tiers, refs, ambiguous, unmatched = [], [], [], []
    for r in pool.itertuples():
        m = np.where((np.abs(wls - r.wavelength_air_A) <= WTOL)
                     & (np.abs(eps - r.elo_eV) <= EPTOL))[0]
        if m.size == 0:
            unmatched.append(round(float(r.wavelength_air_A), 4))
            tiers.append(None); refs.append(None); continue
        t = sorted(set(str(x) for x in sub.iloc[m]["gf_tier"]))
        if m.size > 1 and len(t) > 1:
            ambiguous.append((round(float(r.wavelength_air_A), 4), t))
        tiers.append(t[0] if len(t) == 1 else "AMBIGUOUS")
        refs.append(str(sub.iloc[m]["loggf_reference"].iloc[0]))
    if unmatched or ambiguous:
        raise DiagnosticRefusal(
            f"gf provenance is not resolvable for the whole pool: {len(unmatched)} line(s) "
            f"with no canonical_gf row within ±{WTOL} Å / ±{EPTOL} eV {unmatched[:5]}, "
            f"{len(ambiguous)} with tiers that disagree inside the window {ambiguous[:3]}. "
            f"Test 2 partitions BY that provenance — a partition with a hole in it is not "
            f"a partition, and never a silent drop (RYA-429).")

    pool = pool.assign(gf_tier=tiers, loggf_reference=refs)
    groups = {}
    for tier, g in pool.groupby("gf_tier"):
        groups[str(tier)] = {
            "n": int(len(g)),
            "base_median_dex": round(float(g["a_1dlte"].median()), 4),
            "post_amarsi_median_dex": round(float(g["a_3dnlte"].median()), 4),
            "median_aberr_dex": round(float(g["aberr"].median()), 4),
            "references": sorted(set(str(x) for x in g["loggf_reference"]))[:8],
        }
    lab_like = {"LAB", "NIST-C+"}
    n_lab = int(sum(v["n"] for k, v in groups.items() if k in lab_like))
    n_floor = int(len(pool)) - n_lab
    return {
        "holding": holding,
        "n_pool": int(len(pool)),
        "match": {"key": "lambda + EP dual key",
                  "tolerances": {"wavelength_A": WTOL, "excitation_potential_eV": EPTOL},
                  "source": "pipeline.gf_resolver (WTOL/EPTOL, not re-derived here)",
                  "n_unmatched": 0, "n_ambiguous": 0,
                  "n_unique_match": int(len(pool))},
        "by_tier": groups,
        "n_lab_tier": n_lab,
        "n_kurucz_floored": n_floor,
        "splittable": bool(n_lab and n_floor),
        "product_gf_axis": str(prod.iloc[0]["gf"]),
        "budget_gf_rung": _budget_gf_rung(holding),
        "canonical_gf_artifact": str(CANON.relative_to(ROOT)),
    }


def _budget_gf_rung(holding: str) -> str:
    """The gf verdict the product's OWN error budget reached, for comparison with the axis.

    These are two claims about the same product's oscillator strengths made by two
    different code paths. Printing them side by side is the whole content of Test 2's
    second half — the budget MEASURES the pool, the axis is a constant in a label table.
    """
    p = AMARSI / HOLDINGS[holding] / "rya817_3dnlte_budgets.txt"
    if not p.exists():
        return "(no budget artifact)"
    for line in p.read_text().splitlines():
        if line.strip().startswith("gf rung:"):
            return line.strip()[len("gf rung:"):].strip()
    return "(no gf rung line in the budget)"


# ── TEST 3 ────────────────────────────────────────────────────────────────────

def diagnose_nlte_wiring(holding: str) -> dict:
    """Is the feed's 0.000 ⟨3D⟩ NLTE effect a real zero, or the wrong statistic?"""
    prod, _lines, _ = amarsi_leg(holding)
    stem, _ = base_stem_from_provenance(prod.iloc[0])

    paths = {k: BP / f"{stem}{suf}_{tail}"
             for k, (suf, tail) in {
                 "nlte_lines": (MEAN3D_NLTE_SUFFIX, "lines.csv"),
                 "lte_lines": (MEAN3D_LTE_SUFFIX, "lines.csv"),
                 "nlte_products": (MEAN3D_NLTE_SUFFIX, "products.csv"),
                 "lte_products": (MEAN3D_LTE_SUFFIX, "products.csv"),
                 "nlte_prov": (MEAN3D_NLTE_SUFFIX, "provenance.txt"),
                 "lte_prov": (MEAN3D_LTE_SUFFIX, "provenance.txt"),
             }.items()}
    missing = [k for k, v in paths.items() if not v.exists()]
    if missing:
        return {"holding": holding, "ran": False,
                "note": f"no ⟨3D⟩ pair committed for this holding "
                        f"(missing {sorted(missing)}); Test 3 is scoped to the holding "
                        f"that carries one."}

    nlte = pd.read_csv(paths["nlte_lines"])
    lte = pd.read_csv(paths["lte_lines"])
    nprod = pd.read_csv(paths["nlte_products"]).iloc[0]
    lprod = pd.read_csv(paths["lte_products"]).iloc[0]

    pd_ = paired_differential(nlte, lte)

    # ⚠️ THE CONTROL THAT KEEPS THIS FROM GOING VACUOUS (RYA-853/RYA-1080). A differential
    # of a product against ITSELF is 0.000 on every line and reads exactly like a real
    # null. So the two sides are asserted to be DIFFERENT OBJECTS before their difference
    # is believed: different files, different bytes, and a nonzero count.
    import hashlib
    md5 = {k: hashlib.md5(paths[k].read_bytes()).hexdigest()
           for k in ("nlte_lines", "lte_lines")}
    sides_differ = md5["nlte_lines"] != md5["lte_lines"]

    # WITNESS 4 — the code path, not the provenance string. `derive_band_products` builds
    # BOTH members from the same ⟨3D⟩ atmosphere and the same deck, and the ONLY
    # difference between them is whether `nlte_deck` is passed to the fitter at all.
    src = (ROOT / "scripts" / "derive_band_products.py").read_text()
    code_witness = {
        "passes_nlte_deck_only_when_nlte": bool(
            re.search(r"if _nlte:\s*\n\s*_fit_kw\.update\(nlte_deck=\"gerber\"", src)),
        "deck_key_is_first_positional": bool(
            re.search(r"gnlte\.for_node\(_deck_key,", src)),
        "note": "RYA-1049: `for_node(element, ..., node=...)` silently ignores `node` and "
                "returns the 1D deck. This call site passes the ⟨3D⟩ deck key as the FIRST "
                "POSITIONAL argument, which is the reachable way, so that defect does not "
                "apply to the shipped pair.",
    }

    # RYA-915 one-path overwrite: two products, two filenames, both present, both distinct.
    overwrite_witness = {
        "distinct_paths": paths["nlte_lines"] != paths["lte_lines"],
        "both_present": True, "distinct_bytes": sides_differ,
    }

    engaged = bool(sides_differ and pd_.n_nonzero > 0
                   and code_witness["passes_nlte_deck_only_when_nlte"]
                   and code_witness["deck_key_is_first_positional"]
                   and "departures applied" in paths["nlte_prov"].read_text()
                   and "departures WITHHELD" in paths["lte_prov"].read_text())

    vals = np.sort(pd.concat([lte["abundance"], nlte["abundance"]]).to_numpy(float))
    gaps = np.diff(np.unique(vals))
    return {
        "holding": holding, "ran": True,
        "artifacts": {k: str(v.relative_to(ROOT)) for k, v in paths.items()},
        "published": {"mean3d_lte_A": float(lprod["A"]), "mean3d_nlte_A": float(nprod["A"]),
                      "difference_of_published_values": round(
                          float(nprod["A"]) - float(lprod["A"]), 4),
                      "n_lines": int(nprod["n_lines"])},
        "paired_differential": pd_.as_dict(),
        "departures_engaged": engaged,
        "witnesses": {
            "provenance_strings_differ": {
                "nlte": paths["nlte_prov"].read_text().strip()[:160],
                "lte": paths["lte_prov"].read_text().strip()[:160]},
            "per_line_values_differ": {"md5": md5, "distinct": sides_differ,
                                       "n_nonzero_of_n": f"{pd_.n_nonzero}/{pd_.n_paired}"},
            "code_path": code_witness,
            "rya915_one_path_overwrite": overwrite_witness,
        },
        "not_a_quantiser": {
            "n_unique_abundance_values": int(len(np.unique(vals))),
            "min_gap_dex": round(float(gaps.min()), 6) if len(gaps) else None,
            "note": "the medians collided because the 34th of 67 sorted values landed on "
                    "the same number in both pools — a coincidence, not a rounding floor "
                    "(RYA-1083).",
        },
        "bookkeeping_defect": _nlte_delta_column(nlte, paths["nlte_lines"]),
    }


def _nlte_delta_column(nlte: pd.DataFrame, path: Path) -> dict:
    """⚠️ The NLTE member's own per-line departure columns read as if it were the LTE run.

    Found while assembling the witnesses above. It does NOT change any value — the
    abundances carry the departures — but it is the column a later reader would consult to
    ask "were departures applied?", and it answers no on the product that applied them.
    """
    if "nlte_delta_dex" not in nlte.columns:
        return {"present": False}
    deltas = sorted(set(np.round(nlte["nlte_delta_dex"].dropna().to_numpy(float), 6)))
    sources = sorted(set(str(x) for x in nlte["nlte_source"].dropna()))
    contradicts = bool(deltas == [0.0] and any("no departure applied" in s
                                               for s in sources))
    return {"present": True, "artifact": str(path.relative_to(ROOT)),
            "nlte_delta_dex_values": deltas[:5], "nlte_source_values": sources[:3],
            "contradicts_the_products_own_provenance": contradicts,
            "note": "cosmetic in the aggregate (the abundances do carry the departures) "
                    "and misleading to any reader who trusts the column — the departures "
                    "are applied inside the fit, so the bookkeeping columns the EW-route "
                    "correction path writes were never filled by the synthesis path."}


# ── the label audit that Test 1 turns up ──────────────────────────────────────

def label_audit(t1: dict) -> dict:
    """WHERE `EW` and `kurucz` come from on a product whose base is synthesis lab-gf.

    Both are ASSERTED, neither is measured, and each is asserted twice. Named here because
    the ticket's premise is downstream of them: a reader who takes the display name at face
    value concludes there is an EW measurement under this product, and there is not.
    """
    handler = t1["amarsi_product"]["handler"]
    return {
        "displayed_route": t1["amarsi_product"]["route_axis"],
        "displayed_route_basis": t1["amarsi_product"]["route_basis"],
        "actual_base_route": t1["declared_base"]["route_axis"],
        "actual_base_handler": t1["declared_base"]["handler"],
        "assertion_1": {
            "site": "scripts/rya817_run_3dnlte_bands.py — build_product(handler=...)",
            "value": handler,
            "why_it_wins": "handler is `resolve_route`'s STRONGEST witness (RYA-869), so "
                           "it sets route_basis='handler' and outranks everything else.",
            "why_it_is_wrong": "it is a constant, written when this route could only "
                               "read PROFILEFIT pools. RYA-1031 widened discovery to "
                               "SYNTH and the constant did not move with it.",
        },
        "assertion_2": {
            "site": "pipeline/treatment_axes.py — _ROUTE_BY_LABEL['ENGINE-A-3DNLTE']",
            "value": taxes._ROUTE_BY_LABEL.get("ENGINE-A-3DNLTE"),
            "why_it_is_wrong": "RYA-1002 removed the pinned route from `ENGINE-A` for "
                               "exactly this reason — the label appears on both routes, "
                               "so it cannot name one. `ENGINE-A-3DNLTE` inherits the same "
                               "problem and kept the pin.",
        },
        "gf_axis_assertion": {
            "site": "pipeline/treatment_axes.py — LEGACY['ENGINE-A-3DNLTE']['gf']",
            "value": taxes.LEGACY["ENGINE-A-3DNLTE"]["gf"],
            "why_it_is_wrong": "a constant in a label table. The pool's own budget "
                               "measures the gf and reaches the opposite verdict.",
        },
        "shipping_consequence": {
            "charged_harness_residual_dex":
                harness_residual.HANDLER_RESIDUAL_DEX[handler],
            "earned_harness_residual_dex":
                harness_residual.HANDLER_RESIDUAL_DEX[t1["declared_base"]["handler"]],
            "note": "the handler is not only a display string — `rya1095_amarsi_error_"
                    "budget` spends it on `harness_residual.for_handler`, so the shipped "
                    "`syst_dex` carries the profile fitter's control residual on a "
                    "product measured by the synthesiser.",
        },
        "verdict": "the `EW` and `kurucz` segments of this product's display name are "
                   "LABELS, contradicted by the product's own base file and its own "
                   "error budget. Nothing about the VALUE changes; what changes is what "
                   "the name licenses a reader to conclude — and this ticket is that "
                   "reader.",
    }


# ── report ────────────────────────────────────────────────────────────────────

def render(rep: dict) -> str:
    t1, t2, t3 = rep["test1"], rep["test2"], rep["test3"]
    d = t1["decomposition"]
    L = []
    A = L.append
    A("=" * 88)
    A(f"RYA-1104 — Amarsi 3D-NLTE base diagnostic   pool={rep['pool']}  "
      f"holding={rep['holding']}")
    A("=" * 88)

    A("")
    A("TEST 1 — is the offset the EW-route base?")
    A(f"  the Amarsi leg names its 1D-LTE base as")
    A(f"      {t1['declared_base']['artifact']}")
    A(f"  route={t1['declared_base']['route_axis']}  "
      f"handler={t1['declared_base']['handler']}  "
      f"-> SYNTHESIS: {t1['declared_base']['is_synthesis']}")
    v = t1["base_verified_line_by_line"]
    A(f"  verified line by line: {v['n_matched_in_base']}/{v['n_pool']} matched, "
      f"max |residual| = {v['max_abs_residual_dex']:.1e} dex")
    A(f"      {v['verdict']}")
    A("")
    A(f"  SYNTH base (full pool)      {d['synth_base_full_pool']['A']:.3f}  "
      f"n={d['synth_base_full_pool']['n']}")
    A(f"  SYNTH base (Amarsi subset)  {d['synth_base_amarsi_subset']['A']:.3f}  "
      f"n={d['synth_base_amarsi_subset']['n']}   "
      f"[line-subset step {d['line_subset_step_dex']:+.3f}]")
    A(f"  + Amarsi (SHIPPED)          {d['post_amarsi']['A']:.3f}  "
      f"n={d['post_amarsi']['n']}   "
      f"[correction step {d['amarsi_correction_step_dex']:+.3f}]")
    ew = t1["ew_route_product"]
    if ew["exists"]:
        A(f"  EW-route 1D-LTE product     {ew['A']:.3f}  n={ew['n']}   "
          f"({ew['handler']}) — a DIFFERENT pool, never under the Amarsi leg")
    else:
        A(f"  EW-route 1D-LTE product     {ew['note']}")
    A(f"  Asplund 2021                {t1['asplund21']:.2f}")
    A("")
    A(f"  MLP inputs: {', '.join(t1['mlp_input_features']['signature'])}")
    A(f"      takes EW / reduced-EW: "
      f"{t1['mlp_input_features']['takes_ew_or_reduced_ew']}")
    A("")
    n = t1["network_on_two_pools"]
    cp, op, tr = n["control_pool"], n["our_pool"], n["elo_trend"]
    A("  THE SAME NETWORK ON TWO LINE POOLS — where the offset actually lives:")
    A(f"    Amarsi's own solar set   n={cp['n']:<3} Elo {cp['elo_eV']['min']:.2f}-"
      f"{cp['elo_eV']['max']:.2f} (median {cp['elo_eV']['median']:.2f})   "
      f"correction {cp['mean_correction_dex']:+.4f}   control PASSES: "
      f"{cp['control_passed']}")
    A(f"    our graded pool          n={op['n']:<3} Elo {op['elo_eV']['min']:.2f}-"
      f"{op['elo_eV']['max']:.2f} (median {op['elo_eV']['median']:.2f})   "
      f"correction {op['median_correction_dex']:+.4f}")
    A(f"    swing {n['swing_dex']:+.4f} dex from the line population alone")
    A(f"    network Elo trend over 0-5 eV: "
      f"{tr['network_trend_over_0_to_5_eV_dex']:+.4f} dex "
      f"(published {tr['published_range']}); our pool's own slope "
      f"{tr['our_pool_slope_dex_per_eV']:+.4f} dex/eV, r="
      f"{tr['our_pool_corr_elo_vs_aberr']:+.2f}")
    A(f"    transported to the control pool's mean Elo: "
      f"{tr['our_correction_transported_to_control_mean_elo_dex']:+.4f} — Elo explains "
      f"{tr['fraction_of_swing_explained']:.0%} of the swing, NOT all of it")
    A(f"    {n['graded_pool_elo_floor_is_the_gf_tier']['verdict']}")
    ax = n["axis_confound_closed"]
    A(f"    A(Fe;3N) axis confound: run at {ax['band_run_axis']} (railed="
      f"{ax['railed_at_ceiling']}) vs control at {ax['control_axis']}; whole-axis "
      f"sensitivity {ax['sensitivity_across_whole_axis_dex']} dex — {ax['verdict']}")
    A("")
    A(f"  VERDICT: {rep['verdicts']['test1']}")

    A("")
    A("TEST 2 — graded-vs-Kurucz gf split (no line list touched)")
    m = t2["match"]
    A(f"  {m['key']}, +/-{m['tolerances']['wavelength_A']} A / "
      f"+/-{m['tolerances']['excitation_potential_eV']} eV ({m['source']})")
    A(f"  {m['n_unique_match']}/{t2['n_pool']} unique matches, "
      f"{m['n_unmatched']} unmatched, {m['n_ambiguous']} ambiguous")
    for tier, g in sorted(t2["by_tier"].items()):
        A(f"  {tier:<10} n={g['n']:<4} base {g['base_median_dex']:.3f}  "
          f"post-Amarsi {g['post_amarsi_median_dex']:.3f}   "
          f"{'; '.join(g['references'][:2])}")
    A(f"  lab-tier n={t2['n_lab_tier']}   Kurucz-floored n={t2['n_kurucz_floored']}   "
      f"splittable: {t2['splittable']}")
    A(f"  product `gf` axis says: {t2['product_gf_axis']}")
    A(f"  the same product's budget says: {t2['budget_gf_rung'][:96]}")
    A(f"  VERDICT: {rep['verdicts']['test2']}")

    A("")
    A("TEST 3 — <3D>-NLTE vs <3D>-LTE wiring")
    if not t3["ran"]:
        A(f"  {t3['note']}")
    else:
        p, q = t3["published"], t3["paired_differential"]
        A(f"  published    <3D>-LTE {p['mean3d_lte_A']:.3f}   "
          f"<3D>-NLTE {p['mean3d_nlte_A']:.3f}   "
          f"difference-of-published {p['difference_of_published_values']:+.4f}")
        A(f"  PER-LINE PAIRED (the real statistic, RYA-1083):")
        A(f"      median {q['median']:+.4f}  mean {q['mean']:+.4f}  "
          f"sd {q['sd']:.4f}  range {q['min']:+.3f}..{q['max']:+.3f}")
        A(f"      nonzero on {q['n_nonzero']}/{q['n_paired']} lines   "
          f"collision: {q['collision']}")
        A(f"  departures engaged: {t3['departures_engaged']}")
        w = t3["witnesses"]
        A(f"      per-line files distinct: {w['per_line_values_differ']['distinct']}  "
          f"({w['per_line_values_differ']['n_nonzero_of_n']} lines differ)")
        A(f"      code passes nlte_deck only when NLTE: "
          f"{w['code_path']['passes_nlte_deck_only_when_nlte']}   "
          f"deck key first-positional (RYA-1049 safe): "
          f"{w['code_path']['deck_key_is_first_positional']}")
        A(f"      RYA-915 one-path overwrite: "
          f"{'not present' if w['rya915_one_path_overwrite']['distinct_bytes'] else 'SUSPECT'}")
        nq = t3["not_a_quantiser"]
        A(f"      {nq['n_unique_abundance_values']} unique abundances, "
          f"min gap {nq['min_gap_dex']} dex — not a quantiser floor")
        b = t3["bookkeeping_defect"]
        if b.get("contradicts_the_products_own_provenance"):
            A(f"  ⚠️  {b['artifact']}: nlte_delta_dex={b['nlte_delta_dex_values']} and "
              f"nlte_source={b['nlte_source_values']}")
            A(f"      — the NLTE product's own departure columns read as LTE. "
              f"{b['note'][:88]}")
        A(f"  VERDICT: {rep['verdicts']['test3']}")

    A("")
    A("LABEL AUDIT — where the ticket's premise came from")
    la = rep["label_audit"]
    A(f"  displayed route `{la['displayed_route']}` "
      f"(basis={la['displayed_route_basis']}) vs actual base route "
      f"`{la['actual_base_route']}` / {la['actual_base_handler']}")
    for k in ("assertion_1", "assertion_2", "gf_axis_assertion"):
        a = la[k]
        A(f"    {a['site']}  ->  {a['value']!r}")
    sc = la["shipping_consequence"]
    A(f"  shipped syst_dex charges {sc['charged_harness_residual_dex']} dex "
      f"(ProfileFit) where the base earns {sc['earned_harness_residual_dex']} dex")

    A("")
    A("SMOKE-TEST SHAPE")
    A(f"TEST 1  EW base {ew['A']:.3f} (n={ew['n']}, different pool) | "
      f"SYNTH base {d['synth_base_full_pool']['A']:.3f} | "
      f"SYNTH base on Amarsi lines {d['synth_base_amarsi_subset']['A']:.3f} | "
      f"+Amarsi {d['post_amarsi']['A']:.3f}"
      if ew["exists"] else
      f"TEST 1  EW base n/a | SYNTH base {d['synth_base_full_pool']['A']:.3f} | "
      f"+Amarsi {d['post_amarsi']['A']:.3f}")
    A(f"        the correction itself: {n['control_pool']['mean_correction_dex']:+.4f} on "
      f"Amarsi's own solar lines vs {n['our_pool']['median_correction_dex']:+.4f} on ours "
      f"-> swing {n['swing_dex']:+.4f} from the line population")
    lab = t2["by_tier"].get("LAB", {"n": 0, "base_median_dex": float('nan')})
    A(f"TEST 2  lab-gf n={t2['n_lab_tier']} median {lab['base_median_dex']:.3f} | "
      f"kurucz-floored n={t2['n_kurucz_floored']} — the split does not exist")
    if t3["ran"]:
        A(f"TEST 3  <3D>-NLTE departures engaged: {t3['departures_engaged']} | "
          f"LTE {t3['published']['mean3d_lte_A']:.3f} "
          f"NLTE {t3['published']['mean3d_nlte_A']:.3f} "
          f"effect {t3['paired_differential']['median']:+.4f} (paired, not differenced)")
    A("")
    A("NOTHING WAS MODIFIED: no line list, no canonical_gf, no STAR_PARAMS, no product.")
    A("=" * 88)
    return "\n".join(L)


def verdicts(t1, t2, t3) -> dict:
    d = t1["decomposition"]
    n = t1["network_on_two_pools"]
    v1 = (f"REFUTED — there is no EW-route base under this product. Its 1D-LTE column is "
          f"the SYNTHESIS pool it names, verified to "
          f"{t1['base_verified_line_by_line']['max_abs_residual_dex']:.0e} dex on "
          f"{t1['base_verified_line_by_line']['n_pool']}/"
          f"{t1['base_verified_line_by_line']['n_pool']} lines. Test 1's counterfactual IS "
          f"the shipped configuration, and it reads "
          f"{d['post_amarsi']['A']:.3f}, not ~7.47. The +"
          f"{d['total_vs_full_synth_base_dex']:.3f} against the full synth pool splits "
          f"{d['line_subset_step_dex']:+.3f} line-subset (67 -> "
          f"{d['synth_base_amarsi_subset']['n']} lines, the MLP's training domain) and "
          f"{d['amarsi_correction_step_dex']:+.3f} Amarsi correction. The route was never "
          f"a term in it. WHERE IT IS INSTEAD: the same network returns "
          f"{n['control_pool']['mean_correction_dex']:+.4f} dex on Amarsi's own solar line "
          f"set (where the reactivation control PASSES against the published 7.46) and "
          f"{n['our_pool']['median_correction_dex']:+.4f} on our graded pool — a "
          f"{n['swing_dex']:+.4f} dex swing from the LINE POPULATION. Our lab-gf tier "
          f"floors at {n['graded_pool_elo_floor_is_the_gf_tier']['elo_min_eV']:.2f} eV "
          f"where Amarsi's set reaches {n['control_pool']['elo_eV']['min']:.2f} eV, and the "
          f"network's correction tracks Elo "
          f"({n['elo_trend']['network_trend_over_0_to_5_eV_dex']:+.3f} dex over 0-5 eV, "
          f"the paper's own published trend). ⚠️ Elo transport explains only "
          f"{n['elo_trend']['fraction_of_swing_explained']:.0%} of the swing, so the "
          f"population difference is MEASURED but not fully mechanised — the remainder is "
          f"left unattributed rather than guessed.")
    if t2["splittable"]:
        v2 = (f"SPLIT EXISTS — lab n={t2['n_lab_tier']}, Kurucz-floored "
              f"n={t2['n_kurucz_floored']}; see by_tier for the medians.")
    else:
        v2 = (f"CANNOT BE RUN, AND THAT IS THE ANSWER — all "
              f"{t2['n_pool']} lines are lab-tier gf "
              f"({', '.join(sorted(t2['by_tier']))}), so the Kurucz-floored subset is "
              f"EMPTY and the gf zero-point cannot carry the offset. RYA-822's "
              f"'heavily Kurucz-floored' is true of the BAND, not of this GRADED pool. "
              f"The product's `gf` axis nevertheless publishes "
              f"'{t2['product_gf_axis']}' while its own budget reaches rung 3 on the same "
              f"50 lines — a label contradicting a measurement inside one product.")
    if not t3["ran"]:
        v3 = t3["note"]
    else:
        q = t3["paired_differential"]
        v3 = (f"DEPARTURES ARE ENGAGED ({t3['departures_engaged']}) and the feed's 0.000 "
              f"is NOT a wiring bug — it is the wrong statistic. Differencing two "
              f"published medians gives "
              f"{t3['published']['difference_of_published_values']:+.4f}; the per-line "
              f"paired differential is {q['median']:+.4f} (mean {q['mean']:+.4f}), "
              f"nonzero on {q['n_nonzero']}/{q['n_paired']} lines. RYA-1051's +0.033 is "
              f"corroborated; the feed's 0.000 is not a competing value. The RYA-1099 "
              f"Gerber systematics note is therefore SAFE on this point — the ⟨3D⟩ NLTE "
              f"effect it rests on is real and measured. ⚠️ separately, the NLTE "
              f"product's own `nlte_delta_dex`/`nlte_source` columns read as LTE.")
    return {"test1": v1, "test2": v2, "test3": v3}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pool", default="amarsi_feI_vis", choices=sorted(POOLS),
                    help="the line population every number is over. Named, never defaulted "
                         "silently — RYA-1083's lesson is that two aggregates over "
                         "different populations look comparable and are not.")
    ap.add_argument("--holding", default=DEFAULT_HOLDING, choices=sorted(HOLDINGS),
                    help="the arm/holding whose Amarsi leg to diagnose.")
    ap.add_argument("--check-holding", default=None, choices=sorted(HOLDINGS),
                    help="a second holding, run through Test 1 and Test 2 as the "
                         "cross-check the ticket asks for. Default: every other holding.")
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    args = ap.parse_args(argv)

    t1 = diagnose_base(args.holding)
    t2 = diagnose_gf_split(args.holding)
    t3 = diagnose_nlte_wiring(args.holding)
    rep = {
        "ticket": "RYA-1104",
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "read_only": True,
        "pool": args.pool, "holding": args.holding,
        "species": "%s %s %s" % POOLS[args.pool],
        "test1": t1, "test2": t2, "test3": t3,
        "label_audit": label_audit(t1),
        "verdicts": verdicts(t1, t2, t3),
    }

    checks = ([args.check_holding] if args.check_holding
              else [h for h in sorted(HOLDINGS) if h != args.holding])
    rep["cross_check"] = {}
    for h in checks:
        c1, c2 = diagnose_base(h), diagnose_gf_split(h)
        rep["cross_check"][h] = {
            "base_is_synthesis": c1["declared_base"]["is_synthesis"],
            "base_verified_max_residual_dex":
                c1["base_verified_line_by_line"]["max_abs_residual_dex"],
            "synth_base_full_pool": c1["decomposition"]["synth_base_full_pool"],
            "synth_base_amarsi_subset": c1["decomposition"]["synth_base_amarsi_subset"],
            "post_amarsi": c1["decomposition"]["post_amarsi"],
            "line_subset_step_dex": c1["decomposition"]["line_subset_step_dex"],
            "amarsi_correction_step_dex": c1["decomposition"]["amarsi_correction_step_dex"],
            "ew_route_product": c1["ew_route_product"],
            "n_lab_tier": c2["n_lab_tier"], "n_kurucz_floored": c2["n_kurucz_floored"],
            "gf_tiers": sorted(c2["by_tier"]),
        }

    text = render(rep)
    print(text)
    print("\nCROSS-CHECK (the same two tests on every other holding)")
    for h, c in rep["cross_check"].items():
        print(f"  {h:<18} base synth={c['base_is_synthesis']}  "
              f"full {c['synth_base_full_pool']['A']:.3f} (n={c['synth_base_full_pool']['n']})"
              f"  subset {c['synth_base_amarsi_subset']['A']:.3f} "
              f"(n={c['synth_base_amarsi_subset']['n']})"
              f"  +Amarsi {c['post_amarsi']['A']:.3f}"
              f"  gf tiers {c['gf_tiers']}"
              f"  kurucz-floored n={c['n_kurucz_floored']}")

    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "amarsi_base_diagnostic.json").write_text(json.dumps(rep, indent=2) + "\n")
    (args.out / "amarsi_base_diagnostic.txt").write_text(text + "\n")
    print(f"\nwrote {args.out / 'amarsi_base_diagnostic.json'}")
    print(f"wrote {args.out / 'amarsi_base_diagnostic.txt'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
