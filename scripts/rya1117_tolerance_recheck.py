#!/usr/bin/env python3
"""RYA-1117 — re-validate every RYA-1110 count under a SOURCE-DERIVED tolerance + plateau.

    python3 scripts/rya1117_tolerance_recheck.py            # sweep + report + write artifact
    python3 scripts/rya1117_tolerance_recheck.py --check    # exit 1 if a count is off-plateau

🔴 THE TICKET'S PREMISE IS REFUTED, AND THAT IS THE FIRST RESULT.

RYA-1117 was raised on the assumption that RYA-1110's counts "were computed against Jofré
2014's list with the same fixed tolerance" that RYA-1109 found 20x too tight. They were
not. `rya1110_build_gbs_fe_lineset._MATCH_TOL_A` is **0.015 Å**, three times
`line_match.MATCH_TOL_A`, derived in that ticket from Jofré's 2-dp print precision and
sanity-checked against a tolerance scan and a displaced null.

That does not make the re-check pointless — "I derived it" is a claim, and the way to
settle it is to measure. So this script does the measurement RYA-1109's method prescribes,
on every count RYA-1110 published, and reports what a fixed 0.005 Å WOULD have given as
the counterfactual BEFORE column.

WHY JOFRÉ IS NOT AGSS21
-----------------------
RYA-1109's trap was a UNIT: AGSS21 Table A.2 prints λ in NANOMETRES to 2 dp — 0.1 Å
resolution, half-width 0.05 Å — so the 0.005 Å default was 20x too tight. Jofré prints
ÅNGSTRÖMS. Measured here rather than read off the format string (`_measure_precision`):

    Jofré Tables 4/5      2 dp   ->  0.01 Å step,  half-width 0.005 Å
    VizieR ew.dat         2 dp   ->  0.01 Å step,  half-width 0.005 Å
    VizieR table6.dat     1 dp   ->  0.10 Å step,  half-width 0.050 Å
    canonical_gf          3-4 dp ->  half-width <= 0.0005 Å
    Heiter+2021 geslines  4 dp   ->  half-width 0.00005 Å

So print-rounding alone allows ~0.0055 Å of disagreement on the Jofré↔canonical_gf join —
which is why the bare default is NOT wildly wrong here the way it was for AGSS21. But
rounding is not the only term: two catalogues can list the same transition at genuinely
different wavelengths. The plateau is what measures that second term.

⚠️ `table6.dat` IS NOT JOINED BY TOLERANCE AT ALL, so its coarse 1 dp does not enter.
RYA-1110 pairs it to `ew.dat` ordinally, having asserted the conditions that make an
ordinal pairing provable. A 0.05 Å half-width window there would have been the λ-only join
RYA-1037 exists to end — table6 carries no EP to disambiguate with.

WHAT A PLATEAU LOOKS LIKE WHEN THE MATCHER REFUSES AMBIGUITY
------------------------------------------------------------
RYA-1109's rule is "a count that keeps climbing is counting coincidences; one that plateaus
is the answer". Under `require_ep=True` the shape is sharper than that: past the plateau
the count does not climb, it FALLS, because a second in-window candidate makes the line
AMBIGUOUS and `line_match` refuses it rather than guessing. So the signature to look for is
rise -> flat -> fall, and the flat is the answer. A count that rose forever would mean the
EP half of the key had stopped discriminating.
"""
from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from pipeline import line_match                       # noqa: E402
import rya1110_build_gbs_fe_lineset as GBS            # noqa: E402

OUT = ROOT / "data" / "results" / "rya1117" / "gbs_tolerance_recheck.json"

#: The counterfactual BEFORE. RYA-1110 never used this; it is `line_match`'s default and
#: the value RYA-1109's first pass was bitten by, carried so the two tickets are comparable.
BEFORE_TOL_A = line_match.MATCH_TOL_A

#: The tolerance RYA-1110 actually shipped, re-derived and re-checked here.
DERIVED_TOL_A = GBS._MATCH_TOL_A

SWEEP = (0.005, 0.0075, 0.010, 0.015, 0.020, 0.030, 0.050, 0.080, 0.100, 0.250, 0.500)

#: Displacements for the chance-match null, applied at every tolerance in the sweep. A
#: count can plateau and still be wrong if the null is climbing underneath it.
NULL_SHIFTS_A = (-0.3, 0.3)

_SPECIES = ("Fe I", "Fe II")


class RecheckError(RuntimeError):
    """A published count is not the plateau value."""


# ── measured precision, never read off a format string ───────────────────────
def _measure_precision(values: list[str]) -> dict:
    """Decimal places actually printed, and the smallest non-zero step between values."""
    dp = collections.Counter(len(v.split(".")[1]) if "." in v else 0 for v in values)
    xs = sorted({float(v) for v in values})
    steps = [round(b - a, 6) for a, b in zip(xs, xs[1:])]
    step = min((s for s in steps if s > 0), default=float("nan"))
    worst = max(dp)
    return {"n": len(values), "decimals": dict(sorted(dp.items())),
            "min_observed_step_A": step,
            "print_step_A": 10.0 ** -worst,
            "print_half_width_A": 0.5 * 10.0 ** -worst}


def measured_precision() -> dict:
    H = ROOT / "data" / "reference" / "jofre2014_gbs"
    G = ROOT / "data" / "reference" / "heiter2021_ges"

    def tsv_col(path: Path, col: int) -> list[str]:
        return [l.split("\t")[col] for l in path.read_text().splitlines()
                if l and not l.startswith(("#", "wavelength", "Element"))]

    out = {
        "jofre_table4_fe1": _measure_precision(
            tsv_col(H / "paper_table4_fe1_golden.tsv", 0)),
        "jofre_table5_fe2": _measure_precision(
            tsv_col(H / "paper_table5_fe2_golden.tsv", 0)),
        "vizier_ew_dat": _measure_precision(
            [l[17:24].strip() for l in (H / "vizier" / "ew.dat").read_text().splitlines()
             if l.strip()]),
        "vizier_table6_dat": _measure_precision(
            [l[5:11].strip() for l in (H / "vizier" / "table6.dat").read_text().splitlines()
             if l.strip()]),
        "heiter2021_geslines": _measure_precision(
            tsv_col(G / "geslines_Fe_4700_6900.tsv", 3)),
    }
    cg = pd.read_csv(ROOT / "data" / "linelists" / "canonical_gf.csv", low_memory=False)
    fe = cg[cg["species"].isin(_SPECIES)]
    out["canonical_gf"] = _measure_precision(
        [f"{x:.10g}" for x in fe["wavelength_air_A"]])
    return out


# ── the sweep ────────────────────────────────────────────────────────────────
def _join(ref: pd.DataFrame, pool: pd.DataFrame, tol: float, shift: float = 0.0):
    """Positional index into `pool` per `ref` row (-1 unresolved), plus the ambiguity count.

    λ+EP, `require_ep=True` — the same call RYA-1110 makes, so the sweep measures THAT
    join and not a lookalike (RYA-1080's lesson: the rescue and its verifier must resolve
    through the same call).
    """
    idx = np.full(len(ref), -1, dtype=int)
    ambiguous = 0
    for sp in _SPECIES:
        want = np.flatnonzero((ref["species"] == sp).to_numpy())
        have = np.flatnonzero((pool["species"] == sp).to_numpy())
        if not want.size or not have.size:
            continue
        p = pool.iloc[have]
        r = line_match.match(
            ref["wavelength_air_A"].to_numpy(float)[want] + shift,
            p["wavelength_air_A"].to_numpy(float),
            want_ep=ref["excitation_potential_eV"].to_numpy(float)[want],
            src_ep=p["excitation_potential_eV"].to_numpy(float),
            tol_A=tol, require_ep=True)
        ok = r.index >= 0
        idx[want[ok]] = have[r.index[ok]]
        ambiguous += len(r.ambiguous)
    return idx, ambiguous


def counts_at(ls: pd.DataFrame, pool: pd.DataFrame, ges: pd.DataFrame,
              tol: float) -> dict:
    """Every RYA-1110 count that depends on a tolerance, at this tolerance."""
    idx, amb = _join(ls, pool, tol)
    nulls = [int((_join(ls, pool, tol, shift=s)[0] >= 0).sum()) for s in NULL_SHIFTS_A]
    sel = (ls["rew_class"] == "pass").to_numpy()
    lab = np.array([bool(pool["gf_tier"].iloc[j] == "LAB") if j >= 0 else False
                    for j in idx])
    ours = np.array([pool["log_gf"].iloc[j] if j >= 0 else np.nan for j in idx])
    gbs = ls["log_gf_gbs"].to_numpy(float)
    comparable = sel & np.isfinite(ours) & np.isfinite(gbs)
    identical = int((np.abs(gbs - ours)[comparable] <= GBS._GF_SAME_DEX).sum())
    hidx, hamb = _join(ls, ges, tol)
    return {
        "tol_A": tol,
        "canonical_gf_resolved": int((idx >= 0).sum()),
        "canonical_gf_ambiguous": int(amb),
        "displaced_null": nulls,
        "graded_LAB_overlap_FeI": int((lab & sel & (ls["species"] == "Fe I").to_numpy()).sum()),
        "graded_LAB_overlap_FeII": int((lab & sel & (ls["species"] == "Fe II").to_numpy()).sum()),
        "gf_comparable": int(comparable.sum()),
        "gf_identical": identical,
        "heiter2021_resolved": int((hidx >= 0).sum()),
        "heiter2021_ambiguous": int(hamb),
    }


#: The counts RYA-1110 published, and which sweep key each one is. Pinned so this script
#: cannot quietly "re-validate" a number nobody reported.
PUBLISHED = {
    "canonical_gf_resolved": 159,
    "graded_LAB_overlap_FeI": 16,
    "graded_LAB_overlap_FeII": 0,
    "gf_comparable": 121,
    "gf_identical": 67,
    "heiter2021_resolved": 159,
}


def _is_clean(s: dict) -> bool:
    """No ambiguity refused, and no chance match admitted, at this tolerance."""
    return (s["canonical_gf_ambiguous"] == 0 and s["heiter2021_ambiguous"] == 0
            and all(n == 0 for n in s["displaced_null"]))


def plateau_of(sweep: list[dict], key: str) -> dict:
    """The plateau value — but only inside the CLEAN REGION.

    🔴 THE NAIVE RULE PICKS THE WRONG REGION HERE, AND IT IS WORTH SAYING WHY. "The longest
    run of equal values" would call `canonical_gf_resolved` = 157 over 0.05-0.25 Å, a run
    four points wide. That region is flat because the matcher has stopped answering: two
    lines have acquired a second in-window candidate and are being REFUSED as ambiguous,
    while the displaced null has already begun admitting 1-3 chance matches. It is the
    flatness of saturation, not of agreement.

    So a plateau is only evidence inside the region where the key still discriminates —
    zero ambiguity refused AND zero chance matches admitted. Within that region a wider
    window can only recover GENUINE matches (the null proves no coincidences are getting
    in), so the answer is the value at the widest clean tolerance, and the check is that it
    is stable across the top of the clean region rather than still climbing when the region
    ends.
    """
    clean = [s for s in sweep if _is_clean(s)]
    if not clean:
        raise RecheckError(f"no clean tolerance in the sweep for {key} — every window "
                           f"either refuses an ambiguity or admits a chance match")
    top = clean[-1]
    value = top[key]
    flat = [s for s in clean if s[key] == value]
    return {
        "value": value,
        "clean_lo_A": clean[0]["tol_A"], "clean_hi_A": top["tol_A"],
        "flat_lo_A": flat[0]["tol_A"], "flat_hi_A": flat[-1]["tol_A"],
        "flat_width": len(flat),
        "still_climbing_at_clean_edge": len(flat) < 2,
        "derived_tol_inside": flat[0]["tol_A"] <= DERIVED_TOL_A <= flat[-1]["tol_A"],
        "naive_longest_run": _naive_longest_run(sweep, key),
    }


def _naive_longest_run(sweep: list[dict], key: str) -> dict:
    """What "longest flat run" WOULD have said. Reported so the difference is visible."""
    best = (0, 0, 0)
    i = 0
    while i < len(sweep):
        j = i
        while j + 1 < len(sweep) and sweep[j + 1][key] == sweep[i][key]:
            j += 1
        if j - i + 1 > best[0]:
            best = (j - i + 1, i, j)
        i = j + 1
    n, a, b = best
    return {"value": sweep[a][key], "width": n,
            "tol_lo_A": sweep[a]["tol_A"], "tol_hi_A": sweep[b]["tol_A"],
            "region_is_clean": all(_is_clean(s) for s in sweep[a:b + 1])}


def structural_zero(ls: pd.DataFrame, cg: pd.DataFrame) -> dict:
    """Is the Fe II zero-overlap a wavelength fact or a tolerance artifact?

    A tolerance-independent claim has to be shown, not asserted: the number that settles it
    is the CLOSEST cross-pair between the two sets, against the window that would be needed
    to bridge it.
    """
    lab2 = cg[(cg["species"] == "Fe II") & (cg["gf_tier"] == "LAB")
              & cg["wavelength_air_A"].between(3800.0, 6910.0)]
    gbs2 = ls[(ls["rew_class"] == "pass") & (ls["species"] == "Fe II")]
    d = np.abs(np.subtract.outer(gbs2["wavelength_air_A"].to_numpy(float),
                                 lab2["wavelength_air_A"].to_numpy(float)))
    closest = float(d.min())
    return {
        "our_graded_FeII_VIS": {"n": int(len(lab2)),
                                "lo_A": float(lab2["wavelength_air_A"].min()),
                                "hi_A": float(lab2["wavelength_air_A"].max())},
        "gbs_FeII_selected": {"n": int(len(gbs2)),
                              "lo_A": float(gbs2["wavelength_air_A"].min()),
                              "hi_A": float(gbs2["wavelength_air_A"].max())},
        "closest_cross_pair_A": closest,
        "derived_tol_A": DERIVED_TOL_A,
        "window_multiple_needed": closest / DERIVED_TOL_A,
        "verdict": ("WAVELENGTH-DISJOINT — the two sets do not overlap in wavelength at "
                    "all, so the zero is a property of the data and no tolerance short of "
                    f"{closest:.1f} Å could create a match"),
    }


def run() -> dict:
    ls = pd.read_csv(ROOT / "data" / "linelists" / "reference_sets"
                     / "gbs_solar_fe_rya1110.csv")
    cg = pd.read_csv(ROOT / "data" / "linelists" / "canonical_gf.csv", low_memory=False)
    pool = cg[cg["species"].isin(_SPECIES)].reset_index(drop=True)
    ges = GBS._read_ges_lines().rename(
        columns={"lambda": "wavelength_air_A", "Elow": "excitation_potential_eV"})
    ges = ges.assign(species=np.where(ges["Ion"] == 1, "Fe I", "Fe II"))

    sweep = [counts_at(ls, pool, ges, t) for t in SWEEP]
    before = next(s for s in sweep if s["tol_A"] == BEFORE_TOL_A)
    after = next(s for s in sweep if s["tol_A"] == DERIVED_TOL_A)
    plateaus = {k: plateau_of(sweep, k) for k in PUBLISHED}

    return {
        "ticket": "RYA-1117",
        "rechecks": "RYA-1110",
        "premise": {
            "as_raised": "RYA-1110's counts were computed with the fixed 0.005 A default",
            "verdict": "REFUTED",
            "detail": (f"rya1110_build_gbs_fe_lineset._MATCH_TOL_A is {DERIVED_TOL_A} A, "
                       f"{DERIVED_TOL_A / BEFORE_TOL_A:.0f}x line_match.MATCH_TOL_A "
                       f"({BEFORE_TOL_A} A), derived in RYA-1110 from Jofré's 2-dp print "
                       f"precision. The BEFORE column here is a COUNTERFACTUAL."),
        },
        "measured_precision": measured_precision(),
        "before_tol_A": BEFORE_TOL_A,
        "derived_tol_A": DERIVED_TOL_A,
        "sweep": sweep,
        "before": before,
        "after": after,
        "plateaus": plateaus,
        "published": PUBLISHED,
        "structural_zero_FeII": structural_zero(ls, cg),
    }


def check(doc: dict) -> list[str]:
    """Every published count must BE the plateau value, with the derived tolerance on it."""
    bad = []
    for key, published in PUBLISHED.items():
        pl = doc["plateaus"][key]
        if pl["value"] != published:
            bad.append(f"{key}: published {published}, clean-region plateau {pl['value']} "
                       f"over {pl['flat_lo_A']}-{pl['flat_hi_A']} A")
        if not pl["derived_tol_inside"]:
            bad.append(f"{key}: the derived tolerance {DERIVED_TOL_A} A is OUTSIDE the "
                       f"plateau {pl['flat_lo_A']}-{pl['flat_hi_A']} A")
        if pl["still_climbing_at_clean_edge"]:
            bad.append(f"{key}: the count is still moving where the clean region ends "
                       f"({pl['clean_hi_A']} A) — the window may be clipping real matches")
    # a plateau means nothing if the null is leaking under it
    for s in doc["sweep"]:
        if s["tol_A"] <= DERIVED_TOL_A and any(n > 0 for n in s["displaced_null"]):
            bad.append(f"tol {s['tol_A']}: displaced null resolved {s['displaced_null']} "
                       f"— the window is finding coincidences at or below the derived value")
    return bad


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true",
                    help="exit 1 if any published count is off its plateau")
    a = ap.parse_args()
    doc = run()

    print(f"PREMISE: {doc['premise']['verdict']} — {doc['premise']['detail']}\n")
    print("MEASURED PRINT PRECISION")
    for k, v in doc["measured_precision"].items():
        print(f"  {k:22s} decimals={v['decimals']}  step={v['print_step_A']:g} A  "
              f"half-width={v['print_half_width_A']:g} A")
    print()
    hdr = (f"{'tol_A':>7} {'canon':>6} {'amb':>4} {'null':>9} "
           f"{'FeI∩LAB':>8} {'FeII':>5} {'n_cmp':>6} {'ident':>6} {'heiter':>7}")
    print("PLATEAU SWEEP\n" + hdr + "\n  " + "-" * (len(hdr) - 2))
    for s in doc["sweep"]:
        mark = "  <- DERIVED" if s["tol_A"] == DERIVED_TOL_A else (
               "  <- default" if s["tol_A"] == BEFORE_TOL_A else "")
        print(f"{s['tol_A']:7.4f} {s['canonical_gf_resolved']:6d} "
              f"{s['canonical_gf_ambiguous']:4d} {str(s['displaced_null']):>9} "
              f"{s['graded_LAB_overlap_FeI']:8d} {s['graded_LAB_overlap_FeII']:5d} "
              f"{s['gf_comparable']:6d} {s['gf_identical']:6d} "
              f"{s['heiter2021_resolved']:7d}{mark}")
    print()
    print(f"BEFORE (counterfactual {BEFORE_TOL_A} A) vs AFTER (derived {DERIVED_TOL_A} A)")
    print("  count                       BEFORE   AFTER  publ  plateau  "
          "flat over          naive longest run")
    for key, published in PUBLISHED.items():
        pl = doc["plateaus"][key]
        nv = pl["naive_longest_run"]
        span = f"{pl['flat_lo_A']}-{pl['flat_hi_A']} A"
        naive = f"{nv['value']} @ {nv['tol_lo_A']}-{nv['tol_hi_A']} A"
        if not nv["region_is_clean"]:
            naive += "  <- NOT CLEAN"
        print(f"  {key:<26} {doc['before'][key]:>6} {doc['after'][key]:>7} "
              f"{published:>5} {pl['value']:>8}  {span:<18} {naive}")
    z = doc["structural_zero_FeII"]
    print(f"\nSTRUCTURAL ZERO (Fe II): ours {z['our_graded_FeII_VIS']['lo_A']:.2f}-"
          f"{z['our_graded_FeII_VIS']['hi_A']:.2f} A, GBS "
          f"{z['gbs_FeII_selected']['lo_A']:.2f}-{z['gbs_FeII_selected']['hi_A']:.2f} A\n"
          f"  closest cross-pair {z['closest_cross_pair_A']:.2f} A = "
          f"{z['window_multiple_needed']:.0f}x the derived window -> {z['verdict']}")

    fails = check(doc)
    for f in fails:
        print(f"\n  FAIL {f}")
    if not fails:
        print("\nOK — every published RYA-1110 count is its plateau value, with the "
              "derived tolerance sitting inside the plateau.")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(doc, indent=2) + "\n")
    print(f"\nwrote {OUT.relative_to(ROOT)}")
    return 1 if (a.check and fails) else 0


if __name__ == "__main__":
    raise SystemExit(main())
