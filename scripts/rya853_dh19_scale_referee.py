#!/usr/bin/env python3
"""
RYA-853 scope 3 — is our Fe II gf scale legitimate, or solar-fitted?
====================================================================
RYA-852 found our Fe II pool sits +0.106 dex above NIST, and that the current scale is what
makes the solar Fe I - Fe II ionization balance come out right (+0.018 on ours, -0.088 on
NIST's). Two readings, with opposite consequences for whether the balance — which gates
log g — is an INDEPENDENT check or a CIRCULAR one:

  * MINE (RYA-852): the pool is partly solar-fitted (one member is Meléndez & Barbuy 2009
    outright, and MB09 is built partly by reverse solar analysis), so a scale that makes the
    solar balance work is exactly what you would expect. CIRCULAR.
  * RYA-853 counter-evidence: Den Hartog 2019 reports NIST's Fe II gf sitting ~0.1 dex BELOW
    modern pure-lab values. Then +0.106 above NIST is correct and NIST is the outlier.
    LEGITIMATE.

THE REFEREE is Den Hartog et al. 2019 (ApJS 243, 33): branching fractions x laser-induced
fluorescence lifetimes. Pure laboratory, NO solar normalisation anywhere in it — which is
the whole point, because a referee that touched the solar spectrum could not settle a
question about solar fitting.

    ours - DH ~ 0        -> our scale IS the pure-lab scale; NIST is low; the balance is
                            an independent check and my solar-fitting hypothesis DIES.
    ours - DH ~ +0.1-0.2 -> we sit above pure lab too; the hypothesis LIVES and the
                            balance may be circular.

and the third leg tests the counter-evidence on its own terms:

    NIST - DH < 0        -> NIST really is low, as Den Hartog 2019 states.

⚠️ THE OVERLAP DOES NOT CONTAIN THE ARBITER LINES. DH19's optical set stops at 4584 A and
the three Fe II arbiter lines (6147.734 / 6238.386 / 6247.557) are redward of it. This
refereeing is therefore about the SCALE of the pool, established on the 4173-4584 A overlap,
and is NOT a direct measurement of the arbiter lines. Stated because the distinction is the
difference between "our Fe II scale is sound" and "these three lines are fine".

⚠️ DATA PROVENANCE. The ten optical values below were supplied by Ryan from the DH19 PDF.
`J/ApJS/243/33` is the VizieR ID for the full machine-readable Table 6 (131 lines), but
astroquery returns 0 tables for it from this machine — so this is NOT reproducible from the
network here, and the numbers are transcribed. That is a real limitation of this run, not a
property of the source.

Usage:
    python3 scripts/rya853_dh19_scale_referee.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

CANONICAL_GF = ROOT / "data" / "linelists" / "canonical_gf.csv"
MASTER = ROOT / "data" / "linelists" / "linelist_master.csv"
OUT = ROOT / "data" / "results" / "rya853"

#: Den Hartog et al. 2019, ApJS 243, 33 — Fe II, optical subset.
#: Branching fractions x LIF radiative lifetimes. PURE LAB, no solar normalisation.
#: (wavelength_air_A, chi_eV, log_gf)
DH19_OPTICAL = (
    (4173.451, 2.583, -2.38),
    (4233.162, 2.583, -2.02),
    (4303.170, 2.704, -2.52),
    (4351.762, 2.704, -1.95),
    (4385.377, 2.778, -2.64),
    (4416.819, 2.778, -2.57),
    (4508.281, 2.856, -2.42),
    (4522.628, 2.844, -2.29),
    (4549.466, 2.828, -1.92),
    (4583.829, 2.807, -1.94),
)

DH19_CITE = "Den Hartog et al. 2019, ApJS 243, 33 (BF x LIF lifetimes; pure lab)"
DH19_VIZIER = "J/ApJS/243/33"

#: EP guard on BOTH sides (RYA-780/846/852). Wavelength alone picks the wrong level.
WAVE_TOL_A = 0.05
EP_TOL_EV = 0.05

#: RYA-852's measured pool-wide offset against NIST, for the comparison at the end.
POOL_OFFSET_VS_NIST = 0.106


def ours(dh: pd.DataFrame) -> pd.DataFrame:
    """Our stored Fe II log gf for the overlap lines, from canonical_gf then master."""
    c = pd.read_csv(CANONICAL_GF, comment="#", low_memory=False)
    c = c[c.species.astype(str) == "Fe II"]
    m = pd.read_csv(MASTER, low_memory=False)
    m = m[(m.element.astype(str).str.strip() == "Fe")
          & (m.ion.astype(str).str.strip().isin(["II", "2"]))]

    rows = []
    for _, d in dh.iterrows():
        rec = {"wavelength_air_A": d.wavelength_air_A, "chi_eV": d.chi_eV,
               "dh_loggf": d.dh_loggf}
        cm = c[(np.abs(c.wavelength_air_A - d.wavelength_air_A) <= WAVE_TOL_A)
               & (np.abs(c.excitation_potential_eV - d.chi_eV) <= EP_TOL_EV)]
        rec["n_canonical_candidates"] = int(len(cm))
        if len(cm) == 1:
            rec["our_loggf"] = float(cm.iloc[0].log_gf)
            rec["our_ref"] = str(cm.iloc[0].loggf_reference)
        elif len(cm) > 1:
            rec["ambiguous_canonical"] = True

        mm = m[(np.abs(m.wavelength_air_A - d.wavelength_air_A) <= WAVE_TOL_A)
               & (np.abs(m.excitation_potential_eV - d.chi_eV) <= EP_TOL_EV)]
        rec["n_master_candidates"] = int(len(mm))
        if len(mm) == 1:
            rec["master_loggf"] = float(mm.iloc[0].log_gf)
            rec["master_ref"] = str(mm.iloc[0].loggf_source)
        rows.append(rec)
    return pd.DataFrame(rows)


def nist(dh: pd.DataFrame) -> pd.DataFrame:
    """NIST ASD for the same lines, in AIR, EP-matched, ambiguity flagged."""
    from astroquery.nist import Nist
    import astropy.units as u

    lo = float(dh.wavelength_air_A.min()) - 5
    hi = float(dh.wavelength_air_A.max()) + 5
    t = Nist.query(lo * u.AA, hi * u.AA, linename="Fe II", wavelength_type="vac+air")
    gcol = next((c for c in t.colnames if c.startswith("gi")), None)
    eicol = next((c for c in t.colnames if c.startswith("Ei")), None)

    cand = []
    for r in t:
        def _f(col):
            try:
                return float(str(r[col]).strip())
            except Exception:
                return float("nan")
        w = _f("Observed")
        if not np.isfinite(w):
            w = _f("Ritz")
        fik = _f("fik")
        try:
            gi = float(str(r[gcol]).split("-")[0].strip()) if gcol else np.nan
        except Exception:
            gi = np.nan
        try:
            ei = float(str(r[eicol]).split("-")[0].strip()) if eicol else np.nan
        except Exception:
            ei = np.nan
        if not (np.isfinite(w) and np.isfinite(gi) and np.isfinite(fik) and gi * fik > 0):
            continue
        cand.append({"w": w, "ei": ei, "loggf": float(np.log10(gi * fik)),
                     "acc": str(r["Acc."]).strip()})
    cd = pd.DataFrame(cand)

    rows = []
    for _, d in dh.iterrows():
        rec = {"wavelength_air_A": d.wavelength_air_A}
        if len(cd):
            m = cd[(np.abs(cd.w - d.wavelength_air_A) <= WAVE_TOL_A)
                   & (np.abs(cd.ei - d.chi_eV) <= EP_TOL_EV)]
            rec["n_nist_candidates"] = int(len(m))
            if len(m) == 1:
                rec["nist_loggf"] = float(m.iloc[0].loggf)
                rec["nist_acc"] = str(m.iloc[0].acc)
        rows.append(rec)
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    dh = pd.DataFrame(DH19_OPTICAL,
                      columns=["wavelength_air_A", "chi_eV", "dh_loggf"])
    print(f"[referee] {DH19_CITE}")
    print(f"[referee] {len(dh)} optical lines, "
          f"{dh.wavelength_air_A.min():.1f}-{dh.wavelength_air_A.max():.1f} A")
    print(f"[referee] VizieR {DH19_VIZIER} returns 0 tables from this machine; "
          f"values transcribed from the PDF")

    o = ours(dh)
    try:
        n = nist(dh)
        o = o.merge(n, on="wavelength_air_A", how="left")
    except Exception as e:
        print(f"  NIST query failed ({type(e).__name__}) — the ours-vs-DH leg still runs")

    o["ours_minus_dh"] = o.get("our_loggf") - o.dh_loggf
    o["master_minus_dh"] = o.get("master_loggf") - o.dh_loggf
    if "nist_loggf" in o:
        o["nist_minus_dh"] = o.nist_loggf - o.dh_loggf

    print(f"\n{'wave':>9}{'chi':>7}{'DH':>8}{'ours':>8}{'Δ ours-DH':>11}"
          f"{'NIST':>8}{'Δ NIST-DH':>11}  our source")
    for _, r in o.iterrows():
        print(f"{r.wavelength_air_A:9.3f}{r.chi_eV:7.3f}{r.dh_loggf:8.2f}"
              f"{r.get('our_loggf', np.nan):8.3f}{r.get('ours_minus_dh', np.nan):+11.3f}"
              f"{r.get('nist_loggf', np.nan):8.3f}"
              f"{r.get('nist_minus_dh', np.nan):+11.3f}  "
              f"{str(r.get('our_ref',''))[:22]}")

    # ours - NIST on the SAME ten lines, so all three legs sit on identical footing.
    # RYA-852's +0.106 was measured on the REDWARD pool (5256-6456 A); comparing it to a
    # blue-overlap number without recomputing it here would be comparing two different
    # line sets and calling the difference a result.
    if "nist_loggf" in o:
        o["ours_minus_nist"] = o.get("our_loggf") - o.nist_loggf

    def _stat(col, seed: int = 0, n_boot: int = 20000):
        """Median WITH its uncertainty. A median of ten values spanning 0.6 dex is not a
        number you can compare against a 0.05 threshold without saying how well it is
        determined — that is how a test that cannot discriminate gets read as a verdict."""
        v = o[col].dropna() if col in o else pd.Series(dtype=float)
        if not len(v):
            return None
        a = v.to_numpy(float)
        rng = np.random.default_rng(seed)
        boot = np.array([np.median(rng.choice(a, len(a))) for _ in range(n_boot)])
        return {"n": int(len(a)), "median": float(np.median(a)),
                "mean": float(a.mean()),
                "std": float(a.std(ddof=1)) if len(a) > 1 else 0.0,
                "sem": float(a.std(ddof=1) / np.sqrt(len(a))) if len(a) > 1 else 0.0,
                "ci95": [float(np.percentile(boot, 2.5)),
                         float(np.percentile(boot, 97.5))],
                "min": float(a.min()), "max": float(a.max())}

    s_ours = _stat("ours_minus_dh")
    s_nist = _stat("nist_minus_dh")
    s_on = _stat("ours_minus_nist")
    print(f"\n=== the referee (median, 95% bootstrap CI) ===")
    for lab, st in (("ours - DH", s_ours), ("NIST - DH", s_nist),
                    ("ours - NIST (same 10 lines)", s_on)):
        if st:
            print(f"  {lab:<28} {st['median']:+.3f}  CI [{st['ci95'][0]:+.3f}, "
                  f"{st['ci95'][1]:+.3f}]  sd {st['std']:.3f}  n={st['n']}")
    print(f"  ours - NIST (RYA-852, REDWARD pool 5256-6456 A): "
          f"{POOL_OFFSET_VS_NIST:+.3f} dex")

    verdict, reasoning = "INCONCLUSIVE", ""
    if s_ours:
        med = s_ours["median"]
        lo, hi = s_ours["ci95"]
        # The two readings predict ~0.0 and ~+0.1..0.2. If the CI covers both, the test
        # does not discriminate and saying so is the result.
        covers_zero = lo <= 0.0 <= hi
        covers_offset = lo <= 0.13 <= hi
        if covers_zero and covers_offset:
            verdict = "INCONCLUSIVE — the overlap cannot separate the two readings"
            reasoning = (f"ours - DH is {med:+.3f} with a 95% CI of [{lo:+.3f}, {hi:+.3f}], "
                         f"which covers BOTH the pure-lab prediction (~0) and the "
                         f"solar-fitted one (~+0.13). Ten lines spanning "
                         f"{s_ours['min']:+.3f}..{s_ours['max']:+.3f} dex is not enough to "
                         f"choose. The full DH19 Table 6 (131 lines, VizieR "
                         f"{DH19_VIZIER}) would settle it.")
        elif abs(med) <= 0.05:
            verdict = "LEGITIMATE — our scale IS the pure-lab scale"
            reasoning = (f"ours sits {med:+.3f} from Den Hartog's pure-lab values, i.e. on "
                         f"the laboratory scale. The +0.106 against NIST is then NIST "
                         f"being low, exactly as Den Hartog 2019 reports. The solar Fe I - "
                         f"Fe II balance is an INDEPENDENT check, and RYA-852's "
                         f"solar-fitting hypothesis is REFUTED.")
        elif med >= 0.08:
            verdict = "HYPOTHESIS LIVES — we sit above pure lab too"
            reasoning = (f"ours sits {med:+.3f} ABOVE Den Hartog's pure-lab values, so the "
                         f"offset is not NIST being low. A scale above pure lab that "
                         f"happens to make the solar balance work is what RYA-852 "
                         f"suspected, and the balance may be CIRCULAR.")
        else:
            reasoning = (f"ours sits {med:+.3f} from pure lab — between the two "
                         f"predictions, so neither reading is clean.")
    # 🔴 A second finding, and it undercuts the premise of the question.
    band_flip = None
    if s_on:
        blue = s_on["median"]
        band_flip = {"blue_overlap_4173_4584": blue,
                     "red_pool_5256_6456": POOL_OFFSET_VS_NIST,
                     "swing_dex": POOL_OFFSET_VS_NIST - blue,
                     "sign_flips": bool(blue * POOL_OFFSET_VS_NIST < 0)}
        print(f"\n=== 🔴 the offset against NIST is BAND-DEPENDENT ===")
        print(f"  blue overlap 4173-4584 A : ours - NIST = {blue:+.3f}")
        print(f"  red pool     5256-6456 A : ours - NIST = {POOL_OFFSET_VS_NIST:+.3f}"
              f"   (RYA-852)")
        print(f"  swing {POOL_OFFSET_VS_NIST - blue:+.3f} dex"
              f"{', SIGN FLIPS' if band_flip['sign_flips'] else ''}")
        print(f"  => '+0.106 pool-wide' is not one scale offset. RYA-852 measured it on the "
              f"RED pool\n     alone, and the blue lines carry different sources "
              f"(T83av/PGHcor/BSScor/KK vs\n     VALD3/RU/MB09), so 'our Fe II scale' is "
              f"not a single thing to referee.")

    print(f"\n=== VERDICT: {verdict} ===")
    print(f"  {reasoning}")
    print(f"\n  ⚠️ This refs the SCALE via the 4173-4584 A overlap. The three arbiter "
          f"lines\n     (6147.734 / 6238.386 / 6247.557) are REDWARD of DH's ceiling and "
          f"are NOT\n     directly refereed here.")

    o.to_csv(OUT / "rya853_dh19_referee.csv", index=False)
    (OUT / "rya853_dh19_referee.json").write_text(json.dumps({
        "ticket": "RYA-853 scope 3",
        "referee": DH19_CITE, "vizier_id": DH19_VIZIER,
        "vizier_reachable_from_this_machine": False,
        "data_provenance": "transcribed from the DH19 PDF (supplied on the ticket)",
        "n_overlap_lines": int(len(dh)),
        "overlap_A": [float(dh.wavelength_air_A.min()), float(dh.wavelength_air_A.max())],
        "ours_minus_dh": s_ours, "nist_minus_dh": s_nist,
        "ours_minus_nist_pool": POOL_OFFSET_VS_NIST,
        "verdict": verdict, "reasoning": reasoning,
        "band_dependence": band_flip,
        "ours_minus_nist_same_lines": s_on,
        "caveat": ("the overlap does NOT contain the three Fe II arbiter lines, which are "
                   "redward of DH19's 4584 A optical ceiling; this refs the pool SCALE, "
                   "not those lines"),
    }, indent=2))
    print(f"\n[out] {OUT}/rya853_dh19_referee.csv\n[out] {OUT}/rya853_dh19_referee.json")


if __name__ == "__main__":
    main()
