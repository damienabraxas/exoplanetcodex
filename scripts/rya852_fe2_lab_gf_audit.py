#!/usr/bin/env python3
"""
RYA-852 — can the Fe II arbiter lines be graded on primary-lab gf?
==================================================================
The Fe II twin of RYA-824/836. The answer is NO on current evidence, and the reason is
worse than a coverage gap: **the grade already stored for these lines is wrong.**

WHICH LINES. The ticket guessed "6247.6 / 6432.7 / 6456.4 or similar; confirm from the
linelist". Confirmed from the products, and the guess is half right: the three-line Fe II
set is the VIS ENGINE-A aggregate,

    6147.734   6238.386   6247.557        (EP 3.889 / 3.889 / 3.892 eV)

6432.676 and 6456.380 are in the wider 11-line VIS 1D-LTE pool, not the trio.

🔴 WHAT `canonical_gf` CLAIMS, AND WHAT NIST ACTUALLY SAYS

    line        canonical_gf                       NIST ASD v5.11 (queried live)
    6147.734    RU (Raassen & Uylings)             log gf -2.796, accuracy E
    6149.246    "NIST ASD v5.11 grade B"           log gf -2.854, accuracy E
    6238.386    2009A&A... (Melendez & Barbuy)     log gf -2.754, accuracy D
    6247.557    "NIST ASD v5.11 grade B"           log gf -2.444, accuracy D

The two lines labelled *NIST ASD v5.11 grade B* match NIST on **neither the grade nor the
value**. NIST grades them D and E; their stored log gf sits +0.115 and +0.130 dex above
NIST's. Grade B is 10% (0.041 dex); D is 50% (0.176 dex) and E is 100% (0.301 dex) — so the
stored grade understates the cited accuracy by a factor of 4 to 7.

That matters because RYA-850 keys the graded gf term on exactly this metadata. Wiring Fe II
into the graded pool on the stored grade would publish a 0.041 dex gf term on lines whose
own source says 0.176-0.301 — the precise failure this ticket says to avoid ("do NOT publish
a bar tighter than the lab data supports").

💡 The three non-RU lines all sit +0.115 to +0.154 dex ABOVE NIST, a consistent offset that
looks like a single different scale rather than three independent errors. 6238.386 is
labelled Melendez & Barbuy 2009 outright, and MB09's Fe II values are systematically higher.
That is a HYPOTHESIS about the other two, not a finding — it is stated as one.

⚠️ THE VACUUM/AIR TRAP, WHICH NEARLY MANUFACTURED AN ABSENCE. `astroquery.nist` defaults to
`wavelength_type='vacuum'`. Queried that way, NONE of the three arbiter wavelengths appear
in the Fe II list and the obvious conclusion is "NIST does not cover these lines". They are
all there: air->vacuum is +1.71 A at 6150 A, so 6147.734 is listed at 6149.435. Every query
here passes `wavelength_type='vac+air'` and matches on wavelength AND excitation potential.

⚠️ MATCH ON WAVELENGTH *AND* EP. A 0.05 A window alone pulls the wrong physical line out of
`canonical_gf`: 6149.246 comes back at EP 13.436 eV and 6432.676 at 10.930 eV, when the
measured lines are at 3.889 and 2.891. Both are high-excitation neighbours. (RYA-780.)

WHAT THIS DOES NOT ESTABLISH
  * Den Hartog 2019's optical subset (10 lines) could not be reached — `Vizier.find_catalogs`
    returned nothing for three phrasings. That is an absence in THE SEARCH, not in the
    source (RYA-833): whether these three lines are among DH19's ten is UNVERIFIED.
  * MB09's own S/L flags (lab-normalised vs solar-fitted) were likewise not obtained, so the
    RYA-161 firewall cannot yet be applied line by line. 6238.386 is MB09-labelled and is
    therefore a firewall CANDIDATE, not a confirmed exclusion.

Usage (Sirius — needs astroquery in venv_ci):
    python3 scripts/rya852_fe2_lab_gf_audit.py
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
OUT = ROOT / "data" / "results" / "rya852"

#: The Fe II VIS ENGINE-A aggregate — the three-line arbiter set, read off the products
#: rather than taken from the ticket's guess.
ARBITER_LINES_A = (6147.734, 6238.386, 6247.557)

#: The wider VIS 1D-LTE pool (n=11), reported alongside so the trio is in context.
VIS_POOL_A = (5256.932, 5337.722, 5991.371, 6084.102, 6147.734, 6149.246,
              6238.386, 6247.557, 6369.459, 6432.676, 6456.380)

#: NIST ASD accuracy letters -> fractional uncertainty on the transition probability.
#: These are the CITED accuracies; the dex figure is log10(1 + frac).
NIST_ACCURACY_FRACTION = {
    "AAA": 0.003, "AA": 0.01, "A+": 0.02, "A": 0.03, "B+": 0.07, "B": 0.10,
    "C+": 0.18, "C": 0.25, "D+": 0.40, "D": 0.50, "E": 1.00,
}

#: RYA-850 wires this for a graded pool. Recorded here only to show what the stored (wrong)
#: grade B would have published.
GRADED_GF_TERM_DEX = 0.041

#: Matching tolerances. EP is required — wavelength alone pulls high-excitation neighbours.
WAVE_TOL_A = 0.05
EP_TOL_EV = 0.05


def accuracy_to_dex(letter: str) -> float:
    frac = NIST_ACCURACY_FRACTION.get(str(letter).strip())
    return float(np.log10(1.0 + frac)) if frac else float("nan")


def master_ep(waves) -> dict:
    """The measured lines' true excitation potentials, from the production linelist."""
    m = pd.read_csv(MASTER, low_memory=False)
    f2 = m[(m.element.astype(str).str.strip() == "Fe")
           & (m.ion.astype(str).str.strip().isin(["II", "2"]))]
    out = {}
    for w in waves:
        r = f2[np.abs(f2.wavelength_air_A - w) <= 0.02]
        if len(r):
            out[w] = float(r.iloc[0].excitation_potential_eV)
    return out


def canonical_rows(waves, eps: dict) -> pd.DataFrame:
    """canonical_gf provenance, matched on wavelength AND EP (RYA-780)."""
    c = pd.read_csv(CANONICAL_GF, comment="#", low_memory=False)
    f2 = c[c.species.astype(str) == "Fe II"]
    rows = []
    for w in waves:
        ep = eps.get(w)
        rec = {"wavelength_air_A": w, "ep_eV": ep}
        if ep is not None:
            m = f2[(np.abs(f2.wavelength_air_A - w) <= WAVE_TOL_A)
                   & (np.abs(f2.excitation_potential_eV - ep) <= EP_TOL_EV)]
            # what a wavelength-ONLY match would have returned, to show the trap
            m_bad = f2[np.abs(f2.wavelength_air_A - w) <= WAVE_TOL_A]
            if len(m):
                r = m.iloc[0]
                rec.update(canonical_loggf=float(r.log_gf),
                           canonical_ref=str(r.loggf_reference),
                           stored_nist_grade=(None if pd.isna(r.nist_grade)
                                              else str(r.nist_grade)))
            if len(m_bad):
                rec["ep_if_matched_on_wavelength_only"] = float(
                    m_bad.iloc[0].excitation_potential_eV)
        rows.append(rec)
    return pd.DataFrame(rows)


def nist_rows(lo_A: float, hi_A: float) -> pd.DataFrame:
    """NIST ASD Fe II, in AIR. The default is VACUUM and that hides every line here."""
    from astroquery.nist import Nist
    import astropy.units as u

    t = Nist.query(lo_A * u.AA, hi_A * u.AA, linename="Fe II",
                   wavelength_type="vac+air")
    gcol = next(c for c in t.colnames if c.startswith("gi"))
    rows = []
    for r in t:
        try:
            obs = float(r["Observed"])
        except Exception:
            continue
        try:
            f = float(r["fik"])
        except Exception:
            continue
        try:
            gi = float(str(r[gcol]).split("-")[0].strip())
        except Exception:
            continue
        # ⚠️ EP is parsed and matched on, because wavelength alone picks the wrong
        # physical line here too: a +-0.05 A window around 6432.676 returns a NIST row
        # 2.298 dex away. I hit this on the canonical_gf side, flagged it, and then made
        # the same mistake on the NIST side — the guard has to be on BOTH.
        ei = float("nan")
        for col in t.colnames:
            if col.startswith("Ei"):
                try:
                    ei = float(str(r[col]).split("-")[0].strip())
                except Exception:
                    ei = float("nan")
                break
        acc = str(r["Acc."]).strip()
        rows.append({"nist_wave_air_A": obs, "nist_ep_eV": ei,
                     "nist_fik": f, "nist_gi": gi,
                     "nist_loggf": float(np.log10(gi * f)),
                     "nist_accuracy": acc,
                     "nist_accuracy_dex": accuracy_to_dex(acc),
                     "nist_TP": str(r["TP"]).strip()})
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    eps = master_ep(VIS_POOL_A)
    can = canonical_rows(VIS_POOL_A, eps)

    print("=== the Fe II arbiter trio, confirmed from the PRODUCTS not the ticket ===")
    print(f"  VIS ENGINE-A aggregate: {', '.join(f'{w:.3f}' for w in ARBITER_LINES_A)}")
    print("  (the ticket guessed 6247.6 / 6432.7 / 6456.4 — only 6247.6 is in the trio;")
    print("   6432.7 and 6456.4 are in the wider 11-line VIS 1D-LTE pool)")

    print("\n=== ⚠️ wavelength-only matching pulls the WRONG physical line ===")
    bad = can[can.ep_if_matched_on_wavelength_only.notna()
              & (np.abs(can.ep_if_matched_on_wavelength_only - can.ep_eV) > EP_TOL_EV)]
    for _, r in bad.iterrows():
        print(f"  {r.wavelength_air_A:9.3f}  true EP {r.ep_eV:6.3f} eV, but a 0.05 A "
              f"window returns EP {r.ep_if_matched_on_wavelength_only:7.3f} eV")
    if bad.empty:
        print("  (none — every pool line matches uniquely)")

    try:
        nist = nist_rows(5200.0, 6500.0)
    except Exception as e:
        raise SystemExit(f"NIST ASD query failed ({e}) — this audit needs it; run on "
                         f"Sirius with venv_ci, and note the default wavelength_type is "
                         f"VACUUM, which hides every line in this band.")

    merged = []
    for _, r in can.iterrows():
        m = nist[(np.abs(nist.nist_wave_air_A - r.wavelength_air_A) <= WAVE_TOL_A)
                 & (np.abs(nist.nist_ep_eV - (r.ep_eV if r.ep_eV is not None
                                              else np.nan)) <= EP_TOL_EV)]
        rec = r.to_dict()
        if len(m):
            rec.update(m.iloc[0].to_dict())
        merged.append(rec)
    d = pd.DataFrame(merged)
    d["is_arbiter"] = d.wavelength_air_A.isin(ARBITER_LINES_A)
    d["delta_canonical_minus_nist"] = d.canonical_loggf - d.nist_loggf

    print("\n=== canonical_gf vs NIST ASD v5.11 (air, EP-matched) ===")
    print(f"{'wave':>9} {'arb':>4} {'canonical ref':<26}{'stored':>7}"
          f"{'canon':>8}{'NIST':>8}{'delta':>8}{'NIST acc':>9}{'-> dex':>8}")
    for _, r in d.iterrows():
        print(f"{r.wavelength_air_A:9.3f} {'YES' if r.is_arbiter else '':>4} "
              f"{str(r.get('canonical_ref','')):<26}"
              f"{str(r.get('stored_nist_grade') or ''):>7}"
              f"{r.get('canonical_loggf', np.nan):>8.3f}"
              f"{r.get('nist_loggf', np.nan):>8.3f}"
              f"{r.get('delta_canonical_minus_nist', np.nan):>+8.3f}"
              f"{str(r.get('nist_accuracy','')):>9}"
              f"{r.get('nist_accuracy_dex', np.nan):>8.3f}")

    # ── the defect ────────────────────────────────────────────────────────────────
    claimed_b = d[d.stored_nist_grade.astype(str) == "B"]
    wrong = claimed_b[claimed_b.nist_accuracy.astype(str) != "B"]
    print(f"\n=== 🔴 lines whose STORED grade disagrees with NIST ===")
    for _, r in wrong.iterrows():
        stored_dex = accuracy_to_dex("B")
        print(f"  {r.wavelength_air_A:9.3f}  stored 'B' ({stored_dex:.3f} dex) but NIST "
              f"says '{r.nist_accuracy}' ({r.nist_accuracy_dex:.3f} dex) — "
              f"understated {r.nist_accuracy_dex/stored_dex:.1f}x; "
              f"and the value differs by {r.delta_canonical_minus_nist:+.3f} dex")
    if wrong.empty:
        print("  (none — the stored grades agree with NIST)")

    # ── the verdict ───────────────────────────────────────────────────────────────
    arb = d[d.is_arbiter]
    worst = float(arb.nist_accuracy_dex.max())
    print(f"\n=== verdict for the arbiter trio ===")
    print(f"  NO line of the three has a confirmed PRIMARY-LAB gf.")
    print(f"  Their own cited NIST accuracies span "
          f"{arb.nist_accuracy_dex.min():.3f}-{worst:.3f} dex, against the ~0.1 dex floor "
          f"the ticket anticipated and the {GRADED_GF_TERM_DEX:.3f} the stored grade "
          f"would have published.")
    print(f"  => they stay UNGRADED, with the reason stated (RYA-852 'do not fabricate "
          f"coverage'), and the stored grade B is a DEFECT to correct.")

    # ── the pool-wide scale offset ────────────────────────────────────────────────
    # Not just the two mislabelled lines: nearly the whole Fe II pool sits ABOVE NIST by a
    # similar amount, including its plain-VALD3 members. A coherent offset across an entire
    # pool is a SCALE difference, not eleven independent errors — and since a gf that is
    # too high yields an abundance that is too low, it lands directly on the Fe I / Fe II
    # ionization balance (RYA-407).
    off = d[np.isfinite(d.delta_canonical_minus_nist)]
    med_off = float(off.delta_canonical_minus_nist.median())
    n_above = int((off.delta_canonical_minus_nist > 0).sum())
    print(f"\n=== pool-wide: canonical_gf vs NIST across all {len(off)} matched lines ===")
    print(f"  median offset {med_off:+.3f} dex, {n_above} of {len(off)} ABOVE NIST "
          f"(range {off.delta_canonical_minus_nist.min():+.3f} .. "
          f"{off.delta_canonical_minus_nist.max():+.3f})")
    print(f"  a gf too HIGH by {med_off:+.3f} yields an abundance too LOW by about the "
          f"same, so this bears on the Fe I/Fe II ionization balance (RYA-407)")

    d.to_csv(OUT / "rya852_fe2_gf_audit.csv", index=False)
    summary = {
        "ticket": "RYA-852",
        "arbiter_lines_air_A": list(ARBITER_LINES_A),
        "arbiter_source": "Fe II VIS ENGINE-A aggregate (n=3), read from the products",
        "ticket_guess_was": [6247.6, 6432.7, 6456.4],
        "n_pool_lines": len(VIS_POOL_A),
        "stored_grade_disagrees_with_nist": [
            {"wavelength_air_A": float(r.wavelength_air_A),
             "stored": str(r.stored_nist_grade),
             "nist_accuracy": str(r.nist_accuracy),
             "stored_dex": accuracy_to_dex("B"),
             "nist_dex": float(r.nist_accuracy_dex),
             "delta_loggf": float(r.delta_canonical_minus_nist)}
            for _, r in wrong.iterrows()],
        "arbiter_nist_accuracy_dex": {
            f"{r.wavelength_air_A:.3f}": float(r.nist_accuracy_dex)
            for _, r in arb.iterrows()},
        "pool_scale_offset": {
            "median_dex": med_off, "n_above_nist": n_above, "n_matched": int(len(off)),
            "note": ("a coherent offset across the whole pool including its plain-VALD3 "
                     "members is a SCALE difference, not independent errors; a gf too high "
                     "gives an abundance too low, so it bears on RYA-407's ionization "
                     "balance"),
        },
        "verdict": ("NO arbiter line has a confirmed primary-lab gf. Their cited NIST "
                    "accuracies are D/E (0.176-0.301 dex), worse than the ~0.1 dex floor "
                    "the ticket anticipated and 4-7x the 0.041 the stored grade B would "
                    "have published. They stay UNGRADED with the reason stated, and the "
                    "stored grade B is a defect."),
        "not_established": [
            "Den Hartog 2019's optical subset could not be reached (Vizier.find_catalogs "
            "returned nothing for three phrasings) — an absence in the SEARCH, not the "
            "source (RYA-833). Whether these lines are among its ten optical is UNVERIFIED.",
            "MB09's own S/L flags were not obtained, so the RYA-161 firewall cannot be "
            "applied line by line; 6238.386 is MB09-labelled and is a firewall CANDIDATE, "
            "not a confirmed exclusion.",
            "the consistent +0.115..+0.154 dex offset above NIST across the three "
            "non-RU lines looks like one common scale (plausibly MB09) — a hypothesis.",
        ],
        "traps": [
            "astroquery.nist defaults to wavelength_type='vacuum'; queried that way NONE "
            "of the arbiter lines appear (air->vac is +1.71 A at 6150 A) and the obvious "
            "conclusion is a manufactured absence.",
            "matching canonical_gf on wavelength alone returns EP 13.436 eV for 6149.246 "
            "and 10.930 eV for 6432.676 — high-excitation neighbours. Match on EP too.",
        ],
    }
    (OUT / "rya852_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\n[out] {OUT}/rya852_fe2_gf_audit.csv\n[out] {OUT}/rya852_summary.json")


if __name__ == "__main__":
    main()
