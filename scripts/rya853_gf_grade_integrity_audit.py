#!/usr/bin/env python3
"""
RYA-853 — do the stored gf grades and values match the source they cite?
========================================================================
RYA-852 found two Fe II lines stored as "NIST ASD v5.11 grade B" that NIST grades E and D,
with log gf +0.130 and +0.115 above NIST's own values. This asks how far that goes.

🔴 IT IS A KNOWN DEFECT CLASS, ALREADY FOUND ONCE AND FIXED FOR TWO ROWS. The header of
`data/linelists/nist_reference.csv` records that RYA-592 re-verified the **Mg I rows only**
and found the same thing:

    "the log_gf values RE-DERIVE EXACTLY from the live source, but two columns were wrong
     ... the accuracy CODE (both rows carried A; ASD reports B+ and B) and aki_s-1 (both
     rows disagreed with their own log_gf). ... The OTHER rows in this file are NOT
     re-verified against v5.12 — the drift found here means they should not be assumed
     current."

The file asks for exactly this audit. The Fe II rows are the next instance of a defect the
repo already knew it had.

WHY IT IS NOT COSMETIC. This file is the "Type B uncertainty anchor". RYA-850 keys its
graded gf term on the stored grade, so a fabricated `B` publishes 0.041 dex on a line whose
source says 0.176-0.301. A grade is a claim about how well a number is known, and a wrong
one is not a labelling nit — it is an understated error bar with a citation attached.

⚠️ MY OWN FIRST DIAGNOSIS WAS WRONG, AND IS RECORDED BECAUSE IT IS THE OBVIOUS ONE.
`rya347_fe2_atomic_data_audit.py:108` matches NIST on wavelength alone within 0.1 A, with
no excitation-potential guard, so a wrong-line match looked like the cause. It is not: the
EP in these files is CORRECT (3.889 / 3.892), the right physical line was matched, and the
stored grade and value are simply not what NIST says. The loose match is a latent risk
worth fixing on its own, but it did not produce this.

THE TWO CROSS-MATCH GUARDS (RYA-846/852), audited across the repo:
  (a) AIR vs VACUUM — `astroquery.nist` defaults to VACUUM. At 6150 A that is +1.71 A, so a
      naive query returns "not covered" for every line in the band.
  (b) EP ON BOTH SIDES — a wavelength-only window returns the wrong level (EP 13.436 eV for
      Fe II 6149.246; a -2.298 dex artifact for 6432.676 until the NIST side matched on EP
      too, which I got wrong in RYA-852 after flagging it).

Usage (Sirius — needs astroquery in venv_ci):
    python3 scripts/rya853_gf_grade_integrity_audit.py
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

NIST_REF = ROOT / "data" / "linelists" / "nist_reference.csv"
NIST_XC = ROOT / "data" / "linelists" / "nist_crosscheck.csv"
OUT = ROOT / "data" / "results" / "rya853"

#: NIST ASD accuracy letters -> fractional uncertainty on the transition probability.
NIST_ACCURACY_FRACTION = {
    "AAA": 0.003, "AA": 0.01, "A+": 0.02, "A": 0.03, "B+": 0.07, "B": 0.10,
    "C+": 0.18, "C": 0.25, "D+": 0.40, "D": 0.50, "E": 1.00,
}

#: Matching tolerances. EP is REQUIRED on both sides (guard b).
WAVE_TOL_A = 0.05
EP_TOL_EV = 0.05

#: A stored value this far from its cited source is a transcription defect, not rounding.
VALUE_TOL_DEX = 0.02


def accuracy_to_dex(letter: str) -> float:
    frac = NIST_ACCURACY_FRACTION.get(str(letter).strip())
    return float(np.log10(1.0 + frac)) if frac else float("nan")


def load_stored() -> pd.DataFrame:
    """Both hand-maintained NIST extracts, tagged by which file they came from."""
    frames = []
    for path, tag in ((NIST_REF, "nist_reference"), (NIST_XC, "nist_crosscheck")):
        if not path.exists():
            continue
        d = pd.read_csv(path, comment="#")
        d["stored_file"] = tag
        frames.append(d)
    if not frames:
        raise SystemExit("neither NIST extract is present")
    return pd.concat(frames, ignore_index=True)


def nist_live(species: str, lo_A: float, hi_A: float) -> pd.DataFrame:
    """NIST ASD in AIR, carrying EP so the caller can match on it (guards a and b)."""
    from astroquery.nist import Nist
    import astropy.units as u

    t = Nist.query(lo_A * u.AA, hi_A * u.AA, linename=species,
                   wavelength_type="vac+air")
    gcol = next((c for c in t.colnames if c.startswith("gi")), None)
    eicol = next((c for c in t.colnames if c.startswith("Ei")), None)
    rows = []
    for r in t:
        def _f(col):
            try:
                return float(str(r[col]).strip())
            except Exception:
                return float("nan")
        obs, ritz = _f("Observed"), _f("Ritz")
        wave = obs if np.isfinite(obs) else ritz
        if not np.isfinite(wave):
            continue
        fik = _f("fik")
        try:
            gi = float(str(r[gcol]).split("-")[0].strip()) if gcol else float("nan")
        except Exception:
            gi = float("nan")
        try:
            ei = float(str(r[eicol]).split("-")[0].strip()) if eicol else float("nan")
        except Exception:
            ei = float("nan")
        rows.append({
            "nist_wave_A": wave, "nist_ep_eV": ei, "nist_fik": fik, "nist_gi": gi,
            "nist_loggf": (float(np.log10(gi * fik))
                           if np.isfinite(gi) and np.isfinite(fik) and gi * fik > 0
                           else float("nan")),
            "nist_accuracy": str(r["Acc."]).strip(),
        })
    return pd.DataFrame(rows)


def audit_cross_match_guards() -> dict:
    """Scope 4: does every NIST cross-match in the repo carry both guards?"""
    findings = []
    for p in sorted(ROOT.glob("scripts/*.py")):
        src = p.read_text(errors="ignore")
        if "Nist.query" not in src:
            continue
        findings.append({
            "file": p.name,
            "air_guard": "wavelength_type" in src and "vac+air" in src,
            "ep_guard": bool(re.search(r"(excitation|Ei|ep_eV|elo)", src)),
        })
    # the local-extract matchers do not call NIST but still cross-match on wavelength
    for p in sorted(ROOT.glob("scripts/*.py")):
        src = p.read_text(errors="ignore")
        if "nist_reference" not in src and "nist_crosscheck" not in src:
            continue
        loose = re.findall(r"wavelength_air_A\s*-\s*\w+\)\.abs\(\)\s*<\s*([0-9.]+)", src)
        if loose:
            findings.append({
                "file": p.name, "air_guard": None,
                "ep_guard": bool(re.search(r"excitation_potential_eV\s*-", src)),
                "wavelength_only_window_A": [float(x) for x in loose],
            })
    return {"sites": findings}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    stored = load_stored()
    print(f"[stored] {len(stored)} rows across "
          f"{stored.stored_file.nunique()} hand-maintained NIST extracts")

    guards = audit_cross_match_guards()
    print(f"\n=== scope 4: cross-match guards ===")
    for s in guards["sites"]:
        bits = []
        if s.get("air_guard") is not None:
            bits.append(f"air/vac={'OK' if s['air_guard'] else 'MISSING'}")
        bits.append(f"EP={'OK' if s['ep_guard'] else 'MISSING'}")
        if s.get("wavelength_only_window_A"):
            bits.append(f"wavelength-only window {s['wavelength_only_window_A']} A")
        print(f"  {s['file']:<44} {'  '.join(bits)}")

    # ── verify every stored row against live NIST ────────────────────────────────
    rows = []
    for (el, ion), grp in stored.groupby(["element", "ion"]):
        species = f"{el} {ion}"
        lo = float(grp.wavelength_air_A.min()) - 5
        hi = float(grp.wavelength_air_A.max()) + 5
        try:
            live = nist_live(species, lo, hi)
        except Exception as e:
            print(f"  [{species}] NIST query FAILED: {type(e).__name__}: {str(e)[:80]}")
            live = pd.DataFrame()
        for _, r in grp.iterrows():
            rec = {
                "stored_file": r.stored_file, "element": el, "ion": ion,
                "wavelength_air_A": float(r.wavelength_air_A),
                "ep_eV": float(r.excitation_potential_eV),
                "stored_loggf": float(r.log_gf),
                "stored_grade": str(r.nist_grade).strip(),
                "stored_grade_dex": accuracy_to_dex(r.nist_grade),
            }
            if len(live):
                m = live[(np.abs(live.nist_wave_A - rec["wavelength_air_A"]) <= WAVE_TOL_A)
                         & (np.abs(live.nist_ep_eV - rec["ep_eV"]) <= EP_TOL_EV)]
                if len(m):
                    mm = m.iloc[0]
                    rec.update(nist_loggf=float(mm.nist_loggf),
                               nist_accuracy=str(mm.nist_accuracy),
                               nist_accuracy_dex=accuracy_to_dex(mm.nist_accuracy))
            rows.append(rec)

    d = pd.DataFrame(rows)
    d["grade_matches"] = d.stored_grade == d.get("nist_accuracy")
    d["delta_loggf"] = d.stored_loggf - d.get("nist_loggf")
    d["value_matches"] = d.delta_loggf.abs() <= VALUE_TOL_DEX

    matched = d[d.get("nist_accuracy").notna()] if "nist_accuracy" in d else d.iloc[0:0]
    print(f"\n=== stored vs live NIST ASD ({len(matched)} of {len(d)} rows matched) ===")
    bad_grade = matched[~matched.grade_matches]
    bad_value = matched[~matched.value_matches]
    print(f"  grade mismatches: {len(bad_grade)} of {len(matched)}")
    print(f"  value mismatches (> {VALUE_TOL_DEX} dex): {len(bad_value)} of {len(matched)}")

    if len(bad_grade):
        print(f"\n=== 🔴 stored grade disagrees with NIST ===")
        print(f"{'elem':<7}{'wave':>10}{'stored':>8}{'NIST':>7}"
              f"{'stored dex':>11}{'NIST dex':>10}{'ratio':>7}  file")
        for _, r in bad_grade.sort_values("wavelength_air_A").iterrows():
            ratio = (r.nist_accuracy_dex / r.stored_grade_dex
                     if r.stored_grade_dex and np.isfinite(r.stored_grade_dex) else np.nan)
            print(f"{r.element+' '+str(r.ion):<7}{r.wavelength_air_A:>10.3f}"
                  f"{r.stored_grade:>8}{str(r.nist_accuracy):>7}"
                  f"{r.stored_grade_dex:>11.3f}{r.nist_accuracy_dex:>10.3f}"
                  f"{ratio:>7.1f}x  {r.stored_file}")

    if len(bad_value):
        print(f"\n=== 🔴 stored log gf disagrees with NIST ===")
        for _, r in bad_value.sort_values("wavelength_air_A").iterrows():
            print(f"  {r.element} {r.ion} {r.wavelength_air_A:9.3f}  "
                  f"stored {r.stored_loggf:+.3f}  NIST {r.nist_loggf:+.3f}  "
                  f"delta {r.delta_loggf:+.3f}  [{r.stored_file}]")

    d.to_csv(OUT / "rya853_grade_integrity.csv", index=False)
    summary = {
        "ticket": "RYA-853",
        "n_stored_rows": int(len(d)),
        "n_matched_to_nist": int(len(matched)),
        "n_grade_mismatch": int(len(bad_grade)),
        "n_value_mismatch": int(len(bad_value)),
        "grade_mismatches": [
            {"element": f"{r.element} {r.ion}",
             "wavelength_air_A": float(r.wavelength_air_A),
             "stored": r.stored_grade, "nist": str(r.nist_accuracy),
             "stored_dex": float(r.stored_grade_dex),
             "nist_dex": float(r.nist_accuracy_dex),
             "understated_by": float(r.nist_accuracy_dex / r.stored_grade_dex)
             if r.stored_grade_dex else None,
             "file": r.stored_file}
            for _, r in bad_grade.iterrows()],
        "value_mismatches": [
            {"element": f"{r.element} {r.ion}",
             "wavelength_air_A": float(r.wavelength_air_A),
             "stored_loggf": float(r.stored_loggf),
             "nist_loggf": float(r.nist_loggf),
             "delta": float(r.delta_loggf), "file": r.stored_file}
            for _, r in bad_value.iterrows()],
        "cross_match_guards": guards,
        "known_precedent": ("nist_reference.csv's own header records that RYA-592 "
                            "re-verified the Mg I rows ONLY and found the accuracy CODE "
                            "wrong there too (stored A, ASD says B+ and B), and warns "
                            "that the other rows 'should not be assumed current'. This "
                            "audit is what that warning asked for."),
        "corrected_diagnosis": ("rya347_fe2_atomic_data_audit.py matches NIST on "
                                "wavelength alone within 0.1 A with no EP guard, which "
                                "looks like the cause and is NOT: the EP in these files "
                                "is correct and the right line was matched. The stored "
                                "grade and value are simply not what NIST says. The loose "
                                "match is a latent risk worth fixing separately."),
    }
    (OUT / "rya853_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\n[out] {OUT}/rya853_grade_integrity.csv\n[out] {OUT}/rya853_summary.json")


if __name__ == "__main__":
    main()
