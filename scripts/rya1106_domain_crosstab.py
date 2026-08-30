#!/usr/bin/env python3
"""RYA-1106 sub-items — WHICH axis puts 19 of AGSS21's 40 Fe I lines outside the Amarsi MLP,
and how that boundary lies against our own Codex Graded (laboratory-gf) pool.

    python3 scripts/rya1106_domain_crosstab.py

TWO QUESTIONS, ONE TABLE (Ryan, 2026-08-29)
-------------------------------------------
(a) The Amarsi network refuses 19 of the 40 AGSS21 Fe I lines. Ryan's hypothesis, from
    Amarsi, Liljegren & Nissen 2022 (A&A 663, A87): the wall is the WEAK-LINE restriction.
    That paper restricted its own solar Fe analysis to `REW < -4.9`, and Amarsi cautions
    that the interpolation routine covers "a rather restricted range of optical lines". If
    the 19 are the STRONG lines, the wall is the field's own weak-line discipline and there
    is nothing to fix -- the 21 ARE the reliable set.

(b) Crossover: is each of the 40 in our Codex Graded pool? Prediction: the 19 refused lines
    are DISJOINT from graded (they would be the low-Elo strong lines our 2.85 eV lab-gf
    floor already excludes), and the 21 accepted lines largely ARE the graded overlap. That
    would put Codex Graded (EW) and Asplund (MLP 3D-NLTE) both in the weak-line regime and
    leave the 19 covered by NEITHER product.

🔴 THE REW HYPOTHESIS IS REFUTED BY THE DATA, AND THE REFUTATION IS THE FINDING.
Every one of the 40 AGSS21 lines already satisfies REW < -4.9 -- the published set spans
-5.934 to -4.914, because Table A.2 IS Asplund's own weak-line selection. A cut every line
already passes cannot separate 19 from 21. The axis that actually refuses them is
`Eup - Elo`, which is `hc/lambda_vac` by construction: a WAVELENGTH statement. The network
was trained on 4787.8-6810.3 A (gt02) / 4994.1-6739.5 A (lt02), and the 19 sit outside that
span. So the wall is the network's WAVELENGTH COVERAGE, not line strength.

⚠️ `in_domain` IS A PROPERTY OF THE LINE, NOT OF THE HOLDING. It is decided from atomic data
(Elo, Eup, and their difference) against the training envelope, so all four VIS holdings
refuse the same 19. This script therefore reads the per-line artifact of a REFERENCE holding
for the domain verdict, and asserts the other three agree -- rather than quietly reporting
one holding's flags as though they were universal.

JOINS ARE lambda+EP DUAL KEY (RYA-1037), at the tolerance DERIVED from AGSS21's printed
precision (Table A.2 prints nm to 2 dp = 0.1 A, so +/-0.05 A), and the count is shipped with
a plateau sweep so the reader can see it is not a function of the window (RYA-1109).
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

from pipeline import reference_lineset as rls   # noqa: E402
from pipeline import line_match                 # noqa: E402

CANON = ROOT / "data" / "linelists" / "canonical_gf.csv"
PERLINE = ROOT / "data" / "results" / "rya1106"
OUT = ROOT / "data" / "results" / "rya1106"

#: Amarsi, Liljegren & Nissen 2022 restricted their solar Fe analysis to weak lines.
#: DISPLAY/TEST ONLY -- nothing here cuts on it; it is the hypothesis under test.
AMARSI_REW_CUT = -4.9

#: The holding whose per-line artifact supplies the domain verdict. kpno_kurucz2005 serves
#: all 40 lines, so it is the one holding that can speak for every line in the set.
REFERENCE_HOLDING = "kpno_kurucz2005"
ALL_HOLDINGS = ["kpno_kurucz2005", "kpno_molecfit", "harps_molecfit", "iag"]


def graded_pool() -> pd.DataFrame:
    """Our Codex Graded pool: Fe I at the LAB (primary-laboratory gf) tier."""
    d = pd.read_csv(CANON, low_memory=False)
    d = d[d["species"].astype(str) == "Fe I"]
    d = d[d["gf_tier"].astype(str) == "LAB"]
    return d.dropna(subset=["wavelength_air_A", "excitation_potential_eV"]).reset_index(
        drop=True)


def _match(ref: pd.DataFrame, pool: pd.DataFrame, tol: float) -> np.ndarray:
    r = line_match.match(ref["wavelength_air_A"].to_numpy(float),
                         pool["wavelength_air_A"].to_numpy(float),
                         want_ep=ref["elo_eV"].to_numpy(float),
                         src_ep=pool["excitation_potential_eV"].to_numpy(float),
                         require_ep=True, tol_A=tol)
    return np.asarray(r.index)


def axis_verdict(reason: str) -> str:
    """Collapse the network's per-line refusal text to the AXIS that did the refusing.

    ⚠️ THE CLAUSES MUST BE ANCHORED. A bare `"Elo " in reason` test also fires on the
    substring inside `"transition energy Eup-Elo "`, which reported an Elo-axis refusal on
    all 19 lines when not one of them is refused on Elo. Each axis is matched at a clause
    boundary (start of string, or after "; ") so one axis cannot be read out of another's
    name.
    """
    if not isinstance(reason, str) or not reason.strip():
        return "in domain"
    axes = []
    if re.search(r"transition energy Eup-Elo\b", reason):
        axes.append("Eup-Elo (= hc/lambda: WAVELENGTH)")
    if re.search(r"(?:^|;\s*)Eup\s+[-\d]", reason):
        axes.append("Eup")
    if re.search(r"(?:^|;\s*)Elo\s+[-\d]", reason):
        axes.append("Elo")
    if "lower level" in reason:
        axes.append("Elo level match")
    return "; ".join(dict.fromkeys(axes)) or "other"


def load_perline(holding: str) -> pd.DataFrame:
    p = PERLINE / holding / "asplund_lines_per_line.csv"
    if not p.exists():
        raise SystemExit(f"missing per-line artifact for {holding}: {p}\n"
                         "  run scripts/rya1106_asplund_replication.py first")
    return pd.read_csv(p)


def main() -> int:
    ref = load_perline(REFERENCE_HOLDING)

    # ── the domain verdict must be a property of the LINE, so prove the holdings agree ──
    agree = {}
    for h in ALL_HOLDINGS:
        p = PERLINE / h / "asplund_lines_per_line.csv"
        if not p.exists():
            agree[h] = "artifact absent"
            continue
        o = pd.read_csv(p)
        # 🔴 lambda+EP DUAL KEY (RYA-1037), not lambda alone. These frames are the same 40
        # published lines so a wavelength join happens to be unique here -- which is
        # exactly the reasoning that makes a lambda-only key survive review and then break
        # on the first blend. The row count is asserted so a key that fails to match is
        # LOUD rather than quietly shrinking the comparison to the rows that did.
        m = o.merge(ref[["wavelength_air_A", "elo_eV", "in_domain"]],
                    on=["wavelength_air_A", "elo_eV"], suffixes=("", "_ref"))
        if len(m) != len(o):
            raise SystemExit(
                f"{h}: lambda+EP join matched {len(m)} of {len(o)} lines against "
                f"{REFERENCE_HOLDING} -- the two runs do not describe the same line set")
        same = bool((m["in_domain"].astype(bool) == m["in_domain_ref"].astype(bool)).all())
        agree[h] = "agrees" if same else "DISAGREES"
    disagreeing = [h for h, v in agree.items() if v == "DISAGREES"]
    if disagreeing:
        raise SystemExit(f"in_domain differs across holdings {disagreeing} -- it is being "
                         "treated as a per-line property and is not one; STOP")

    # ── the graded join, on the dual key at the derived tolerance ─────────────────────
    pool = graded_pool()
    tol = rls.SETS["asplund"].match_tol_A
    idx = _match(ref, pool, tol)
    ref = ref.copy()
    ref["in_graded"] = idx >= 0
    ref["graded_gf"] = [None if i < 0 else float(pool["log_gf"].iloc[i]) for i in idx]
    ref["in_domain"] = ref["in_domain"].astype(bool)
    ref["axis"] = [axis_verdict(r) for r in ref["domain_reason"]]

    plateau = {str(t): int((_match(ref, pool, t) >= 0).sum())
               for t in (0.005, 0.01, 0.02, 0.03, 0.05, 0.08, 0.10, 0.25, 0.50)}

    out = ref[~ref["in_domain"]].sort_values("wavelength_air_A")
    ind = ref[ref["in_domain"]].sort_values("wavelength_air_A")

    L = []; A = L.append
    A("=" * 100)
    A("RYA-1106 sub-items — the MLP domain wall, and its crossover with Codex Graded")
    A("=" * 100)
    A(f"  line set   : AGSS21 Table A.2 Fe I, n={len(ref)}")
    A(f"  graded pool: canonical_gf Fe I at gf_tier=LAB, n={len(pool)}")
    A(f"  join       : lambda+EP dual key (RYA-1037), tol +/-{tol} A (AGSS21 prints nm to 2 dp)")
    A(f"  plateau    : {plateau}")
    A(f"  domain     : identical across all four holdings -> {agree}")
    A("")

    # ── (a) the 19, and the REW hypothesis ───────────────────────────────────────────
    A("-" * 100)
    A(f"(a) THE {len(out)} OUT-OF-DOMAIN LINES — which axis refuses them")
    A("-" * 100)
    A(f"  {'lambda_air':>11}{'Elo':>7}{'Eup':>7}{'Eup-Elo':>9}{'log gf':>8}{'REW':>8}"
      f"{'graded':>8}   axis")
    for _, r in out.iterrows():
        A(f"  {r['wavelength_air_A']:>11.3f}{r['elo_eV']:>7.3f}{r['eup_eV']:>7.3f}"
          f"{r['eup_eV'] - r['elo_eV']:>9.4f}{r['loggf_asplund']:>8.3f}"
          f"{r['rew_agss21']:>8.3f}{('YES' if r['in_graded'] else 'no'):>8}   {r['axis']}")
    A("")

    n_viol = int((ref["rew_agss21"] >= AMARSI_REW_CUT).sum())
    A(f"  REW TEST vs Amarsi's own REW < {AMARSI_REW_CUT} restriction:")
    A(f"    all 40 AGSS21 lines span REW {ref['rew_agss21'].min():.3f} .. "
      f"{ref['rew_agss21'].max():.3f}   -> {n_viol} of 40 violate REW < {AMARSI_REW_CUT}")
    A(f"    out-of-domain ({len(out)}): REW {out['rew_agss21'].min():.3f} .. "
      f"{out['rew_agss21'].max():.3f}   median {out['rew_agss21'].median():.3f}")
    A(f"    in-domain     ({len(ind)}): REW {ind['rew_agss21'].min():.3f} .. "
      f"{ind['rew_agss21'].max():.3f}   median {ind['rew_agss21'].median():.3f}")
    A(f"    out-of-domain Elo {out['elo_eV'].min():.3f}..{out['elo_eV'].max():.3f} "
      f"(median {out['elo_eV'].median():.3f});  in-domain Elo "
      f"{ind['elo_eV'].min():.3f}..{ind['elo_eV'].max():.3f} "
      f"(median {ind['elo_eV'].median():.3f})")
    A("")
    if n_viol == 0:
        A(f"  🔴 VERDICT: the weak-line hypothesis is REFUTED. Every one of the 40 lines")
        A(f"     ALREADY satisfies REW < {AMARSI_REW_CUT} -- Table A.2 IS Asplund's weak-line")
        A(f"     selection -- so that cut cannot be what separates the {len(out)} from the {len(ind)}.")
    else:
        A(f"  VERDICT: {n_viol} of 40 violate REW < {AMARSI_REW_CUT}; the cut is live on this set.")
    axis_counts = out["axis"].value_counts()
    A(f"     The axis that DOES refuse them:")
    for k, v in axis_counts.items():
        A(f"       {v:>3} x  {k}")
    A(f"     Eup-Elo is hc/lambda_vac by construction, so this is a WAVELENGTH statement:")
    A(f"     the network saw a restricted optical span and these lines fall outside it.")
    A(f"     out-of-domain lambda span {out['wavelength_air_A'].min():.1f}.."
      f"{out['wavelength_air_A'].max():.1f} A;  in-domain "
      f"{ind['wavelength_air_A'].min():.1f}..{ind['wavelength_air_A'].max():.1f} A")
    A("")

    # ── (b) the 2x2 ─────────────────────────────────────────────────────────────────
    A("-" * 100)
    A("(b) 2x2 CROSS-TAB over all 40 AGSS21 Fe I lines — {in Codex Graded} x {MLP in-domain}")
    A("-" * 100)
    cells = {}
    for g in (True, False):
        for dm in (True, False):
            cells[(g, dm)] = ref[(ref["in_graded"] == g) & (ref["in_domain"] == dm)]
    A(f"  {'':<22}{'MLP in-domain':>16}{'MLP out-of-domain':>20}{'total':>8}")
    for g in (True, False):
        lab = "in Codex Graded" if g else "NOT in Graded"
        a, b = len(cells[(g, True)]), len(cells[(g, False)])
        A(f"  {lab:<22}{a:>16}{b:>20}{a + b:>8}")
    ta = len(cells[(True, True)]) + len(cells[(False, True)])
    tb = len(cells[(True, False)]) + len(cells[(False, False)])
    A(f"  {'total':<22}{ta:>16}{tb:>20}{ta + tb:>8}")
    A("")

    strong_and_graded = cells[(True, False)]
    if len(strong_and_graded) == 0:
        A("  ✅ PREDICTION HOLDS: the out-of-domain lines are DISJOINT from Codex Graded.")
        A("     Codex Graded (EW, lab-gf) and Asplund (MLP 3D-NLTE) both live inside the")
        A(f"     network's optical span; the {len(out)} refused lines are covered by NEITHER")
        A("     product -- they are Bride-only.")
    else:
        A(f"  🔴 FINDING — {len(strong_and_graded)} line(s) are out-of-domain AND graded, which the")
        A("     prediction said should be empty. Named in full:")
        for _, r in strong_and_graded.iterrows():
            A(f"       {r['wavelength_air_A']:.3f}  Elo {r['elo_eV']:.3f}  "
              f"REW {r['rew_agss21']:.3f}  axis: {r['axis']}")
    A("")
    A("  per-line detail, all 40:")
    A(f"  {'lambda_air':>11}{'Elo':>7}{'log gf':>8}{'REW':>8}{'graded':>8}{'domain':>9}")
    for _, r in ref.sort_values("wavelength_air_A").iterrows():
        A(f"  {r['wavelength_air_A']:>11.3f}{r['elo_eV']:>7.3f}{r['loggf_asplund']:>8.3f}"
          f"{r['rew_agss21']:>8.3f}{('YES' if r['in_graded'] else 'no'):>8}"
          f"{('in' if r['in_domain'] else 'OUT'):>9}")
    A("=" * 100)

    text = "\n".join(L)
    print(text)

    OUT.mkdir(parents=True, exist_ok=True)
    ref_out = ref[["wavelength_air_A", "elo_eV", "eup_eV", "loggf_asplund", "rew_agss21",
                   "in_graded", "graded_gf", "in_domain", "axis", "domain_reason"]]
    ref_out.to_csv(OUT / "rya1106_domain_crosstab.csv", index=False)
    doc = {
        "generated_by": "scripts/rya1106_domain_crosstab.py",
        "line_set": "asplund_agss21 Fe I (AGSS21 Table A.2)",
        "n_lines": int(len(ref)),
        "graded_pool": {"source": "data/linelists/canonical_gf.csv",
                        "filter": "species == 'Fe I' and gf_tier == 'LAB'",
                        "n": int(len(pool))},
        "join": {"key": "lambda+EP dual (RYA-1037)", "tol_A": tol, "plateau": plateau},
        "domain_agreement_across_holdings": agree,
        "rew_test": {
            "amarsi_cut": AMARSI_REW_CUT,
            "n_violating_cut": n_viol,
            "verdict": ("REFUTED -- every line already satisfies the cut"
                        if n_viol == 0 else "cut is live on this set"),
            "all": [round(float(ref["rew_agss21"].min()), 4),
                    round(float(ref["rew_agss21"].max()), 4)],
            "out_of_domain": [round(float(out["rew_agss21"].min()), 4),
                              round(float(out["rew_agss21"].max()), 4)],
            "in_domain": [round(float(ind["rew_agss21"].min()), 4),
                          round(float(ind["rew_agss21"].max()), 4)],
        },
        "axis_counts": {str(k): int(v) for k, v in axis_counts.items()},
        "crosstab": {"graded_and_in_domain": len(cells[(True, True)]),
                     "graded_and_out_of_domain": len(cells[(True, False)]),
                     "not_graded_and_in_domain": len(cells[(False, True)]),
                     "not_graded_and_out_of_domain": len(cells[(False, False)])},
        "prediction_holds": bool(len(strong_and_graded) == 0),
    }
    (OUT / "rya1106_domain_crosstab.json").write_text(json.dumps(doc, indent=2) + "\n")
    (OUT / "rya1106_domain_crosstab.txt").write_text(text + "\n")
    print(f"\nwrote {OUT / 'rya1106_domain_crosstab.json'}")
    print(f"wrote {OUT / 'rya1106_domain_crosstab.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
