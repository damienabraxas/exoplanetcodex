#!/usr/bin/env python3
"""RYA-1214 Step 1 — adjudicate the NIST ASD CNO gf into `canonical_gf.csv`.

WHAT THIS CLOSES
----------------
RYA-1160 pulled 5,119 graded NIST ASD transitions for C I/II, N I/II, O I/II, validated
the pull against the store (values reproduce to <=0.001 dex on every uniquely matched
control row) and then stopped: `adjudicated: false`, `canonical_gf_modified: false`.
So the gf sat on disk in `data/linelists/primary_gf/` while every CNO line in
`canonical_gf.csv` still carried an UNGRADED Kurucz or VALD value. This script is the
bridge.

WHAT A CNO gf CAN AND CANNOT BE (RYA-1172)
------------------------------------------
CNO has ZERO primary-laboratory gf in this repo and there is no lab table to register --
`pipeline.gf_grades.LAB_TABLES` holds Fe I, Fe II and Al I and nothing else. The
community standard for the light elements is CRITICALLY-EVALUATED THEORY:

    C I, C II, N I, N II   MCHF (Froese Fischer & Tachiev), NIST partial update 2006
    O I, O II              Opacity Project, Wiese Fuhr & Deters 1996 Monograph 7

Both are calculations. `gf_tier` here is therefore `NIST-C+` -- the compilation rung of
the RYA-946 ladder (lab > NIST-C+ > Kurucz > VALD) -- and NEVER `LAB`. The per-species
pedigree that distinguishes a C grade from an O grade is `pipeline.cno_gf_pedigree`, and
this script refuses to adjudicate any species that module does not state a pedigree for
(C III-C V, N III-N V, O III-O VI): an unknown pedigree stays unknown rather than being
extended from a neighbouring species (RYA-1072).

THE MATCH, AND WHY ITS TOLERANCES ARE MEASURED AND NOT CHOSEN
--------------------------------------------------------------
Wavelength AND excitation potential, unique or nothing (RYA-780). The tolerances are read
off a bimodality in the residuals, not picked. Taking, for every CNO row, the NIST
candidate nearest in wavelength among those agreeing in EP to 0.001 eV:

    |dlambda|          n     median |d log gf|     fraction > 0.1 dex
    0.000-0.005 A    529          0.001                  0.08
    0.005-0.010 A    101          0.002                  0.07
    0.010-0.020 A    113          0.003                  0.06
    0.020-0.050 A    134          0.002                  0.07
    0.050-0.100 A     90          0.004                  0.10
    > 0.100 A       1151          1.046                  0.87

A real identification reproduces the other catalogue's log gf; a coincidence does not.
The break is at ~0.1 A and it is a factor of ~300 in the median, so WTOL_A = 0.05 sits
inside the flat region with a wide margin and is nowhere near the wall. That is wider
than RYA-822's 0.006 A for near-UV Fe I -- deliberately, and for a measured reason: this
band is ~0.1 CNO lines per Angstrom where the near-UV Fe forest is ~5.6, so the
coincidence rate that forced the tight Fe window does not exist here. `--null-trials`
measures it rather than asserting it.

EP is the load-bearing coordinate and it is TIGHT (0.001 eV). `canonical_gf` quotes EP to
3 decimals, so 0.001 eV is the quantisation itself; it separates multiplets sharing a
wavelength, which is the failure mode a loose EP invites.

WHAT IS DELIBERATELY NOT TOUCHED
--------------------------------
* **The Tachiev/MCHF O I values.** RYA-1160 derived Tachiev log gf for the 777 triplet
  and the 8446 triplet from Tachiev's own A-values and found a SYSTEMATIC -0.016 dex
  (777) and -0.005 dex (8446) offset against NIST. AGSS21 names Tachiev for N i ONLY and
  states no gf source for O i, so which source AGSS21 used for oxygen is OPEN. This
  script REPORTS both values side by side in `oi_tachiev_vs_nist.csv` and adopts
  neither over the other -- Ryan decides. Silently resolving it would pick the dominant
  oxygen indicator's scale by accident.
* **Rows already adjudicated to NIST.** `adjudicated_rya367` ([O I] 6300) and the three
  `nist` rows from RYA-1171 keep their history; where this pull agrees that is recorded
  as corroboration, not rewritten.
* **N I is adjudicated but gains nothing independent.** NIST's own N I source IS Tachiev
  (TP code T7370 on all five AGSS21 lines), so the pull and the paper are one source seen
  twice. Recorded in the provenance so the agreement is never read as confirmation.

WRITING
-------
Line surgery: the file is read as raw lines, only changed rows are re-serialised, and
every untouched line stays byte-identical. A pandas round-trip of `canonical_gf.csv`
rewrites float formatting on rows nobody meant to touch (RYA-853), so it is not used.
The round-trip of the UNCHANGED lines is asserted byte-identical before anything is
written, which is what makes the re-serialisation of the changed ones safe.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.cno_gf_pedigree import pedigree_for, UNESTABLISHED  # noqa: E402
from pipeline.gf_grades import nist_sigma_dex, NIST_ACC_PCT  # noqa: E402

CANON = ROOT / "data" / "linelists" / "canonical_gf.csv"
PULLS = ROOT / "data" / "linelists" / "primary_gf"
TACHIEV = ROOT / "data" / "audit" / "rya1160_cno_nist_gf" / "tachiev_a_to_gf.csv"
OUT = ROOT / "data" / "audit" / "rya1214_cno_products"

#: The six species a held source states a pedigree for. Higher ionisation stages are
#: NOT adjudicated -- `cno_gf_pedigree.UNESTABLISHED` is the honest answer for them.
SPECIES_FILES = {
    "C I": "CI", "C II": "CII", "N I": "NI", "N II": "NII", "O I": "OI", "O II": "OII",
}
BANDS = ("900_3000", "3000_25000")

#: Read off the bimodality in the module docstring. Not chosen, and swept by --null-trials.
WTOL_A = 0.05
EPTOL_EV = 0.001

#: The compilation rung of RYA-946's ladder. Never LAB: this is evaluated theory.
TIER = "NIST-C+"
STATUS = "nist_rya1214"

#: Statuses that record a decision already taken about this line's gf. This pull does not
#: overwrite them; where it agrees, that is corroboration and is reported as such.
PROTECTED_STATUS = {"adjudicated_rya367", "nist", "nist_rya822", "nist_rya945",
                    "lab_rya945", "lab_rya834", "lab_rya1047", "isotope_rya1075"}

#: Above this the NIST class is a WIDER bar than the RYA-161 Kurucz systematic (0.20 dex).
#: Adopted anyway -- the rung is about pedigree, and a measured-but-wide bar beats an
#: assumed one (RYA-799) -- but counted and reported, never hidden.
K07_SYSTEMATIC_DEX = 0.20

#: 🔴 THE PHYSICALITY GATE, AND IT IS NOT A THRESHOLD I CHOSE.
#:
#: 39 rows of the RYA-1160 pull carry an oscillator strength of 8e5 to 5e7. A bound-bound
#: absorption f-value is bounded by the oscillator-strength sum rule and is never above
#: ~2, so `log gf = log10(gi * fik)` came out at +6.6 to +7.8 on those rows -- gf > 4
#: MILLION. Three of them (C II 5039.639, 5207.982, 5717.984) sit on a real canonical line
#: and would have been ADOPTED, moving the store's gf by +8.6 to +9.7 dex. RYA-1160
#: validated the pull on 5 control rows, all of which reproduced, so nothing caught it.
#:
#: The discriminator is the Ladenburg identity, which the pull's own columns over-determine:
#:
#:      gi * fik  ==  1.4992e-16 * lambda[A]^2 * gk * Aki
#:
#: Solving it for gk must return a STATISTICAL WEIGHT. Measured over all 5,119 graded rows:
#:
#:      5,080 rows   gk_implied in 0.89 - 12.02, clustering on the integers 1..12
#:         39 rows   gk_implied in 8.2e8 - 9.9e10
#:
#: Eight orders of magnitude of separation, and the bad set is exactly the 39. They are all
#: C II, all NIST accuracy class E, and all carry the SAME transition-probability reference
#: code `T10528` -- so this is one upstream delivery, not scattered noise.
#:
#: ⚠️ REFUSED, NOT REPAIRED. gk is not in the tidy pull, so the identity bounds the row
#: without recovering it; inventing a gk to back out a "corrected" gf would be fabricating
#: the datum. The rows are quarantined, counted, and reported (RYA-711: quarantine, never
#: cull), and the canonical value they would have overwritten stands untouched.
GK_MIN, GK_MAX = 0.5, 50.0


def load_nist(species: str) -> pd.DataFrame:
    tag = SPECIES_FILES[species]
    frames = []
    for band in BANDS:
        p = PULLS / f"nist_asd_{tag}_{band}.tsv"
        if not p.exists():
            raise SystemExit(f"missing NIST pull {p} — re-run scripts/rya1160_pull_nist_cno.py")
        frames.append(pd.read_csv(p, sep="\t"))
    n = pd.concat(frames, ignore_index=True)
    n = n[n.log_gf.notna() & n.ei_eV.notna() & n.wavelength_A.notna()]
    # A graded row is one carrying a NIST accuracy class we can price. An unmapped class
    # would become a NaN sigma wearing a tier, which is worse than staying ungraded.
    n = n[n.nist_grade.isin(NIST_ACC_PCT)]
    n = n.reset_index(drop=True)
    n["gk_implied"] = (n.gi * n.fik
                       / (1.4992e-16 * n.wavelength_A ** 2 * n["aki_s-1"]))
    return n


def split_nonphysical(n: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Ladenburg gate — see GK_MIN/GK_MAX. Returns (usable, quarantined)."""
    ok = n.gk_implied.between(GK_MIN, GK_MAX)
    return n[ok].reset_index(drop=True), n[~ok].reset_index(drop=True)


def match(canon_sp: pd.DataFrame, nist: pd.DataFrame) -> pd.DataFrame:
    """(wavelength, EP) match, admitted ONLY when the MULTIPLICITIES AGREE.

    🔴 A UNIQUE MATCH IS NOT ENOUGH, AND ASSUMING IT WAS WOULD HAVE BIASED THE STORE LOW.

    `canonical_gf` is one row per PHYSICAL line: `cluster_physical_lines` (RYA-353)
    collapses the VALD components of one feature onto a gf-weighted centroid and records
    how many it collapsed in `hfs_n_components`. NIST ASD serves those components
    SEPARATELY. So the two catalogues disagree about what a row is, and the match has to
    reconcile that rather than ignore it:

      * `hfs_n_components == n_candidates`  the same set of transitions on both sides.
        Adopt `log10(sum 10**log_gf)` -- the gf of a blend is the SUM of its components,
        which is the same convention the canonical row was built with.
      * anything else                        MULTIPLICITY MISMATCH. Refused.

    The rule is stated on multiplicity alone, and the AGREEMENT is then the test of it:
    where it fires with n > 1 it reproduces the store's existing value to a median
    0.0006 dex (max 0.033) over 45 rows. Nothing selects a candidate by how close its gf
    is to the value being replaced -- that would be choosing the answer.

    What the rule refuses is exactly where the damage was: 13 rows have a SINGLE NIST
    candidate against `hfs_n_components >= 2`, i.e. NIST served one component of a feature
    we hold as a sum. Their deltas are systematically NEGATIVE -- 12 of 13 below zero,
    the hfs=3 case by -0.346 dex -- because a component is not a blend. Adopting them on
    uniqueness alone would have written a systematic under-count into the oxygen and
    carbon pools with full provenance attached.
    """
    nw = nist.wavelength_A.to_numpy(float)
    ne = nist.ei_eV.to_numpy(float)
    ng = nist.log_gf.to_numpy(float)
    out = []
    for idx, r in canon_sp.iterrows():
        w, ep = float(r.wavelength_air_A), float(r.excitation_potential_eV)
        hfs = int(r.hfs_n_components) if pd.notna(r.hfs_n_components) else 1
        cand = np.where((np.abs(nw - w) <= WTOL_A) & (np.abs(ne - ep) <= EPTOL_EV))[0]
        rec = {"row": idx, "line_id": r.line_id, "species": r.species,
               "wavelength_air_A": w, "excitation_potential_eV": ep,
               "hfs_n_components": hfs,
               "log_gf_before": float(r.log_gf),
               "loggf_reference_before": r.loggf_reference,
               "gf_tier_before": r.gf_tier,
               "adjudication_status_before": r.adjudication_status,
               "n_candidates": int(cand.size),
               "verdict": ""}
        if cand.size == 0:
            rec["verdict"] = "NO_NIST_ROW"
        elif cand.size != hfs:
            rec["verdict"] = f"MULTIPLICITY_MISMATCH_{cand.size}_vs_hfs_{hfs}"
        else:
            gf_sum = float(np.log10(np.sum(10.0 ** ng[cand])))
            grades = [str(nist.nist_grade[int(j)]) for j in cand]
            # The blend is only as well known as its WORST component. NIST_ACC_PCT is the
            # full ladder including the '+' tiers, so "worst" is by percentage, never by
            # letter order (RYA-592's inverted ladder).
            worst = max(grades, key=lambda g: NIST_ACC_PCT[g])
            rec.update({
                "nist_log_gf": gf_sum,
                "nist_grade": worst,
                "nist_grades_all": "|".join(grades),
                "nist_acc_pct": NIST_ACC_PCT[worst],
                "nist_tp_code": "|".join(sorted({str(nist.ref_transition_probability[int(j)])
                                                 for j in cand})),
                "nist_wavelength_A": float(nw[cand].mean()),
                "d_wavelength_A": float(nw[cand].mean() - w),
                "d_ep_eV": float(ne[cand].mean() - ep),
                "delta_nist_minus_store": gf_sum - float(r.log_gf),
                "verdict": "ADOPT" if cand.size == 1 else f"ADOPT_SUM_{cand.size}",
            })
        out.append(rec)
    d = pd.DataFrame(out)
    d["matched"] = d.verdict.str.startswith("ADOPT")
    return d


def null_rate(canon_sp: pd.DataFrame, nist: pd.DataFrame, trials: int, seed: int) -> float:
    """How often does THIS matcher fire on randomised wavelengths?

    EP is kept real and only the wavelength is randomised. That isolates exactly the
    coincidence being worried about and is the conservative direction -- randomising both
    would report a lower null than the real hazard.
    """
    nw, ne = nist.wavelength_A.to_numpy(float), nist.ei_eV.to_numpy(float)
    ep = canon_sp.excitation_potential_eV.to_numpy(float)
    lo, hi = float(canon_sp.wavelength_air_A.min()), float(canon_sp.wavelength_air_A.max())
    if not np.isfinite(lo) or hi <= lo:
        return float("nan")
    rng = np.random.default_rng(seed)
    rates = []
    for _ in range(trials):
        w = rng.uniform(lo, hi, len(ep))
        hits = sum(1 for k in range(len(w))
                   if np.count_nonzero((np.abs(nw - w[k]) <= WTOL_A)
                                       & (np.abs(ne - ep[k]) <= EPTOL_EV)) == 1)
        rates.append(hits / len(w))
    return float(np.mean(rates))


# --------------------------------------------------------------------------- writing


def _row_to_line(fields: list[str]) -> str:
    buf = io.StringIO()
    csv.writer(buf, lineterminator="\n").writerow(fields)
    return buf.getvalue()


def apply_line_surgery(edits: dict[int, dict[str, str]]) -> tuple[int, int]:
    """Rewrite ONLY the data lines in `edits` (0-based data-row index -> {col: value}).

    Asserts first that re-serialising an UNCHANGED line reproduces it byte for byte. If
    that fails the writer's quoting convention differs from the file's and no edit is
    safe, so it raises rather than writing a file that differs on rows nobody touched.
    """
    with CANON.open("r", newline="") as fh:
        raw = fh.read().split("\n")
    trailing_newline = raw[-1] == ""
    if trailing_newline:
        raw = raw[:-1]
    header = raw[0]
    cols = next(csv.reader([header]))
    n_data = len(raw) - 1

    # Control: the round-trip must be the identity on lines we are NOT changing.
    probe = [i for i in range(n_data) if i not in edits][:2000]
    for i in probe:
        line = raw[i + 1]
        if _row_to_line(next(csv.reader([line]))).rstrip("\n") != line:
            raise SystemExit(
                f"STOP: csv round-trip is not the identity on data row {i} — the writer's "
                f"quoting differs from the file's, so a surgical edit would silently "
                f"rewrite rows nobody asked to change (RYA-853).")

    changed = 0
    for i, patch in edits.items():
        fields = next(csv.reader([raw[i + 1]]))
        for col, val in patch.items():
            fields[cols.index(col)] = val
        new = _row_to_line(fields).rstrip("\n")
        if new != raw[i + 1]:
            raw[i + 1] = new
            changed += 1
    text = "\n".join(raw) + ("\n" if trailing_newline else "")
    with CANON.open("w", newline="") as fh:
        fh.write(text)
    return changed, n_data


def _fmt(x: float) -> str:
    """Match the file's own convention: plain repr, no scientific notation, no padding."""
    s = f"{x:.6f}".rstrip("0").rstrip(".")
    return s if s not in ("", "-0") else "0"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true", help="write canonical_gf.csv")
    ap.add_argument("--null-trials", type=int, default=20)
    a = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    canon = pd.read_csv(CANON, low_memory=False)

    all_matches, per_species, edits, quarantined = [], [], {}, []
    for sp in SPECIES_FILES:
        ped = pedigree_for(sp)
        if ped.key == UNESTABLISHED.key:      # structurally unreachable for these six
            raise SystemExit(f"{sp} has no stated pedigree — refusing to adjudicate it")
        c = canon[canon.species == sp]
        if c.empty:
            continue
        nist_all = load_nist(sp)
        nist, bad = split_nonphysical(nist_all)
        if len(bad):
            bad = bad.assign(species=sp)
            quarantined.append(bad)
        m = match(c, nist)
        m["gf_authority"] = ped.key
        m["gf_method"] = ped.method
        m["gf_vintage"] = ped.vintage
        m["is_laboratory"] = ped.is_laboratory
        null = null_rate(c, nist, a.null_trials, seed=1214)
        hit = float(m.matched.mean())
        n_amb = int(m.verdict.str.startswith("MULTIPLICITY_MISMATCH").sum())
        n_hit = int(m.matched.sum())

        protected = m.matched & m.adjudication_status_before.isin(PROTECTED_STATUS)
        adopt = m.matched & ~protected
        for _, r in m[adopt].iterrows():
            edits[int(r.row)] = {
                "log_gf": _fmt(r.nist_log_gf),
                "loggf_reference": f"NIST ASD v5.11 grade {r.nist_grade} ({r.nist_tp_code})",
                "nist_grade": str(r.nist_grade),
                "adjudication_status": STATUS,
                "gf_sigma_dex": _fmt(nist_sigma_dex(r.nist_grade)),
                "gf_tier": TIER,
            }
        sig = m.loc[adopt, "nist_grade"].map(nist_sigma_dex)
        per_species.append({
            "species": sp,
            "gf_authority": ped.key,
            "method": ped.method.split(" --")[0].split(" —")[0][:60],
            "vintage": ped.vintage,
            "is_laboratory": ped.is_laboratory,
            "canonical_rows": int(len(c)),
            "nist_graded_rows_in_pull": int(len(nist_all)),
            "nist_rows_nonphysical_quarantined": int(len(nist_all) - len(nist)),
            "nist_rows_usable": int(len(nist)),
            "matched": n_hit,
            "adopt_single": int((m.verdict == "ADOPT").sum()),
            "adopt_summed_blend": int(m.verdict.str.startswith("ADOPT_SUM").sum()),
            "multiplicity_mismatch_refused": n_amb,
            "no_nist_row": int((m.verdict == "NO_NIST_ROW").sum()),
            "protected_not_overwritten": int(protected.sum()),
            "adopted": int(adopt.sum()),
            "match_rate": round(hit, 4),
            "randomised_null_rate": round(null, 4) if np.isfinite(null) else None,
            "grades": json.dumps({k: int(v) for k, v
                                  in m.loc[adopt, "nist_grade"].value_counts().items()}),
            "sigma_median_dex": round(float(sig.median()), 4) if len(sig) else None,
            "n_sigma_wider_than_K07": int((sig > K07_SYSTEMATIC_DEX).sum()),
            "delta_nist_minus_store_median": round(float(m.loc[adopt, "delta_nist_minus_store"].median()), 4) if adopt.any() else None,
            "delta_nist_minus_store_absmax": round(float(m.loc[adopt, "delta_nist_minus_store"].abs().max()), 4) if adopt.any() else None,
            "n_moved_gt_0p1_dex": int((m.loc[adopt, "delta_nist_minus_store"].abs() > 0.1).sum()) if adopt.any() else 0,
        })
        all_matches.append(m)
        if np.isfinite(null) and hit <= null:
            raise SystemExit(
                f"STOP ({sp}): match rate {hit:.3f} is at or below the randomised null "
                f"{null:.3f} — this matcher is finding coincidences, not lines.")

    md = pd.concat(all_matches, ignore_index=True)
    md.to_csv(OUT / "cno_gf_adjudication.csv", index=False)
    ps = pd.DataFrame(per_species)
    ps.to_csv(OUT / "cno_gf_adjudication_by_species.csv", index=False)

    q = (pd.concat(quarantined, ignore_index=True) if quarantined
         else pd.DataFrame(columns=["species"]))
    q.to_csv(OUT / "nist_pull_nonphysical_rows.csv", index=False)

    # ------------------------------------------------------------------ O I: state it
    tach = pd.read_csv(TACHIEV)
    oi = tach[(tach.species == "O I") & tach.tachiev_log_gf.notna()].copy()
    store = canon[canon.species == "O I"]
    rows = []
    for _, t in oi.iterrows():
        s = store[(store.wavelength_air_A - float(t.wavelength_A)).abs() <= WTOL_A]
        rows.append({
            "wavelength_A": t.wavelength_A,
            "lower": t.lower, "upper": t.upper,
            "nist_log_gf": t["nist_log_gf"], "nist_grade": t.nist_grade,
            "tachiev_log_gf": t.tachiev_log_gf,
            "tachiev_rel_unc_dex": t.tachiev_rel_unc_dex,
            "tachiev_minus_nist_dex": t.delta_dex,
            "store_line_id": s.line_id.iloc[0] if len(s) == 1 else "",
            "store_log_gf": float(s.log_gf.iloc[0]) if len(s) == 1 else np.nan,
            "adopted_here": "NIST (this script) — Tachiev NOT adopted",
            "decision": "OPEN — referred to Ryan (RYA-1214 EOS)",
        })
    pd.DataFrame(rows).to_csv(OUT / "oi_tachiev_vs_nist.csv", index=False)

    # ------------------------------- what the PROTECTED rows look like against this pull
    #
    # 🔴 A PROTECTED ROW IS NOT A VERIFIED ROW. `PROTECTED_STATUS` stops this script
    # overwriting a decision somebody already took; it says nothing about whether the
    # current NIST ASD still agrees with it. Writing them out is how the one that does
    # NOT agree becomes visible instead of being protected into silence.
    #
    # It is [O I] 6300.304 — the single most important oxygen indicator in the programme.
    # The store carries -9.717 at grade **A**, cited "NIST ASD v5.11 grade A (Storey &
    # Zeippen 2000)" and adjudicated under RYA-367. This pull's NIST ASD rows for that
    # line are the M1 transition at -9.776 grade **B+** and an E2 partner at -12.20,
    # 2.4 dex weaker and negligible in the sum. So the value differs by 0.059 dex AND the
    # grade is one class better than the source now publishes. The sibling [O I] 6363.776
    # reproduces NIST exactly (-10.2580 vs -10.2581), which makes this line-specific
    # rather than a scale offset and rules out the easy explanation.
    #
    # Same shape as RYA-1171 (the 777 triplet carried A+ where ASD publishes A), one line
    # over, and nobody re-checked 6300 because RYA-1160's control had five rows. REPORTED,
    # NOT CHANGED: which source AGSS21 used for [O i] is the same open question as the
    # 777 triplet's, and resolving the dominant oxygen indicator's scale as a side effect
    # of a bulk ingest is exactly what this ticket's firewall forbids.
    prot = md[md.matched & md.adjudication_status_before.isin(PROTECTED_STATUS)].copy()
    if len(prot):
        prot["agrees_within_0p02_dex"] = prot.delta_nist_minus_store.abs() <= 0.02
        prot["action"] = "NOT OVERWRITTEN — prior adjudication stands; reported for Ryan"
        prot[["line_id", "species", "wavelength_air_A", "adjudication_status_before",
              "loggf_reference_before", "log_gf_before", "nist_log_gf", "nist_grade",
              "nist_tp_code", "delta_nist_minus_store", "agrees_within_0p02_dex",
              "action"]].to_csv(OUT / "protected_rows_vs_this_pull.csv", index=False)

    # ------------------------------------------------------------------------- report
    print("=== RYA-1214 Step 1 — CNO gf adjudication ===")
    print(f"  match window: |dlambda| <= {WTOL_A} A AND |dEP| <= {EPTOL_EV} eV, UNIQUE")
    print(ps[["species", "canonical_rows", "nist_rows_usable",
              "nist_rows_nonphysical_quarantined", "matched",
              "adopt_single", "adopt_summed_blend",
              "multiplicity_mismatch_refused", "protected_not_overwritten", "adopted",
              "match_rate", "randomised_null_rate", "sigma_median_dex",
              "n_sigma_wider_than_K07", "delta_nist_minus_store_absmax",
              "n_moved_gt_0p1_dex"]].to_string(index=False))
    print(f"\n  total rows to adopt: {len(edits)}")
    if len(q):
        print(f"\n  🔴 NON-PHYSICAL NIST ROWS QUARANTINED: {len(q)} "
              f"(grades {sorted(set(q.nist_grade))}, TP codes "
              f"{sorted(set(q.ref_transition_probability.astype(str)))}) — "
              f"implied gk {q.gk_implied.min():.3g}-{q.gk_implied.max():.3g}")
    if len(prot):
        dis = prot[~prot.agrees_within_0p02_dex]
        print(f"\n  PROTECTED rows checked against this pull: {len(prot)}, "
              f"disagreeing: {len(dis)}")
        if len(dis):
            print(dis[["line_id", "species", "wavelength_air_A",
                       "adjudication_status_before", "log_gf_before", "nist_log_gf",
                       "nist_grade", "delta_nist_minus_store"]].to_string(index=False))
    print("\n  O I 777 / 8446 — Tachiev vs NIST, STATED not resolved:")
    print(pd.DataFrame(rows)[["wavelength_A", "nist_log_gf", "tachiev_log_gf",
                              "tachiev_minus_nist_dex", "nist_grade"]].to_string(index=False))

    prov = {
        "ticket": "RYA-1214",
        "step": "1 — adjudicate the RYA-1160 NIST ASD CNO pull into canonical_gf",
        "species_adjudicated": list(SPECIES_FILES),
        "species_refused": "C III-C V, N III-N V, O III-O VI — cno_gf_pedigree states no "
                           "authority for those stages (UNESTABLISHED); an unknown "
                           "pedigree is recorded as unknown, never extended from a "
                           "neighbouring species (RYA-1072)",
        "gf_tier_written": TIER,
        "never_LAB": "C/N/O have ZERO primary-laboratory gf in this repo and no lab table "
                     "in gf_grades.LAB_TABLES. Opacity Project and MCHF are CALCULATIONS. "
                     "Tiering them LAB would repeat RYA-1005.",
        "match_tolerances": {"wavelength_A": WTOL_A, "excitation_potential_eV": EPTOL_EV,
                             "basis": "read off the |d log gf| bimodality — flat at "
                                      "~0.002 dex out to 0.1 A, median 1.05 dex beyond"},
        "null_trials": a.null_trials,
        "per_species": per_species,
        "rows_adopted": len(edits),
        "nist_pull_defect": {
            "n_rows": int(len(q)),
            "what": "oscillator strengths of 8e5-5e7 in the RYA-1160 pull; log gf came "
                    "out at +6.6 to +7.8. The Ladenburg identity gi*fik = 1.4992e-16 * "
                    "lambda^2 * gk * Aki returns gk = 8.2e8-9.9e10 on them and gk = "
                    "0.89-12.02 on the other 5,080 graded rows.",
            "confined_to": "C II only, NIST accuracy class E only, transition-probability "
                           "reference code T10528 only — one upstream delivery.",
            "would_have_been_adopted": 3,
            "impact_avoided_dex": "+8.60 to +9.74 on C II 5039.639 / 5207.982 / 5717.984",
            "action": "QUARANTINED, not repaired. gk is not in the tidy pull, so the "
                      "identity bounds the row without recovering it. "
                      "nist_pull_nonphysical_rows.csv carries all of them.",
        },
        "oi_777_decision": "OPEN. Tachiev and NIST differ systematically by -0.016 dex on "
                           "the 777 triplet and -0.005 dex on 8446. AGSS21 names Tachiev "
                           "for N i only and states no O i gf source. Both values are "
                           "reported in oi_tachiev_vs_nist.csv; NIST is what the store "
                           "already carried and this script does not change it. Ryan "
                           "decides.",
        "n_i_caveat": "NIST's N I values ARE Tachiev (TP T7370). Adjudicating N I against "
                      "this pull is one source seen twice, not independent confirmation.",
        "canonical_gf_modified": bool(a.apply),
    }
    (OUT / "cno_gf_adjudication.prov.json").write_text(json.dumps(prov, indent=2))

    if a.apply:
        changed, n_data = apply_line_surgery(edits)
        print(f"\nWROTE {CANON.relative_to(ROOT)}: {changed} of {n_data} data lines changed")
        if changed != len(edits):
            print(f"  NOTE: {len(edits) - changed} targeted row(s) were already identical")
    else:
        print("\n[dry-run] re-run with --apply to write canonical_gf.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
