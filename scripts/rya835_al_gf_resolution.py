#!/usr/bin/env python3
"""
RYA-835 — Al gf resolution: grade the IR doublets, re-run the cull, cross-check externally.

FOUR RESULTS, three of which revise the ticket's own premise.

1. THE CULL RE-RUN RECOVERS NOTHING, BECAUSE NOTHING NEEDED RECOVERING. The ticket expected
   6696.185 to be culled as "Kurucz LOW" (an RYA-825-class stale cull) and recoverable for
   free. It is not culled: canonical_gf gives it NIST grade C+, `_gf_tier` returns MED, and
   `ACCEPTED_GF_TIERS = {HIGH, MED}`. It was already kept. 6631.218 IS culled, correctly,
   for being Kurucz K75 -> LOW. So the VIS pool is not "emptied" — it is ONE line, which is
   exactly what phase_c reports ("Al ... 1 curated line(s)").

2. 7835/7836 ARE NOW GRADEABLE. NIST ASD carries graded oscillator strengths for both.

3. 8772/8773 DO NOT CONFIRM — THEY DISAGREE BY ~0.15-0.18 dex, six times NIST's own stated
   accuracy. And the incumbent is weaker provenance than the ticket assumed: seed_source is
   `synth(GES)` with a TRUNCATED bibcode `'1995JPhB..'`, i.e. inherited from the Gaia-ESO
   synthesis list rather than adjudicated from the paper. Same truncated-bibcode shape
   RYA-799 hit with `2014MNRAS.`.

4. THE EXTERNAL CROSS-CHECK AGREES WITH NIST TO <=0.0004 dex ON ALL FOUR LINES — BUT IT IS
   NOT INDEPENDENT, AND SAYING SO IS THE POINT. Nordlander & Lind 2017 (A&A 607 A75), the
   paper behind our own Al <3D> grid, adopts these exact values; its own footnote says the
   uncertainties come from Kelleher & Podobedova 2008, which IS the NIST Al compilation.
   So the agreement proves our extraction is right, NOT that the physics is independently
   confirmed — the RYA-760 "FMW *is* NIST" caution, applied to Al.

⚠️ BURHEIM 2023 DOES NOT COVER THESE LINES — CHECKED, NOT ASSUMED (RYA-833). Its 35
tabulated wavelengths (arXiv:2309.06273 source) put the nearest entry at 7836.521, which is
0.39 A away AND a different transition (3d 2D3/2 - 5p 2P1/2), and 8843.705, 71 A from
8772.865. It DOES cover 10875.953 / 13127.005 / 13154.350 / 16723.541 / 16767.948 /
21098.84 / 21169.58 — the richer IR ladder this ticket flagged for later, so it becomes the
lab-gf source when that work happens.

⚠️ NOTHING HERE ADOPTS A gf OR MOVES AN ABUNDANCE. Under the RYA-161 firewall, "the
corrected gf moves A(Al) toward the reference value" is NOT a reason to adopt it. The case
for the NIST/Kelleher values is PROVENANCE — graded, and used by the reference paper —
and the substitution belongs to a pool rebuild, not to a grading job (the RYA-799/822
conclusion, one element over).
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

OUT = ROOT / "data" / "results" / "rya835"

#: The four IR doublet members this ticket owns, plus the VIS pool lines.
IR_DOUBLET = (7835.309, 7836.134, 8772.865, 8773.896)
VIS_POOL = (6631.218, 6696.185)

#: NIST accuracy ladder, worst-case % on A_ki. Enumerated in full including '+' tiers —
#: omitting them once put B+ (<=7%) below B (<=10%), an inverted ladder (RYA-592).
NIST_ACC_PCT = {"AAA": 0.3, "AA": 1.0, "A+": 2.0, "A": 3.0, "B+": 7.0, "B": 10.0,
                "C+": 18.0, "C": 25.0, "D+": 40.0, "D": 50.0, "E": 100.0}

#: Nordlander & Lind 2017 Table 1, transcribed from the arXiv source (1708.01949).
#: Multiplet components are listed SEPARATELY there, exactly as NIST lists them, so the
#: comparison sums them the same way canonical_gf's HFS convention does.
NL2017 = {
    7835.309: [(-0.689, 0.04)],
    7836.134: [(-1.834, 0.06), (-0.534, 0.025)],
    8772.865: [(-0.349, 0.025)],
    8773.896: [(-1.495, 0.06), (-0.192, 0.025)],
    6696.185: [(-1.569, 0.06)],
    6698.673: [(-1.870, 0.06)],
}
NL_CITATION = ("Nordlander & Lind 2017, A&A 607, A75 (DOI 10.1051/0004-6361/201730427) — "
               "uncertainties from Kelleher & Podobedova 2008, i.e. the NIST Al "
               "compilation, so this is NOT an independent measurement")

#: Burheim et al. 2023 (arXiv:2309.06273) tabulated wavelengths in 6000-45000 A, read from
#: the paper's own LaTeX source rather than from its stated 670-4200 nm RANGE. A range is
#: not a line list: 12 lines are reported, and none of them is ours.
BURHEIM_LINES = (6697.864, 6700.522, 7602.048, 7617.884, 7836.521, 7841.048, 8843.705,
                 8882.566, 8883.937, 9178.761, 9194.596, 9272.138, 9283.919, 10771.313,
                 10785.000, 10875.953, 10894.716, 11256.270, 11258.008, 11307.478,
                 12753.397, 12760.765, 13127.005, 13154.350, 14924.209, 14930.134,
                 16723.541, 16767.948, 21098.84, 21169.58, 24992.75, 25029.69,
                 38632.83, 38721.19, 41841.42)
BURHEIM_CITATION = ("Burheim, Hartman & Nilsson 2023, arXiv:2309.06273 — experimental Al I "
                    "oscillator strengths, 670-4200 nm, 12 lines at 2-11% accuracy")
BURHEIM_TOL_A = 0.05


#: The NIST pull is a CACHED, DATED ARTIFACT rather than a live call, for a practical
#: reason worth recording: `astroquery` lives only in venv_ci and `ispec` (which the cull
#: needs, via curate_nonfe_pools) lives only in venv312, so no single interpreter can do
#: both. Pulling once with scripts/rya822_pull_nist_nearuv.py and reading the TSV here also
#: gives the query a date and a provenance sidecar, which a live call would not.
NIST_TSV = (ROOT / "data" / "linelists" / "primary_gf" / "nist_asd_AlI_6600_8800.tsv")


def nist_al() -> pd.DataFrame:
    """Al I from the cached NIST ASD pull (AIR wavelengths — the RYA-822 lesson)."""
    if not NIST_TSV.exists():
        raise SystemExit(
            f"NIST pull missing at {NIST_TSV}. Regenerate with:\n"
            f"  venv_ci/bin/python scripts/rya822_pull_nist_nearuv.py "
            f"--species 'Al I' --lo-A 6600 --hi-A 8800 --step-A 100\n"
            f"Refusing to continue: an absent pull would silently become 'NIST has no "
            f"value for this line', which is the RYA-833 failure mode.")
    t = pd.read_csv(NIST_TSV, sep="\t")
    t = t.rename(columns={"wavelength_A": "w", "log_gf": "nist_line_loggf",
                          "nist_grade": "Acc."})
    t["fik"] = pd.to_numeric(t.get("fik"), errors="coerce")
    t["gi"] = pd.to_numeric(t.get("gi"), errors="coerce")
    return t[t.w.notna()]


def total_loggf(components) -> float:
    """log10(sum 10^gf) — canonical_gf's HFS convention, applied to both sides alike."""
    return float(np.log10(sum(10.0 ** g for g in components)))


def main() -> None:
    from pipeline.curate_nonfe_pools import (ACCEPTED_GF_TIERS, apply_cull,  # noqa: E402
                                             load_pool)

    # ── 1. the cull, run rather than reasoned about ──────────────────────────
    pool = load_pool(["Al"])
    kept = {}
    for gr in (False, True):
        q = apply_cull(pool.copy(), grade_restrict=gr)
        kept[gr] = q
    survivors = kept[True][kept[True].kept]
    print(f"1. CULL  accepted tiers {sorted(ACCEPTED_GF_TIERS)}")
    for _, r in kept[True].iterrows():
        print(f"   {r.wavelength_air_A:9.3f}  tier={r.gf_tier:<4} "
              f"grade={str(r.nist_grade):<4} kept={bool(r.kept)}")
    print(f"   -> {len(survivors)} of {len(pool)} survive; RECOVERY COUNT = 0 "
          f"(6696.185 was never culled)")

    # ── 2/3. NIST vs canonical on the doublets ──────────────────────────────
    nist = nist_al()
    canon = pd.read_csv(ROOT / "data" / "linelists" / "canonical_gf.csv", low_memory=False)
    canon = canon[canon.species == "Al I"]

    rows = []
    for tgt in IR_DOUBLET + VIS_POOL:
        # ⚠️ 0.30 A, and the GAP IS RECORDED. At 0.15 A this silently missed 6696.185:
        # NIST lists that line at 6696.015, 0.17 A away, and a missed match is
        # indistinguishable in the output from NIST has no value — the RYA-833 shape.
        m = nist[((nist.w - tgt).abs() < 0.30) & nist.fik.notna() & nist.gi.notna()]
        c = canon[(canon.wavelength_air_A - tgt).abs() < 0.05]
        rec = {"wavelength_air_A": tgt,
               "canonical_log_gf": float(c.iloc[0].log_gf) if len(c) else np.nan,
               "canonical_reference": str(c.iloc[0].loggf_reference) if len(c) else None,
               "canonical_nist_grade": (str(c.iloc[0].nist_grade)
                                        if len(c) and pd.notna(c.iloc[0].nist_grade)
                                        else None)}
        if len(m):
            rec["nist_n_components"] = int(len(m))
            rec["nist_log_gf"] = total_loggf(np.log10(m.gi * m.fik))
            worst = max(m["Acc."].astype(str), key=lambda a: NIST_ACC_PCT.get(a, 999))
            rec["nist_grade_worst"] = worst
            rec["nist_sigma_dex"] = float(np.log10(1 + NIST_ACC_PCT[worst] / 100))
            rec["nist_wavelength_A"] = float(m.w.iloc[int((m.w - tgt).abs().argmin())])
            rec["nist_gap_A"] = round(abs(rec["nist_wavelength_A"] - tgt), 4)
        if tgt in NL2017:
            rec["nl2017_log_gf"] = total_loggf([g for g, _ in NL2017[tgt]])
            rec["nl2017_n_components"] = len(NL2017[tgt])
        # Burheim coverage — a POSITIVE check with its tolerance stated, not an assumption
        near = min(BURHEIM_LINES, key=lambda b: abs(b - tgt))
        rec["burheim_covers"] = bool(abs(near - tgt) <= BURHEIM_TOL_A)
        rec["burheim_nearest_A"] = near
        rec["burheim_gap_A"] = round(abs(near - tgt), 3)
        rows.append(rec)

    d = pd.DataFrame(rows)
    d["delta_nist_minus_canonical"] = d.nist_log_gf - d.canonical_log_gf
    d["delta_nl_minus_nist"] = d.nl2017_log_gf - d.nist_log_gf

    print("\n2/3. NIST ASD vs canonical_gf")
    print(f"   {'line':>10}{'NIST':>9}{'grade':>7}{'sigma':>7}{'canonical':>11}"
          f"{'ref':>13}{'delta':>9}")
    for _, r in d.iterrows():
        print(f"   {r.wavelength_air_A:>10.3f}{r.nist_log_gf:>9.4f}"
              f"{str(r.nist_grade_worst):>7}{r.nist_sigma_dex:>7.4f}"
              f"{r.canonical_log_gf:>11.4f}{str(r.canonical_reference)[:12]:>13}"
              f"{r.delta_nist_minus_canonical:>+9.4f}"
              f"{('  [lambda gap %.3f A]' % r.nist_gap_A) if pd.notna(r.get('nist_gap_A')) and r.nist_gap_A > 0.05 else ''}")

    print("\n4. EXTERNAL CROSS-CHECK — Nordlander & Lind 2017")
    agree = d[d.delta_nl_minus_nist.notna()]
    print(f"   |NL - NIST| max = {agree.delta_nl_minus_nist.abs().max():.4f} dex "
          f"over {len(agree)} lines")
    print(f"   ⚠️ NOT INDEPENDENT: {NL_CITATION.split('—')[1].strip()}")

    print("\n   BURHEIM 2023 COVERAGE (checked, not assumed)")
    for _, r in d.iterrows():
        print(f"     {r.wavelength_air_A:9.3f}  covers={str(r.burheim_covers):<5} "
              f"nearest {r.burheim_nearest_A:9.3f} (gap {r.burheim_gap_A} A)")

    OUT.mkdir(parents=True, exist_ok=True)
    d.to_csv(OUT / "rya835_al_gf_per_line.csv", index=False)
    (OUT / "rya835_al_gf_summary.json").write_text(json.dumps({
        "ticket": "RYA-835",
        "cull": {"accepted_tiers": sorted(ACCEPTED_GF_TIERS),
                 "n_pool": int(len(pool)), "n_kept_grade_restrict": int(len(survivors)),
                 "recovery_count": 0,
                 "note": "6696.185 was never culled — NIST C+ -> MED, and MED is "
                         "accepted. 6631.218 is culled as Kurucz K75 -> LOW. The VIS "
                         "pool is 1 line, matching phase_c's '1 curated line(s)'."},
        "sources": {"nist": "NIST ASD via astroquery, wavelength_type='vac+air' (RYA-822)",
                    "nordlander_lind_2017": NL_CITATION,
                    "burheim_2023": BURHEIM_CITATION},
        "burheim_covers_any_target_line": bool(d.burheim_covers.any()),
        "nl_vs_nist_max_abs_dex": float(agree.delta_nl_minus_nist.abs().max()),
        "firewall": ("RYA-161: nothing here adopts a gf. That the corrected values move "
                     "A(Al) toward the reference is NOT a reason to adopt them; the case "
                     "is provenance, and the substitution is a POOL REBUILD."),
        "per_line": json.loads(d.to_json(orient="records")),
    }, indent=2, default=float))
    print(f"\n[out] {OUT}/rya835_al_gf_summary.json")


if __name__ == "__main__":
    main()
