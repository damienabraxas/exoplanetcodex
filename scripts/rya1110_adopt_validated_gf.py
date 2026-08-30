#!/usr/bin/env python3
"""RYA-1110 - adopt the 12 VALIDATED gf for GBS lines Jofre never published. Ryan's call.

    python3 scripts/rya1110_adopt_validated_gf.py [--check]

Jofre+2014 Tables 4/5 publish no log gf for 21 of the GBS Fe set's 159 lines. RYA-1111's
options report scored Heiter et al. 2021 -- the PUBLISHED version of the GES-v3 list Jofre
used -- against the 138 lines Jofre DID publish, per reference code. Ryan ratified adopting
the 12 whose code reproduces Jofre EXACTLY on >= 5 held-out lines (2026-08-30).

    ADOPTED   12   MRW (11) + RU (1); both codes 57/57 and 5/5 exact on the held-out set
    REFUSED    9   1 THIN (BKK, exact but only n=2) + 8 RISKY (deviations to 0.180 dex)

🔴 THE ADOPTED VALUE IS NOT WRITTEN INTO `log_gf_gbs`, AND THAT IS THE WHOLE DESIGN.
`log_gf_gbs` means "what Jofre published". These 12 lines have no such value and never
will. Writing Heiter's number there would make the file assert something false about a
publication, and no later reader could tell the twelve apart from the 138. So the adoption
lands in ITS OWN artifact with its own columns, joined at read time:

    log_gf_adopted          the value
    gf_adopted_source       whose it is
    gf_adopted_basis        WHY it was adoptable -- the measured held-out rate
    gf_adopted_ticket       who decided

⚠️ AND RYA-1110's ARTIFACT IS NOT TOUCHED. `gbs_solar_fe_rya1110.csv` stays byte-identical;
this writes a SIDECAR. That keeps "what the paper published" and "what we chose to adopt"
physically separate objects -- separately reviewable, separately revertible, and impossible
to conflate by reading one file. (RYA-1084's lesson about round-tripping a data file
applies here too: the safest edit to that CSV is none.)

🔴 THE SELECTION IS READ FROM THE REPORT, NEVER RE-DERIVED. The 12 come from
`data/results/rya1111/gbs_missing_gf_options.json` by `verdict == "ADOPTABLE"`. Recomputing
the grading here would put a second copy of the threshold in the tree, free to drift from
the one Ryan actually approved -- and the whole point of the grading is that it is the
audited basis for the decision.

⚠️ WHAT THIS DOES NOT DO. It does not touch canonical_gf, any product, or any measured
value. It does not adopt the 9 refused lines. It changes no gf that Jofre DID publish: the
138 keep their published values and are not consulted here except as the validation set the
report already scored.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

REPORT = ROOT / "data" / "results" / "rya1111" / "gbs_missing_gf_options.json"
GBS = ROOT / "data" / "linelists" / "reference_sets" / "gbs_solar_fe_rya1110.csv"
OUT = ROOT / "data" / "linelists" / "reference_sets" / "gbs_gf_adoption_rya1110.csv"

ADOPT_VERDICT = "ADOPTABLE"
DECISION = ("Ryan, 2026-08-30, on RYA-1111's options report: adopt the 12 whose reference "
            "code reproduces Jofre's published gf EXACTLY on at least 5 held-out lines; "
            "leave the 1 THIN and 8 RISKY refused.")
SOURCE = ("Heiter et al. 2021, A&A 645, A106 -- the published version of the GES line "
          "list Jofre+2014 used as GES-v3 ('Heiter et al., in prep')")

COLS = ["species", "wavelength_air_A", "elo_eV", "log_gf_adopted", "gf_adopted_source",
        "gf_adopted_reference_code", "gf_adopted_basis", "gf_adopted_ticket",
        "gf_adopted_decision"]


def rows_to_adopt() -> list[dict]:
    rep = json.loads(REPORT.read_text())
    out = []
    for r in rep["per_line"]:
        if r["verdict"] != ADOPT_VERDICT:
            continue
        s = r["code_stats"]
        out.append({
            "species": r["species"],
            "wavelength_air_A": f'{r["wavelength_air_A"]:.2f}',
            "elo_eV": f'{r["elo_eV"]:.4f}',
            "log_gf_adopted": f'{r["heiter2021_log_gf"]:.3f}',
            "gf_adopted_source": SOURCE,
            "gf_adopted_reference_code": r["reference_code"],
            "gf_adopted_basis": (
                f'reference code {r["reference_code"]} reproduces Jofre+2014\'s PUBLISHED '
                f'gf exactly on {s["exact"]}/{s["n"]} held-out lines '
                f'(max deviation {s["max_dev_dex"]} dex)'),
            "gf_adopted_ticket": "RYA-1110 (decision) / RYA-1111 (validation)",
            "gf_adopted_decision": DECISION,
        })
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="verify the committed sidecar matches the report; write nothing")
    a = ap.parse_args()

    rows = rows_to_adopt()
    rep = json.loads(REPORT.read_text())
    tally = rep["tally"]

    print(f"GBS gf ADOPTION -- {len(rows)} of {rep['per_line'].__len__()} unpublished lines")
    print(f"  report tally : {tally}")
    print(f"  adopting     : {ADOPT_VERDICT} only")
    print(f"  refusing     : {tally.get('THIN', 0)} THIN + {tally.get('RISKY', 0)} RISKY "
          f"= {tally.get('THIN', 0) + tally.get('RISKY', 0)} lines, unchanged")
    for r in rows:
        print(f"    {r['species']:6} {float(r['wavelength_air_A']):9.2f} A  "
              f"Elo {r['elo_eV']}  code {r['gf_adopted_reference_code']:6} "
              f"log gf {r['log_gf_adopted']:>7}")

    # 🔴 THE VALUES MUST BE THE SOURCE'S, NOT A ROUNDED ECHO OF THEM. Re-read Heiter's own
    # column from the GBS artifact and require an exact match, so a formatting slip in this
    # script cannot quietly change an adopted gf.
    held = {}
    for r in csv.DictReader(GBS.open()):
        held[(r["species"], round(float(r["wavelength_air_A"]), 2))] = r
    bad = []
    for r in rows:
        k = (r["species"], round(float(r["wavelength_air_A"]), 2))
        src = held.get(k)
        if src is None:
            bad.append(f'{k} is not in the GBS set at all')
            continue
        if src["log_gf_gbs"].strip():
            bad.append(f'{k} HAS a published Jofre gf -- it must not be adopted over')
        if abs(float(src["heiter2021_log_gf"]) - float(r["log_gf_adopted"])) > 5e-4:
            bad.append(f'{k} adopted {r["log_gf_adopted"]} != Heiter '
                       f'{src["heiter2021_log_gf"]}')
    if bad:
        print("\nREFUSING -- the adoption does not match its source:")
        for b in bad:
            print(f"  [FAIL] {b}")
        return 1
    print("\n  every adopted value re-read from Heiter's own column and matched;")
    print("  none of the 12 has a published Jofre gf to override")

    if a.check:
        if not OUT.exists():
            print(f"\n[FAIL] {OUT.name} missing")
            return 1
        have = list(csv.DictReader(OUT.open()))
        same = ([{k: r[k] for k in COLS} for r in have] == rows)
        print(f"\n--check: committed sidecar {'MATCHES' if same else 'DIFFERS FROM'} "
              f"the report")
        return 0 if same else 1

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS)
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {OUT.relative_to(ROOT)}  ({len(rows)} adopted)")
    print(f"⚠️ {GBS.name} is NOT modified -- the adoption is a sidecar, joined at read time")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
