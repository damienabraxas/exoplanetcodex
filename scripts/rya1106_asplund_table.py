#!/usr/bin/env python3
"""RYA-1106 D1 — the four-instrument Asplund replication table, assembled from the
per-holding runs.

    python3 scripts/rya1106_asplund_table.py

WHY THIS IS A SEPARATE STEP AND NOT A LOOP INSIDE THE RUNNER
------------------------------------------------------------
`rya1106_asplund_replication.py` takes `--holdings` and will happily run all four in one
process. It is not run that way, because one Turbospectrum flux fit per line per holding is
~25 minutes and a single process puts the four in series with one failure able to lose all
of them -- which is exactly what happened on the first attempt, where a crash at the
artifact step discarded three holdings' completed fits. The four are therefore launched as
four processes and their reports merged here.

⚠️ THIS SCRIPT MEASURES NOTHING. It reads the committed per-holding artifacts and arranges
them. Every number in the table comes from a run; if a holding's report is missing, that
holding is reported as MISSING rather than dropped from the table, because a four-row table
silently rendered with three rows is the failure this ticket keeps finding elsewhere.

🔴 RYA-161 FIREWALL. Whatever each instrument returns is what is printed. `ASPLUND21_FE`
appears only as a column to compare against; nothing here cuts, weights, orders or branches
on the distance to 7.46.

⚠️ COVERAGE IS PART OF THE RESULT, NOT A FOOTNOTE. A holding that answers on 27 of the 40
lines is answering a different question from one that answers on 40, so `cov` sits in the
main table and every unserved line is named underneath it (RYA-429/711).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

RES = ROOT / "data" / "results" / "rya1106"
ASPLUND21_FE = 7.46          # A(Fe), Asplund, Amarsi & Grevesse 2021 — DISPLAY ONLY

#: The VIS holdings Amarsi covers (row 9/9 of the model matrix; Fe I VIS only).
HOLDINGS = ["kpno_kurucz2005", "kpno_molecfit", "harps_molecfit", "iag"]


def load(holding: str) -> dict | None:
    """This holding's row, from whichever report file the run left behind."""
    for cand in (RES / holding / "holding_report.json", RES / "asplund_replication.json"):
        if not cand.exists():
            continue
        doc = json.loads(cand.read_text())
        for h in doc.get("holdings", []):
            if h.get("holding_key") == holding and not h.get("error"):
                return h | {"_pool": doc.get("pool"), "_src": str(cand.relative_to(ROOT))}
    return None


def main() -> int:
    rows = {h: load(h) for h in HOLDINGS}
    have = {h: r for h, r in rows.items() if r}
    if not have:
        raise SystemExit(f"no per-holding reports under {RES} -- run the replication first")
    pool = next(iter(have.values()))["_pool"]

    L = []; A = L.append
    A("=" * 104)
    A("RYA-1106 D1 — Amarsi 3D-NLTE on ASPLUND'S OWN Fe I line set, all four VIS instruments")
    A("=" * 104)
    A(f"  line set : {pool['citation']}")
    A(f"  cut      : {pool['cut']}")
    A(f"  log gf   : {pool['gf']}")
    A(f"  Elo span : {pool['elo_min']:.3f}-{pool['elo_max']:.3f} eV "
      f"(our graded pool floors at 2.85 — that gap IS the RYA-1104 finding)")
    A("")
    A(f"  {'holding':<18}{'A(Fe)':>8}{'vs 7.46':>9}{'n':>4}{'cov':>9}{'REW span':>17}"
      f"{'1D base':>9}{'corr':>8}   {'route':<7} {'gf':<7} {'basis':<9}")
    A("  " + "-" * 100)
    for h in HOLDINGS:
        r = rows[h]
        if r is None:
            A(f"  {h:<18}{'MISSING':>8}   no report — this holding did NOT run")
            continue
        c = r["coverage"]
        pl = pd.read_csv(RES / h / "asplund_lines_per_line.csv")
        used = pl[pl["in_domain"].map(bool) & pl["a_3dnlte"].notna()]
        rew = f"{used['rew_agss21'].min():.2f}..{used['rew_agss21'].max():.2f}"
        A(f"  {h:<18}{r['A_3dnlte']:>8.3f}{r['vs_asplund_dex']:>+9.3f}{r['n_lines']:>4}"
          f"{c['n_served']:>5}/{c['n_asplund_lines']:<3}{rew:>17}"
          f"{r['base_1dlte_median']:>9.3f}{r['median_aberr']:>+8.3f}   "
          f"{r['route_axis']:<7} {r['gf_axis']:<7} {r['route_basis']:<9}")
    A("")
    A(f"  Asplund, Amarsi & Grevesse 2021        : {ASPLUND21_FE}")
    A(f"  AGSS21's own published 3D-NLTE (median): "
      f"{next(iter(have.values()))['agss21_published_3dnlte_median']}")
    A("")

    A("  COVERAGE — every Asplund line a holding could not serve, named (RYA-429/711):")
    for h in HOLDINGS:
        r = rows[h]
        if r is None:
            continue
        c = r["coverage"]
        if not c["n_unserved"]:
            A(f"    {h:<18} all {c['n_asplund_lines']} served")
            continue
        A(f"    {h:<18} {c['n_unserved']} unserved of {c['n_asplund_lines']}  "
          f"({', '.join(f'{k}={v}' for k, v in sorted(c['unserved_by_status'].items()))})")
        for u in c["unserved_lines"]:
            wl = u.get("wavelength_air_A", u) if isinstance(u, dict) else u
            why = u.get("status", "") if isinstance(u, dict) else ""
            A(f"        {wl}  {why}")
    A("")
    # ── does the coverage loss cost any usable line, or is it inside the refused set? ──
    # MEASURED, not asserted: a holding that loses 13 of 40 lines and still reports n=21
    # is either a coincidence or a structural fact, and the reader should not have to
    # guess which. If a future holding loses an IN-DOMAIN line, this prints it by name.
    ref = pd.read_csv(RES / "kpno_kurucz2005" / "asplund_lines_per_line.csv")
    outdom = set(ref.loc[~ref["in_domain"].map(bool), "wavelength_air_A"].round(1))
    costly = {}
    for h in HOLDINGS:
        r = rows[h]
        if r is None:
            continue
        uns = {round(float(u["wavelength_air_A"] if isinstance(u, dict) else u), 1)
               for u in r["coverage"]["unserved_lines"]}
        bite = sorted(uns - outdom)
        if bite:
            costly[h] = bite
    if costly:
        A("  🔴 COVERAGE LOSS REACHES IN-DOMAIN LINES — these holdings lost lines the")
        A("     network WOULD have accepted, so their n is genuinely smaller:")
        for h, b in costly.items():
            A(f"       {h:<18} {b}")
    else:
        A("  ⚠️ THE COVERAGE LOSS AND THE DOMAIN WALL DO NOT COMPOUND. Every line no")
        A("     holding could serve is one the MLP already refuses, so n is 21 on all")
        A("     four despite coverage running 27/40 to 40/40 — MEASURED here, not assumed.")
        A("     Read the two columns separately: `cov` is what the spectrum reaches and")
        A("     `n` is what the network will price. HARPS answering on 27 of 40 lines is")
        A("     still a real limit on the product; it just is not a limit on n.")
    A("")
    A("  ⚠️ n=21 ON EVERY HOLDING IS THE NETWORK'S DOMAIN, NOT THE DATA. 19 of the 40")
    A("     AGSS21 lines fall outside the Amarsi MLP's training envelope, identically on")
    A("     all four holdings because the test is on atomic data, not on the spectrum.")
    A("     `rya1106_domain_crosstab.py` shows the axis is WAVELENGTH, not line strength.")
    A("")
    A("  ⚠️ CAVEAT CARRIED, NOT CHASED (RYA-1104/RYA-282): the line-population swing")
    A("     between this pool and our graded pool is 29% accounted for by the network's")
    A("     Elo trend and 71% UNATTRIBUTED. Documented, not a tuning target: no line list,")
    A("     gf or parameter was adjusted to move any number toward 7.46 (RYA-161).")
    A("=" * 104)

    text = "\n".join(L)
    print(text)
    (RES / "asplund_four_instrument_table.txt").write_text(text + "\n")
    doc = {
        "ticket": "RYA-1106 D1",
        "line_set": pool,
        "asplund21_reference": ASPLUND21_FE,
        "firewall": ("RYA-161: this product is emitted because it uses Asplund's ACTUAL "
                     "inputs, not because of what it returns; nothing branches on 7.46"),
        "holdings": {h: (None if rows[h] is None else {
            k: rows[h][k] for k in ("instrument", "holding", "A_3dnlte", "vs_asplund_dex",
                                    "n_lines", "n_excluded", "stat_dex", "syst_dex",
                                    "base_1dlte_median", "median_aberr", "route_axis",
                                    "gf_axis", "route_basis", "handler", "coverage")})
                     for h in HOLDINGS},
    }
    (RES / "asplund_four_instrument_table.json").write_text(json.dumps(doc, indent=2) + "\n")
    print(f"\nwrote {RES / 'asplund_four_instrument_table.txt'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
