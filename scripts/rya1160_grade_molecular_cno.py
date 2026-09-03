#!/usr/bin/env python3
"""
RYA-1160 — grade the 7,800 molecular CNO rows against the acquired primary archives.

🔴 WHY THIS WRITES A LEDGER AND NOT canonical_gf. Those 7,800 rows are themselves the
RYA-1149 defect: molecular transitions sitting inside the ATOMIC canonical store, which
RYA-1130 exists to prevent. Stamping grades onto them in place would deepen a schema
violation that is still open. The grades are emitted as a standalone artifact for the
RYA-1130 molecular store to consume once it exists.

Sources are the RYA-1136 primary holdings, whose parsers are REUSED rather than
reimplemented: Brooke 2013 (C2), Masseron 2014 (CH), Brooke 2014 (CN), Brooke 2015 (NH),
Brooke 2016 (OH).

⚠️ FRAME. canonical_gf is AIR above 2000 A; the primary parsers yield VACUUM. Canonical
wavelengths are converted air->vac with the same Morton/IAU formula RYA-1136 uses, so
both sides are compared in one frame. Skipping this is a ~1.5 A error at 5000 A -- larger
than the match tolerance, so it would silently produce zero matches and read as "the
primary source does not cover these lines".

Matching is wavelength + lower-energy, single-match-or-refuse. No argmin (RYA-1144).
"""
from __future__ import annotations

import csv, importlib.util, json, sys, collections
from bisect import bisect_left, bisect_right
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/audit/rya1160_cno_nist_gf"
INGEST = ROOT / "scripts/ingest_cno_molecular_primary_rya1136.py"
MOL = {"C2", "CH", "CN", "NH", "OH"}
HC_EV_CM = 1.0 / 8065.544005
W_TOL, E_TOL = 0.05, 0.05      # Angstrom (vac), eV -- E is loose: zero-points differ


def air_to_vac(w: float) -> float:
    s2 = (1e4 / w) ** 2
    n = 1 + 1e-8 * (8342.13 + 2406030 / (130 - s2) + 15997 / (38.9 - s2))
    return w * n


def parse_ax() -> dict:
    """The NH and OH A-X lists RYA-1158 found acquired-but-unread.

    🔴 THESE ARE WHY NH AND OH CAME BACK 100% ABSENT ON THE FIRST PASS. Our canonical
    NH/OH rows are NUV/VIS A-X transitions; the RYA-1136 ingest reads only the X-X
    (infrared) members of the same archives, so the correct primary data was in the tree
    and never opened. Reading it here grades those rows and closes that half of RYA-1158.

    NH-A-X-linelist.csv publishes Position(angair), f-value, Eupper, Elower, A, Branch,
    v', v", J', J" -- richer identity than the X-X lists the intake does parse.
    """
    import math
    base = ROOT / "data/reference/cno_molecular_primary"
    out = collections.defaultdict(list)
    nh = base / "nh_brooke2014/NH-A-X-linelist.csv"
    if nh.exists():
        with nh.open(encoding="utf-8-sig", errors="replace") as fh:
            for row in csv.DictReader(fh):
                try:
                    wair = float(row["Position(angair)"]); f = float(row["f-value"])
                    elo = float(row["Elower"]); jl = float(row['J"'])
                except (TypeError, ValueError, KeyError):
                    continue
                if f <= 0:
                    continue
                gf = (2 * jl + 1) * f
                out["NH"].append((air_to_vac(wair), elo * HC_EV_CM, math.log10(gf),
                                  str(row.get("Branch", "")).strip(),
                                  str(nh.relative_to(ROOT))))
    oh = base / "oh_brooke2016/OH-A-X-linelist-final.csv"
    if oh.exists():
        with oh.open(encoding="utf-8-sig", errors="replace") as fh:
            hdr = None
            for row in csv.reader(fh):
                cells = [c.strip() for c in row]
                if hdr is None:
                    low = [c.lower() for c in cells]
                    if any("angair" in c or "position" in c for c in low):
                        hdr = {c.lower(): i for i, c in enumerate(cells)}
                    continue
                def g(*keys):
                    for k in keys:
                        for h, i in hdr.items():
                            if k in h and i < len(cells):
                                try: return float(cells[i])
                                except ValueError: return None
                    return None
                wair, f, elo, jl = g("angair"), g("f-value"), g("elower"), g('j"')
                if None in (wair, f, elo, jl) or f <= 0:
                    continue
                gf = (2 * jl + 1) * f
                out["OH"].append((air_to_vac(wair), elo * HC_EV_CM, math.log10(gf),
                                  "", str(oh.relative_to(ROOT))))
    for k in out:
        out[k].sort(key=lambda t: t[0])
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    spec = importlib.util.spec_from_file_location("ing1136", INGEST)
    m = importlib.util.module_from_spec(spec); sys.modules["ing1136"] = m
    spec.loader.exec_module(m)

    idx = collections.defaultdict(list)
    for tr in m.inventory():
        idx[tr.species].append((tr.wavelength_vac_A, tr.lower_energy_eV, tr.loggf,
                                tr.label, tr.source))
    for sp, extra in parse_ax().items():
        idx[sp].extend(extra)
        print(f"  A-X list read: {sp} +{len(extra)} transitions", flush=True)
    for k in idx:
        idx[k].sort(key=lambda t: t[0])

    with (ROOT / "data/linelists/canonical_gf.csv").open(newline="") as fh:
        rows = [r for r in csv.DictReader(fh) if r["species"] in MOL]

    ledger, stats = [], collections.Counter()
    for r in rows:
        sp = r["species"]
        src = idx.get(sp, [])
        wv = air_to_vac(float(r["wavelength_air_A"]))
        ep = float(r["excitation_potential_eV"] or 0)
        waves = [t[0] for t in src]
        lo, hi = bisect_left(waves, wv - W_TOL), bisect_right(waves, wv + W_TOL)
        cand = [t for t in src[lo:hi] if abs(t[1] - ep) <= E_TOL]
        rec = {"line_id": r["line_id"], "species": sp,
               "wavelength_air_A": r["wavelength_air_A"],
               "wavelength_vac_A": f"{wv:.4f}", "EP_eV": r["excitation_potential_eV"],
               "our_log_gf": r["log_gf"], "our_tier": r["gf_tier"]}
        if len(cand) != 1:
            stats["AMBIGUOUS" if cand else "ABSENT"] += 1
            rec |= {"match": f"AMBIGUOUS({len(cand)})" if cand else "ABSENT",
                    "primary_log_gf": "", "delta_dex": "", "primary_label": "",
                    "primary_source": "", "graded": "NO"}
        else:
            w2, e2, lg, lab, srcf = cand[0]
            stats["UNIQUE"] += 1
            rec |= {"match": "UNIQUE", "primary_log_gf": f"{lg:.4f}",
                    "delta_dex": f"{float(r['log_gf']) - lg:+.4f}" if r["log_gf"] else "",
                    "primary_label": lab, "primary_source": srcf, "graded": "YES"}
        ledger.append(rec)

    with (OUT / "molecular_grade_ledger.csv").open("w", newline="") as fh:
        w_ = csv.DictWriter(fh, fieldnames=list(ledger[0]), lineterminator="\n")
        w_.writeheader(); w_.writerows(ledger)

    per = collections.defaultdict(collections.Counter)
    for r in ledger:
        per[r["species"]][r["match"].split("(")[0]] += 1
    d = [abs(float(r["delta_dex"])) for r in ledger if r["delta_dex"]]
    verdict = {"ticket": "RYA-1160", "rows": len(ledger), **dict(stats),
               "per_species": {k: dict(v) for k, v in per.items()},
               "value_deltas": len(d),
               "median_abs_delta_dex": round(sorted(d)[len(d)//2], 4) if d else None,
               "max_abs_delta_dex": round(max(d), 4) if d else None,
               "canonical_gf_modified": False,
               "why_not_written": ("these 7,800 rows are the RYA-1149 defect -- molecular "
                                   "data inside the ATOMIC store. Grades are emitted for "
                                   "the RYA-1130 molecular store, not stamped in place.")}
    (OUT / "molecular_grade_verdict.json").write_text(json.dumps(verdict, indent=2) + "\n")
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()
