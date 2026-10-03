#!/usr/bin/env python3
"""RYA-1233 -- the published Si reference line sets, one file per (set, species), plus
the registry the orchestrator reads (governing process step 7).

Ryan, 2026-10-03: graded lines are the best lines, as close to Asplund's as possible --
the published sets, with the gf they used. RYA-1218 defined Si Reference Grade as
external published line-set membership and transcribed the sets; this turns them into
what derive_band_products' `--lines-from-set NAME=CSV` consumes, so the orchestrator can
build Reference pools without an agent:

  SI_AGSS21        Asplund 2021 -> Amarsi & Asplund 2017 Table 1 (VIS / red-optical),
                   the lines they USED; gf Garz 1973 + 0.097 dex, Si II 6371 lab mean.
                   Source: data/reference/si_agss21/si_agss21_lines.csv (RYA-1169).
  SI_ELGUETA2026   Elgueta et al. 2026 (A&A 710, A111) -- the Si I lines robust (Rob=Y)
                   on the Sun/G-dwarf block, Y/J/H, 9800-18000 A; atomic data VALD3
                   ("laboratory atomic data", their sec. 2). Source: RYA-1058's
                   normalised CDS table data/audit/rya1058_elgueta/normalized_lines.csv.
  SI_BERGEMANN2013 Bergemann et al. 2013 (ApJ 764, 115) Table 1, the four J-band Si I
                   lines. Source: RYA-1218 data/audit/rya1218_si_protocol/si_nir_line_pool.csv.

`match_tol_A` is the window a set line is matched into the synthesis list with: half a
unit of the source's last printed decimal, floored at RYA-1171's 0.005 A dual-key
tolerance. One species per file -- the set loader filters by wavelength only.
"""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "reference" / "line_sets"
REGISTRY = OUT / "REGISTRY.csv"
FLOOR_A = 0.005

SETS = {
    "SI_AGSS21": {"bib": "amarsi2017_si", "doi": "10.1093/mnras/stw2445",
                  "source": "Asplund 2021 -> Amarsi & Asplund 2017 (MNRAS 464, 264) Table 1",
                  "gf_basis": "the gf the authors used: Garz 1973 + 0.097 dex (Si I); "
                              "Schulz-Gulde/Blanco/Matheron mean (Si II)"},
    "SI_ELGUETA2026": {"bib": "elgueta2026", "doi": "10.1051/0004-6361/202659148",
                       "source": "Elgueta et al. 2026 (A&A 710, A111), Sun/G-dwarf Rob=Y",
                       "gf_basis": "VALD3 (2015-01-21, v968) laboratory atomic data, as the "
                                   "authors used"},
    "SI_BERGEMANN2013": {"bib": "bergemann2013_si_jband", "doi": "10.1088/0004-637X/764/2/115",
                         "source": "Bergemann et al. 2013 (ApJ 764, 115) Table 1",
                         "gf_basis": "the authors' Table 1 log gf"},
}


def _tol(text: str) -> float:
    dec = len(text.split(".")[1]) if "." in text else 0
    return max(FLOOR_A, 0.5 * 10 ** (-dec))


def _rows():
    out = []
    with (ROOT / "data/reference/si_agss21/si_agss21_lines.csv").open() as fh:
        for r in csv.DictReader(fh):
            if r["reference_status"] != "used":
                continue
            out.append(("SI_AGSS21", r["species"], r["wavelength_air_A"], r["ep_eV"],
                        r["published_loggf"], r["band"]))
    with (ROOT / "data/audit/rya1058_elgueta/normalized_lines.csv").open() as fh:
        for r in csv.DictReader(fh):
            if r["species"] == "SiI" and r["sun_Gd_rob"].strip() == "Y":
                out.append(("SI_ELGUETA2026", "Si I", r["wavelength_A"], r["lower_ep_eV"],
                            r["elgueta_loggf"], r["band"]))
    with (ROOT / "data/audit/rya1218_si_protocol/si_nir_line_pool.csv").open() as fh:
        for r in csv.DictReader(fh):
            if r["source"] == "Bergemann2013":
                out.append(("SI_BERGEMANN2013", "Si I", r["wavelength_air_A"], r["ep_eV"],
                            r["loggf"], r["band"]))
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    groups: dict = {}
    for row in _rows():
        groups.setdefault((row[0], row[1]), []).append(row)
    reg = []
    for (name, species), rows in sorted(groups.items()):
        # The MOST precise row sets the window: a csv round-trip drops trailing zeros
        # (Table 1's 7034.90 reads back as 7034.9), which would otherwise widen it 10x.
        tol = min(_tol(r[2]) for r in rows)
        fn = OUT / f"{name.lower()}_{species.replace(' ', '')}.csv"
        with fn.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["line_set", "species", "wavelength_air_A", "ep_eV", "loggf_published",
                        "band_published", "match_tol_A", "source", "doi"])
            for r in sorted(rows, key=lambda r: float(r[2])):
                w.writerow([r[0], r[1], r[2], r[3], r[4], r[5], tol,
                            SETS[name]["source"], SETS[name]["doi"]])
        el, ion = species.split()
        ws = [float(r[2]) for r in rows]
        reg.append({"element": el, "ion": ion, "set_name": name,
                    "csv": str(fn.relative_to(ROOT)), "n_lines": len(rows),
                    "lo_A": min(ws), "hi_A": max(ws), "bibliography_key": SETS[name]["bib"],
                    "doi": SETS[name]["doi"], "source": SETS[name]["source"],
                    "gf_basis": SETS[name]["gf_basis"], "ticket": "RYA-1233"})
        print(f"  {name:<17} {species:<6} {len(rows):3d} lines  {min(ws):8.1f}-{max(ws):8.1f} A"
              f"  tol {tol} -> {fn.relative_to(ROOT)}")
    with REGISTRY.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(reg[0]))
        w.writeheader()
        w.writerows(reg)
    print(f"registry: {REGISTRY.relative_to(ROOT)} ({len(reg)} sets)")


if __name__ == "__main__":
    main()
