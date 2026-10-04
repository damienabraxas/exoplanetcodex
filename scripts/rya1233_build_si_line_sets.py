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
  SI_DESHMUKH2022  Deshmukh et al. 2022 (A&A 668, A48) Table 2, the IR Si I lines of their
                   chosen subsample (10288.94, 12390.15, 12395.83), gf Pehlivan Rhodin --
                   the IR lines Lodders 2025 cites as the modern Si gf scale. Their optical
                   subsample is NOT registered: it was measured on Pehlivan Rhodin gf, the
                   canonical optical gf is Asplund's Garz, and one line carries one gf.
                   Source: data/reference/si_deshmukh2022/table2_subsample.csv.

STEP 7 -- GRADED. A published set is a line SELECTION; whether each line is graded is a
property of the gf the synthesis will use. Every set file carries `gf_graded` (the
canonical row has a published uncertainty: a stored gf_sigma_dex -- lab, Garz, Pehlivan
Rhodin -- or a NIST accuracy grade), and `<file>_graded.csv` holds only those lines. The
orchestrator dispatches the GRADED file (`graded_csv` in the registry); a line nobody has
published an uncertainty for is never measured, so it never meets a budget it cannot pass.
(The CRIRES+ Si run of 2026-10-03 measured Elgueta's whole set on unpriced VALD/Kurucz gf
-- the step-7 gap this closes.)

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
    "SI_DESHMUKH2022": {"bib": "deshmukh2022_si", "doi": "10.1051/0004-6361/202142072",
                        "source": "Deshmukh et al. 2022 (A&A 668, A48) Table 2, IR subsample",
                        "gf_basis": "Pehlivan Rhodin (2018/2024), as the authors used"},
}


def _graded(species: str, wave: float, cg: list) -> bool:
    """Does the canonical row this line synthesises with carry a published gf uncertainty?"""
    m = [r for r in cg if r["species"] == species and abs(float(r["wavelength_air_A"]) - wave) < 0.01]
    if len(m) != 1:
        return False
    r = m[0]
    sig = r.get("gf_sigma_dex", "")
    return (sig not in ("", "nan") and float(sig) == float(sig)) or bool(r.get("nist_grade", "").strip())


def _tol(text: str) -> float:
    dec = len(text.split(".")[1]) if "." in text else 0
    return max(FLOOR_A, 0.5 * 10 ** (-dec))


def _canonical_rows() -> list:
    with (ROOT / "data/linelists/canonical_gf.csv").open(newline="") as fh:
        return [r for r in csv.DictReader(fh) if r["species"].startswith("Si ")]


def _canonical_waves() -> dict:
    """canonical line id -> its wavelength in the synthesis-facing list."""
    with (ROOT / "data/linelists/canonical_gf.csv").open(newline="") as fh:
        return {r["line_id"]: r["wavelength_air_A"] for r in csv.DictReader(fh)}


def _rows():
    """(set, species, wavelength to MATCH on, ep, published log gf, band, published wavelength).

    The match wavelength is the CANONICAL one wherever a prior ticket already matched the
    published line to its canonical id (RYA-1169 for AGSS21, RYA-1058 for Elgueta): a
    source printed to 0.01 A sits up to 0.012 A from the list (Amarsi's 6741.64 vs
    6741.628), which a printed-precision window misses. The published wavelength is kept
    beside it for provenance."""
    cw = _canonical_waves()
    out = []
    with (ROOT / "data/reference/si_agss21/si_agss21_lines.csv").open() as fh:
        for r in csv.DictReader(fh):
            if r["reference_status"] != "used":
                continue
            out.append(("SI_AGSS21", r["species"],
                        cw.get(r["canonical_line_id"], r["wavelength_air_A"]), r["ep_eV"],
                        r["published_loggf"], r["band"], r["wavelength_air_A"]))
    with (ROOT / "data/audit/rya1058_elgueta/normalized_lines.csv").open() as fh:
        for r in csv.DictReader(fh):
            if r["species"] == "SiI" and r["sun_Gd_rob"].strip() == "Y":
                out.append(("SI_ELGUETA2026", "Si I",
                            cw.get(r["canonical_line_id"], r["wavelength_A"]), r["lower_ep_eV"],
                            r["elgueta_loggf"], r["band"], r["wavelength_A"]))
    with (ROOT / "data/audit/rya1218_si_protocol/si_nir_line_pool.csv").open() as fh:
        for r in csv.DictReader(fh):
            if r["source"] == "Bergemann2013":
                out.append(("SI_BERGEMANN2013", "Si I", r["wavelength_air_A"], r["ep_eV"],
                            r["loggf"], r["band"], r["wavelength_air_A"]))
    cg = _canonical_rows()
    with (ROOT / "data/reference/si_deshmukh2022/table2_subsample.csv").open() as fh:
        for r in csv.DictReader(l for l in fh if not l.startswith("#")):
            w = float(r["wavelength_air_A"])
            if w < 9000:
                continue
            # Match to the canonical wavelength (Table 2 prints 0.01 A); EP from E_low.
            m = [c for c in cg if c["species"] == r["species"]
                 and abs(float(c["wavelength_air_A"]) - w) <= 0.006]
            assert len(m) == 1, (w, len(m))
            out.append(("SI_DESHMUKH2022", r["species"], m[0]["wavelength_air_A"],
                        f"{float(r['elow_cm']) / 8065.544:.3f}", r["loggf_new"], "NIR",
                        r["wavelength_air_A"]))
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    groups: dict = {}
    for row in _rows():
        groups.setdefault((row[0], row[1]), []).append(row)
    reg = []
    cg = _canonical_rows()
    for (name, species), rows in sorted(groups.items()):
        # The MOST precise row sets the window: a csv round-trip drops trailing zeros
        # (Table 1's 7034.90 reads back as 7034.9), which would otherwise widen it 10x.
        tol = min(_tol(r[6]) for r in rows)
        fn = OUT / f"{name.lower()}_{species.replace(' ', '')}.csv"
        gfn = OUT / f"{name.lower()}_{species.replace(' ', '')}_graded.csv"
        head = ["line_set", "species", "wavelength_air_A", "wavelength_published_A",
                "ep_eV", "loggf_published", "band_published", "match_tol_A",
                "source", "doi", "gf_graded"]
        n_graded = 0
        with fn.open("w", newline="") as fh, gfn.open("w", newline="") as gh:
            w, gw = csv.writer(fh), csv.writer(gh)
            w.writerow(head)
            gw.writerow(head)
            for r in sorted(rows, key=lambda r: float(r[2])):
                g = _graded(species, float(r[2]), cg)
                row = [r[0], r[1], r[2], r[6], r[3], r[4], r[5], tol,
                       SETS[name]["source"], SETS[name]["doi"], g]
                w.writerow(row)
                if g:
                    gw.writerow(row)
                    n_graded += 1
        el, ion = species.split()
        ws = [float(r[2]) for r in rows]
        reg.append({"element": el, "ion": ion, "set_name": name,
                    "csv": str(fn.relative_to(ROOT)), "n_lines": len(rows),
                    "graded_csv": str(gfn.relative_to(ROOT)), "n_graded": n_graded,
                    "lo_A": min(ws), "hi_A": max(ws), "bibliography_key": SETS[name]["bib"],
                    "doi": SETS[name]["doi"], "source": SETS[name]["source"],
                    "gf_basis": SETS[name]["gf_basis"], "ticket": "RYA-1233"})
        print(f"  {name:<17} {species:<6} {len(rows):3d} lines ({n_graded} graded)  {min(ws):8.1f}-{max(ws):8.1f} A"
              f"  tol {tol} -> {fn.relative_to(ROOT)}")
    with REGISTRY.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(reg[0]))
        w.writeheader()
        w.writerows(reg)
    print(f"registry: {REGISTRY.relative_to(ROOT)} ({len(reg)} sets)")


if __name__ == "__main__":
    main()
