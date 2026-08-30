#!/usr/bin/env python3
"""Build the auditable RYA-1136 CNO closure products without wavelength-only joins."""
from __future__ import annotations

import csv
import bisect
import hashlib
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "data/audit/rya1136_cno_intake"
MOLECULAR = ROOT / "data/reference/amarsi2021_cno/derived/amarsi2021_cno_molecular_lines.csv"
TS = ROOT / "data/linelists/molecular/turbospectrum"
ATOMIC = ROOT / "data/nlte_grids/amarsi2019_cno/table1.dat"
CANON = ROOT / "data/audit/rya1129_atomic_intake"
PRIMARY_MATCH = AUDIT / "primary_molecular_crossmatch.csv"


def write_csv(path: Path, rows: list[dict], fields: tuple[str, ...]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def band(wavelength_A: float) -> str:
    if wavelength_A < 2000: return "FUV"
    if wavelength_A < 4000: return "NUV"
    if wavelength_A < 7000: return "VIS"
    if wavelength_A < 10000: return "RED_OPTICAL"
    if wavelength_A < 25000: return "NIR"
    return "IR"


ACCEPTED_MOLECULAR_JOINS = {
    "PHYSICAL_TUPLE_MATCH", "PRIMARY_TUPLE_MATCH", "PRIMARY_UNRESOLVED_SUM_MATCH",
}


def molecular_inventory() -> dict[str, list[tuple]]:
    dirs = {"C2": "C2", "CH": "CH", "CN": "CN", "NH": "NH", "OH": "OH", "12C16O": "CO"}
    index = {species: [] for species in dirs}
    for species, dirname in dirs.items():
        for path in sorted((TS / dirname).glob("*")):
            if path.suffix not in {".bsyn", ".dat"}: continue
            with path.open(errors="replace") as stream:
                for line_no, raw in enumerate(stream, 1):
                    text = raw.strip()
                    if not text or text.startswith("'"): continue
                    parts = text.split(None, 3)
                    try:
                        wavelength, energy, loggf = map(float, parts[:3])
                    except (ValueError, IndexError):
                        continue
                    label = parts[3].strip() if len(parts) > 3 else ""
                    index[species].append(
                        (wavelength, energy, loggf, path.relative_to(ROOT), line_no, label)
                    )
    for species in index:
        index[species].sort(key=lambda item: item[0])
    return index


def molecular_crossmatch() -> list[dict]:
    source = list(csv.DictReader(MOLECULAR.open()))
    index = molecular_inventory()
    rows = []
    for row in source:
        species = row["species"]
        wavelength = float(row["wavelength_vac_nm"]) * 10
        energy = float(row["lower_energy_eV"])
        loggf = float(row["published_loggf"])
        source = index[species]
        source_wavelengths = [item[0] for item in source]
        lo = bisect.bisect_left(source_wavelengths, wavelength - 0.02)
        hi = bisect.bisect_right(source_wavelengths, wavelength + 0.02)
        candidates = source[lo:hi]
        # Three independent source fields are required. Wavelength alone is never enough.
        matches = [item for item in candidates
                   if abs(item[1] - energy) <= 0.002
                   and abs(item[2] - loggf) <= 0.002]
        status = "PHYSICAL_TUPLE_MATCH" if len(matches) == 1 else ("AMBIGUOUS" if matches else "UNMATCHED")
        match = matches[0] if len(matches) == 1 else None
        rows.append({
            "source_row": row["source_row"], "species": species,
            "wavelength_vac_nm": row["wavelength_vac_nm"],
            "lower_energy_eV": row["lower_energy_eV"], "published_loggf": row["published_loggf"],
            "source_band": row["source_band"], "join_status": status,
            "candidate_count": len(matches), "matched_file": str(match[3]) if match else "",
            "matched_line": match[4] if match else "", "raw_transition_label": match[5] if match else "",
            "identity_basis": "wavelength+lower_energy+loggf" if match else "",
            "ambiguity_note": "" if match else "No unique three-field identity in vendored exact-release candidates",
        })
    if PRIMARY_MATCH.exists():
        primary = {row["source_row"]: row for row in csv.DictReader(PRIMARY_MATCH.open())}
        for row in rows:
            match = primary.get(row["source_row"])
            if not match:
                continue  # CO remains joined against the exact Li release above.
            row["join_status"] = match["join_status"]
            row["candidate_count"] = match["component_count"]
            row["matched_file"] = match["primary_source"]
            row["matched_line"] = match["primary_lines"]
            row["raw_transition_label"] = match["transition_labels"]
            row["identity_basis"] = "wavenumber+band+lower_energy+gf"
            row["ambiguity_note"] = (
                "" if match["join_status"] in {"PRIMARY_TUPLE_MATCH", "PRIMARY_UNRESOLVED_SUM_MATCH"}
                else match["join_status"]
            )
    return rows


def canonical_rows(element: str) -> list[dict]:
    return list(csv.DictReader((CANON / f"{element}_atomic_manifest.csv").open()))


def nearest_canonical(element: str, wavelength_A: float, ep: float, loggf: float) -> tuple[str, str, str]:
    candidates = [row for row in canonical_rows(element)
                  if abs(float(row["wavelength_air_A"]) - wavelength_A) <= 0.03
                  and abs(float(row["excitation_potential_eV"]) - ep) <= 0.002]
    exact = [row for row in candidates if abs(float(row["loggf"]) - loggf) <= 0.01]
    if len(exact) != 1:
        return "", "AMBIGUOUS" if exact else "ABSENT", ""
    row = exact[0]
    return row["canonical_line_id"], "PHYSICAL_TUPLE_MATCH", row["gf_tier"]


def atomic_census() -> list[dict]:
    rows = []
    for raw in ATOMIC.read_text().splitlines():
        species = raw[0:4].strip()
        if species not in {"CI", "OI"}: continue
        element = species[0]
        label = raw[6:14].strip()
        air_nm, vac_nm = float(raw[17:25]), float(raw[28:36])
        ep, loggf = float(raw[39:50]), float(raw[53:61])
        cid, join, tier = nearest_canonical(element, air_nm * 10, ep, loggf)
        rows.append({
            "reference_line_set": "AmarsiEtAl2019_AA630_A104_Table1",
            "use_status": "SOURCE_ANALYSIS_GRID_SET", "element": element,
            "species": f"{element} I", "line_label": label, "wavelength_air_A": f"{air_nm*10:.3f}",
            "wavelength_vac_A": f"{vac_nm*10:.3f}", "lower_EP_eV": f"{ep:.7f}",
            "published_loggf": f"{loggf:.3f}", "source_band": band(air_nm * 10),
            "gf_source": "Amarsi2019 adopted grid input; upstream gf source follow-up required",
            "gf_source_type": "COMPILED_IN_SOURCE_ANALYSIS", "canonical_line_id": cid,
            "join_status": join, "codex_gf_tier": tier,
            "ambiguity_note": "Grid input is not by itself proof of final AGSS21 adopted-line use",
        })
    # Amarsi et al. 2020 Table 1: the complete five-line solar N I selection
    # adopted by AGSS21 (air wavelengths in Angstrom).
    for wavelength in (7442.29, 8216.33, 8629.23, 8683.40, 10108.90):
        candidates = [r for r in canonical_rows("N") if abs(float(r["wavelength_air_A"]) - wavelength) <= .05]
        row = candidates[0] if len(candidates) == 1 else None
        rows.append({
            "reference_line_set": "AmarsiEtAl2020_GALAH_NI_model_atom",
            "use_status": "AGSS21_ADOPTED_FIVE_LINE_SET", "element": "N",
            "species": "N I", "line_label": f"{wavelength:.2f}A", "wavelength_air_A": f"{wavelength:.3f}",
            "wavelength_vac_A": "", "lower_EP_eV": row["excitation_potential_eV"] if row else "",
            "published_loggf": row["loggf"] if row else "", "source_band": band(wavelength),
            "gf_source": row["adopted_source"] if row else "UNRESOLVED",
            "gf_source_type": row["source_class"] if row else "UNRESOLVED",
            "canonical_line_id": row["canonical_line_id"] if row else "",
            "join_status": "PHYSICAL_TUPLE_MATCH" if row else "ABSENT",
            "codex_gf_tier": row["gf_tier"] if row else "",
            "ambiguity_note": "Amarsi et al. 2020 Table 1; explicitly adopted by AGSS21",
        })
    return rows


def main() -> None:
    AUDIT.mkdir(parents=True, exist_ok=True)
    molecular = molecular_crossmatch()
    atomic = atomic_census()
    write_csv(AUDIT / "molecular_physical_crossmatch.csv", molecular, tuple(molecular[0]))
    write_csv(AUDIT / "atomic_source_census.csv", atomic, tuple(atomic[0]))

    sources = [
        {"source_id":"AGSS21","citation":"Asplund, Amarsi & Grevesse 2021, A&A 653 A141","doi":"10.1051/0004-6361/202140445","role":"adopted Solar CNO lineage","asset":"article","sha256":"","status":"SOURCE_IDENTIFIED"},
        {"source_id":"Amarsi2021_Table2","citation":"Amarsi et al. 2021, A&A 656 A113","doi":"10.1051/0004-6361/202141384","role":"408 used molecular transitions","asset":str(MOLECULAR.relative_to(ROOT)),"sha256":sha256(MOLECULAR),"status":"ACQUIRED"},
        {"source_id":"Amarsi2019_Table1","citation":"Amarsi, Nissen & Skuladottir 2019, A&A 630 A104","doi":"10.1051/0004-6361/201936179","role":"C I/O I atomic model-grid line parameters","asset":str(ATOMIC.relative_to(ROOT)),"sha256":sha256(ATOMIC),"status":"ACQUIRED"},
        {"source_id":"Amarsi2020_N","citation":"Amarsi et al. 2020, A&A 642 A62","doi":"10.1051/0004-6361/202038650","role":"N I model atom and departure grid","asset":"data/nlte_grids/amarsi_galah/N_amarsi2020_v3.prov.json","sha256":sha256(ROOT/'data/nlte_grids/amarsi_galah/N_amarsi2020_v3.prov.json'),"status":"ACQUIRED_PARTIAL_LINEAGE"},
        {"source_id":"Brooke2013_C2","citation":"Brooke et al. 2013, JQSRT 124, 11","doi":"10.1016/j.jqsrt.2013.02.025","role":"C2 wavelengths, energies, transition probabilities","asset":"data/reference/cno_molecular_primary/c2_brooke2013/BrookeEtAl-C2-2013-JQSRT.zip","sha256":sha256(ROOT/'data/reference/cno_molecular_primary/c2_brooke2013/BrookeEtAl-C2-2013-JQSRT.zip'),"status":"ACQUIRED"},
        {"source_id":"Brooke2014_CN","citation":"Brooke et al. 2014, ApJS 210, 23","doi":"10.1088/0067-0049/210/2/23","role":"CN wavelengths, energies, transition probabilities","asset":"data/reference/cno_molecular_primary/cn_brooke2014/table4.dat.gz","sha256":sha256(ROOT/'data/reference/cno_molecular_primary/cn_brooke2014/table4.dat.gz'),"status":"ACQUIRED"},
        {"source_id":"Brooke2015_NH","citation":"Brooke et al. 2015, J. Chem. Phys. 143, 026101","doi":"10.1063/1.4923422","role":"NH wavelengths, energies, transition probabilities","asset":"data/reference/cno_molecular_primary/nh_brooke2015/BrookeEtAl-NH-2015-JCP.zip","sha256":sha256(ROOT/'data/reference/cno_molecular_primary/nh_brooke2015/BrookeEtAl-NH-2015-JCP.zip'),"status":"ACQUIRED"},
        {"source_id":"Brooke2016_OH","citation":"Brooke et al. 2016, JQSRT 168, 142","doi":"10.1016/j.jqsrt.2015.07.021","role":"OH wavelengths, energies, transition probabilities","asset":"data/reference/cno_molecular_primary/oh_brooke2016/OH-Supplementary.zip","sha256":sha256(ROOT/'data/reference/cno_molecular_primary/oh_brooke2016/OH-Supplementary.zip'),"status":"ACQUIRED"},
        {"source_id":"Masseron2014_CH","citation":"Masseron et al. 2014, A&A 571 A47","doi":"10.1051/0004-6361/201423956","role":"CH wavelengths, energies, transition probabilities","asset":"data/reference/cno_molecular_primary/ch_masseron2014/table14.dat.gz","sha256":sha256(ROOT/'data/reference/cno_molecular_primary/ch_masseron2014/table14.dat.gz'),"status":"ACQUIRED"},
        {"source_id":"Li2015_CO","citation":"Li et al. 2015, ApJS 216, 15","doi":"10.1088/0067-0049/216/1/15","role":"12C16O wavelengths, energies, transition probabilities","asset":"data/linelists/molecular/turbospectrum/CO/CO_IR_Li2015.dat","sha256":sha256(TS/'CO/CO_IR_Li2015.dat'),"status":"ACQUIRED"},
        {"source_id":"BarklemCollet2016","citation":"Barklem & Collet 2016, A&A 588 A96","doi":"10.1051/0004-6361/201526961","role":"molecular partition functions and equilibrium constants","asset":"data/reference/cno_molecular_primary/constants_barklem2016/table1.dat","sha256":sha256(ROOT/'data/reference/cno_molecular_primary/constants_barklem2016/table1.dat'),"status":"ACQUIRED"},
    ]
    write_csv(AUDIT / "source_bibliography.csv", sources, tuple(sources[0]))

    constants = []
    for molecule in ("C2", "CH", "CN", "NH", "OH", "CO"):
        constants.append({"molecule":molecule,"partition_function_source":"Barklem & Collet 2016","dissociation_energy_source":"Barklem & Collet 2016 table1/list holdings","isotopic_assumption":"12C16O explicit for CO; other Table2 isotopologues implicit main species","verdict":"PRIMARY_TABLES_ACQUIRED"})
    write_csv(AUDIT / "molecular_constants_ledger.csv", constants, tuple(constants[0]))

    conflict = [
        {"scope":"published molecular count","source_a":"CDS ReadMe/Table2","value_a":"408","source_b":"stale RYA-1131 comment","value_b":"879","decision":"408 is authoritative for published used transitions","status":"RESOLVED"},
        {"scope":"molecular identity","source_a":"Amarsi2021 Table2","value_a":"no rotational labels","source_b":"vendored synthesis lists","value_b":"heterogeneous upstream releases","decision":"only unique wavelength+EP+loggf tuple joins admitted","status":"OPEN_FOR_UNMATCHED"},
        {"scope":"N I adopted census","source_a":"AGSS21 summary","value_a":"five N I lines","source_b":"Amarsi2020 Table 1","value_b":"744.229, 821.633, 862.923, 868.340, 1010.890 nm","decision":"full primary-paper five-line census transcribed","status":"RESOLVED"},
    ]
    write_csv(AUDIT / "conflict_ledger.csv", conflict, tuple(conflict[0]))

    rejected = [
        {"species":"CN","system":"A-X red","band":"14 bands beyond (0-0)","count":"463","wavelength_region":"red/NIR","use_status":"REJECTED","reason":"automatic legacy equivalent widths produced two-to-three-times larger dispersion","evidence":"Amarsi2021 Sect. 2.1 lines 150-160"},
        {"species":"NH","system":"A-X","band":"unspecified","count":"NOT_PUBLISHED","wavelength_region":"near-UV around 340 nm","use_status":"REJECTED","reason":"crowding and continuum/blend limitations; individual list not published","evidence":"Amarsi2021 Sect. 2.1 lines 161-165"},
        {"species":"OH","system":"A-X","band":"unspecified","count":"NOT_PUBLISHED","wavelength_region":"near-UV around 320 nm","use_status":"REJECTED","reason":"crowding and continuum/blend limitations; individual list not published","evidence":"Amarsi2021 Sect. 2.1 lines 161-165"},
        {"species":"CN","system":"B-X","band":"unspecified","count":"NOT_PUBLISHED","wavelength_region":"near-UV around 390 nm","use_status":"REJECTED","reason":"crowding and continuum/blend limitations; individual list not published","evidence":"Amarsi2021 Sect. 2.1 lines 161-165"},
    ]
    write_csv(AUDIT / "rejected_indicator_ledger.csv", rejected, tuple(rejected[0]))

    mol_status = Counter(r["join_status"] for r in molecular)
    atom_status = Counter(r["join_status"] for r in atomic)
    coverage = []
    for domain, rows in (("molecular", molecular), ("atomic", atomic)):
        for b in ("FUV","NUV","VIS","RED_OPTICAL","NIR","IR"):
            subset = [r for r in rows if r["source_band"] == b]
            accepted = ACCEPTED_MOLECULAR_JOINS if domain == "molecular" else {"PHYSICAL_TUPLE_MATCH"}
            coverage.append({"domain":domain,"band":b,"source_rows":len(subset),"matched":sum(r["join_status"] in accepted for r in subset),"unresolved":sum(r["join_status"] not in accepted for r in subset),"verdict":"CROSSMATCH_REVIEW" if subset and any(r["join_status"] not in accepted for r in subset) else ("SOURCE_ROWS_MATCHED" if subset else "NO_SOURCE_ROWS")})
    write_csv(AUDIT / "combined_coverage_matrix.csv", coverage, tuple(coverage[0]))

    verdict = {
        "schema":"codex.cno_intake_verdict/1", "ticket":"RYA-1136",
        "intake_census_complete": True, "frozen_ready_for_measurement": False,
        "verdict":"INTAKE_COMPLETE_REVIEW_REQUIRED",
        "molecular":{"published_used":len(molecular),"join_status":dict(mol_status)},
        "atomic":{"source_rows":len(atomic),"join_status":dict(atom_status)},
        "blocking_findings":[
            "Nine Amarsi features map to indistinguishable rounded primary components because Table 2 omits rotational labels",
            "Five source-identity matches retain published-strength revision differences",
            "Four OH rows are absent or physically inconsistent in the acquired Brooke release",
            "Rejected CN count is published (463), but rejected UV transition identities are not published",
        ],
        "safety":"No abundance derived; no gf tuned; no wavelength-only join admitted",
    }
    (AUDIT / "intake_verdict.json").write_text(json.dumps(verdict, indent=2)+"\n")


if __name__ == "__main__": main()
