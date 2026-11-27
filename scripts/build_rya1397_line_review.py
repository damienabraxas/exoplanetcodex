"""Build the Solar baseline and literature candidate comparison for RYA-1397."""

import argparse
import bisect
import csv
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline.band_policy import BandPolicyError, resolve

AUDIT = ROOT / "data/audit/rya1397"
REFERENCE = ROOT / "data/reference"
LAB = REFERENCE / "lab_gf_tables"
PRODUCTS = ROOT / "data/products/solar"
FIELDS = (
    "target species band source source_path source_row source_kind source_doi "
    "wavelength wavelength_unit wavelength_medium lower_energy_eV upper_energy_eV "
    "lower_level upper_level lower_J upper_J hfs_components source_loggf "
    "source_sigma_dex source_A_s-1 source_A_uncertainty_s-1 gauge_difference theoretical_method "
    "canonical_candidate_ids canonical_loggf baseline_presence identity_status decision note"
).split()


def rows(path):
    with path.open() as stream:
        yield from csv.DictReader(line for line in stream if not line.startswith("#"))


def write(name, fields, records):
    with (AUDIT / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def band(wavelength):
    try:
        return resolve(float(wavelength)).name
    except BandPolicyError:
        return "OUTSIDE_POLICY"


def product_ids(baseline):
    accepted = defaultdict(set)
    coverage = []
    for target, entry in baseline.items():
        path = ROOT / entry["baseline_path"]
        data = path.read_bytes()
        if hashlib.sha256(data).hexdigest() != entry["baseline_sha256"]:
            raise ValueError(f"Baseline changed: {path}")
        products = json.loads(data)["products"]
        if target in {"Fe", "Fe II"}:
            ion = "II" if target == "Fe II" else "I"
            products = [product for product in products if product["ion"] == ion]
        missing = sum(not product.get("uncertainty_indicator_ids") for product in products)
        for product in products:
            for indicator in product.get("uncertainty_indicator_ids") or ():
                accepted[target].add(indicator)
        coverage.append({
            "target": target, "product_rows": len(products),
            "rows_missing_indicator_ids": missing,
            "reported_indicator_ids": len(accepted[target]),
        })
    return accepted, coverage


def target_for(species):
    return "Fe II" if species == "Fe II" else species.strip("[]").split()[0]


def locate(index, species, wavelength, ep, tolerance):
    wavelengths, candidates = index.get(species, ((), ()))
    if not candidates or ep == "":
        return []
    wave = float(wavelength)
    energy = float(ep)
    start = bisect.bisect_left(wavelengths, wave - tolerance)
    stop = bisect.bisect_right(wavelengths, wave + tolerance)
    return [
        row for row in candidates[start:stop]
        if row["excitation_potential_eV"] and
        abs(float(row["excitation_potential_eV"]) - energy) <= 0.02
    ]


def compare(row, matches, accepted):
    ids = {entry["physical_id"] for entry in matches}
    strengths = {entry["log_gf"] for entry in matches if entry["log_gf"]}
    prior = ids & accepted[row["target"]]
    if not row["lower_energy_eV"]:
        identity, decision = "LOWER_ENERGY_MISSING", "DEFER_IDENTITY"
    elif len(ids) > 1:
        identity, decision = "AMBIGUOUS", "DEFER_IDENTITY"
    elif not ids:
        identity, decision = "NO_WAVELENGTH_ENERGY_LINK", "REVIEW_NEW_CANDIDATE"
    elif row["source_loggf"] and any(
        abs(float(value) - float(row["source_loggf"])) > 0.015 for value in strengths
    ):
        identity, decision = "PROVISIONAL_WAVELENGTH_ENERGY", "REVIEW_GF_DIFFERENCE"
    else:
        identity, decision = "PROVISIONAL_WAVELENGTH_ENERGY", "REVIEW_EXISTING_ID"
    row.update(
        canonical_candidate_ids=";".join(sorted(ids)), canonical_loggf=";".join(sorted(strengths)),
        baseline_presence="REPORTED_ID" if prior else "NO_REPORTED_ID_MATCH",
        identity_status=identity, decision=decision,
    )


def legacy_fe_rows(baseline, index, accepted):
    for product in json.loads((ROOT / baseline["Fe"]["baseline_path"]).read_bytes())["products"]:
        if (product["ion"] != "I" or product.get("uncertainty_indicator_ids")
                or product["treatment"] != "ENGINE-A-3DNLTE"):
            continue
        summary = Path(product["provenance"]["copied_to"])
        path = summary.with_name(summary.name.replace("_products.csv", "_per_line.csv"))
        pool = []
        for ordinal, entry in enumerate(rows(ROOT / path), 1):
            if entry["in_domain"] != "True" or not math.isfinite(float(entry["a_3dnlte"] or "nan")):
                continue
            if entry.get("in_aggregate_1dlte", "True") != "True":
                continue
            wavelength = entry["wavelength_air_A"]
            ep = entry["elo_eV"]
            row = {
                "target": "Fe", "species": "Fe I", "band": product["band"],
                "source": "Legacy Fe I product measurement", "source_kind": "LEGACY_BASELINE",
                "source_path": path, "source_row": ordinal,
                "wavelength": wavelength, "wavelength_unit": "Å", "wavelength_medium": "air",
                "lower_energy_eV": ep, "upper_energy_eV": entry["eup_eV"],
                "source_loggf": entry.get("loggf") or entry.get("loggf_asplund", ""),
            }
            compare(row, locate(index, "Fe I", wavelength, ep, 0.03), accepted)
            row.update(
                baseline_presence="RECOVERED_PRODUCT_POOL", decision="PRESERVE_BASELINE",
                note="Recovered using the original acceptance predicate and a verified product row count.",
            )
            pool.append(row)
        if len(pool) != product["n_lines"]:
            raise ValueError(f"Recovered pool does not match recorded product count: {path}")
        yield from pool


def reference_rows(targets):
    directory = REFERENCE / "line_sets"
    for entry in rows(directory / "REGISTRY.csv"):
        path = ROOT / entry["csv"]
        for ordinal, row in enumerate(rows(path), 1):
            species = row["species"]
            if target_for(species) not in targets or not row["wavelength_air_A"]:
                continue
            yield path, ordinal, {
                "species": species, "source": row["source"], "source_kind": "PUBLISHED_SET",
                "source_doi": row["doi"], "wavelength": row["wavelength_air_A"],
                "lower_energy_eV": row["ep_eV"], "source_loggf": row["loggf_published"],
            }, float(row["match_tol_A"] or 0.01)


def lab_rows(targets):
    for path in (LAB / "lab_gf_lines.csv", *sorted((LAB / "transcribed").glob("*.csv"))):
        for ordinal, row in enumerate(rows(path), 1):
            species = row["species"]
            if target_for(species) not in targets or not row["wavelength_air_A"]:
                continue
            yield path, ordinal, {
                "species": species, "source": row["paper"], "source_kind": "LAB_TABLE",
                "source_doi": row["doi"], "wavelength": row["wavelength_air_A"],
                "lower_energy_eV": row["ep_eV"], "upper_energy_eV": row["ep_upper_eV"],
                "lower_J": row["j_lower"], "source_loggf": row["loggf"],
                "source_sigma_dex": row["loggf_sigma_dex"],
                "source_A_s-1": f"{float(row['A_1e6_s']) * 1e6:.9g}" if row["A_1e6_s"] else "",
                "source_A_uncertainty_s-1": f"{float(row['dA_1e6_s']) * 1e6:.9g}" if row["dA_1e6_s"] else "",
            }, 0.03


def molecular_rows():
    path = REFERENCE / "amarsi2021_cno/derived/amarsi2021_cno_molecular_lines.csv"
    for row in rows(path):
        yield {
            "target": row["element_parameter"].removeprefix("logeps"), "species": row["species"],
            "band": row["source_band"], "source": row["reference_line_set"],
            "source_path": path.relative_to(ROOT), "source_row": row["source_row"],
            "source_kind": "MOLECULAR_SOLAR_INDICATOR", "source_doi": "10.1051/0004-6361/202141384",
            "wavelength": row["wavelength_vac_nm"], "wavelength_unit": "nm",
            "wavelength_medium": "vacuum", "lower_energy_eV": row["lower_energy_eV"],
            "source_loggf": row["published_loggf"], "identity_status": row["transition_identity_status"],
            "decision": "DEFER_MOLECULAR_IDENTITY",
        }


def inspected_rows(index, accepted):
    path = AUDIT / "inspected_source_rows.csv"
    for ordinal, source in enumerate(rows(path), 1):
        row = {key: value for key, value in source.items() if key in FIELDS}
        target = target_for(row["species"])
        row.update(
            target=target, source_path=path.relative_to(ROOT), source_row=ordinal,
            source_kind="INSPECTED_PAPER", source=f"{source['source']} ({source['table']})",
        )
        wave = float(row["wavelength"]) * (10 if row["wavelength_unit"] == "nm" else 1)
        if row["wavelength_medium"] == "air":
            row["band"] = band(wave)
            compare(row, locate(index, row["species"].strip("[]"), wave, row["lower_energy_eV"], 0.03), accepted)
        else:
            row.update(identity_status="VERIFY_WAVELENGTH_MEDIUM", decision="DEFER_IDENTITY")
        row["decision"] = source["disposition"] or row["decision"]
        yield row


def oxygen_rows(accepted):
    path = AUDIT / "oxygen_pereira_table2.csv"
    for ordinal, source in enumerate(rows(path), 1):
        indicator = source["baseline_candidate_id"]
        yield {
            "target": "O", "species": source["species"],
            "band": band(source["baseline_candidate_air_A"]),
            "source": "Pereira et al. 2009 (Table 2)", "source_path": path.relative_to(ROOT),
            "source_row": ordinal, "source_kind": "INSPECTED_PAPER", "source_doi": source["source_doi"],
            "wavelength": source["source_wavelength_nm"], "wavelength_unit": "nm",
            "wavelength_medium": "as_published", "lower_energy_eV": source["source_lower_excitation_eV"],
            "source_loggf": source["source_log_gf"], "canonical_candidate_ids": indicator,
            "canonical_loggf": source["baseline_log_gf"],
            "baseline_presence": "REPORTED_ID" if indicator in accepted["O"] else "NO_REPORTED_ID_MATCH",
            "identity_status": "PROVISIONAL_LEVEL_CHECK_OWED",
            "decision": source["decision"], "note": source["rationale"],
        }


def local_supplement_rows(library, index, accepted):
    acquisition = ROOT / "data/audit/rya1396"
    local_files = {}
    for name in ("aluminum_downloads.csv", "silicon_downloads.csv", "supplement_downloads.csv"):
        local_files.update({row["local_file"]: row["sha256"] for row in rows(acquisition / name)})

    def checked(relative):
        path = library / relative
        if hashlib.sha256(path.read_bytes()).hexdigest() != local_files[relative]:
            raise ValueError(f"Source differs from RYA-1396 receipt: {path}")
        return path

    relative = "Iron_Tables/Den_Hartog_2019_ApJS_243_33_table6.csv"
    for ordinal, entry in enumerate(rows(checked(relative)), 1):
        wave = entry["wavelength_air_A"]
        ep = f"{float(entry['lower_level_cm-1']) / 8065.544005:.7f}"
        row = {
            "target": "Fe II", "species": "Fe II", "band": band(wave),
            "source": "Den Hartog et al. 2019 (Table 6)", "source_path": relative,
            "source_row": ordinal, "source_kind": "ACQUIRED_LAB_TABLE",
            "source_doi": "10.3847/1538-4365/ab322e", "wavelength": wave,
            "wavelength_unit": "Å", "wavelength_medium": "air", "lower_energy_eV": ep,
            "upper_energy_eV": f"{float(entry['upper_level_cm-1']) / 8065.544005:.7f}",
            "lower_J": entry["lower_J"], "upper_J": entry["upper_J"],
            "source_loggf": entry["log_gf"], "source_sigma_dex": entry["log_gf_uncertainty"],
            "source_A_s-1": f"{float(entry['transition_probability_1e6_s-1']) * 1e6:.9g}",
            "source_A_uncertainty_s-1": f"{float(entry['transition_probability_uncertainty_1e6_s-1']) * 1e6:.9g}",
        }
        compare(row, locate(index, "Fe II", wave, ep, 0.03), accepted)
        if 3170 <= float(wave) <= 3227.99:
            row.update(decision="DEFER_SOURCE_SOLAR_EXCLUSION", note="Falls in the near-UV interval excluded by the source as strong or severely blended.")
        elif row["band"] == "OUTSIDE_POLICY":
            row["decision"] = "DEFER_OUTSIDE_BANDS"
        yield row

    relative = "Aluminum_VizieR/J_AA_621_A16/all_tables.tsv"
    for ordinal, entry in vizier_rows(checked(relative), "recno\tEl\tUpper"):
        species = {"AlI": "Al I", "AlII": "Al II"}[entry["El"]]
        yield {
            "target": "Al", "species": species, "source": "Papoulia et al. 2019 (transition table)",
            "source_kind": "THEORY_CANDIDATE", "source_path": relative, "source_row": ordinal,
            "source_doi": "10.1051/0004-6361/201833764", "wavelength": entry["lambda"],
            "wavelength_unit": "Å", "wavelength_medium": "as_published",
            "lower_level": entry["Lower"], "upper_level": entry["Upper"],
            "source_loggf": loggf(entry["gf"]), "source_A_s-1": entry["A"],
            "theoretical_method": "MCDHF/RCI",
            "identity_status": "LEVEL_JOIN_OWED", "decision": "DEFER_THEORY_REVIEW",
        }

    for table, species in (("b3", "Si I"), ("b4", "Si I"), ("b7", "Si II"), ("b8", "Si II")):
        relative = f"Silicon_VizieR/J_AA_682_A184/table{table}.tsv"
        header = "recno\tULevel" if table in {"b3", "b4"} else "ULevel\t"
        for ordinal, entry in vizier_rows(checked(relative), header):
            experimental = table == "b3"
            row = {
                "target": "Si", "species": species, "source": f"Pehlivan Rhodin et al. 2024 (Table {table.upper()})",
                "source_path": relative, "source_row": ordinal,
                "source_doi": "10.1051/0004-6361/202245686", "wavelength": entry["lambda"],
                "wavelength_unit": "nm", "wavelength_medium": "vacuum" if species == "Si II" else "as_published",
                "lower_level": entry["LLevel"], "upper_level": entry["ULevel"],
                "source_loggf": entry.get("loggfExp", "") if experimental else loggf(entry.get("gf", "")),
                "source_sigma_dex": entry.get("e_loggfExp", "") if experimental else "",
                "source_A_s-1": entry.get("Aul" if experimental else "A", ""),
                "theoretical_method": "" if experimental else "MCDHF/RCI",
                "source_kind": "ACQUIRED_LAB_TABLE" if experimental else "THEORY_CANDIDATE",
                "note": ("Experimental log gf from column loggfExp." if experimental else
                         "Original calculated A and gf columns. Rescaled alternatives remain in the source."),
                "identity_status": "LEVEL_JOIN_OWED",
                "decision": "REVIEW_LAB_LEVELS" if experimental else "DEFER_THEORY_REVIEW",
            }
            if entry.get("LEnergy"):
                row["lower_energy_eV"] = f"{float(entry['LEnergy']) / 8065.544005:.7f}"
                row["upper_energy_eV"] = f"{float(entry['UEnergy']) / 8065.544005:.7f}"
            yield row


def vanadium_mirror_rows(library, index, accepted):
    directory = Path("Vanadium_Kurucz/2301_VII")
    biemont = directory / "gf2301bgfm.lab"
    wood = directory / "gf2301wlds.lab"

    def read_levels(relative):
        for ordinal, line in enumerate((library / relative).read_text().splitlines(), 1):
            if line[18:24].strip() != "23.01":
                continue
            yield ordinal, (
                float(line[:11]), line[11:18].strip(),
                float(line[24:36]), float(line[36:41]), line[41:52].strip(),
                float(line[52:64]), float(line[64:69]), line[69:80].strip(),
            )

    wood_levels = [fields for _, fields in read_levels(wood)]
    for ordinal, fields in read_levels(biemont):
        wave_nm, gf, lower, lower_j, lower_term, upper, upper_j, upper_term = fields
        wave = f"{wave_nm * 10:.3f}"
        if any(
            lower_j == candidate[3] and upper_j == candidate[6]
            and abs(lower - candidate[2]) < 0.1 and abs(upper - candidate[5]) < 0.1
            for candidate in wood_levels
        ):
            continue
        ep = f"{lower / 8065.544005:.7f}"
        row = {
            "target": "V", "species": "V II", "band": band(wave),
            "source": "Biémont et al. 1989, Kurucz mirror",
            "source_path": biemont, "source_row": ordinal, "source_kind": "LAB_MIRROR",
            "wavelength": wave, "wavelength_unit": "Å", "wavelength_medium": "air",
            "lower_energy_eV": ep, "upper_energy_eV": f"{upper / 8065.544005:.7f}",
            "lower_level": lower_term, "upper_level": upper_term,
            "lower_J": lower_j, "upper_J": upper_j, "source_loggf": gf,
        }
        compare(row, locate(index, "V II", wave, ep, 0.03), accepted)
        row.update(
            decision="DEFER_PRIMARY_SCAN_CHECK",
            note="Unmatched in the held Wood mirror at 0.1 cm-1 per level and exact J. Verify against the Biémont primary scan.",
        )
        yield row


def cno_theory_rows():
    path = REFERENCE / "cno_atomic_primary/derived/agss21_cno_atomic_gf.csv"
    for row in rows(path):
        if row["species"] != "C I" or row["wavelength_usable"] != "YES":
            continue
        wavelength = float(row["wavelength_A"])
        if band(wavelength) == "OUTSIDE_POLICY":
            continue
        yield {
            "target": "C", "species": "C I", "source": "Li et al. 2021 (Table A2)",
            "source_kind": "THEORY_CANDIDATE", "source_path": row["source_file"],
            "source_row": row["source_line"], "source_doi": "10.1093/mnras/stab214",
            "wavelength": row["wavelength_A"], "wavelength_unit": "Å", "wavelength_medium": "as_published",
            "lower_level": row["lower_level"], "upper_level": row["upper_level"],
            "source_loggf": row["log_gf"], "source_A_s-1": row["A_s-1"],
            "gauge_difference": row["dT_gauge_difference"], "theoretical_method": "MCDHF/RCI",
            "identity_status": "LEVEL_JOIN_OWED", "decision": "DEFER_THEORY_REVIEW",
        }

    path = REFERENCE / "nitrogen_method_rya1220/li2023_gf_candidates.csv"
    for ordinal, row in enumerate(rows(path), 1):
        yield {
            "target": "N", "species": "N I", "band": band(row["wavelength_air_A"]),
            "source": "Li et al. 2023 (Table 7)", "source_kind": "THEORY_CANDIDATE",
            "source_path": path.relative_to(ROOT), "source_row": ordinal,
            "source_doi": "10.3847/1538-4365/acb705", "wavelength": row["wavelength_air_A"],
            "wavelength_unit": "Å", "wavelength_medium": "air",
            "source_loggf": row["loggf_babushkin"], "theoretical_method": "MCDHF/RCI",
            "identity_status": "LEVEL_JOIN_OWED", "decision": "DEFER_THEORY_REVIEW",
            "note": f"Babushkin log gf shown. Coulomb log gf {row['loggf_coulomb']}.",
        }


def loggf(value):
    return f"{math.log10(float(value)):.6f}" if value and float(value) > 0 else ""


def vizier_rows(path, header_prefix):
    with path.open() as stream:
        for line in stream:
            if line.startswith(header_prefix):
                fields = line.rstrip("\n").split("\t")
                break
        else:
            raise ValueError(f"Missing table header: {path}")
        next(stream)
        next(stream)
        for ordinal, line in enumerate(stream, 1):
            if line.startswith("#") or not line.strip():
                break
            values = [value.strip() for value in line.rstrip("\n").split("\t")]
            if len(values) != len(fields):
                raise ValueError(f"Malformed row {ordinal}: {path}")
            yield ordinal, dict(zip(fields, values))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-documents", type=Path, required=True)
    library = parser.parse_args().reference_documents
    baseline = {row["target"]: row for row in rows(AUDIT / "solar_product_baseline.csv")
                if row["baseline_path"]}
    accepted, coverage = product_ids(baseline)
    index = defaultdict(list)
    canonical = {}
    for row in rows(ROOT / "data/linelists/canonical_gf.csv"):
        if target_for(row["species"]) in baseline:
            canonical[row["physical_id"]] = row
            index[row["species"]].append((float(row["wavelength_air_A"]), row))
    for species, entries in index.items():
        entries.sort(key=lambda item: item[0])
        index[species] = tuple(zip(*entries))

    result = defaultdict(list)
    for target, ids in accepted.items():
        for indicator in sorted(ids):
            row = canonical.get(indicator, {})
            result[target].append({
                "target": target, "species": row.get("species", ""),
                "band": band(row["wavelength_air_A"]) if row else "",
                "source": "Committed product uncertainty indicator",
                "source_path": baseline[target]["baseline_path"], "source_kind": "BASELINE",
                "source_doi": row.get("gf_source_doi", ""),
                "wavelength": row.get("wavelength_air_A", ""), "wavelength_unit": "Å" if row else "",
                "wavelength_medium": "air" if row else "",
                "lower_energy_eV": row.get("excitation_potential_eV", ""),
                "hfs_components": row.get("hfs_n_components", ""),
                "canonical_candidate_ids": indicator, "canonical_loggf": row.get("log_gf", ""),
                "baseline_presence": "REPORTED_ID",
                "identity_status": "CANONICAL_ID" if row else "PRODUCT_ID_ONLY",
                "decision": "PRESERVE_BASELINE",
            })
    seen = set()
    for source_rows in (reference_rows(baseline), lab_rows(baseline)):
        for path, ordinal, row, tolerance in source_rows:
            key = tuple(row.values())
            if key in seen:
                continue
            seen.add(key)
            wavelength = row["wavelength"]
            row.update(
                target=target_for(row["species"]), band=band(wavelength),
                source_path=path.relative_to(ROOT), source_row=ordinal,
                wavelength_unit="Å", wavelength_medium="air",
            )
            compare(row, locate(index, row["species"], wavelength, row["lower_energy_eV"], tolerance), accepted)
            result[row["target"]].append(row)
    for source_rows in (molecular_rows(), inspected_rows(index, accepted), oxygen_rows(accepted), cno_theory_rows(),
                        legacy_fe_rows(baseline, index, accepted),
                        local_supplement_rows(library, index, accepted),
                        vanadium_mirror_rows(library, index, accepted)):
        for row in source_rows:
            result[row["target"]].append(row)

    for entry in coverage:
        target = entry["target"]
        counts = Counter(row["source_kind"] for row in result[target])
        entry.update({kind: counts[kind] for kind in ("PUBLISHED_SET", "LAB_TABLE", "MOLECULAR_SOLAR_INDICATOR", "INSPECTED_PAPER",
                                                     "ACQUIRED_LAB_TABLE", "THEORY_CANDIDATE", "LAB_MIRROR", "LEGACY_BASELINE")})
    write("line_review_coverage.csv", list(coverage[0]), coverage)
    write("solar_line_review.csv", FIELDS, (row for target in baseline for row in result[target]))


if __name__ == "__main__":
    main()
