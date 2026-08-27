#!/usr/bin/env python3
"""Normalize and audit the immutable Chiappino et al. (2026) CRIRES+ line table."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/reference/chiappino2026/raw"
OUT = ROOT / "data/audit/rya1059_chiappino"
CANONICAL = ROOT / "data/linelists/canonical_gf.csv"
DOI = "10.3847/1538-4357/ae7de8"
ARXIV = "2606.11329"
EXPECTED = {
    "apjae7de8t1_ascii.txt": "5c3e97e8450798d6b642ace8ac9a8bbe9e78f163acd5fb6519e33cb3c382170d",
    "apjae7de8t2_mrt.txt": "3b6c04aa5762543d5b3b8ffd929cb660aa9af30b573cf3be0578196565414c3e",
    "apjae7de8t3_ascii.txt": "b176b81c892e37dae5764a8e37748fc54253d4a8aee2d07711fe4d7eca67f1e3",
    "apjae7de8t4_ascii.txt": "06963750ec85b6f4e4baaecf6ed7432209fe2c31efc58775afa7ee3010429b26",
}
MOLECULES = {"C12O", "C13O", "OH", "CN"}
ODD_HFS = {"Al", "Sc", "V", "Co", "Cu"}

sys.path.insert(0, str(ROOT))
from pipeline.line_match import match  # noqa: E402


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def band(wave: float) -> str:
    if 11160 <= wave <= 13560: return "J"
    if 14840 <= wave <= 18540: return "H"
    if 19210 <= wave <= 24720: return "K"
    return "OUTSIDE_PAPER_SETUPS"


def parse_source() -> list[dict[str, object]]:
    rows = []
    started = False
    for raw in (RAW / "apjae7de8t2_mrt.txt").read_text().splitlines():
        if raw.startswith("Fe I  "):
            started = True
        if not started or not raw.strip():
            continue
        m = re.match(r"^(\S+(?:\s+[IV]+)?)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s*$", raw)
        if not m:
            raise ValueError(f"unparsed publisher line: {raw!r}")
        label, w, ep, gf = m.groups()
        molecular = label.replace(" ", "") in MOLECULES
        if molecular:
            species, element, ion = label.replace(" ", ""), "CNO", "molecule"
        else:
            species = label
            element, ion = label.split()
        rows.append({"source_species": label, "species": species, "element": element,
                     "ion_or_molecule": ion, "atomic_or_molecular": "molecular" if molecular else "atomic",
                     "wavelength_air_A": float(w), "wavelength_system": "air",
                     "excitation_potential_eV": float(ep), "loggf_used": float(gf),
                     "band": band(float(w)), "crires_setup": {"J": "J1226", "H": "H1582", "K": "K2166"}.get(band(float(w)), ""),
                     "stated_gf_source": "Plez online compilation" if molecular else "VALD3",
                     "traced_primary_source": "UNTRACED_MOLECULAR_TRANSITION" if molecular else "UNTRACED_BEYOND_VALD3",
                     "gf_grade": "MOLECULAR_SOURCE_UNTRACED" if molecular else "VALD3_FALLBACK",
                     "gf_sigma": "", "hfs_isotope_flag": "isotope-sensitive" if species in {"C12O", "C13O"} else ("HFS_REQUIRED" if element in ODD_HFS else ""),
                     "diagnostic_role": {"C12O": "C abundance; 12C/13C", "C13O": "12C/13C", "OH": "O abundance", "CN": "N abundance"}.get(species, "atomic abundance line"),
                     "source_ref": DOI})
    return rows


def canonical_index() -> dict[str, list[dict[str, str]]]:
    with CANONICAL.open() as h:
        rows = list(csv.DictReader(h))
    out: dict[str, list[dict[str, str]]] = defaultdict(list)
    for r in rows:
        if r.get("species") and r.get("wavelength_air_A") and r.get("excitation_potential_eV"):
            out[r["species"]].append(r)
    return out


def reconcile(rows: list[dict[str, object]]) -> None:
    can = canonical_index()
    for row in rows:
        candidates = can.get(str(row["species"]), []) if row["atomic_or_molecular"] == "atomic" else []
        if candidates:
            result = match([float(row["wavelength_air_A"])], [float(x["wavelength_air_A"]) for x in candidates],
                           want_ep=[float(row["excitation_potential_eV"])],
                           src_ep=[float(x["excitation_potential_eV"]) for x in candidates],
                           tol_A=0.05, ep_tol_eV=0.02)
            unique = result.n_resolved == 1 and not result.ambiguous
            hit = candidates[int(result.index[0])] if unique else {}
            status = "PHYSICAL_KEY_UNIQUE" if unique else ("AMBIGUOUS_PHYSICAL_KEY" if result.ambiguous else "NO_PHYSICAL_KEY_MATCH")
        else:
            hit, status = {}, "MOLECULAR_SCHEMA_SEPARATE" if row["atomic_or_molecular"] == "molecular" else "NO_PHYSICAL_KEY_MATCH"
        row.update({"canonical_match_status": status, "canonical_line_id": hit.get("line_id", ""),
                    "canonical_loggf": hit.get("log_gf", ""), "canonical_gf_tier": hit.get("gf_tier", ""),
                    "canonical_gf_source": hit.get("loggf_reference", ""),
                    "delta_loggf": round(float(row["loggf_used"]) - float(hit["log_gf"]), 6) if hit else "",
                    "promotion_allowed": False,
                    "action": "TRACE_MOLECULAR_LINE_LIST" if row["atomic_or_molecular"] == "molecular" else "TRACE_PRIMARY_GF_BEFORE_PROMOTION"})


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields); w.writeheader(); w.writerows(rows)


def build(out: Path) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    manifest = []
    for name, expected in EXPECTED.items():
        actual = sha256(RAW / name)
        if actual != expected:
            raise SystemExit(f"SHA-256 mismatch for {name}: {actual} != {expected}")
        manifest.append({"filename": name, "sha256": actual, "bytes": (RAW / name).stat().st_size,
                         "source_url": "https://iopscience.iop.org/article/10.3847/1538-4357/ae7de8/data",
                         "retrieved_utc": "2026-08-27", "license_note": "publisher supplementary data; preserve unchanged"})
    write_csv(RAW.parent / "SHA256SUMS.csv", manifest)
    (RAW.parent / "metadata.json").write_text(json.dumps({"title": "CRIRES+ Reveals the Chemistry of the Stellar Subpopulations in the Bulge Fossil Fragment Liller 1", "doi": DOI, "arxiv": ARXIV, "journal": "ApJ 1006, 32", "year": 2026, "source_table": "apjae7de8t2_mrt.txt", "raw_policy": "immutable"}, indent=2) + "\n")

    rows = parse_source(); reconcile(rows)
    write_csv(out / "normalized_lines.csv", rows)
    grouped = []
    for species, b in sorted({(str(r["species"]), str(r["band"])) for r in rows}):
        g = [r for r in rows if r["species"] == species and r["band"] == b]
        grouped.append({"species": species, "band": b, "chiappino_lines": len(g),
                        "canonical_unique_matches": sum(r["canonical_match_status"] == "PHYSICAL_KEY_UNIQUE" for r in g),
                        "new_candidate_transitions": sum(r["canonical_match_status"] == "NO_PHYSICAL_KEY_MATCH" for r in g),
                        "primary_lab_gradeable_now": sum(r["canonical_gf_tier"] == "LAB" for r in g),
                        "weaker_source_only": sum(r["gf_grade"] == "VALD3_FALLBACK" for r in g),
                        "current_crires_holding_path": "alpha_cen_a_crires_plus;alpha_cen_b_crires_plus" if b in {"J", "H", "K"} else "none",
                        "telluric_risk": "REQUIRES_PER_SPECTRUM_VERIFICATION"})
    write_csv(out / "species_band_reverse_index.csv", grouped)

    fe = [x for x in grouped if x["species"] == "Fe I"]
    write_csv(out / "fe_jhk_delta.csv", fe)
    al = [r for r in rows if r["species"] == "Al I"]
    for r in al:
        r["al_classification"] = "ALREADY_COVERED" if r["canonical_match_status"] == "PHYSICAL_KEY_UNIQUE" else "EMPIRICALLY_ATTRACTIVE_WEAKER_GF"
        r["route_to"] = "RYA-716;RYA-1001"
    write_csv(out / "al_completeness_delta.csv", al)
    cno = [r for r in rows if r["species"] in {"C I", "C12O", "C13O", "CN", "OH"}]
    for r in cno:
        r["transferability"] = "GIANT_SPECIFIC_USABILITY_ONLY_PENDING_FGK_TEST" if r["atomic_or_molecular"] == "molecular" else "REQUIRES_FGK_SYNTHESIS_TEST"
        r["route_to"] = "RYA-719;RYA-720;RYA-721"
    write_csv(out / "cno_diagnostic_source_map.csv", cno)
    rejected = [r for r in rows if r["traced_primary_source"].startswith("UNTRACED")]
    write_csv(out / "empirical_lines_failing_provenance.csv", rejected)

    memo = """# Chiappino et al. 2026 methodology memo

- CRIRES+ settings: J1226 (1116–1356 nm), H1582 (1484–1854 nm), K2166 (1921–2472 nm).
- Observations: 0.4 arcsec slit, R approximately 50,000, final S/N per resolution element at least 40.
- Reduction: CR2RES v1.4.1 (dark/flat correction, nod-pair sky subtraction, arc wavelength calibration, optimal 1D extraction).
- Spectral synthesis: TURBOSPECTRUM in LTE with MARCS atmospheres.
- Atomic transitions: VALD3, including HFS entries for odd elements; the published table does not give per-line underlying references.
- Molecules: B. Plez online compilation. C uses 12CO, N uses CN, O uses OH; 13CO is retained separately for 12C/13C.
- Continuum/photon noise: continuum-placement uncertainty is stated as 1–2%. Multi-line random uncertainty is standard deviation / sqrt(N); single-line species receive 0.10 dex.
- Systematics: perturbations of ±50–100 K Teff, ±0.2 dex log(g), ±0.3 km/s microturbulence; reported abundance responses are generally <0.15 dex.
- Transfer firewall: the Liller 1 stars are cool RGB stars. Line usability, blends, continuum, LTE behavior, and telluric exposure do not transfer automatically to FGK dwarfs.
- Telluric handling and normalization are not specified with enough operational detail in the source text to inherit; every Codex spectrum still requires its own verified correction/provenance.
"""
    (out / "methodology.md").write_text(memo)
    summary = {"ticket": "RYA-1059", "lines": len(rows), "atomic_lines": sum(r["atomic_or_molecular"] == "atomic" for r in rows), "molecular_lines": sum(r["atomic_or_molecular"] == "molecular" for r in rows), "bands": dict(Counter(str(r["band"]) for r in rows)), "species": len(set(str(r["species"]) for r in rows)), "physical_key_unique": sum(r["canonical_match_status"] == "PHYSICAL_KEY_UNIQUE" for r in rows), "provenance_blocked": len(rejected), "mutation_policy": "DISCOVERY_REFERENCE_ONLY"}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (out / "README.md").write_text("# RYA-1059 Chiappino 2026 CRIRES+ digest\n\nOfficial publisher tables are preserved under `data/reference/chiappino2026/raw/` and SHA-256 verified before every build. Derived tables treat line usage as empirical usability evidence only; no Chiappino/VALD3 value is promoted as gf authority. Reproduce with `python3 scripts/rya1059_digest_chiappino.py`.\n")
    return summary


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__); p.add_argument("--out", type=Path, default=OUT); a = p.parse_args()
    print(json.dumps(build(a.out), indent=2))


if __name__ == "__main__": main()
