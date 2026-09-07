#!/usr/bin/env python3
"""Build the auditable RYA-1136 CNO closure products without wavelength-only joins."""
from __future__ import annotations

import csv
import os
import re
import functools
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


#: The AGSS21-adopted N I five-line set: (air wavelength in Angstrom, published lower-level
#: excitation potential in eV). The EP is what makes the join physical rather than
#: wavelength-only (RYA-1143); without it these five rows carried an unearned
#: PHYSICAL_TUPLE_MATCH.
N_I_ADOPTED_SET = (
    (7442.29, 10.3301), (8216.33, 10.3359), (8629.23, 10.6900),
    (8683.40, 10.3301), (10108.90, 11.7529),
)

#: 🔴 RYA-1144. `AMBIGUOUS_SUM_MATCH` is deliberately ABSENT: a feature whose published
#: loggf is reproduced by more than one subset of primary components has not been
#: identified, only fitted, and must not be counted as coverage.
#: `PRIMARY_UNRESOLVED_SUM_MATCH` now means exactly ONE viable subset.
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
            # 🔴 RYA-1181 (A3). `system` and the VIBRATIONAL band are published by the
            # source and were parsed into primary_molecular_crossmatch, then dropped in
            # the merge to here -- the artifact everything downstream actually reads.
            # `source_band` is the REPORTING bin (VIS/NIR/IR); it is not the band.
            "system": row.get("system", ""),
            "vibrational_band": row.get("band", ""),
            "candidate_count": len(matches), "matched_file": str(match[3]) if match else "",
            "matched_line": match[4] if match else "", "raw_transition_label": match[5] if match else "",
            "identity_basis": "wavelength+lower_energy+loggf" if match else "",
            **split_transition_label(species, str(match[5]) if match else ""),
            # 🔴 RYA-1180: a match is only as good as what it matched AGAINST, and that was
            # nowhere on the row. Derived from the matched file's PATH, not hardcoded per
            # species, so a future source is graded by where it actually came from.
            "source_provenance_grade": provenance_grade(str(match[3]) if match else ""),
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
            row["source_provenance_grade"] = provenance_grade(match["primary_source"])
            # RYA-1181: the primary side publishes these and this merge used to discard them.
            row["system"] = match.get("system", row.get("system", ""))
            row["vibrational_band"] = match.get("band", row.get("vibrational_band", ""))
            row["primary_j_lower"] = match.get("primary_j_lower", "")
            row.update(split_transition_label(row["species"], match["transition_labels"]))
            row["ambiguity_note"] = (
                "" if match["join_status"] in {"PRIMARY_TUPLE_MATCH", "PRIMARY_UNRESOLVED_SUM_MATCH"}
                else match["join_status"]
            )
    return rows


def canonical_rows(element: str) -> list[dict]:
    return list(csv.DictReader((CANON / f"{element}_atomic_manifest.csv").open()))


def nearest_canonical(element: str, wavelength_A: float, ep: float,
                      loggf: float | None) -> tuple[str, str, str]:
    """Wavelength AND excitation potential, then loggf when the caller has one.

    🔴 `loggf=None` means "the published gf is not in hand". The join then rests on
    wavelength + EP, which is still a two-field physical identity and still refuses
    ambiguity -- it is NOT a wavelength-only key. It must never be called with
    loggf=None merely to make a stubborn line match.
    """
    candidates = [row for row in canonical_rows(element)
                  if abs(float(row["wavelength_air_A"]) - wavelength_A) <= 0.03
                  and abs(float(row["excitation_potential_eV"]) - ep) <= 0.002]
    exact = (candidates if loggf is None else
             [row for row in candidates if abs(float(row["loggf"]) - loggf) <= 0.01])
    if len(exact) != 1:
        return "", "AMBIGUOUS" if exact else "ABSENT", ""
    row = exact[0]
    return row["canonical_line_id"], "PHYSICAL_TUPLE_MATCH", row["gf_tier"]


#: Where a matched line actually came from. A join status says the three fields agreed;
#: it says nothing about whether the thing they agreed with is a primary table or a
#: redistribution of one -- and RYA-1142 A4 found the intake's ONLY clean-match class
#: (all 80 12C16O PHYSICAL_TUPLE_MATCH rows) resting entirely on an ExoMol redistribution
#: converted by a script external to this repo. That has to be visible on the row.
PROVENANCE_GRADES = {
    "data/reference/cno_molecular_primary/": "PRIMARY_PUBLISHED_TABLE",
    "data/linelists/molecular/turbospectrum/": "REDISTRIBUTION_VENDORED_SYNTHESIS_LIST",
}


#: RYA-854's single source of truth for every reference the project cites. RYA-1170: an
#: intake must POINT at it, never copy a value out of it -- a copied DOI drifted into a
#: different paper and stayed wrong through a full QA pass.
SSOT_BIBLIOGRAPHY = ROOT / "data" / "refs" / "bibliography.csv"

#: `local_file` is stored relative to the CODEX ROOT, one level above the reference
#: library, and that library is a working directory rather than a repo artifact -- it
#: differs per machine and is legitimately absent on CI. `generate_sources_page`
#: audit_library says so explicitly and treats absence as non-fatal; this follows it.
#: Override for a machine whose library sits elsewhere.
CODEX_LIBRARY_ROOT = Path(
    os.environ.get("CODEX_LIBRARY_ROOT", Path.home() / "Documents" / "Exoplanet Codex"))


@functools.lru_cache(maxsize=1)
def ssot_bibliography() -> dict[str, dict]:
    with SSOT_BIBLIOGRAPHY.open(newline="") as fh:
        return {r["key"]: r for r in csv.DictReader(fh)}


def ssot_local_file(key: str) -> str:
    """The SSOT's own local_file path for this key, verbatim. '' if it declares none."""
    return ssot_bibliography()[key]["local_file"].strip()


def ssot_local_sha(key: str, source_id: str = "") -> str:
    """sha256 of the held document -- COMPUTED where the library is reachable, otherwise
    CARRIED FORWARD from what this artifact already records.

    🔴 A CHECKSUM IS A RECORDED FACT, NOT A PER-RUN MEASUREMENT. The reference library
    lives above the repo and is absent on CI, and `test_cno_closure_rya1136` re-runs this
    build there. Recomputing-or-blanking would therefore erase the checksum on every CI
    run and then trip RYA-1170's own "empty sha256 on a HELD key is a broken link" rule --
    a guard defeated by the machine it runs on. So on a machine without the library the
    previously recorded value stands, unchanged, and is neither re-derived nor faked.
    """
    rel = ssot_local_file(key)
    if not rel:
        return ""
    path = CODEX_LIBRARY_ROOT / rel
    if path.is_file():
        return sha256(path)
    return _recorded_sha(source_id or key)


def _recorded_sha(source_id: str) -> str:
    """The sha256 this artifact already carries for a row, or '' if it has none yet."""
    out = AUDIT / "source_bibliography.csv"
    if not out.exists():
        return ""
    with out.open(newline="") as fh:
        for row in csv.DictReader(fh):
            if row.get("source_id") == source_id:
                return (row.get("sha256") or "").strip()
    return ""


def split_transition_label(species: str, label: str) -> dict:
    """🔴 RYA-1181 (A3) — a label column that is an UNSPLIT REMAINDER is not identity.

    Two of the intake's "transition label" columns were whatever the parser had left over:

      CO  the tail of the Turbospectrum record --
          "0.000 199.0 2.41E+01 'X' 'X' 99.0 98.0 'X' 'X' 0.0 0.0 'v1-0_J99-98_Li2015'"
          which carries J', J'', v' and v'' inside a string nothing can join on.
      CH  branch and the observed-minus-calculated residual glued together, "R  0.00274".

    Both are split into named columns here. Fields absent for a species stay EMPTY rather
    than being invented -- an absent field and a fabricated one must never look alike.
    """
    out = {"branch": "", "o_minus_c": "", "j_upper": "", "j_lower_label": "",
           "v_upper": "", "v_lower": ""}
    lab = (label or "").strip()
    if not lab:
        return out
    first = lab.split(";")[0].strip()
    if species == "12C16O":
        # the Li2015 tag is the only self-describing part: 'v1-0_J99-98_Li2015'
        m = re.search(r"'v(\d+)-(\d+)_J(\d+)-(\d+)_", first)
        if m:
            out["v_upper"], out["v_lower"] = m.group(1), m.group(2)
            out["j_upper"], out["j_lower_label"] = m.group(3), m.group(4)
        return out
    parts = first.split()
    if parts and re.fullmatch(r"[A-Za-z]{1,3}\d*", parts[0]):
        out["branch"] = parts[0]
        if len(parts) > 1:
            try:
                float(parts[1])
                out["o_minus_c"] = parts[1]
            except ValueError:
                pass
    return out


def provenance_grade(matched_file: str) -> str:
    """Grade a match by WHERE it matched, from the path. '' for an unmatched row.

    Unrecognised is UNCLASSIFIED, never a default grade: silently calling an unknown
    source primary is the defect this column exists to stop (RYA-1072 -- an allow-list
    for the true case must not launder the case it does not recognise).
    """
    if not matched_file:
        return ""
    for prefix, grade in PROVENANCE_GRADES.items():
        if matched_file.startswith(prefix):
            return grade
    return "UNCLASSIFIED_SOURCE"


def barklem_adopted_de(path: Path) -> dict[str, tuple[str, str]]:
    """{molecule: (De, e_De)} from Barklem & Collet 2016 table1.dat, read POSITIONALLY.

    The CDS byte-by-byte description defines four dissociation-energy columns; only the
    last is the paper's own adopted value:

        18- 27  HH    Huber & Herzberg 1979      (comparison)
        40- 49  Luo   Luo 2007                   (comparison)
        62- 71  G2    G2 theory                  (comparison)
        84- 93  De    "Dissociation energy adopted"   <- this one
        94-103  e_De  error in the adopted value

    Splitting on whitespace and taking a value that "looks right" is how RYA-853 captured
    another author's comparison column as the paper's own number. Byte offsets, or nothing.
    """
    out: dict[str, tuple[str, str]] = {}
    for line in path.read_text().splitlines():
        mol = line[0:5].strip()
        if not mol:
            continue
        out[mol] = (line[83:93].strip(), line[93:103].strip())
    return out


#: 🔴 RYA-1183 (RYA-1142 A6b) — THE Ni I BLEND AT [O I] 6300 IS A PHYSICAL COMPONENT.
#:
#: [O I] 6300.30 is the single most-used solar oxygen diagnostic and it is NOT usable
#: alone: Ni I 6300.34 sits inside it. `pipeline.cno_synthesis` already treats the feature
#: as `kind='forbidden_blend'` requiring a joint synthesis with A(Ni) pinned -- but the
#: atomic census had no blend column and no Ni I row, so the best-known contaminant of the
#: diagnostic was retained NOWHERE in the intake. A census that lists the [O I] line and
#: not its blend partner describes a line that does not exist in isolation.
#:
#: Every value here is read from our own store (canonical_gf line gf_101075), not typed
#: from memory, and the row is marked BLEND_COMPONENT so it can never be mistaken for an
#: AGSS21 source line.
NI_6300_BLEND_ROW = {
    "reference_line_set": "JohanssonEtAl2003_ApJ584_L107",
    "use_status": "BLEND_COMPONENT_NOT_A_SOURCE_LINE", "element": "Ni",
    "species": "Ni I", "line_label": "6300.34A",
    "wavelength_air_A": "6300.342", "wavelength_vac_A": "",
    "lower_EP_eV": "4.2660000", "published_loggf": "-2.110",
    "source_band": "RED_OPTICAL",
    "gf_source": "Johansson, Litzen, Lundberg & Zhang 2003, ApJ 584, L107 "
                 "(arXiv:astro-ph/0301382)",
    "gf_source_type": "PRIMARY_LABORATORY_BLEND_PARTNER",
    "canonical_line_id": "gf_101075", "join_status": "CARRIED_AS_BLEND_COMPONENT",
    "codex_gf_tier": "OTHER", "codex_loggf": "-2.11", "codex_EP_eV": "4.266",
    "ambiguity_note": ("Not an AGSS21 adopted line. Retained because [O I] 6300.300 cannot "
                       "be measured without it -- pipeline.cno_synthesis registers the "
                       "feature as kind='forbidden_blend', '[O I] 6300.30 + Ni I 6300.34 "
                       "joint synthesis (A(Ni) pinned)'. RYA-1183 / RYA-1142 A6b."),
    "blend_role": "CONTAMINANT_COMPONENT",
    "blend_partner_A": "6300.300",
    "blend_note": "Ni I 6300.34 contaminates [O I] 6300.30; joint synthesis required.",
}


def blend_fields(element: str, wavelength_A: float) -> dict:
    """Blend/component columns. Only [O I] 6300 is populated, because it is the only
    feature in this census the pipeline registers as a blend -- an empty column here means
    'no blend recorded', never 'checked and clean'."""
    if element == "O" and abs(float(wavelength_A) - 6300.30) <= 0.05:
        return {"blend_role": "PRIMARY_DIAGNOSTIC_IN_A_BLEND",
                "blend_partner_A": "6300.342",
                "blend_note": ("Ni I 6300.34 sits inside this feature; "
                               "pipeline.cno_synthesis requires a joint synthesis with "
                               "A(Ni) pinned. Component carried in this census (RYA-1183).")}
    return {"blend_role": "", "blend_partner_A": "", "blend_note": ""}


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
            "codex_loggf": "", "codex_EP_eV": "",
            "ambiguity_note": "Grid input is not by itself proof of final AGSS21 adopted-line use",
            **blend_fields(element, air_nm * 10),
        })
    # Amarsi et al. 2020 Table 1: the complete five-line solar N I selection adopted by
    # AGSS21 (air wavelengths in Angstrom), with the PUBLISHED excitation potential so
    # the join has a second physical field to stand on.
    #
    # 🔴 RYA-1143. This loop previously selected on `abs(wavelength - w) <= .05` and
    # nothing else, then stamped the result PHYSICAL_TUPLE_MATCH -- a wavelength-only
    # key wearing a physical-identity label, which is exactly the RYA-1034 defect. It
    # also reported `lower_EP_eV` and `published_loggf` READ OUT OF the matched
    # canonical row, so our own stored value round-tripped back under a column named
    # `published_loggf` (the RYA-1035 vendor-echo defect).
    #
    # Both are fixed here: the join is EP-aware via nearest_canonical(), and the
    # canonical readback is reported in explicitly-named codex_* columns while the
    # published_* columns carry the SOURCE's values or stay empty.
    for wavelength, published_ep in N_I_ADOPTED_SET:
        cid, join, tier = nearest_canonical("N", wavelength, published_ep, None)
        row = next((r for r in canonical_rows("N") if r["canonical_line_id"] == cid), None)
        rows.append({
            "reference_line_set": "AmarsiEtAl2020_GALAH_NI_model_atom",
            "use_status": "AGSS21_ADOPTED_FIVE_LINE_SET", "element": "N",
            "species": "N I", "line_label": f"{wavelength:.2f}A", "wavelength_air_A": f"{wavelength:.3f}",
            "wavelength_vac_A": "", "lower_EP_eV": f"{published_ep:.4f}",
            "published_loggf": "", "source_band": band(wavelength),
            "gf_source": "Amarsi et al. 2020 Table 1 (published loggf not transcribed)",
            "gf_source_type": "PUBLISHED_SOURCE_NOT_YET_TRANSCRIBED",
            "canonical_line_id": cid,
            "join_status": join,
            "codex_gf_tier": tier,
            "codex_loggf": row["loggf"] if row else "",
            "codex_EP_eV": row["excitation_potential_eV"] if row else "",
            "ambiguity_note": ("Amarsi et al. 2020 Table 1; explicitly adopted by AGSS21. "
                               "EP-aware join (RYA-1143). published_loggf is EMPTY because "
                               "the source value has not been transcribed -- the codex_* "
                               "columns are OUR store, not the paper's."),
            **blend_fields("N", wavelength),
        })
    rows.append(NI_6300_BLEND_ROW)
    return rows


def main() -> None:
    AUDIT.mkdir(parents=True, exist_ok=True)
    molecular = molecular_crossmatch()
    atomic = atomic_census()
    write_csv(AUDIT / "molecular_physical_crossmatch.csv", molecular, tuple(molecular[0]))
    write_csv(AUDIT / "atomic_source_census.csv", atomic, tuple(atomic[0]))

    # 🔴 RYA-1170 -- THE DOI IS READ FROM THE SSOT, NOT COPIED BESIDE IT.
    #
    # `Amarsi2019_Table1` carried 10.1051/0004-6361/201936179, which Crossref resolves to
    # Curran & Moss, "Quasi-stellar object redshift estimates", A&A 629 -- a different
    # paper, different authors, different volume -- as the cited source for the ENTIRE
    # C I / O I atomic census. data/refs/bibliography.csv held the right DOI
    # (10.1051/0004-6361/201936265) the whole time. That is RYA-355's single-source
    # defect: a value duplicated into a second file instead of pointed at one copy, and
    # the duplicate drifted. Note the row's own CITATION string was right ("A&A 630 A104")
    # and disagreed with its own DOI.
    #
    # So a row that names an `ssot_key` no longer states a DOI at all -- `ssot_doi()`
    # fetches it, and a divergence is now impossible rather than merely detected. Rows
    # whose source bibliography.csv does not hold keep their own DOI and say so with an
    # empty ssot_key; a guard reports them so the absence stays visible.
    def ssot_doi(key: str) -> str:
        row = ssot_bibliography()[key]
        return row["doi"].strip()

    sources = [
        {"source_id":"AGSS21","source_type":"article","citation":"Asplund, Amarsi & Grevesse 2021, A&A 653 A141","ssot_key":"asplund2021","doi":ssot_doi("asplund2021"),"role":"adopted Solar CNO lineage","asset":ssot_local_file("asplund2021"),"sha256":ssot_local_sha("asplund2021","AGSS21"),"status":"ACQUIRED_REFERENCE_LIBRARY","provenance_note":"RYA-1170 amendment: this row carried asset='article' with an EMPTY sha256 while bibliography.csv key `asplund2021` records the paper as HELD (verified=extracted). That is a broken link, not a missing source, and it is what left the A1/A7 reconciliation with no acquired referent. Now points at the SSOT's own local_file. NOTE the reference library lives one level above the repo and is a working directory, not a repo artifact (see generate_sources_page.audit_library) -- so the checksum is verifiable only where the library is present, and is recorded rather than re-derived elsewhere."},
        {"source_id":"Amarsi2021_Table2","source_type":"primary_published_table","citation":"Amarsi et al. 2021, A&A 656 A113","ssot_key":"","doi":"10.1051/0004-6361/202141384","role":"408 used molecular transitions","asset":str(MOLECULAR.relative_to(ROOT)),"sha256":sha256(MOLECULAR),"status":"ACQUIRED","provenance_note":""},
        {"source_id":"Amarsi2019_Table1","source_type":"primary_published_table","citation":"Amarsi, Nissen & Skuladottir 2019, A&A 630 A104","ssot_key":"amarsi2019","doi":ssot_doi("amarsi2019"),"role":"C I/O I atomic model-grid line parameters","asset":str(ATOMIC.relative_to(ROOT)),"sha256":sha256(ATOMIC),"status":"ACQUIRED","provenance_note":""},
        {"source_id":"Amarsi2020_N","source_type":"derived_grid","citation":"Amarsi et al. 2020, A&A 642 A62","ssot_key":"","doi":"10.1051/0004-6361/202038650","role":"N I model atom and departure grid","asset":"data/nlte_grids/amarsi_galah/N_amarsi2020_v3.prov.json","sha256":sha256(ROOT/'data/nlte_grids/amarsi_galah/N_amarsi2020_v3.prov.json'),"status":"ACQUIRED_PARTIAL_LINEAGE","provenance_note":""},
        {"source_id":"Brooke2013_C2","source_type":"primary_published_table","citation":"Brooke et al. 2013, JQSRT 124, 11","ssot_key":"","doi":"10.1016/j.jqsrt.2013.02.025","role":"C2 wavelengths, energies, transition probabilities","asset":"data/reference/cno_molecular_primary/c2_brooke2013/BrookeEtAl-C2-2013-JQSRT.zip","sha256":sha256(ROOT/'data/reference/cno_molecular_primary/c2_brooke2013/BrookeEtAl-C2-2013-JQSRT.zip'),"status":"ACQUIRED","provenance_note":""},
        {"source_id":"Brooke2014_CN","source_type":"primary_published_table","citation":"Brooke et al. 2014, ApJS 210, 23","ssot_key":"","doi":"10.1088/0067-0049/210/2/23","role":"CN wavelengths, energies, transition probabilities","asset":"data/reference/cno_molecular_primary/cn_brooke2014/table4.dat.gz","sha256":sha256(ROOT/'data/reference/cno_molecular_primary/cn_brooke2014/table4.dat.gz'),"status":"ACQUIRED","provenance_note":""},
        {"source_id":"Brooke2015_NH","source_type":"primary_published_table","citation":"Brooke et al. 2015, J. Chem. Phys. 143, 026101","ssot_key":"","doi":"10.1063/1.4923422","role":"NH wavelengths, energies, transition probabilities","asset":"data/reference/cno_molecular_primary/nh_brooke2015/BrookeEtAl-NH-2015-JCP.zip","sha256":sha256(ROOT/'data/reference/cno_molecular_primary/nh_brooke2015/BrookeEtAl-NH-2015-JCP.zip'),"status":"ACQUIRED","provenance_note":""},
        {"source_id":"Brooke2016_OH","source_type":"primary_published_table","citation":"Brooke et al. 2016, JQSRT 168, 142","ssot_key":"","doi":"10.1016/j.jqsrt.2015.07.021","role":"OH wavelengths, energies, transition probabilities","asset":"data/reference/cno_molecular_primary/oh_brooke2016/OH-Supplementary.zip","sha256":sha256(ROOT/'data/reference/cno_molecular_primary/oh_brooke2016/OH-Supplementary.zip'),"status":"ACQUIRED","provenance_note":""},
        {"source_id":"Masseron2014_CH","source_type":"primary_published_table","citation":"Masseron et al. 2014, A&A 571 A47","ssot_key":"","doi":"10.1051/0004-6361/201423956","role":"CH wavelengths, energies, transition probabilities","asset":"data/reference/cno_molecular_primary/ch_masseron2014/table14.dat.gz","sha256":sha256(ROOT/'data/reference/cno_molecular_primary/ch_masseron2014/table14.dat.gz'),"status":"ACQUIRED","provenance_note":""},
        {"source_id":"Li2015_CO","citation":"Li et al. 2015, ApJS 216, 15","doi":"10.1088/0067-0049/216/1/15","source_type":"redistribution","role":"12C16O transition data, TWICE DERIVED: ExoMol redistribution of Li2015, then converted to Turbospectrum babsma format","asset":"data/linelists/molecular/turbospectrum/CO/CO_IR_Li2015.dat","sha256":sha256(TS/'CO/CO_IR_Li2015.dat'),"status":"ACQUIRED_AS_REDISTRIBUTION","provenance_note":"NOT a primary Li 2015 ApJS 216,15 table. The asset's own second line reads 'ExoMol Li2015'; MOLECULAR_MANIFEST records it as an RYA-236 conversion of the ExoMol Li2015 CO list to Turbospectrum babsma (species code 0608.012016) by a converter EXTERNAL to this repo, so the conversion cannot be re-run or verified here. No Li 2015 primary table has been acquired: data/reference/cno_molecular_primary/ has no CO directory. This is the sole source behind all 80 12C16O PHYSICAL_TUPLE_MATCH rows -- the intake's entire clean-match class (RYA-1142 A4, RYA-1180)."},
        {"source_id":"BarklemCollet2016","citation":"Barklem & Collet 2016, A&A 588 A96","doi":"10.1051/0004-6361/201526961","source_type":"primary_published_table","role":"adopted molecular DISSOCIATION ENERGIES (table1.dat) -- partition functions and equilibrium constants NOT acquired","asset":"data/reference/cno_molecular_primary/constants_barklem2016/table1.dat","sha256":sha256(ROOT/'data/reference/cno_molecular_primary/constants_barklem2016/table1.dat'),"status":"ACQUIRED_PARTIAL","provenance_note":"The role previously read 'molecular partition functions and equilibrium constants' while the asset was table1.dat, which is the DISSOCIATION ENERGY table. Per the holding's own ReadMe: table1=dissociation energies (ON DISK), table6=partition functions (ABSENT), table7=equilibrium constants (ABSENT), table2/*=per-molecule constants (ABSENT; list.dat is only its filename index). See molecular_constants_ledger.csv for what is actually acquired (RYA-1142 A5, RYA-1180)."},
    ]
    write_csv(AUDIT / "source_bibliography.csv", sources, tuple(sources[0]))

    # 🔴 RYA-1180 (RYA-1142 A5). This ledger used to emit six BYTE-IDENTICAL rows all
    # claiming `PRIMARY_TABLES_ACQUIRED` with `partition_function_source = Barklem &
    # Collet 2016` -- for tables that are NOT ON DISK. Per the holding's own ReadMe:
    #   table1.dat  dissociation energies              ON DISK   (the only one)
    #   table6.dat  partition functions for molecules  ABSENT
    #   table7.dat  equilibrium constants              ABSENT
    #   table2/*    per-molecule constants             ABSENT (list.dat is its index)
    # So the ledger asserted acquisition of exactly what was missing and recorded none of
    # what was present. It now carries the six ADOPTED De values PARSED FROM the acquired
    # file -- never retyped, so the ledger cannot drift from its source (RYA-1101) -- with
    # the table id, the row label and the file checksum, per molecule.
    #
    # ⚠️ THE ADOPTED COLUMN IS NOT THE FIRST ONE. table1.dat carries four De columns:
    # HH (Huber & Herzberg 1979), Luo (2007), G2 theory, and `De` = "Dissociation energy
    # adopted" at bytes 84-93. They disagree: NH is HH 3.470, Luo 3.470, G2 3.378, and
    # ADOPTED 3.419 -- a value equal to none of its inputs. Read POSITIONALLY by the
    # ReadMe's byte offsets, never by whitespace order (RYA-853's comparison-column defect).
    constants = []
    de_table = ROOT / "data/reference/cno_molecular_primary/constants_barklem2016/table1.dat"
    de_sha = sha256(de_table)
    adopted = barklem_adopted_de(de_table)
    for molecule in ("C2", "CH", "CN", "NH", "OH", "CO"):
        de, e_de = adopted[molecule]
        constants.append({
            "molecule": molecule,
            "dissociation_energy_eV": de,
            "dissociation_energy_unc_eV": e_de,
            "dissociation_energy_source": "Barklem & Collet 2016, table1.dat col 'De' (bytes 84-93, 'Dissociation energy adopted')",
            "de_source_asset": "data/reference/cno_molecular_primary/constants_barklem2016/table1.dat",
            "de_source_sha256": de_sha,
            "de_source_row_label": molecule,
            "partition_function_source": "Barklem & Collet 2016 table6.dat",
            "partition_function_status": "NOT_ACQUIRED",
            "equilibrium_constant_status": "NOT_ACQUIRED (table7.dat)",
            "per_molecule_constants_status": "NOT_ACQUIRED (table2/* absent; list.dat is only its filename index)",
            "isotopic_assumption": "12C16O explicit for CO; other Table2 isotopologues implicit main species",
            "verdict": "DISSOCIATION_ENERGY_ACQUIRED; PARTITION_FUNCTIONS_AND_EQUILIBRIUM_CONSTANTS_NOT_ACQUIRED",
        })
    write_csv(AUDIT / "molecular_constants_ledger.csv", constants, tuple(constants[0]))

    conflict = [
        {"scope":"published molecular count","source_a":"CDS ReadMe/Table2","value_a":"408","source_b":"stale RYA-1131 comment","value_b":"879","decision":"408 is authoritative for published used transitions","status":"RESOLVED"},
        {"scope":"molecular identity","source_a":"Amarsi2021 Table2","value_a":"no rotational labels","source_b":"vendored synthesis lists","value_b":"heterogeneous upstream releases","decision":"only unique wavelength+EP+loggf tuple joins admitted","status":"OPEN_FOR_UNMATCHED"},
        {"scope":"N I adopted census","source_a":"AGSS21 summary","value_a":"five N I lines","source_b":"Amarsi2020 Table 1","value_b":"744.229, 821.633, 862.923, 868.340, 1010.890 nm","decision":"full primary-paper five-line census transcribed","status":"RESOLVED"},
    ]
    write_csv(AUDIT / "conflict_ledger.csv", conflict, tuple(conflict[0]))

    rejected = [
        {"species":"CN","system":"A-X red","band":"14 bands beyond (0-0)","count":"463","wavelength_region":"red/NIR","use_status":"REJECTED","reason":"automatic legacy equivalent widths produced two-to-three-times larger dispersion","evidence":"Amarsi2021 Sect. 2.1 lines 150-160","transitions_held_at":"","held_count_basis":""},
        {"species":"NH","system":"A-X","band":"unspecified","count":"6653","wavelength_region":"near-UV around 340 nm","use_status":"NOT_SELECTED_BY_SOURCE","reason":"Amarsi rejected these for crowding and continuum/blend limitations and did not publish WHICH lines were used. The SELECTION is unpublished; the TRANSITIONS are held in this repo and are counted here (RYA-1181/RYA-1142 A9).","evidence":"Amarsi2021 Sect. 2.1 lines 161-165","transitions_held_at":"data/reference/cno_molecular_primary/nh_brooke2014/NH-A-X-linelist.csv","held_count_basis":"RYA-1181 measured, near-UV 2000-4000 A"},
        {"species":"OH","system":"A-X","band":"unspecified","count":"5400","wavelength_region":"near-UV around 320 nm","use_status":"NOT_SELECTED_BY_SOURCE","reason":"Amarsi rejected these for crowding and continuum/blend limitations and did not publish WHICH lines were used. The SELECTION is unpublished; the TRANSITIONS are held in this repo and are counted here (RYA-1181/RYA-1142 A9).","evidence":"Amarsi2021 Sect. 2.1 lines 161-165","transitions_held_at":"data/reference/cno_molecular_primary/oh_brooke2016/OH-A-X-linelist-final.csv","held_count_basis":"RYA-1181 measured, near-UV 2000-4000 A"},
        {"species":"CN","system":"B-X","band":"unspecified","count":"45786","wavelength_region":"near-UV around 390 nm","use_status":"NOT_SELECTED_BY_SOURCE","reason":"Amarsi rejected these for crowding and continuum/blend limitations and did not publish WHICH lines were used. The SELECTION is unpublished; the TRANSITIONS are held in this repo and are counted here (RYA-1181/RYA-1142 A9).","evidence":"Amarsi2021 Sect. 2.1 lines 161-165","transitions_held_at":"data/reference/cno_molecular_primary/cn_brooke2014/table4.dat.gz","held_count_basis":"RYA-1181 measured, near-UV 2000-4000 A"},
    ]
    write_csv(AUDIT / "rejected_indicator_ledger.csv", rejected, tuple(rejected[0]))

    # 🔴 RYA-1181 (A9 / scope truth). RYA-1131 is titled "across FUV/NUV/IR" and RYA-1136
    # "UV-IR". The DELIVERED inventory is neither: both domains match ZERO rows in FUV and
    # ZERO in NUV, and the atomic census does not reach below 5052 A at all. Empty bins in
    # a coverage matrix read as "looked and found nothing"; here nothing was looked at.
    # Stated as its own record so the claim and the delivery can be compared without
    # reading two ticket titles. (`intake_verdict.json` is RYA-1183's to correct.)
    mol = list(csv.DictReader((AUDIT / "molecular_physical_crossmatch.csv").open()))
    atm = list(csv.DictReader((AUDIT / "atomic_source_census.csv").open()))
    mol_A = [float(r["wavelength_vac_nm"]) * 10 for r in mol]
    atm_A = [float(r["wavelength_air_A"]) for r in atm if r.get("wavelength_air_A", "").strip()]
    (AUDIT / "delivered_scope_rya1181.json").write_text(json.dumps({
        "ticket": "RYA-1181",
        "claimed_by_ticket_titles": "RYA-1131 'across FUV/NUV/IR'; RYA-1136 'UV-IR'",
        "delivered": "VIS to IR",
        "atomic_span_A": [round(min(atm_A), 1), round(max(atm_A), 1)],
        "molecular_span_A": [round(min(mol_A), 1), round(max(mol_A), 1)],
        "matched_rows_by_band": {
            "atomic": dict(Counter(r["source_band"] for r in atm)),
            "molecular": dict(Counter(r["source_band"] for r in mol)),
        },
        "fuv_matched": 0, "nuv_matched": 0,
        "why_the_uv_bins_are_empty": (
            "Not a null result. Amarsi's UV SELECTION is unpublished, so no UV row enters "
            "Table 2 and nothing is matched there. The UV TRANSITIONS are held in this "
            "repo and are counted in rejected_indicator_ledger.csv and "
            "held_uv_transitions_rya1181.csv (RYA-1181 A9). UV is a FOLLOW-ON, not a gap "
            "that was searched and came up empty."),
    }, indent=2) + "\n")

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
            # RYA-1144: these were previously argmin-resolved and counted as coverage.
            f"{mol_status.get('AMBIGUOUS_SUM_MATCH', 0)} features are reproduced by MORE THAN ONE "
            "subset of primary components; they are fitted, not identified, and are no longer "
            "counted as matched coverage",
            # RYA-1146: the only clean-match class rests on a converted redistribution.
            "All 80 CO PHYSICAL_TUPLE_MATCH rows resolve against an ExoMol->Turbospectrum "
            "conversion whose converter is not in this repo; the Li 2015 primary was never acquired",
            # RYA-1158: the scope claim and the delivered inventory disagree.
            "Titled UV-IR but delivered VIS-IR: zero FUV and zero NUV rows in both domains, and "
            "the NH A-X / OH A-X / CN B-X ultraviolet transitions are acquired but unread",
        ],
        "safety":("No abundance derived; no gf tuned; no wavelength-only join admitted "
                  "(RYA-1143: the N I five-line join is EP-aware as of this build, and the "
                  "canonical readback is reported in codex_* columns rather than as "
                  "published_*)"),
    }
    # 🔴 RYA-1183 B3(2) — A BLOCKER THAT HAS BEEN WORKED ON IS NOT A BLOCKER THAT IS GONE.
    #
    # RYA-1179/1180/1181 have all landed since these findings were written, and it would be
    # easy to read that as "resolved" and drop them. Two of the three are NOT resolved:
    # those tickets made the defect VISIBLE and RECORDED, which is a different thing from
    # removing the condition that blocks the freeze. Each blocker now carries its
    # disposition so the distinction survives the next reader.
    verdict["blocker_disposition"] = {
        "sum_matches_fitted_not_identified": {
            "addressed_by": "RYA-1150",
            "state": "RESOLVED_IN_ACCOUNTING",
            "detail": "The 26 ambiguous sum-matches are excluded from accepted coverage; "
                      "they are reported, not counted.",
        },
        "exomol_co_redistribution": {
            "addressed_by": "RYA-1180",
            "state": "STILL_BLOCKING",
            "detail": "Now labelled a redistribution and graded on all 80 rows "
                      "(source_provenance_grade=REDISTRIBUTION_VENDORED_SYNTHESIS_LIST), "
                      "but no Li 2015 primary table has been acquired. The intake's only "
                      "clean-match class still rests on a twice-derived list.",
        },
        "n_i_wavelength_only_admission": {
            "addressed_by": "RYA-1143, confirmed by RYA-1179",
            "state": "RESOLVED",
            "detail": "The N I join is EP-aware (nearest_canonical gates on wavelength "
                      "<=0.03 A AND EP <=0.002 eV in one conjunction). RYA-1179's "
                      "scope-aware guard does not flag it, so the safety line's "
                      "'no wavelength-only join admitted' is now TRUE rather than "
                      "aspirational.",
        },
        "uv_scope": {
            "addressed_by": "RYA-1181",
            "state": "STILL_BLOCKING",
            "detail": "Scope now stated as VIS-to-IR and ~105,858 held UV transitions are "
                      "counted, but zero FUV/NUV rows are delivered and the held UV is "
                      "staged rather than matched.",
        },
    }
    (AUDIT / "intake_verdict.json").write_text(json.dumps(verdict, indent=2)+"\n")

    # 🔴 RYA-1183 (RYA-1142 A6b-1) — AN EMPTY BIN IS NOT A FINDING UNTIL IT SAYS WHY.
    #
    # The atomic census is neutrals only, and nothing recorded whether that is the SOURCE's
    # shape or ours. Measured from the raw table rather than asserted, and the two causes
    # turn out to be different:
    #
    #   C II / N II / O II   ABSENT FROM THE SOURCE. table1.dat carries only CI, OI and
    #                        FeII species codes -- there is no ionised CNO to admit, so
    #                        this is not a filter and no ionised row may be fabricated.
    #   Fe II                PRESENT IN THE SOURCE AND DELIBERATELY EXCLUDED. 142 Fe II
    #                        rows sit in table1.dat and `atomic_census` drops them
    #                        (`if species not in {"CI","OI"}`). That IS a filter, and it is
    #                        recorded as one rather than left looking like an absence.
    #   N I                  NOT FROM THIS SOURCE AT ALL. table1.dat has no N; the five
    #                        N I rows come from Amarsi et al. 2020 Table 1.
    import collections as _c
    _raw = _c.Counter(l[0:4].strip() for l in ATOMIC.read_text().splitlines() if l[0:4].strip())
    (AUDIT / "atomic_census_scope_rya1183.json").write_text(json.dumps({
        "ticket": "RYA-1183",
        "scope": "NEUTRALS ONLY",
        "species_codes_present_in_source": dict(sorted(_raw.items())),
        "source_file": str(ATOMIC.relative_to(ROOT)),
        "ionised_cno": {
            "present_in_census": 0,
            "cause": "ABSENT_FROM_SOURCE",
            "statement": ("Amarsi 2019 Table 1 contains no C II / N II / O II — its only "
                          "species codes are CI, OI and FeII. The empty ionised-CNO bins "
                          "are the source's shape, NOT a filter we applied, and no "
                          "ionised row is fabricated to fill them."),
        },
        "fe_ii": {
            "present_in_source": _raw.get("FeII", 0),
            "present_in_census": 0,
            "cause": "EXCLUDED_BY_THIS_BUILDER",
            "statement": ("Fe II IS in the source and is deliberately dropped by "
                          "atomic_census() because this is a CNO intake. Recorded so an "
                          "absence-by-filter is never read as an absence-by-source."),
        },
        "n_i": {
            "present_in_census": 5,
            "cause": "DIFFERENT_SOURCE",
            "statement": ("table1.dat has no nitrogen at all. The five N I rows come from "
                          "Amarsi et al. 2020 Table 1, the AGSS21-adopted N I set."),
        },
    }, indent=2) + "\n")

    # 🔴 RYA-1150. summary.json is written by the INGEST script and was never revisited,
    # so it still advertised canonical_matched=0 / crossmatch_review=408 while this file
    # reported 390 in accepted classes -- two artifacts in one directory 390 rows apart,
    # and the stale one read clean. The builder now owns it, and a consistency assert
    # makes a future divergence fail loudly instead of sitting there.
    accepted = sum(v for k, v in mol_status.items() if k in ACCEPTED_MOLECULAR_JOINS)
    summary_path = AUDIT / "summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
    summary.update({
        "canonical_matched": accepted,
        "crossmatch_review": len(molecular) - accepted,
        "verdict": verdict["verdict"],
        "frozen_ready": verdict["frozen_ready_for_measurement"],
        "reconciled_with": "intake_verdict.json (RYA-1150)",
        # RYA-1183 B3(1): the status label is taken FROM intake_verdict.json above, so the
        # two artifacts cannot state different verdicts. `BLOCKED_MOLECULAR_DATA` appears
        # in no artifact and never did -- the recorded status is
        # INTAKE_COMPLETE_REVIEW_REQUIRED, which says something different and narrower:
        # the census IS complete, and what is blocked is the freeze, not the intake.
        "verdict_label_note": (
            "INTAKE_COMPLETE_REVIEW_REQUIRED, not BLOCKED_MOLECULAR_DATA. The census is "
            "complete (intake_census_complete=true); what is withheld is "
            "frozen_ready_for_measurement, on the blocking_findings listed in "
            "intake_verdict.json. Both artifacts take this string from one writer."),
    })
    summary_path.write_text(json.dumps(summary, indent=2)+"\n")
    assert summary["canonical_matched"] + summary["crossmatch_review"] == len(molecular), \
        "summary.json and intake_verdict.json disagree on the molecular total"


if __name__ == "__main__": main()
