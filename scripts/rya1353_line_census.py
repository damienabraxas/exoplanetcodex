#!/usr/bin/env python3
"""Read-only O I/II and Fe I/II holdings census; never grades or adopts a transition.

Run from the repository with its Python environment. Output is an exclusive new
directory. Source rows remain separate: a tuple candidate is not a physical-ID proof.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import platform
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd

from config.synth_bands import SYNTH_BANDS
from pipeline.gf_grades import CITATIONS, is_gf_graded
from pipeline.line_match import EP_TOL_EV, MATCH_TOL_A, match
from pipeline.si_protocol_audit import holding_specs

SPECIES = ("O I", "O II", "Fe I", "Fe II")
_MATCH_INPUTS = {}


def number(value):
    if value is None or str(value).strip() in ("", "nan", "None"):
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Nonfinite input: {value!r}")
    return result


def read(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader((s for s in stream if not s.startswith("#")),
                                   delimiter="\t" if path.suffix == ".tsv" else ","))


def write(path, rows):
    if not rows:
        raise ValueError(f"Refusing empty table without a schema: {path}")
    fields = list(dict.fromkeys(k for r in rows for k in r))
    payload = io.StringIO(newline="")
    writer = csv.DictWriter(payload, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    data = payload.getvalue().encode("utf-8")
    if len(data) > 1_000_000:
        path = path.with_suffix(path.suffix + ".gz")
        data = gzip.compress(data, mtime=0)
    with path.open("xb") as stream:
        stream.write(data)


def band(wave):
    for name, spec in SYNTH_BANDS.items():
        if spec.lo_A <= wave < spec.hi_A:
            return name
    # Explicit census extensions, not new synthesis bands or engine support.
    if wave < min(s.lo_A for s in SYNTH_BANDS.values()):
        return "UV_below_synthesis_floor"
    return "outside_configured_synthesis_bands"


def medium(wave):
    # Canonical ASD pulls declare vac+air: vacuum below 2000 A, air above.
    return "vacuum" if wave < 2000 else "air"


def tuple_candidates(wave, ep, rows, wave_key, ep_key, tol=MATCH_TOL_A):
    """Strict shared matcher; unknown EP and duplicate wavelengths never become joins.

    Rows must already be restricted to species and a common wavelength medium.
    No physical membership, source grade or component identity is transferred.
    """
    if wave is None or ep is None:
        return "MISSING_IDENTITY_HOLD", None
    key = (id(rows), wave_key, ep_key)
    if key not in _MATCH_INPUTS:
        eligible = [r for r in rows if number(r.get(wave_key)) is not None
                    and number(r.get(ep_key)) is not None]
        waves = np.array([number(r[wave_key]) for r in eligible])
        eps = np.array([number(r[ep_key]) for r in eligible])
        _MATCH_INPUTS[key] = (rows, eligible, waves, eps, Counter(waves))
    _, eligible, waves, eps, counts = _MATCH_INPUTS[key]
    result = match([wave], waves, want_ep=[ep],
                   src_ep=eps,
                   tol_A=tol, require_ep=True)
    if result.ambiguous:
        return "AMBIGUOUS_HOLD", None
    if not result.resolved[0]:
        return "NO_TUPLE_CANDIDATE_IN_SNAPSHOT", None
    selected = eligible[int(result.index[0])]
    # The shared matcher collapses equal-wavelength rows. A census cannot assume
    # identical wavelength means identical upper level/component or identical gf.
    if counts[number(selected[wave_key])] > 1:
        return "DUPLICATE_WAVELENGTH_HOLD", None
    return "TUPLE_CANDIDATE_COMPONENT_REVIEW", selected


def queue(row):
    """Ordinal REVIEW order, never measurement eligibility or solver feasibility."""
    if row["species"].startswith("O"):
        family = "O_FIRST"
    else:
        family = "FE_STRESS_TEST"
    if row["gf_availability"] == "NO_FINITE_GF_IN_CANONICAL":
        return family + ":GF_DATA_REVIEW"
    if row["experimental_gf_state"] == "RECORDED_LAB":
        return family + ":WITHHELD_LAB_BENCHMARK"
    if row["predicted_depth_state"] in ("PREDICTED_WEAK", "PREDICTED_STRONG"):
        return family + ":OBSERVABILITY_REVIEW"
    return family + ":IDENTITY_AND_OBSERVABILITY_HOLD"


def build(out):
    _MATCH_INPUTS.clear()
    out.mkdir(parents=True, exist_ok=False)
    inputs = set()

    def take(rel):
        path = ROOT / rel
        inputs.add(path)
        return read(path)

    canonical = [r for r in take("data/linelists/canonical_gf.csv") if r["species"] in SPECIES]
    if not canonical or len({r["line_id"] for r in canonical}) != len(canonical):
        raise ValueError("Empty canonical scope or duplicate canonical line IDs")
    by_species = {sp: [r for r in canonical if r["species"] == sp] for sp in SPECIES}
    current_mask = is_gf_graded(pd.DataFrame(canonical)).tolist()
    raw = take("data/linelists/linelist_solar.csv")
    raw_scope = [r for r in raw if f"{r['element']} {r['ion']}" in SPECIES]
    raw_species = {sp: [r for r in raw_scope if f"{r['element']} {r['ion']}" == sp] for sp in SPECIES}
    write(out / "solar_intake_rows.csv", [dict(source_row_id=f"solar:{i}", **r)
                                           for i, r in enumerate(raw_scope, 2)])

    # A normalized holding table preserves actual declared file spans independently
    # of catalog capability. It does not open spectra, infer SNR or correct tellurics.
    holdings = take("data/catalog/holdings_manifest_registry.csv")
    instruments = {r["instrument_id"]: r for r in take("data/catalog/instrument_catalog.csv")}
    solar_holdings = [h for h in holdings if h["system_id"] == "solar"]
    # Reuse the existing AST adapter; importing the measurement driver would load
    # compute inputs. Exact runtime holding IDs can carry their own registry state.
    runtime_path = ROOT / 'scripts/measure_band_ew.py'
    inputs.update((runtime_path, ROOT / 'pipeline/si_protocol_audit.py'))
    runtime = holding_specs(runtime_path)
    spans = []
    for h in solar_holdings:
        path = ROOT / h["manifest_path"]
        if path.exists() and path.suffix in (".csv", ".tsv"):
            rows = take(h["manifest_path"])
            rows = [r for r in rows if r.get("instrument_id") == h["instrument_id"] and "loader" in r]
        else:
            rows = []
        for i, r in enumerate(rows):
            spans.append(dict(span_id=f"{h['holding_id']}:{i}", holding_id=h["holding_id"],
                              instrument_id=h["instrument_id"],
                              wave_min_A=r.get("wave_min_A", ""), wave_max_A=r.get("wave_max_A", ""),
                              span_status=r.get("span_status", ""),
                              association="MANIFEST_INSTRUMENT_SCOPE_ONLY_EXACT_HOLDING_MAPPING_PENDING",
                              telluric_applied="UNKNOWN_EXACT_PRODUCT_MAPPING_PENDING",
                              normalization_state="UNKNOWN_EXACT_PRODUCT_MAPPING_PENDING",
                              manifest_path=h["manifest_path"], raw_manifest_row=json.dumps(r, sort_keys=True)))
        spec = runtime.get(h['holding_id'])
        if spec and spec['span'] is not None:
            if spec['instrument_id'] != h['instrument_id']:
                raise ValueError(f"Runtime/registry instrument mismatch: {h['holding_id']}")
            spans.append(dict(span_id=f"{h['holding_id']}:runtime", holding_id=h['holding_id'],
                instrument_id=h['instrument_id'], wave_min_A=spec['span'][0], wave_max_A=spec['span'][1],
                span_status='RUNTIME_DECLARATION_NOT_PIXEL_VERIFICATION',
                association='EXACT_RUNTIME_HOLDING_ID', telluric_applied=h['telluric_applied'],
                normalization_state=h['normalization_state'], manifest_path='scripts/measure_band_ew.py',
                raw_manifest_row=json.dumps(spec, sort_keys=True)))
    write(out / "holding_spans.csv", spans)
    write(out / "solar_holdings.csv", [dict(**h, catalog_capability=json.dumps(instruments[h['instrument_id']], sort_keys=True),
                 n_manifest_instrument_spans=sum(s['holding_id'] == h['holding_id'] for s in spans),
                 runtime_holding_spec=json.dumps(runtime.get(h['holding_id']), sort_keys=True),
                 span_mapping_state=('EXACT_RUNTIME_DECLARED_SPAN_PIXEL_REVIEW_PENDING'
                     if h['holding_id'] in runtime and runtime[h['holding_id']]['span'] is not None
                     else 'HOLD_DYNAMIC_LOADER_OR_MANIFEST_ASSOCIATION'),
                 coverage_interpretation="CATALOG_CAPABILITY_IS_NOT_ACTUAL_LINE_COVERAGE") for h in solar_holdings])

    sources = []
    for sp, stem in (("O I", "OI"), ("O II", "OII"), ("Fe I", "FeI"), ("Fe II", "FeII")):
        for path in sorted((ROOT / "data/linelists/primary_gf").glob(f"nist_asd_{stem}_*.tsv")):
            rows = take(str(path.relative_to(ROOT)))
            prov = path.with_suffix(".prov.json")
            if prov.exists():
                inputs.add(prov)
            for i, r in enumerate(rows, 2):
                sources.append(dict(species=sp, source_kind="NIST_EVALUATED_SNAPSHOT",
                     source_file=str(path.relative_to(ROOT)), source_row=i,
                     wavelength_A=r.get("wavelength_A", ""), lower_EP_eV=r.get("ei_eV", ""),
                     source_loggf=r.get("log_gf", ""), source_sigma_dex="", source_doi="10.18434/T4W30F",
                     raw_record=json.dumps(r, sort_keys=True)))
    for sp, file in (("Fe I", "fe1_lab_loggf.csv"), ("Fe I", "fe1_lab_loggf_ruffoni2013.csv"),
                     ("Fe II", "fe2_lab_loggf_dh19.csv")):
        rel = "data/reference/fe_gf_lab/" + file
        for i, r in enumerate(take(rel), 2):
            sig = number(r["e_loggf_dex"])
            if sig is None or sig <= 0:
                raise ValueError(f"Invalid primary lab uncertainty: {rel}:{i}")
            sources.append(dict(species=sp, source_kind="PRIMARY_EXPERIMENTAL_TABLE",
                   source_file=rel, source_row=i, wavelength_A=r["wavelength_air_A"],
                   lower_EP_eV=r["elo_eV"], source_loggf=r["loggf"], source_sigma_dex=r["e_loggf_dex"],
                   source_doi=CITATIONS[r['source']][1], raw_record=json.dumps(r, sort_keys=True)))
        prov = (ROOT / rel).with_suffix(".prov.json")
        if prov.exists():
            inputs.add(prov)
    for i, r in enumerate(take("data/audit/rya1136_cno_intake/atomic_source_census.csv"), 2):
        if r["species"] in SPECIES:
            sources.append(dict(species=r["species"], source_kind="PUBLISHED_O_SOURCE_CENSUS",
                 source_file="data/audit/rya1136_cno_intake/atomic_source_census.csv", source_row=i,
                 wavelength_A=r["wavelength_air_A"] or r["wavelength_vac_A"],
                 lower_EP_eV=r["lower_EP_eV"], source_loggf=r["published_loggf"],
                 source_sigma_dex="", source_doi="", raw_record=json.dumps(r, sort_keys=True)))
    for sp in ("Fe I", "Fe II"):
        rel = f"data/reference/line_sets/fe_agss21_{sp.replace(' ', '')}.csv"
        for i, r in enumerate(take(rel), 2):
            sources.append(dict(species=sp, source_kind="PUBLISHED_FE_REFERENCE_SET",
                 source_file=rel, source_row=i, wavelength_A=r["wavelength_air_A"],
                 lower_EP_eV=r["ep_eV"], source_loggf=r["loggf_published"],
                 source_sigma_dex="", source_doi=r["doi"], raw_record=json.dumps(r, sort_keys=True)))

    alternatives = defaultdict(list)
    for i, r in enumerate(sources):
        r["source_record_id"] = f"source:{i}"
        w, ep = number(r["wavelength_A"]), number(r["lower_EP_eV"])
        r["wavelength_medium"] = medium(w) if w is not None else "UNKNOWN"
        status, hit = tuple_candidates(w, ep, by_species[r["species"]],
                                      "wavelength_air_A", "excitation_potential_eV")
        r["tuple_status"] = status
        r["candidate_canonical_line_id"] = hit["line_id"] if hit else ""
        r["identity_gate"] = "HOLD_UPPER_LEVEL_J_MULTIPLET_ISOTOPE_HFS_NOT_PROVEN"
        r["gf_state_in_source"] = "GF_PRESENT" if number(r["source_loggf"]) is not None else "NO_GF_IN_THIS_SOURCE_ROW"
        r["delta_loggf_source_minus_canonical"] = (
            number(r["source_loggf"]) - number(hit["log_gf"])
            if hit and number(r["source_loggf"]) is not None and number(hit["log_gf"]) is not None else "")
        if hit:
            alternatives[hit["line_id"]].append(r)
    write(out / "source_records.csv", sources)

    exclusions = take("data/catalog/line_curation_exclusions.csv")
    problems = take("data/registry/problem_children.csv")
    write(out / "exclusion_evidence.csv", [dict(kind="CURATED_EXCLUSION", **r) for r in exclusions if r['species'] in SPECIES]
          + [dict(kind="PROBLEM_CHILD_SCOPE_NOT_IDENTITY_JOIN", **r) for r in problems if r['species'] in SPECIES])
    write(out / "source_rejected_indicators.csv", take("data/audit/rya1136_cno_intake/rejected_indicator_ledger.csv"))
    # Preserve all Table 6 rows, including missing Ladenburg entries and upper
    # limits, and explicitly firewall the two solar-refined alternative columns.
    table6 = ROOT / 'data/reference/fe_gf_lab/ruffoni2013/apj485713t6_ascii.txt'
    inputs.add(table6)
    excluded_alternatives = []
    upper = ''
    for n, line in enumerate(table6.read_text().splitlines(), 1):
        fields = line.split('\t')
        if len(fields) < 12:
            continue
        try:
            wavelength_nm = float(fields[2])
        except ValueError:
            continue  # Published headings/residual-BF rows, not transitions.
        upper = fields[0] or upper
        for label, col in (('BF_AND_EFFECTIVE_LIFETIME', 7), ('RECOMMENDED', 11)):
            excluded_alternatives.append(dict(species='Fe I', source_file=str(table6.relative_to(ROOT)),
                source_line=n, upper_term=upper, lower_term=fields[1], wavelength_air_nm=wavelength_nm,
                alternative_column=label, value_as_printed=fields[col],
                disposition='EXCLUDED_FROM_INDEPENDENT_SOLAR_BENCHMARK_RYA161',
                reason='SOLAR_REFINED_LIFETIME_OR_RECOMMENDED_SCALE_NOT_INDEPENDENT_EXPERIMENT',
                raw_source_line_json=json.dumps(line)))
    write(out / 'excluded_gf_alternatives.csv', excluded_alternatives)
    records = []
    for r, graded in zip(canonical, current_mask):
        row = dict(r)
        w, ep = number(r["wavelength_air_A"]), number(r["excitation_potential_eV"])
        if w is None or w <= 0:
            raise ValueError(f"Invalid canonical wavelength: {r['line_id']}")
        row["wavelength_medium"] = medium(w)
        row["medium_basis"] = "CANONICAL_COLUMN_NAME_AIR_ASD_VAC_BELOW_2000_CONVENTION_REQUIRES_UPSTREAM_REVIEW"
        row["band"] = band(w)
        row["current_code_gf_graded"] = graded
        row["experimental_gf_state"] = "RECORDED_LAB" if "LAB" in r["gf_tier"] else "NOT_RECORDED_AS_LAB_NOT_PROOF_NONE_EXISTS"
        row["gf_availability"] = "FINITE_CANONICAL_GF" if number(r["log_gf"]) is not None else "NO_FINITE_GF_IN_CANONICAL"
        status, intake = tuple_candidates(w, ep, raw_species[r['species']], "wavelength_air_A", "excitation_potential_eV")
        row["solar_intake_tuple_status"] = status
        depth = number(intake["central_depth"]) if intake else None
        row["predicted_central_depth"] = depth
        row["predicted_depth_state"] = ("UNKNOWN" if depth is None else "BELOW_PREDICTED_FLOOR" if depth < .05
                                        else "PREDICTED_WEAK" if depth <= .60 else "PREDICTED_STRONG")
        row["observed_feature_depth"] = "UNKNOWN_NOT_MEASURED_BY_CENSUS"
        row["codex_deep_membership"] = "HOLD_OBSERVED_DEPTH_AND_GRADE_POLICY"
        row["blend_flag_intake"] = intake["blend_flag"] if intake else "UNKNOWN"
        row["alternative_record_ids"] = ";".join(s["source_record_id"] for s in alternatives[r['line_id']])
        row["primary_lab_tuple_candidate_ids"] = ";".join(s["source_record_id"] for s in alternatives[r['line_id']]
                                                            if s["source_kind"] == "PRIMARY_EXPERIMENTAL_TABLE")
        row["published_reference_tuple_candidate_ids"] = ";".join(s["source_record_id"] for s in alternatives[r['line_id']]
                                                            if s["source_kind"].startswith("PUBLISHED"))
        row["identity_missing"] = "upper_level;J_lower;J_upper;multiplet;isotope;HFS_component_ids"
        row["scientific_use"] = "HOLD_IDENTITY_AND_HOLDING_CONDITIONING"
        row["theory_tractability"] = "O_ATOMIC_STRUCTURE_REVIEW" if r['species'].startswith('O') else "FE_OPEN_SHELL_MIXING_REQUIRES_BENCHMARK"
        row["diagnostic_value"] = "REFERENCE_TUPLE_CANDIDATE" if row["published_reference_tuple_candidate_ids"] else "NOT_ASSESSED"
        row["mdwarf_suitability"] = "HOLD_MOLECULAR_BLENDS_PRESSURE_BROADENING_MAGNETISM_AND_TARGET_DEPTH"
        row["review_queue"] = queue(row)
        records.append(row)
    write(out / "canonical_census.csv", records)
    by_id = {r['line_id']: r for r in records}
    priorities = []
    for r in sources:
        c = by_id.get(r['candidate_canonical_line_id'])
        experimental = r['source_kind'] == 'PRIMARY_EXPERIMENTAL_TABLE'
        # Ordinal triage, defined before any calculation: O I first; existing
        # experimental data are benchmark/ingest work, not missing-lab targets.
        group = ('BENCHMARK_OR_INGEST_EXISTING_EXPERIMENT' if experimental else
                 'O_I_THEORY_REVIEW' if r['species'] == 'O I' else
                 'FE_THEORY_STRESS_TEST_REVIEW' if r['species'].startswith('Fe') else
                 'O_II_TARGET_APPLICABILITY_REVIEW')
        priorities.append(dict(source_record_id=r['source_record_id'], species=r['species'],
              wavelength_A=r['wavelength_A'], wavelength_medium=r['wavelength_medium'],
              lower_EP_eV=r['lower_EP_eV'], candidate_canonical_line_id=r['candidate_canonical_line_id'],
              review_group=group, source_gf_state=r['gf_state_in_source'],
              experimental_status='EXPERIMENTAL_TABLE_HELD' if experimental else 'NOT_ESTABLISHED_BY_THIS_SOURCE',
              identity_status=r['tuple_status'],
              observability='PENDING_EXACT_HOLDING_REVIEW',
              predicted_depth_state=c['predicted_depth_state'] if c else 'UNKNOWN',
              diagnostic_value=c['diagnostic_value'] if c else 'SOURCE_ONLY_REQUIRES_DIAGNOSTIC_REVIEW',
              tractability='O_STRUCTURE_BENCHMARK_REQUIRED' if r['species'].startswith('O') else 'FE_LEVEL_MIXING_CONVERGENCE_REQUIRED',
              mdwarf_suitability='HOLD_TARGET_MOLECULAR_BLEND_PRESSURE_MAGNETIC_NLTE_REVIEW',
              exclusion_reason='NOT_ADMITTED_IDENTITY_OR_OBSERVABILITY_OR_BENCHMARK_UNVALIDATED'))
    write(out / "candidate_review.csv", priorities)

    # Compact normalized eligibility matrix: explicit line membership per holding,
    # no per-row claim that range overlap means an observed usable line.
    matrix = []
    for s in spans:
        lo, hi = number(s['wave_min_A']), number(s['wave_max_A'])
        if lo is None or hi is None or lo > hi:
            raise ValueError(f"Invalid holding span {s}")
        for sp in SPECIES:
            candidates = [r for r in records if r['species'] == sp and lo <= number(r['wavelength_air_A']) <= hi]
            matrix.append(dict(span_id=s['span_id'], species=sp, canonical_range_candidate_count=len(candidates),
                       line_ids=";".join(r['line_id'] for r in candidates),
                       coverage="DECLARED_RANGE_CANDIDATES_NOT_GAP_OR_EDGE_VERIFIED",
                       association=s['association'],
                       eligibility="HOLD_EXACT_FRAME_GAPS_SNR_BLEND_CONTINUUM_TELLURIC_AND_MODEL_CHECK",
                       telluric_applied=s['telluric_applied'], normalization_state=s['normalization_state']))
    write(out / "holding_line_matrix.csv", matrix)

    engine = take("data/catalog/engine_coverage.csv")
    write(out / "engine_scope.csv", [dict(**r, census_interpretation="DECLARED_ELEMENT_BAND_SCOPE_NOT_TRANSITION_VALIDATION")
                                         for r in engine if r.get('element') in ('O', 'Fe')])
    products = []
    for rel in ("data/products/solar/O.json", "data/products/solar/Fe.json"):
        path = ROOT / rel
        inputs.add(path)
        obj = json.loads(path.read_text())
        # Products describe current published usage; do not reverse infer line IDs.
        for section in ('products', 'archive', 'quarantine'):
            for i, product in enumerate(obj.get(section, [])):
                if not isinstance(product, dict):
                    raise ValueError(f"Unexpected product structure in {rel}:{section}")
                products.append(dict(source_file=rel, source_section=section, source_index=i,
                    **{k: product.get(k, '') for k in ('element', 'ion', 'band', 'instrument', 'holding',
                        'tier', 'selector', 'treatment', 'route', 'grade', 'n_lines', 'n_excluded')},
                    census_interpretation="EXISTING_PRODUCT_SELECTION_NOT_LINE_IDENTITY_PROOF"))
    write(out / "product_usage.csv", products)
    bib = take("data/refs/bibliography.csv")
    doi_to_keys = defaultdict(list)
    for r in bib:
        if r['doi']:
            doi_to_keys[r['doi'].lower()].append(r['key'])
    write(out / "canonical_provenance_resolution.csv", [dict(line_id=r['line_id'],
          doi=r['gf_source_doi'], loggf_reference=r['loggf_reference'],
          bibliography_keys=";".join(doi_to_keys[r['gf_source_doi'].lower()]),
          resolution="DOI_RESOLVED" if doi_to_keys[r['gf_source_doi'].lower()] else "HOLD_REFERENCE_CODE_OR_MISSING_DOI_REQUIRES_SOURCE_REVIEW") for r in canonical])
    for path in [ROOT / "config/synth_bands.yaml", ROOT / "pipeline/gf_grades.py", ROOT / "pipeline/line_match.py", Path(__file__)]:
        inputs.add(path)
    write(out / "input_manifest.csv", [dict(path=str(p.relative_to(ROOT)), bytes=p.stat().st_size,
          sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(inputs)])
    grouped = Counter((r['species'], r['band'], r['gf_tier'], r['current_code_gf_graded']) for r in records)
    write(out / "coverage_counts.csv", [dict(species=k[0], band=k[1], gf_tier=k[2], current_code_gf_graded=k[3],
                                           canonical_rows=v) for k, v in sorted(grouped.items())])
    summary = dict(schema="rya1353.census/1", scope="COMMITTED_HOLDINGS_SNAPSHOT_NOT_EXHAUSTIVE_ATOMIC_UNIVERSE",
        base_sha=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        python=sys.version, platform=platform.platform(), canonical_rows=len(records), solar_intake_rows=len(raw_scope),
        source_rows=len(sources), species_counts=dict(Counter(r['species'] for r in records)),
        source_match_states=dict(Counter(r['tuple_status'] for r in sources)),
        review_queue_counts=dict(Counter(r['review_queue'] for r in records)),
        all_canonical_rows_retained=True, all_source_rows_retained=True,
        measured_or_published=False, identity_confirmed_matches=0,
        wave_tolerance_A=MATCH_TOL_A, ep_tolerance_eV=EP_TOL_EV,
        current_grade_policy_conflict="RYA-1351_PRIMARY_LAB_VS_RYA-1233_CURRENT_CODE_REQUIRES_RYA-1352_REVIEW")
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(f"Read {len(inputs)} files; wrote {len(list(out.iterdir()))} census artifacts to {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    build(parser.parse_args().out)
