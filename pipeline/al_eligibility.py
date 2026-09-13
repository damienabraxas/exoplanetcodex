"""RYA-1134 Al line/holding and engine coverage, with explicit non-run cells.

Catalog or declared loader extent is not real-pixel coverage. This planning
matrix preserves that distinction and cannot authorize a measurement by itself.
"""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd

from pipeline import band_policy, telluric_policy
from pipeline.gerber_nlte import DECKS
from pipeline.model_registry import LINE_SETS

POOL_AXES = {"reference": "reference_membership", "our-graded": "codex_membership",
             "our-deep-graded": "deep_membership", "asplund-al": "replication_membership"}


def matrices(ledger: pd.DataFrame, root: Path, *, holding_specs=None):
    if holding_specs is None:
        sys.path.insert(0, str(root / "scripts"))
        from measure_band_ew import _INSTRUMENT_HOLDINGS
        holding_specs = {h.holding_id: h for values in _INSTRUMENT_HOLDINGS.values() for h in values}
    holdings = pd.read_csv(root / "data/catalog/holdings_manifest_registry.csv").fillna("")
    holdings = holdings[holdings.system_id.eq("solar")]
    instruments = pd.read_csv(root / "data/catalog/instrument_catalog.csv").fillna("").set_index("instrument_id")
    specs = holding_specs
    manifest_columns = pd.read_csv(root / "data/audit/rya1132_al_intake/al_line_manifest.csv", nrows=0).columns
    missing_axes = [axis for axis in ("line_set", "telluric_applied", "normalization_state", "observed_conditioning")
                    if axis not in manifest_columns]
    schema_reason = "RYA1176_MISSING:" + ",".join(missing_axes) if missing_axes else ""
    rows = []
    for r in ledger.to_dict("records"):
        wavelength = r["wavelength_air"]
        try:
            policy = band_policy.resolve(wavelength)
            band = policy.name
        except band_policy.BandPolicyError:
            policy, band = None, "NO_DECLARED_BAND"
        for h in holdings.to_dict("records"):
            inst = instruments.loc[h["instrument_id"]]
            spec = specs.get(h["holding_id"])
            # The instrument extent is used only to prove OUTSIDE. A catalog
            # cannot prove a usable pixel exists inside its envelope.
            lo, hi = float(inst.wavelength_min_nm)*10, float(inst.wavelength_max_nm)*10
            extent_basis = "instrument catalog only"
            if spec and spec.span_A:
                lo, hi = spec.span_A
                extent_basis = "exact HoldingSpec declared span"
            coverage = "OUTSIDE_DECLARED_RANGE" if not lo <= wavelength <= hi else "IN_DECLARED_RANGE_PIXELS_UNVERIFIED"
            if spec and spec.span_A and coverage != "OUTSIDE_DECLARED_RANGE" and not spec.covers(wavelength, 1.):
                coverage = "EDGE_WINDOW_INCOMPLETE"
            state, evidence = telluric_policy.verified_band_state(h["holding_id"], wavelength)
            risks = []
            if r["atomic_handoff_status"] != "READY":
                risks.append("ATOMIC_HANDOFF_HOLD")
            if not spec:
                risks.append("NO_HOLDING_READER")
            if policy is None:
                risks.append("NO_DECLARED_BAND_POLICY")
            if not h["normalization_state"] or h["normalization_state"] == "unknown":
                risks.append("NORMALIZATION_UNKNOWN")
            susceptible = telluric_policy.in_telluric_band(wavelength)
            if (policy and policy.telluric_required) or susceptible:
                if state != "clean":
                    risks.append("EXACT_WINDOW_TELLURIC_VERIFICATION_REQUIRED")
            # RYA-1176 remains unmerged: unknown must not become native-as-delivered.
            risks.append("OBSERVED_CONDITIONING_UNESTABLISHED")
            if schema_reason:
                risks.append(schema_reason)
            if coverage == "OUTSIDE_DECLARED_RANGE":
                disposition, reason = "N/A", coverage
            else:
                disposition, reason = "HOLD", "|".join(risks + [coverage])
            rows.append(dict(canonical_line_id=r["canonical_line_id"], species=r["species"],
                wavelength_air=wavelength, band=band, holding_id=h["holding_id"], instrument_id=h["instrument_id"],
                reference_membership=r["reference_membership"], codex_membership=r["codex_membership"],
                deep_membership=r["deep_membership"], replication_membership=r["replication_membership"],
                verified_source_record_id=r["verified_source_record_id"], verified_source_class=r["verified_source_class"],
                gf_uncertainty_bound_dex=r["gf_uncertainty_bound_dex"], physical_identity_status=r["physical_identity_status"],
                hfs_component_status=r["hfs_component_status"], feature_depth=r.get("feature_depth", np.nan),
                blend_flag=r.get("blend_flag", "unknown"), coverage=coverage, coverage_basis=extent_basis,
                telluric_risk="SUSCEPTIBLE" if susceptible else "NOT_IN_REGISTERED_COMPLEX",
                telluric_applied=h["telluric_applied"] or "unknown", telluric_window_state=state,
                telluric_window_evidence=evidence, normalization_state=h["normalization_state"] or "unknown",
                normalization_verification="REGISTRY_DECLARATION_NOT_REMEASURED",
                observed_conditioning="unknown", route="NEITHER_UNTIL_GATES_PASS", status=disposition, reason=reason))
    line_matrix = pd.DataFrame(rows)
    models = pd.read_csv(root / "data/catalog/model_registry.csv").fillna("").to_dict("records")
    model_citations = {1: "kurucz2014_atlas;sneden1973;gerber2023",
        2: "gerber2023;gustafsson2008", 3: "amarsi2020_galah;wehrhahn2023_pysme;nordlander_lind2017",
        4: "gerber2023;gustafsson2008;nordlander_lind2017",
        5: "gerber2023;magic2013;nordlander_lind2017",
        6: "gerber2023;magic2013;nordlander_lind2017"}
    bibliography = pd.read_csv(root / "data/refs/bibliography.csv").set_index("key")
    for chain in model_citations.values():
        if any(key not in bibliography.index for key in chain.split(";")):
            raise ValueError(f"unresolved Al model bibliography chain: {chain}")
    # Models 5/6 already own the mean-3D treatment axes. Their Al binding is
    # Al@mean3D in derive_band_products, not an additional duplicate model.
    engine_rows = []
    for (band, holding, instrument), cell in line_matrix.groupby(["band", "holding_id", "instrument_id"], sort=True):
        for pool, axis in POOL_AXES.items():
            if pool not in LINE_SETS:
                raise ValueError(f"unregistered line_set {pool}")
            members = cell[cell[axis].eq("MEMBER")]
            in_range = members[~members.coverage.eq("OUTSIDE_DECLARED_RANGE")]
            for model in models:
                for route in ("EW", "synth"):
                    status, reason = "HOLD", ";".join(x for x in (schema_reason, "EXACT_HOLDING_INTAKE_REQUIRED") if x)
                    al_deck = "Al@mean3D" if model["model_id"] in (5, 6) else "Al" if model["model_id"] in (2, 4) else ""
                    if in_range.empty:
                        status, reason = "N/A", "NO_POOL_LINES_IN_DECLARED_HOLDING_RANGE"
                    elif route == "EW" and (pool == "our-deep-graded" or str(model["route"]) == "synth" or al_deck):
                        status, reason = "N/A", "SYNTHESIS_ONLY_POOL_OR_MODEL"
                    elif route == "EW" and band == "near-UV":
                        status, reason = "N/A", "BAND_POLICY_SYNTHESIS_ONLY_RYA713_RYA1133"
                    elif model["status"] in ("in-dev", "not-emitted"):
                        status, reason = "HOLD", "EXPERIMENTAL_OR_NO_EMITTER"
                    elif str(model["scale"]) == "3D-NLTE":
                        status, reason = "N/A", "NO_VALIDATED_FULL_3D_NLTE_AL_IMPLEMENTATION"
                    elif al_deck and al_deck not in DECKS:
                        status, reason = "HOLD", "AL_DECK_NOT_REGISTERED"
                    elif model["model_id"] == 3:
                        reason += ";AL_CORRECTION_GRID_PHYSICAL_IDENTITY_REVALIDATION_REQUIRED"
                    engine_rows.append(dict(band=band, holding_id=holding, instrument_id=instrument,
                        line_set=pool, route=route, model_id=model["model_id"], engine=model["display_name"],
                        al_implementation="PySME-derived additive 1D-NLTE corrections" if model["model_id"] == 3 else "Gerber Al synthesis" if al_deck else "REGISTRY_ROUTE_REQUIRES_AL_VALIDATION",
                        bibliography_chain=model_citations.get(model["model_id"], "NOT_ADOPTED_AL_MODEL"),
                        treatment=model["stored_token"], scale=model["scale"],
                        al_deck_key=al_deck,
                        model_binding_status="REGISTERED_NOT_EXECUTION_VALIDATED" if al_deck in DECKS else "AL_BINDING_REVIEW_REQUIRED",
                        n_pool_lines=len(members), n_in_declared_range=len(in_range),
                        n_verified_eligible=0, status=status, reason=reason,
                        applicability_evidence="current model_registry.csv; gerber_nlte.DECKS; no stale HAVE/DISK_ONLY inference"))
    return line_matrix, pd.DataFrame(engine_rows)
