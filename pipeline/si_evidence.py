"""RYA-1218/754: audit existing Si results without promoting partial budgets.

Atomic grade, published-list membership, measurement acceptance and uncertainty
completeness are independent. This module changes none of the input abundances.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

from pipeline import line_match
from pipeline import uncertainty_contract as uc


def attach_identity(lines, canonical):
    """Resolve each species separately with wavelength AND excitation energy."""
    output = []
    for (element, ion), group in lines.groupby(["element", "ion"], sort=False):
        source = canonical[canonical.species == f"{element} {ion}"].reset_index(drop=True)
        if source.empty:
            raise line_match.LineMatchError(f"No atomic rows for {element} {ion}")
        if source.duplicated(["wavelength_air_A", "excitation_potential_eV"]).any():
            raise ValueError("duplicate atomic identities require adjudication")
        result = line_match.match(
            group.wavelength_air_A, source.wavelength_air_A,
            want_ep=group.ep_eV, src_ep=source.excitation_potential_eV,
            require_ep=True)
        indexes = line_match.require_resolved(result, what="Si result audit", species=f"{element} {ion}")
        for (_, row), idx in zip(group.iterrows(), indexes):
            atom = source.iloc[idx]
            output.append({**row.to_dict(), "line_id": atom.line_id,
                           "snapshot_loggf": float(atom.log_gf),
                           "snapshot_gf_source": atom.loggf_reference})
    return pd.DataFrame(output)


def accepted(lines):
    mask = lines.in_aggregate.astype(str).str.lower().eq("true")
    result = lines[mask].copy()
    if not np.isfinite(result.abundance.to_numpy(float)).all():
        raise ValueError("accepted line has no finite abundance")
    if result.line_id.duplicated().any():
        raise ValueError("accepted pool repeats a physical line")
    return result


def diagnostic_statistics(lines, reported):
    """Show estimator choice; SE(mean) is not relabelled as sigma(median)."""
    values = accepted(lines).abundance.to_numpy(float)
    if len(values) != int(reported["n_lines"]):
        raise ValueError("product count disagrees with its accepted physical pool")
    if not len(values):
        return {"n": 0, "median": None, "mean": None, "scatter": None, "sem_mean": None}
    median, mean = float(np.median(values)), float(np.mean(values))
    if not math.isclose(median, float(reported["A"]), abs_tol=0.00051):
        raise ValueError("reported abundance is not the archived pool median")
    scatter = float(np.std(values, ddof=1)) if len(values) > 1 else None
    return {"n": len(values), "median": median, "mean": mean, "scatter": scatter,
            "sem_mean": scatter / math.sqrt(len(values)) if scatter is not None else None,
            "reported_estimator": "median", "reported_stat_dex": float(reported["stat_dex"]),
            "statistical_admission": "HOLD_estimator_and_independence_assessment"}


def paired_difference(left, right):
    """Diagnostic only: difference actual matched lines, never unmatched totals."""
    a, b = accepted(left).set_index("line_id"), accepted(right).set_index("line_id")
    common = sorted(set(a.index) & set(b.index))
    if not common:
        raise ValueError("no common accepted physical lines")
    delta = a.loc[common, "abundance"] - b.loc[common, "abundance"]
    return {"matched_ids": common, "n": len(common),
            "left_only": sorted(set(a.index) - set(common)),
            "right_only": sorted(set(b.index) - set(common)),
            "mean_delta": float(delta.mean()), "median_delta": float(delta.median()),
            "per_line_delta": dict(zip(common, map(float, delta))),
            "sigma_differential": None,
            "status": "DIAGNOSTIC_joint_covariance_not_supplied"}


def hold_budget(product, ids, *, source, statistics):
    """Use RYA-587 unchanged; archived placeholders cannot populate real terms."""
    reasons = {
        "measurement": "Archived aggregate is a median; its scatter/sqrt(N) is an SE of the mean. Independence and estimator-specific uncertainty are not established.",
        "transition_data": "Archived 0.17 dex is an explicitly unmeasured gf placeholder. Need the used transition scale, per-line evidence and shared-source covariance.",
        "stellar.teff": "No matched-pool adopted-uncertainty response supplied.",
        "stellar.logg": "No matched-pool adopted-uncertainty response supplied.",
        "stellar.xi": "1D synthesis has xi; no matched-pool response supplied.",
        "stellar.metallicity": "No matched-pool adopted-uncertainty response supplied.",
        "continuum": "A shipped continuum is not a measurement of its uncertainty.",
        "profile_ew": "Fit curvature and fit-quality diagnostics are not a validated noise/correlation likelihood.",
        "pseudo_continuum": "Requires assessed applicability and overlap with continuum/profile terms.",
        "telluric": "Correction state alone does not quantify residual uncertainty or establish N/A at each line.",
        "holding_instrument": "Matched holding differences available diagnostically; no adopted covariance model.",
        "nlte": "Partial line corrections do not establish model uncertainty on the entire pool.",
        "model_atmosphere": "Labeled engine comparisons are not automatically independent atmosphere realizations.",
        "hfs_isotopes": "Applicability must be supported for the physical Si transitions.",
        "blends": "Accepted minima do not establish absence or uncertainty of blends.",
        "molecular_coupling": "Atomic Si still requires assessment of molecular background/blend coupling.",
    }
    terms = [uc.component(name, None, state="HOLD", source=source,
                          evidence={"reason": reasons[name],
                                    **({"archived_statistics": statistics} if name == "measurement" else {})})
             for name in uc.COMPONENTS]
    scope = uc.product_scope(product, star="solar", indicator_ids=ids)
    return uc.assemble(scope, terms, covariance=[], covariance_source=source,
                       assumptions="No missing covariance or component is inferred as zero.")
