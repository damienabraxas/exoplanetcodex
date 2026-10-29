"""RYA-587 publication contract, extending error_budget's measured/HOLD/N/A terms.

Legacy ErrorBudget and uncertainty_stack results remain diagnostics until their
evidence satisfies this contract. No element or treatment can bypass it.
All numeric components are already propagated to abundance units (dex).
"""
from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping

import numpy as np

from pipeline.error_budget import Term

SCHEMA = "rya587.uncertainty.v1"
COMPONENTS = (
    "measurement", "transition_data", "stellar.teff", "stellar.logg",
    "stellar.xi", "stellar.metallicity", "continuum", "profile_ew",
    "pseudo_continuum", "telluric", "holding_instrument", "nlte",
    "model_atmosphere", "hfs_isotopes", "blends", "molecular_coupling",
)
STATES = {"MEASURED", "DEFINED", "N/A", "HOLD"}


class UncertaintyError(ValueError):
    """Missing or inconsistent evidence; never a zero uncertainty."""


def _finite(value, label, *, nonnegative=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise UncertaintyError(f"{label}: expected a numeric value")
    if not math.isfinite(value) or (nonnegative and value < 0):
        raise UncertaintyError(f"{label}: invalid value {value!r}")
    return float(value)


def _source(value, label):
    if not isinstance(value, str) or not value.strip():
        raise UncertaintyError(f"{label}: provenance is required")


def pool_digest(indicator_ids):
    """Stable physical transition/indicator IDs, not rounded wavelengths alone."""
    ids = list(indicator_ids)
    if not ids or any(not isinstance(x, str) or not x.strip() for x in ids):
        raise UncertaintyError("a nonempty physical indicator pool is required")
    if len(set(ids)) != len(ids):
        raise UncertaintyError("duplicate indicator IDs would count evidence twice")
    return hashlib.sha256(json.dumps(sorted(ids), separators=(",", ":")).encode()).hexdigest()


def product_scope(product, *, star, indicator_ids):
    from pipeline.product_eligibility import key_of
    _source(star, "star")
    return {"star": star, "product_key": key_of(product),
            "pool_sha256": pool_digest(indicator_ids)}


def component(name, sigma_dex, *, state, source, evidence=None):
    """Use the existing Term semantics; values here must not be divided by N again."""
    if state not in STATES:
        raise UncertaintyError(f"{name}: unknown applicability state {state!r}")
    _source(source, name)
    if state in {"HOLD", "N/A"}:
        if sigma_dex is not None:
            raise UncertaintyError(f"{name}: {state} cannot carry a numeric value")
    else:
        sigma_dex = _finite(sigma_dex, name, nonnegative=True)
    term = Term(name, sigma_dex, False, source, applicable=state != "N/A")
    return {"name": term.name, "state": state, "sigma_dex": term.dex,
            "source": term.source, "evidence": dict(evidence or {})}


def covariance_variance(jacobian, covariance):
    """J C Jᵀ, with finite, symmetric, positive-semidefinite covariance required."""
    j = np.asarray(jacobian, dtype=float)
    c = np.asarray(covariance, dtype=float)
    if j.ndim != 1 or c.shape != (len(j), len(j)) or not len(j):
        raise UncertaintyError("covariance dimensions do not match the Jacobian")
    if not np.isfinite(j).all() or not np.isfinite(c).all():
        raise UncertaintyError("nonfinite covariance or sensitivity")
    if not np.allclose(c, c.T, rtol=1e-12, atol=1e-15):
        raise UncertaintyError("covariance is not symmetric")
    tolerance = max(float(np.max(np.abs(c))), 1e-30) * 1e-12
    if float(np.linalg.eigvalsh(c).min()) < -tolerance:
        raise UncertaintyError("covariance is not positive semidefinite")
    return max(float(j @ c @ j), 0.0)


def scatter_measurement(values, *, indicator_ids, source, independent):
    """SE of vetted independent lines; correlated or single-line pools remain HOLD."""
    ids, vals = list(indicator_ids), list(values)
    digest = pool_digest(ids)
    if len(ids) != len(vals):
        raise UncertaintyError("line IDs and abundance values have different lengths")
    vals = [_finite(x, "line abundance") for x in vals]
    n = len(vals)
    raw = float(np.std(vals, ddof=1)) if n > 1 else None
    evidence = {"method": "line_scatter", "raw_sigma": raw, "n_lines": n,
                "pool_sha256": digest, "independent": independent,
                "small_n": n < 4}
    valid = independent is True and n > 1
    return component("measurement", raw / math.sqrt(n) if valid else None,
                     state="MEASURED" if valid else "HOLD", source=source,
                     evidence=evidence)


def fit_measurement(sigma_dex, *, source, likelihood, correlation_treatment, n_pixels):
    """Fit/profile likelihood sigma; pixel count is never a line-count divisor."""
    _source(likelihood, "likelihood")
    _source(correlation_treatment, "pixel correlation treatment")
    if not isinstance(n_pixels, int) or isinstance(n_pixels, bool) or n_pixels < 1:
        raise UncertaintyError("positive pixel count required")
    return component("measurement", sigma_dex, state="MEASURED", source=source,
                     evidence={"method": "profile_likelihood", "likelihood": likelihood,
                               "correlation_treatment": correlation_treatment,
                               "n_pixels": n_pixels})


def transition_data(indicator_ids, sigma_dex, weights, *, covariance, sources,
                    covariance_source):
    """Propagate line-specific transition errors, including common-source covariance."""
    ids = list(indicator_ids)
    digest = pool_digest(ids)
    sigmas = [_finite(x, "transition sigma", nonnegative=True) for x in sigma_dex]
    w = [_finite(x, "line response weight") for x in weights]
    if not (len(ids) == len(sigmas) == len(w) == len(sources)):
        raise UncertaintyError("transition data arrays have different lengths")
    for s in sources:
        _source(s, "line-specific transition data")
    _source(covariance_source, "transition covariance")
    variance = covariance_variance(w, covariance)
    if not np.allclose(np.diag(covariance), np.square(sigmas), rtol=1e-10, atol=1e-15):
        raise UncertaintyError("transition covariance disagrees with line-specific sigmas")
    return component("transition_data", math.sqrt(variance), state="MEASURED",
                     source=covariance_source,
                     evidence={"indicator_ids": ids, "pool_sha256": digest,
                               "line_sigma_dex": sigmas, "weights": w,
                               "line_sources": list(sources), "covariance": covariance})


def paired_response(nominal, minus, plus, *, delta, source, parameter_source):
    """Symmetric ±1σ responses on the exact nominal indicator pool.

    Preserve both sides to expose nonlinearity. A moved/missing pool is HOLD at
    the caller, not an automatic restriction to whichever lines survived.
    """
    _source(source, "perturbation runs")
    _source(parameter_source, "stellar parameter uncertainty")
    delta = _finite(delta, "adopted parameter sigma", nonnegative=True)
    if delta <= 0:
        raise UncertaintyError("zero-sigma parameters need a sourced DEFINED zero")
    if set(nominal) != set(minus) or set(nominal) != set(plus):
        raise UncertaintyError("perturbed indicator pools differ from the nominal pool")
    digest = pool_digest(nominal)
    lo, hi = [], []
    for key in sorted(nominal):
        a = _finite(nominal[key], "nominal abundance")
        lo.append(_finite(minus[key], "minus abundance") - a)
        hi.append(_finite(plus[key], "plus abundance") - a)
    # Linearized signed one-sigma response. Side responses remain available for
    # adjudication; no automatic claim that asymmetric responses are adequate.
    response = float(np.mean(np.subtract(hi, lo)) / 2)
    return {"signed_response_dex": response, "delta_minus_dex": float(np.mean(lo)),
            "delta_plus_dex": float(np.mean(hi)), "delta_parameter": delta,
            "jacobian": response / delta, "pool_sha256": digest,
            "per_indicator_minus": lo, "per_indicator_plus": hi,
            "indicator_ids": sorted(nominal), "source": source,
            "parameter_source": parameter_source}


#: RYA-1251 (Ryan 2026-10-28): the PUBLISHED sigma is Asplund, Amarsi & Grevesse 2021's
#: (Sect. 2.1): line-to-line statistics (+) half the 3D effect (+) the NLTE term (half the
#: correction, floor 0.03; none for molecules, as Amarsi+2021). Atomic data and the
#: observational terms are NOT in it -- Asplund: "in general ... not attempted". The other
#: components stay on the document as evidence; they no longer widen the headline bar.
ASPLUND_COMPONENTS = ("measurement", "model_atmosphere", "nlte")
ASPLUND_BASIS = ("Asplund, Amarsi & Grevesse 2021 Sect. 2.1: statistical (+) half the 3D "
                 "effect (+) NLTE (half the correction, floor 0.03); molecules: Amarsi+2021 "
                 "compositional term instead of NLTE, trend term not priced; atomic data and "
                 "observational terms not included")


#: Asplund's NLTE floor. It applies to every atomic result, including one whose NLTE is inside
#: the synthesis (ENGINE-B) or not applied at all (an LTE product has no smaller NLTE error).
ASPLUND_NLTE_FLOOR = 0.03
_NUMERIC = {"MEASURED", "DEFINED"}


def asplund_terms(document, molecule=None, indicator=None):
    """(statistical, systematic) in the Asplund/Amarsi convention, or None.

    statistical = the measurement term (uncorrelated between lines/groups); systematic =
    half the 3D effect (+) NLTE (half the correction, floor 0.03), and for a molecule
    Amarsi+2021's compositional term (our `molecular_coupling`) in place of NLTE. Amarsi
    Eq. 3 treats the systematic part as perfectly correlated between groups, so a combiner
    averages it instead of shrinking it. Amarsi's "trend" term is not priced here.
    For a named molecule the 3D term is at least half Amarsi+2021 Table 2's 3D - MARCS shift
    (Asplund's "half the 3D effect"); a larger measured term is kept.
    RYA-1317: when `indicator` sits on one of Amarsi's groups (CH A-X for the G band, not the
    CH species row) the group's own shift is used, and the term is also at least Amarsi Eq.
    1's atmosphere term (the standard error over 3D, <3D>, ATMO, HM74) for that group.
    None when the measurement is not numeric or any term used is on HOLD.
    """
    terms = {c["name"]: c for c in (document or {}).get("components", [])}
    if terms.get("measurement", {}).get("state") not in _NUMERIC:
        return None
    molecular = terms.get("molecular_coupling", {}).get("state") in _NUMERIC
    used = ASPLUND_COMPONENTS + (("molecular_coupling",) if molecular else ())
    if any(terms.get(n, {}).get("state") == "HOLD" for n in used):
        return None
    def value(name):
        c = terms.get(name, {})
        return c["sigma_dex"] if c.get("state") in _NUMERIC else 0.0
    third = value("molecular_coupling") if molecular else max(value("nlte"), ASPLUND_NLTE_FLOOR)
    atmosphere = value("model_atmosphere")
    evidence = terms.get("model_atmosphere", {}).get("evidence") or {}
    if molecular and evidence.get("bound"):
        # A RYA-587 bound charges the FULL published 1D->3D shift; Asplund charges half of
        # the 3D effect, so the measured grid leg stands and the half-shift floor follows.
        atmosphere = abs(float(evidence.get("measured_leg_dex") or 0.0))
    group = None
    if molecular and indicator:
        from pipeline import molecular_model_form
        group = molecular_model_form.group_terms(indicator)
    if group is not None and group["atmos_sem_dex"] is not None:
        atmosphere = max(atmosphere, abs(group["shift_3d_dex"]) / 2.0, group["atmos_sem_dex"])
    elif molecular and molecule:
        from pipeline import molecular_model_form
        try:
            shift = float(molecular_model_form._row(molecule).model_form_3D_minus_MARCS)
            atmosphere = max(atmosphere, abs(shift) / 2.0)
        except molecular_model_form.ModelFormError:
            pass
    return value("measurement"), math.sqrt(atmosphere ** 2 + third ** 2)


def molecule_of(product):
    """The molecule a product measures, from its first indicator id (e.g. 'CH_Gband:...')."""
    ids = product.get("uncertainty_indicator_ids") or [""]
    return str(ids[0]).split(":")[0].split("_")[0] or None


ASPLUND_STAMPS = ("sigma_asplund", "sigma_asplund_stat", "sigma_asplund_syst",
                  "sigma_asplund_basis")
ASPLUND_LEGACY_3DNLTE_BASIS = (
    "Asplund, Amarsi & Grevesse 2021 Sect. 2.1 on a 3D-NLTE product with no term-by-term "
    "budget: the product's measured statistical term (+) half the mean per-line 3D-NLTE "
    "correction from its own per-line file (+) the 0.03 NLTE floor; atomic data and "
    "observational terms not included")


def legacy_3dnlte_asplund_terms(product, root=None):
    """(stat, syst) for a legacy 3D-NLTE row, from the row and its own per-line file.

    Legacy rows carry no component budget, but a 3D-NLTE row's artifact records each
    line's 3D-NLTE correction (`aberr`). The engine returns 3D and NLTE as ONE correction,
    so half of its mean stands for the model term and the NLTE term takes the 0.03 floor,
    the same split the budgeted ENGINE-A rows carry. None when anything is missing.
    """
    import csv
    from pathlib import Path
    if "uncertainty" in product or product.get("treatment") != "ENGINE-A-3DNLTE":
        return None
    stat = product.get("sigma_stat")
    path = (product.get("provenance") or {}).get("copied_to") or ""
    if stat is None or not path.endswith("_products.csv"):
        return None
    per_line = Path(root or Path(__file__).resolve().parents[1]) / path.replace(
        "_products.csv", "_per_line.csv")
    try:
        rows = list(csv.DictReader(per_line.open()))
    except OSError:
        return None
    corrections = [float(r["aberr"]) for r in rows
                   if r.get("in_domain") == "True" and r.get("status", "ok") == "ok"
                   and r.get("a_3dnlte") not in ("", None) and r.get("aberr") not in ("", None)]
    if not corrections:
        return None
    half = abs(sum(corrections) / len(corrections)) / 2.0
    return float(stat), math.hypot(half, ASPLUND_NLTE_FLOOR)


def legacy_asplund_stamp_problems(product, legacy_row, *, root=None):
    """None unless the row differs from its predecessor ONLY in the Asplund stamps.

    Then [] when every stamp equals the recompute from the row's own artifacts, else the
    mismatches. The stamps are derived, never asserted, so a recompute is their evidence.
    """
    changed = {k for k in set(product) | set(legacy_row)
               if k not in LEGACY_RIDE_ALONG_STAMPS and product.get(k) != legacy_row.get(k)}
    if not changed or not changed <= set(ASPLUND_STAMPS):
        return None
    terms = legacy_3dnlte_asplund_terms(product, root)
    if terms is None:
        return ["no per-line 3D-NLTE artifact to recompute the Asplund stamps from"]
    want = {"sigma_asplund": math.hypot(*terms), "sigma_asplund_stat": terms[0],
            "sigma_asplund_syst": terms[1], "sigma_asplund_basis": ASPLUND_LEGACY_3DNLTE_BASIS}
    return [f"{k} {product.get(k)!r} != recomputed {v!r}"
            for k, v in want.items() if product.get(k) != v]


def asplund_sigma(document, molecule=None, indicator=None):
    """The Asplund-convention sigma: statistical (+) systematic from `asplund_terms`, or None."""
    got = asplund_terms(document, molecule, indicator)
    return None if got is None else math.hypot(*got)


def assemble(scope, components, *, covariance, covariance_source, assumptions,
             required_extra=()):
    """One authoritative total. A required HOLD makes sigma_reported null.

    covariance is the full component covariance in dex², in the order of the
    numeric components. Explicit diagonal matrices document independence rather
    than guessing it. Signed stellar/shared-scale responses belong in this matrix.
    """
    _source(covariance_source, "covariance")
    _source(assumptions, "covariance/independence assumptions")
    records = [component(**dict(c)) for c in components]
    names = [c["name"] for c in records]
    if len(set(names)) != len(names):
        raise UncertaintyError("duplicate components would double-count uncertainty")
    missing = sorted(set(COMPONENTS).union(required_extra) - set(names))
    for name in missing:
        records.append(component(name, None, state="HOLD",
                                 source="RYA-587: required evidence not supplied"))
    for key in ("star", "product_key", "pool_sha256"):
        _source(scope.get(key), f"scope.{key}")
    numeric = [c for c in records if c["state"] in {"MEASURED", "DEFINED"}]
    if numeric:
        variance = covariance_variance(np.ones(len(numeric)), covariance)
        diagonal = np.diag(np.asarray(covariance, dtype=float))
        if not np.allclose(diagonal, [c["sigma_dex"] ** 2 for c in numeric],
                           rtol=1e-10, atol=1e-15):
            raise UncertaintyError("covariance diagonal disagrees with component variances")
    else:
        if covariance != []:
            raise UncertaintyError("covariance supplied without numeric components")
        variance = None
    holds = [c["name"] for c in records if c["state"] == "HOLD"]
    return {"schema": SCHEMA, "scope": dict(scope), "components": records,
            "required_extra": list(required_extra),
            "covariance": covariance, "covariance_order": [c["name"] for c in numeric],
            "covariance_source": covariance_source, "assumptions": assumptions,
            "state": "HOLD" if holds else "COMPLETE", "holds": holds,
            "sigma_reported": math.sqrt(variance) if not holds and variance is not None else None,
            "sigma_measured_partial": math.sqrt(variance) if variance is not None else None}


def validate(document, *, scope):
    """Recompute totals and verify identity; never trust a caller's COMPLETE flag."""
    if not isinstance(document, Mapping) or document.get("schema") != SCHEMA:
        raise UncertaintyError("missing canonical RYA-587 uncertainty document")
    if document.get("scope") != scope:
        raise UncertaintyError("uncertainty belongs to a different star/product/indicator pool")
    rebuilt = assemble(scope, document["components"], covariance=document["covariance"],
                       covariance_source=document["covariance_source"],
                       assumptions=document["assumptions"],
                       required_extra=document.get("required_extra", ()))
    for field in ("state", "holds", "covariance_order"):
        if document.get(field) != rebuilt[field]:
            raise UncertaintyError(f"stale or inconsistent {field}")
    # RYA-1233: a recomputed quadrature sum differs in the last bit between runs (BLAS
    # summation order), so the total is compared to rounding, not bit-for-bit.
    got, want = document.get("sigma_reported"), rebuilt["sigma_reported"]
    if (got is None) != (want is None) or (
            want is not None and not math.isclose(got, want, rel_tol=1e-12, abs_tol=1e-15)):
        raise UncertaintyError("stale or inconsistent sigma_reported")
    if rebuilt["holds"]:
        raise UncertaintyError("unresolved components: " + ", ".join(rebuilt["holds"]))
    terms = {c["name"]: c for c in rebuilt["components"]}
    for name in ("measurement", "transition_data", "stellar.teff", "stellar.logg"):
        if terms[name]["state"] == "N/A":
            raise UncertaintyError(f"{name} cannot be omitted from an abundance budget")
    measurement = terms["measurement"]["evidence"]
    method = measurement.get("method")
    if method == "line_scatter":
        if measurement.get("independent") is not True or measurement.get("n_lines", 0) < 2:
            raise UncertaintyError("scatter SE needs vetted independent lines and N > 1")
        if measurement.get("pool_sha256") != scope["pool_sha256"]:
            raise UncertaintyError("measurement uncertainty uses a different pool")
        expected = _finite(measurement.get("raw_sigma"), "raw scatter", nonnegative=True)
        expected /= math.sqrt(measurement["n_lines"])
        if not math.isclose(expected, terms["measurement"]["sigma_dex"], rel_tol=1e-10):
            raise UncertaintyError("measurement sigma is not raw scatter / sqrt(N)")
    elif method == "profile_likelihood":
        _source(measurement.get("correlation_treatment"), "pixel correlation treatment")
        _source(measurement.get("likelihood"), "profile likelihood")
        n_pixels = measurement.get("n_pixels")
        if not isinstance(n_pixels, int) or isinstance(n_pixels, bool) or n_pixels < 1:
            raise UncertaintyError("positive pixel count required")
    else:
        raise UncertaintyError("measurement method is not declared")
    data = terms["transition_data"]["evidence"]
    recomputed = transition_data(data.get("indicator_ids", []), data.get("line_sigma_dex", []),
                                data.get("weights", []), covariance=data.get("covariance", []),
                                sources=data.get("line_sources", []),
                                covariance_source=terms["transition_data"]["source"])
    if (data.get("pool_sha256") != scope["pool_sha256"]
            or recomputed["evidence"]["pool_sha256"] != scope["pool_sha256"]):
        raise UncertaintyError("transition-data uncertainty uses a different indicator pool")
    if not math.isclose(recomputed["sigma_dex"], terms["transition_data"]["sigma_dex"], rel_tol=1e-10):
        raise UncertaintyError("transition-data sigma is inconsistent with its line evidence")
    if method == "line_scatter":
        # RYA-1232: components closer than 0.5 A are ONE resolved feature (fitted in windows
        # that each see the same blend). The scatter is then over features (n_lines) and its
        # n_components must be exactly the physical pool; without n_components the scatter
        # line count must equal the pool, as before.
        n_pool = measurement.get("n_components", measurement["n_lines"])
        if n_pool != len(data["indicator_ids"]) or measurement["n_lines"] > n_pool:
            raise UncertaintyError("scatter line count disagrees with the physical pool")
    for name in ("stellar.teff", "stellar.logg", "stellar.xi", "stellar.metallicity"):
        c = terms[name]
        if c["state"] == "MEASURED":
            ev = c["evidence"]
            if ev.get("pool_sha256") != scope["pool_sha256"]:
                raise UncertaintyError(f"{name}: sensitivity is not measured on this pool")
            _source(ev.get("parameter_source"), f"{name}: parameter source")
            _source(ev.get("response_assessment"), f"{name}: asymmetric response assessment")
            response = _finite(ev.get("signed_response_dex"), f"{name}: response")
            minus = _finite(ev.get("delta_minus_dex"), f"{name}: minus response")
            plus = _finite(ev.get("delta_plus_dex"), f"{name}: plus response")
            if not math.isclose(response, (plus - minus) / 2, rel_tol=1e-10, abs_tol=1e-14):
                raise UncertaintyError(f"{name}: signed response disagrees with perturbations")
            if not math.isclose(abs(response), c["sigma_dex"], rel_tol=1e-10):
                raise UncertaintyError(f"{name}: sigma disagrees with its signed response")
            if _finite(ev.get("delta_parameter"), f"{name}: adopted sigma") <= 0:
                raise UncertaintyError(f"{name}: invalid adopted parameter uncertainty")
        elif c["state"] == "DEFINED" and c["sigma_dex"] != 0:
            raise UncertaintyError(f"{name}: nonzero parameter terms need perturbation evidence")
    return rebuilt


def differential_uncertainty(target, reference, *, matched_ids, target_scope,
                             reference_scope, cross_covariance, source, assumptions):
    """Matched-route target minus reference: propagate the JOINT covariance.

    Cross covariance can cancel a shared gf scale. An explicit zero block is an
    independence assumption, never the default. Absolute budgets remain attached.
    """
    a, b = validate(target, scope=target_scope), validate(reference, scope=reference_scope)
    digest = pool_digest(matched_ids)
    if target_scope["star"] == reference_scope["star"]:
        raise UncertaintyError("target-star differential needs a distinct reference star")
    if any(s["pool_sha256"] != digest for s in (target_scope, reference_scope)):
        raise UncertaintyError("differential requires budgets rederived on the matched pool")
    # Identity contains route, treatment, line_set; holdings/instruments may differ.
    from pipeline.product_eligibility import KEY_FIELDS
    ak, bk = target_scope["product_key"].split("|"), reference_scope["product_key"].split("|")
    for field in ("element", "ion", "route", "treatment", "line_set"):
        idx = KEY_FIELDS.index(field)
        if len(ak) != len(KEY_FIELDS) or len(bk) != len(KEY_FIELDS) or ak[idx] != bk[idx]:
            raise UncertaintyError(f"differential has unmatched {field}")
    _source(source, "cross covariance")
    _source(assumptions, "differential assumptions")
    ac, bc = np.asarray(a["covariance"]), np.asarray(b["covariance"])
    cross = np.asarray(cross_covariance, dtype=float)
    if cross.shape != (len(ac), len(bc)):
        raise UncertaintyError("cross covariance dimensions do not match the absolute budgets")
    joint = np.block([[ac, cross], [cross.T, bc]])
    variance = covariance_variance(np.r_[np.ones(len(ac)), -np.ones(len(bc))], joint)
    return {"schema": SCHEMA, "quantity": "target_minus_reference_dex",
            "target": a, "reference": b, "matched_pool_sha256": digest,
            "cross_covariance": cross.tolist(), "source": source,
            "assumptions": assumptions, "sigma_differential": math.sqrt(variance)}


def ratio_uncertainty(log_abundances, covariance, *, numerator, denominator,
                      source, assumptions, linear_ratio=False):
    """Propagate a joint abundance covariance to [X/Fe] or C/O.

    For [X/Fe], inputs must already be matched-reference differential abundances
    and their joint covariance. For C/O, use absolute A(C), A(O). The shared
    reference and molecular-coupling covariance must remain in the input matrix.
    """
    _source(source, "ratio covariance provenance")
    _source(assumptions, "ratio covariance assumptions")
    values = np.asarray(log_abundances, dtype=float)
    if not (0 <= numerator < len(values) and 0 <= denominator < len(values)) or numerator == denominator:
        raise UncertaintyError("ratio needs two distinct abundance indices")
    if not np.isfinite(values).all():
        raise UncertaintyError("nonfinite abundance ratio input")
    j = np.zeros(len(values))
    j[numerator], j[denominator] = 1.0, -1.0
    var = covariance_variance(j, covariance)
    difference = float(values[numerator] - values[denominator])
    value = 10 ** difference if linear_ratio else difference
    sigma = math.sqrt(var) * (math.log(10) * value if linear_ratio else 1)
    return {"value": value, "sigma": sigma, "jacobian": j.tolist(),
            "covariance": covariance, "source": source, "assumptions": assumptions,
            "basis": "linearized_number_ratio" if linear_ratio else "log_abundance_difference"}


def publication_problems(product):
    """Gate used by every element publisher. Missing evidence is explicit HOLD."""
    try:
        star = product.get("star")
        ids = product.get("uncertainty_indicator_ids")
        scope = product_scope(product, star=star, indicator_ids=ids or [])
        budget = validate(product.get("uncertainty"), scope=scope)
        terms = {c["name"]: c for c in budget["components"]}
        if product.get("band") in {"NIR", "H", "K"} and terms["telluric"]["state"] == "N/A":
            raise UncertaintyError("IR products require quantified telluric residual uncertainty")
        if terms["stellar.xi"]["state"] == "N/A":
            evidence = terms["stellar.xi"]["evidence"]
            if evidence.get("parameter_exists") is not False:
                raise UncertaintyError("xi N/A requires evidence that this route has no xi parameter")
        if str(product.get("selector", "")).startswith("MOL-") and terms["molecular_coupling"]["state"] == "N/A":
            raise UncertaintyError("molecular products require coupling assessment")
        if product.get("sigma_reported") is None or not math.isclose(
                _finite(product["sigma_reported"], "sigma_reported", nonnegative=True),
                budget["sigma_reported"], rel_tol=1e-10, abs_tol=0.00005):
            raise UncertaintyError("published sigma_reported does not match the canonical total")
        if product.get("X_H") is not None or product.get("sigma_differential") is not None:
            differential = product.get("differential_uncertainty")
            if not isinstance(differential, Mapping):
                raise UncertaintyError("target-star abundance needs matched-reference differential covariance")
            reference = differential["reference"]
            result = differential_uncertainty(budget, reference, matched_ids=ids,
                target_scope=scope, reference_scope=reference["scope"],
                cross_covariance=differential["cross_covariance"],
                source=differential["source"], assumptions=differential["assumptions"])
            if differential.get("target") != budget:
                raise UncertaintyError("differential contains a stale target absolute budget")
            if differential.get("sigma_differential") != result["sigma_differential"]:
                raise UncertaintyError("differential total disagrees with joint covariance")
            if product.get("sigma_differential") != result["sigma_differential"]:
                raise UncertaintyError("published differential sigma disagrees with canonical result")
        if any(product.get(k) is not None for k in ("X_Fe", "C_O")):
            raise UncertaintyError("multi-element ratio requires joint component evidence; ratio_uncertainty alone is not publication evidence")
    except (ValueError, TypeError, KeyError) as exc:
        return [str(exc)]
    return []


#: 🔴 RYA-1224 -- THE ONLY TWO FIELDS THE LEGACY COMPARISON IGNORES, AND THE LIST IS
#: CLOSED. Both are readings taken AT EMIT TIME about the emit itself, not statements about
#: the science: `generated_at` is a wall clock and `code_commit` is `git rev-parse HEAD`.
#: `rya1178_emit_fe_schema.enrich` restamps both on every product on every run, so under a
#: whole-dict comparison EVERY row lost its exemption to a ride-along and the feed became
#: unwritable by its own emitter -- measured at all 160 Fe rows, with those two fields the
#: ONLY difference from the committed feed (RYA-1224 EOS section 8).
#:
#: ⚠️ NOTHING ELSE MAY JOIN THIS TUPLE. The exemption's whole job is to refuse a CHANGED
#: product without evidence, so every field that could carry a change -- `A`, `n_lines`,
#: `n_excluded`, any sigma, `xi_state`, any label or note -- must stay inside the
#: comparison. Adding a third name here would blind the guard to a real move, which is the
#: opposite of the defect being fixed. `test_rya1224_xi_floor` asserts the tuple's exact
#: contents and mutation-tests a 0.001 dex nudge to A, so a widening fails loudly.
LEGACY_RIDE_ALONG_STAMPS = ("generated_at", "code_commit")


def _legacy_comparable(product):
    """A product's content with the emit-time stamps dropped, for the legacy comparison."""
    return {k: v for k, v in product.items() if k not in LEGACY_RIDE_ALONG_STAMPS}


#: 🔴 RYA-1224 (Ryan's ruling, 2026-09-18) -- THE xi-LAYER CORRECTION ROUTE. A legacy row may
#: be RE-PUBLISHED with a corrected xi layer, in place of the full canonical budget it has
#: never carried, when the change is provably CONFINED to that layer and the asserting side
#: names the artifact it was measured on.
#:
#: The problem this solves: RYA-1224 made the min_paired floor a property of the rule rather
#: than of the artifacts, which corrected `xi_state`/`sigma_xi` on real products. Those are
#: exactly the fields the legacy exemption must refuse a change to without evidence -- so
#: the more honest number could not be published, while the demonstrably WRONG one (a
#: derivative borrowed from a pool the product does not have) stayed live because it was
#: already in the file. The gate was preserving a known defect.
#:
#: ⚠️ THIS IS NOT A WEAKER GATE, IT IS A NARROWER ONE, AND EVERY CLAUSE IS LOAD-BEARING:
#:   * the changed set must be a SUBSET of the xi layer -- one field outside it and the
#:     route does not apply at all, so `A`, `n_lines`, `n_excluded`, `provenance`, every
#:     label and every other sigma are as protected as before;
#:   * `A`, `n_lines`, `n_excluded` are ALSO checked by name, because they are the RYA-161
#:     firewall and a subset test that silently changed meaning would take them with it;
#:   * a row that ASSERTS a derivative must name an artifact that EXISTS ON DISK -- resolved,
#:     not merely non-empty (RYA-1170: point at the SSOT, never take its word for it);
#:   * a row that WITHDRAWS one needs no artifact, because nothing is being claimed -- but it
#:     must leave no readable derivative behind, which is the RYA-1224 defect itself;
#:   * the sigma arithmetic is RECOMPUTED here. This route skips `validate`, so it cannot
#:     also skip the recompute that makes `sigma_reported` mean something.
XI_LAYER_FIELDS = (
    "xi_state", "sigma_xi", "xi_note", "xi_source", "dA_dxi_dex_per_kms",
    "sigma_syst_complete", "sigma_syst_components",
    "sigma_reported", "sigma_reported_basis", "sigma_reported_caveat",
)

#: The RYA-161 firewall fields, named so a change to them can never ride this route.
XI_IMMUTABLE_FIELDS = ("A", "n_lines", "n_excluded")

#: States in which no derivative is asserted, so no artifact is owed.
XI_HOLD_STATES = ("UNMEASURED", "NOT_IN_CAMPAIGN")


def _quadrature(*terms):
    got = [t for t in terms if t]
    return round(math.sqrt(sum(t * t for t in got)), 6) if got else None


def xi_layer_correction_problems(product, legacy_row, *, root=None):
    """Whether this re-publication qualifies as an evidenced xi-layer correction.

    Returns ``None`` when the route does not apply (the change is not confined to the xi
    layer), ``[]`` when the row qualifies, or a list of reasons when it applies and fails.
    """
    changed = {k for k in set(product) | set(legacy_row)
               if k not in LEGACY_RIDE_ALONG_STAMPS
               and product.get(k) != legacy_row.get(k)}
    if not changed or not changed <= set(XI_LAYER_FIELDS):
        return None                      # not an xi-only change: this route has no opinion

    problems = []
    for field in XI_IMMUTABLE_FIELDS:
        if product.get(field) != legacy_row.get(field):
            problems.append(f"{field} moved on an xi-layer correction -- refused (RYA-161)")

    state, sigma_xi = product.get("xi_state"), product.get("sigma_xi")
    if sigma_xi is not None:
        if state != "MEASURED":
            problems.append(f"sigma_xi is published but xi_state is {state!r}")
        source = product.get("xi_source")
        if not isinstance(source, str) or not source.strip():
            problems.append("an asserted xi derivative must name its artifact in xi_source")
        else:
            from pathlib import Path
            base = Path(root) if root else Path(__file__).resolve().parents[1]
            if not (base / source).is_file():
                problems.append(f"xi_source does not resolve to a file: {source}")
        if product.get("dA_dxi_dex_per_kms") is None:
            problems.append("a published sigma_xi with no readable derivative")
    else:
        if state not in XI_HOLD_STATES:
            problems.append(f"no sigma_xi but xi_state is {state!r}, not a hold")
        #: 🔴 THE RYA-1224 DEFECT ITSELF -- the honest label with the borrowed number still
        #: attached. A hold that leaves the float readable is not a hold.
        if "dA_dxi_dex_per_kms" in product:
            problems.append("a held pool still publishes a readable dA_dxi_dex_per_kms")

    complete = _quadrature(product.get("sigma_syst"), sigma_xi)
    if product.get("sigma_syst_complete") != complete:
        problems.append(f"sigma_syst_complete is {product.get('sigma_syst_complete')!r}, "
                        f"not quadrature(sigma_syst, sigma_xi) = {complete!r}")
    reported = _quadrature(product.get("sigma_stat"), product.get("sigma_syst_complete"))
    if product.get("sigma_reported") != reported:
        problems.append(f"sigma_reported is {product.get('sigma_reported')!r}, not "
                        f"quadrature(sigma_stat, sigma_syst_complete) = {reported!r}")
    return problems


def assert_publication_feed(document, *, previous=None):
    """Require evidence for additions/changes; preserve exact existing legacy rows.

    ``previous`` is the destination read before writing, never a candidate-provided
    exemption. Retention does not certify legacy uncertainty completeness.

    The comparison is on SCIENTIFIC CONTENT: the two emit-time stamps in
    ``LEGACY_RIDE_ALONG_STAMPS`` are excluded, everything else is compared exactly. A clock
    tick is not a product change; anything else is (RYA-1224).
    """
    legacy = (previous or {}).get("products", [])
    #: ⚠️ DICT EQUALITY, NOT A CANONICAL JSON DUMP. Round-tripping through `json.dumps` to
    #: get a hashable key would normalise the very things a change can hide in -- float
    #: repr, int/float identity, key order in nested provenance -- and RYA-1037 was bitten
    #: by exactly that class of normalisation. n is 160; an O(n^2) exact comparison is
    #: cheap and cannot launder a difference.
    legacy_comparable = [_legacy_comparable(q) for q in legacy]
    from pipeline.product_eligibility import key_of
    #: RYA-1224: the xi-layer route needs THIS row's predecessor, not "some equal row", so
    #: legacy is indexed on the RYA-1127 identity key. A duplicated key is left out of the
    #: index entirely -- with two candidates there is no unique predecessor, and guessing
    #: one would be the neighbour-matching error RYA-1206 was built to refuse.
    by_key = {}
    for q in legacy:
        by_key.setdefault(key_of(q), []).append(q)
    legacy_by_key = {k: v[0] for k, v in by_key.items() if len(v) == 1}
    failures = []
    for product in document.get("products", []):
        if "uncertainty" not in product and _legacy_comparable(product) in legacy_comparable:
            continue
        #: RYA-1224 (Ryan, 2026-09-18): a legacy row whose change is CONFINED to the xi
        #: layer, with the asserting side naming a resolvable artifact and the sigma
        #: arithmetic recomputed here. `None` means the route has no opinion and the full
        #: contract decides, exactly as before.
        predecessor = legacy_by_key.get(key_of(product))
        if "uncertainty" not in product and predecessor is not None:
            stamped = legacy_asplund_stamp_problems(product, predecessor)
            if stamped == []:
                continue
            if stamped:
                failures.append(f"{key_of(product)}: Asplund stamp refused: "
                                f"{'; '.join(stamped)}")
                continue
            verdict = xi_layer_correction_problems(product, predecessor)
            if verdict == []:
                continue
            if verdict:
                failures.append(f"{key_of(product)}: xi-layer correction refused: "
                                f"{'; '.join(verdict)}")
                continue
        problems = publication_problems(product)
        if problems:
            failures.append(f"{key_of(product)}: {'; '.join(problems)}")
    if failures:
        raise UncertaintyError("RYA-587 refuses incomplete live products:\n" + "\n".join(failures))
