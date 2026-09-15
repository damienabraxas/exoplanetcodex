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
    for field in ("state", "holds", "sigma_reported", "covariance_order"):
        if document.get(field) != rebuilt[field]:
            raise UncertaintyError(f"stale or inconsistent {field}")
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
    else:
        raise UncertaintyError("measurement method is not declared")
    data = terms["transition_data"]["evidence"]
    recomputed = transition_data(data.get("indicator_ids", []), data.get("line_sigma_dex", []),
                                data.get("weights", []), covariance=data.get("covariance", []),
                                sources=data.get("line_sources", []),
                                covariance_source=terms["transition_data"]["source"])
    if data.get("pool_sha256") != scope["pool_sha256"]:
        raise UncertaintyError("transition-data uncertainty uses a different indicator pool")
    if not math.isclose(recomputed["sigma_dex"], terms["transition_data"]["sigma_dex"], rel_tol=1e-10):
        raise UncertaintyError("transition-data sigma is inconsistent with its line evidence")
    for name in ("stellar.teff", "stellar.logg", "stellar.xi", "stellar.metallicity"):
        c = terms[name]
        if c["state"] == "MEASURED":
            ev = c["evidence"]
            if ev.get("pool_sha256") != scope["pool_sha256"]:
                raise UncertaintyError(f"{name}: sensitivity is not measured on this pool")
            _source(ev.get("parameter_source"), f"{name}: parameter source")
            _source(ev.get("response_assessment"), f"{name}: asymmetric response assessment")
            response = _finite(ev.get("signed_response_dex"), f"{name}: response")
            _finite(ev.get("delta_minus_dex"), f"{name}: minus response")
            _finite(ev.get("delta_plus_dex"), f"{name}: plus response")
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


def assert_publication_feed(document):
    """Final write-boundary guard, including metadata updates and schema restamps."""
    failures = []
    for product in document.get("products", []):
        problems = publication_problems(product)
        if problems:
            from pipeline.product_eligibility import key_of
            failures.append(f"{key_of(product)}: {'; '.join(problems)}")
    if failures:
        raise UncertaintyError("RYA-587 refuses incomplete live products:\n" + "\n".join(failures))
