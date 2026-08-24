"""
RYA-985 — the abundance harness follows the star.

The chain was solar-only: `derive_band_products`, `control_synthesis_handler`,
`line_width` and `lte_identity_probe` each read `STAR_PARAMS["solar"]` and there was no
`--star` anywhere. Measuring tau Ceti would have produced a plausible number computed on the
Sun's atmosphere and the Sun's NLTE grid node.

What these pin:

  * every star reaches ITS OWN grid node — the `feh=0.0` that fetched solar departure
    coefficients whatever star was asked for is the plausible-looking-wrong-number failure;
  * `--star solar` is BIT-IDENTICAL to the old hardcoded behaviour;
  * an unknown star LOUD-FAILS — no silent solar fallback;
  * a width bound is never computed from another star's parameters, and a MISSING term is
    refused rather than dropped;
  * a SOLVED parameter is never presented as a fundamental (RYA-957 / Heiter+2015).
"""
from __future__ import annotations

import sys
import types
from pathlib import Path
from unittest import mock

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from config.constants import STAR_PARAMS, get_star_params        # noqa: E402
from pipeline import line_width as LW                            # noqa: E402


def _nodes_for(star: str):
    """Capture the MPIA grid node `engine_a_delta` would request for `star`."""
    import importlib
    mod = importlib.import_module("scripts.derive_band_products")
    seen = []

    def fake_submit(chunk, code, node):
        seen.append(node[0])
        return {node[0]["name"]: {}}

    fake = types.ModuleType("scripts.build_nlte_grids_mpia")
    fake._submit_batch = fake_submit
    with mock.patch.dict(sys.modules, {"scripts.build_nlte_grids_mpia": fake}):
        mod.engine_a_delta("Fe", "I", np.array([5000.0]), star=star)
    return seen[0]


# ── the grid node follows the star ───────────────────────────────────────────────────
@pytest.mark.parametrize("star", ["tau_ceti", "eps_eri", "alpha_cen_a", "alpha_cen_b"])
def test_a_non_solar_star_reaches_its_own_grid_node(star):
    """🔴 The failure this ticket exists to prevent: `feh` was pinned to 0.0, so every
    request fetched the SOLAR departure coefficients whatever star was being measured."""
    n = _nodes_for(star)
    p = get_star_params(star)
    assert n["teff_K"] == int(round(float(p["teff"])))
    assert n["logg"] == pytest.approx(float(p["logg"]))
    assert n["feh"] == pytest.approx(float(p["feh_ref"]))
    assert n["name"] == star, "the batch result is keyed by node name; 'sun' would mis-read it"


def test_tau_ceti_does_not_request_the_solar_node():
    """[Fe/H] = -0.49 is not a rounding difference from 0.00."""
    n = _nodes_for("tau_ceti")
    assert n["feh"] == pytest.approx(-0.49)
    assert n["feh"] != 0.0
    assert n["name"] != "sun"


def test_solar_is_bit_identical_to_the_old_hardcoded_node():
    """The default must change nothing. The old code built exactly this dict."""
    n = _nodes_for("solar")
    p = STAR_PARAMS["solar"]
    assert n == {"name": "sun", "teff_K": int(round(float(p["teff"]))),
                 "logg": float(p["logg"]), "feh": 0.0}


def test_an_unknown_star_loud_fails_with_no_solar_fallback():
    with pytest.raises(KeyError, match="No STAR_PARAMS record"):
        _nodes_for("hd_does_not_exist")


# ── the CLI ──────────────────────────────────────────────────────────────────────────
@pytest.mark.parametrize("script", ["derive_band_products", "control_synthesis_handler"])
def test_the_driver_exposes_star_defaulting_to_solar(script):
    src = (ROOT / "scripts" / f"{script}.py").read_text(encoding="utf-8")
    assert '"--star"' in src, f"{script} has no --star"
    i = src.index('"--star"')
    assert 'default="solar"' in src[i:i + 200], f"{script}'s --star does not default to solar"


def test_no_hardcoded_solar_remains_in_the_measurement_chain():
    for f in ("scripts/derive_band_products.py", "scripts/control_synthesis_handler.py",
              "pipeline/line_width.py", "scripts/lte_identity_probe.py"):
        code = "\n".join(l for l in (ROOT / f).read_text(encoding="utf-8").splitlines()
                         if not l.lstrip().startswith("#"))
        assert 'STAR_PARAMS["solar"]' not in code, f"{f} still hardcodes the solar row"


# ── widths follow the star, and refuse rather than substitute ───────────────────────
def test_the_width_floor_follows_the_star():
    assert LW.irreducible_sigma_kms("tau_ceti") != LW.irreducible_sigma_kms("solar")
    assert LW.irreducible_sigma_kms() == LW.irreducible_sigma_kms("solar")


def test_a_missing_broadening_term_is_refused_not_dropped():
    """🔴 A ceiling with a missing term is an UNDERESTIMATE, and rejecting lines with it
    rejects them for a reason that is not real.

    RYA-988 re-pointed this at the INVARIANT instead of one example star. tau_ceti was
    the original subject; it now carries a cited vmac (Bruntt+2010 Table B1), which made
    the old premise assert fire. Ask stars.yaml which stars actually lack the term rather
    than naming one that the next adoption ticket can quietly fix.
    """
    missing = [s for s in STAR_PARAMS if "vmac" not in get_star_params(s)]
    assert missing, ("no star lacks vmac — build a fixture record rather than deleting "
                     "the guard; the refusal is the behaviour under test")
    for star in missing:
        with pytest.raises(LW.MissingBroadeningTerm, match="vmac"):
            LW.max_stellar_sigma_kms(star)
    LW.max_stellar_sigma_kms("solar")          # the Sun has one, so it still works


def test_a_solved_microturbulence_is_refused_for_the_floor():
    """55 Cnc A carries xi_init/xi_xcheck because it SOLVES xi — not a config constant."""
    assert "xi" not in get_star_params("55cnc_a")
    with pytest.raises(LW.MissingBroadeningTerm, match="xi"):
        LW.irreducible_sigma_kms("55cnc_a")


def test_the_ceiling_uses_the_stars_own_parameters():
    a = LW.physical_ceiling_sigma_A(5500.0, 0.01, "solar")
    b = LW.physical_ceiling_sigma_A(5500.0, 0.01, "alpha_cen_b")
    assert a != b, "the ceiling ignored the star"


# ── pin vs solve ─────────────────────────────────────────────────────────────────────
def test_tau_cetis_logg_is_never_presented_as_a_fundamental():
    """RYA-957: Heiter+2015 brackets it as 'not to be used as a reference'."""
    r = LW.star_parameter_provenance("tau_ceti")
    assert r["logg_is_fundamental"] is False
    assert "logg" in r["solved"] and "logg" not in r["pinned"]
    assert "starting point" in r["note"]


@pytest.mark.parametrize("star", ["solar", "alpha_cen_a", "alpha_cen_b", "eps_eri"])
def test_a_pinned_logg_is_reported_as_fundamental(star):
    r = LW.star_parameter_provenance(star)
    assert r["logg_is_fundamental"] is True
    assert "PINNED fundamental" in r["note"]


def test_the_provenance_carries_the_cited_source():
    r = LW.star_parameter_provenance("tau_ceti")
    assert "Heiter" in r["source"] or "GBS" in r["source"]
    assert r["logg_basis"]
