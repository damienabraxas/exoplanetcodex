"""The ⟨3D⟩ atmosphere's interface contract with iSpec's Turbospectrum path.

The register has carried "the ⟨3D⟩ band product is still owed, and the remaining step is a
PHYSICS one, not plumbing" since v117, for Al and then Fe. It is neither: it is an
INTERFACE CONTRACT, and every clause is readable in iSpec's own source. What made it look
like physics is that the obvious route — let iSpec write a MARCS file — is genuinely
impossible, because a mul23 ⟨3D⟩ model has five columns and the MARCS writer needs Depth,
Pₑ, P_g and P_rad, none of which it carries.

These tests pin the three numbers in `ispec/synth/turbospectrum.py` that the contract turns
on. Each one is a SWITCH whose wrong value produces a well-formed spectrum for a structure
nobody chose:

    :114  is_marcs_model = len(atmosphere_layers[0]) == 11   -> babsma's MARCS-FILE flag
    :159  interpolated_tau = atmosphere_layers[:, 7]         -> the departure file's τ
    :331  if remove_tmp_atm_file: os.remove(...)             -> a caller's file is SAFE

They run against the real committed ⟨3D⟩ model. The deck binary is Sirius-only; the
deck-vs-atmosphere τ agreement quoted below was measured against it.
"""
from pathlib import Path

import numpy as np
import pytest

from pipeline import mean3d_atmosphere as M

MOD = (Path(__file__).resolve().parents[1] / "data" / "atmospheres"
       / "stagger_avg3d_rya442" / "sun_avg3d_stagger.mod")


def _numeric_rows(text: list[str]) -> list[int]:
    """Indices of the real depth rows.

    ⚠️ NOT "lines with five tokens" — the model's own credit line,
    `* Ekaterina Semenova, July 2020.`, has exactly five. The parser is immune because it
    float-converts, but a test fixture built on the token count alone silently slices the
    `TAU5000 SCALE` header away and then asserts against the wrong refusal."""
    out = []
    for i, ln in enumerate(text):
        parts = ln.split()
        if len(parts) != 5:
            continue
        try:
            [float(x) for x in parts]
        except ValueError:
            continue
        out.append(i)
    return out


@pytest.fixture(scope="module")
def loaded():
    return M.load(MOD)


# ── the model parses, and refuses rather than guesses ────────────────────────

def test_the_committed_model_reads(loaded):
    layers, model = loaded
    assert model["model_id"] == "t5777g44m00"
    assert (model["teff"], model["logg"]) == (5777.0, 4.44)
    assert model["feh"] == 0.0
    assert model["ndep"] == 101
    assert model["data"].shape == (101, 5)


def test_the_declared_NDEP_is_checked_against_what_was_read(tmp_path):
    """A truncated model still synthesises — it just has fewer layers than it says. That
    is why the header's own count is a check and not decoration."""
    text = MOD.read_text().splitlines()
    rows = _numeric_rows(text)
    short = tmp_path / "short.mod"
    short.write_text("\n".join(text[:rows[0]] + [text[i] for i in rows[:50]]) + "\n")
    with pytest.raises(M.Mean3DAtmosphereError, match="NDEP=101 but 50 depth rows"):
        M.read_mul23(short)


def test_a_rosseland_scaled_model_is_REFUSED(tmp_path):
    """🔴 τ_ROSSELAND AND τ₅₀₀ ARE DIFFERENT DEPTH VARIABLES. RYA-1013 measured how far
    apart they run on this very grid: the STAGGER cube trimmed on Rosseland has a median
    column top of −7.28 against τ₅₀₀'s −2.64. Pairing τ₅₀₀ departures with a Rosseland
    atmosphere misplaces every layer and still returns a spectrum."""
    swapped = tmp_path / "ross.mod"
    swapped.write_text(MOD.read_text().replace("TAU5000 SCALE", "TAUROSS SCALE"))
    with pytest.raises(M.Mean3DAtmosphereError, match="not 'TAU5000 SCALE'"):
        M.read_mul23(swapped)


def test_non_monotonic_depth_is_REFUSED(tmp_path):
    """Every downstream index — the departure pairing above all — assumes the ordering."""
    text = MOD.read_text().splitlines()
    rows = _numeric_rows(text)
    text[rows[5]], text[rows[6]] = text[rows[6]], text[rows[5]]
    bad = tmp_path / "unsorted.mod"
    bad.write_text("\n".join(text) + "\n")
    with pytest.raises(M.Mean3DAtmosphereError, match="not strictly increasing"):
        M.read_mul23(bad)


# ── the three iSpec switches ─────────────────────────────────────────────────

def test_the_array_must_not_look_like_a_MARCS_model(loaded):
    """🔴 ELEVEN COLUMNS IS A SWITCH, NOT A SIZE. `is_marcs_model = len(layers[0]) == 11`
    becomes babsma's `MARCS-FILE` flag. RYA-442 measured `.true.` on a TAU5000 model as a
    LOUD failure ("This model is probably not a MARCS model!") — but loud in babsma's log,
    several layers below anything that reports a result."""
    layers, _ = loaded
    assert layers.shape[1] != 11
    M.assert_not_marcs_shaped(layers)
    with pytest.raises(M.Mean3DAtmosphereError, match="11 columns"):
        M.assert_not_marcs_shaped(np.zeros((101, 11)))


def test_column_7_carries_log_tau500(loaded):
    """iSpec writes `atmosphere_layers[:, 7]` into the departure file as its τ
    (turbospectrum.py:159), matching the MARCS writer's own map (atmospheres.py:274)."""
    layers, model = loaded
    assert M.TAU_COLUMN == 7
    assert np.allclose(layers[:, 7], model["data"][:, 0])
    assert layers[0, 7] == -5.0 and layers[-1, 7] == 5.0


def test_too_few_columns_is_REFUSED():
    """A 5-column array would be the honest shape of the data and an INDEX ERROR at
    line 159 — or worse, a silent wrap on a different container type."""
    with pytest.raises(M.Mean3DAtmosphereError, match="at least 8"):
        M.assert_not_marcs_shaped(np.zeros((101, 5)))


def test_the_model_is_plane_parallel(loaded):
    """`radius = layers[0][-1]`, and spherical needs `nvalues == 11 and radius > 2.0`. With
    ncols ≠ 11 the model is plane-parallel by construction, which is right for ⟨3D⟩ solar —
    but the last column is set explicitly rather than left to chance."""
    layers, _ = loaded
    assert layers[0, -1] == 0.0


def test_every_column_iSpec_does_not_read_is_NaN(loaded):
    """🔴 THE ARRAY IS NOT THE ATMOSPHERE. babsma reads the real structure from the mul23
    FILE; this array exists only to satisfy three reads. A zero in the unused columns would
    look like a physical value — NaN makes the first caller that mistakes it fail loudly
    instead of quietly using a temperature of 0."""
    layers, _ = loaded
    unused = [c for c in range(layers.shape[1])
              if c not in (M.TAU_COLUMN, layers.shape[1] - 1)]
    assert np.isnan(layers[:, unused]).all()


# ── the gate that makes iSpec's τ overwrite honest ───────────────────────────

def test_tau_consistency_passes_on_the_models_own_tau(loaded):
    """The identity case: a deck whose τ IS this model's. Measured against the real Fe
    ⟨3D⟩ deck on Sirius, max |Δ log τ| = 5.0e-5 — the model file's four-decimal print
    precision (deck −3.09995 written as −3.0999), not a physical difference."""
    layers, model = loaded
    same = {"tau": model["data"][:, 0], "ndep": model["ndep"]}
    assert M.assert_tau_consistent(same, model) == 0.0


def test_tau_consistency_tolerates_the_files_print_precision(loaded):
    """⚠️ THE TOLERANCE IS DERIVED, NOT CHOSEN. Four decimals in the file against float64
    in the deck means a half-digit (5e-5) disagreement is expected and means nothing."""
    _, model = loaded
    tau = model["data"][:, 0].copy()
    tau[19] += 5.0e-5                      # exactly the residual measured on the real deck
    assert M.assert_tau_consistent({"tau": tau, "ndep": len(tau)}, model) == pytest.approx(5e-5)


def test_a_REGRIDDED_deck_is_caught(loaded):
    """🔴 WHAT THE GATE IS ACTUALLY FOR. The deck's spacing is 0.1 dex; a deck built on a
    shifted or different grid disagrees by ~that, three orders above the print precision.
    Without the gate iSpec's overwrite would apply those departures at the wrong depths and
    return a perfectly well-formed spectrum."""
    _, model = loaded
    tau = model["data"][:, 0] + 0.1
    with pytest.raises(M.Mean3DAtmosphereError, match="DIFFERENT depth scales"):
        M.assert_tau_consistent({"tau": tau, "ndep": len(tau)}, model)


def test_a_depth_count_mismatch_is_caught_by_name(loaded):
    """iSpec pairs departures and atmosphere by INDEX, so unequal lengths are not a
    truncation, they are a relabelling of every layer."""
    _, model = loaded
    tau = model["data"][:50, 0]
    with pytest.raises(M.Mean3DAtmosphereError, match="depth count differs"):
        M.assert_tau_consistent({"tau": tau, "ndep": 50}, model)


def test_the_loaded_array_satisfies_the_gerber_depth_check(loaded):
    """The two guards must agree: `gerber_nlte.assert_depth_match` compares the deck's ndep
    against `len(atmosphere)`, and the Fe/Al ⟨3D⟩ decks are 101."""
    from pipeline import gerber_nlte as G
    layers, model = loaded
    G.assert_depth_match({"ndep": 101}, layers)
    with pytest.raises(G.GerberDeckError, match="depth mismatch"):
        G.assert_depth_match({"ndep": 56}, layers)      # a MARCS deck against ⟨3D⟩ layers
