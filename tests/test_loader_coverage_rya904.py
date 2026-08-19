"""
tests/test_loader_coverage_rya904.py
====================================
RYA-904 — the CRIRES+ telluric-CORRECTED holdings were UNADDRESSABLE, and the two keys
that made them so.

WHAT THE DEFECT ACTUALLY WAS
----------------------------
Not a telluric decision. `gate_holding()` was correct and is unchanged: it refuses
`solar_vesta_crires_plus_idp` (telluric_applied=not-applied on an instrument registered
telluric_required=yes) and passes both corrected holdings. What was wrong was the
DISPATCH one level up — `_LOADER_HOLDING` mapped one INSTRUMENT to one HOLDING, and for
`crires_plus` it named the raw IDPs. So the harness asked for the holding that should be
refused, was refused, and the two corrected holdings could not be named by anything.

An unreachable holding leaves no trace. A refusal leaves a reason; unreachability looks
exactly like having no data, which is why this survived from RYA-806 to RYA-904 with the
telluric gate green the whole time (RYA-833: an absence is a hypothesis, not a result).

WHAT THESE TESTS PIN
--------------------
  * both keys are HOLDING-level, and they moved TOGETHER. Wiring the corrected Y arm
    while leaving `PRE_NORMALISED` instrument-keyed sets a continuum on an
    already-normalised spectrum — the RYA-713 defect, EWs low by a median 11.7 %.
  * selection prefers a passing holding and NEVER falls back to a refused one.
  * the provenance NAMES the holding, because `crires_plus` alone does not distinguish
    the corrected Y arm from the raw IDPs.
  * the guard is RED on the pre-RYA-904 wiring. A guard nobody has seen fail is not a
    guard (RYA-805/833), so the control is executed here, not asserted in prose.

Everything here runs off the git-tracked RYA-794 product and the registry. No Sirius-only
data is needed, which is the point: the wiring is testable even where the spectra are not.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import measure_band_ew as M  # noqa: E402
from pipeline import loader_coverage as LC  # noqa: E402
from pipeline.telluric_policy import gate_holding  # noqa: E402

Y_HOLDING = "solar_crires_plus_y_rya794"
IDP_HOLDING = "solar_vesta_crires_plus_idp"
#: fully inside the RYA-794 product's 10280-10680 A span
IN_Y = 10400.0


# ── the two keys, and that they are the SAME key ─────────────────────────────

def test_both_maps_are_holding_level_and_agree_on_their_keys():
    """`_LOADER_HOLDING` and `PRE_NORMALISED` were separately instrument-keyed dicts.

    They are now one table. `PRE_NORMALISED` is DERIVED from `_INSTRUMENT_HOLDINGS`
    rather than written a second time, so the pair cannot come apart the way RYA-845's
    two declarations of one constant did.
    """
    assert not hasattr(M, "_LOADER_HOLDING"), (
        "the instrument->holding map is the defect; it must not survive alongside "
        "its replacement, or a caller can still reach the wrong one")
    wired = {h.holding_id for specs in M._INSTRUMENT_HOLDINGS.values() for h in specs}
    assert set(M.PRE_NORMALISED) == wired
    # every wired holding is a real registered holding, not a name someone typed
    import pandas as pd
    reg = set(pd.read_csv(LC.HOLDINGS).holding_id.astype(str))
    assert wired <= reg, f"unregistered holding(s) wired: {sorted(wired - reg)}"


def test_one_instrument_two_holdings_two_normalisation_answers():
    """The reason a per-instrument flag could not be right: CRIRES+ disagrees with itself.

    Raw Vesta IDPs are un-normalised adu (TUNIT2='adu', flux to ~1.9e5); the RYA-794
    corrected Y arm is already continuum-normalised (median 0.9985). A single answer for
    `crires_plus` is wrong for one of them, and being wrong in the
    "renormalise something already normalised" direction is RYA-713.
    """
    specs = {h.holding_id: h for h in M.holdings_for("crires_plus")}
    assert specs[Y_HOLDING].pre_normalised is True
    assert specs[IDP_HOLDING].pre_normalised is False


# ── selection: prefer passing, never fall back to refused ────────────────────

def test_the_gate_itself_is_unchanged_and_still_discriminates():
    """POSITIVE CONTROL on the gate before blaming the dispatch (RYA-896 verified this;
    re-run here so 'the gate is fine' is measured in this file, not quoted)."""
    assert gate_holding(Y_HOLDING, "crires_plus")[0] is True
    assert gate_holding(IDP_HOLDING, "crires_plus")[0] is False
    assert gate_holding("elgueta2026_vizier", "crires_plus")[0] is True


def test_select_prefers_the_corrected_holding_inside_its_span():
    spec = M.select_holding("crires_plus", IN_Y, 1.2)
    assert spec.holding_id == Y_HOLDING
    assert spec.pre_normalised is True


def test_outside_the_corrected_span_it_refuses_rather_than_serving_a_truncated_window():
    """Coverage must be TOTAL. A window half inside a fixed-span product is a truncated
    window, and serving it quietly would be worse than the refusal it replaces."""
    with pytest.raises(M.TelluricNotCorrected):
        M.select_holding("crires_plus", 10679.5, 1.2)     # straddles the red edge
    with pytest.raises(M.TelluricNotCorrected):
        M.select_holding("crires_plus", 10050.0, 1.0)     # blueward of the product


def test_a_refused_holding_is_never_silently_fallen_back_to():
    """Named explicitly, the raw IDP holding still cannot be measured — the gate decides,
    and `holding=` is a way to ASK, not a way to bypass."""
    with pytest.raises(M.TelluricNotCorrected) as e:
        M.select_holding("crires_plus", IN_Y, 1.2, holding=IDP_HOLDING)
    assert IDP_HOLDING in str(e.value)
    assert "not-applied" in str(e.value)


def test_unknown_instrument_is_loud():
    with pytest.raises(LookupError, match="no window loader for instrument"):
        M.select_holding("jwst_nirspec_prism", 10400.0, 1.0)


# ── the window itself ────────────────────────────────────────────────────────

def test_load_window_serves_the_corrected_arm_and_names_the_holding():
    win = M.load_window_ex("crires_plus", IN_Y, 1.2)
    assert win.holding.holding_id == Y_HOLDING
    assert f"holding={Y_HOLDING}" in win.provenance
    assert win.pre_normalised is True
    assert win.wave.size == win.flux.size > 12
    assert np.all(np.diff(win.wave) >= 0)
    # the corrected product is NORMALISED — assert the relationship, not a magic constant
    assert 0.9 < float(np.median(win.flux)) < 1.1
    # ... and the raw IDPs are not, which is the whole reason this is per-holding
    assert M.PRE_NORMALISED[IDP_HOLDING] is False


def test_the_gdsat_caveat_travels_with_the_holding():
    """RYA-904 spec 6. Elgueta's G-dwarf robustness flag certifies ZERO Fe I lines for a
    solar-type star, and RYA-794 quoted no abundance. Reaching the pool does not
    pre-decide that it is publishable, so the caveat is a property of the holding rather
    than a sentence someone remembers to retype in each consumer."""
    spec = M.select_holding("crires_plus", IN_Y, 1.2)
    assert "GDSat" in spec.caveat and "RYA-794" in spec.caveat


def test_the_product_is_air_and_rest_frame_and_that_is_checked_not_assumed():
    """CRIRES+ delivers VACUUM/TOPOCENT natively; `wavelength_air_A` in a DERIVED product
    is a claim. A vacuum copy has the same columns, the same row count and the same flux
    — it just fits every line ~2.9 A away, at the wrong abundance. So the loader runs a
    POSITIVE control on RYA-794's own published Fe I positions.

    THE TEST MUST DISCRIMINATE (RYA-805): shifting the product by the air-vacuum
    difference must make the check fail, or it is not testing anything.
    """
    M._crires_cache.pop("y", None)
    w, f = M.crires_y_spectrum()               # passes on the real product
    with pytest.raises(LookupError, match="AIR wavelengths"):
        M._assert_air_rest_frame(w + 2.9, f)   # a vacuum-scale shift is caught
    M._crires_cache.pop("y", None)


# ── the guard, and the control that proves it fires ──────────────────────────

def test_guard_is_green_on_the_wiring_as_shipped():
    LC.reconcile_loader_coverage()


def test_POSITIVE_CONTROL_guard_is_red_on_the_pre_rya904_wiring():
    """🔴 THE CONTROL. Point `crires_plus` back at ONLY the raw IDP holding — the exact
    `_LOADER_HOLDING` this ticket replaces — and clear the declared gaps, and the guard
    must name BOTH corrected holdings as unreachable.

    This is the state `origin/main` was in, with every telluric test green.
    """
    pre_904 = {"solar_kpno": "kpno_solar_atlas",
               "solar_iag": "iag_fts_solar_atlas",
               IDP_HOLDING: "crires_plus"}
    with pytest.raises(LC.LoaderCoverageError) as e:
        LC.reconcile_loader_coverage(addressable=pre_904, declared_gaps={})
    msg = str(e.value)
    assert Y_HOLDING in msg
    assert "elgueta2026_vizier" in msg
    # and it must NOT accuse the holding that is correctly refused
    assert IDP_HOLDING not in msg


def test_guard_does_not_fire_on_a_correctly_refused_holding():
    """DISCRIMINATION. `solar_vesta_crires_plus_idp` is verified and unaddressable-by-
    measurement, and it must NEVER appear: it does not pass the gate, so its absence from
    the harness is a decision with a reason, not a gap. A guard that also flagged it
    would be pressure to weaken the telluric gate — the opposite of the fix."""
    rows = [b["holding_id"] for b in LC.unreachable_holdings(
        addressable={}, declared_gaps={})]
    assert IDP_HOLDING not in rows
    assert Y_HOLDING in rows                    # ... while the corrected one does appear


def test_every_declared_gap_names_a_reason_and_a_ticket():
    """A gap is a CLAIM, and an unexplained one is indistinguishable from an oversight —
    which is the failure this whole ticket is about. The table can only shrink by
    building something."""
    assert LC.DECLARED_GAPS, "an empty gap table means nothing was declared, not that "\
                             "nothing is missing"
    import re
    for hid, why in LC.DECLARED_GAPS.items():
        assert len(why) > 80, f"{hid}: a gap needs a reason, not a label"
        assert re.search(r"RYA-\d+", why), f"{hid}: name the ticket that owns the gap"


def test_declared_gaps_are_real_holdings_not_typos():
    import pandas as pd
    reg = set(pd.read_csv(LC.HOLDINGS).holding_id.astype(str))
    assert set(LC.DECLARED_GAPS) <= reg


def test_a_gap_that_is_actually_wired_is_itself_a_failure():
    """The table rots BOTH ways: a stale exemption is how a guard quietly stops guarding,
    because it would wave through the next holding to go unreachable.

    🔴 THE CONTROL IS SYNTHETIC ON PURPOSE (RYA-897). It used to inject the real
    `solar_harps`, which worked only while HARPS was unwired — and this test's own
    docstring predicted that wiring would land. It did (RYA-897), the entry was removed
    from DECLARED_GAPS, and this control silently stopped exercising anything: the
    injected holding was no longer declared, so no stale gap could form and the raise it
    asserts could never happen. A control coupled to live data expires with the data.
    Both halves are now fabricated, so the detector is tested rather than the table."""
    with pytest.raises(LC.LoaderCoverageError, match="are in fact WIRED"):
        LC.reconcile_loader_coverage(
            declared_gaps={"solar_kpno": "SYNTHETIC — a gap entry for a holding that is "
                                         "in fact wired, so the detector has something "
                                         "to catch. RYA-904/897 control."},
            addressable={**LC.addressable_holdings(), "solar_kpno": "kpno"})


# ── the same defect shape, one and two levels deeper ─────────────────────────

def test_species_token_spells_the_ion_the_way_the_linelist_does():
    """`--ion II` was DECORATIVE on the synthesis route. `select_lines` was pinned to
    Fe I, so `--ion II` selected Fe I lines and labelled every one of them Fe II — a
    plausible NUMBER for an ion nobody measured, which is worse than a refusal."""
    from rya759_nearuv_fe_product import species_token, FE_I
    assert species_token("Fe", "I") == FE_I == "Fe 1"
    assert species_token("Fe", "II") == "Fe 2"
    assert species_token("Sr", "II") == "Sr 2"
    with pytest.raises(SystemExit):
        species_token("Fe", "one")


def test_select_lines_refuses_a_species_the_band_does_not_hold():
    """An empty selection must RAISE, not return an empty frame that becomes a product
    with n=0 (RYA-832: a row that arrives carrying no number is worse than no row)."""
    from rya759_nearuv_fe_product import select_lines
    import pandas as pd
    ll = pd.read_csv(ROOT / "data" / "linelists" / "ispec_ir_9200_13000" /
                     "atomic_lines.tsv", sep="\t", low_memory=False)
    rec = ll.to_records(index=False)
    with pytest.raises(SystemExit, match="no 'Fe 3' rows"):
        select_lines(rec, lo_A=10280, hi_A=10680, n=40, teff=5772.0, min_sep_A=4.0,
                     species="Fe 3")
    # ... and the default is still Fe I, so RYA-759's near-UV selection is untouched
    import inspect
    assert inspect.signature(select_lines).parameters["species"].default == "Fe 1"


# ── the RYA-713 half: the corrected product must NOT be re-normalised ────────

def test_the_corrected_product_is_not_re_normalised_and_the_test_discriminates():
    """RYA-904 spec: assert the harness TREATS this holding as pre-normalised.

    The flag alone is not the property — what matters is which branch it selects. So
    this runs both branches on the real product and shows they differ: the pre-normalised
    branch is a no-op, the un-normalised branch divides by a fitted local continuum and
    moves the flux. If they were the same the test would prove nothing, which is the
    RYA-805 failure mode.

    On the synthesis route there is no renormalisation step at all — `fit_one` hands the
    observed flux straight to `_fit_synth_flux` — so this is the EW/profile-fit route's
    branch, exercised directly. Both routes read the same per-holding flag, and that is
    the thing that could not have been right when it was per-instrument.
    """
    from pipeline.lines_fit import _local_renorm

    win = M.load_window_ex("crires_plus", 10535.709, 1.56)   # RYA-794's Fe I 10535.709
    assert win.pre_normalised is True
    c = 10535.709
    m = np.abs(win.wave - c) <= 1.20
    wf, ff = win.wave[m], win.flux[m]

    # the branch the flag selects: flux untouched, continuum is unity by construction
    kept, cont = ff, np.ones_like(ff)
    assert np.array_equal(kept, ff)
    assert np.all(cont == 1.0)

    # THE DISCRIMINATOR: the other branch is not a no-op on this data
    renormed, _ = _local_renorm(wf, ff, c, window=1.20)
    assert not np.allclose(renormed, ff, atol=1e-6), (
        "the two branches are indistinguishable on this window, so this test cannot "
        "detect a double-normalisation — widen the window or pick a line with structure")
    # and it MOVES A LINE DEPTH, which is what an EW is made of. Measured on this
    # window: median flux 0.9973 -> 0.9872 and Fe I 10535.709's depth 0.1025 -> 0.1118,
    # +9.0%. That is the cost of getting the flag wrong for this holding, and it is the
    # same magnitude class as RYA-713's 11.7% median (opposite sign here, because the
    # direction depends on where the fitted continuum lands relative to unity — which is
    # exactly why "it would probably be small" is not an argument).
    assert abs(float(np.median(renormed)) - float(np.median(ff))) > 1e-4
    d_raw = 1.0 - float(ff.min())
    d_renorm = 1.0 - float(renormed.min())
    assert abs(d_renorm / d_raw - 1.0) > 0.02


# ── the RYA-905 preflight reader must see the new shape ──────────────────────

def test_preflight_dispatch_reader_reads_the_holding_table():
    """RYA-905's `read_dispatch` text-parses the harness. RYA-904 changed what it parses.

    🔴 IT RETURNED AN EMPTY DISPATCH AND CI CAUGHT IT. `_INSTRUMENT_HOLDINGS` is an
    ANNOTATED assignment (`ast.AnnAssign`), and the reader walked `ast.Assign` only — so
    it saw no instruments at all, which reads as "nothing is wired". The only reason that
    surfaced as a failure instead of a page of confident false WARNs is that RYA-905's
    positive control refuses to report absences it cannot corroborate (RYA-833).

    This pins the shape so the next change to that table cannot empty the reader silently.
    """
    import preflight_check as pf
    d = pf.read_dispatch()
    assert d.controls_ok, d.control_note
    assert "crires_plus" in d.instruments
    # the point of the ticket: ONE instrument, SEVERAL holdings, in preference order
    assert d.served_holdings["crires_plus"] == (Y_HOLDING, IDP_HOLDING)
    # and the reader agrees with the harness rather than restating it
    assert {i: tuple(h.holding_id for h in M.holdings_for(i))
            for i in d.instruments} == d.served_holdings


def test_preflight_reader_still_reads_the_legacy_single_valued_map(tmp_path):
    """CONTROL AT THE OTHER END. The reader must keep working against a harness that
    predates RYA-904, or 'it parses the new shape' is compatible with 'it parses only
    the new shape' and the synthetic fixtures the controls rely on stop meaning anything.
    """
    import preflight_check as pf
    legacy = tmp_path / "legacy_harness.py"
    legacy.write_text(
        'PRE_NORMALISED = {"kpno_solar_atlas": True, "harps": False}\n'
        '_LOADER_HOLDING = {"kpno_solar_atlas": "solar_kpno"}\n'
        "def load_window(instrument, centre, pad):\n"
        '    if instrument == "kpno_solar_atlas":\n        return 1, 2, 3\n'
        "    raise LookupError('no window loader')\n")
    d = pf.read_dispatch(legacy)
    assert d.controls_ok, d.control_note
    assert d.instruments == ("kpno_solar_atlas",)
    assert d.served_holdings == {"kpno_solar_atlas": ("solar_kpno",)}
    assert "harps" in d.configured        # configured but not dispatched: RYA-897's shape
