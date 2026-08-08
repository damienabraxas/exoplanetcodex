"""
tests/test_isotope_gf_convention_rya684.py
==========================================
RYA-684 guards for the isotope-fraction double-application.

Three layers, because the finding has three separable failure modes:

  A. The MEASUREMENT can rot — the audit script must still be able to tell a
     fraction-folded list from a fraction-free one.  Tested on synthetic lists
     with the answer built in, so it runs anywhere, Sirius or Mac.

  B. The RECORD can drift — the committed audit record is what the guard reads,
     so its verdicts are pinned.  If a re-vendored line list flips a species,
     this fails and the flip has to be looked at, not absorbed.

  C. The CONVENTION can be reintroduced — a future harness that fits a target
     species against a double-folded block must be stopped at preflight.

The physics constant that anchors all of it: Eu II's two isotopes are 0.478 /
0.522, so a doubly-applied fraction weakens the blended feature by
-log10(0.478^2 + 0.522^2) = 0.3002 dex, and RYA-565's two-leg VALD-vs-GES
comparison measured +0.300.
"""
import json
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

from pipeline import isotope_gf_convention as igc  # noqa: E402
import rya684_isotope_gf_audit as audit  # noqa: E402

RECORD = ROOT / 'data' / 'audit' / 'rya684_isotope_gf_audit.json'

# Turbospectrum's makeabund.f defaults for the species this ticket turns on.
EU = {151: 0.478, 153: 0.522}
BA = {130: 0.00106, 132: 0.00101, 134: 0.02417, 135: 0.06592,
      136: 0.07854, 137: 0.1123, 138: 0.7170}


def _write_list(path, Z, ion, per_iso_loggf, wave=6645.100, ep=1.380):
    """Write a minimal Turbospectrum line list: one feature, one line per isotope."""
    with open(path, 'w') as fh:
        for iso, loggf in sorted(per_iso_loggf.items()):
            code = f"{Z}.{iso:03d}" if iso else f"{Z}.000"
            fh.write(f"'  {code:<18s}' {ion:4d} {1:9d}\n")
            fh.write(f"'Test {Z} {ion}'\n")
            fh.write(f"  {wave:.4f}  {ep:.4f} {loggf:6.3f}   -7.500    4.0  1.00E+08 "
                     f"'x' 'x'   0.0    1.0 'TEST LS:a LS:b'\n")
    return str(path)


# ── A. the measurement itself ────────────────────────────────────────────────

def test_audit_calls_a_fraction_folded_list_folded(tmp_path):
    """gf carrying the isotope fraction reconstructs to one physical gf."""
    phys = -1.000
    p = _write_list(tmp_path / 'folded.list', 63, 2,
                    {a: phys + math.log10(f) for a, f in EU.items()})
    res = audit.audit_species(audit.read_ts_lists([p]), {(63, a): f for a, f in EU.items()},
                              'test')
    assert len(res) == 1
    assert res[0]['verdict'] == 'FOLDED'
    assert res[0]['components_colocated'] is True
    assert res[0]['offset_if_folded_dex'] == pytest.approx(0.3002, abs=5e-4)


def test_audit_calls_a_fraction_free_list_fraction_free(tmp_path):
    """Each isotope carrying the FULL gf is the correct, engine-applies-it form."""
    phys = -1.000
    p = _write_list(tmp_path / 'free.list', 63, 2, {a: phys for a in EU})
    res = audit.audit_species(audit.read_ts_lists([p]), {(63, a): f for a, f in EU.items()},
                              'test')
    assert res[0]['verdict'] == 'FRACTION_FREE'


def test_offset_is_minus_log10_sum_f_squared():
    """The predicted offset is the sum of SQUARED fractions, not log10(2) by luck.

    Eu's isotopes are 0.478/0.522, not 0.5/0.5 — log10(2) = 0.30103 while the
    real answer is 0.30020. They agree to 3 decimals only because Eu happens to
    split nearly evenly; Ba's seven isotopes do not, and there the two differ by
    0.03 dex.
    """
    assert igc.colocated_offset(EU.values()) == pytest.approx(0.30020, abs=1e-5)
    assert igc.colocated_offset(EU.values()) != pytest.approx(math.log10(2), abs=1e-4)
    assert igc.colocated_offset(BA.values()) == pytest.approx(0.26938, abs=1e-5)


def test_mono_isotopic_species_cannot_be_exposed():
    """Mn/Co/Sc are HFS-split but single-isotope: sum f^2 == 1 => zero offset.

    This is the ticket's central premise correction — HFS is not isotope
    structure, and no amount of hyperfine splitting creates this defect.
    """
    for mono in ([1.0], [1.00], [1.0]):
        assert igc.colocated_offset(mono) == pytest.approx(0.0, abs=1e-12)


# ── B. the committed record ──────────────────────────────────────────────────

def test_committed_record_pins_the_exposed_species():
    """The five VALD-surface species RYA-684 measured as folded, and only those."""
    rec = igc.load_audit_record(RECORD)
    exposed = igc.folded_species(igc.VALD_FOR_GRID_SURFACE, rec)
    atomic = {(Z, ion) for (Z, ion), r in exposed.items() if r['kind'] == 'atomic'}
    assert atomic == {(3, 1), (20, 2), (29, 1), (56, 2), (63, 2)}, (
        "the set of fraction-folded species on the shipped VALD for-grid lists "
        "changed — re-run the audit and adjudicate before updating this pin")


def test_committed_record_pins_the_eu_offset():
    """Eu II is the anchor: +0.3002 predicted against RYA-565's measured +0.300."""
    rec = igc.load_audit_record(RECORD)
    assert igc.double_application_offset(63, 2, record=rec) == pytest.approx(0.3002, abs=5e-4)


def test_ges_and_engine_b_surfaces_are_fraction_free():
    """The two surfaces that feed live values must stay form (A).

    Every reported abundance routes through either the iSpec/GES path or the
    Gerber Engine-B deck.  If either ever ships folded gf, live values move.
    """
    rec = igc.load_audit_record(RECORD)
    for surface in ('ges(iSpec path)', 'ts-nlte-ges(Engine-B)'):
        folded = igc.folded_species(surface, rec)
        assert not folded, (
            f"{surface} now ships fraction-folded gf for {sorted(folded)} — every "
            f"abundance measured on that surface is biased high and must be re-derived")


def test_record_still_names_the_engine_application_site():
    """The finding rests on bsyn.f multiplying by isotopfrac; keep the citation live."""
    rec = igc.load_audit_record(RECORD)
    codes = [s['code'] for s in rec['engine_application_site']['sites']]
    assert any('ntot' in c and 'isotopfrac' in c for c in codes)


# ── C. the preflight guard ───────────────────────────────────────────────────

def test_guard_rejects_a_double_folded_target(tmp_path):
    phys = -1.000
    p = _write_list(tmp_path / 'target.list', 63, 2,
                    {a: phys + math.log10(f) for a, f in EU.items()})
    with pytest.raises(SystemExit) as e:
        igc.assert_target_convention(p, 63, 2)
    assert 'ISOTOPE CONVENTION VIOLATION' in str(e.value)


def test_guard_accepts_the_uncoded_target_pattern(tmp_path):
    """The RYA-581 Ba pattern: write Z.000, fold the fraction, engine applies 1.0."""
    p = _write_list(tmp_path / 'ok.list', 56, 2, {0: -1.000})
    out = igc.assert_target_convention(p, 56, 2)
    assert out['target_isotope_coded'] is False
    assert out['blend_blocks_exposed'] == []


def test_guard_records_but_tolerates_exposed_blend_blocks(tmp_path):
    """A folded block belonging to a BLEND species is recorded, never fatal."""
    p = tmp_path / 'blend.list'
    with open(p, 'w') as fh:
        fh.write("'  56.000            '    2         1\n'Ba II'\n")
        fh.write("  5853.6680  0.6040 -1.000   -7.578    2.0  1.00E+08 'x' 'x'   0.0    1.0 'x'\n")
        fh.write("'  63.151            '    2         1\n'Eu II'\n")
        fh.write("  5853.9000  1.2300 -1.320   -7.500    4.0  1.00E+08 'x' 'x'   0.0    1.0 'x'\n")
    out = igc.assert_target_convention(p, 56, 2)
    assert out['blend_blocks_exposed'] == [dict(Z=63, ion=2, isotope=151, n_lines=1)]


def test_guard_loud_fails_without_the_audit_record(tmp_path):
    """No silent fallback to a hardcoded species list."""
    with pytest.raises(SystemExit) as e:
        igc.load_audit_record(tmp_path / 'nope.json')
    assert 'audit record missing' in str(e.value)


# ── the routing conclusion this ticket reported ──────────────────────────────

def test_no_live_value_routes_through_an_exposed_target_block():
    """Every harness window that hits its OWN exposed target is Eu or Ba.

    Eu II is owed-no-value and its VALD leg was a deliberate control; Ba II is
    fitted from a Z.000 block the harness writes itself. If a THIRD target ever
    shows up here, a reported abundance is being fitted against double-folded gf.
    """
    rec = igc.load_audit_record(RECORD)
    own = {w['target'] for w in rec['harness_windows_hit']
           for c in w['contaminants'] if c['is_the_target_species']}
    assert own <= {'Eu II', 'Ba II'}, (
        f"a new target species is being fitted against its own exposed block: "
        f"{sorted(own - {'Eu II', 'Ba II'})}")
