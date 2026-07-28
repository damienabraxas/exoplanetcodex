"""
tests/test_vac_air_single_source_rya501.py
==========================================
RYA-501 — single-source guard for the vacuum↔air converter (same guard class as the
RYA-355 single-source gf check). The Birch & Downs 1994 refractive-index formula must
be DEFINED in exactly one place — ``pipeline/wavelength_util.py`` — and every other
user imports it. A second copy is a defect the moment it exists (it can drift on the
next edit), so this fails loudly if the formula's constants reappear anywhere else.

This closes the RYA-264 straggler: two historical intake scripts
(``intake_ir_atlases_rya390``, ``intake_solar_atlases_rya459``) once carried a local
B&D copy; they now import the shared converter.
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CANONICAL = ROOT / 'pipeline' / 'wavelength_util.py'
GUARD = Path(__file__).resolve()

# Distinctive refractive-index coefficients — matches BOTH the Birch & Downs 1994
# form (2.406147e-2 / 0.02406147; 8.34254e-5 / 0.0000834254) AND the Morton 2000
# variant (0.02408926 / 0.00008336) so a reintroduction of either is caught.
_FORMULA_TOKENS = re.compile(r'2\.406147|02406147|02408926|8\.34254|0000834254|00008336')


def _repo_py_files():
    for p in ROOT.rglob('*.py'):
        if '.git' in p.parts:
            continue
        yield p


def test_formula_defined_in_exactly_one_place():
    hits = []
    for p in _repo_py_files():
        if p == GUARD:                      # the guard names the constants to search for
            continue
        text = p.read_text(encoding='utf-8', errors='ignore')
        if _FORMULA_TOKENS.search(text):
            hits.append(p.relative_to(ROOT).as_posix())
    assert hits == ['pipeline/wavelength_util.py'], (
        "vac↔air refractive-index formula must be defined ONLY in "
        "pipeline/wavelength_util.py (single source of truth, RYA-264/501). "
        f"Found the constants in: {hits}. Import vac_to_air/air_to_vac from the "
        "shared util instead of re-deriving Birch & Downs / Morton locally.")


def test_canonical_actually_defines_it():
    # guard against the guard: the one allowed locus really does hold the formula
    assert _FORMULA_TOKENS.search(CANONICAL.read_text())


def test_intake_scripts_import_the_shared_converter():
    for name in ('intake_ir_atlases_rya390.py', 'intake_solar_atlases_rya459.py'):
        src = (ROOT / 'scripts' / name).read_text()
        assert 'from pipeline.wavelength_util import' in src, name
        # and no local converter function survives
        assert 'def vac_to_air_A' not in src and 'def air_to_vac_A' not in src, name


def test_pipeline_users_delegate_not_redefine():
    # the other in-pipeline users import/delegate; none redefines the constants
    for mod in ('uv_conditioning.py', 'uv_line_selection.py', 'crires_telluric.py'):
        src = (ROOT / 'pipeline' / mod).read_text()
        assert 'wavelength_util' in src, mod
        assert not _FORMULA_TOKENS.search(src), f"{mod} redefines the formula"
