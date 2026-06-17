"""
pipeline/species.py
===================
Canonical species/ion encoding normalizer (RYA-345).

Species/ion identity is encoded inconsistently across the project's data files,
and every cross-file matcher used to hand-handle a different format — getting it
wrong SILENTLY broke matching (a 0-match looks identical to "no such line"). This
cost real debugging time in RYA-339 (MPIA grid 'I'/'II' strings) and RYA-343 (GES
regions 'Fe 2' string + int32 ion + 'T'/'F' molecule flag).

This module collapses every encoding to ONE canonical key:

    atoms     → (Z:int, ion:int)         e.g. Fe II → (26, 2)
    molecules → ('mol', name:str)        e.g. CH    → ('mol', 'CH')

Encodings normalized (the full project inventory — see RYA-345):

    encoding            example            source
    ----------------------------------------------------------------------
    'Symbol ion_int'    'Fe 2'             GES synth line list (element col)
    Symbol + roman      'Fe', 'II'         linelist_solar.csv / solar_ew
    Symbol + int        'Fe', 2            GES regions (ion int32 {1,2,3})
    roman ion           'II'               MPIA NLTE grid (ion col)
    int ion             2                  (generic)
    MOOG/TS code        '26.1' / '26.01'   spectrum_moog_species / turbospectrum
    Z + ion             26, 2              (generic numeric)
    iSpec region note   'Fe 1'             iSpec line_regions 'note' field
    molecule flag       molecule='T'       GES synth list (molecule col)

The element symbol ↔ atomic-number mapping is read from iSpec's
chemical_elements_symbols.dat (the project's existing source of truth) — no
hardcoded periodic table here.

NO silent fallbacks: a species whose ion cannot be resolved RAISES; an
inconsistent (embedded vs explicit) ion RAISES.
"""
from __future__ import annotations

import re
from pathlib import Path

from config.constants import ISPEC_DIR

# ── element symbol ↔ atomic number (from source, cached) ──────────────────────
_CHEM_FILE = ISPEC_DIR / 'input' / 'abundances' / 'chemical_elements_symbols.dat'
_sym2z: dict[str, int] = {}
_z2sym: dict[int, str] = {}

# Roman numeral ↔ int for ionization stages (I=neutral=1 … VIII).
_ROMAN2INT = {'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5, 'VI': 6, 'VII': 7, 'VIII': 8}
_INT2ROMAN = {v: k for k, v in _ROMAN2INT.items()}

# A bare numeric species code: "26", "26.1", "26.01", "26.0".
_NUM_CODE_RE = re.compile(r'^\d+(\.\d*)?$')

MOLECULE = 'mol'  # sentinel head for molecular keys → ('mol', name)


def _load_elements() -> None:
    """Populate symbol↔Z maps from chemical_elements_symbols.dat (TSV)."""
    if _sym2z:
        return
    if not _CHEM_FILE.exists():
        raise FileNotFoundError(
            f"Element table not found: {_CHEM_FILE}\n"
            "Set ISPEC_DIR or run: cd ispec && tar -xzf input.tar.gz"
        )
    with open(_CHEM_FILE) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        zi, si = header.index('atomic_num'), header.index('symbol')
        for line in fh:
            cols = line.rstrip('\n').split('\t')
            if len(cols) <= max(zi, si):
                continue
            z, sym = int(cols[zi]), cols[si].strip()
            _sym2z[sym] = z
            _z2sym[z] = sym


def element_z(symbol: str) -> int:
    """Atomic number for an element symbol. Raises on unknown symbol."""
    _load_elements()
    sym = str(symbol).strip()
    if sym not in _sym2z:
        raise ValueError(f"Unknown element symbol: {symbol!r}")
    return _sym2z[sym]


def z_symbol(z: int) -> str:
    """Element symbol for an atomic number. Raises on unknown Z."""
    _load_elements()
    if int(z) not in _z2sym:
        raise ValueError(f"Unknown atomic number: {z!r}")
    return _z2sym[int(z)]


def parse_ion(ion) -> int:
    """Normalize an ionization stage (int, '2', 'II') → int ≥ 1. Raises on junk."""
    if isinstance(ion, bool):
        raise ValueError(f"Invalid ion: {ion!r}")
    if isinstance(ion, (int,)) or (isinstance(ion, float) and float(ion).is_integer()):
        v = int(ion)
    else:
        s = str(ion).strip().upper()
        if s in _ROMAN2INT:
            v = _ROMAN2INT[s]
        elif s.isdigit():
            v = int(s)
        else:
            raise ValueError(f"Unparseable ion: {ion!r}")
    if v < 1:
        raise ValueError(f"Ion stage must be ≥ 1 (neutral = 1), got {ion!r}")
    return v


def _is_molecule_flag(flag) -> bool:
    if flag is None:
        return False
    if isinstance(flag, bool):
        return flag
    return str(flag).strip().upper() in ('T', 'TRUE', 'MOL', 'Y', 'YES', '1')


def _molecule_key(name) -> tuple:
    """Canonical molecular key — ('mol', cleaned-name). First whitespace token,
    so 'MgH 1' and 'MgH' collapse together."""
    sym = str(name).strip().split()[0]
    return (MOLECULE, sym)


def _ion_from_code(frac: str) -> int:
    """Ionization from a MOOG/Turbospectrum code's fractional part.
    The digits after the decimal are (ion − 1): '26.0'→I, '26.1'→II, '26.01'→II."""
    return (int(frac) if frac != '' else 0) + 1


def _from_numeric_code(code, ion) -> tuple:
    """Resolve a numeric species code ('26.1', '26', 26, 26.1) → key."""
    s = str(code).strip()
    if '.' in s:
        zpart, frac = s.split('.', 1)
    else:
        zpart, frac = s, None
    z = int(zpart)
    if z > 99:                      # molecular MOOG/TS codes (e.g. 106 CH, 822 TiO)
        return (MOLECULE, s)
    code_ion = _ion_from_code(frac) if frac is not None else None
    if ion is not None:
        ion_int = parse_ion(ion)
        if code_ion is not None and code_ion != ion_int:
            raise ValueError(
                f"Inconsistent ion: code {code!r} implies {code_ion}, "
                f"explicit ion is {ion_int}")
        return (z, ion_int)
    if code_ion is None:
        raise ValueError(f"Numeric code {code!r} has no ionization; pass ion=")
    return (z, code_ion)


def species_key(element, ion=None, molecule=None) -> tuple:
    """Map any project species/ion encoding to the canonical key.

    Atoms → (Z, ion);  molecules → ('mol', name).

    Accepts: 'Fe 2', ('Fe', 'II'), ('Fe', 2), ('Fe 2', 2), '26.1', '26.01',
    (26, 2), iSpec note 'Fe 1', molecule flag 'T'/'F'.

    Raises (never silently mis-keys) on unknown symbol or unresolvable/
    inconsistent ion.
    """
    if _is_molecule_flag(molecule):
        return _molecule_key(element)

    # numeric Z or MOOG/TS code
    if isinstance(element, (int, float)) and not isinstance(element, bool):
        return _from_numeric_code(element, ion)

    s = str(element).strip()
    if s == '':
        raise ValueError("Empty species specifier")
    if _NUM_CODE_RE.match(s):
        return _from_numeric_code(s, ion)

    # 'Symbol [embedded-ion]'  e.g. 'Fe', 'Fe 2', 'Fe II', 'Fe 1'
    parts = s.split()
    sym = parts[0]
    embedded = parts[1] if len(parts) > 1 else None

    if sym not in _peek_symbols():
        return _molecule_key(s)     # not an atomic symbol → molecule (CH, MgH, …)

    z = element_z(sym)
    explicit_ion = parse_ion(ion) if ion is not None else None
    embedded_ion = parse_ion(embedded) if embedded is not None else None
    if explicit_ion is not None and embedded_ion is not None and explicit_ion != embedded_ion:
        raise ValueError(
            f"Inconsistent ion for {element!r}: embedded {embedded_ion} "
            f"vs explicit {explicit_ion}")
    ion_int = explicit_ion if explicit_ion is not None else embedded_ion
    if ion_int is None:
        raise ValueError(f"No ionization stage for {element!r}; pass ion=")
    return (z, ion_int)


def _peek_symbols() -> set:
    _load_elements()
    return _sym2z.keys()


def species_note(element, ion=None) -> str:
    """Render the canonical iSpec line-region note 'Symbol ion_int' (e.g. 'Fe 2').
    Replaces the old `_ours_to_ispec_note` and fixes its ion>II bug."""
    key = species_key(element, ion)
    if key[0] == MOLECULE:
        return key[1]
    return f"{z_symbol(key[0])} {key[1]}"
