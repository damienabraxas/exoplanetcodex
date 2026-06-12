# =============================================================================
# THE EXOPLANET CODEX
# exoplanetcodex.org  |  github.com/damienabraxas/exoplanetcodex
# =============================================================================
#
# File:         vald_parse.py
# Module:       data/linelists (shared linelist tooling)
# Description:  Shared VALD3 "Extract Stellar" Long Format parsing utilities.
#               One parser for all stars — refactored from the RYA-223 Procyon
#               inspection logic so every unpack/merge script shares the same
#               data-line identification and field extraction.
#
# Author:       Ryan Schmitt
# Contributors: Claude (Anthropic) via Claude Code
# Created:      2026-06-12
# Last modified: 2026-06-12
# Linear issue: RYA-269 — BUILD: 55 Cnc A VALD UV + NIR unpack/verify/merge
#
# -----------------------------------------------------------------------------
# KEY REFERENCES
# -----------------------------------------------------------------------------
# Ryabchikova et al. 2015 — Phys. Scr. 90, 054005 — VALD3 description
#
# -----------------------------------------------------------------------------
# DEPENDENCIES
# -----------------------------------------------------------------------------
# External: none (stdlib only)
# =============================================================================

import re
from pathlib import Path

# VALD web extractions are silently capped at 100,000 output transitions;
# a capped delivery carries this warning on line 1 (RYA-64 post-mortem).
TRUNCATION_WARNING = 'WARNING: Output was truncated to 100000 lines'

_ROMAN = {1: 'I', 2: 'II', 3: 'III', 4: 'IV', 5: 'V'}


def read_vald_header(path):
    """
    Read the VALD metadata line and detect web-cap truncation.

    The metadata line reads:
      ' 5000.00000, 30000.00000, NNNNN, MMMMMMM, 0.9 Wavelength region, ...'
    On truncated deliveries it is preceded by the TRUNCATION_WARNING line.

    Returns dict: wl_start, wl_end, n_selected, n_processed, vmicro, truncated
    """
    with open(path, errors='replace') as f:
        line1 = f.readline().strip()
        truncated = line1.startswith(TRUNCATION_WARNING.split(':')[0])
        meta = f.readline().strip() if truncated else line1
    fields = [p.strip() for p in meta.split(',')]
    return {
        'wl_start'   : float(fields[0]),
        'wl_end'     : float(fields[1]),
        'n_selected' : int(fields[2]),
        'n_processed': int(fields[3]),
        'vmicro'     : float(fields[4].split()[0]),
        'truncated'  : truncated,
    }


# Transition species token: bare symbol + ionisation stage, e.g. 'Fe 1',
# 'MgH 1', 'C2 1'. Distinguishes transitions from the trailing model-atmosphere
# block ('castelli_...krz', 'H :  0.92', ...) whose quoted lines otherwise
# pass the structural test.
_SPECIES_RE = re.compile(r'[A-Za-z][A-Za-z0-9]{0,3} [1-5]')


def is_vald_data_line(line):
    """
    Identify a Long Format transition data line by structure (RYA-223):
      Data:   'Fe 1',  3780.26, ...   → text after the closing quote starts ','
      Config: '  LS   ...'            → nothing after the closing quote
      Ref:    '_  Kurucz...'          → nothing after the closing quote
      Footer: 'castelli_...krz', 'H :  0.92', ...  → species token malformed
    Returns the (species, remainder) pair for data lines, else None.
    """
    if not line.startswith("'"):
        return None
    quote_parts = line.split("'")
    if len(quote_parts) < 3 or not quote_parts[2].startswith(','):
        return None
    species = quote_parts[1].strip()
    if not _SPECIES_RE.fullmatch(species):
        return None
    return species, quote_parts[2]


def parse_vald_long(path, max_examples=5):
    """
    Parse a VALD3 Long Format extraction into transition records.

    Each transition spans 4 physical lines; only the leading data line is
    consumed. Fields per the Long Format column header:
      WL_air(A), log gf, E_low(eV), J lo, E_up(eV), J up,
      Lande lower/upper/mean, Rad., Stark, Waals damping, central depth

    NOTE: VALD labels the column WL_air(A) but delivers vacuum wavelengths
    below 2000 Å. No conversion is performed here — callers must record the
    convention per line (RYA-269 spec).

    Returns (records, failures):
      records  — list of dicts: species, element, ion (Roman), wavelength,
                 log_gf, e_low_eV, damping_rad, damping_stark, damping_vdW,
                 central_depth
      failures — list of (line_number, line_text, error) up to max_examples,
                 plus n_failures total as the last element's count; callers
                 must treat a non-empty list as reportable (no silent drops).
    """
    records = []
    failures = []
    n_failures = 0

    with open(path, errors='replace') as f:
        for i, line in enumerate(f, start=1):
            hit = is_vald_data_line(line)
            if hit is None:
                continue
            species, rest = hit
            try:
                fields = rest.split(',')
                parts = species.split()
                if len(parts) != 2:
                    raise ValueError(f'unexpected species token {species!r}')
                element = parts[0]
                ion = _ROMAN[int(parts[1])]
                records.append({
                    'species'       : species,
                    'element'       : element,
                    'ion'           : ion,
                    'wavelength'    : float(fields[1]),
                    'log_gf'        : float(fields[2]),
                    'e_low_eV'      : float(fields[3]),
                    'damping_rad'   : float(fields[10]),
                    'damping_stark' : float(fields[11]),
                    'damping_vdW'   : float(fields[12]),
                    'central_depth' : float(fields[13]),
                })
            except (ValueError, IndexError, KeyError) as e:
                n_failures += 1
                if len(failures) < max_examples:
                    failures.append((i, line.rstrip()[:120], str(e)[:80]))

    return records, {'n_failures': n_failures, 'examples': failures,
                     'n_parsed': len(records)}
