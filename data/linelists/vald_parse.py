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


# =============================================================================
# RYA-321 — metallicity intake gate
# -----------------------------------------------------------------------------
# A metal-rich star extracted at the solar default under-selects weak metal lines
# (an incomplete pool, worst for trace / neutron-capture species), because VALD
# selects lines from the stellar atmosphere at the requested COMPOSITION. The
# composition is set via the "Extract Stellar" form's free-text Chemical
# composition field ('M/H: <value>').
#
# WHERE TO READ THE APPLIED [M/H] (corrected 2026-06-15 from the job.019562 email
# + delivery): NOT the Castelli model filename. VALD splits the request into two
# things:
#   * model-atmosphere STRUCTURE — a Castelli/Kurucz grid node. When VALD lacks
#     the exact metal-rich structure it substitutes the solar-structure node and
#     stamps 'castelli_ap00k2_...' in the footer (with a "VALD does not have the
#     exact model, will use ... instead" line in the job email) EVEN WHEN M/H was
#     applied. So the filename is the WRONG signal — it false-rejects correct
#     metal-rich extractions. (Earlier RYA-321 audit read the filename and wrongly
#     concluded α Cen / 55 Cnc were solar-default; they were not.)
#   * COMPOSITION — the element abundance block that follows the filename, where
#     ALL metals are uniformly offset from VALD's solar scale by the requested
#     [M/H] (H, He fixed). This IS what drove line selection. We read it here.
#
# VALD's solar reference abundances (log N_X/N_tot, the [M/H]=0 block, taken from
# vald_solar_raw.txt). Applied [M/H] = delivered abundance - this reference,
# averaged over these metals (they agree exactly under uniform scaling). This is a
# property of VALD's line-selection solar scale (like TRUNCATION_WARNING), NOT a
# per-star metallicity, so it does not violate the no-hardcoded-metallicity rule
# (that rule binds feh_ref, still read live from STAR_PARAMS).
_VALD_SOLAR_ABUND = {
    'Na': -5.71, 'Mg': -4.46, 'Si': -4.49, 'Ca': -5.68,
    'Ti': -7.02, 'Cr': -6.37, 'Fe': -4.54, 'Ni': -5.79,
}

# Quoted 'El: -x.xx' tuples of the abundance block (e.g. 'Fe: -4.54', 'V : -8.04').
# Transition records ('Fe 1', 3780,...) have no colon after the element, so they
# do not match.
_ABUND_RE = re.compile(r"'([A-Z][a-z]?)\s*:\s*(-?\d+\.\d+)'")


def applied_metallicity_from_lines(lines, max_spread=0.03):
    """Core: compute applied [M/H] from an iterable of text lines (so a caller can
    stream a .gz without decompressing to disk). See parse_vald_applied_metallicity
    for semantics. Returns the mean metal offset, or None."""
    found = {}
    for line in lines:
        for el, val in _ABUND_RE.findall(line):
            if el in _VALD_SOLAR_ABUND and el not in found:
                found[el] = float(val)
        if len(found) == len(_VALD_SOLAR_ABUND):
            break
    if not found:
        return None
    offsets = [found[el] - _VALD_SOLAR_ABUND[el] for el in found]
    if max(offsets) - min(offsets) > max_spread:
        return None
    return sum(offsets) / len(offsets)


def parse_vald_applied_metallicity(path, max_spread=0.03):
    """Return the [M/H] VALD applied to the composition (from the footer abundance
    block), or None if it cannot be read. Computed as the mean offset of the
    reference metals (_VALD_SOLAR_ABUND) from VALD's solar scale. The per-element
    offsets must agree to within `max_spread` (VALD scales all metals uniformly);
    a larger spread means a parse error or a changed solar reference, so we return
    None (fail loud — do not trust a half-read block)."""
    with open(path, errors='replace') as f:
        return applied_metallicity_from_lines(f, max_spread)


def parse_vald_structure_grid(path):
    """Informational: the [M/H] of the Castelli STRUCTURE model VALD used (footer
    'castelli_apNNk2_...' filename), or None. This is the atmosphere structure
    node, NOT the applied composition (see module note) — used only to report when
    VALD substituted a solar-structure model. Do NOT gate on this."""
    m = None
    rx = re.compile(r'castelli_a([pm])(\d{2})k\d', re.IGNORECASE)
    with open(path, errors='replace') as f:
        for line in f:
            m = rx.search(line)
            if m:
                sign = -1.0 if m.group(1).lower() == 'm' else 1.0
                return sign * int(m.group(2)) / 10.0
    return None


def verify_metallicity(path, star_id, star_params, tol=0.05):
    """Intake gate: the [M/H] VALD APPLIED to the composition must match the star's
    catalog [Fe/H]. Source of truth is STAR_PARAMS[star_id]['feh_ref'] (RYA-292) —
    never hardcoded. A missing record (no STAR_PARAMS entry, or no abundance block
    in the delivery) REJECTs loudly: cannot verify == do not trust, no silent
    fallback. Returns (verdict, message).

    tol defaults to 0.05 dex (~1σ on the GBS [Fe/H]); the composition is applied
    continuously so this only needs to absorb reasonable reference differences
    (e.g. a +0.35 vs catalog +0.31 [Fe/H] choice), while still catching a true
    solar-default extraction (|Δ| ≈ the star's [Fe/H], ≥0.2 for our metal-rich set).

    NB: the spec drafted the lookup as STAR_PARAMS[...]['feh'], but the live record
    key is 'feh_ref' (the catalog [Fe/H]; [Fe/H] is solved, not pinned).
    """
    rec = star_params.get(star_id)
    if not rec or 'feh_ref' not in rec:
        return ('REJECT',
                f"{path}: no STAR_PARAMS[{star_id}].feh_ref — cannot verify "
                f"metallicity. Wire the star into STAR_PARAMS (RYA-292) first.")
    catalog = float(rec['feh_ref'])
    found = parse_vald_applied_metallicity(path)
    if found is None:
        return ('REJECT',
                f"{path}: no composition (abundance block) readable — cannot "
                f"confirm applied [M/H] (catalog [Fe/H]={catalog:+.2f}). Re-extract "
                f"with 'M/H: {catalog:+.2f}' in the Chemical composition field.")
    if abs(found - catalog) > tol:
        return ('REJECT',
                f"{path}: METALLICITY MISMATCH — VALD applied [M/H]={found:+.2f}, "
                f"catalog [Fe/H]={catalog:+.2f} (|Δ|={abs(found - catalog):.2f} > "
                f"{tol}). Re-extract with 'M/H: {catalog:+.2f}' in Chemical "
                f"composition.")
    return ('ACCEPT',
            f"{path}: metallicity OK (VALD applied [M/H]={found:+.2f} vs catalog "
            f"{catalog:+.2f}, |Δ|={abs(found - catalog):.2f}).")
