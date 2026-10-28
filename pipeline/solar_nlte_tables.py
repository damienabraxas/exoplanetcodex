"""Declared, solar-only same-EW corrections from published model columns.

No extrapolation to other stars. A cited 3D-NLTE minus MARCS-LTE correction
added to our LTE fit is labelled as such, never as in-transfer 3D synthesis.
Canonical IDs resolve rounded published wavelengths without fuzzy feature joins.
"""
from pathlib import Path
import csv
import math
import yaml

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / 'data/reference/solar_nlte_tables.yaml'


def declarations(element, ion, star):
    rows = yaml.safe_load(REGISTRY.read_text())['tables']
    selected = [r for r in rows if r['element'] == element and str(ion) in r['ions']
                and r['star'] == star]
    # RYA-1251: an element may have more than one cited table (Na: AGSS21 Table A.1 for
    # its five lines, Canocchi+2026 3D for the others) -- but never two for one line;
    # `deltas` refuses that as ambiguous.
    return selected


def input_paths(element, ion, star):
    specs = declarations(element, ion, star)
    return ([REGISTRY, Path(__file__)] + [ROOT / r['source'] for r in specs]) if specs else []


def deltas(element, ion, waves, leg, *, star):
    out = {}
    specs = declarations(element, ion, star)
    if not specs:
        return out
    canonical = {r['line_id']: r for r in csv.DictReader(
        (ROOT / 'data/linelists/canonical_gf.csv').open())}
    for spec in specs:
        if leg not in spec['legs']:
            continue            # a table may cite one leg only (Canocchi+2026: 3D)
        for row in csv.DictReader((ROOT / spec['source']).open()):
            species = f'{element} {ion}'
            if row[spec['species_column']] != species:
                continue
            key = row.get(spec.get('canonical_id_column', ''), '').strip()
            if key:
                c = canonical[key]
                if c['species'] != species:
                    raise ValueError(f'Wrong-ion canonical join: {key}')
                w = float(c['wavelength_air_A'])
            else:
                w = float(row[spec['wavelength_column']])
            delta = float(row[spec['legs'][leg]]) - float(row[spec['baseline']])
            if not math.isfinite(delta):
                raise ValueError(f'Non-finite cited correction for {species} {w}')
            for requested in waves:
                requested = float(requested)
                if abs(requested - w) <= spec['match_tol_A']:
                    if requested in out:
                        raise ValueError(f'Ambiguous cited correction: {species} {requested}')
                    out[requested] = round(delta, 4)
    return out


def source(element, ion, leg, *, star):
    return '; '.join(f"{s['citation']}; CITED solar {leg}-NLTE minus 1D-LTE per-line "
                     "additive correction; our atmosphere may differ; not in-transfer 3D"
                     for s in declarations(element, ion, star) if leg in s['legs'])


def applied_source(element, ion, applied, *, star, native_source):
    """Name only sources that actually supplied a finite correction in this pool."""
    waves = [float(w) for w, delta in applied.items() if math.isfinite(float(delta))]
    cited = deltas(element, ion, waves, "1D", star=star)
    parts = []
    if any(w not in cited for w in waves):
        parts.append(native_source)
    if cited:
        parts.append(source(element, ion, "1D", star=star))
    return '; '.join(parts) or 'no finite NLTE correction supplied'


def model(element, ion, *, star):
    models = {s['model'] for s in declarations(element, ion, star)}
    if len(models) != 1:
        raise ValueError(f'Ambiguous cited model axis: {models}')
    return models.pop()
