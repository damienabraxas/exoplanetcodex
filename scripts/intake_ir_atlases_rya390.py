#!/usr/bin/env python3
"""
scripts/intake_ir_atlases_rya390.py
===================================
RYA-390 Part A — intake the three IR reference atlases for the K-band CO arm and
extract the CO segment (4255–4367 cm⁻¹ ≈ 2.290–2.350 µm) from each.

Three references, three roles (RYA-373 three-way validation):
  1. ACE-FTS solar atlas        — SPACE, telluric-FREE solar truth (primary).
                                   Hase, Wallace et al. 2010, JQSRT 111, 521.
  2. NSO Kitt Peak photatl       — TERRESTRIAL solar (residual telluric), the
                                   cross-instrument reference. Livingston &
                                   Wallace 1991, NSO Tech. Report #91-001.
  3. Wallace telluric near-IR    — pure TELLURIC transmission (validates the
                                   molecfit telluric MODEL directly). Wallace,
                                   Hinkle & Livingston, NOAO/NSO.

All three are FTS products in **vacuum wavenumber (cm⁻¹)**. We carry the native
wavenumber and add vacuum + air wavelength (Å) per segment so downstream RYA-373
can match either convention. Air conversion: Birch & Downs (1994) / Edlén — the
formula VALD3/iSpec use.

Raw atlases live OUTSIDE the repo (large), in the external spectra tree; this
script reads them, writes the small CO segments + provenance into the repo.

Usage:  python scripts/intake_ir_atlases_rya390.py
Out:    data/solar_reference/ir_atlases/{*_co_*.csv, ir_atlases_provenance_rya390.json}
"""
import hashlib
import json
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = (ROOT.parent / 'data' / 'spectra' / 'exoplanetcodex-data' /
       'Solar Calibration' / 'IR Reference Atlases')
OUT = ROOT / 'data' / 'solar_reference' / 'ir_atlases'

# K-band CO segment (CO Δv=2 first-overtone region used by the ¹³C/¹⁸O arm).
CO_LO_CM, CO_HI_CM = 4255.0, 4367.0


def vac_to_air_A(wl_vac_A: np.ndarray) -> np.ndarray:
    """Birch & Downs (1994) / Edlén vacuum→air wavelength (Å), the VALD3/iSpec
    convention. s = 1e4/λ_vac(Å) in µm⁻¹. Valid across the optical–near-IR."""
    s2 = (1.0e4 / wl_vac_A) ** 2
    n = (1.0 + 0.0000834254 + 0.02406147 / (130.0 - s2)
         + 0.00015998 / (38.9 - s2))
    return wl_vac_A / n


def add_wavelengths(df: pd.DataFrame, wn_col: str = 'wavenumber_cm-1') -> pd.DataFrame:
    """Attach vacuum + air wavelength (Å) from the FTS vacuum wavenumber."""
    wl_vac = 1.0e8 / df[wn_col].to_numpy(dtype=float)   # cm⁻¹ → Å (vacuum)
    df.insert(1, 'wavelength_vac_A', np.round(wl_vac, 4))
    df.insert(2, 'wavelength_air_A', np.round(vac_to_air_A(wl_vac), 4))
    return df


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def in_band(wn: np.ndarray) -> np.ndarray:
    return (wn >= CO_LO_CM) & (wn <= CO_HI_CM)


# ── parsers ───────────────────────────────────────────────────────────────────

def parse_ace(path: Path) -> pd.DataFrame:
    """ACE-FTS CSV: 'MRCA,' header then wavenumber_cm-1,intensity (telluric-free)."""
    df = pd.read_csv(path, header=None, names=['wavenumber_cm-1', 'intensity'],
                     skiprows=1)
    df = df[in_band(df['wavenumber_cm-1'].to_numpy(dtype=float))].reset_index(drop=True)
    return add_wavelengths(df)


def parse_photatl(files: list[Path]) -> pd.DataFrame:
    """photatl wnNNNN: whitespace wavenumber, solar, atmospheric, total
    (Livingston & Wallace). Files overlap 2 cm⁻¹ each end — concat, band-clip,
    drop duplicate wavenumbers from the overlap."""
    frames = []
    for p in files:
        frames.append(pd.read_csv(
            p, sep=r'\s+', header=None,
            names=['wavenumber_cm-1', 'solar', 'atmospheric', 'total']))
    df = pd.concat(frames, ignore_index=True)
    df = df[in_band(df['wavenumber_cm-1'].to_numpy(dtype=float))]
    df = (df.drop_duplicates(subset='wavenumber_cm-1')
            .sort_values('wavenumber_cm-1').reset_index(drop=True))
    return add_wavelengths(df)


def parse_wallace_ratio(path: Path) -> pd.DataFrame:
    """Wallace ratio04300.txt: whitespace wavenumber, telluric_ratio (transmission)."""
    df = pd.read_csv(path, sep=r'\s+', header=None,
                     names=['wavenumber_cm-1', 'telluric_ratio'])
    df = df[in_band(df['wavenumber_cm-1'].to_numpy(dtype=float))].reset_index(drop=True)
    return add_wavelengths(df)


def seg_meta(df: pd.DataFrame) -> dict:
    wn = df['wavenumber_cm-1']
    return {'n_points': int(len(df)),
            'wn_cm-1_range': [round(float(wn.min()), 4), round(float(wn.max()), 4)],
            'wl_air_A_range': [round(float(df['wavelength_air_A'].min()), 3),
                               round(float(df['wavelength_air_A'].max()), 3)]}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sources = []

    # 1. ACE-FTS — telluric-free solar truth
    ace_raw = RAW / 'ACE-FTS' / 'ace-solar-spectrum.txt'
    ace = parse_ace(ace_raw)
    ace.to_csv(OUT / 'ace_fts_solar_co_4255_4367.csv', index=False)
    sources.append({
        'name': 'ACE-FTS solar atlas', 'role': 'solar truth (telluric-free, PRIMARY)',
        'citation': 'Hase, Wallace, McLeod, Harrison & Bernath 2010, JQSRT 111, 521 '
                    '("A complete solar spectrum based on ACE data")',
        'source_url': 'https://databace.scisat.ca/solarspectrum/',
        'raw_file': ace_raw.name, 'raw_bytes': ace_raw.stat().st_size,
        'raw_sha256': sha256(ace_raw),
        'native_units': 'vacuum wavenumber (cm⁻¹), 0.005 cm⁻¹ step; normalized intensity',
        'full_range_cm-1': [700.0, 4430.0], 'resolution_cm-1': 0.02, 'snr': '~400',
        'telluric_status': 'telluric-FREE (high-sun ACE occultation, space-based)',
        'disk_integration': 'roughly disk-integrated (occultation) — matches reflected-Vesta',
        'segment_file': 'ace_fts_solar_co_4255_4367.csv',
        'segment_columns': ['wavenumber_cm-1', 'wavelength_vac_A', 'wavelength_air_A', 'intensity'],
        **seg_meta(ace)})

    # 2. NSO photatl — terrestrial solar (solar + atmospheric + total columns)
    ph_files = [RAW / 'NSO_photatl' / f for f in
                ('wn4250', 'wn4275', 'wn4300', 'wn4325', 'wn4350')]
    ph = parse_photatl(ph_files)
    ph.to_csv(OUT / 'nso_photatl_co_4255_4367.csv', index=False)
    sources.append({
        'name': 'NSO Kitt Peak IR photosphere atlas (photatl)',
        'role': 'terrestrial solar cross-reference (independent ground-based reduction)',
        'citation': 'Livingston & Wallace 1991, "An Atlas of the Solar Spectrum in the '
                    'Infrared from 1850 to 9000 cm⁻¹ (1.1 to 5.4 µm)", NSO Tech. Report '
                    '#91-001. Acknowledgement: "NSO/Kitt Peak FTS data used here were '
                    'produced by NSF/NOAO."',
        'source_url': 'https://nispdata.nso.edu/ftp/pub/atlas/photatl/',
        'raw_file': '+'.join(f.name for f in ph_files),
        'raw_bytes': sum(f.stat().st_size for f in ph_files),
        'raw_sha256': {f.name: sha256(f) for f in ph_files},
        'native_units': 'vacuum wavenumber (cm⁻¹); cols = solar, atmospheric, total',
        'telluric_status': 'TERRESTRIAL (disk-center): col solar=telluric-corrected, '
                           'col atmospheric=residual telluric, col total=observed. '
                           'NB gaps in strong-telluric regions are LINEAR INTERPOLATIONS '
                           '(README); intensities not continuous file-to-file.',
        'disk_integration': 'disk-CENTER (not disk-integrated) — caveat vs reflected-Vesta',
        'segment_file': 'nso_photatl_co_4255_4367.csv',
        'segment_columns': ['wavenumber_cm-1', 'wavelength_vac_A', 'wavelength_air_A',
                            'solar', 'atmospheric', 'total'],
        **seg_meta(ph)})

    # 3. Wallace telluric near-IR — pure telluric transmission
    wr_raw = RAW / 'Wallace_telluric' / 'ratio04300.txt'
    wr = parse_wallace_ratio(wr_raw)
    wr.to_csv(OUT / 'wallace_telluric_co_ratio.csv', index=False)
    sources.append({
        'name': 'Wallace Telluric Near-IR Atlas',
        'role': 'pure-telluric reference (validates the molecfit telluric MODEL)',
        'citation': 'Wallace, Hinkle & Livingston, "The Absorption Spectrum of the '
                    "Earth's Atmosphere from 0.578 to 5.43 µm (1800 to 17400 cm⁻¹)\", "
                    'NOAO/NSO. Telluric = ratio of 5-airmass / 2-airmass FTS scans '
                    '(1990/12/18 #4 & #5, center disc) → solar cancels.',
        'source_url': 'https://nispdata.nso.edu/ftp/pub/Wallace_telluric_near_ir_atlas/',
        'raw_file': wr_raw.name, 'raw_bytes': wr_raw.stat().st_size,
        'raw_sha256': sha256(wr_raw),
        'native_units': 'vacuum wavenumber (cm⁻¹); telluric_ratio = transmission (0–1)',
        'telluric_status': 'PURE TELLURIC transmission (Earth atmosphere)',
        'disk_integration': 'N/A (telluric, not solar)',
        'segment_file': 'wallace_telluric_co_ratio.csv',
        'segment_columns': ['wavenumber_cm-1', 'wavelength_vac_A', 'wavelength_air_A',
                            'telluric_ratio'],
        'COVERAGE_CAVEAT': 'the ASCII telluric ratio (ratio04300.txt) spans only '
                           '4299.8–4338.6 cm⁻¹ — the MIDDLE of the 4255–4367 CO band. '
                           'Full-band telluric: use the photatl atmospheric column '
                           '(4248–4377) or the Wallace eps plots (atl04240–atl04360). '
                           'Telluric line IDs: linelist_TOTAL/H2O/CH4_ext.txt.',
        **seg_meta(wr)})

    manifest = {
        'ticket': 'RYA-390 Part A', 'date': str(date.today()),
        'co_band_cm-1': [CO_LO_CM, CO_HI_CM],
        'co_band_air_A': [round(float(vac_to_air_A(np.array([1e8 / CO_HI_CM]))[0]), 1),
                          round(float(vac_to_air_A(np.array([1e8 / CO_LO_CM]))[0]), 1)],
        'convention': 'native = FTS vacuum wavenumber (cm⁻¹); wavelength_vac_A = 1e8/wn; '
                      'wavelength_air_A via Birch & Downs (1994)/Edlén (VALD3/iSpec). '
                      'Vesta (reflected solar) = integrated-disk.',
        'raw_store': str(RAW),
        'sources': sources,
    }
    (OUT / 'ir_atlases_provenance_rya390.json').write_text(json.dumps(manifest, indent=2))

    print(f"\n=== RYA-390 IR atlas intake — CO band {CO_LO_CM}-{CO_HI_CM} cm⁻¹ "
          f"({manifest['co_band_air_A'][0]}-{manifest['co_band_air_A'][1]} Å air) ===")
    for s in sources:
        print(f"\n  {s['name']}  [{s['role']}]")
        print(f"    segment: {s['segment_file']}  ({s['n_points']} pts, "
              f"{s['wn_cm-1_range'][0]}-{s['wn_cm-1_range'][1]} cm⁻¹, "
              f"{s['wl_air_A_range'][0]}-{s['wl_air_A_range'][1]} Å air)")
        print(f"    telluric: {s['telluric_status'][:70]}")
    print(f"\n  [out] {OUT}/  (+ ir_atlases_provenance_rya390.json)")


if __name__ == '__main__':
    main()
