#!/usr/bin/env python3
"""
scripts/intake_solar_atlases_rya459.py
======================================
RYA-459 (under the RYA-162 epic) — intake the Tier-1 solar reference atlases that
reach the lines HARPS-VIS (380-690 nm) cannot, and wire them with HONEST provenance.

Two Tier-1 sources, two provenance classes:

  1. Kitt Peak Solar Flux Atlas  (Kurucz, Furenlid, Brault & Testerman 1984,
     NSO Atlas No. 1, "Solar Flux Atlas from 296 to 1300 nm") — provenance=MEASURED.
     FTS, McMath/Pierce, NSO/Kitt Peak. 251 segments lm0296..lm1296 (4 nm each).
     Cols: air wavelength (nm), pseudo-residual (normalized) flux, observed
     irradiance (uW/cm^2/nm). THE measured anchor for solar N (N I red + NH 3360),
     CN violet, [O I] 6300, O I 777, and the K/Co/Sc/P diagnostics.

  2. CALSPEC solar reference composite (Colina, Bohlin & Castelli 1996, AJ 112, 307;
     Bohlin, Dickinson & Calzetti 2001, AJ 122, 2118 / sun_reference_stis_002.fits)
     — provenance=CITED-COMPOSITE. Hubble CANNOT observe the Sun directly, so the
     UV is a composite of Woods+1996 (UV) / Neckel&Labs 1984 (vis) / Arvesen+1969 /
     Castelli model (IR). Vacuum Angstrom, FLAM (erg/s/cm^2/A), dlambda ~20 A (low
     resolution). The cited UV reference + an absolute-flux cross-check; NEVER a
     direct solar measurement. Every value carries provenance=cited-composite.

The raw atlases live OUTSIDE the repo (external store, large). This script reads
them, extracts the small diagnostic-region segments + the composite into the repo,
and records full provenance (URL, retrieval date, bytes, sha256) per source.

Usage:  python scripts/intake_solar_atlases_rya459.py
Out:    data/solar_reference/kpno_flux_atlas/{*.csv, kpno_provenance_rya459.json, README.md}
        data/solar_reference/uv_composite/{sun_calspec_composite.csv, uv_provenance_rya459.json}
"""
import hashlib
import json
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))          # run from anywhere; import the pipeline SSOT
# RYA-501: air↔vac come from the single shared converter (Birch & Downs 1994) — no
# local formula copy. NOTE: the canonical vac_to_air keeps λ<2000 Å in VACUUM (IAU/VALD
# convention, air undefined there; RYA-303/426), which the old local copy did not. The
# KPNO path (air→vac, ≥3350 Å) is unaffected; only the UV-composite vac→air dips below
# 2000 Å (see extract_uv_composite). Committed atlas outputs are NOT regenerated here.
from pipeline.wavelength_util import vac_to_air, air_to_vac   # noqa: E402

STORE = (ROOT.parent / 'data' / 'spectra' / 'exoplanetcodex-data' / 'Solar Calibration')
KP_RAW = STORE / 'Kitt Peak Flux Atlas'
UV_RAW = STORE / 'UV Composite'
OUT_KP = ROOT / 'data' / 'solar_reference' / 'kpno_flux_atlas'
OUT_UV = ROOT / 'data' / 'solar_reference' / 'uv_composite'

KP_URL = 'https://nispdata.nso.edu/ftp/pub/atlas/fluxatl/'
UV_URL = ('https://archive.stsci.edu/hlsps/reference-atlases/cdbs/current_calspec/'
          'sun_reference_stis_002.fits')

# Diagnostic regions (air Angstrom) — the lines that matter for the RYA-369 N
# strategy + the CNO arms + the P/K/Co/Sc DATA-GAP elements. window = +/- margin.
DIAGNOSTICS = [
    # name,                element, line_A,  lo_A,   hi_A,   consumer
    ('NH_3360',            'N', 3360.0, 3350.0, 3375.0, 'RYA-369 N (NH UV band head)'),
    ('CN_violet_3883',     'N', 3883.0, 3870.0, 3895.0, 'RYA-369 N cross-check (CN B-X)'),
    ('CoI_3845',           'Co', 3845.46, 3843.0, 3848.0, 'Co DATA-GAP probe'),
    ('ScII_4246',          'Sc', 4246.82, 4244.0, 4249.0, 'Sc DATA-GAP probe'),
    ('OI_6300',            'O', 6300.30, 6295.0, 6305.0, 'O forbidden (cross-check, RYA-455)'),
    ('NI_7442_7468',       'N', 7455.0, 7438.0, 7472.0, 'RYA-369 N I red (7442/7468)'),
    ('KI_7665_7699',       'K', 7682.0, 7662.0, 7702.0, 'K DATA-GAP (resonance doublet)'),
    ('OI_777_triplet',     'O', 7773.0, 7768.0, 7778.0, 'O I 777 (PRIMARY O, RYA-455)'),
    ('NI_8216_8223',       'N', 8219.5, 8210.0, 8228.0, 'RYA-369 N I red (8216/8223)'),
    ('NI_8680_8718',       'N', 8699.0, 8676.0, 8722.0, 'RYA-369 N I red (8680-8718 multiplet)'),
    ('PI_10581_10596',     'P', 10589.0, 10575.0, 10602.0, 'P near-IR multiplet (alt to FUV)'),
]

KP_COLS = ['wavelength_air_A', 'residual_flux', 'irradiance_uW_cm2_nm']


def _digest(path: Path, algo: str) -> str:
    h = hashlib.new(algo)
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def sha256(path: Path) -> str:
    return _digest(path, 'sha256')


def md5(path: Path) -> str:
    return _digest(path, 'md5')


def _load_kp_segment(seg_path: Path) -> pd.DataFrame:
    """One lm#### file: 3 whitespace cols (nm air, residual, irradiance)."""
    arr = np.loadtxt(seg_path)
    df = pd.DataFrame({
        'wavelength_air_A': arr[:, 0] * 10.0,   # nm -> Angstrom
        'residual_flux':    arr[:, 1],
        'irradiance_uW_cm2_nm': arr[:, 2],
    })
    return df


def _kp_inventory():
    return sorted(p for p in KP_RAW.glob('lm[0-9]*') if p.is_file())


def extract_kp_diagnostics():
    OUT_KP.mkdir(parents=True, exist_ok=True)
    segs = _kp_inventory()
    if not segs:
        raise FileNotFoundError(
            f"No Kitt Peak segments in {KP_RAW} — fetch first from {KP_URL} "
            f"(lm0296..lm1296). See README_SOURCE for steps.")
    # segment start nm from filename lm0628 -> 628.0 nm
    seg_start = {p: int(p.name[2:]) for p in segs}

    extracted = []
    for name, el, line_A, lo_A, hi_A, consumer in DIAGNOSTICS:
        lo_nm, hi_nm = lo_A / 10.0, hi_A / 10.0
        # any segment whose [start, start+4.1] overlaps the window
        relevant = [p for p in segs if not (seg_start[p] + 4.11 < lo_nm or seg_start[p] > hi_nm)]
        if not relevant:
            extracted.append({'name': name, 'covered': False, 'line_A': line_A,
                              'window_A': [lo_A, hi_A], 'consumer': consumer})
            continue
        parts = [_load_kp_segment(p) for p in relevant]
        df = pd.concat(parts, ignore_index=True).sort_values('wavelength_air_A')
        df = df.drop_duplicates('wavelength_air_A')
        df = df[(df['wavelength_air_A'] >= lo_A) & (df['wavelength_air_A'] <= hi_A)].reset_index(drop=True)
        df['wavelength_vac_A'] = np.round(air_to_vac(df['wavelength_air_A'].values), 4)
        df = df[['wavelength_air_A', 'wavelength_vac_A', 'residual_flux', 'irradiance_uW_cm2_nm']]
        out = OUT_KP / f'kpno_{name}.csv'
        df.to_csv(out, index=False, float_format='%.5f')
        # median sampling -> nominal resolution proxy
        dl = float(np.median(np.diff(df['wavelength_air_A'].values))) if len(df) > 1 else float('nan')
        extracted.append({
            'name': name, 'element': el, 'line_A': line_A, 'window_A': [lo_A, hi_A],
            'consumer': consumer, 'covered': bool(len(df) > 0),
            'segment_file': out.name, 'segment_bytes': out.stat().st_size,
            'segment_md5': md5(out), 'n_points': int(len(df)),
            'median_dlambda_A': round(dl, 5),
            'wl_air_A_range': [round(float(df['wavelength_air_A'].min()), 3),
                               round(float(df['wavelength_air_A'].max()), 3)] if len(df) else None,
            'flux_min': round(float(df['residual_flux'].min()), 4) if len(df) else None,
        })
    return segs, seg_start, extracted


def write_kp_provenance(segs, seg_start, extracted):
    inv = []
    for p in segs:
        inv.append({'file': p.name, 'start_nm': seg_start[p], 'bytes': p.stat().st_size})
    # manifest md5 over the sorted (name,bytes) tuples — a single fingerprint for the
    # 251-segment raw set without 251 per-file hashes.
    manifest = '\n'.join(f"{s['file']}:{s['bytes']}" for s in inv)
    manifest_md5 = hashlib.md5(manifest.encode()).hexdigest()
    full_lo = min(seg_start.values()) if seg_start else None
    full_hi = (max(seg_start.values()) + 4.1) if seg_start else None
    prov = {
        'ticket': 'RYA-459 (under RYA-162)',
        'source_name': 'Kitt Peak Solar Flux Atlas',
        'citation': ('Kurucz, Furenlid, Brault & Testerman 1984, "Solar Flux Atlas from '
                     '296 to 1300 nm", National Solar Observatory Atlas No. 1. '
                     'Acknowledgement: NSO/Kitt Peak FTS data used here were produced by NSF/NOAO.'),
        'source_url': KP_URL,
        'retrieval_date': date.today().isoformat(),
        'provenance': 'measured',
        'instrument': 'FTS at McMath/Pierce Solar Telescope, NSO/Kitt Peak',
        'native_units': ('col1 = air wavelength (nm, uneven grid); col2 = pseudo-residual '
                         '(continuum-normalized) flux; col3 = observed irradiance (uW/cm^2/nm)'),
        'normalization_state': 'residual_flux = normalized (0-1); irradiance = absolute flux',
        'full_coverage_nm': [full_lo, full_hi],
        'n_segments': len(segs),
        'segment_bytes_total': int(sum(s['bytes'] for s in inv)),
        'raw_manifest_md5': manifest_md5,
        'raw_store': str(KP_RAW),
        'raw_inventory': inv,
        'extracted_diagnostics': extracted,
    }
    (OUT_KP / 'kpno_provenance_rya459.json').write_text(json.dumps(prov, indent=2))
    return prov


def extract_uv_composite():
    OUT_UV.mkdir(parents=True, exist_ok=True)
    fpath = UV_RAW / 'sun_reference_stis_002.fits'
    if not fpath.exists():
        raise FileNotFoundError(
            f"CALSPEC composite missing at {fpath} — fetch from {UV_URL}")
    with fits.open(fpath) as h:
        d = h[1].data
        hdr0 = h[0].header
    w_vac = np.asarray(d['WAVELENGTH'], dtype=float)
    flux = np.asarray(d['FLUX'], dtype=float)
    syserr = np.asarray(d['SYSERROR'], dtype=float)
    fwhm = np.asarray(d['FWHM'], dtype=float)
    df = pd.DataFrame({
        'wavelength_vac_A': np.round(w_vac, 4),
        # vac→air via the shared SSOT: λ<2000 Å stay vacuum (IAU/VALD; the old local
        # copy converted them). Committed output is NOT regenerated by this dedup.
        'wavelength_air_A': np.round(vac_to_air(w_vac), 4),
        'flux_erg_s_cm2_A': flux,
        'syserror': syserr,
        'fwhm_A': fwhm,
        'provenance': 'cited-composite',          # propagates to any downstream value
    })
    out = OUT_UV / 'sun_calspec_composite.csv'
    df.to_csv(out, index=False)
    # the HISTORY block documents which sub-source covers which range — record it
    sub_sources = [
        {'range_A': [1195, 4100], 'source': 'Woods et al. 1996 (UV)', 'kind': 'cited-UV-composite'},
        {'range_A': [4100, 8700], 'source': 'Neckel & Labs 1984', 'kind': 'measured-ground'},
        {'range_A': [8700, 9600], 'source': 'Arvesen et al. 1969', 'kind': 'measured-ground'},
        {'range_A': [9600, 26950], 'source': 'Castelli model', 'kind': 'model'},
    ]
    prov = {
        'ticket': 'RYA-459 (under RYA-162)',
        'source_name': 'CALSPEC solar reference composite (sun_reference_stis_002)',
        'citation': ('Colina, Bohlin & Castelli 1996, AJ 112, 307; '
                     'Bohlin, Dickinson & Calzetti 2001, AJ 122, 2118'),
        'source_url': UV_URL,
        'retrieval_date': date.today().isoformat(),
        'provenance': 'cited-composite',
        'critical_note': ('Hubble cannot observe the Sun directly; this is a COMPOSITE of '
                          'literature UV + ground-based visible + a model IR tail. The UV '
                          '(<4100 A) is a CITED composite, never a direct solar measurement. '
                          'Tagged cited-composite in full so any downstream UV abundance '
                          'inherits the cited (not measured) flag, per the RYA-455 discipline.'),
        'wavelength_frame': 'vacuum (native); air column added (Birch & Downs 1994)',
        'flux_units': 'FLAM = erg/s/cm^2/A (absolute)',
        'resolution': 'low (median dlambda ~20 A; R ~ 150-300) — a flux composite, not a line atlas',
        'n_points': int(len(df)),
        'wl_vac_A_range': [round(float(w_vac.min()), 1), round(float(w_vac.max()), 1)],
        'raw_file': fpath.name,
        'raw_bytes': fpath.stat().st_size,
        'raw_md5': md5(fpath),
        'raw_sha256': sha256(fpath),
        'composite_sub_sources': sub_sources,
        'segment_file': out.name,
        'segment_bytes': out.stat().st_size,
        'segment_md5': md5(out),
    }
    (OUT_UV / 'uv_provenance_rya459.json').write_text(json.dumps(prov, indent=2))
    return prov


def main():
    print(f"RYA-459 solar-atlas intake\n  KP raw : {KP_RAW}\n  UV raw : {UV_RAW}\n")
    segs, seg_start, extracted = extract_kp_diagnostics()
    kp = write_kp_provenance(segs, seg_start, extracted)
    print(f"Kitt Peak: {len(segs)} raw segments ({kp['full_coverage_nm']} nm); "
          f"{sum(e['covered'] for e in extracted)}/{len(extracted)} diagnostics extracted -> {OUT_KP}")
    uv = extract_uv_composite()
    print(f"CALSPEC UV composite: {uv['n_points']} pts {uv['wl_vac_A_range']} A "
          f"[{uv['provenance']}] -> {OUT_UV}")
    print("\nDone. Run: python -m pipeline.audit_solar_reference --verify")


if __name__ == '__main__':
    main()
