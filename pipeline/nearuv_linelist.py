"""
pipeline/nearuv_linelist.py — the near-UV VALD holdings, as an iSpec line list.
==============================================================================
RYA-759 Move 2. Brings `vald_solar_nearuv_2000_3780_hfson_raw.txt` into SYNTHESIS
by writing it in iSpec's own atomic-linelist format, so that everything below the
production synth core is unchanged: `ispec.generate_spectrum(code='turbospectrum')`
writes the Turbospectrum file itself, from this array, exactly as it does for the
optical GES list.

WHY NOT `scripts/vald_to_turbospectrum.py`
------------------------------------------
That script wrote the Turbospectrum `.lte` file BY HAND, for raw `babsma_lu`/`bsyn_lu`.
Move 1 established two things about that route:

  * the 0-byte spectrum was never the list — it was `babsma` dying on a model/flag
    mismatch, so the hand-written list has never actually been read by `bsyn`; and
  * the raw-binary path is the Gerber NLTE-anchor validation gate's plumbing. A
    validation harness's plumbing must not become a production abundance path.

So the format problem is handed back to iSpec, which already owns it
(`ispec.lines.__turbospectrum_write_atomic_linelist`). Every trap the hand-rolled
converter hit — species-code alignment, the LTE-vs-NLTE row layout — is a trap in
code we no longer write. The VALD extraction is not wasted: it is the input here.

WHAT COMES FROM WHERE — no second parser, no invented physics
-------------------------------------------------------------
The VALD walk is `data.linelists.vald_parse.parse_vald_long`, the project's shared,
no-silent-drop intake (RYA-223) — the same reader the per-star intakes use, extended
in RYA-759 to carry the upper-level and Landé fields it already read. Element → Z is
iSpec's own chemical-elements table, not a table written here.

Two VALD columns are deliberately NOT translated, and both are stated rather than
guessed (the ticket's standing instruction):

  * **Stark.** VALD's Stark parameter is carried into `stark` for the record, but
    Turbospectrum's LTE line format has no Stark field at all — iSpec's writer never
    emits it. Nothing is assumed about matching conventions.
  * **van der Waals = 0.** 13,029 of the 3000–3780 Å lines carry vdW = 0 because VALD
    supplied none. That is passed through as 0.0, which is Turbospectrum's documented
    "use your own default" value. It is a stated choice; `band_stats()` reports the
    count so the choice is visible in every run rather than buried.

Orbital types are 'X' where VALD gives none. That is not an invention: 'X' is already
present 2,947 times in the production GES list iSpec ships, so it is a value the
downstream reader is known to accept.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from data.linelists.vald_parse import parse_vald_long  # noqa: E402

#: eV → cm⁻¹ (CODATA 2018). iSpec's format carries both; VALD delivers only eV.
EV_TO_CM1 = 8065.543937

#: Default near-UV holdings + band. The band's blue edge is 3000 Å because that is
#: where the Kitt Peak flux atlas starts (2960 Å) and where our extraction has depth;
#: the red edge is the VALD delivery's own 3780 Å boundary.
DEFAULT_RAW = ROOT / 'data' / 'linelists' / 'vald_solar_nearuv_2000_3780_hfson_raw.txt'
DEFAULT_LO_A, DEFAULT_HI_A = 3000.0, 3780.0

#: Roman → integer ionization stage. iSpec's `ion` column is 1 for neutral.
_ION_INT = {'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5, 'VI': 6,
            'VII': 7, 'VIII': 8, 'IX': 9, 'X': 10}


class NearUVLinelistError(RuntimeError):
    """Raised when the near-UV list cannot be built. Never returns an empty list."""


def _ensure_ispec_on_path() -> None:
    """iSpec is a source tree, not a pip install — the repo's own convention."""
    from config.constants import ISPEC_DIR
    if str(ISPEC_DIR) not in sys.path:
        sys.path.insert(0, str(ISPEC_DIR))


def _ispec_atomic_dtype():
    """iSpec's own atomic-linelist dtype — read from iSpec, never restated here.

    Restating the 43 fields would create a second definition that silently rots the
    day iSpec changes one. If iSpec is not importable this raises, which is correct:
    without iSpec there is nothing to build a list *for*.
    """
    _ensure_ispec_on_path()
    try:
        from ispec.lines import _get_atomic_linelist_definition
    except ImportError as exc:      # pragma: no cover — environment, not logic
        raise NearUVLinelistError(
            f"iSpec is not importable, so the atomic-linelist definition cannot be "
            f"read from it: {exc}. This builder deliberately has no local copy of the "
            f"format — see the module docstring."
        ) from exc
    return _get_atomic_linelist_definition()


def _element_z_map(chem_elements=None) -> dict[str, int]:
    """symbol → atomic number, from iSpec's chemical-elements table.

    `chem_elements` is the array `ispec.read_chemical_elements` returns (the same one
    the production synth path already loads); pass it in to avoid a second read.
    """
    if chem_elements is None:
        _ensure_ispec_on_path()
        import ispec
        from pipeline.abundances_derive import _SYNTH_CHEM_FILE
        chem_elements = ispec.read_chemical_elements(_SYNTH_CHEM_FILE)
    return {str(r['symbol']): int(r['atomic_num']) for r in chem_elements}


def read_band(raw_path=DEFAULT_RAW, lo_A: float = DEFAULT_LO_A,
              hi_A: float = DEFAULT_HI_A) -> tuple[list[dict], dict]:
    """Parse the VALD extract and keep [lo_A, hi_A). Returns (records, report).

    The parser's failure report is passed straight back — a caller that ignores it is
    ignoring dropped lines, which RYA-429's no-silent-drop rule forbids.
    """
    raw_path = Path(raw_path)
    if not raw_path.exists():
        raise NearUVLinelistError(f"VALD extract not found: {raw_path}")
    records, report = parse_vald_long(str(raw_path))
    band = [r for r in records if lo_A <= r['wavelength'] < hi_A]
    report = dict(report, n_in_band=len(band), lo_A=lo_A, hi_A=hi_A,
                  source=str(raw_path))
    if not band:
        raise NearUVLinelistError(
            f"{raw_path.name} yielded 0 lines in {lo_A:.0f}-{hi_A:.0f} A "
            f"({report['n_parsed']} parsed overall). An empty line list synthesises a "
            f"flat spectrum at every abundance, which reads as a fit rather than a "
            f"failure — refusing to return one.")
    return band, report


def to_ispec_array(records: list[dict], *, chem_elements=None,
                   gf_sources: dict[tuple[str, int, float], str] | None = None,
                   reference_code: str = 'VALD3') -> np.ndarray:
    """VALD records → an iSpec atomic-linelist structured array.

    `gf_sources` optionally maps (element, ion_int, wavelength_A) → the VALD gf tag
    (e.g. 'K14', 'BWL'), as produced by `scripts/ingest_vald_references.py`. It is
    recorded in `reference_code` so the gf provenance travels with the line; it never
    changes a `loggf`.
    """
    dtype = _ispec_atomic_dtype()
    zmap = _element_z_map(chem_elements)

    unknown = sorted({r['element'] for r in records if r['element'] not in zmap})
    if unknown:
        raise NearUVLinelistError(
            f"{len(unknown)} species in the extract are absent from iSpec's "
            f"chemical-elements table and cannot be given an atomic number: "
            f"{', '.join(unknown[:10])}. Dropping them would silently thin the "
            f"blend forest the near-UV synthesis depends on.")

    arr = np.zeros(len(records), dtype=dtype)
    for i, r in enumerate(records):
        z = zmap[r['element']]
        ion_i = _ION_INT.get(r['ion'])
        if ion_i is None:
            raise NearUVLinelistError(
                f"unmapped ionization stage {r['ion']!r} for {r['element']} at "
                f"{r['wavelength']:.4f} A")
        wl = float(r['wavelength'])
        e_lo, e_up = float(r['e_low_eV']), float(r['e_up_eV'])
        tag = None
        if gf_sources is not None:
            tag = gf_sources.get((r['element'], ion_i, round(wl, 4)))

        arr['element'][i] = f"{r['element']} {ion_i}"
        arr['wave_A'][i] = wl
        arr['wave_nm'][i] = wl / 10.0
        arr['loggf'][i] = float(r['log_gf'])
        arr['lower_state_eV'][i] = e_lo
        arr['lower_state_cm1'][i] = e_lo * EV_TO_CM1
        arr['lower_j'][i] = float(r['j_low'])
        arr['upper_state_eV'][i] = e_up
        arr['upper_state_cm1'][i] = e_up * EV_TO_CM1
        arr['upper_j'][i] = float(r['j_up'])
        # Turbospectrum takes the upper statistical weight, not J. VALD gives J.
        arr['upper_g'][i] = 2.0 * float(r['j_up']) + 1.0
        arr['lande_lower'][i] = float(r['lande_lower'])
        arr['lande_upper'][i] = float(r['lande_upper'])
        arr['spectrum_transition_type'][i] = 'GA'
        # VALD's Rad. column is log10(gamma_rad / s^-1); Turbospectrum wants the
        # linear rate. 0.0 means VALD supplied none -> leave 0.0 (TS default), never
        # 10**0 = 1 s^-1, which would be a physically absurd radiative width.
        rad_log = float(r['damping_rad'])
        arr['rad'][i] = rad_log
        arr['turbospectrum_rad'][i] = 10.0 ** rad_log if rad_log > 0.0 else 0.0
        arr['stark'][i] = float(r['damping_stark'])
        vdw = float(r['damping_vdW'])
        arr['waals'][i] = vdw
        arr['waals_single_gamma_format'][i] = vdw
        arr['turbospectrum_fdamp'][i] = vdw
        arr['spectrum_fudge_factor'][i] = 1.0
        arr['theoretical_depth'][i] = float(r['central_depth'])
        arr['theoretical_ew'][i] = 0.0
        arr['lower_orbital_type'][i] = 'X'
        arr['upper_orbital_type'][i] = 'X'
        arr['molecule'][i] = 'F'
        arr['spectrum_synthe_isotope'][i] = 0
        arr['ion'][i] = ion_i
        arr['spectrum_moog_species'][i] = f"{z}.{ion_i - 1}"
        arr['turbospectrum_species'][i] = f"{z}.000000"
        arr['width_species'][i] = f"{z}.{ion_i - 1:02d}"
        arr['reference_code'][i] = (tag or reference_code)[:10]
        arr['nlte'][i] = 'F'
        arr['nlte_level_low'][i] = 0
        arr['nlte_level_up'][i] = 0
        arr['nlte_label_low'][i] = 'none'
        arr['nlte_label_up'][i] = 'none'
        # Built and validated for Turbospectrum only. Claiming support for codes this
        # list has never been run through would be a claim we cannot back; the other
        # writers filter on their own flag and will correctly see an empty list.
        arr['turbospectrum_support'][i] = 'T'
        for f in ('spectrum_support', 'moog_support', 'width_support',
                  'synthe_support', 'sme_support'):
            arr[f][i] = 'F'

    arr.sort(order='wave_A')
    return arr


def band_stats(arr: np.ndarray) -> dict:
    """The facts a reader needs to judge the list, including the ones that look bad."""
    els, counts = np.unique(arr['element'], return_counts=True)
    order = np.argsort(-counts)
    return {
        'n_lines': int(len(arr)),
        'n_species': int(len(els)),
        'lo_A': float(arr['wave_A'].min()),
        'hi_A': float(arr['wave_A'].max()),
        'n_vdw_zero': int(np.sum(arr['turbospectrum_fdamp'] == 0.0)),
        'n_vdw_positive': int(np.sum(arr['turbospectrum_fdamp'] > 0.0)),
        'n_rad_zero': int(np.sum(arr['turbospectrum_rad'] == 0.0)),
        'top_species': [(str(els[i]), int(counts[i])) for i in order[:8]],
    }


def load_gf_sources(path) -> dict[tuple[str, int, float], str]:
    """Read a `vald_gf_sources_*.csv` (scripts/ingest_vald_references.py) into the
    key `to_ispec_array` expects. Absent file → empty map, and the caller records
    that the lines carry the bare 'VALD3' tag."""
    import pandas as pd
    p = Path(path)
    if not p.exists():
        return {}
    df = pd.read_csv(p)
    return {(str(r.element), int(r.ion), round(float(r.wavelength_air_A), 4)):
            str(r.vald_gf_source) for r in df.itertuples() if str(r.vald_gf_source)}


def write(arr: np.ndarray, out_path) -> str:
    """Write via iSpec's own writer, then read it back and check it round-trips.

    The read-back is the point: a file iSpec cannot re-read is a file iSpec cannot
    synthesise from, and that is precisely the failure this whole ticket has been
    chasing. Cheap here, invisible later.
    """
    _ensure_ispec_on_path()
    import ispec
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ispec.write_atomic_linelist(arr, linelist_filename=str(out_path))
    back = ispec.read_atomic_linelist(str(out_path))
    if len(back) != len(arr):
        raise NearUVLinelistError(
            f"round-trip lost lines: wrote {len(arr)}, read back {len(back)} from "
            f"{out_path}")
    dw = float(np.max(np.abs(np.sort(back['wave_A']) - np.sort(arr['wave_A']))))
    if dw > 1e-3:
        raise NearUVLinelistError(
            f"round-trip moved wavelengths by up to {dw:.4g} A in {out_path}")
    return str(out_path)


def build(raw_path=DEFAULT_RAW, lo_A: float = DEFAULT_LO_A, hi_A: float = DEFAULT_HI_A,
          out_path=None, *, gf_sources_csv=None, chem_elements=None) -> tuple[str, dict]:
    """Parse → convert → write → verify. Returns (path, report)."""
    records, report = read_band(raw_path, lo_A, hi_A)
    gf = load_gf_sources(gf_sources_csv) if gf_sources_csv else {}
    arr = to_ispec_array(records, chem_elements=chem_elements, gf_sources=gf or None)
    if out_path is None:
        out_path = (ROOT / 'data' / 'linelists' /
                    f'ispec_nearuv_{int(lo_A)}_{int(hi_A)}' / 'atomic_lines.tsv')
    path = write(arr, out_path)
    report.update(band_stats(arr))
    report['gf_sources_matched'] = int(sum(
        1 for i in range(len(arr)) if str(arr['reference_code'][i]) != 'VALD3'))
    report['path'] = path
    return path, report
