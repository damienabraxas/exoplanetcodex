"""
pipeline/nearuv_synth.py — near-UV LTE synthesis on the PRODUCTION atmosphere path.
==================================================================================
RYA-759 Move 2. The near-UV synthesis runs on the same model atmosphere as the optical
Fe production run — ATLAS9.Castelli, interpolated by iSpec to the Sun's STAR_PARAMS —
with Turbospectrum driven THROUGH iSpec. Not because the raw-binary route is broken
(Move 1 named its fault exactly), but because a near-UV-vs-optical comparison made on
two different atmospheres measures the atmospheres as much as the band.

    optical production            near-UV (this module)
    ------------------            ---------------------
    ATLAS9.Castelli               ATLAS9.Castelli          <- same grid
    STAR_PARAMS['solar']          STAR_PARAMS['solar']     <- same params
    ispec.generate_spectrum       ispec.generate_spectrum  <- same driver
    GES v6, 420-920 nm            near-UV VALD, 300-378 nm <- the ONLY difference

The Gerber raw `babsma_lu`/`bsyn_lu` deck on a MARCS grid node remains the NLTE-anchor
validation gate (`scripts/ts_gerber_gate.py`) and nothing else.

THE THREE GUARDS, AND WHY THESE THREE
-------------------------------------
Every one is a failure this project has actually shipped, silently:

1. **All-zero / empty flux** — RYA-506 and RYA-682: iSpec writes a zero-row artifact
   and exits 0. A zero-flux window fits an abundance from nothing.
2. **Zero lines in band** — RYA-713 at 4065.381 Å: below the GES list's blue limit the
   synthesis is IDENTICAL at every abundance, so chi2 is flat and the "fit" returns
   noise that looks like a measurement. This is the near-UV's most likely failure, since
   the whole point is that the usual list does not reach here.
3. **Model interpolation** — Move 1: the model atmosphere is where the 0-byte spectrum
   actually came from. It is checked before anything downstream is believed.

Each RAISES. The RYA-765 tracer records it as well, so a run leaves a timeline; the
tracer is observability, never the enforcement — a `trace_check` that only warns is the
warn-only degrade RYA-765 already caught being swallowed.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.intake_debug import trace_asset, trace_check, trace_fallback  # noqa: E402

#: Peak |Δflux| over a ±0.5 dex abundance swing below which the band is not usably
#: represented in the line list — i.e. the synthesis is not actually using it.
#:
#: RYA-772: the per-line measurement handler on the RYA-713 line carries the same
#: constant under the same name. It is defined here rather than imported because that
#: module is not on main and this ticket is near-UV only. When RYA-713 lands, the two
#: must be reconciled to one definition — a duplicated threshold is exactly the drift
#: the no-copy-paste rule exists to stop, and this note is the marker for it.
MIN_SENSITIVITY = 0.002

#: canonical_gf.csv spans 3780.0383–9199.897 Å. Below its blue edge there is no
#: canonical entry to single-source against, and `apply_to_synth_array` raises rather
#: than defaulting (RYA-353). Measured from the file, stated here so the near-UV path
#: can refuse to *pretend* it had canonical gf.
CANONICAL_GF_LO_A = 3780.0


#: iSpec returns exactly 1e-10 at the FIRST and LAST pixel of every synthesis window.
#: Measured, and not a near-UV quirk: the optical GES path at 5500-5502 A does the same
#: (RYA-759). It is a boundary sentinel, not flux. Left in, it makes the deepest point of
#: every window the window's own edge, and any integration over the full window carries
#: two spurious zero-flux pixels. Trimmed here, counted, and reported.
EDGE_SENTINEL = 1e-8


class NearUVSynthesisError(RuntimeError):
    """A near-UV synthesis that cannot be believed. Never returns a quiet zero."""


def trim_edge_sentinels(wave_A, flux):
    """Drop iSpec's leading/trailing 1e-10 boundary pixels. Returns (wave, flux, n)."""
    f = np.asarray(flux, dtype=float)
    w = np.asarray(wave_A, dtype=float)
    lo, hi = 0, len(f)
    while lo < hi and (not np.isfinite(f[lo]) or abs(f[lo]) <= EDGE_SENTINEL):
        lo += 1
    while hi > lo and (not np.isfinite(f[hi - 1]) or abs(f[hi - 1]) <= EDGE_SENTINEL):
        hi -= 1
    return w[lo:hi], f[lo:hi], int(lo + (len(f) - hi))


def assert_atmosphere(atmosphere, *, teff: float, logg: float, feh: float,
                      vturb: float, model_grid: str) -> None:
    """Guard 3. An atmosphere that is None/empty means the interpolation failed."""
    ok = atmosphere is not None and len(np.atleast_1d(atmosphere)) > 0
    trace_asset('model_atmosphere', ok, path=model_grid,
                detail=f"{model_grid} at Teff={teff:.0f} logg={logg:.3f} "
                       f"[Fe/H]={feh:.2f} xi={vturb:.2f}")
    if not ok:
        raise NearUVSynthesisError(
            f"model atmosphere interpolation returned nothing for {model_grid} at "
            f"Teff={teff:.0f}, logg={logg:.3f}, [Fe/H]={feh:.2f}. Synthesising against "
            f"an empty atmosphere is the Move 1 failure in a different costume.")
    n = len(np.atleast_1d(atmosphere))
    trace_check('model_atmosphere_layers', n > 1, detail=f"{n} layers", n_layers=n)
    if n <= 1:
        raise NearUVSynthesisError(
            f"{model_grid} interpolated to {n} layer(s) at Teff={teff:.0f}, "
            f"logg={logg:.3f} — not an atmosphere.")


def assert_linelist_covers(linelist, lo_A: float, hi_A: float, *,
                           element: str | None = None) -> int:
    """Guard 2, static half: does the list contain anything in the band at all?

    Returns the in-band line count. Zero raises — a synthesis over a band the list
    does not reach produces the same spectrum at every abundance.
    """
    path_note = f"{lo_A:.1f}-{hi_A:.1f} A"
    n = 0
    if linelist is not None and len(linelist):
        w = np.asarray(linelist['wave_A'], dtype=float)
        m = (w >= lo_A) & (w <= hi_A)
        n = int(m.sum())
    trace_asset('synthesis_linelist', n > 0, path=path_note,
                detail=f"{n} lines in band")
    if n == 0:
        raise NearUVSynthesisError(
            f"the synthesis line list contains 0 lines in {path_note}. Every trial "
            f"abundance would synthesise an IDENTICAL spectrum, so a fit here returns "
            f"noise shaped like a measurement (RYA-713 measured exactly this at "
            f"4065.381 A, below the GES list's 4200 A blue limit). Coverage gap.")
    n_el = None
    if element is not None and len(linelist):
        pref = f"{element} "
        n_el = int(np.sum([str(e).startswith(pref) for e in linelist['element'][
            (np.asarray(linelist['wave_A'], dtype=float) >= lo_A) &
            (np.asarray(linelist['wave_A'], dtype=float) <= hi_A)]]))
        trace_check(f'linelist_{element}_in_band', n_el > 0,
                    detail=f"{n_el} {element} lines in {path_note}", n_lines=n_el)
        if n_el == 0:
            raise NearUVSynthesisError(
                f"the line list covers {path_note} ({n} lines) but contains no "
                f"{element} line there, so no {element} abundance can be fitted from "
                f"this window.")
    return n


def assert_usable_flux(flux, *, where: str) -> np.ndarray:
    """Guard 1. Empty, all-zero, or non-finite flux is a failure, not a spectrum."""
    f = np.asarray(flux, dtype=float)
    if f.size == 0:
        trace_fallback('synth_empty', f"{where}: generate_spectrum returned 0 points",
                       severity='ERROR')
        raise NearUVSynthesisError(
            f"{where}: the synthesis returned an EMPTY spectrum. This is the 0-byte "
            f"failure shape — it must raise, never be averaged.")
    n_finite = int(np.sum(np.isfinite(f)))
    if n_finite == 0:
        trace_fallback('synth_all_nan', f"{where}: {f.size} points, none finite",
                       severity='ERROR')
        raise NearUVSynthesisError(
            f"{where}: all {f.size} synthesised flux points are non-finite.")
    if not np.any(f[np.isfinite(f)] != 0.0):
        trace_fallback('synth_all_zero',
                       f"{where}: {f.size} points, every one zero (RYA-506 class)",
                       severity='ERROR')
        raise NearUVSynthesisError(
            f"{where}: the synthesis returned an ALL-ZERO spectrum over {f.size} "
            f"points. RYA-506/682: iSpec writes a zero artifact and exits 0, so this "
            f"is checked rather than trusted.")
    trace_check('synth_flux_usable', True,
                detail=f"{where}: {f.size} pts, {n_finite} finite, "
                       f"min={np.nanmin(f):.4f} max={np.nanmax(f):.4f}")
    return f


def assert_sensitive(flux_lo, flux_hi, *, where: str,
                     min_sensitivity: float = MIN_SENSITIVITY) -> float:
    """Guard 2, dynamic half: does the band RESPOND to the abundance?

    The static check proves lines exist in the list; this proves the synthesiser
    actually used them. A band that does not move under a ±0.5 dex swing cannot carry
    an abundance no matter how many lines the list claims.
    """
    lo = np.asarray(flux_lo, dtype=float)
    hi = np.asarray(flux_hi, dtype=float)
    sens = float(np.nanmax(np.abs(hi - lo))) if lo.size and lo.size == hi.size else 0.0
    ok = sens >= min_sensitivity
    trace_check('abundance_sensitivity', ok,
                detail=f"{where}: peak |dflux| = {sens:.5f} over a +/-0.5 dex swing "
                       f"(floor {min_sensitivity})", sensitivity=sens)
    if not ok:
        raise NearUVSynthesisError(
            f"{where}: a +/-0.5 dex abundance swing moves the synthetic flux by at "
            f"most {sens:.5f} (floor {min_sensitivity}). The band is present in the "
            f"list but not synthesised from it — a flat chi2, which fits to noise.")
    return sens


def gf_provenance(lo_A: float, hi_A: float) -> dict:
    """State the gf single-sourcing status for a band instead of inferring it.

    Returns the flag `_load_synth_resources` must be called with, and the reason —
    so the reason travels into the run's provenance rather than living in a comment.
    """
    below = lo_A < CANONICAL_GF_LO_A
    if below:
        detail = (f"canonical_gf.csv starts at {CANONICAL_GF_LO_A:.1f} A; this band "
                  f"starts at {lo_A:.1f} A, so canonical gf single-sourcing (RYA-353) "
                  f"is NOT available and the list's own VALD log gf is used as "
                  f"delivered. Per-line VALD gf SOURCES are recorded, but no canonical "
                  f"adjudication exists below {CANONICAL_GF_LO_A:.1f} A.")
        trace_fallback('canonical_gf_unavailable', detail, severity='WARN',
                       lo_A=lo_A, canonical_lo_A=CANONICAL_GF_LO_A)
    else:
        detail = "canonical gf single-sourcing (RYA-353) applied"
    return {'apply_canonical_gf': not below, 'detail': detail}


def build_solar_context(element: str, resolving_power: float, *,
                        linelist_file: str = None,
                        apply_canonical_gf: bool = True,
                        tmp_dir: str = '/tmp/ispec_codex_synth') -> dict:
    """Assemble the synthesis context from the pipeline's own loaders — never invented.

    Every value comes from an existing production source:
      * stellar parameters from `STAR_PARAMS` via `get_star_params` (config/stars.yaml);
      * the model atmosphere from `_load_atmosphere` on **ATLAS9.Castelli**, named
        explicitly rather than left to a default — a model silently mismatched to its
        flag is what produced this ticket's 0-byte spectra in the first place;
      * broadening from `_resolve_broadening`, which FAILS LOUD rather than inheriting
        the Sun's profile (RYA-288);
      * the linelist through `_load_synth_resources`, so the near-UV list and the
        optical GES list are loaded by one function.

    Two refusals are deliberate and are the reason this is a function rather than a
    dict literal:

    1. **`xi` is never defaulted.** Microturbulence sets the saturation regime and
       therefore biases A(X) itself, not merely its error bar. A star without one
       raises.
    2. **`a_start` is NOT our banked answer.** It is seeded from iSpec's own internal
       solar composition, which knows nothing about the value any control is testing.
       Seeding from the banked number makes a flat χ² return that number exactly and
       score as a perfect result.
    """
    from config.constants import get_star_params
    from pipeline.abundances_derive import (_ATLAS9, _ISPEC_SOLAR_ABUND_FILE,
                                            _load_atmosphere, _load_synth_resources,
                                            _resolve_broadening)
    from pipeline.cno_synthesis import _atom_codes, _solar_A as _ispec_solar_A_map
    import ispec

    p = get_star_params('solar')
    if 'xi' not in p:
        raise NearUVSynthesisError(
            "no microturbulence ('xi') in STAR_PARAMS for the Sun. Refusing to default "
            "it — xi sets the saturation regime and biases the fitted abundance.")
    teff, logg = float(p['teff']), float(p['logg'])
    feh = float(p.get('feh', p.get('feh_ref', 0.0)))
    vturb = float(p['xi'])

    # R is DISCARDED from the resolver on purpose: resolving power is an instrument
    # constant and the caller says which instrument this run is against.
    _r_unused, vmac, vsini, fit_vmac = _resolve_broadening('solar')
    if fit_vmac:
        raise NearUVSynthesisError(
            "STAR_PARAMS asks for vmac='fit' on the Sun. Fitting vmac by full-profile "
            "chi2 rails to the upper bound at high SNR (RYA-309); this path will not "
            "do it. Set a fixed solar vmac or adopt a literature RT value.")

    linelist, isotopes, chem = _load_synth_resources(
        linelist_file=linelist_file, apply_canonical_gf=apply_canonical_gf)
    solar_abund = ispec.read_solar_abundances(_ISPEC_SOLAR_ABUND_FILE)
    atm = _load_atmosphere(teff, logg, feh, vturb, model_grid=_ATLAS9)
    codes = _atom_codes([element], chem, solar_abund)

    return dict(atmosphere=atm, teff=teff, logg=logg, feh=feh, vturb=vturb,
                linelist=linelist, isotopes=isotopes, solar_abund=solar_abund,
                chem_elements=chem, atom_code=int(codes[element]),
                resolving_power=float(resolving_power),
                macroturbulence=float(vmac), vsini=float(vsini),
                model_grid=_ATLAS9, tmp_dir=tmp_dir,
                solar_A=float(_ispec_solar_A_map([element], chem, solar_abund)[element]))


def synthesize_band(context: dict, lo_A: float, hi_A: float, *, element: str,
                    trial_A: float, step_A: float = 0.02,
                    check_sensitivity: bool = True) -> dict:
    """Synthesise [lo_A, hi_A] at `trial_A` through the production core, guarded.

    `context` is the dict `scripts.control_synthesis_handler.build_context` returns —
    the same context the optical control uses, so this shares its atmosphere, its
    broadening resolution and its refusal to default `xi`.

    Returns {'wave_A', 'flux', 'n_lines_in_band', 'sensitivity'}.
    """
    from pipeline.abundances_derive import _synth_flux_at_abund

    assert_atmosphere(context['atmosphere'], teff=context['teff'],
                      logg=context['logg'], feh=context['feh'],
                      vturb=context['vturb'],
                      model_grid=context.get('model_grid', 'ATLAS9.Castelli'))
    n_band = assert_linelist_covers(context['linelist'], lo_A, hi_A, element=element)

    wave_A = np.arange(lo_A, hi_A + 0.5 * step_A, step_A)
    kw = dict(atmosphere=context['atmosphere'], teff=context['teff'],
              logg=context['logg'], feh=context['feh'], vturb=context['vturb'],
              linelist=context['linelist'], isotopes=context['isotopes'],
              solar_abund=context['solar_abund'], element=element,
              atom_code=context['atom_code'],
              R=float(context['resolving_power']),
              macroturbulence=float(context['macroturbulence']),
              vsini=float(context['vsini']),
              tmp_dir=context.get('tmp_dir', '/tmp/ispec_codex_synth'))

    raw = assert_usable_flux(
        _synth_flux_at_abund(wave_A / 10.0, trial_A=float(trial_A), **kw),
        where=f"{element} {lo_A:.1f}-{hi_A:.1f} A @ A={trial_A:.3f}")
    wave_A, flux, n_trim = trim_edge_sentinels(wave_A, raw)
    trace_check('edge_sentinels_trimmed', True,
                detail=f"dropped {n_trim} boundary pixel(s) at iSpec's {EDGE_SENTINEL:g} "
                       f"sentinel", n_trimmed=n_trim)
    if flux.size == 0:
        raise NearUVSynthesisError(
            f"{element} {lo_A:.1f}-{hi_A:.1f} A: every pixel was a boundary sentinel — "
            f"nothing was synthesised.")

    sens = float('nan')
    if check_sensitivity:
        w_nm = wave_A / 10.0
        f_lo = _synth_flux_at_abund(w_nm, trial_A=float(trial_A) - 0.5, **kw)
        f_hi = _synth_flux_at_abund(w_nm, trial_A=float(trial_A) + 0.5, **kw)
        # Trim by INDEX, not by value: a genuinely saturated core can sit near zero, and
        # comparing differently-trimmed arrays would compare different wavelengths.
        sens = assert_sensitive(f_lo[1:-1], f_hi[1:-1],
                                where=f"{element} {lo_A:.1f}-{hi_A:.1f} A")

    return {'wave_A': wave_A, 'flux': flux, 'n_lines_in_band': n_band,
            'sensitivity': sens, 'n_edge_trimmed': n_trim}
