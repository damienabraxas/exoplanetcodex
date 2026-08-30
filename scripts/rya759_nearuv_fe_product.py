#!/usr/bin/env python3
"""
RYA-759 — the near-UV Fe I 1D-LTE product: A(Fe) from flux-fit synthesis, 3000-3780 A.

This is the abundance step the smoke harness never had. Every other route is closed in
this band: interval integration is falsified everywhere, and profile fitting was TESTED
and FALSIFIED here (901 Fe I candidates -> 0 measurable, because the median line gap is
0.146 A, smaller than a strong line's own wings). Synthesis is what is left, and it works
because it never needs an isolated line -- it fits the whole crowded window at once.

METHOD. For each Fe I line, minimise reduced chi2 over a single free A(Fe) with
`_fit_synth_flux` -- the SAME production core that reproduces the optical A(Fe I)=7.520
exactly -- against the Kitt Peak flux atlas, on ATLAS9.Castelli at STAR_PARAMS solar.
Blends are synthesised in-window, never floored.

FOUR THINGS THAT MAKE THIS BAND DIFFERENT, each handled explicitly rather than inherited:

1. NO MEASURED EW EXISTS, so the production wing-wide window rule (half-width from EW)
   CANNOT be applied -- there is nothing to key it on. A fixed half-width is used and
   `--half-width-A` sweeps it, because a window choice that changes the answer is a
   systematic that has to be reported, not a parameter to tune.
2. NO gf SINGLE-SOURCING BELOW 3780 A. `canonical_gf.csv` spans 3780.04-9199.90 A, so
   `apply_canonical_gf` is OFF and the list's own VALD gf is used as delivered. That is a
   declared provenance gap, not a setting.
3. NO MOLECULAR OPACITY in-band from either source (VALD gives no isotopologue code;
   iSpec's molecular lists start at 400 nm). Where that bites, the synthesis sits ABOVE
   the observed pseudo-continuum -- so this records the per-window synth/obs continuum
   ratio as a first-class diagnostic and reports lines whose window is poorly reproduced.
4. THE PSEUDO-CONTINUUM SYSTEMATIC (0.10 dex, RYA-713) DOES NOT AVERAGE DOWN. It is
   carried into the budget as a term, never into the scatter.

⚠️ HYDROGEN IS SYNTHESISED AND THIS BAND IS NOT A HOLE. The record said the merging Balmer
series 3646-3771 A was absent because H I is excluded from our list. The exclusion is real;
the conclusion was wrong. `ispec/synth/turbospectrum.py:283` appends TS's own
`DATA/Hlinedata` on EVERY call, and line 133 drops hydrogen from the atomic list on purpose
("usually encoded into radiative transfer codes"). Measured: the H 13 / H 12 windows come
out at synth/obs = 0.889 / 0.927 against a 1.050 control, with the OBSERVED atlas showing
the same depression (0.36 vs 0.67). So the full 3000-3780 A is in play.

⚠️ THE CONTROL THAT GATES THIS NUMBER PASSES ON COMPENSATING ERRORS. Recorded verdict:
abundance -0.0100 dex (inside +/-0.05) while the EW ratio is 1.433, MAD 0.359. Its own
words: "ABUNDANCE AGREES, EW DOES NOT ... Treat with MORE suspicion than a clean failure,
because it looks like success." That caveat is stamped into this artifact so it cannot be
read off without it.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

from pipeline.nearuv_synth import build_solar_context, gf_provenance   # noqa: E402
# No copy-paste (the Ba->Al lesson): the atlas readers already exist on this ticket.
from rya759_nearuv_synth import _kp_segments, _load_kp_window          # noqa: E402

LO_A, HI_A = 3000.0, 3780.0
NEARUV_LINELIST = ROOT / 'data' / 'linelists' / 'ispec_nearuv_3000_3780' / 'atomic_lines.tsv'
OUT = ROOT / 'data' / 'audit' / 'nearuv_fe_product'

#: RYA-713, ratified. A systematic, so it is a BUDGET TERM and never enters the scatter.
PSEUDO_CONTINUUM_DEX = 0.10

#: A window the synthesis cannot reproduce cannot be trusted to yield an abundance from it.
#: Not a tuned threshold -- it is "the model is within 25% of the observed blanketing",
#: which is loose on purpose, and every line is reported with its ratio either way.
CONTINUUM_RATIO_LO, CONTINUUM_RATIO_HI = 0.75, 1.25


#: Neutral iron as the list spells it. iSpec writes 'Fe 1' (space, arabic), NOT 'Fe I' --
#: my first pass matched the roman form and selected nothing. It failed loudly rather than
#: returning an empty product, which is the only reason that was a two-minute correction.
FE_I = 'Fe 1'

#: Theoretical central depth for candidates. The floor keeps lines actually visible above
#: this band's crowding; the ceiling drops cores so black that chi2 is nearly flat in A(X)
#: -- those return a number without really measuring one.
DEPTH_FLOOR, DEPTH_CEIL = 0.15, 0.90


def species_token(element: str, ion: str) -> str:
    """('Fe', 'I') -> 'Fe 1', as iSpec's linelist actually spells it — RYA-904.

    Roman numerals are what the CLI, the accounting file and every ticket use; 'Fe 1'
    (space, arabic) is what the list uses. Matching the roman form selects NOTHING, which
    RYA-759 hit and caught only because it fails loudly.
    """
    roman = {"I": 1, "II": 2, "III": 3, "IV": 4, "V": 5, "VI": 6}
    ion = str(ion).strip()
    if ion in roman:
        return f"{element} {roman[ion]}"
    if ion.isdigit():
        return f"{element} {int(ion)}"
    raise SystemExit(f"cannot spell ion {ion!r} as a linelist species token; add it to "
                     f"`species_token` rather than guessing.")


def select_lines(linelist: np.ndarray, *, lo_A: float, hi_A: float, n: int,
                 teff: float, min_sep_A: float, species: str = FE_I) -> pd.DataFrame:
    """Candidates for ONE species, strongest-first by the list's THEORETICAL CENTRAL DEPTH.

    🔴 `species` — RYA-904. THIS WAS PINNED TO Fe I AND `--ion` WAS DECORATIVE.
    `derive_band_products.py --ion II` built a Fe context, called this, got Fe I lines,
    and labelled every one of them `Fe II`. Nothing downstream could tell: the
    wavelengths are real, the fits converge, and the product carries an ion it did not
    measure. That is worse than the loader defect this ticket is about, because it
    produces a plausible NUMBER rather than a refusal. The default is unchanged, so
    RYA-759's near-UV selection is bit-for-bit what it was.

    ⚠️ `theoretical_ew` IS ALL ZERO in this list -- checked, not assumed: 4,364 Fe I rows,
    every one 0.0, because our VALD converter does not populate it. `theoretical_depth` IS
    populated (0.001-0.995, median 0.433), so depth is the strength measure here.

    That is also WHY the window below is fixed rather than wing-wide: the production rule
    `_wingwide_window_nm` keys its half-width on an EW, and there is neither a measured one
    (profile fitting was tested and falsified in this band) nor a theoretical one. A fixed
    half-width that gets SWEPT is the honest substitute for a rule that cannot run.

    Minimum separation keeps the set spread across the band instead of piling into one
    crowded complex, so neighbouring fits are not the same photons counted twice.
    """
    names = linelist.dtype.names
    w_A = np.asarray(linelist['wave_A'] if 'wave_A' in names
                     else linelist['wave_nm'] * 10.0, dtype=float)
    el = np.asarray([str(x).strip() for x in linelist['element']])
    m = (w_A >= lo_A) & (w_A <= hi_A) & (el == species)
    if not m.any():                      # never guess a column value silently
        raise SystemExit(f"no {species!r} rows in {lo_A}-{hi_A} A; element values look "
                         f"like {sorted(set(el))[:8]}")
    df = pd.DataFrame({
        'wave_A': w_A[m],
        'loggf': np.asarray(linelist['loggf'], dtype=float)[m],
        'ep_eV': np.asarray(linelist['lower_state_eV'], dtype=float)[m],
        'theo_depth': np.asarray(linelist['theoretical_depth'], dtype=float)[m],
    })
    theta = 5040.0 / float(teff)
    df['strength'] = df['loggf'] - df['ep_eV'] * theta
    df = df[df['theo_depth'].between(DEPTH_FLOOR, DEPTH_CEIL)]
    if df.empty:
        raise SystemExit(f'no {species} line in {lo_A}-{hi_A} A has theoretical depth in '
                         f'[{DEPTH_FLOOR}, {DEPTH_CEIL}] — check the column')
    df = df.sort_values('theo_depth', ascending=False).reset_index(drop=True)

    kept: list[int] = []
    taken: list[float] = []
    for i, wv in enumerate(df['wave_A'].to_numpy()):
        if all(abs(wv - t) >= min_sep_A for t in taken):
            kept.append(i)
            taken.append(float(wv))
        if len(kept) >= n:
            break
    return df.iloc[kept].sort_values('wave_A').reset_index(drop=True)


def fit_one(ctx: dict, segs, wave_A: float, hw_A: float, tmp_dir: str,
            load=None, *, nlte_deck=None, nlte_deck_key=None,
            atmosphere_layers_file=None, atmosphere=None) -> dict:
    """Flux-fit A(Fe) in one window, plus the continuum diagnostic for that window.

    🔴 `load` — RYA-904. THE OBSERVED SPECTRUM WAS HARD-PINNED TO KITT PEAK HERE.
    `derive_band_products.synthesis_route` takes an `--instrument` argument, tags the
    product with it, and then called this function, which read the Kitt Peak atlas
    whatever the argument said. So `--instrument crires_plus` would have produced a
    product LABELLED CRIRES+ and MEASURED ON KITT PEAK — the same defect the ticket
    exists to fix, one level deeper, and silent because the KPNO atlas does cover the
    IR windows and the fit would have converged perfectly happily.

    `load(centre, pad) -> (wave_A, flux, provenance)` supplies the observed window
    instead. The DEFAULT IS UNCHANGED: `None` reads Kitt Peak through this module's own
    reader exactly as before, so RYA-759's published near-UV values cannot move by way
    of this argument.

    🔴 RYA-1044 — `nlte_deck` / `nlte_deck_key` / `atmosphere_layers_file` / `atmosphere`
    exist so the SYNTHESIS route can run an Engine-B leg over these same lines, WITHOUT a
    second copy of this function. RYA-701 is the reason: one Ba->Al copy of a fitting
    routine produced thirteen defects, and a copy here would additionally let the
    Engine-B leg drift from the 1D-LTE leg it is supposed to be differenced against --
    which would silently corrupt the very quantity the pair exists to measure.

    ALL FOUR DEFAULT TO None AND CHANGE NOTHING WHEN UNSET. With them unset this is
    character-for-character the call RYA-759 published against, and `atmosphere` falls
    back to `ctx['atmosphere']` exactly as before. They are keyword-only so no positional
    caller can acquire one by accident. (Measured, not assumed: over 30 probes spanning the near-UV,
    red-optical and NIR, including six segment seams, this reader and
    `measure_band_ew.load_kp_window` return bit-identical arrays — so the holding
    dispatch that now supplies `load` is the same data by a different door.)
    """
    from pipeline.abundances_derive import _fit_synth_flux

    lo_A, hi_A = wave_A - hw_A, wave_A + hw_A
    obs_source = "kitt peak atlas segments (default reader)"
    try:
        if load is None:
            ow_A, of = _load_kp_window(segs, wave_A, hw_A + 0.4)
        else:
            ow_A, of, obs_source = load(wave_A, hw_A + 0.4)
    except LookupError as e:
        return {'status': 'no_atlas', 'reason': str(e), 'a_synth': np.nan,
                'red_chi2': np.nan, 'cont_ratio': np.nan, 'obs_source': obs_source}

    a_solar = float(ctx['solar_A'])
    # RYA-1044: only present when the caller supplied them, so an unset call reaches
    # `_fit_synth_flux` with exactly the arguments it reached before they existed.
    _extra = {k: v for k, v in (
        ("nlte_deck", nlte_deck), ("nlte_deck_key", nlte_deck_key),
        ("atmosphere_layers_file", atmosphere_layers_file)) if v is not None}
    r = _fit_synth_flux(
        ow_A / 10.0, np.asarray(of, dtype=float),
        ctx['atmosphere'] if atmosphere is None else atmosphere,
        ctx['teff'], ctx['logg'], ctx['feh'], ctx['vturb'],
        ctx['linelist'], ctx['isotopes'], ctx['solar_abund'], ctx.get('element', 'Fe'),
        int(ctx['atom_code']), lo_A / 10.0, hi_A / 10.0,
        max(a_solar - 3.0, 1.0), a_solar + 5.0,
        float(ctx['resolving_power']), float(ctx['macroturbulence']),
        # 🔴 RYA-1043/1040: THE DECK WAS NEVER FORWARDED. `_fit_synth_flux` has taken
        # `nlte_deck` all along, but this call omitted it, so every synthesis-route fit
        # was LTE whatever `--engine-b-deck` said. Combined with `main()` returning into
        # `synthesis_route` BEFORE the EW route's Engine-B block, the flag parsed, was
        # accepted, and did nothing — which is why the Gerber column is empty on every
        # graded and deep-graded product we hold.
        # RYA-1050 residual: forward the stage when the context knows one. Today this
        # script has no ion axis at all -- it is the near-UV Fe I product and `ctx` is
        # built by `build_solar_context('Fe', ...)`, which sets no 'ion' -- so this is
        # None and the call is unchanged. It is wired anyway because this script calls
        # `_fit_synth_flux` DIRECTLY and runs no NLTE guard of its own: if an ion axis is
        # ever added here, the stage must reach the guard or an Fe II leg would clear on
        # Fe I's labels.
        float(ctx['vsini']), tmp_dir=tmp_dir, ion=ctx.get('ion'), **_extra)

    # How well is this window's blanketing reproduced at all? A window the model cannot
    # reproduce is not a window to take an abundance from, and the ratio says which.
    cont = np.nan
    try:
        from pipeline.nearuv_synth import synthesize_band
        s = synthesize_band(ctx, lo_A, hi_A, element='Fe',
                            trial_A=(r['A_X'] if np.isfinite(r.get('A_X', np.nan))
                                     else a_solar),
                            step_A=0.02, check_sensitivity=False)
        obs_i = np.interp(s['wave_A'], ow_A, of)
        cont = float(np.median(s['flux']) / np.median(obs_i))
    except Exception:
        pass

    return {'status': r['status'], 'reason': r.get('reason', ''),
            # RYA-904 — WHICH SPECTRUM this abundance was fitted against. It travels out
            # of the fitter because the fitter is the only thing that knows it, and a
            # product that cannot name its own observed source is exactly what the
            # instrument-decorative bug above looked like from the outside.
            'obs_source': obs_source,
            'a_synth': float(r.get('A_X', np.nan)),
            'red_chi2': float(r.get('red_chi2', np.nan)),
            'n_pix': int(r.get('n_pix', 0)), 'cont_ratio': cont,
            # RYA-847: pass the constraint metrics straight through. `_fit_synth_flux`
            # measures them; this harness is the only thing between the fitter and the
            # band product, and anything it does not forward is a quantity nobody
            # downstream can gate on -- which is exactly what happened to `red_chi2`
            # before RYA-843 went looking for it.
            'sigma_A': float(r.get('sigma_A', np.nan)),
            'frac_rise_weaker': float(r.get('frac_rise_weaker', np.nan)),
            'edge_distance_dex': float(r.get('edge_distance_dex', np.nan))}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--limit', type=int, default=40, help='max Fe I lines to fit')
    ap.add_argument('--half-width-A', type=float, default=0.40,
                    help='constant fit half-width in A. The production EW-keyed rule '
                         'cannot run here (no measured OR theoretical EW exists), so this '
                         'is fixed and MUST be swept — it is a systematic, not a knob.')
    ap.add_argument('--min-sep-A', type=float, default=4.0)
    ap.add_argument('--lo-A', type=float, default=LO_A)
    ap.add_argument('--hi-A', type=float, default=HI_A)
    ap.add_argument('--resolving-power', type=float, default=None)
    ap.add_argument('--tag', default='main')
    ap.add_argument('--out', default=str(OUT))
    a = ap.parse_args()

    if a.resolving_power is None:
        cat = pd.read_csv(ROOT / 'data' / 'catalog' / 'instrument_catalog.csv')
        row = cat[cat.iloc[:, 0].astype(str) == 'kpno_solar_atlas']
        if row.empty:
            raise SystemExit('kpno_solar_atlas absent from data/catalog/instrument_catalog.csv')
        a.resolving_power = float(row.iloc[0]['resolving_power_max'])

    if not NEARUV_LINELIST.exists():
        raise SystemExit(f'near-UV linelist missing: {NEARUV_LINELIST}\n'
                         f'build it with pipeline.nearuv_linelist.build()')

    prov = gf_provenance(a.lo_A, a.hi_A)
    print(f'[gf]  {prov["detail"]}')
    ctx = build_solar_context('Fe', a.resolving_power,
                              linelist_file=str(NEARUV_LINELIST),
                              apply_canonical_gf=prov['apply_canonical_gf'])
    print(f'[atm] ATLAS9.Castelli Teff={ctx["teff"]:.0f} logg={ctx["logg"]:.3f} '
          f'[Fe/H]={ctx["feh"]:.2f} xi={ctx["vturb"]:.2f}')
    print(f'[brd] R={ctx["resolving_power"]:.0f} vmac={ctx["macroturbulence"]} '
          f'vsini={ctx["vsini"]}  a_start(iSpec internal)={ctx["solar_A"]:.3f}')

    segs = _kp_segments()
    print(f'[kp]  {len(segs)} atlas segments')

    cand = select_lines(ctx['linelist'], lo_A=a.lo_A, hi_A=a.hi_A, n=a.limit,
                        teff=ctx['teff'], min_sep_A=a.min_sep_A)
    hw_note = f'constant +/-{a.half_width_A:.2f} A (no EW exists to key the production rule)'
    print(f'[sel] {len(cand)} Fe I lines, {cand["wave_A"].min():.1f}-'
          f'{cand["wave_A"].max():.1f} A, theoretical depth '
          f'{cand["theo_depth"].min():.3f}-{cand["theo_depth"].max():.3f}')
    print(f'[win] {hw_note}\n')

    # A private tmp dir: the shared default is a race when two synthesis jobs overlap.
    # It must be CREATED -- the fitter writes into it and does not mkdir, so a missing
    # directory surfaces as "synthesis error: No such file or directory" on every line.
    tmp_dir = f'/tmp/ispec_nearuv_{a.tag}_{os.getpid()}'
    os.makedirs(tmp_dir, exist_ok=True)
    rows = []
    for i, rec in cand.iterrows():
        wv = float(rec['wave_A'])
        out = fit_one(ctx, segs, wv, a.half_width_A, tmp_dir)
        rows.append({'wave_A': wv, 'loggf': float(rec['loggf']),
                     'ep_eV': float(rec['ep_eV']), 'strength': float(rec['strength']),
                     'theo_depth': float(rec['theo_depth']),
                     'half_width_A': float(a.half_width_A), **out})
        flag = '' if CONTINUUM_RATIO_LO <= (out['cont_ratio'] or 0) <= CONTINUUM_RATIO_HI \
            else '  [continuum off]'
        print(f'  {i + 1:>3}/{len(cand)}  {wv:9.3f}  {out["status"]:<12} '
              f'A={out["a_synth"]:.3f}  chi2r={out["red_chi2"]:.2f}  '
              f'cont={out["cont_ratio"]:.3f}{flag}', flush=True)

    df = pd.DataFrame(rows)
    ok = df[(df['status'] == 'ok') & np.isfinite(df['a_synth'])].copy()
    clean = ok[ok['cont_ratio'].between(CONTINUUM_RATIO_LO, CONTINUUM_RATIO_HI)]

    def _agg(d: pd.DataFrame) -> dict:
        if d.empty:
            return {'n': 0, 'median': None, 'mad': None, 'sem': None}
        v = d['a_synth'].to_numpy(float)
        med = float(np.median(v))
        mad = float(1.4826 * np.median(np.abs(v - med)))
        return {'n': int(v.size), 'median': med, 'mad': mad,
                'sem': float(mad / np.sqrt(v.size)) if v.size else None}

    all_stat, clean_stat = _agg(ok), _agg(clean)
    print(f'\n{"":<26}{"n":>5}{"A(Fe I)":>10}{"scatter":>10}{"SEM":>8}')
    for name, s in (('all converged', all_stat), ('continuum-clean', clean_stat)):
        if s['n']:
            print(f'  {name:<24}{s["n"]:>5}{s["median"]:>10.3f}'
                  f'{s["mad"]:>10.3f}{s["sem"]:>8.3f}')

    art = {
        'ticket': 'RYA-759', 'product': 'Fe I near-UV 1D-LTE',
        'band_A': [a.lo_A, a.hi_A], 'instrument': 'kpno_solar_atlas',
        'engine': 'Turbospectrum via iSpec', 'model_grid': ctx['model_grid'],
        'method': 'flux-fit synthesis (_fit_synth_flux), blend-aware, LTE',
        'resolving_power': ctx['resolving_power'],
        'half_width_mode': hw_note,
        'gf_provenance': prov['detail'],
        'pseudo_continuum_systematic_dex': PSEUDO_CONTINUUM_DEX,
        'all_converged': all_stat, 'continuum_clean': clean_stat,
        'n_attempted': int(len(df)),
        'status_counts': {k: int(v) for k, v in df['status'].value_counts().items()},
        'caveats': [
            'GATE PASSES ON COMPENSATING ERRORS: the SynthesisHandler optical control is '
            '-0.0100 dex on abundance but 1.433 on EW ratio (MAD 0.359). Its own verdict: '
            '"ABUNDANCE AGREES, EW DOES NOT ... Treat with MORE suspicion than a clean '
            'failure, because it looks like success." This number inherits that.',
            'NO gf single-sourcing below 3780 A (canonical_gf.csv starts at 3780.04 A); '
            'VALD gf used as delivered, accuracy unadjudicated for 127 source tags.',
            'NO molecular opacity in-band from either source (VALD has no isotopologue '
            'code; iSpec molecular lists start at 400 nm).',
            f'Pseudo-continuum systematic {PSEUDO_CONTINUUM_DEX} dex does NOT average '
            'down and is NOT in the scatter above.',
            'No measured EW exists in this band, so the production EW-keyed window rule '
            'cannot apply; half-width is fixed and must be swept.',
            'LTE only. UV Fe I is heavily over-ionized, so the missing NLTE correction is '
            'large and POSITIVE -- a low value here is expected physics, not a bug.',
        ],
    }
    outdir = Path(a.out)
    outdir.mkdir(parents=True, exist_ok=True)
    df.to_csv(outdir / f'nearuv_fe_per_line_{a.tag}.csv', index=False)
    (outdir / f'nearuv_fe_product_{a.tag}.json').write_text(json.dumps(art, indent=2))
    print(f'\n[out] {outdir}/nearuv_fe_product_{a.tag}.json')


if __name__ == '__main__':
    main()
