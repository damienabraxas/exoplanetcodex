"""
scripts/rya346_fe2_core_absorber_diag.py
=========================================
RYA-346 — Fe II core-residual diagnostic. RYA-344 was an honest STOP: the obs−synth
residual for the 5 χ²ᵣ-failing solar Fe II lines is centred on the Fe II line CORE
(not a separable nearby blend), and the flux fit floats Fe II HIGH yet the core
deficit persists. An abundance-insensitive core deficit is a SHAPE problem: either a
non-negligible absorber sitting on the core that is absent from our extraction, or
the Fe II line's own loggf/damping.

This script tests the first hypothesis with 5 targeted VALD "Extract All" pulls
(±0.2 Å, low threshold), one per Fe II line, filed at data/linelists/raw/fe2_diag_<wl>.txt.

Per window:
  1. RYA-342/344 flux-space fit → χ²ᵣ, A(Fe II), and the residual peak λ vs the core.
  2. From the VALD pull, every transition within ±0.05 (→±0.10) Å of the core that is
     ABSENT from atomic_lines.tsv (species+λ+EP match via the RYA-345 normalizer),
     ranked by VALD solar central_depth.
  3. Decisive test: add the strongest absent line to an IN-MEMORY scratch copy of the
     synth line list (atomic_lines.tsv on disk is NEVER touched) and re-fit. χ²ᵣ and
     A(Fe II) before/after.

Acceptance (a candidate "explains" the residual) — BOTH must hold:
  • χ²ᵣ < 10 after the addition, AND
  • A(Fe II) falls back toward the ~7.50 consensus (no longer propped high).
A candidate that drops χ²ᵣ but leaves Fe II high just absorbed the error differently.

DIAGNOSTIC ONLY: canonical atomic_lines.tsv untouched; no guessed absorber; every
dropped/unmatched line logged. Params from constants.py (STAR_PARAMS['solar']). Å throughout.
"""
import sys, glob, os
from pathlib import Path
import numpy as np, pandas as pd

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO))
sys.path.insert(0, str(_REPO / 'data' / 'linelists'))
sys.path.insert(0, str(_REPO / 'scripts'))

import vald_parse
from pipeline.species import species_key
import pipeline.abundances_derive as ad
from rya344_recover_fe2_blends import _solar_context, _fit_line, CLEAN_POOL_A, _TMP

_RAW = _REPO / 'data' / 'linelists' / 'raw'
_OUT = _REPO / 'data' / 'results'; _OUT.mkdir(parents=True, exist_ok=True)

# RYA-343 recoverable set → pull file + observed EW (mÅ) for the fit window sizing
TARGETS = [  # (wave_air_A, ew_obs_mA, pull file)
    (5234.623, 120.7, 'fe2_diag_5234.txt'),
    (5991.371,  44.3, 'fe2_diag_5991.txt'),
    (6084.102,  28.5, 'fe2_diag_6084.txt'),
    (6247.557,  62.9, 'fe2_diag_6247.txt'),
    (6456.380,  85.4, 'fe2_diag_6456.txt'),
]
CORE_TOL_TIGHT = 0.05   # Å — first-pass core tolerance
CORE_TOL_WIDE  = 0.10   # Å — widen if the tight pass finds nothing
TRIVIAL_CD     = 0.01   # central_depth below which a line cannot fill an abundance-
                        # insensitive core deficit (flagged, but still scratch-tested
                        # as a control if it is the strongest absent line)


def intake_table():
    print("\n=== INTAKE TABLE (VALD Extract All, ±0.2 Å, filed to data/linelists/raw/) ===")
    print(f"{'file':<20}{'region (Å)':<26}{'sel/proc':<12}{'parse':<14}{'coverage':<10}verdict")
    rows = []
    for wave_A, _ew, fn in TARGETS:
        p = _RAW / fn
        hdr = open(p, errors='replace').readline().strip()
        reg = hdr.split('Wavelength')[0].strip()
        lo, hi = float(reg.split(',')[0]), float(reg.split(',')[1])
        sel, proc = reg.split(',')[2].strip(), reg.split(',')[3].strip()
        recs, info = vald_parse.parse_vald_long(str(p))
        not_trunc = hdr.split(',')[0].strip()[:4].isdigit()      # line 1 = VALD region header
        # window centred on the VALD catalog λ (~±0.003 Å off the RYA-343 target); the
        # diagnostic only needs the ±0.18 Å core region — require that, not exactly ±0.2.
        covers = (lo <= wave_A - 0.18 + 1e-6) and (hi >= wave_A + 0.18 - 1e-6)
        ok = (info['n_failures'] == 0) and not_trunc and covers
        _region = f"{lo:.2f}-{hi:.2f}"
        _selproc = f"{sel}/{proc}"
        _parse = f"{info['n_parsed']} ok, {info['n_failures']} fail"
        _cov = 'spans ±0.2' if covers else 'SHORT'
        print(f"{fn:<20}{_region:<26}{_selproc:<12}{_parse:<14}{_cov:<10}{'OK' if ok else 'CHECK'}")
        rows.append(dict(file=fn, target=wave_A, region_lo=lo, region_hi=hi,
                         n_selected=int(sel), n_processed=int(proc),
                         n_parsed=info['n_parsed'], n_fail=info['n_failures'],
                         not_truncated=not_trunc, covers_window=covers, verdict='OK' if ok else 'CHECK'))
    pd.DataFrame(rows).to_csv(_OUT / 'rya346_intake.csv', index=False)
    print(f"(wrote {(_OUT/'rya346_intake.csv').relative_to(_REPO)})")


def _present_in_synth(ctx, sp_key, w, e, wtol, etol=0.15):
    """Is a transition of this species+EP already in atomic_lines.tsv within wtol Å?"""
    ll = ctx['linelist']
    el = np.asarray([str(x) for x in ll['element']])
    wv = np.asarray(ll['wave_A'], dtype=float)
    ep = np.asarray(ll['lower_state_eV'], dtype=float)
    for i in np.where(np.abs(wv - w) <= wtol)[0]:
        try:
            if species_key(el[i]) == sp_key and abs(ep[i] - e) <= etol:
                return True, round(float(wv[i]), 3)
        except ValueError:
            continue
    return False, None


def _scratch_linelist(ctx, vald_rec):
    """In-memory scratch copy of the synth line list with ONE VALD transition appended
    (atomic_lines.tsv on disk untouched). Clones a same-species row for the format,
    overwrites the physical fields from the VALD record."""
    ll = ctx['linelist']
    sp_key = species_key(vald_rec['element'], vald_rec['ion'])
    el = np.asarray([str(x) for x in ll['element']])
    donor = None
    for i in range(len(ll)):
        try:
            if species_key(el[i]) == sp_key:
                donor = i; break
        except ValueError:
            continue
    if donor is None:
        donor = 0                      # fall back to any row for the dtype/format
    row = ll[donor:donor + 1].copy()
    w = float(vald_rec['wavelength'])
    row['wave_A'][0] = w
    if 'wave_nm' in ll.dtype.names: row['wave_nm'][0] = w / 10.0
    row['loggf'][0] = float(vald_rec['log_gf'])
    row['lower_state_eV'][0] = float(vald_rec['e_low_eV'])
    for nm, key in (('rad', 'damping_rad'), ('stark', 'damping_stark'), ('waals', 'damping_vdW')):
        if nm in ll.dtype.names:
            row[nm][0] = float(vald_rec[key])
    return np.concatenate([ll, row])


def diagnose(ctx):
    rows = []
    for wave_A, ew_obs, fn in TARGETS:
        print(f"\n{'='*72}\n  WINDOW {wave_A:.3f} Å  (pull {fn})\n{'='*72}")
        # 1. fit + residual peak
        a0, chi0, status, resid = _fit_line(ctx, wave_A, ew_obs)
        peak_w, peak_d = np.nan, 0.0
        if resid is not None:
            w, d = resid['wave_A'], resid['dflux']
            core = (np.abs(w - wave_A) <= 0.18) & (d < 0)
            if core.any():
                j = np.argmin(d[core]); peak_w = float(w[core][j]); peak_d = float(-d[core][j])
        print(f"  fit: A(FeII)={a0:.3f}  χ²ᵣ={chi0:.1f}  | residual peak {peak_w:.3f} Å "
              f"(Δcore {peak_w-wave_A:+.3f} Å, depth {peak_d:.3f})")

        # 2. candidates from the pull: absent from synth list, near the core
        recs, _ = vald_parse.parse_vald_long(str(_RAW / fn))
        cand = []
        for r in recs:
            if r['species'] == 'Fe 2' and abs(r['wavelength'] - wave_A) < 0.05:
                continue                                   # the target itself
            sp = species_key(r['element'], r['ion'])
            pres, w_in = _present_in_synth(ctx, sp, r['wavelength'], r['e_low_eV'],
                                           wtol=CORE_TOL_WIDE)
            dcore = abs(r['wavelength'] - wave_A)
            cand.append(dict(species=r['species'], wave=round(r['wavelength'], 3),
                             ep=r['e_low_eV'], loggf=r['log_gf'], cd=r['central_depth'],
                             in_synth=pres, synth_wave=w_in, dcore=round(dcore, 3)))
        cdf = pd.DataFrame(cand)
        # core-candidate set: absent from synth, within core tolerance (tight→wide)
        def core_absent(tol):
            if cdf.empty: return cdf
            return cdf[(~cdf.in_synth) & (cdf.dcore <= tol)]
        core_set = core_absent(CORE_TOL_TIGHT)
        used_tol = CORE_TOL_TIGHT
        if core_set.empty:
            core_set = core_absent(CORE_TOL_WIDE); used_tol = CORE_TOL_WIDE
        if cdf.empty:
            print("  pull transitions (excl. target): NONE — only the Fe II line in ±0.2 Å")
        else:
            print(f"  pull transitions (excl. target): {len(cdf)} "
                  f"({int((~cdf.in_synth).sum())} absent from synth list)")
            for _, c in cdf.sort_values('cd', ascending=False).head(6).iterrows():
                print(f"     {c.species:6} {c.wave:.3f}  EP={c.ep:.2f} loggf={c.loggf:+.2f} "
                      f"cd={c.cd:.4f}  Δcore={c.dcore:+.3f}  in_synth={c.in_synth}")
        print(f"  CORE-candidate set (absent & within ±{used_tol} Å of core): {len(core_set)}")

        # 3. decisive scratch test — strongest absent line anywhere in the pull, as the
        #    best-case control (if even this can't close it, no absorber explains the residual)
        a1 = chi1 = np.nan; tested = None
        absent_any = cdf[~cdf.in_synth] if not cdf.empty else cdf
        if not absent_any.empty:
            top = absent_any.sort_values('cd', ascending=False).iloc[0]
            tested = dict(species=top.species, wave=float(top.wave), ep=float(top.ep),
                          log_gf=float(top.loggf), e_low_eV=float(top.ep),
                          element=top.species.split()[0],
                          ion={'1':'I','2':'II','3':'III'}[top.species.split()[1]],
                          wavelength=float(top.wave),
                          damping_rad=0.0, damping_stark=0.0, damping_vdW=0.0)
            # carry the VALD damping for the tested record
            rec = next(r for r in recs if abs(r['wavelength']-top.wave) < 1e-3 and r['species']==top.species)
            tested.update(damping_rad=rec['damping_rad'], damping_stark=rec['damping_stark'],
                          damping_vdW=rec['damping_vdW'])
            scratch_ll = _scratch_linelist(ctx, tested)
            ctx2 = dict(ctx); ctx2['linelist'] = scratch_ll
            a1, chi1, _s, _r = _fit_line(ctx2, wave_A, ew_obs)
            label = "CONTROL" if (top.dcore > CORE_TOL_WIDE or top.cd < TRIVIAL_CD) else "candidate"
            print(f"  SCRATCH TEST [{label}]: +{top.species} {top.wave:.3f} (cd {top.cd:.4f}, "
                  f"Δcore {top.dcore:+.3f}) → A(FeII) {a0:.3f}→{a1:.3f}  χ²ᵣ {chi0:.1f}→{chi1:.1f}")
        else:
            print("  SCRATCH TEST: skipped — no line absent from the synth list to add "
                  "(VALD returns only already-modelled lines / the Fe II target)")

        chi_pass = np.isfinite(chi1) and chi1 < 10.0
        ab_back  = np.isfinite(a1) and abs(a1 - CLEAN_POOL_A) <= 0.10 and a1 < a0 - 1e-3
        absorber_found = bool(chi_pass and ab_back)
        verdict = 'absorber-found → BUILD' if absorber_found else 'no-absorber → Fe II atomic-data audit'
        print(f"  VERDICT: {verdict}")
        rows.append(dict(target=wave_A, a_before=a0, chi2r_before=round(chi0, 1),
                         resid_peak=round(peak_w, 3), resid_dcore=round(peak_w - wave_A, 3),
                         n_pull_lines=len(cdf), n_core_candidates=len(core_set),
                         tested_line=(f"{tested['species']} {tested['wave']:.3f}" if tested else None),
                         a_after=round(a1, 3) if np.isfinite(a1) else None,
                         chi2r_after=round(chi1, 1) if np.isfinite(chi1) else None,
                         absorber_found=absorber_found, verdict=verdict))
    out = pd.DataFrame(rows)
    out.to_csv(_OUT / 'rya346_diag.csv', index=False)
    nfound = int(out.absorber_found.sum())
    print(f"\n{'='*72}\n  TALLY: {nfound}/5 absorber-found, {5-nfound}/5 no-absorber")
    print(f"  Clean Fe II pool: 5 → {5 + nfound}  (RYA-341 balance "
          f"{'re-crowns' if nfound else 'UNCHANGED at −0.015'})")
    print(f"  wrote {(_OUT/'rya346_diag.csv').relative_to(_REPO)}")


def main():
    Path(_TMP).mkdir(parents=True, exist_ok=True)
    intake_table()
    ctx = _solar_context()
    diagnose(ctx)


if __name__ == '__main__':
    main()
