"""
scripts/rya347_fe2_atomic_data_audit.py
=======================================
RYA-347 — Fe II atomic-data audit on the 5 χ²ᵣ-failing solar Fe II lines
(5234.623, 5991.371, 6084.102, 6247.557, 6456.380). Lineage: RYA-343 (blend char)
→ RYA-344 (0/5 from existing extraction) → RYA-346 (no coincident absorber, deficit
opacity-insensitive) → here: is the residual the Fe II line's own gf and/or van der
Waals damping?

KEY PHYSICS the audit is built around:
  • gf↔abundance degeneracy. For an isolated line whose abundance is FLOATED by the
    flux fit, changing log gf by Δ is exactly a −Δ shift in the best-fit A(Fe II);
    the profile SHAPE (hence χ²ᵣ) is invariant. So a persistent high χ²ᵣ is a SHAPE
    problem (damping wings / profile), NOT gf. gf matters only for whether the fitted
    ABUNDANCE is biased. The fork therefore reads: χ²ᵣ ⇒ damping; abundance ⇒ gf.

SCRATCH ONLY — the canonical synth line list, linelist_solar.csv, EWs, spectra and
quarantine are NOT modified. Scratch edits are in-memory structured-array rows; the
flux fit (RYA-342 machinery) is re-run per condition. Sourced values only (NIST/
linelist_solar with citation); no guessed gf or vdW. χ²ᵣ<10 gate unchanged (RYA-342).
Matching via the RYA-345 species normalizer. Å throughout. Outputs to data/results/.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO))
sys.path.insert(0, str(_REPO / 'data' / 'linelists'))
sys.path.insert(0, str(_REPO / 'scripts'))

import pipeline.abundances_derive as ad
from pipeline.species import species_key
from rya344_recover_fe2_blends import _solar_context, _fit_line, CLEAN_POOL_A, _TMP

_OUT = _REPO / 'data' / 'results'; _OUT.mkdir(parents=True, exist_ok=True)
_SYNTH = Path(ad._SYNTH_LINELIST_FILE)
_LL_SOLAR = _REPO / 'data' / 'linelists' / 'linelist_solar.csv'
_NIST = _REPO / 'data' / 'linelists' / 'nist_reference.csv'
_NIST_XC = _REPO / 'data' / 'linelists' / 'nist_crosscheck.csv'

TARGETS = [  # (wave_A, ew_obs_mA)
    (5234.623, 120.7), (5991.371, 44.3), (6084.102, 28.5),
    (6247.557, 62.9), (6456.380, 85.4),
]
_FE_II = (26, 2)
ABUND_TOL = 0.10   # "unbiased" = within this of the RYA-341 clean-pool consensus (CLEAN_POOL_A)


# ── Step 0: locate + confirm atomic data ─────────────────────────────────────
def _synth_row_index(ctx, wave_A):
    ll = ctx['linelist']
    el = np.asarray([str(x) for x in ll['element']]); wv = np.asarray(ll['wave_A'], float)
    cand = [i for i in np.where(np.abs(wv - wave_A) <= 0.05)[0]
            if _safe_key(el[i]) == _FE_II]
    if not cand:
        return None
    return cand[int(np.argmin(np.abs(wv[cand] - wave_A)))]


def _safe_key(e):
    try: return species_key(e)
    except ValueError: return None


def step0(ctx):
    ll = ctx['linelist']
    sol = pd.read_csv(_LL_SOLAR, low_memory=False)
    sol = sol[(sol.element == 'Fe') & (sol.ion == 'II')]
    print(f"\n=== STEP 0 — synth list located: {_SYNTH}")
    print(f"{'line':>9} | {'synth gf':>8} {'LLsolar gf':>10} {'Δgf':>6} | "
          f"{'synth stark':>11} {'synth waals':>11} {'tt':>3} {'wsgl':>7} | LLsolar(stark,vdW) ref")
    rows = []
    for wave_A, _ew in TARGETS:
        i = _synth_row_index(ctx, wave_A)
        s = ll[i]
        m = sol.iloc[(sol.wavelength_air_A - wave_A).abs().argmin()]
        dgf = float(s['loggf']) - float(m.log_gf)
        rows.append(dict(line=wave_A, synth_gf=float(s['loggf']), llsolar_gf=float(m.log_gf),
                         # RYA-853 scope 4: step1's NIST match needs the EP, on BOTH sides.
                         # Without it the guard added there degenerates to "never match",
                         # which reads as "no graded source" and is a silent false absence.
                         ep_eV=float(m.excitation_potential_eV),
                         dgf=round(dgf, 3), synth_stark=float(s['stark']),
                         synth_waals=float(s['waals']), transition_type=str(s['spectrum_transition_type']),
                         waals_single_gamma=float(s['waals_single_gamma_format']),
                         llsolar_stark=float(m.damping_stark), llsolar_vdw=float(m.damping_vdW)))
        print(f"{wave_A:9.3f} | {float(s['loggf']):8.3f} {float(m.log_gf):10.3f} {dgf:6.3f} | "
              f"{float(s['stark']):11.2f} {float(s['waals']):11.3f} {str(s['spectrum_transition_type']):>3} "
              f"{float(s['waals_single_gamma_format']):7.3f} | ({float(m.damping_stark):.2f},{float(m.damping_vdW):.3f})")
    df = pd.DataFrame(rows)
    print(f"\n  synth Stark all identical? {df.synth_stark.nunique()==1} ({df.synth_stark.iloc[0]})  "
          f"| transition_type all 'AO' (ABO active)? {(df.transition_type=='AO').all()}")
    print(f"  LLsolar Stark all identical? {df.llsolar_stark.nunique()==1} ({df.llsolar_stark.iloc[0]})  "
          f"(catalog-default Stark — the RYA-347 'default damping' evidence; lives in the EW-path file)")
    print(f"  gf mismatch synth vs LLsolar? {(df.dgf.abs()>0.001).any()} "
          f"(max |Δgf| = {df.dgf.abs().max():.3f} dex) → DUPLICATED-VALUE DEFECT (fit consumes synth gf)")
    df.to_csv(_OUT / 'rya347_step0.csv', index=False)
    return df


# ── Step 1: gf cross-check ───────────────────────────────────────────────────
def step1(ctx, s0):
    nist = pd.read_csv(_NIST, comment='#'); nist = nist[(nist.element=='Fe')&(nist.ion=='II')]
    nxc = pd.read_csv(_NIST_XC, comment='#'); nxc = nxc[(nxc.element=='Fe')&(nxc.ion=='II')]
    print("\n=== STEP 1 — gf cross-check (synth=GES v6 vs sourced refs) ===")
    print(f"{'line':>9} {'synth(GESv6)':>12} {'LLsolar(VALD3)':>15} {'NIST ASD':>10} {'grade':>6}  source/flag")
    out = []
    # 🔴 RYA-853 scope 4: EP guard on BOTH sides, and ambiguity refused.
    # This matched on wavelength alone within +/-0.1 A and took .iloc[0] -- a 0.1 A window
    # at 6150 A spans several levels, and file order decided the winner. RYA-853 measured
    # the consequence directly: matched on wavelength alone, 6149.246 picks up an EP of
    # 13.436 eV instead of 3.889. The window is now +/-0.05 A AND the EP must agree.
    WAVE_TOL_A, EP_TOL_EV = 0.05, 0.05

    def _pick(tbl, w, ep):
        """The unique row matching on BOTH axes, or None. Never argmin: two rows the
        data cannot separate must not be separated by proximity."""
        if not np.isfinite(ep):
            return None
        m = tbl[((tbl.wavelength_air_A - w).abs() <= WAVE_TOL_A)
                & ((tbl.excitation_potential_eV - ep).abs() <= EP_TOL_EV)]
        return m.iloc[0] if len(m) == 1 else None

    for _, r in s0.iterrows():
        w = r.line
        ep = float(getattr(r, 'ep_eV', np.nan))
        hit = _pick(nist, w, ep)
        if hit is None:
            hit = _pick(nxc, w, ep)
        nist_gf = float(hit.log_gf) if hit is not None else np.nan
        grade = hit.nist_grade if hit is not None else '—'
        flag = ('NIST ASD v5.11' if np.isfinite(nist_gf)
                else 'NO EP-matched graded source → flag for manual NIST/M&B 2009 pull')
        print(f"{w:9.3f} {r.synth_gf:12.3f} {r.llsolar_gf:15.3f} "
              f"{(f'{nist_gf:.3f}' if np.isfinite(nist_gf) else '—'):>10} {str(grade):>6}  {flag}")
        out.append(dict(line=w, synth_gf=r.synth_gf, llsolar_gf=r.llsolar_gf,
                        nist_gf=nist_gf, nist_grade=grade,
                        best_gf=nist_gf if np.isfinite(nist_gf) else np.nan, source=flag))
    df = pd.DataFrame(out); df.to_csv(_OUT / 'rya347_step1_gf.csv', index=False)
    return df


# ── Step 3: scratch test ─────────────────────────────────────────────────────
# vdW NOTE: the synth list already uses ABO (transition_type 'AO', turbospectrum_fdamp
# = packed σ.α). The classical single-gamma value the ticket cited lives only in
# linelist_solar (the EW path). Empirically, switching the TS path to classical Unsold
# is a no-op here (Unsold ≈ the ABO σ for these levels). So condition B is the decisive
# damping-FIXABILITY test: a sensitivity scan over physical ABO σ (0.5–4× baseline,
# abundance refit). If no σ in that range reaches χ²ᵣ<10, damping cannot be the fault
# (the residual is on the core; damping reshapes the wings). This DIAGNOSES fixability;
# it does NOT adopt an unsourced value — a real correction still needs a sourced Barklem σ.
_SIGMA_MULTS = [0.5, 2.0, 4.0]


def _scratch(ctx, wave_A, loggf=None, sigma_mult=None):
    """In-memory scratch synth list with one Fe II row modified (ABO σ scaled by
    sigma_mult; gf overwritten if given). atomic_lines.tsv on disk untouched."""
    ll = ctx['linelist'].copy()
    i = _synth_row_index(ctx, wave_A)
    if loggf is not None:
        ll['loggf'][i] = float(loggf)
    if sigma_mult is not None:
        base = ll['turbospectrum_fdamp'][i]
        sig = int(base); alpha = base % 1
        new = max(20, int(sig * sigma_mult)) + alpha
        ll['turbospectrum_fdamp'][i] = new
        if 'waals' in ll.dtype.names: ll['waals'][i] = new
        if 'spectrum_transition_type' in ll.dtype.names: ll['spectrum_transition_type'][i] = 'AO'
    return ll


def step3(ctx, s1):
    print("\n=== STEP 3 — scratch flux-fit: per line × {baseline, A:gf, B:vdW-scan, C:both} ===")
    print("  A = sourced gf (NIST where available); B = best of an ABO σ scan (0.5–4×);")
    print("  C = both. χ²ᵣ unchanged under A confirms gf↔abundance degeneracy.")
    print(f"{'line':>9} {'cond':>9} {'A(FeII)':>8} {'χ²ᵣ':>7}  note")
    rows = []
    for wave_A, ew in TARGETS:
        bestgf = s1[s1.line == wave_A].best_gf.iloc[0]
        gf_alt = bestgf if np.isfinite(bestgf) else None
        # baseline
        a0, chi0, st0, _ = _fit_line(ctx, wave_A, ew)
        print(f"{wave_A:9.3f} {'baseline':>9} {a0:8.3f} {chi0:7.1f}  {st0}")
        rows.append(dict(line=wave_A, cond='baseline', a_synth=round(a0,3), red_chi2=round(chi0,1), note=st0))
        # A: gf
        if gf_alt is not None:
            ctxA = dict(ctx); ctxA['linelist'] = _scratch(ctx, wave_A, loggf=gf_alt)
            aA, chiA, stA, _ = _fit_line(ctxA, wave_A, ew)
            print(f"{wave_A:9.3f} {'A:gf':>9} {aA:8.3f} {chiA:7.1f}  {stA}")
            rows.append(dict(line=wave_A, cond='A:gf', a_synth=round(aA,3), red_chi2=round(chiA,1), note=stA))
        else:
            print(f"{wave_A:9.3f} {'A:gf':>9} {'—':>8} {'—':>7}  no sourced gf alt → manual NIST/M&B pull")
            rows.append(dict(line=wave_A, cond='A:gf', a_synth=np.nan, red_chi2=np.nan, note='no sourced gf'))
        # B: ABO σ scan → keep the best (min χ²ᵣ)
        best = None
        for m in _SIGMA_MULTS:
            ctxB = dict(ctx); ctxB['linelist'] = _scratch(ctx, wave_A, sigma_mult=m)
            aB, chiB, _, _ = _fit_line(ctxB, wave_A, ew)
            if best is None or chiB < best[1]: best = (aB, chiB, m)
        aB, chiB, mB = best
        print(f"{wave_A:9.3f} {'B:vdW':>9} {aB:8.3f} {chiB:7.1f}  best σ×{mB} of scan {_SIGMA_MULTS}")
        rows.append(dict(line=wave_A, cond='B:vdW', a_synth=round(aB,3), red_chi2=round(chiB,1),
                         note=f'best σ×{mB}'))
        # C: both (best σ + gf)
        if gf_alt is not None:
            ctxC = dict(ctx); ctxC['linelist'] = _scratch(ctx, wave_A, loggf=gf_alt, sigma_mult=mB)
            aC, chiC, stC, _ = _fit_line(ctxC, wave_A, ew)
            print(f"{wave_A:9.3f} {'C:both':>9} {aC:8.3f} {chiC:7.1f}  gf+σ×{mB}")
            rows.append(dict(line=wave_A, cond='C:both', a_synth=round(aC,3), red_chi2=round(chiC,1), note=f'gf+σ×{mB}'))
        else:
            rows.append(dict(line=wave_A, cond='C:both', a_synth=np.nan, red_chi2=np.nan, note='no sourced gf'))
    df = pd.DataFrame(rows); df.to_csv(_OUT / 'rya347_step3_scratch.csv', index=False)
    return df


# ── Step 4: fork ─────────────────────────────────────────────────────────────
def step4(s3, s1):
    print("\n=== STEP 4 — per-line fork ===")
    piv = s3.pivot_table(index='line', columns='cond', values=['red_chi2', 'a_synth'], aggfunc='first')
    verdicts = []
    for wave_A in [t[0] for t in TARGETS]:
        def gv(metric, cond):
            try: return piv.loc[wave_A, (metric, cond)]
            except KeyError: return np.nan
        chi_b = gv('red_chi2', 'baseline'); a_b = gv('a_synth', 'baseline')
        chi_A = gv('red_chi2', 'A:gf');    a_A = gv('a_synth', 'A:gf')
        chi_B = gv('red_chi2', 'B:vdW');   a_B = gv('a_synth', 'B:vdW')
        chi_C = gv('red_chi2', 'C:both');  a_C = gv('a_synth', 'C:both')
        def ok(chi, a): return np.isfinite(chi) and chi < 10.0 and abs(a - CLEAN_POOL_A) <= ABUND_TOL
        if ok(chi_A, a_A):      v = 'gf-fault, fixable → KEEP (differential-safe)'
        elif ok(chi_B, a_B):    v = 'damping-fault, fixable → KEEP (sourced ABO needed)'
        elif ok(chi_C, a_C):    v = 'both-needed, fixable → KEEP'
        else:
            nist_ok = np.isfinite(s1[s1.line==wave_A].nist_gf.iloc[0])
            v = ('unfixable with sourced values → RETIRE from solar clean pool; '
                 'distrust on metal-rich (55 Cnc)'
                 + ('' if nist_ok else ' [also needs manual NIST/M&B gf + Barklem ABO pulls to be definitive]'))
        verdicts.append(dict(line=wave_A, chi2_base=round(chi_b,1), chi2_A=round(chi_A,1) if np.isfinite(chi_A) else None,
                             chi2_B=round(chi_B,1) if np.isfinite(chi_B) else None,
                             chi2_C=round(chi_C,1) if np.isfinite(chi_C) else None, verdict=v))
        print(f"  {wave_A:.3f}: χ²ᵣ base={chi_b:.1f} A={chi_A if not np.isfinite(chi_A) else round(chi_A,1)} "
              f"B={chi_B if not np.isfinite(chi_B) else round(chi_B,1)} "
              f"C={chi_C if not np.isfinite(chi_C) else round(chi_C,1)} → {v}")
    df = pd.DataFrame(verdicts); df.to_csv(_OUT / 'rya347_step4_fork.csv', index=False)
    nkeep = sum('KEEP' in v['verdict'] for v in verdicts)
    print(f"\n  HEADLINE: {nkeep}/5 fixable→keep, {5-nkeep}/5 retire/needs-manual")
    print(f"  Projected clean Fe II pool: 5 → {5 + nkeep}  "
          f"(RYA-341 {'re-crowns via follow-up BUILD' if nkeep else 'UNCHANGED −0.015'})")
    return df


def main():
    Path(_TMP).mkdir(parents=True, exist_ok=True)
    ctx = _solar_context()
    s0 = step0(ctx)
    s1 = step1(ctx, s0)
    s3 = step3(ctx, s1)
    step4(s3, s1)


if __name__ == '__main__':
    main()
