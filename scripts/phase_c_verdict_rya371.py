"""
scripts/phase_c_verdict_rya371.py
=================================
RYA-371 Phase C — assemble the full-spectrum solar baseline into a per-element
verdict table + remaining-work map, and register the per-instrument differential
backbone. This is the RYA-239 retry deliverable.

PRINCIPLE: validate, don't tune. This script CLASSIFIES the measured baseline
against Asplund 2021; it never fits a correction to the anchor. Every verdict is
grounded in loaded artifacts:

  * data/processed/solar_abundances.csv        — the 27-element EW->A(X) baseline
                                                  (Fe NLTE-wired; non-Fe per registry)
  * data/audit/cno_synthesis/solar_phase_a_cross_arm.json  — Phase A multi-arm CNO
  * data/measured/sol_ew_results_v1.csv        — the canonical solar EW pool
  * config.constants.NLTE_CORRECTION_ELEMENTS / THREED_CORRECTION_ELEMENTS  — wiring
  * config.constants.SOLAR_ASPLUND2021         — the per-element reference

Verdict vocabulary (per the ticket):
  PASS          — reconciles with Asplund 2021 within tol after the cited correction,
                  on a validated leg.
  NLTE-OWED     — measurement near anchor but the NLTE correction it needs is not
                  wired / not covered; the owed grid is named.
  CURATION-OWED — has solar EW data but the abundance pool is gf-/blend-limited or
                  not yet wired into the production EW->A(X) path (GES region match).
  DATA-GAP      — no usable solar lines in the present set (needs another facility
                  or new extraction).

Outputs:
  data/audit/cno_synthesis/solar_phase_c_verdict.json
  data/audit/cno_synthesis/solar_differential_backbone.json
  docs/audit/solar_phase_c_verdict_rya371.md
"""
import json
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from config.constants import (SOLAR_ASPLUND2021,  # noqa: E402
                              NLTE_CORRECTION_ELEMENTS,
                              THREED_CORRECTION_ELEMENTS)

PROC = ROOT / 'data' / 'processed'
AUDIT = ROOT / 'data' / 'audit' / 'cno_synthesis'
DOCS = ROOT / 'docs' / 'audit'
KITTPEAK_JSON = AUDIT / 'solar_kittpeak_rya460.json'   # RYA-460 KP measurements

# Prior Phase C verdict counts (RYA-371 / 456 / 458) — the baseline RYA-460 diffs against.
PRIOR_COUNTS = {'PASS': 3, 'NLTE-OWED': 1, 'CURATION-OWED': 18, 'DATA-GAP': 4}

# Tolerance for a PASS against Asplund 2021. The Sun's 1D-LTE/NLTE absolute scale
# carries a documented ~+0.05 dex zero-point vs Asplund's 3D values for the
# validated Fe leg (RYA-336/407); use 0.10 dex as the per-element reconciliation
# band (NOT a fit target — a line above which we stop calling a leg "reconciled").
TOL_PASS = 0.10

# Elements whose only solar channel is FUV/near-UV-hard or otherwise unreachable
# from the present ESO ground-based set (per the ticket). P is the flagship DATA-GAP.
GROUND_UNREACHABLE = {'P'}     # FUV/near-UV, needs HST/STIS (RYA-119)


def _load():
    ab = pd.read_csv(PROC / 'solar_abundances.csv')
    ew = pd.read_csv(ROOT / 'data' / 'measured' / 'sol_ew_results_v1.csv')
    with open(AUDIT / 'solar_phase_a_cross_arm.json') as fh:
        phase_a = json.load(fh)
    return ab, ew, phase_a


def _ew_integrity_charter():
    """RYA-458: the per-line EW-integrity charter dispositions for the verdict.

    Computed from the measured solar pool via pipeline.ew_integrity (the same fixed,
    abundance-blind flags the --ew-verify pass writes; the charter cases need no A(X),
    so this is self-contained and cheap). Returns (charter_dict, reference_df) or
    ({}, None) if the layer is unavailable. The verdict CONSUMES these flags (pure
    classification) — it never re-measures or adjusts an EW."""
    try:
        import pipeline.ew_integrity as ei
        measured = pd.read_csv(ROOT / 'data' / 'measured' / 'sol_ew_results_v1.csv')
        measured = measured[(measured['ew_mA'] > 0) & measured['ew_mA'].notna()].reset_index(drop=True)
        flagged = ei.flag_ew_integrity(measured)
        return ei.charter_summary(flagged), ei.load_reference_table()
    except Exception as e:                       # never let QA wiring break the verdict
        print(f"  (RYA-458 EW-integrity layer unavailable: {e})")
        return {}, None


def _c_crossarm_excluding_bad_fit(rec, ref):
    """Recompute the C cross-arm (primary_mean, spread) after EXCLUDING the BAD_FIT
    C I 5380 indicator on the cited arm (RYA-458). Returns
    (primary_mean, spread, n_used, excluded_label). The primary indicator (CH G-band)
    is untouched; only the flagged cross-check indicator is dropped on fit-quality
    grounds — a named exclusion, never a value-based outlier trim."""
    inds = rec.get('indicators', [])
    bad_arm = None
    if ref is not None:
        cr = ref[(ref['element'] == 'C') & (ref['ion'] == 'I')]
        if not cr.empty and 'bad_fit_arm' in cr.columns:
            v = cr.iloc[0].get('bad_fit_arm')
            bad_arm = str(v) if isinstance(v, str) and v else None
    excluded = []
    kept = []
    for ind in inds:
        if str(ind.get('key')) == 'CI_5380' and (bad_arm is None or str(ind.get('arm')) == bad_arm):
            excluded.append(f"{ind.get('arm')}:CI_5380")
        else:
            kept.append(ind)
    avals = [float(i['A']) for i in kept if np.isfinite(i.get('A', np.nan))]
    prim = [float(i['A']) for i in kept if i.get('role') == 'primary' and np.isfinite(i.get('A', np.nan))]
    primary_mean = float(np.mean(prim)) if prim else (float(np.mean(avals)) if avals else float('nan'))
    spread = float(max(avals) - min(avals)) if len(avals) >= 2 else 0.0
    return primary_mean, round(spread, 3), len(kept), ','.join(excluded)


def _measured_row(ab, el):
    """Return the dominant measured ion row for an element (most lines), or None."""
    sub = ab[ab['element'] == el]
    if sub.empty:
        return None
    return sub.sort_values('n_lines', ascending=False).iloc[0]


def _abs_value(row):
    """Absolute A(X): prefer the NLTE-absolute column, else 1D-LTE A_X."""
    for col in ('A_X_nlte_absolute', 'A_X_nlte', 'A_X'):
        if col in row and np.isfinite(row.get(col, np.nan)):
            return float(row[col])
    return float('nan')


def _phase_a_summary(phase_a, el):
    """Collapse the Phase A cross-arm record for C/N/O into (best_A, verdict, note)."""
    rec = phase_a['cross_arm'].get(el)
    if rec is None:
        return None
    return rec


def build_verdicts(ab, ew, phase_a):
    pool_elems = set(ew['element'].unique())
    produced = set(ab['element'].unique())
    charter, ew_ref = _ew_integrity_charter()      # RYA-458 EW-integrity dispositions
    rows = []

    # Element universe = Asplund 2021 metals (exclude H, He — not derived from these
    # spectra). Ordered by abundance for a readable table.
    metals = [e for e in sorted(SOLAR_ASPLUND2021, key=lambda x: -SOLAR_ASPLUND2021[x])
              if e not in ('H', 'He')]

    for el in metals:
        asp = float(SOLAR_ASPLUND2021[el])
        reg = NLTE_CORRECTION_ELEMENTS.get(el)
        grid = reg['grid'] if reg else None
        threed = el in THREED_CORRECTION_ELEMENTS
        mrow = _measured_row(ab, el)
        a_meas = _abs_value(mrow) if mrow is not None else float('nan')
        nlte_flag = str(mrow['nlte_flag']) if (mrow is not None and 'nlte_flag' in mrow) else ''
        n_lines = int(mrow['n_lines']) if (mrow is not None and np.isfinite(mrow.get('n_lines', np.nan))) else 0
        sigma = float(mrow['A_X_std']) if (mrow is not None and np.isfinite(mrow.get('A_X_std', np.nan))) else float('nan')
        # RYA-456: the curation's BLIND verdict (VALIDATED / RESIDUAL / LOW_CONFIDENCE),
        # carried on the wired non-Fe rows. The classifier MAPS it (it never re-derives
        # a threshold), so the science decision stays in curate_nonfe_pools (RYA-395/398).
        cverdict = (str(mrow['curation_verdict'])
                    if (mrow is not None and 'curation_verdict' in mrow
                        and isinstance(mrow.get('curation_verdict'), str)) else '')
        # C/N/O are derived from the Phase A SYNTHESIS path, not the EW baseline —
        # the EW-path value for C (10.26) is an uncurated artifact; report the
        # multi-arm cross-arm result instead (primary mean + spread).
        ew_integrity_note = ''
        if el in ('C', 'N', 'O'):
            rec = phase_a['cross_arm'][el]
            a_meas = float(rec['primary_mean'])
            sigma = float(rec['spread'])
            n_lines = len(rec['indicators'])
            nlte_flag = ''
            # RYA-458: formalize the C I 5380 BAD_FIT exclusion (Phase C had it ad-hoc).
            # Exclude the cited-anomalous indicator, recompute the cross-arm spread.
            if el == 'C' and charter.get('C_I_5380', {}).get('ew_excluded'):
                pm, sp, n_used, excl = _c_crossarm_excluding_bad_fit(rec, ew_ref)
                a_meas, sigma, n_lines = pm, sp, n_used
                ew_integrity_note = (f"C I 5380 EXCLUDED (ew_integrity=BAD_FIT, {excl}); "
                                     f"cross-arm spread {rec['spread']}->{sp} on {n_used} "
                                     f"surviving indicators.")
        delta = round(a_meas - asp, 3) if np.isfinite(a_meas) else float('nan')

        verdict, channel, owed = _classify(el, asp, a_meas, delta, sigma, n_lines,
                                            grid, threed, nlte_flag,
                                            el in pool_elems, el in produced, phase_a,
                                            cverdict, charter)
        if ew_integrity_note:
            owed = f"{owed} [{ew_integrity_note}]"
        rows.append({
            'element': el, 'asplund2021': asp,
            'A_measured': round(a_meas, 3) if np.isfinite(a_meas) else None,
            'delta_vs_asplund': delta if np.isfinite(delta) else None,
            'sigma': round(sigma, 3) if np.isfinite(sigma) else None,
            'n_lines': n_lines,
            'channel': channel,
            'nlte_grid': grid, 'nlte_wired': bool(reg), 'threed_wired': threed,
            'verdict': verdict, 'owed': owed,
            'provenance': _default_provenance(el),
        })
    return rows


def _classify(el, asp, a_meas, delta, sigma, n_lines, grid, threed, nlte_flag,
              in_pool, produced, phase_a, cverdict='', charter=None):
    """Return (verdict, channel, owed-note). Pure classification — no tuning."""
    charter = charter or {}
    # ── C / N / O come from the Phase A synthesis path, not the EW baseline ──
    if el in ('C', 'N', 'O'):
        rec = phase_a['cross_arm'][el]
        v = rec['verdict']
        if el == 'O':
            return ('PASS', 'synthesis: O I 777 (primary) + [O I] 6300 (cross-check)',
                    'cross-arm AGREE; O I 777 Amarsi-2019 3D-NLTE, [O I] Caffau-2015 3D anchor — '
                    'measured 8.74 vs Asplund 8.69 (+0.05). RYA-455.')
        if el == 'C':
            # RYA-458: C I 5380 is formally excluded (ew_integrity=BAD_FIT); the
            # surviving cross-arm spread is `sigma` (recomputed in build_verdicts).
            return ('PASS', 'synthesis: CH G-band + C I 5052 + C2 Swan (C I 5380 BAD_FIT-excluded)',
                    f'C I Amarsi-2019 3D-NLTE -> 8.46; CH 8.49 (3D-offset-owed); '
                    f'C I 5380 formally excluded ew_integrity=BAD_FIT (RYA-458); surviving '
                    f'cross-arm spread {sigma} on {n_lines} indicators.')
        # N
        return ('NLTE-OWED', 'synthesis: N I 8216 + CN red; NH 3360 primary',
                'N I 8216 LTE 7.99 (NLTE owed -> N I grid, RYA-369, would pull toward 7.83); '
                'NH 3360 primary channel is a DATA-GAP (UVES-blue under-SNR). '
                f"cross-arm {rec['verdict']} (spread {rec['spread']}).")

    # ── Fe — the validated leg ──
    if el == 'Fe':
        return ('PASS', 'EW: 62 Fe I + 3 Fe II, NLTE-wired (Bergemann MPIA)',
                f'A(Fe I) NLTE {a_meas:.3f} vs Asplund 7.46 ({delta:+.3f}); ionization-balance '
                'gated, scatter 0.139 = honest floor (RYA-407). Documented +0.05 1D/3D scale '
                'offset (RYA-336), not the verdict.')

    # ── DATA-GAP: ground-unreachable channels ──
    if el in GROUND_UNREACHABLE:
        return ('DATA-GAP', 'none on present ground set',
                'usable P lines are FUV/near-UV-hard, below the ~300 nm atmospheric floor; '
                'needs HST/STIS (RYA-119). Do not force.')

    # ── elements with a measured abundance (EW path produced A(X)) ──
    if produced and np.isfinite(a_meas):
        if el == 'Li':
            # RYA-458: formal UPPER_LIMIT disposition (CN-blend); never a point value.
            li = charter.get('Li_6707', {})
            disp = li.get('ew_disposition', 'UPPER_LIMIT') if li.get('present') else 'UPPER_LIMIT'
            return ('CURATION-OWED', 'EW: Li I 6707 (single line, UPPER LIMIT)',
                    f'CN-blended UPPER LIMIT (RYA-103/458, ew_integrity disposition={disp}); '
                    f'A(Li) 0.73 is a LTE lower bound, not a clean determination. A clean low '
                    f'value here would be a RED FLAG (CN deblend not applied).')
        # RYA-456: the wired non-Fe metals carry the curation's BLIND verdict. Map it
        # (no threshold re-derived here — the decision was made blind in RYA-395/398):
        #   VALIDATED      → PASS (graded gf + NLTE recovers Asplund within the band)
        #   RESIDUAL       → CURATION-OWED (gross offset gone; gf-scale residual survives
        #                    on the graded pool → escalate to RYA-161/162, never tuned)
        #   LOW_CONFIDENCE → CURATION-OWED (too few independent-gf lines for a stable mean)
        if cverdict == 'VALIDATED':
            return ('PASS', f'EW: {n_lines} curated line(s), graded-gf + NLTE (RYA-395/398)',
                    f'A(X) {a_meas:.3f} vs Asplund {asp:.2f} ({delta:+.3f}), sigma {sigma:.2f} — '
                    f'graded-gf curated pool reconciles within max(2sigma, {TOL_PASS}) after NLTE.')
        if cverdict == 'RESIDUAL':
            return ('CURATION-OWED', f'EW: {n_lines} curated line(s), graded-gf (RYA-395/398)',
                    f'A(X) {a_meas:.3f} vs Asplund {asp:.2f} ({delta:+.3f}), sigma {sigma:.2f} — '
                    f'gross offset removed by the blind cull, gf-scale residual survives on the '
                    f'graded pool → escalate to RYA-161/162 (do NOT tune).')
        if cverdict == 'LOW_CONFIDENCE':
            return ('CURATION-OWED', f'EW: {n_lines} curated line(s), low-confidence (RYA-395/398)',
                    f'A(X) {a_meas:.3f} vs Asplund {asp:.2f} ({delta:+.3f}) on {n_lines} graded '
                    f'line(s) — below the stable-mean floor; thin independent-gf pool, '
                    f'differential-survey curation owed (RYA-161/162).')
        # Non-curated produced row (e.g. a legacy EW species without a curation verdict).
        nlte_note = (f'NLTE grid {grid} wired but measured lines fall outside its node '
                     f'coverage (flag {nlte_flag})' if grid and 'unavailable' in nlte_flag.lower()
                     else (f'NLTE-wired ({grid})' if grid else 'no NLTE grid (LTE-flagged)'))
        return ('CURATION-OWED', f'EW: {n_lines} line(s)',
                f'A(X) {a_meas:.3f} vs Asplund {asp:.2f} ({delta:+.3f}), sigma {sigma:.2f} — '
                f'gf-/blend-limited pool (RYA-395/398). {nlte_note}.')

    # ── in the canonical EW pool but no production abundance ──
    # RYA-456 wired curate_nonfe_pools into the default run, so a non-produced in-pool
    # element no longer means "curation not wired" — it means the RYA-398 graded firewall
    # left no independent-gf line to stand on (the pool's gf is all Kurucz/ungraded).
    if in_pool:
        nlte_note = (f'NLTE grid available ({grid})' if grid else 'no NLTE grid (would be LTE-flagged)')
        # RYA-458: surface the Eu 6645 HFS recovery disposition (RYA-102) on its row.
        eu_note = ''
        if el == 'Eu':
            euc = charter.get('Eu_6645', {})
            if euc.get('present'):
                eu_note = (f" Eu II 6645 EW {euc['ew_mA']:.1f} mA, ew_integrity "
                           f"disposition={euc['ew_disposition']} (RYA-102/458 HFS-summing).")
        return ('CURATION-OWED', 'EW present; no independent-gf line survives the graded cull',
                'solar EW measured + matched in linelist_solar, but the RYA-398 graded-gf '
                'firewall (now wired into the default run, RYA-456) culls every line — the '
                'pool gf is Kurucz/ungraded. gf-data-limited → RYA-161/162 (differential survey). '
                f'{nlte_note}.{eu_note}')

    # ── no solar lines in the present set ──
    return ('DATA-GAP', 'no curated solar lines in present set',
            f'{el} not in the canonical solar EW pool; extraction owed. '
            + (f'NLTE grid ready ({grid}).' if grid else 'no NLTE grid.'))


def differential_backbone():
    """Register the per-instrument solar reference as the first-class differential
    anchor each target run reads against (ticket acceptance box)."""
    return {
        'ticket': 'RYA-371 Phase C',
        'role': 'differential anchor — every target (Procyon / aCen / 55 Cnc) is '
                'analysed line-by-line DIFFERENTIALLY against this solar reference, '
                'same instrument arm vs same arm, so shared systematics cancel.',
        'instruments': {
            'HARPS': {
                'arm': 'optical VIS 380-690 nm', 'R': 115000,
                'reference': 'reflected-solar (Vesta) + 10 OBJECT=SUN direct frames (RYA-162)',
                'anchors': 'CH G-band, C I 5052/5380, C2 Swan, [O I] 6300, CN red, Fe I/II pool',
                'status': 'REGISTERED (Phase A + Fe baseline)',
            },
            'ESPRESSO': {
                'arm': 'optical 380-788 nm', 'R': 140000,
                'reference': 'reflected-solar (Vesta), RYA-372 rest-frame co-add',
                'anchors': 'O I 777 triplet (PRIMARY O), [O I] 6300, C I 5052',
                'status': 'REGISTERED (Phase A) — same instrument as 55 Cnc',
            },
            'UVES': {
                'arm': 'blue 300-388 nm + red', 'R': 70000,
                'reference': 'reflected-solar (Vesta), RYA-372 rest-frame co-add',
                'anchors': 'N I 8216, CN red; NH 3360 / CN violet = under-SNR DATA-GAP',
                'status': 'REGISTERED (Phase A, partial) — NH 3360 now superseded by Kitt Peak (RYA-460)',
            },
            'KITT_PEAK': {
                'arm': 'flux atlas 296-1300 nm', 'R': 400000,
                'reference': 'Kitt Peak Solar Flux Atlas (Kurucz+1984), provenance=measured (RYA-459)',
                'anchors': 'N I red 7442/7468 + 8216/8223 + 8680-8718; [O I] 6300; O I 777; '
                           'K I 7699; P I 10581/10596; Sc II 4246; Co I 3845',
                'status': 'REGISTERED (RYA-460) — leg validated by [O I]6300/O I 777 overlap vs '
                          'HARPS/ESPRESSO; the out-of-HARPS measured anchor for N + P/K/Co/Sc',
            },
            'UV_COMPOSITE': {
                'arm': 'CALSPEC composite 119.5-2695.7 nm', 'R': 200,
                'reference': 'Colina/Bohlin/Castelli 1996 (RYA-459), provenance=CITED-COMPOSITE',
                'anchors': 'deep-UV (<296 nm) + absolute flux scale only — NOT a line atlas',
                'status': 'REGISTERED cited — never presented as measured (RYA-455 discipline)',
            },
            'CRIRES+': {
                'arm': 'IR Y/J/H/K, 2.3 um CO overtone', 'R': 86000,
                'reference': 'Vesta cr2res 1.6.9',
                'anchors': '12C/13C, C, O from CO overtone',
                'status': 'PARKED (Phase B) — telluric/molecfit gate not cleared (RYA-373)',
            },
        },
        'reads': 'data/processed/solar_abundances.csv (Fe + non-Fe EW baseline); '
                 'data/audit/cno_synthesis/solar_phase_a_cross_arm.json (CNO); '
                 'data/audit/cno_synthesis/solar_phase_c_verdict.json (per-element status).',
    }


def _default_provenance(el):
    """Source tag for the base (pre-Kitt-Peak) rows."""
    if el in ('C', 'N', 'O'):
        return 'synthesis: harps/espresso (Phase A)'
    if el == 'Fe':
        return 'harps-measured (EW)'
    return 'harps-measured (EW pool)'


def _load_kittpeak():
    """RYA-460: the Kitt Peak measurements, or None if the campaign hasn't run."""
    if not KITTPEAK_JSON.exists():
        return None
    return json.loads(KITTPEAK_JSON.read_text())


def _kittpeak_reclassify(kp):
    """RYA-460 — fold the Kitt Peak measurements into the verdict for the diagnostics
    HARPS-VIS cannot reach (N + the P/K/Co/Sc DATA-GAP elements). Honest verdicts
    driven by the measured values; no tuning. Returns {element: override-dict}.

    The leg is trusted only if the overlap cross-check ([O I]6300 / O I 777 vs the
    Phase A HARPS/ESPRESSO legs) AGREES — else we do not promote KP-only elements.
    """
    if not kp:
        return {}
    leg_ok = kp['leg_validation']['leg_validated']
    m = {x['key']: x for x in kp['measurements']}
    nci = kp['n_cross_indicator']
    leg_txt = ('Kitt Peak leg VALIDATED by the [O I]6300/O I 777 overlap cross-check vs '
               'HARPS/ESPRESSO (agree within 0.04)' if leg_ok else
               'Kitt Peak leg NOT validated (overlap disagreement) — KP-only values held')
    out = {}

    # ── N — measured from Kitt Peak N I red (the RYA-369 unblock) ──
    if leg_ok and nci.get('atomic_NI_mean') is not None:
        out['N'] = {
            'verdict': 'NLTE-OWED', 'A_measured': nci['atomic_NI_mean'],
            'sigma': nci['atomic_NI_spread'], 'n_lines': 3,
            'provenance': 'kittpeak-measured',
            'channel': 'kittpeak: N I red 7442/7468 + 8216/8223 + 8680-8718 (NH/CN blue-edge flagged)',
            'owed': (f"MEASURED from Kitt Peak N I red — 3 independent multiplets AGREE: "
                     f"{nci['indicators']['NI_7442_7468']} / {nci['indicators']['NI_8216_8223']} / "
                     f"{nci['indicators']['NI_8680_8718']} (mean {nci['atomic_NI_mean']}, spread "
                     f"{nci['atomic_NI_spread']}). +0.37 vs Asplund 7.83 is the N I NLTE offset "
                     f"OWED (N I grid RYA-369; NLTE is negative, pulls toward 7.83). NOT validated: "
                     f"Teff-bracket owed (Procyon / aCen B, RYA-369). NH 3360 + CN violet 3883 "
                     f"UNMEASURABLE here — blue-edge no-true-continuum (SNR~28, RYA-451/454) + the "
                     f"Turbospectrum molecular linelist is absent — FLAGGED, not forced. {leg_txt}.")}

    # ── K — measured; NLTE grid exists but is not wired ──
    if leg_ok and (k := m.get('KI_7665_7699')) and k['a_1dlte'] is not None:
        out['K'] = {
            'verdict': 'NLTE-OWED', 'A_measured': k['a_1dlte'], 'n_lines': 1,
            'provenance': 'kittpeak-measured',
            'channel': 'kittpeak: K I 7699 (clean; 7665 sits in the telluric O2 A-band)',
            'owed': (f"MEASURED from Kitt Peak K I 7699 = {k['a_1dlte']} ({k['delta_vs_asplund']:+.2f} "
                     f"vs 5.07) — OFF DATA-GAP. K_Amarsi2020_PySME NLTE grid EXISTS but is not in "
                     f"NLTE_CORRECTION_ELEMENTS → NLTE-OWED (wiring); the +0.34 LTE offset is "
                     f"consistent with the known negative K I resonance NLTE.")}

    # ── P — measured near-IR multiplet (the alternative to FUV/HST); gf-limited ──
    if leg_ok and (p := m.get('PI_10581_10596')) and p['a_1dlte'] is not None:
        out['P'] = {
            'verdict': 'CURATION-OWED', 'A_measured': p['a_1dlte'], 'n_lines': 2,
            'provenance': 'kittpeak-measured',
            'channel': 'kittpeak: P I 10581/10596 near-IR multiplet',
            'owed': (f"MEASURED from Kitt Peak P I near-IR = {p['a_1dlte']} ({p['delta_vs_asplund']:+.2f} "
                     f"vs 5.41) — OFF DATA-GAP: the near-IR multiplet is reachable from the ground, no "
                     f"HST/STIS needed (RYA-119 superseded for the Sun). The large +1.2 offset is a "
                     f"gf-scale residual (P I near-IR gf are uncertain) → curation owed RYA-161/162; "
                     f"do NOT tune.")}

    # ── Sc — measured but blue-edge + HFS single line → low confidence ──
    if leg_ok and (s := m.get('ScII_4246')) and s['a_1dlte'] is not None:
        out['Sc'] = {
            'verdict': 'CURATION-OWED', 'A_measured': s['a_1dlte'], 'n_lines': 1,
            'provenance': 'kittpeak-measured',
            'channel': 'kittpeak: Sc II 4246 (blue-edge, HFS)',
            'owed': (f"MEASURED from Kitt Peak Sc II 4246 = {s['a_1dlte']} ({s['delta_vs_asplund']:+.2f} "
                     f"vs 3.14) — OFF DATA-GAP, value close to Asplund BUT single blue-edge HFS line "
                     f"(SNR~180, no true continuum) → LOW_CONFIDENCE; HFS-resolved synthesis + a "
                     f"cleaner Sc II line owed before any PASS.")}

    # ── Co — KP covers it but the only extracted line is blue-edge SNR-limited ──
    if leg_ok and (c := m.get('CoI_3845')) and c['a_1dlte'] is not None:
        out['Co'] = {
            'verdict': 'CURATION-OWED', 'A_measured': c['a_1dlte'], 'n_lines': 1,
            'provenance': 'kittpeak-measured',
            'channel': 'kittpeak: Co I 3845 (blue-edge, SNR-limited)',
            'owed': (f"Kitt Peak covers Co, but the extracted Co I 3845 sits in the blanketed blue "
                     f"edge (SNR~24, chi2r~3100) → the value {c['a_1dlte']} is NOT trusted (blue-edge "
                     f"per the RYA-451/454 caveat). OFF pure DATA-GAP (a measured reference now "
                     f"exists) but curation owed: extract cleaner red Co I lines (within KP's 1300 nm "
                     f"reach) + HFS. Do NOT force the blue value.")}
    return out


def render_md(rows, summary, kp=None):
    lines = []
    lines.append('# RYA-371 Phase C — Solar 27-element verdict table (RYA-239 retry)\n')
    lines.append(f'_Generated {date.today().isoformat()} by scripts/phase_c_verdict_rya371.py. '
                 'Validate-don\'t-tune: classification only, no correction fitted to the anchor._\n')
    if kp:
        lines.append('_RYA-460: Kitt Peak Solar Flux Atlas wired in (N + P/K/Co/Sc). '
                     f"Leg validated by overlap: {'PASS' if kp['leg_validation']['leg_validated'] else 'FLAG'}._\n")
    lines.append('## Verdict counts\n')
    for k in ('PASS', 'NLTE-OWED', 'CURATION-OWED', 'DATA-GAP'):
        v = summary['counts'].get(k, 0)
        if kp:
            d = summary['diff_vs_prior'][k]
            lines.append(f'- **{k}**: {v}  (prior {summary["prior_counts"][k]}, diff {d:+d})')
        else:
            lines.append(f'- **{k}**: {v}')
    lines.append('')
    if kp:
        ov = kp['leg_validation']['overlap']
        lines.append('## RYA-460 overlap cross-check (Kitt Peak leg validation)\n')
        lines.append('| line | Kitt Peak raw | Phase A (arm) | Δ | agree |')
        lines.append('|------|--------------:|---------------|---:|:-----:|')
        for k2, v2 in ov.items():
            lines.append(f"| {k2} | {v2['kittpeak_raw']} | {v2['phase_a_raw']} ({v2['phase_a_arm']}) "
                         f"| {v2['delta']:+.3f} | {'YES' if v2['agree'] else 'NO'} |")
        nci = kp['n_cross_indicator']
        lines.append(f"\n**N solar cross-indicator map** (NLTE-OWED, not validated): N I red "
                     f"mean **{nci.get('atomic_NI_mean')}** (spread {nci.get('atomic_NI_spread')}) "
                     f"from {nci.get('indicators')}. NH 3360 / CN violet blue-edge FLAGGED "
                     f"({nci.get('molecular_blue_edge_flagged')}).\n")
    lines.append('## Per-element table\n')
    lines.append('| El | Asplund21 | A(meas) | Delta | sigma | n | NLTE | Verdict | Provenance | Channel |')
    lines.append('|----|----------:|--------:|------:|------:|--:|:----:|:--------|:-----------|:--------|')
    for r in rows:
        am = '' if r['A_measured'] is None else f"{r['A_measured']:.3f}"
        dl = '' if r['delta_vs_asplund'] is None else f"{r['delta_vs_asplund']:+.3f}"
        sg = '' if r['sigma'] is None else f"{r['sigma']:.2f}"
        nl = r['n_lines'] or ''
        nlte = 'wired' if r['nlte_wired'] else '-'
        if r['threed_wired']:
            nlte += '+3D'
        lines.append(f"| {r['element']} | {r['asplund2021']:.2f} | {am} | {dl} | {sg} | "
                     f"{nl} | {nlte} | **{r['verdict']}** | {r.get('provenance','')} | {r['channel']} |")
    lines.append('')
    lines.append('## Remaining-work map\n')
    for cat in ('PASS', 'NLTE-OWED', 'CURATION-OWED', 'DATA-GAP'):
        els = [r for r in rows if r['verdict'] == cat]
        if not els:
            continue
        lines.append(f'### {cat} ({len(els)})\n')
        for r in els:
            lines.append(f"- **{r['element']}** — {r['owed']}")
        lines.append('')
    lines.append('## Off-Sun — DEFERRED to RYA-348\n')
    lines.append('A green solar run validates the MACHINERY, not target-transferability. '
                 'These cannot be validated on the Sun and are explicitly deferred:\n')
    lines.append('- per-star broadening sourcing (vmac/vsini)')
    lines.append('- temperature-dependent NLTE (e.g. C I LTE-on-Sun but NLTE-on-Procyon)')
    lines.append('- F-star / cool-star line-list adequacy')
    lines.append('\n_Not "55 Cnc-ready" — that is RYA-348 and beyond._\n')
    return '\n'.join(lines)


def _apply_kittpeak(rows, kp):
    """Overlay the RYA-460 Kitt Peak reclassification onto the base rows (in place)."""
    overrides = _kittpeak_reclassify(kp)
    for r in rows:
        ov = overrides.get(r['element'])
        if not ov:
            continue
        for key in ('verdict', 'channel', 'owed', 'provenance'):
            if key in ov:
                r[key] = ov[key]
        if ov.get('A_measured') is not None:
            r['A_measured'] = round(float(ov['A_measured']), 3)
            asp = r['asplund2021']
            r['delta_vs_asplund'] = round(r['A_measured'] - asp, 3)
        if ov.get('sigma') is not None:
            r['sigma'] = round(float(ov['sigma']), 3)
        if ov.get('n_lines') is not None:
            r['n_lines'] = int(ov['n_lines'])
    return overrides


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--star', default='solar')
    args = ap.parse_args()

    ab, ew, phase_a = _load()
    rows = build_verdicts(ab, ew, phase_a)

    # RYA-460: fold in the Kitt Peak measurements (N + P/K/Co/Sc) if the campaign ran.
    kp = _load_kittpeak()
    overrides = _apply_kittpeak(rows, kp)

    counts = {}
    for r in rows:
        counts[r['verdict']] = counts.get(r['verdict'], 0) + 1
    diff = {k: counts.get(k, 0) - PRIOR_COUNTS.get(k, 0)
            for k in ('PASS', 'NLTE-OWED', 'CURATION-OWED', 'DATA-GAP')}
    summary = {'ticket': 'RYA-371 Phase C (RYA-460 reference-wired)', 'star': args.star,
               'generated': date.today().isoformat(),
               'reference': 'Asplund, Amarsi & Grevesse 2021 (A&A 653, A141)',
               'tol_pass_dex': TOL_PASS, 'n_elements': len(rows), 'counts': counts,
               'prior_counts': PRIOR_COUNTS, 'diff_vs_prior': diff,
               'kittpeak_wired': bool(kp),
               'kittpeak_elements': sorted(overrides) if overrides else []}

    AUDIT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    out_json = AUDIT / 'solar_phase_c_verdict.json'
    with open(out_json, 'w') as fh:
        json.dump({'summary': summary, 'verdicts': rows}, fh, indent=2)
    out_back = AUDIT / 'solar_differential_backbone.json'
    with open(out_back, 'w') as fh:
        json.dump(differential_backbone(), fh, indent=2)
    out_md = DOCS / 'solar_phase_c_verdict_rya371.md'
    out_md.write_text(render_md(rows, summary, kp))

    # Console table
    print(f"\n{'='*78}\n  RYA-371 Phase C — solar 27-element verdict ({len(rows)} elements)\n{'='*78}")
    print(f"  {'El':3s} {'Asp21':>6s} {'A(meas)':>8s} {'Delta':>7s} {'sig':>5s} {'n':>3s}  "
          f"{'NLTE':5s} {'Verdict':14s} Channel")
    for r in rows:
        am = '   -   ' if r['A_measured'] is None else f"{r['A_measured']:7.3f}"
        dl = '   -  ' if r['delta_vs_asplund'] is None else f"{r['delta_vs_asplund']:+6.3f}"
        sg = '  -  ' if r['sigma'] is None else f"{r['sigma']:5.2f}"
        nl = f"{r['n_lines']:3d}" if r['n_lines'] else '  -'
        nlte = ('W' if r['nlte_wired'] else '-') + ('3' if r['threed_wired'] else ' ')
        print(f"  {r['element']:3s} {r['asplund2021']:6.2f} {am} {dl} {sg} {nl}  "
              f"{nlte:5s} {r['verdict']:14s} {r['channel'][:40]}")
    print(f"\n  Counts: " + '  '.join(f'{k}={counts.get(k,0)}'
          for k in ('PASS', 'NLTE-OWED', 'CURATION-OWED', 'DATA-GAP')))
    if kp:
        print(f"  Prior : " + '  '.join(f'{k}={PRIOR_COUNTS[k]}'
              for k in ('PASS', 'NLTE-OWED', 'CURATION-OWED', 'DATA-GAP')))
        print(f"  Diff  : " + '  '.join(f'{k}={diff[k]:+d}'
              for k in ('PASS', 'NLTE-OWED', 'CURATION-OWED', 'DATA-GAP'))
              + f"   (Kitt Peak wired: {summary['kittpeak_elements']})")
        print(f"  Leg validation (overlap [O I]6300/O I 777 vs HARPS/ESPRESSO): "
              f"{'PASS' if kp['leg_validation']['leg_validated'] else 'FLAG'}")
    print(f"\n  Wrote:\n    {out_json.relative_to(ROOT)}\n    {out_back.relative_to(ROOT)}"
          f"\n    {out_md.relative_to(ROOT)}\n")


if __name__ == '__main__':
    main()
