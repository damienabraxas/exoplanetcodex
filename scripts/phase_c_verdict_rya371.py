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
from pipeline import data_namespace as ns  # noqa: E402  RYA-469 gold solar reference
from pipeline.provenance_honesty import (  # noqa: E402  RYA-596 tripwire, shared (RYA-653)
    ZERO_SURVIVOR_CHANNEL, assert_blank_cause_is_honest)
from pipeline.solar_scale_provenance import (  # noqa: E402  RYA-681 value-keyed idempotency
    apply_reported_scale_correction)
from pipeline.ratified_constraints import (  # noqa: E402  RYA-674 emission-time gate
    assert_ratified_constraints_satisfied)

PROC = ROOT / 'data' / 'processed'
AUDIT = ROOT / 'data' / 'audit' / 'cno_synthesis'
DOCS = ROOT / 'docs' / 'audit'
KITTPEAK_JSON = AUDIT / 'solar_kittpeak_rya460.json'   # RYA-460 KP measurements
CU_V_SYNTH_JSON = (ROOT / 'data' / 'audit' / 'cu_v_hfs_synthesis' /
                   'solar_cu_v_hfs_synthesis_rya466.json')   # RYA-466 HFS-synthesis Cu/V
MN_SYNTH_JSON = (ROOT / 'data' / 'audit' / 'mn_hfs_synthesis' /
                 'solar_mn_hfs_synthesis_rya473.json')       # RYA-473 HFS-synthesis Mn
S_SYNTH_JSON = ROOT / 'data' / 'results' / 'solar_s_costasilva_rya492.json'  # RYA-492 CS-gf S
BA_SYNTH_JSON = ROOT / 'data' / 'results' / 'solar_ba_synthesis_rya559.json'  # RYA-559 Ba II 5853 synth
BA_DEBLEND_JSON = ROOT / 'data' / 'results' / 'solar_ba_deblend_rya581.json'  # RYA-581 in-window deblend
CO_SYNTH_JSON = ROOT / 'data' / 'results' / 'co_synthesis_rya564.json'   # RYA-564 Co I red-line synth

# Prior Phase C verdict counts — the immediate baseline this run diffs against.
# RYA-462 diffs against RYA-460's reference-wired verdict (3 / 2 / 21 / 0); RYA-460 in
# turn moved it off the pre-Kitt-Peak 3 / 1 / 18 / 4 (DATA-GAP eliminated).
PRIOR_COUNTS = {'PASS': 3, 'NLTE-OWED': 2, 'CURATION-OWED': 21, 'DATA-GAP': 0}

# Tolerance for a PASS against Asplund 2021. The Sun's 1D-LTE/NLTE absolute scale
# carries a documented ~+0.05 dex zero-point vs Asplund's 3D values for the
# validated Fe leg (RYA-336/407); use 0.10 dex as the per-element reconciliation
# band (NOT a fit target — a line above which we stop calling a leg "reconciled").
TOL_PASS = 0.10

# Elements whose only solar channel is FUV/near-UV-hard or otherwise unreachable
# from the present ESO ground-based set (per the ticket). P is the flagship DATA-GAP.
GROUND_UNREACHABLE = {'P'}     # FUV/near-UV, needs HST/STIS (RYA-119)


def _load(gold_version: str = 'CURRENT'):
    """Load the verdict's inputs. Returns (gold_df, ew_pool, phase_a, resolved_version).

    RYA-469: the solar verdict is computed from a FROZEN gold reference, not a
    regenerable working file — a perturbed regen can never become the baseline.

    RYA-674: WHICH frozen reference is now a PARAMETER, defaulting to CURRENT so
    nothing else changes. It was hardcoded, and once RYA-681 made gold v3's
    self-contradictory Fe row a loud refusal to load, that hardcoding meant the verdict
    artifact could not be regenerated AT ALL — and since the RYA-654 element status
    tracker is generated from that committed artifact, the tracker could not be updated
    after any element was fixed. Naming a different frozen input is auditable (the
    resolved version is stamped into the emitted summary and every artifact header);
    a flag that skipped the scale guard, or an in-memory label repair, would be the
    silent-correction pattern RYA-681 just spent a session removing.
    """
    ab, resolved = ns.read_solar_reference(gold_version)
    ew = pd.read_csv(ROOT / 'data' / 'measured' / 'sol_ew_results_v1.csv')
    with open(AUDIT / 'solar_phase_a_cross_arm.json') as fh:
        phase_a = json.load(fh)
    return ab, ew, phase_a, resolved


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


FE1_RESID_RYA407 = (ROOT / 'data' / 'audit' / 'fe1_scatter' /
                    'fe1_per_line_residuals_rya407.csv')     # RYA-407 honest-floor pool


def _fe1_raw_scatter():
    """Fe I line-to-line raw scatter from the canonical RYA-407 honest-floor pool.

    The frozen gold reference carries no per-line Fe scatter column, so a regenerated
    verdict would otherwise report sigma=null and the RYA-166 A4 gate could not read
    it. Restore the scatter from its canonical source, computed exactly as the pipeline
    does (np.nanstd, numpy default ddof=0 — the same 0.139 that A_X_std produced when
    the gold was first built), which sits at/under the 0.1398 ddof=1 floor
    (ACCEPTANCE_PROFILES['G']['fe1_scatter_max'], RYA-407/446). Returns NaN if absent."""
    if not FE1_RESID_RYA407.exists():
        return float('nan')
    pool = pd.read_csv(FE1_RESID_RYA407)['a_1dlte'].astype(float).values
    return float(np.nanstd(pool)) if pool.size > 1 else float('nan')


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

        # ── RYA-553: apply the tabulated Magic-2013 1D→3D solar Fe correction ──
        # The reported Fe I anchor is the value the RYA-166 gate and the RYA-527
        # gold-v3 freeze both read. Our NLTE grids emit on the 1D-NLTE scale
        # (~+0.05 above 3D-true); add the NEGATIVE tabulated offset, AFTER NLTE, to
        # move the reported anchor onto the true 7.46 scale. SOLAR ONLY (off-solar
        # per-Teff/[Fe/H] generalisation is owed, RYA-550). Never silent: logged +
        # recorded pre/post on the row.
        #
        # RYA-681 — idempotency is keyed on the VALUE, not on a prose label. The old
        # guard was `if '3D' not in mrow['method_scale'].upper()`, i.e. idempotent
        # with respect to a string that the gold BUILDER writes (and hardcoded to
        # '1D-NLTE (Fe I)'), not with respect to the number it labels. RYA-665 froze
        # the post-correction 7.466 under that pre-correction label and the guard
        # re-armed → RYA-669 measured 7.416, with every gate green. The decision now
        # lives in pipeline.solar_scale_provenance, which reads the gold row's
        # explicit `scale_state` declaration, CROSS-CHECKS it against where the value
        # actually sits between the two published scale centres, and RAISES on
        # contradiction instead of guessing. A desynchronised freeze is now a loud
        # load failure, not a silent second subtraction.
        fe_1d3d = None
        if el == 'Fe' and np.isfinite(a_meas):
            a_meas, fe_1d3d = apply_reported_scale_correction(el, a_meas, mrow)
            if fe_1d3d['applied']:
                print(f"  RYA-553 Fe 1D→3D: A(Fe I) {fe_1d3d['a_1dnlte_pre']:.3f} (1D-NLTE) "
                      f"{fe_1d3d['correction_dex']:+.3f} -> {a_meas:.3f} (3D-NLTE, Magic 2013)"
                      f"  [scale from {fe_1d3d['gold_scale_source']}, value-corroborated]")
            else:
                print(f"  RYA-553 Fe 1D→3D: SKIPPED — {fe_1d3d['reason']}; "
                      f"A(Fe I) stays {a_meas:.3f}")
            # RYA-407/446: the frozen gold has no per-line Fe scatter — restore it.
            if not np.isfinite(sigma):
                sigma = _fe1_raw_scatter()

        # RYA-456: the curation's BLIND verdict (VALIDATED / RESIDUAL / LOW_CONFIDENCE),
        # carried on the wired non-Fe rows. The classifier MAPS it (it never re-derives
        # a threshold), so the science decision stays in curate_nonfe_pools (RYA-395/398).
        cverdict = (str(mrow['curation_verdict'])
                    if (mrow is not None and 'curation_verdict' in mrow
                        and isinstance(mrow.get('curation_verdict'), str)) else '')
        # RYA-596: the gold row's ratified confidence TIER (RYA-522). An `owed`-tier
        # row freezes NO value BY DESIGN even when the curation produced one — so a
        # blank A_measured here means "held by the tier", NOT "the graded cull left
        # nothing". The classifier must not conflate the two (see _classify).
        gold_tier = (str(mrow['confidence'])
                     if (mrow is not None and 'confidence' in mrow
                         and isinstance(mrow.get('confidence'), str)) else '')
        gold_note = (str(mrow['note'])
                     if (mrow is not None and 'note' in mrow
                         and isinstance(mrow.get('note'), str)) else '')
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
                                            cverdict, charter, gold_tier, gold_note)
        if ew_integrity_note:
            owed = f"{owed} [{ew_integrity_note}]"
        _assert_blank_cause_is_honest(el, channel, n_lines, a_measured=a_meas,
                                      site='phase_c verdict (RYA-371)')
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
            'fe_1d3d_correction': fe_1d3d,   # RYA-553: solar Fe 1D→3D pre/post (None for non-Fe)
        })
    # RYA-674 §2C: the emission-time gate on ratified constraints. Runs on the BASE
    # assembly here and again in main() after the dedicated-channel overlays, because
    # an overlay is itself an emission — the Li 1.409 leak was an overlay-shaped move.
    assert_ratified_constraints_satisfied(rows, 'phase_c verdict generator (RYA-371)')
    return rows


# RYA-653: the tripwire and the claim it guards now live in ONE place
# (pipeline/provenance_honesty.py) so the gold reference builder enforces the
# same invariant from the same code, not from a forked copy. Re-exported under
# their original names — this module's callers and tests are unchanged.
_assert_blank_cause_is_honest = assert_blank_cause_is_honest


def _classify(el, asp, a_meas, delta, sigma, n_lines, grid, threed, nlte_flag,
              in_pool, produced, phase_a, cverdict='', charter=None,
              gold_tier='', gold_note=''):
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
                f'A(Fe I) 3D-NLTE {a_meas:.3f} vs Asplund 7.46 ({delta:+.3f}); ionization-balance '
                'gated, scatter 0.139 = honest floor (RYA-407). The +0.05 1D→3D solar offset '
                '(Magic 2013) is now APPLIED at the reported layer (RYA-553), so the anchor sits '
                'on the true 3D scale; FE_GATE [7.41,7.51] governs.')

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
    if in_pool:
        nlte_note = (f'NLTE grid available ({grid})' if grid else 'no NLTE grid (would be LTE-flagged)')
        # RYA-458: surface the Eu 6645 HFS recovery disposition (RYA-102) on its row.
        eu_note = ''
        if el == 'Eu':
            euc = charter.get('Eu_6645', {})
            if euc.get('present'):
                eu_note = (f" Eu II 6645 EW {euc['ew_mA']:.1f} mA, ew_integrity "
                           f"disposition={euc['ew_disposition']} (RYA-102/458 HFS-summing).")

        # ── RYA-596: TWO distinct blank causes, never conflated ──
        # This branch used to assert one cause unconditionally — "the RYA-398 graded-gf
        # firewall culls every line" — for ANY in-pool row with a blank A(X). That claim
        # was FALSE for Ca/Ti/Ni/Na/Al/Sr and was self-contradicted by the `n_lines` on
        # the very same row (2/10/2/2/1/1 graded survivors). The real cause there is the
        # RYA-522 gold tiering: an `owed`-tier row freezes NO value BY RATIFIED DESIGN
        # (Ryan 2026-07-05, build_solar_reference_v2_rya522.py — "suspect → held, not
        # immortalised"), and this verdict READS the frozen gold back in (RYA-469). The
        # blank is a deliberate HOLD round-tripping through the gold, not a cull.
        # n_lines here IS the curated graded-survivor count carried on the gold row, so
        # it is the evidence that decides which cause we may state.
        if n_lines > 0:
            held = (f'A(X) was produced by the RYA-395/398 graded cull on {n_lines} surviving '
                    f'independent-gf line(s) and is HELD UNFROZEN at the ratified gold tier '
                    f"'{gold_tier or 'owed'}' (RYA-522: a value is frozen only if we would stake "
                    f'a differential on it). This is a deliberate hold, NOT a zero-survivor cull '
                    f'— the graded firewall did leave lines to stand on. To surface a number the '
                    f'tier must be re-ratified (RYA-522), or the pool broadened via RYA-161/162 '
                    f'(differential survey).')
            if gold_note:
                held += f' Gold note: "{gold_note}".'
            return ('CURATION-OWED',
                    f'EW: {n_lines} curated line(s); value HELD at gold tier '
                    f"'{gold_tier or 'owed'}' (RYA-522) — not a graded-cull blank",
                    f'{held} {nlte_note}.{eu_note}')

        # n_lines == 0 → the graded firewall genuinely left nothing (Mg, Y, Zr, Eu).
        return ('CURATION-OWED', ZERO_SURVIVOR_CHANNEL,
                'solar EW measured + matched in linelist_solar, but the RYA-398 graded-gf '
                'firewall (wired into the default run, RYA-456) culls every line — zero '
                'independent-gf survivors (pool gf is Kurucz/ungraded, or every line fails '
                'SAT/HIERR/BLEND). gf-data-limited → RYA-161/162 (differential survey). '
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


def _apply_nlte_grid_delta(element, wave_A, a_1dlte, star=None):
    """RYA-462: apply the vendored NLTE grid delta to a (Kitt Peak) 1D-LTE measurement
    through the EXISTING interpolation subsystem (pipeline.nlte_corrections) — the same
    path Ca/Ti/Cr/Na travel. Returns (a_nlte, delta, flag). If the element/line is out
    of grid coverage it returns (a_1dlte, nan, 'NLTE_unavailable ...') and NEVER silently
    corrects. Validate-don't-tune: delta is READ from the grid, never fitted to Asplund."""
    star = star or {'teff': 5772.0, 'logg': 4.44, 'feh': 0.0}
    try:
        from pipeline import nlte_corrections as N
        if element not in NLTE_CORRECTION_ELEMENTS:
            return float(a_1dlte), float('nan'), 'NLTE_unavailable (element not registered)'
        if not N.element_grid_in_bounds(element, star['teff'], star['logg'], star['feh']):
            return float(a_1dlte), float('nan'), 'NLTE_unavailable (star out of grid hull)'
        d = N._mpia_element_delta(element, wave_A, star['teff'], star['logg'], star['feh'])
        if d is None or not np.isfinite(d):
            return float(a_1dlte), float('nan'), 'NLTE_unavailable (no grid node within tol)'
        flag = NLTE_CORRECTION_ELEMENTS[element].get('flag', 'NLTE_1D')
        return float(a_1dlte) + float(d), float(d), flag
    except Exception as e:                            # never let the wiring break the verdict
        return float(a_1dlte), float('nan'), f'NLTE_unavailable ({e})'


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

    # ── N — measured from Kitt Peak N I red; NLTE grid now WIRED (RYA-556) ──
    # The N_Amarsi2020_PySME grid was PRESENT-but-UNWIRED here: registered in
    # NLTE_CORRECTION_ELEMENTS (RYA-369/526) but this KP channel emitted the 1D-LTE
    # mean and left N in NLTE-OWED (the merged-not-wired gap RYA-527 surfaced). Apply
    # its vendored solar delta to each of the 3 KP N I multiplets through the SAME
    # interpolation subsystem the other registry elements use (Ca/Ti/Cr/Na/K), and
    # average — the near-LTE correction (per-line ~-0.011/-0.015/-0.015) clears the
    # NLTE debt. Validate-don't-tune: the delta is READ from the grid, never fitted.
    # N does NOT become PASS: the +0.36 residual vs Asplund is a KP red-multiplet
    # gf/data-channel floor (RYA-161), routed to curation, NOT an NLTE debt.
    N_GRID_NODES = {'NI_7442_7468': 7468.31, 'NI_8216_8223': 8216.34, 'NI_8680_8718': 8683.4}
    if leg_ok and nci.get('atomic_NI_mean') is not None:
        a_lte = float(nci['atomic_NI_mean'])
        asp_n = float(SOLAR_ASPLUND2021.get('N', 7.83))
        per_line = [(_apply_nlte_grid_delta('N', w, a_lte), key)
                    for key, w in N_GRID_NODES.items()]
        deltas = [pl[0][1] for pl in per_line]           # (a_nlte, delta, flag)
        nflag = per_line[0][0][2]
        if all(np.isfinite(d) for d in deltas):          # grid resolved for every multiplet
            nd = float(np.mean(deltas))
            a_nlte = round(a_lte + nd, 3)
            print(f"  RYA-556 N I NLTE: A(N) {a_lte:.3f} (1D-LTE) {nd:+.4f} "
                  f"-> {a_nlte:.3f} (N_Amarsi2020_PySME, RYA-369/526) — off NLTE-OWED")
            out['N'] = {
                'verdict': 'CURATION-OWED', 'A_measured': a_nlte,
                'sigma': nci['atomic_NI_spread'], 'n_lines': 3, 'provenance': 'kittpeak-measured',
                'channel': 'kittpeak: N I red 7468/8216/8683 — NLTE-wired (N_Amarsi2020_PySME, RYA-369/526)',
                'owed': (f"MEASURED from Kitt Peak N I red — 3 independent multiplets AGREE "
                         f"(1D-LTE mean {a_lte:.3f}, spread {nci['atomic_NI_spread']}). N I NLTE "
                         f"delta APPLIED via the registered N_Amarsi2020_PySME grid (RYA-369/526) "
                         f"through the existing interpolation subsystem (RYA-556 wiring; {nflag}; "
                         f"validate-don't-tune): per-line {deltas[0]:+.4f}/{deltas[1]:+.4f}/"
                         f"{deltas[2]:+.4f}, mean {nd:+.4f} -> A(N) {a_nlte:.3f} ({a_nlte - asp_n:+.3f} "
                         f"vs Asplund {asp_n:.2f}). The NLTE debt is now CLEARED (off NLTE-OWED). The "
                         f"remaining +{a_nlte - asp_n:.2f} is a KP red-multiplet gf/data-channel floor "
                         f"(RYA-161) — curation owed, NOT an NLTE debt; do NOT tune. NOT validated: "
                         f"Teff-bracket owed (Procyon / aCen B, RYA-369). NH 3360 + CN violet 3883 "
                         f"UNMEASURABLE here — blue-edge no-true-continuum (SNR~28, RYA-451/454) + the "
                         f"Turbospectrum molecular linelist is absent — FLAGGED, not forced. {leg_txt}.")}
        else:
            # Grid unavailable at runtime -> LOUD, held NLTE-OWED, never silently LTE.
            print(f"  RYA-556 N I NLTE: grid unavailable ({nflag}) — N held NLTE-OWED, "
                  f"never silently LTE (delta not applied).")
            out['N'] = {
                'verdict': 'NLTE-OWED', 'A_measured': a_lte,
                'sigma': nci['atomic_NI_spread'], 'n_lines': 3, 'provenance': 'kittpeak-measured',
                'channel': 'kittpeak: N I red 7442/7468 + 8216/8223 + 8680-8718 (NH/CN blue-edge flagged)',
                'owed': (f"MEASURED from Kitt Peak N I red (1D-LTE mean {a_lte:.3f}, spread "
                         f"{nci['atomic_NI_spread']}). N_Amarsi2020_PySME grid registered "
                         f"(RYA-369/526) but the delta could NOT be interpolated at runtime "
                         f"({nflag}) -> held NLTE-OWED, never silently LTE. {leg_txt}.")}

    # ── K — measured; NLTE grid now WIRED (RYA-462) ──
    # The K_Amarsi2020_PySME grid was PRESENT-but-UNWIRED (grid + full PySME machinery
    # existed, K just wasn't in NLTE_CORRECTION_ELEMENTS). RYA-462 registers it; here we
    # apply its vendored solar delta to the Kitt Peak K I 7699 1D-LTE value through the
    # same interpolation subsystem the other registry elements use. Validate-don't-tune:
    # the delta is read from the grid (solar 7699 ~ -0.31), not fitted to Asplund.
    if leg_ok and (k := m.get('KI_7665_7699')) and k['a_1dlte'] is not None:
        asp_k = float(k.get('asplund2021', SOLAR_ASPLUND2021.get('K', 5.07)))
        a_nlte, kd, kflag = _apply_nlte_grid_delta('K', 7698.964, k['a_1dlte'])
        if np.isfinite(kd):
            reconciled = abs(a_nlte - asp_k) <= TOL_PASS
            out['K'] = {
                'verdict': 'PASS' if reconciled else 'NLTE-OWED',
                'A_measured': a_nlte, 'n_lines': 1, 'provenance': 'kittpeak-measured',
                'channel': 'kittpeak: K I 7699 (clean; 7665 in the telluric O2 A-band) — NLTE-wired',
                'owed': (f"MEASURED Kitt Peak K I 7699 = {k['a_1dlte']} (1D-LTE, {k['delta_vs_asplund']:+.2f} "
                         f"vs {asp_k:.2f}). K_Amarsi2020_PySME NLTE delta {kd:+.3f} APPLIED via the existing "
                         f"interpolation subsystem (RYA-462 wiring; {kflag}; validate-don't-tune) -> A(K) "
                         f"{a_nlte:.3f} ({a_nlte - asp_k:+.3f} vs Asplund {asp_k:.2f}). "
                         + ("Reconciles within TOL after the cited NLTE correction — the severe negative "
                            "K I resonance NLTE is real, not tuned."
                            if reconciled else
                            "An NLTE residual survives -> still owed (do NOT tune)."))}
        else:
            out['K'] = {
                'verdict': 'NLTE-OWED', 'A_measured': k['a_1dlte'], 'n_lines': 1,
                'provenance': 'kittpeak-measured',
                'channel': 'kittpeak: K I 7699 (clean; 7665 sits in the telluric O2 A-band)',
                'owed': (f"MEASURED Kitt Peak K I 7699 = {k['a_1dlte']} ({k['delta_vs_asplund']:+.2f} vs "
                         f"{asp_k:.2f}). K_Amarsi2020_PySME grid registered (RYA-462) but the delta could "
                         f"not be interpolated here ({kflag}) -> held NLTE-OWED, never silently LTE.")}

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

    # ── Co — the blue-edge 3845 line is DIAGNOSTIC-ONLY (RYA-564) ──
    # RYA-460 originally reported A(Co) from the Kitt Peak Co I 3845 extraction. That value
    # (6.128, +1.188 vs Asplund) is an artifact of the blanketed blue edge (SNR~24,
    # chi2r~3100) — the same class as the Sr I +2.13 the RYA-524 audit was built to catch —
    # and the registry already marked it "value NOT trusted". RYA-564 DEMOTES it here, at
    # source and UNCONDITIONALLY: the line stays on the record as a diagnostic, but it can
    # never be the reported A(Co) again, whether or not the red-line synthesis has run.
    # The reportable value comes from _co_reclassify (clean red Co I lines); if that has not
    # run, or no red line clears reliability, Co honestly reports NO VALUE.
    if leg_ok and (c := m.get('CoI_3845')) and c['a_1dlte'] is not None:
        out['Co'] = {
            'verdict': 'CURATION-OWED', 'A_measured_blank': True, 'n_lines': 0,
            'provenance': 'diagnostic-only (no reported value)',
            'channel': 'kittpeak: Co I 3845 — DIAGNOSTIC-ONLY (blue-edge artifact, demoted RYA-564)',
            'owed': (f"NO VALUE REPORTED from Kitt Peak. The extracted Co I 3845 sits in the "
                     f"blanketed blue edge (SNR~24, chi2r~3100); its A(Co)={c['a_1dlte']} "
                     f"({c['delta_vs_asplund']:+.3f} vs 4.94) is NOT credible solar Co — the same "
                     f"artifact class as the Sr I +2.13 (RYA-524). DEMOTED to diagnostic-only by "
                     f"RYA-564; it must never enter the freeze. The measurement route is "
                     f"HFS-resolved synthesis on the clean RED Co I lines (RYA-564); until that "
                     f"lands a reliable line, Co is owed WITHOUT a value — never the blue value.")}
    return out


def _load_cu_v_synthesis():
    """RYA-466: the HFS-resolved synthesis Cu/V measurement, or None if it hasn't run."""
    if not CU_V_SYNTH_JSON.exists():
        return None
    return json.loads(CU_V_SYNTH_JSON.read_text())


def _cu_v_reclassify(data):
    """RYA-466 — fold the HFS-resolved synthesis measurement of Cu/V into the verdict.
    These elements are UNMEASURABLE by the EW path (hyperfine-split → the single-profile
    EW fit can't reach them; they never enter the gf-curation pool — the RYA-354 finding).
    Synthesis on the GES HFS-resolved line list measures them. The verdicts are driven by
    the measured values, not tuned: a synthesised value that still sits high is a finding,
    not a PASS. Returns {element: override-dict}."""
    if not data:
        return {}
    out = {}
    cu = data.get('Cu', {})
    if cu.get('A_nlte') is not None:
        a = float(cu['A_nlte']); asp = float(cu['asplund2021']); d = a - asp
        nlte = ('live RYA-402 b-factor' if cu.get('nlte_live')
                else 'vendored RYA-402 b-factor (live .grd offline)')
        reconciled = abs(d) <= TOL_PASS
        out['Cu'] = {
            'verdict': 'PASS' if reconciled else 'CURATION-OWED',
            'A_measured': a, 'sigma': cu.get('scatter'), 'n_lines': cu.get('n_lines'),
            'provenance': 'synthesis: HFS-resolved (RYA-466)',
            'channel': f"HFS synthesis: Cu I {cu.get('n_lines')} lines (5105/5218/5220/5700/5782), "
                       f"gf=Kock&Richter, NLTE {nlte}",
            'owed': (f"MEASURED via HFS-resolved synthesis — the EW path could not reach Cu "
                     f"(all 5 lines hyperfine-split, RYA-354 finding); synthesis on the GES "
                     f"HFS line list (4-10 components/feature) measures it. gf adjudicated to "
                     f"GES=Kock&Richter (KR), VALD3 superseded ({cu.get('gf_provenance','')}). "
                     f"A(Cu)_LTE {cu.get('A_lte_median')} + RYA-402 b-factor NLTE "
                     f"{cu.get('nlte_delta'):+.3f} ({nlte}) = {a:.3f} ({d:+.3f} vs Asplund "
                     f"{asp:.2f}; σ {cu.get('scatter')}, n={cu.get('n_lines')}). "
                     + ("Reconciles within TOL — the HFS-synthesis path closes the Cu gap."
                        if reconciled else
                        "The +offset survives the small Cu NLTE → 1D-LTE optical Cu I sits high; "
                        "3D / fuller-NLTE curation owed (RYA-161/162). The blocker was the "
                        "MEASUREMENT TOOL, now fixed; the residual is a finding, do NOT tune."))}
    v = data.get('V', {})
    if v.get('A_lte_median') is not None:
        a = float(v['A_lte_median']); asp = float(v['asplund2021']); d = a - asp
        out['V'] = {
            'verdict': 'CURATION-OWED',
            'A_measured': a, 'sigma': v.get('scatter'), 'n_lines': v.get('n_lines'),
            'provenance': 'synthesis: HFS-resolved LTE (RYA-466)',
            'channel': f"HFS synthesis LTE: V I {v.get('n_lines')} lines — NLTE-VOID (no model atom)",
            'owed': (f"MEASURED via HFS-resolved LTE synthesis — V was unmeasurable by the EW path "
                     f"(HFS, RYA-354). A(V)_LTE {a:.3f} ({d:+.3f} vs Asplund {asp:.2f}; σ "
                     f"{v.get('scatter')}, n={v.get('n_lines')}) — tight, AGREES with Asplund, BUT "
                     f"V is the NLTE-VOID (no model atom anywhere, RYA-463) → LTE-only, LOWER "
                     f"CONFIDENCE; a V NLTE model + graded gf (V I Lawler+2014, V II Wood+2014) "
                     f"owed before any PASS. Off no-value, do NOT certify on LTE agreement alone.")}
    return out


def _load_mn_synthesis():
    """RYA-473: the HFS-resolved synthesis Mn measurement, or None if it hasn't run."""
    if not MN_SYNTH_JSON.exists():
        return None
    return json.loads(MN_SYNTH_JSON.read_text())


def _mn_reclassify(data):
    """RYA-473 — fold the HFS-resolved synthesis measurement of Mn into the verdict.
    Mn was no-value: its gf is graded (Den Hartog, RYA-468) but the EW path SAT-culls the
    e6S→z6P triplet (REW −4.78..−4.82 over the −4.90 knee, HFS-split). Synthesis on the GES
    HFS-resolved line list (6 components/feature = Den Hartog Table 4) measures it. Driven
    by the measured value, not tuned: a synthesised Mn that still sits high is a finding.
    Returns {'Mn': override-dict}."""
    if not data:
        return {}
    mn = data.get('Mn', {})
    a_lte = mn.get('A_lte_median')
    if a_lte is None:
        return {}
    a = float(mn['A_nlte']) if mn.get('A_nlte') is not None else float(a_lte)
    asp = float(mn['asplund2021']); d = a - asp
    has_nlte = mn.get('nlte_delta') is not None
    nlte = ('live Amarsi HFS-resolved' if mn.get('nlte_live')
            else 'vendored MPIA/Bergemann grid (live Amarsi .grd offline)')
    reconciled = abs(d) <= TOL_PASS
    # Mn carries the RYA-411 NLTE caveat: the MPIA grid δ is its high-EP reference value,
    # the low-EP triplet is not a grid node → do NOT certify PASS on the vendored δ alone.
    vendored_caveat = has_nlte and not mn.get('nlte_live')
    verdict = 'PASS' if (reconciled and not vendored_caveat) else 'CURATION-OWED'
    nlte_txt = (f"+ Mn NLTE {mn.get('nlte_delta'):+.3f} ({nlte}) = {a:.3f}"
                if has_nlte else f"= A(Mn)_LTE {a:.3f} (NLTE unavailable)")
    return {'Mn': {
        'verdict': verdict,
        'A_measured': a, 'sigma': mn.get('scatter'), 'n_lines': mn.get('n_lines'),
        'provenance': 'synthesis: HFS-resolved (RYA-473)',
        'channel': f"HFS synthesis: Mn I {mn.get('n_lines')} lines "
                   f"(6013/6016/6021, Den Hartog e6S→z6P), gf=Den Hartog+2011 (MED), "
                   f"NLTE {nlte}",
        'owed': (f"MEASURED via HFS-resolved synthesis — the EW path SAT-culls the Den Hartog "
                 f"triplet (REW −4.78..−4.82 over the −4.90 knee, HFS-split hfs_n=6; RYA-468 "
                 f"finding: gf graded but saturation is the blocker). Synthesis on the GES HFS "
                 f"line list (6 components/feature = Den Hartog Table 4, cited) measures it. "
                 f"A(Mn)_LTE {a_lte} {nlte_txt} ({d:+.3f} vs Asplund {asp:.2f}; σ "
                 f"{mn.get('scatter')}, n={mn.get('n_lines')}). "
                 + ("Reconciles within TOL — the HFS-synthesis path closes the Mn gap."
                    if (reconciled and not vendored_caveat) else
                    ("The NLTE δ is the MPIA grid's HIGH-EP reference value (+0.107); the "
                     "low-EP triplet is NOT a grid node and its HFS-resolved Amarsi δ differs "
                     "(RYA-411, .grd offline) → do NOT certify PASS on the vendored δ. The "
                     "MEASUREMENT-TOOL blocker is fixed (off no-value); the line-exact NLTE is "
                     "owed (RYA-411 Amarsi grid)." if vendored_caveat else
                     "The offset survives → 1D Mn sits off Asplund; fuller NLTE/3D curation "
                     "owed. The blocker was the MEASUREMENT TOOL, now fixed; the residual is a "
                     "finding, do NOT tune.")))}}


def _apply_mn_synthesis(rows, data):
    """Overlay the RYA-473 Mn HFS-synthesis reclassification onto the base rows (in place)."""
    overrides = _mn_reclassify(data)
    for r in rows:
        ov = overrides.get(r['element'])
        if not ov:
            continue
        for key in ('verdict', 'channel', 'owed', 'provenance'):
            if key in ov:
                r[key] = ov[key]
        if ov.get('A_measured') is not None:
            r['A_measured'] = round(float(ov['A_measured']), 3)
            r['delta_vs_asplund'] = round(r['A_measured'] - r['asplund2021'], 3)
        if ov.get('sigma') is not None:
            r['sigma'] = round(float(ov['sigma']), 3)
        if ov.get('n_lines') is not None:
            r['n_lines'] = int(ov['n_lines'])
    return overrides


def _load_s_synthesis():
    """RYA-492: the Costa-Silva-gf synthesis S measurement, or None if it hasn't run."""
    if not S_SYNTH_JSON.exists():
        return None
    return json.loads(S_SYNTH_JSON.read_text())


def _s_reclassify(data):
    """RYA-557 — repoint solar S to the RYA-492 Costa-Silva gf synthesis value (7.486),
    the single cited source for S. The generic EW path leaves S owed-no-value: the EW
    sanity/gf-grade cull rejects all solar S I lines but 6757.15, and the 6743 line that
    carries the CS gf (canonical_gf, RYA-492) is a SYNTHESIS line absent from the EW pool.
    So the CS-gf value lives on the synthesis channel — fold it in like Cu/V/Mn.

    S stays CURATION-OWED: the +0.37 vs Asplund 7.12 is a gf-scale floor (RYA-161), NOT
    closed by the gf (adopting CS gf moved it only -0.03, 7.516->7.486). Validate-don't-
    tune: the value is READ from the committed measurement, never fitted toward Asplund.
    Returns {'S': override-dict}."""
    if not data:
        return {}
    cs = data.get('costa_silva', {})
    a = cs.get('a_nlte')
    if a is None:
        return {}
    a = float(a)
    asp = float(SOLAR_ASPLUND2021.get('S', 7.12))
    d = a - asp
    ctrl = (data.get('control_ges', {}) or {}).get('a_nlte')
    reconciled = abs(d) <= TOL_PASS          # False for S (gf-floor) — never PASS
    return {'S': {
        'verdict': 'PASS' if reconciled else 'CURATION-OWED',
        'A_measured': a, 'sigma': cs.get('sigma'), 'n_lines': len(cs.get('per_window', [])) or 2,
        'provenance': 'synthesis: Costa-Silva gf (RYA-492)',
        'channel': 'synthesis: S I 6743.53 + 6757.15 windows, gf=Costa Silva+2020 '
                   '(A&A 634 A136) Table1, NLTE Amarsi 2025 (RYA-492)',
        'owed': (f"REPOINTED to the RYA-492 Costa-Silva-2020 atlas-tuned S I gf (single "
                 f"cited source; canonical_gf S I 6743 -0.6103->-0.5476). Synthesis A(S)_NLTE "
                 f"{a:.3f} (σ {cs.get('sigma')}; GES-gf control {ctrl}) — the CS gf moved it "
                 f"only ~-0.03 (7.516->{a:.3f}). {d:+.3f} vs Asplund {asp:.2f} is a gf-SCALE "
                 f"floor (RYA-161), NOT a line-ID error and NOT closed by the gf; stays "
                 f"CURATION-OWED, do NOT tune. PROVENANCE of the other S numbers: the RYA-527 "
                 f"two-engine Engine-A 7.369 is the EW path — the EW cull keeps only S I "
                 f"6757.15 (A_LTE 7.386) + NLTE delta -0.017 = 7.369 (single blend-limited "
                 f"line, no 6743, does NOT use the CS gf); the frozen gold v1 7.753 is the "
                 f"older EW cull (n=2). These are distinct channels, now reconciled: the "
                 f"reported verdict value is the CS-gf synthesis {a:.3f}.")}}


def _apply_s_synthesis(rows, data):
    """Overlay the RYA-492 Costa-Silva-gf S synthesis reclassification onto the rows (RYA-557)."""
    overrides = _s_reclassify(data)
    for r in rows:
        ov = overrides.get(r['element'])
        if not ov:
            continue
        for key in ('verdict', 'channel', 'owed', 'provenance'):
            if key in ov:
                r[key] = ov[key]
        if ov.get('A_measured') is not None:
            r['A_measured'] = round(float(ov['A_measured']), 3)
            r['delta_vs_asplund'] = round(r['A_measured'] - r['asplund2021'], 3)
        if ov.get('sigma') is not None:
            r['sigma'] = round(float(ov['sigma']), 3)
        if ov.get('n_lines') is not None:
            r['n_lines'] = int(ov['n_lines'])
    return overrides


def _load_ba_synthesis():
    """The Ba II 5853 measurement, or None if it hasn't run.

    RYA-581 SUPERSEDES RYA-559: prefer the in-window deblend (profile fit with the
    blends modelled) over the EW->COG inversion whenever it is present. The 559 record
    stays as the fallback so the fold still works on a checkout that predates the
    deblend. Both are routed through _ba_reclassify, which branches on `ticket`."""
    if BA_DEBLEND_JSON.exists():
        return json.loads(BA_DEBLEND_JSON.read_text())
    if not BA_SYNTH_JSON.exists():
        return None
    return json.loads(BA_SYNTH_JSON.read_text())


def _ba_reclassify_deblend(data):
    """RYA-581 — fold the Ba II 5853 IN-WINDOW DEBLEND into the verdict, superseding the
    RYA-559 EW->COG value.

    RYA-559 landed A(Ba)_NLTE 2.410 by inverting the OBSERVED pool EW (74.62 mA) through
    an HFS-resolved curve of growth. That EW carries blend_flag=True — ~10 mA over the
    clean solar Ba II 5853 (~64 mA) — and an EW inversion cannot deblend: one scalar
    cannot separate barium from the rest of the absorption in the integration window, so
    the neighbours were charged to Ba. RYA-559 said so and routed the debt here.

    RYA-581 fits the PROFILE instead (RYA-551 Sr II pattern): the full VALD3 in-window
    block is synthesised alongside the Ba II HFS/isotope components and A(Ba) is fitted
    by chi2, with the Engine-A Korotin2015 delta read at the solar node as before.

    Validate-don't-tune: the profile, the canonical gf and the Korotin delta are all
    READ; A is whatever chi2 returns. Nothing is fitted toward Asplund 2.27. The verdict
    follows the measured reconciliation — if the deblended A lands inside TOL_PASS then
    Ba PASSes, and RYA-581 says explicitly to report that honestly rather than force-hold
    it owed. The single-line caveat rides in the note, flagged for the RYA-527 freeze."""
    a = float(data['A_nlte'])
    asp = float(data.get('asplund2021', SOLAR_ASPLUND2021.get('Ba', 2.27)))
    d = a - asp
    dk = data.get('engineA_korotin_delta')
    bm = data.get('blend_model', {}) or {}
    ev = data.get('deblend_evidence', {}) or {}
    cb = data.get('correction_budget_dex', {}) or {}
    hfs = data.get('hfs', {}) or {}
    top = ', '.join(f"{k} {v} mA" for k, v in
                    list((bm.get('per_species_core_EW_mA') or {}).items())[:4])
    reconciled = abs(d) <= TOL_PASS
    # A profile fit that railed, ran away in chi2, or lost its sensitivity to A is not
    # allowed to buy a PASS on arithmetic alone.
    reliable = bool(data.get('reliable'))
    verdict = 'PASS' if (reconciled and reliable) else 'CURATION-OWED'
    return {'Ba': {
        'verdict': verdict,
        'A_measured': a, 'sigma': data.get('sigma'),
        'n_lines': 1, 'provenance': 'synthesis: Ba II 5853 in-window deblend (RYA-581)',
        'channel': 'synthesis: Ba II 5853.668 in-window blend fit (Turbospectrum, HFS + '
                   'full VALD3 in-window block, chi2 profile fit) + Engine-A Korotin2015 '
                   '1D-NLTE delta',
        'owed': (
            f"MEASURED by in-window blend fit — SUPERSEDES the RYA-559 EW->COG value 2.410. "
            f"RYA-559 inverted the blend_flag=True pool EW (74.62 mA vs the clean solar line "
            f"~64 mA); an EW inversion cannot deblend, so the neighbours were charged to Ba. "
            f"RYA-581 synthesises the full VALD3 in-window block ({bm.get('n_rows')} rows, "
            f"{bm.get('n_species')} species over {bm.get('window_A')} A) alongside the "
            f"{hfs.get('n_components')} Ba II HFS/isotope components and fits A(Ba) to the "
            f"observed profile by chi2 (RYA-551 pattern). Modelled blend in the core "
            f"+/-{bm.get('core_hw_A')} A = {bm.get('blend_core_EW_mA')} mA, dominated by "
            f"{top}. A(Ba)_LTE {data.get('A_lte')} + Engine-A Korotin2015 delta "
            f"{dk:+.4f} = A(Ba) {a:.3f} ({d:+.3f} vs Asplund {asp:.2f}). Deblend is "
            f"demonstrated, not assumed: red_chi2 {ev.get('red_chi2_ba_alone')} (Ba alone) "
            f"-> {ev.get('red_chi2_blends_modelled')} (blends modelled), "
            f"{ev.get('chi2_improvement_factor')}x better. Correction budget vs 2.410: "
            f"{cb.get('from_dropping_the_EW_inversion')} dex from abandoning the EW "
            f"inversion + {cb.get('from_modelling_the_in_window_blend')} dex from the "
            f"in-window blend model; the fitted synthetic core EW "
            f"{cb.get('fitted_core_EW_mA')} mA reproduces the RYA-559 calibration "
            f"(A=2.27 -> 66.5 mA) and the literature clean line (~64-66 mA), confirming the "
            f"CLEAN Ba II 5853 was recovered. Validate-don't-tune: profile, canonical gf and "
            f"Korotin delta all READ, never fitted toward Asplund; do NOT tune. "
            f"sigma {data.get('sigma')} is the HARPS-vs-IAG arm scatter — a precision "
            f"floor, NOT a total error budget (no gf / Korotin-delta / 1D-MARCS / "
            f"single-line terms). "
            + (f"RECONCILED within TOL_PASS {TOL_PASS} -> PASS, but it rests on ONE line: "
               f"flagged for the RYA-527 freeze review, and a second clean Ba II line "
               f"(6141.713 / 6496.897) is the confirming follow-up. Note the two-engine "
               f"record is still UNEVALUABLE for Ba — the synthesis route is not wired into "
               f"_dedicated_engine_B(), so gate 3 cannot see this value (RYA-669)."
               if verdict == 'PASS' else
               f"NOT PASS: |{d:+.3f}| exceeds TOL_PASS {TOL_PASS}"
               + ("" if reliable else " / the profile fit did not meet the reliability floor")
               + "."))}}


def _ba_reclassify(data):
    """RYA-559 — fold the Ba II 5853 synthesis measurement into the verdict. Ba was
    owed-NO-value: the EW path SAT-culls the strong Ba II 5853 line (REW -4.90, at the
    saturation knee, +HFS/isotope), so it never enters the produced pool. The registry
    route is SYNTHESIS of Ba II + the wired Korotin2015 NLTE grid (Ba is the majority
    ion). This lands the value: an HFS-resolved LTE curve-of-growth (Turbospectrum bsyn,
    23 VALD3 components, total loggf -1.000) inverts the OBSERVED solar EW (74.62 mA) ->
    A(Ba)_LTE, then the production Engine-A Korotin2015 1D-NLTE delta (solar 5853 node,
    -0.0285) is applied. Driven by the measured value, not tuned: a synthesised Ba that
    sits high is a finding, not a PASS.

    Ba stays CURATION-OWED: the +0.14 vs Asplund 2.27 is driven by the blend-inflated
    pool EW (blend_flag=True; ~10 mA over the clean solar Ba II 5853 ~64 mA) — a gf/blend
    floor routed to RYA-161/162, NOT closed here. Validate-don't-tune: the EW and the
    Korotin delta are READ, never fitted toward Asplund. Returns {'Ba': override-dict}."""
    if not data:
        return {}
    a = data.get('A_nlte')
    if a is None:
        return {}
    if data.get('ticket') == 'RYA-581':
        return _ba_reclassify_deblend(data)
    a = float(a)
    asp = float(data.get('asplund2021', SOLAR_ASPLUND2021.get('Ba', 2.27)))
    d = a - asp
    dk = data.get('engineA_korotin_delta')
    a_lte = data.get('A_lte')
    clean = (data.get('clean_ew_crosscheck', {}) or {}).get('values', {})
    clean_txt = '; '.join(f"{ew} mA->{v.get('A_nlte')}" for ew, v in clean.items())
    reconciled = abs(d) <= TOL_PASS          # False for Ba (blend/gf floor) — never PASS
    return {'Ba': {
        'verdict': 'PASS' if reconciled else 'CURATION-OWED',
        'A_measured': a, 'sigma': data.get('sigma_ew_dex'),
        'n_lines': 1, 'provenance': 'synthesis: Ba II 5853 HFS (RYA-559)',
        'channel': 'synthesis: Ba II 5853.668 HFS-resolved LTE COG (Turbospectrum, '
                   'VALD3 gf loggf -1.000) + Engine-A Korotin2015 1D-NLTE delta',
        'owed': (f"MEASURED via HFS-resolved synthesis — the EW path SAT-culls the strong "
                 f"Ba II 5853 (REW -4.90, at the knee, +HFS/isotope; registry owed). "
                 f"A(Ba)_LTE {a_lte} (bsyn COG inverting the observed solar EW "
                 f"{data.get('ew_obs_mA')} mA, 23 VALD3 HFS components, total loggf "
                 f"{data.get('total_loggf')}) + Engine-A Korotin2015 delta {dk:+.4f} "
                 f"(solar 5853 node; validate-don't-tune) = A(Ba) {a:.3f} ({d:+.3f} vs "
                 f"Asplund {asp:.2f}). Off no-value. Engine-B (Gerber atom.ba111 TS-native "
                 f"NLTE) corroborates at the departure level (Ba II 4554 delta -0.018, "
                 f"RYA-534) — 5853 is absent from the GES level-ID block so Engine-A drives "
                 f"the delta; both agree Ba II NLTE is small-negative, so A is engine-"
                 f"insensitive. NOT PASS: the pool EW carries blend_flag=True (~10 mA over "
                 f"the clean solar line ~64 mA), inflating A by ~+0.15 — a blend/gf floor "
                 f"(RYA-161/162), NOT real Ba enhancement. Clean-EW cross-check ({clean_txt}) "
                 f"reconciles with Asplund, confirming the offset is the blend. Deblend owed; "
                 f"do NOT tune.")}}


def _apply_ba_synthesis(rows, data):
    """Overlay the RYA-559 Ba II 5853 synthesis reclassification onto the rows (in place)."""
    overrides = _ba_reclassify(data)
    for r in rows:
        ov = overrides.get(r['element'])
        if not ov:
            continue
        for key in ('verdict', 'channel', 'owed', 'provenance'):
            if key in ov:
                r[key] = ov[key]
        if ov.get('A_measured') is not None:
            r['A_measured'] = round(float(ov['A_measured']), 3)
            r['delta_vs_asplund'] = round(r['A_measured'] - r['asplund2021'], 3)
        if ov.get('sigma') is not None:
            r['sigma'] = round(float(ov['sigma']), 3)
        if ov.get('n_lines') is not None:
            r['n_lines'] = int(ov['n_lines'])
    return overrides


def _load_co_synthesis():
    """RYA-564: the red-line Co I HFS-synthesis measurement, or None if it hasn't run."""
    if not CO_SYNTH_JSON.exists():
        return None
    return json.loads(CO_SYNTH_JSON.read_text())


def _co_reclassify(data):
    """RYA-564 — fold the RED-line Co I measurement into the verdict, REPLACING the demoted
    blue-edge 3845 artifact (+1.188, KP SNR~24, chi2r~3100; see _kittpeak_reclassify).

    Co is hyperfine-split, so the EW path cannot reach it (the RYA-354/466 finding). The
    measurement is an HFS-resolved Turbospectrum flux fit on clean red Co I lines — the GES
    NLTE line list supplies the HFS components, the VALD in-window block supplies the blends
    — plus a PER-LINE 1D-NLTE delta read from the RYA-534-validated Gerber grid (Engine-B
    TS-native; the deck raises rather than running silent-LTE). Validate-don't-tune: both the
    gf (canonical_gf single source) and the NLTE delta are READ, never fitted.

    If NO red line clears the reliability floor the element reports NO VALUE — the ticket's
    CRITICAL condition. Falling back to 3845 is never an option. Returns {'Co': override}."""
    if not data:
        return {}
    s = data.get('_summary', {})
    a = s.get('A_Co')
    med_delta = s.get('median_nlte_delta')
    anchor_ok = s.get('nlte_anchor_consistent')
    asp = float(SOLAR_ASPLUND2021.get('Co', 4.94))
    nlte_txt = (f"per-line 1D-NLTE from the RYA-534-validated Gerber TS-native grid "
                f"(median {med_delta:+.4f} vs the Bergemann+2010 anchor +0.100 ± 0.12 — "
                f"{'CONSISTENT' if anchor_ok else 'CHECK'}; read, never fitted)")

    if a is None:
        # Reliability floor not cleared anywhere -> honest no-value. NEVER the 3845 artifact.
        return {'Co': {
            'verdict': 'CURATION-OWED', 'A_measured_blank': True, 'n_lines': 0,
            'provenance': 'synthesis: Co I red HFS (RYA-564) — no line cleared reliability',
            'channel': 'synthesis: Co I red HFS-resolved (RYA-564) — measurable-owed, no value',
            'owed': (f"NO VALUE. The red Co I HFS synthesis ran ({nlte_txt}) but no line cleared "
                     f"the reliability floor: {s.get('reason')}. Co is measurable-owed. The "
                     f"blue-edge 3845 artifact stays DIAGNOSTIC-ONLY — it is not a fallback.")}}

    a = float(a)
    d = a - asp
    n = int(s.get('n_reliable', 0))
    lines = s.get('reliable_lines', {}) or {}
    line_txt = ', '.join(f"{w}={v:.3f}" for w, v in sorted(lines.items()))
    po = s.get('primary_only') or {}
    iag = s.get('iag_crosscheck') or {}
    reconciled = abs(d) <= TOL_PASS
    checks = []
    if po:
        checks.append(f"gf-'agreed' lines only: {po['median']:.3f} (n={po['n']}, "
                      f"scatter {po['scatter']:.3f})")
    if iag:
        checks.append(f"IAG FTS arm: {iag['median']:.3f} (n={iag['n']}, "
                      f"scatter {iag['scatter']:.3f})")
    check_txt = ('; '.join(checks)) if checks else 'none'
    return {'Co': {
        'verdict': 'PASS' if reconciled else 'CURATION-OWED',
        'A_measured': a, 'sigma': s.get('scatter'), 'n_lines': n,
        'provenance': 'synthesis: Co I red HFS (RYA-564)',
        'channel': (f'synthesis: Co I red HFS-resolved flux fit ({n} lines, HARPS; blends '
                    f'modelled) + per-line Gerber 1D-NLTE (RYA-534 deck)'),
        'owed': (f"MEASURED on clean RED Co I lines via HFS-resolved synthesis — the EW path "
                 f"cannot reach Co (hyperfine-split, RYA-354/466) and the previous number came "
                 f"from the blue-edge Co I 3845, now DEMOTED to diagnostic-only (that +1.188 was "
                 f"an SNR~24/chi2r~3100 artifact, not solar Co). Per-line A(Co) [{line_txt}] -> "
                 f"A(Co) {a:.3f} ({d:+.3f} vs Asplund {asp:.2f}), scatter {s.get('scatter')}; "
                 f"{nlte_txt}. Cross-checks — {check_txt}. "
                 + ("RECONCILES within tol on a validated leg: the red-line value replaces the "
                    "artifact and Co clears. Scale caveat (RYA-561/593): this is a 1D-NLTE value "
                    "against a 3D-NLTE Asplund reference — the un-applied 3D term is folded into "
                    "the offset; the class-wide 3D-metals correction is post-Beta (RYA-593)."
                    if reconciled else
                    "A residual survives tol -> curation owed (gf/blend floor, RYA-161/162); "
                    "do NOT tune."))}}


def _apply_co_synthesis(rows, data):
    """Overlay the RYA-564 red-line Co synthesis reclassification onto the rows (in place).
    Runs AFTER _apply_kittpeak so it replaces the demoted 3845 row."""
    overrides = _co_reclassify(data)
    for r in rows:
        ov = overrides.get(r['element'])
        if not ov:
            continue
        for key in ('verdict', 'channel', 'owed', 'provenance'):
            if key in ov:
                r[key] = ov[key]
        if ov.get('A_measured_blank'):
            r['A_measured'] = None
            r['delta_vs_asplund'] = None
            r['sigma'] = None
        if ov.get('A_measured') is not None:
            r['A_measured'] = round(float(ov['A_measured']), 3)
            r['delta_vs_asplund'] = round(r['A_measured'] - r['asplund2021'], 3)
        if ov.get('sigma') is not None:
            r['sigma'] = round(float(ov['sigma']), 3)
        if ov.get('n_lines') is not None:
            r['n_lines'] = int(ov['n_lines'])
    return overrides


def _apply_cu_v_synthesis(rows, data):
    """Overlay the RYA-466 Cu/V HFS-synthesis reclassification onto the base rows (in place)."""
    overrides = _cu_v_reclassify(data)
    for r in rows:
        ov = overrides.get(r['element'])
        if not ov:
            continue
        for key in ('verdict', 'channel', 'owed', 'provenance'):
            if key in ov:
                r[key] = ov[key]
        if ov.get('A_measured') is not None:
            r['A_measured'] = round(float(ov['A_measured']), 3)
            r['delta_vs_asplund'] = round(r['A_measured'] - r['asplund2021'], 3)
        if ov.get('sigma') is not None:
            r['sigma'] = round(float(ov['sigma']), 3)
        if ov.get('n_lines') is not None:
            r['n_lines'] = int(ov['n_lines'])
    return overrides


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
        # RYA-564: an explicit demotion — the channel measured something but it is NOT
        # reportable (Co I 3845). Blank the value rather than leaving a stale one standing;
        # a demoted artifact must never survive as the element's number.
        if ov.get('A_measured_blank'):
            r['A_measured'] = None
            r['delta_vs_asplund'] = None
            r['sigma'] = None
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
    # RYA-674: the gold source is a NAMED INPUT, not a hardcoded 'CURRENT'. Default is
    # unchanged, so an unflagged run behaves exactly as before; naming a frozen version
    # is how the verdict is regenerated while CURRENT carries a row that the RYA-681
    # scale guard refuses to load. The resolved version is stamped into the summary.
    ap.add_argument('--gold-version', default='CURRENT',
                    help="frozen gold solar reference to compute the verdict against "
                         "('CURRENT' or a 'vN' token, e.g. v2). RYA-469/674.")
    args = ap.parse_args()

    ab, ew, phase_a, gold_version = _load(args.gold_version)
    rows = build_verdicts(ab, ew, phase_a)

    # RYA-460: fold in the Kitt Peak measurements (N + P/K/Co/Sc) if the campaign ran.
    kp = _load_kittpeak()
    overrides = _apply_kittpeak(rows, kp)

    # RYA-466: fold in the HFS-resolved synthesis measurement of Cu / V (the elements the
    # EW path cannot reach — hyperfine-split). Moves Cu/V off "no value" → measured.
    cuv = _load_cu_v_synthesis()
    cuv_overrides = _apply_cu_v_synthesis(rows, cuv)

    # RYA-473: fold in the HFS-resolved synthesis measurement of Mn (the Den Hartog e6S→z6P
    # triplet the EW path SAT-culls). Moves Mn off "no value" → measured.
    mn = _load_mn_synthesis()
    mn_overrides = _apply_mn_synthesis(rows, mn)

    # RYA-557: repoint S to the RYA-492 Costa-Silva-gf synthesis value (7.486). The EW path
    # leaves S owed-no-value (all but 6757.15 culled; the CS-gf 6743 line is synthesis-only).
    s_synth = _load_s_synthesis()
    s_overrides = _apply_s_synthesis(rows, s_synth)

    # RYA-559: fold in the Ba II 5853 synthesis value (2.41). The EW path SAT-culls the
    # strong Ba II 5853 (REW -4.90, HFS/blend); synthesis + Korotin NLTE lands it (owed).
    ba_synth = _load_ba_synthesis()
    ba_overrides = _apply_ba_synthesis(rows, ba_synth)

    # RYA-564: fold in the RED-line Co I HFS synthesis. MUST run after _apply_kittpeak — it
    # replaces the blue-edge Co I 3845 artifact (+1.188) that KP demoted to diagnostic-only.
    co_synth = _load_co_synthesis()
    co_overrides = _apply_co_synthesis(rows, co_synth)

    counts = {}
    for r in rows:
        counts[r['verdict']] = counts.get(r['verdict'], 0) + 1
    diff = {k: counts.get(k, 0) - PRIOR_COUNTS.get(k, 0)
            for k in ('PASS', 'NLTE-OWED', 'CURATION-OWED', 'DATA-GAP')}
    summary = {'ticket': 'RYA-371 Phase C (RYA-462 NLTE-grid-wired: K)', 'star': args.star,
               'generated': date.today().isoformat(),
               'reference': 'Asplund, Amarsi & Grevesse 2021 (A&A 653, A141)',
               # RYA-674: the gold version this run ACTUALLY read, resolved from the
               # named input — not `current_version()`, which is a second lookup that a
               # freeze between load and write would silently desynchronise. Metadata
               # about a number must come from the number's own provenance.
               'solar_ref_version': gold_version,          # RYA-469 gold reference frozen
               'tol_pass_dex': TOL_PASS, 'n_elements': len(rows), 'counts': counts,
               'prior_counts': PRIOR_COUNTS, 'diff_vs_prior': diff,
               'kittpeak_wired': bool(kp),
               'kittpeak_elements': sorted(overrides) if overrides else [],
               'cu_v_synthesis_wired': bool(cuv),
               'cu_v_synthesis_elements': sorted(cuv_overrides) if cuv_overrides else [],
               'mn_synthesis_wired': bool(mn),
               'mn_synthesis_elements': sorted(mn_overrides) if mn_overrides else [],
               's_synthesis_wired': bool(s_synth),
               's_synthesis_elements': sorted(s_overrides) if s_overrides else [],
               'ba_synthesis_wired': bool(ba_synth),
               'ba_synthesis_elements': sorted(ba_overrides) if ba_overrides else [],
               'co_synthesis_wired': bool(co_synth),
               'co_synthesis_elements': sorted(co_overrides) if co_overrides else [],
               # RYA-564: the blue-edge Co I 3845 is demoted unconditionally, whether or not
               # the red-line synthesis has run — record it so the freeze can assert it.
               'demoted_diagnostic_only': ['Co I 3845 (blue-edge artifact, RYA-564)']}

    # RYA-674 §2C: re-gate AFTER every dedicated-channel overlay. The overlays above
    # replace reported values (Kitt Peak, Cu/V, Mn, S, Ba, Co), and an overlay that
    # substitutes a vetoed or excluded value is exactly the RYA-669 shape.
    assert_ratified_constraints_satisfied(
        rows, 'phase_c verdict generator (RYA-371, post-overlay)')

    AUDIT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    payload = {'summary': summary, 'verdicts': rows}
    out_json = AUDIT / 'solar_phase_c_verdict.json'   # cno_synthesis audit trail (star-prefixed)
    with open(out_json, 'w') as fh:
        json.dump(payload, fh, indent=2)
    # RYA-469: the canonical per-star verdict is namespaced (data/outputs/{star}/{star}_verdict.json)
    out_ns = ns.output_path(args.star, 'verdict.json')
    with open(out_ns, 'w') as fh:
        json.dump(payload, fh, indent=2)
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
