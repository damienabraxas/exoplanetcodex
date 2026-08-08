#!/usr/bin/env python3
"""
scripts/rya527_two_engine_run.py
================================
RYA-527 REAL two-engine solar run (RYA-525 floor) — NOT an overlay.

Drives BOTH engines over real solar data per line and calls the RYA-525 selector:
  Engine A = 1D-NLTE  = EW->A(X) (iSpec/MOOG, absolute) + production NLTE grid delta
                        (nlte_corrections / nlte_cno), per line.
  Engine B = synthesis = Turbospectrum flux-fit (abundances_derive synthesis-v2),
                        per line; + Gerber TS-native NLTE delta for the 11 Family-A
                        grids that PREFLIGHT reconciled on Sirius (RYA-534).
Per line -> LineEngines -> engine_selection.select_element -> ElementRecord
(reported value, engineA/engineB aggregates, cross-engine delta, mix flag).

LOUD-FAIL (RYA-525): a synthesis-required element (problem_children) with NO
Engine-B value RAISES — never a silent single-engine fall-back.

Engine-B sources, labelled per element:
  - fresh synth-v2 per-line (this run's data/outputs/solar/solar_per_line_synth_v2.csv);
  - for synthesis-required HFS/Sr elements SUPPRESSED from the EW pool (so absent
    from synth-v2), the dedicated Engine-B synthesis measurement (Mn RYA-473,
    Cu/V RYA-466, Sr II RYA-551) — an Engine-B output, injected as a single synth
    line, clearly sourced. A fully-fresh re-run of those harnesses is the Sirius
    step; here they carry their committed synthesis value.

Preflight (scripts/rya527_preflight_reconciliation.py) MUST be green first.
"""
import argparse
import json
import sys
import warnings
from pathlib import Path

warnings.filterwarnings('ignore')
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))   # repo root before pipeline

from pipeline import _runtime as _rt  # noqa: F401,E402  BLAS/fork pins before numpy
import numpy as np
import pandas as pd

import pipeline.abundances_derive as ad
from pipeline import nlte_corrections as nc
from pipeline.engine_selection import (LineEngines, select_element, ENGINE_A, ENGINE_B,
                                       TwoEngineError, RATIFIED_EXCLUDED_SPECIES,
                                       exclusion_reason, is_upper_limit_disposition)
from pipeline.ratified_constraints import (  # RYA-674 emission-time gate
    RowKind, assert_ratified_constraints_satisfied)
from config.constants import (get_star_params, TARGET_ELEMENTS, NLTE_CORRECTION_ELEMENTS,
                              SOLAR_ASPLUND2021)
import pipeline.problem_children as pc
from pipeline import two_engine_inputs as tei

ROOT = Path(__file__).resolve().parent.parent
# RYA-682: the canonical (RYA-469 namespaced) location comes from data_namespace via
# two_engine_inputs — never hand-built here, so the driver and the generator cannot
# drift apart. data/outputs/ is gitignored: this input is GENERATED, not committed.
ENGINE_B_PL = tei.engine_b_per_line_path('solar')
OUT_DIR = ROOT / 'data' / 'audit' / 'rya527_two_engine'

# Gerber TS-native NLTE delta per element (Engine-B NLTE), RYA-534 anchor-validated;
# preflight-reconciled grids on Sirius. Ti ships atom.ti503b (xfail RYA-548) -> its
# Engine-B NLTE is a cross-engine DIAGNOSTIC only, not applied to the reported value.
GERBER_NLTE_DELTA = {'O': -0.105, 'Mg': -0.023, 'Si': -0.034, 'Ca': -0.009, 'Ni': +0.018,
                     'Na': -0.068, 'Co': +0.099, 'Ba': -0.018, 'Mn': +0.043, 'Sr': -0.013,
                     'Ti': +0.221}
GERBER_XFAIL = {'Ti'}   # RYA-548

# Dedicated Engine-B synthesis measurements for the synthesis-required elements
# (CNO via nlte_cno; HFS metals + Sr II via the HFS/Sr synthesis harnesses). These
# are Turbospectrum synthesis-harness outputs (Engine-B), NOT verdict re-labels.
CNO_PHASE_A = ROOT / 'data' / 'audit' / 'cno_synthesis' / 'solar_phase_a_cross_arm.json'
MN_JSON  = ROOT / 'data' / 'audit' / 'mn_hfs_synthesis' / 'solar_mn_hfs_synthesis_rya473.json'
CUV_JSON = ROOT / 'data' / 'audit' / 'cu_v_hfs_synthesis' / 'solar_cu_v_hfs_synthesis_rya466.json'
SR2_JSON = ROOT / 'data' / 'results' / 'sr2_synthesis_rya551.json'
ZR2_JSON = ROOT / 'data' / 'results' / 'zr2_synthesis_rya560.json'   # RYA-560 Zr II LTE synth
ZR2_DEBLEND_JSON = ROOT / 'data' / 'results' / 'zr2_deblend_rya585.json'  # RYA-585 deblend refit
MG5528_JSON = ROOT / 'data' / 'results' / 'mg_5528_synthesis_rya592.json'  # RYA-592 Mg 2nd line


def _solar_params():
    r = get_star_params('solar')
    return {'teff_K': float(r['teff']), 'logg': float(r['logg']), 'feh': 0.0,
            'vturb_kms': float(r.get('xi', 1.0))}


def _nlte_delta_A(el, wave, p):
    """Production Engine-A NLTE grid delta for one line. (delta, a_in_hull).
    No grid -> LTE (delta 0, in-hull True). Out of hull -> (nan, False)."""
    if el not in NLTE_CORRECTION_ELEMENTS:
        return 0.0, True
    if not nc.element_grid_in_bounds(el, p['teff_K'], p['logg'], p['feh']):
        return float('nan'), False
    d = nc._mpia_element_delta(el, wave, p['teff_K'], p['logg'], p['feh'])
    return (float(d) if d is not None and np.isfinite(d) else 0.0), True


def _engine_A_perline(p):
    """Fresh Engine-A (1D-NLTE) per line for every EW-pool element. spec_abund is
    absolute A(X) (hydrogen=12), same scale as Engine-B a_synth."""
    ew_df = ad._load_solar_ews()
    atm = ad._load_atmosphere(p['teff_K'], p['logg'], p['feh'], p['vturb_kms'])
    lm, a_abs, _, _ = ad._ew_to_abundance(ew_df, p, atm)
    out = {}
    for i in range(len(lm)):
        tok = str(lm['note'][i]).split()
        if len(tok) < 2:
            continue
        el = tok[0]
        ion = 'I' if tok[1] == '1' else ('II' if tok[1] == '2' else tok[1])
        a = float(a_abs[i])
        if not np.isfinite(a):
            continue
        wave = float(lm['wave_A'][i])
        d, in_hull = _nlte_delta_A(el, wave, p)
        a_nlte = a + d if np.isfinite(d) else a
        out.setdefault((el, ion), {})[round(wave, 1)] = dict(
            wave=wave, ew=float(lm['ew'][i]), a=a_nlte, in_hull=in_hull)
    return out


def _engine_B_perline():
    """Fresh Engine-B (synthesis-v2) per line + Gerber TS-native NLTE delta.

    The artifact is validated by the preflight (RYA-682) before any compute runs;
    this re-assert keeps the function safe to call on its own.
    """
    tei.assert_engine_b_artifact('solar')
    df = pd.read_csv(ENGINE_B_PL)
    out = {}
    for _, r in df.iterrows():
        el, ion = str(r['element']), str(r['ion'])
        a = r.get('a_synth')
        if a is None or not np.isfinite(a) or str(r.get('status')) != 'ok':
            continue
        b = float(a)
        if el in GERBER_NLTE_DELTA and el not in GERBER_XFAIL:
            b += GERBER_NLTE_DELTA[el]     # Engine-B TS-native NLTE (reconciled grid)
        # b_chi2=None on purpose: synth-v2 red_chi2 uses a 0.01 model-adequacy floor
        # (median ~105), NOT comparable to the two-engine synth_chi2_gate (10) which
        # expects a noise-normalised ~1. The 'status'==ok filter above is the
        # catastrophic-failure gate; regime routing does the quality selection.
        out.setdefault((el, ion), {})[round(float(r['wavelength_air_A']), 1)] = dict(
            wave=float(r['wavelength_air_A']), b=b,
            gerber=(el in GERBER_NLTE_DELTA and el not in GERBER_XFAIL))
    return out


def _dedicated_engine_B():
    """Committed Engine-B SYNTHESIS-HARNESS measurements for the synthesis-required
    elements (CNO nlte_cno primary indicator; Mn/Cu/V HFS synth; Sr II synth)."""
    out = {}
    if CNO_PHASE_A.exists():
        ca = json.loads(CNO_PHASE_A.read_text()).get('cross_arm', {})
        for el in ('C', 'O'):
            prim = next((i for i in ca.get(el, {}).get('indicators', [])
                         if i.get('role') == 'primary'), None)
            if prim and prim.get('A') is not None:
                out[(el, 'I')] = (float(prim['A']), f"nlte_cno synthesis {prim.get('key')} (RYA-491/237)")
    if MN_JSON.exists():
        m = json.loads(MN_JSON.read_text()).get('Mn', {})
        v = m.get('A_nlte') or m.get('A_lte_median')
        if v is not None:
            out[('Mn', 'I')] = (float(v), 'RYA-473 HFS synth')
    if CUV_JSON.exists():
        d = json.loads(CUV_JSON.read_text())
        for el in ('Cu', 'V'):
            e = d.get(el, {})
            v = e.get('A_nlte') or e.get('A_lte_median')
            if v is not None:
                out[(el, 'I')] = (float(v), 'RYA-466 HFS synth')
    if SR2_JSON.exists():
        v = json.loads(SR2_JSON.read_text()).get('4077.709', {}).get('harps', {}).get('A_NLTE')
        if v is not None:
            out[('Sr', 'II')] = (float(v), 'RYA-551 Sr II synth')
    # Zr II — the majority ion -> LTE-robust (registry 279/458, the Sr II/V II
    # precedent); A_LTE IS the value, no NLTE grid. RELIABILITY-GATED throughout:
    # emit only a line that cleared the dEW/dA floor and is not railed.
    #
    # Two sources, tried best-first. RYA-585 (deblend) supersedes RYA-560 for the
    # three strong lines because it re-fits the SAME syntheses with the blends
    # modelled in-window and a blend-pixel continuum, and additionally gates on a
    # sane red_chi2. RYA-560 remains the fallback so the original measurement stays
    # reproducible and wired if the deblend artifact is absent.
    #
    # As of the RYA-585 Sirius run BOTH are silent and Zr stays MEASURABLE-OWED.
    # The deblend fixed what it set out to fix — red_chi2 collapsed from 41-91 to
    # <=1.7, confirming the blend/continuum systematic was real — but every line
    # still sits below the sensitivity floor (best dEW/dA 36.5 < 40) because these
    # three cores are saturated (sat_index 0.36-0.69). That is an intrinsic property
    # of the line set, not a modelling defect, so refitting cannot rescue it; the
    # next lever is cleaner blue Zr II lines (RYA-458). Never a silent sub-floor
    # value. When a reliable Zr II line lands, it flows through here unchanged.
    for _src, _path, _tag in ((585, ZR2_DEBLEND_JSON, 'RYA-585 Zr II deblend LTE'),
                              (560, ZR2_JSON, 'RYA-560 Zr II synth LTE')):
        if not _path.exists():
            continue
        zr = json.loads(_path.read_text())
        rel = [d['harps']['A_LTE'] for w, d in zr.items()
               if isinstance(d, dict) and isinstance(d.get('harps'), dict)
               and d['harps'].get('reliable') and d['harps'].get('A_LTE') is not None]
        if rel:
            out[('Zr', 'II')] = (float(np.mean(rel)),
                                 f"{_tag} (n={len(rel)} reliable)")
            break
    if MG5528_JSON.exists():
        # RYA-592: the SECOND clean Mg I line (5528.405), measured by in-window blend-fit
        # synthesis so Mg could stop being single-line. CONCORDANCE-GATED, and as of the
        # RYA-592 Sirius run the gate is CLOSED: 5528 is reliable (dEW/dA 130 mA/dex,
        # red_chi2 1.4, 288 blend components modelled) but lands 0.21-0.23 dex BELOW the
        # same harness's 5711 -- outside the 0.10 band. Emitting it would silently average
        # two measurements that disagree, which is exactly what the RYA-525 floor forbids
        # ("never silently average two disagreeing scales"), and it would move Mg's reported
        # value on evidence that is itself contested. So while DISCORDANT this contributes
        # NOTHING and Mg stays single-line CURATION-OWED with the reason recorded (the
        # RYA-560 Zr pattern: wired, gated, currently silent). When the discordance is
        # adjudicated (element_status_tracker_drift.md section E), the line flows through
        # here without further wiring. Note this can only ever add an ENGINE-B line: 5528's
        # EW is 3.4x the ratified saturation knee, so it has no Engine-A EW route and
        # cannot create the dCE that RYA-561 gate 3 requires.
        mg = json.loads(MG5528_JSON.read_text())
        v = mg.get('_verdict', {})
        if v.get('lines_concordant') and v.get('target_reliable'):
            out[('Mg', 'I')] = (float(v['target_A_NLTE_engineB']),
                                'RYA-592 Mg I 5528 in-window blend-fit synth (concordant)')
        else:
            print(f"[two-engine] RYA-592 Mg I 5528 HELD OUT: reliable="
                  f"{v.get('target_reliable')}, concordant={v.get('lines_concordant')} "
                  f"(|d| {v.get('concordance_worst_abs_dex')} vs band "
                  f"{v.get('concordance_band')}) -> Mg stays single-line")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--star', default='solar')
    # RYA-669: Phase 2 re-emits into its OWN directory so the 2026-07-18 record stays
    # readable beside the fresh one. Overwriting it would destroy the only evidence of
    # what the pre-v3 floor actually produced, which is what the diff is read against.
    ap.add_argument('--out-dir', default=None,
                    help='repo-relative output dir (default data/audit/rya527_two_engine)')
    args = ap.parse_args()

    # ── RYA-682 preflight: every input, checked BEFORE any compute ───────────
    # The Engine-A leg below costs minutes (GES linelist load, EW triage, MOOG
    # baseline). Discovering a missing input after paying for that is a defect in
    # its own right, and an input that is present-but-unusable never surfaced at
    # all. Both are settled here, first, in one place.
    if args.star != 'solar':
        raise SystemExit(
            f"RYA-682: --star={args.star!r} is not supported. Every dedicated Engine-B "
            f"input this driver reads is solar-specific (CNO cross-arm, Mn/Cu/V HFS, "
            f"Sr II, Zr II, Mg 5528), and the Engine-B per-line path was hardwired to "
            f"solar regardless of this flag — so a non-solar run would silently have "
            f"produced a SOLAR result under another star's name. Use 'solar'.")
    tei.assert_committed_inputs({
        'CNO cross-arm (RYA-491/237)': CNO_PHASE_A,
        'Mn HFS synth (RYA-473)': MN_JSON,
        'Cu/V HFS synth (RYA-466)': CUV_JSON,
        'Sr II synth (RYA-551)': SR2_JSON,
        'Zr II synth (RYA-560)': ZR2_JSON,
        'Zr II deblend (RYA-585)': ZR2_DEBLEND_JSON,
        'Mg I 5528 synth (RYA-592)': MG5528_JSON,
    })
    eb = tei.assert_engine_b_artifact(args.star)
    print(f"[two-engine] preflight OK — Engine-B per-line {eb['path']} "
          f"({eb['usable_rows']} usable lines); {tei.env_summary()}")

    out_dir = OUT_DIR if args.out_dir is None else (ROOT / args.out_dir)
    p = _solar_params()
    print(f"[two-engine] solar params {p}")

    a_pl = _engine_A_perline(p)
    b_pl = _engine_B_perline()
    ded_b = _dedicated_engine_B()
    print(f"[two-engine] Engine-A species {sorted(a_pl)}")
    print(f"[two-engine] Engine-B(synth-v2) species {sorted(b_pl)}; dedicated Engine-B {sorted(ded_b)}")

    # RYA-520 class: synthesis-required elements never offer raw-EW as Engine A.
    synth_req = {'C', 'N', 'O'} | {
        el for el in TARGET_ELEMENTS
        if (dd := pc.disposition_for(el)) and dd.get('required_treatment') in ('synthesis', 'HFS_sum')}

    species = sorted(set(a_pl) | set(b_pl) | set(ded_b),
                     key=lambda k: (-SOLAR_ASPLUND2021.get(k[0], -9), k[0]))
    records, loud = [], []
    for (el, ion) in species:
        disp = pc.disposition_for(el)
        is_synth_req = el in synth_req
        a_lines = {} if is_synth_req else a_pl.get((el, ion), {})   # RYA-520 suppression
        # Engine-B: dedicated synthesis-harness output (CNO/HFS/Sr) takes precedence
        # for the synthesis-required elements; else the fresh synth-v2 per-line.
        use_dedicated = (el, ion) in ded_b and (is_synth_req or (el, ion) not in b_pl)
        lines = []
        if use_dedicated:
            bv, src = ded_b[(el, ion)]
            lines.append(LineEngines(wavelength=0.0, species=f"{el} {ion}",
                                     a_value=None, a_err=None, b_value=bv, b_err=0.05,
                                     b_chi2=None, ew_mA=None, blend_flag=False,
                                     is_problem_child=True))
        else:
            b_lines = b_pl.get((el, ion), {})
            for w in sorted(set(a_lines) | set(b_lines)):
                a = a_lines.get(w); b = b_lines.get(w)
                lines.append(LineEngines(
                    wavelength=(a or b)['wave'], species=f"{el} {ion}",
                    a_value=(a['a'] if a else None), a_err=0.05,
                    a_in_hull=(a['in_hull'] if a else True),
                    b_value=(b['b'] if b else None), b_err=0.05, b_chi2=None,
                    ew_mA=(a['ew'] if a else None), blend_flag=False,
                    is_problem_child=bool(disp)))
        if not lines:
            continue
        # RYA-525 loud-fail: synthesis-required element with NO Engine-B anywhere.
        synth_required = is_synth_req or (disp and disp.get('required_treatment') in ('synthesis', 'HFS_sum'))
        has_B = any(l.b_value is not None for l in lines)
        if synth_required and not has_B:
            loud.append(f"{el} {ion}: synthesis-required (problem_children "
                        f"{disp['required_treatment']}) but NO Engine-B value")
            continue
        try:
            rec = select_element(f"{el} {ion}", lines)
        except TwoEngineError as e:
            loud.append(f"{el} {ion}: select_element raised: {e}")
            continue
        asp = SOLAR_ASPLUND2021.get(el)
        # ── RYA-674 §2C: the ratified DEMOTION, applied where the record is made ──
        # RYA-558 (Cr II) and RYA-563 (Li / upper_limit) both ratified the same
        # treatment: the floor may compute such a species, but it must be carried as a
        # cross-engine DIAGNOSTIC, never as a reported value. The July artifact recorded
        # Cr II 5.676 and Li I 1.409 in `reported`, and the RYA-527 re-emit adopted both
        # (RYA-669). Demoting here rather than at each consumer is the point of the
        # ticket: the next consumer does not have to remember.
        # Scope is the EXPLICIT, cited exclusion list, not the registry-ion rule that
        # `is_ratified_excluded_species` additionally applies — Ti II / Si II are real
        # exclusions in the SELECTOR but Ryan has ratified no emission-time constraint
        # on them, and RYA-674 adds the three known ones with no interpretation.
        veto = None
        if f"{el} {ion}" in RATIFIED_EXCLUDED_SPECIES:
            veto = f"ratified-excluded species — {exclusion_reason(f'{el} {ion}')} (RYA-240/558)"
        elif is_upper_limit_disposition(el):
            veto = (f"{el} carries the registry upper_limit disposition — the floor may "
                    f"not emit a synthesis point value for it (RYA-563/103/458)")
        records.append(dict(
            element=el, ion=ion, asplund2021=asp,
            reported=(None if veto else round(rec.value, 3)),
            diagnostic_only=bool(veto),
            diagnostic_value=(round(rec.value, 3) if veto else None),
            diagnostic_reason=veto,
            err=round(rec.err, 3), n_lines=rec.n_lines,
            delta_vs_asplund=(round(rec.value - asp, 3) if asp is not None else None),
            engineA=(round(rec.engineA_value, 3) if rec.engineA_value is not None else None),
            engineB=(round(rec.engineB_value, 3) if rec.engineB_value is not None else None),
            selected_engines=list(rec.selected_engines),
            cross_engine_mix=rec.cross_engine_mix, mix_flagged=rec.mix_flagged,
            mean_cross_engine_delta=(round(rec.mean_cross_engine_delta, 3)
                                     if rec.mean_cross_engine_delta is not None else None)))

    if loud:
        raise SystemExit("RYA-525 TWO-ENGINE LOUD-FAIL (synthesis-required missing Engine-B):\n  - "
                         + "\n  - ".join(loud))

    # RYA-674 §2C: gate the per-species records before writing. These are
    # SPECIES_RECORD rows — a diagnostic table may legitimately carry species we would
    # never report — so what is checked is that every vetoed one is marked
    # `diagnostic_only` rather than sitting in `reported` for a consumer to adopt.
    assert_ratified_constraints_satisfied(
        records, 'two-engine floor record emitter (RYA-527/525)', RowKind.SPECIES_RECORD)

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / 'solar_two_engine_records.json').write_text(json.dumps(
        dict(ticket='RYA-527 real two-engine run (RYA-525 floor)',
             gerber_nlte_delta=GERBER_NLTE_DELTA, gerber_xfail=sorted(GERBER_XFAIL),
             records=records), indent=2))
    print(f"\n  el  ion  Asplund  reported  d      engines            engineA engineB  mix")
    for r in records:
        shown = (f"{r['reported']:7.3f}" if r['reported'] is not None
                 else f"[{r['diagnostic_value']:.3f}]")   # RYA-674 demoted to diagnostic
        print(f"  {r['element']:>3s} {r['ion']:<3s} {str(r['asplund2021']):>6s}  "
              f"{shown:>7s}  {str(r['delta_vs_asplund']):>6s}  "
              f"{','.join(e.replace('engine','') for e in r['selected_engines']):<16s} "
              f"{str(r['engineA']):>6s}  {str(r['engineB']):>6s}  "
              f"{'MIX*' if r['mix_flagged'] else ('mix' if r['cross_engine_mix'] else '')}")
    print(f"\n  wrote {out_dir.relative_to(ROOT)}/solar_two_engine_records.json "
          f"({len(records)} species)")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
