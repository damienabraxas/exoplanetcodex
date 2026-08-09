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
    Cu/V RYA-466, Sr II RYA-551, Co I RYA-564, Ba II RYA-581) — an Engine-B output,
    injected as a single synth line, clearly sourced. A fully-fresh re-run of those
    harnesses is the Sirius step; here they carry their committed synthesis value.

Preflight (scripts/rya527_preflight_reconciliation.py) MUST be green first.

RYA-680 — Co I and Ba II are wired (they were not)
--------------------------------------------------
RYA-673's Engine-B wiring audit classified both `NO_HARNESS_INVOCATION`: the
ratified synthesis result was committed in the repo and `_dedicated_engine_B()`
never read it. The consequence was not a wrong number, it was NO number — neither
species entered the species loop at all, so RYA-525's own loud-fail could not see
them either (it iterates the union of the three coverage sources; a species absent
from all three never enters the loop). Gate 3 read UNEVALUABLE for both, and no
amount of measurement quality could move that. Both are read here now.

RYA-691 — the `reliable` contract, honoured at every read
---------------------------------------------------------
Each dedicated read now returns a RELIABILITY BASIS alongside its value, and that
basis is written into the record. Three states, and they are not interchangeable:

  * the artifact carries the RYA-679 flag and it is True  -> gated, value admitted;
  * the artifact carries the flag and it is False         -> RAISE (never silent);
  * the artifact carries no flag at all                   -> admitted, but recorded
    as UNGATED with the reason, so "nothing gated this" can never be read as
    "this was gated".

The third state is not a loophole, it is a fact about four of these artifacts: the
CNO cross-arm (RYA-491/237) and the Mn / Cu / V HFS harnesses (RYA-473/466) are not
in-window profile fits and predate the RYA-679 rule, so they have no `dEW_dA` and no
`railed` to compute it from. Fabricating a uniform check over them would assert a
gate that does not exist. What is forbidden is SILENCE about which state applies.
"""
import argparse
import json
import sys
import warnings
from collections import Counter
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
from pipeline.reliability_contract import reliability_basis  # RYA-691, ratified RYA-699
from config.constants import (get_star_params, TARGET_ELEMENTS, NLTE_CORRECTION_ELEMENTS,
                              SOLAR_ASPLUND2021)
import pipeline.problem_children as pc
from pipeline import two_engine_inputs as tei
from pipeline import kittpeak_engine_b as keb   # RYA-695 Engine-B route for N/K/P/Sc

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
CO_JSON = ROOT / 'data' / 'results' / 'co_synthesis_rya564.json'          # RYA-564 Co I red HFS
# RYA-680: Ba reads the RYA-581 IN-WINDOW DEBLEND, never the RYA-559 EW->COG file.
# `data/results/solar_ba_synthesis_rya559.json` is still in the tree (it is the
# measurement RYA-581 supersedes and its cross-checks are the evidence for the
# supersession) and it holds A(Ba) 2.410 — a value RYA-581 demonstrated is inflated
# by a blend_flag=True pool EW that an EW inversion cannot deblend. Wiring that file
# here would silently undo a merged result, so this driver does not know its path at
# all: there is nothing to fall back to and nothing to get wrong. phase_c makes the
# same choice (`_ba_reclassify` dispatches to `_ba_reclassify_deblend` on ticket
# RYA-581), so the floor and the verdict read the SAME barium.
BA_DEBLEND_JSON = ROOT / 'data' / 'results' / 'solar_ba_deblend_rya581.json'


class DedicatedEngineBError(RuntimeError):
    """A dedicated Engine-B route could not produce a usable record.

    RYA-680/518: raised, never swallowed. The failure mode this replaces is the one
    that made this ticket necessary — a route that returns nothing, an element that
    silently leaves the species set, and a gate that reads UNEVALUABLE with no trace
    of why. Every message here names the artifact, the species and the reason.
    """


def _reliability_basis(obj, what, absent_reason, key='reliable'):
    """RYA-691: THE `reliable` contract for one dedicated Engine-B read.

    Returns the basis string to record with the value. Raises when the flag is
    present and false — an unreliable measurement must never reach the map, and it
    must never reach it quietly either (RYA-518, rule ratified RYA-679:
    ``reliable = (not railed) AND dEW_dA >= 40.0``, produced by
    ``scripts/solar_profile_fit.assess_reliability``; this is the CONSUMER side and
    deliberately re-derives nothing).

    `absent_reason` is required, not optional: an artifact with no flag is a real
    and legitimate state, but only if the record says so in words.

    RYA-699 moved the two basis strings into `pipeline.reliability_contract` so the
    emission-time gate can READ back what this WRITES. The behaviour here is
    unchanged; what changed is that a module which never calls this helper can no
    longer emit an unreadable basis quietly.
    """
    if key not in obj:
        return reliability_basis(None, key=key, absent_reason=absent_reason)
    if bool(obj.get(key)):
        return reliability_basis(True, key=key)
    raise DedicatedEngineBError(
        f"RYA-691: {what} is marked {key}=False and MUST NOT be emitted. The RYA-679 "
        f"rule (not railed AND dEW_dA >= 40.0) demoted this measurement; a consumer "
        f"that used it anyway would carry a demotion the artifact recorded and the "
        f"record hid. Fix the measurement or hold the element owed — do not exempt "
        f"this read.")


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


def _nlte_or_lte(rec, what, tag):
    """RYA-691 §2/§3C — resolve the ENGINE, explicitly, and say which one it was.

    Replaces ``v = rec.get('A_nlte') or rec.get('A_lte_median')``. That construct had
    two defects and one of them is live:

      * it substituted the LTE median for a missing NLTE value under the SAME
        provenance tag, so a consumer could not tell an NLTE result from an LTE
        fall-back. **This fires today for V I** — `data/audit/cu_v_hfs_synthesis/`
        records `nlte_void: true` (no V I grid, RYA-466), so V's emitted value has
        always been A_lte_median wearing 'RYA-466 HFS synth'. The VALUE is right and
        does not change here; what changes is that the record now says LTE.
      * `or` is falsy-triggered, so `A_nlte = 0.0` would also fall through. On the
        A(X) = 12 + log(N_X/N_H) scale a literal 0.0 means N_X/N_H = 1e-12 — nine
        orders below the rarest element measured here, and no artifact in the repo
        carries it — so the case is NOT reachable with today's inputs. It is fixed
        anyway: a numeric test is `is None`, never truthiness, and "unreachable
        today" is a property of the data, not of the code.

    Returns (value, tag-with-engine).
    """
    v = rec.get('A_nlte')
    if v is not None:
        return float(v), f"{tag} (NLTE)"
    v = rec.get('A_lte_median')
    if v is not None:
        why = 'nlte_void — no grid (RYA-466)' if rec.get('nlte_void') else 'no A_nlte in artifact'
        return float(v), f"{tag} (LTE FALL-BACK: {why})"
    raise DedicatedEngineBError(
        f"RYA-691: {what} carries neither 'A_nlte' nor 'A_lte_median' — there is no "
        f"value to emit. An Engine-B route that produces nothing must say so, not "
        f"drop the element (RYA-680).")


def _dedicated_engine_B():
    """Committed Engine-B SYNTHESIS-HARNESS measurements for the synthesis-required
    elements (CNO nlte_cno primary indicator; Mn/Cu/V HFS synth; Sr II synth; Zr II;
    Mg I 5528; RYA-680: Co I red HFS and Ba II 5853 in-window deblend).

    Returns {(element, ion): (value, source_tag, reliability_basis)}. The basis is
    carried, not discarded: `source_tag` used to be unpacked into a local and thrown
    away, so the emitted record said nothing at all about where an Engine-B value came
    from or what cleared it.

    Every route here is LOUD (RYA-680/518). A route whose artifact is present but
    unusable RAISES with the reason; it never returns empty and leaves the element to
    vanish from the species set — which is precisely how Co and Ba came to read
    UNEVALUABLE while their measurements sat committed in the tree.
    """
    out = {}
    if CNO_PHASE_A.exists():
        cross_arm = json.loads(CNO_PHASE_A.read_text()).get('cross_arm', {})
        for el in ('C', 'O'):
            ca = cross_arm.get(el, {})
            prims = [i for i in ca.get('indicators', [])
                     if i.get('role') == 'primary' and i.get('A') is not None]
            if not prims:
                raise DedicatedEngineBError(
                    f"RYA-680: the CNO cross-arm artifact {CNO_PHASE_A.name} has no "
                    f"primary indicator carrying a value for {el} — Engine B cannot be "
                    f"formed. {el} is synthesis-required (RYA-520), so this is a broken "
                    f"input, not an element to skip.")
            prim = prims[0]
            # RYA-691 §3A — C/O are NOT brought under the flag, and this is the reason.
            # `role == 'primary'` is a SELECTION rule (which indicator speaks for the
            # element), not a reliability rule, and the RYA-491/237 cross-arm artifact
            # carries no `reliable` key on any indicator: it is a multi-indicator
            # cross-arm reconciliation, not an in-window profile fit, so it has neither
            # `dEW_dA` nor `railed` and the RYA-679 rule is not computable over it. Its
            # own quality statement is the cross-arm verdict + primary spread, which is
            # what gets recorded. Asserting a gate here would be asserting one that does
            # not exist.
            n_prim = len(prims)
            basis = reliability_basis(None, absent_reason=(
                f"selection rule role='primary' ({prim.get('key')}"
                + (f", FIRST of {n_prim} primary indicators in artifact order"
                   if n_prim > 1 else "")
                + f"); the RYA-491/237 cross-arm artifact carries no RYA-679 "
                  f"reliability flag (not a profile fit: no dEW_dA / railed). Its "
                  f"own quality statement: verdict={ca.get('verdict')!r}, "
                  f"primary spread {ca.get('spread')}"))
            out[(el, 'I')] = (float(prim['A']),
                              f"nlte_cno synthesis {prim.get('key')} (RYA-491/237)", basis)
    if MN_JSON.exists():
        m = json.loads(MN_JSON.read_text()).get('Mn', {})
        v, tag = _nlte_or_lte(m, 'Mn I (RYA-473 HFS synth)', 'RYA-473 HFS synth')
        n_ok = sum(1 for l in m.get('per_line', []) if str(l.get('status')) == 'ok')
        out[('Mn', 'I')] = (v, tag, _reliability_basis(
            m, 'Mn I (RYA-473 HFS synth)',
            f"the RYA-473 HFS flux-fit artifact carries no RYA-679 reliability flag "
            f"(not an in-window profile fit); its own quality statement is per-line "
            f"status {n_ok}/{len(m.get('per_line', []))} ok, "
            f"n_lines={m.get('n_lines')}, scatter={m.get('scatter')}"))
    if CUV_JSON.exists():
        d = json.loads(CUV_JSON.read_text())
        for el in ('Cu', 'V'):
            e = d.get(el, {})
            if not e:
                raise DedicatedEngineBError(
                    f"RYA-680: {CUV_JSON.name} has no '{el}' block — Engine B cannot be "
                    f"formed for {el} from a route the preflight declared present.")
            v, tag = _nlte_or_lte(e, f'{el} I (RYA-466 HFS synth)', 'RYA-466 HFS synth')
            n_ok = sum(1 for l in e.get('per_line', []) if str(l.get('status')) == 'ok')
            out[(el, 'I')] = (v, tag, _reliability_basis(
                e, f'{el} I (RYA-466 HFS synth)',
                f"the RYA-466 HFS flux-fit artifact carries no RYA-679 reliability flag "
                f"(not an in-window profile fit); its own quality statement is per-line "
                f"status {n_ok}/{len(e.get('per_line', []))} ok, "
                f"n_lines={e.get('n_lines')}, scatter={e.get('scatter')}"))
    if SR2_JSON.exists():
        # RYA-691: Sr II is the ONE read of the eight that had a `reliable` flag sitting
        # in the artifact and did not consult it — six lines above the Zr block that
        # says "RELIABILITY-GATED throughout". 4077.709 is the RYA-551 primary line;
        # its HARPS fit is reliable=True (dEW_dA 203.5, not railed), so gating it
        # changes nothing today. Under any future demotion it now raises instead of
        # carrying a value the artifact had already marked unusable.
        sr = json.loads(SR2_JSON.read_text()).get('4077.709', {}).get('harps', {})
        v = sr.get('A_NLTE')
        if v is None:
            raise DedicatedEngineBError(
                f"RYA-680: {SR2_JSON.name} has no A_NLTE on the 4077.709 HARPS fit — "
                f"the RYA-551 primary Sr II line. Sr II is synthesis-required; a "
                f"missing primary is a broken artifact, not a silent skip.")
        out[('Sr', 'II')] = (float(v), 'RYA-551 Sr II synth (4077.709 HARPS)',
                             _reliability_basis(sr, 'Sr II 4077.709 HARPS (RYA-551)',
                                                'no reliable key on the 4077.709 HARPS fit'))
    if CO_JSON.exists():
        # ── RYA-680: Co I — wired. Previously NO_HARNESS_INVOCATION (RYA-673) ──────
        # RYA-564 measured A(Co) on clean RED Co I lines by HFS-resolved Turbospectrum
        # flux fit with a PER-LINE Gerber TS-native 1D-NLTE delta, after demoting the
        # untrusted blue-edge 3845 artifact (+1.188 dex, KP SNR~24) to diagnostic-only.
        # The value emitted here is `_summary.A_Co` — the median over the lines that
        # cleared the RYA-679 floor — which is the SAME field phase_c's `_co_reclassify`
        # reads, so the floor and the verdict cannot drift onto different Co.
        #
        # NOTE the NLTE delta is already inside these numbers (per line, from the
        # RYA-534-validated grid). GERBER_NLTE_DELTA['Co'] = +0.099 is the synth-v2
        # leg's element-level delta and is deliberately NOT applied on top — that would
        # double-count the same physics.
        co = json.loads(CO_JSON.read_text())
        s = co.get('_summary', {})
        a, lines_ok = s.get('A_Co'), (s.get('reliable_lines') or {})
        if a is None or not lines_ok:
            raise DedicatedEngineBError(
                f"RYA-680: Co I (RYA-564) produced no reportable value — "
                f"A_Co={a!r}, n_reliable={s.get('n_reliable')!r}, "
                f"reason={s.get('reason')!r}. RYA-564's ratified rule is that if no red "
                f"line clears the reliability floor the element reports NO VALUE and the "
                f"blue-edge 3845 artifact is NEVER a fall-back — so this raises rather "
                f"than reaching for it. Co stays measurable-owed.")
        # The summary's reliable set must agree with the per-line flags it was built
        # from; two views of one fact that could disagree is the RYA-669 defect shape.
        for w in lines_ok:
            if not bool((co.get(w, {}).get('harps') or {}).get('reliable')):
                raise DedicatedEngineBError(
                    f"RYA-691: Co I {w} is listed in _summary.reliable_lines but its "
                    f"HARPS fit is not marked reliable — {CO_JSON.name} contradicts "
                    f"itself and the value it summarises cannot be trusted.")
        out[('Co', 'I')] = (
            float(a),
            f"RYA-564 Co I red HFS synth (NLTE; median of {len(lines_ok)} reliable "
            f"lines {', '.join(sorted(lines_ok))})",
            reliability_basis(True, detail=(
                f"{len(lines_ok)} of "
                f"{len([k for k in co if not k.startswith('_')])} fitted lines cleared "
                f"(not railed AND dEW_dA >= 40.0)")))
    if BA_DEBLEND_JSON.exists():
        # ── RYA-680: Ba II — wired to the RYA-581 DEBLEND. Read this before editing ──
        # A(Ba) here is 2.237 (RYA-581 in-window blend-fit), NOT 2.410 (RYA-559
        # EW->COG). RYA-673's audit map points `HARNESS_RESULTS['Ba']` at the RYA-559
        # file; following that map would wire the superseded, blend-inflated value back
        # in and silently undo a merged result. The ticket asserts the ticket: this
        # route refuses any artifact that is not RYA-581.
        ba = json.loads(BA_DEBLEND_JSON.read_text())
        if str(ba.get('ticket')) != 'RYA-581':
            raise DedicatedEngineBError(
                f"RYA-680: {BA_DEBLEND_JSON.name} declares ticket "
                f"{ba.get('ticket')!r}, not 'RYA-581'. The two-engine floor reads the "
                f"in-window DEBLEND only — the RYA-559 EW->COG value 2.410 is "
                f"superseded (an EW inversion cannot deblend: the blend_flag=True pool "
                f"EW 74.62 mA vs the clean line ~64 mA was charged entirely to Ba) and "
                f"must never be substituted for it.")
        v = ba.get('A_nlte')
        if v is None:
            raise DedicatedEngineBError(
                f"RYA-680: {BA_DEBLEND_JSON.name} carries no A_nlte — the RYA-581 "
                f"profile fit produced no NLTE value. Ba II is synthesis-required; "
                f"there is no second route and no fall-back.")
        out[('Ba', 'II')] = (
            float(v),
            f"RYA-581 Ba II 5853.668 in-window deblend synth + Korotin2015 NLTE "
            f"(delta {ba.get('engineA_korotin_delta')}; supersedes RYA-559 2.410)",
            _reliability_basis(ba, 'Ba II 5853.668 (RYA-581 deblend)',
                               'no reliable key on the RYA-581 deblend artifact'))
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
            out[('Zr', 'II')] = (float(np.mean(rel)), f"{_tag} (n={len(rel)} reliable)",
                                 reliability_basis(True, detail=(
                                     f"{len(rel)} line(s) cleared "
                                     f"(not railed AND dEW_dA >= 40.0)")))
            break
    if ('Zr', 'II') not in out:
        # RYA-691: gated-shut is a legitimate outcome; gated-shut IN SILENCE is not.
        # Mg prints its hold-out below and Zr did not, so the only trace of Zr's was a
        # comment in this file. Zr then vanishes from the species set entirely — it is
        # absent from the EW pool and from synth-v2 too — and RYA-525's loud-fail never
        # sees it, because that guard iterates the union of the three coverage sources
        # and a species in none of them never enters the loop (RYA-673 §2). That hole is
        # NOT closed here (closing it aborts the run on six more species and is a
        # science decision, not a refactor); it is made audible.
        print(f"[two-engine] RYA-585/560 Zr II HELD OUT: no line cleared the RYA-679 "
              f"floor (not railed AND dEW_dA >= 40.0) in either artifact -> Zr II "
              f"emits NO Engine-B value and does not enter the species set")
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
                                'RYA-592 Mg I 5528 in-window blend-fit synth (concordant)',
                                reliability_basis(True, detail=(
                                    'target_reliable=True, AND RYA-592 '
                                    'concordance-gated: lines_concordant=True')))
        else:
            print(f"[two-engine] RYA-592 Mg I 5528 HELD OUT: reliable="
                  f"{v.get('target_reliable')}, concordant={v.get('lines_concordant')} "
                  f"(|d| {v.get('concordance_worst_abs_dex')} vs band "
                  f"{v.get('concordance_band')}) -> Mg stays single-line")
    # ── RYA-695: N I / K I / P I / Sc II — the Kitt Peak synthesis channel ────────
    # These four read `neither`-wired in RYA-673 and MEASURED OFF-ORCHESTRATOR: each
    # carries a real number in the verdict produced by a channel this floor could not
    # see. The audit's framing was wrong in a way that mattered — the Kitt Peak channel
    # is ALREADY an Engine-B synthesis measurement. `wire_reference_atlases_rya460.py`
    # fits every window with `cno_synthesis._fit_element`, the same Turbospectrum
    # flux-fit engine this driver already accepts as Engine B for C and O. So these
    # were never missing an engine, only an INVOCATION — the NO_HARNESS_INVOCATION
    # class RYA-680 closed for Co I and Ba II.
    #
    # They become B_only, not both, and that is the honest outcome: their lines lie
    # redward of HARPS-VIS (which is why RYA-459 acquired the atlas), so Engine A is
    # out of spectral range rather than unwired, and the Gerber-2023 deck ships no N,
    # K, P or Sc atom. See pipeline/kittpeak_engine_b.SINGLE_ENGINE_REASON.
    kp_records = keb.kittpeak_engine_B()
    if kp_records:
        # One measurement, two consumers: the floor's value must equal the verdict's
        # (RYA-680's "the floor and the verdict read the SAME barium", generalised).
        keb.assert_matches_phase_c(kp_records)
        for _k, _v in kp_records.items():
            if _k in out:
                raise DedicatedEngineBError(
                    f"RYA-695: the Kitt Peak route and an existing dedicated route both "
                    f"claim {_k[0]} {_k[1]}. Two Engine-B values for one species must be "
                    f"adjudicated, never silently overwritten by dict order.")
            out[_k] = _v
        print(f"[two-engine] RYA-695 Kitt Peak Engine-B wired: "
              f"{', '.join(f'{e} {i}' for e, i in sorted(kp_records))}")
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
        'Co I red HFS synth (RYA-564)': CO_JSON,          # RYA-680
        'Ba II 5853 deblend (RYA-581)': BA_DEBLEND_JSON,  # RYA-680 — NOT the RYA-559 file
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
        b_source = b_reliability = None
        if use_dedicated:
            # RYA-691: `src` used to be unpacked and dropped on the floor, so the record
            # said nothing about which harness produced the Engine-B value or what (if
            # anything) cleared it. Both now travel with the number.
            bv, b_source, b_reliability = ded_b[(el, ion)]
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
            # RYA-691: provenance of a DEDICATED Engine-B value — which harness, and
            # what cleared it. None for the synth-v2 leg, whose provenance is the
            # per-line artifact named in the preflight line.
            engineB_source=b_source, engineB_reliability=b_reliability,
            cross_engine_mix=rec.cross_engine_mix, mix_flagged=rec.mix_flagged,
            mean_cross_engine_delta=(round(rec.mean_cross_engine_delta, 3)
                                     if rec.mean_cross_engine_delta is not None else None),
            # RYA-695: WHY each engine won, not just which one did.
            #
            # `engine_selection.select_line` has always attached a `reason` and a
            # `regime` to every LineWinner — "clean-weak line -> 1D-NLTE", "hard line
            # (blend/saturation/HFS) -> synthesis", "only Engine-B eligible",
            # "indeterminate regime -> lower line-scatter sigma" — and this emitter
            # dropped all of it, keeping only the engine names. So the artifact could
            # say Fe I mixed both engines but not that 62 lines went to A as clean-weak
            # while 19 went to B as blends, which is the actual selection story and the
            # thing a reader needs to judge whether the mix is sound.
            #
            # Emitted as a COUNT per distinct reason (and per regime) rather than a
            # per-line list: the per-line detail is already reconstructible from the
            # inputs, and an 81-entry array per element would bury the summary that
            # makes the table readable. Sorted by descending count, then text, so the
            # artifact is diffable across runs.
            selection_reasons=[
                {'reason': r, 'engine': e, 'n_lines': n}
                for (r, e), n in sorted(
                    Counter((w.reason, w.engine) for w in rec.per_line).items(),
                    key=lambda kv: (-kv[1], kv[0]))],
            line_regimes=[
                {'regime': g, 'n_lines': n}
                for g, n in sorted(Counter(w.regime for w in rec.per_line).items(),
                                   key=lambda kv: (-kv[1], kv[0]))]))

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
