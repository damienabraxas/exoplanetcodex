"""RYA-695 — the Kitt Peak synthesis channel, as an Engine-B route.

WHY THIS EXISTS
===============
RYA-673 reported N I, K I, P I and Sc II as `neither`-wired: no Engine A, no
Engine B, and — worse — MEASURED OFF-ORCHESTRATOR. Each carries a real number in
`solar_phase_c_verdict.json` (N 8.188, K 5.099, P 6.61, Sc 3.203) produced by a
channel the two-engine floor cannot see, so the floor could never cross-check
them and the RYA-525 loud-fail never even reached them (a species absent from all
three coverage sources never enters the guard's loop — RYA-673 §2).

That framing contained a factual error worth stating plainly, because it sent
four elements to the wrong follow-on ticket:

    **The Kitt Peak channel is already an Engine-B synthesis measurement.**

`scripts/wire_reference_atlases_rya460.py` measures every one of these windows
with ``pipeline.cno_synthesis._fit_element`` — the *same* Turbospectrum flux-fit
engine the orchestrator already accepts as Engine B for C and O via the RYA-491/237
cross-arm artifact. It is a profile fit against a normalized observed spectrum, not
an EW inversion: the artifact carries ``red_chi2`` and ``fit_s`` per window. So
these four were never missing an engine. They were missing an INVOCATION — exactly
the `NO_HARNESS_INVOCATION` class RYA-680 closed for Co I and Ba II, and the same
shape of defect: a ratified measurement sitting committed in the tree that the
orchestrator never reads.

WHY THE OTHER ENGINE STILL DOES NOT EXIST FOR THESE (evidence, not assumption)
==============================================================================
Wiring this route makes these species `B_only`, not `both`, and that is the honest
answer rather than a shortfall:

* **N I, K I, P I** — their diagnostic lines (N I 7468/8216/8683, K I 7699,
  P I 10581/10596) lie REDWARD of HARPS-VIS (380–690 nm). That is the entire
  reason RYA-459/460 acquired the Kitt Peak atlas. No HARPS line ⇒ no curated EW
  pool entry ⇒ Engine A has nothing to invert, and the synth-v2 leg (which fits
  the same HARPS pool) has nothing to fit either. Engine A is not unwired for
  these; it is out of spectral range.
* **Sc II 4246** — in range, but Sc is `HFS_sum` treatment, so RYA-520 suppresses
  its raw-EW leg by design.
* No second NLTE engine is available for any of them: the Gerber-2023 TS-native
  deck staged on Sirius (`/mnt/codex-data/grids/nlte/gerber_ts`) ships exactly
  ELEVEN model atoms — Ba, Ca, Co, Mg, Mn, Na, Ni, O, Si, Sr, Ti. There is no
  N, K, P or Sc atom to build an Engine-B NLTE leg from.

THE VALUE MUST BE THE VERDICT'S VALUE
=====================================
The KP artifact stores 1D-LTE (`a_1dlte`); the NLTE delta is applied downstream by
`phase_c_verdict_rya371._kittpeak_reclassify` through `pipeline.nlte_corrections`.
This module applies the SAME delta through the SAME subsystem, then
`assert_matches_phase_c()` checks the result against the committed verdict. Two
readers of one measurement that could silently disagree is the RYA-669 defect
shape; RYA-680's Ba comment states the rule ("the floor and the verdict read the
SAME barium") and this enforces it for the KP four.

Everything here is LOUD (RYA-518/680): a present-but-unusable artifact raises with
the reason. It never returns empty and lets an element vanish from the species set.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from config.constants import NLTE_CORRECTION_ELEMENTS

ROOT = Path(__file__).resolve().parent.parent
KITTPEAK_JSON = ROOT / 'data' / 'audit' / 'cno_synthesis' / 'solar_kittpeak_rya460.json'
PHASE_C_JSON = ROOT / 'data' / 'audit' / 'cno_synthesis' / 'solar_phase_c_verdict.json'

#: Agreement required between this route and the phase_c verdict, in dex. The two
#: compute the same quantity from the same artifact through the same grid, so the
#: only expected difference is phase_c's 3-decimal rounding of the N mean.
PHASE_C_TOL = 0.0011


class KittPeakEngineBError(RuntimeError):
    """A Kitt Peak Engine-B route could not produce a usable record (RYA-518)."""


#: element -> (ion, measurement keys, NLTE interpolation nodes in Å, ticket)
#:
#: The NLTE nodes mirror `_kittpeak_reclassify` EXACTLY — they are the wavelengths at
#: which that function interpolates the registered grid, not the window centres. An
#: element with `None` has no registered grid and is reported 1D-LTE, flagged.
#:
#: Co I 3845 is deliberately ABSENT. RYA-564 demoted that blue-edge extraction to
#: diagnostic-only (+1.188 dex on KP SNR~24 / chi2r~3100) and Co's Engine B is the
#: RYA-564 red-line HFS harness, already wired. Reading 3845 here would resurrect a
#: value the Codex has ratified as an artifact.
KITTPEAK_ROUTES = {
    'N':  ('I',  ('NI_7442_7468', 'NI_8216_8223', 'NI_8680_8718'),
           (7468.31, 8216.34, 8683.4), 'RYA-460/556'),
    'K':  ('I',  ('KI_7665_7699',), (7698.964,), 'RYA-460/462'),
    'P':  ('I',  ('PI_10581_10596',), None, 'RYA-460'),
    'Sc': ('II', ('ScII_4246',), None, 'RYA-460'),
}

#: Why no second engine exists, per species. Cited so the audit and the tracker can
#: report a RATIFIED single-engine disposition instead of an open acquisition task.
SINGLE_ENGINE_REASON = {
    'N':  ('N I red 7468/8216/8683 lies redward of HARPS-VIS (380-690 nm) — the reason '
           'the RYA-459 Kitt Peak atlas was acquired. No HARPS line ⇒ no curated EW pool '
           'entry ⇒ no Engine A and no synth-v2 line. The Gerber-2023 TS deck ships no '
           'N atom, so no Engine-B NLTE alternative exists either.'),
    'K':  ('K I 7699 lies redward of HARPS-VIS (7665 additionally sits in the telluric '
           'O2 A-band). No HARPS line ⇒ no Engine A and no synth-v2 line. The '
           'Gerber-2023 TS deck ships no K atom.'),
    'P':  ('P I 10581/10596 is a near-IR multiplet, far redward of HARPS-VIS. Ratified '
           'LTE-only-by-design (RYA-460): no P I NLTE grid is published for either '
           'engine and P I is near-LTE at the Sun.'),
    'Sc': ('Sc II 4246 is in range but Sc is HFS_sum treatment, so RYA-520 suppresses '
           'the raw-EW leg by design. Ratified LTE-only-by-design (RYA-460): no Sc NLTE '
           'grid exists for either engine.'),
}


def _solar_star():
    """The solar node the KP channel interpolates its grids at (phase_c's own default)."""
    return {'teff': 5772.0, 'logg': 4.44, 'feh': 0.0}


def _grid_delta(element: str, wave_A: float, star=None):
    """NLTE delta for one line through `pipeline.nlte_corrections` — the subsystem
    Ca/Ti/Cr/Na/K/N already travel. Returns (delta, flag). Never silently corrects:
    an unregistered element or an out-of-hull star returns (nan, reason)."""
    star = star or _solar_star()
    from pipeline import nlte_corrections as N
    if element not in NLTE_CORRECTION_ELEMENTS:
        return float('nan'), 'NLTE_unavailable (element not registered)'
    if not N.element_grid_in_bounds(element, star['teff'], star['logg'], star['feh']):
        return float('nan'), 'NLTE_unavailable (star out of grid hull)'
    d = N._mpia_element_delta(element, wave_A, star['teff'], star['logg'], star['feh'])
    if d is None or not np.isfinite(d):
        return float('nan'), 'NLTE_unavailable (no grid node within tol)'
    return float(d), NLTE_CORRECTION_ELEMENTS[element].get('flag', 'NLTE_1D')


def load_kittpeak():
    """The RYA-460 artifact, or None if the campaign has not run."""
    if not KITTPEAK_JSON.exists():
        return None
    return json.loads(KITTPEAK_JSON.read_text(encoding='utf-8'))


def kittpeak_engine_B(kp=None):
    """Kitt Peak Engine-B measurements -> {(element, ion): (value, source_tag, basis)}.

    Shape-compatible with `_dedicated_engine_B()` in the two-engine orchestrator, so
    it merges without special-casing at the call site.

    The leg-validation gate is honoured: RYA-460 promotes KP-only elements ONLY if
    the [O I] 6300 / O I 777 overlap cross-check against the Phase-A HARPS/ESPRESSO
    legs agrees. A KP leg that failed its own overlap check must not feed the floor.
    """
    kp = kp if kp is not None else load_kittpeak()
    if not kp:
        return {}
    if not kp.get('leg_validation', {}).get('leg_validated'):
        # Not a silent skip: the leg exists and was REFUSED, and the floor must be able
        # to say which of the two it was (RYA-518).
        print('[two-engine] RYA-460 Kitt Peak leg NOT validated by the [O I]6300 / '
              'O I 777 overlap cross-check -> N/K/P/Sc emit NO Engine-B value')
        return {}

    meas = {m['key']: m for m in kp.get('measurements', [])}
    out = {}
    for el, (ion, keys, nodes, ticket) in KITTPEAK_ROUTES.items():
        rows = [meas[k] for k in keys if k in meas]
        if not rows:
            raise KittPeakEngineBError(
                f"RYA-695: {KITTPEAK_JSON.name} carries none of the {el} windows "
                f"{list(keys)}. The KP campaign declared this element measurable; a "
                f"missing window is a broken artifact, not an element to skip.")
        vals = [r.get('a_1dlte') for r in rows]
        if any(v is None for v in vals):
            raise KittPeakEngineBError(
                f"RYA-695: {el} has a Kitt Peak window with no a_1dlte in "
                f"{KITTPEAK_JSON.name} (windows {list(keys)}, values {vals}). An "
                f"Engine-B route that produced nothing must say so, not emit a partial "
                f"mean over the windows that happened to fit.")

        # N is the one multi-window element, and phase_c averages the three multiplets
        # BEFORE applying NLTE (`n_cross_indicator.atomic_NI_mean`). Reproduced exactly:
        # averaging after correction would differ in the last decimal and trip the
        # phase_c cross-check for no scientific reason.
        if el == 'N':
            nci = kp.get('n_cross_indicator', {})
            a_lte = nci.get('atomic_NI_mean')
            if a_lte is None:
                raise KittPeakEngineBError(
                    f"RYA-695: {KITTPEAK_JSON.name} has no n_cross_indicator."
                    f"atomic_NI_mean — the N I cross-indicator mean phase_c reports "
                    f"from. Engine B cannot be formed for N.")
            a_lte = float(a_lte)
        else:
            a_lte = float(np.mean([float(v) for v in vals]))

        if nodes:
            per = [_grid_delta(el, w) for w in nodes]
            deltas, flag = [p[0] for p in per], per[0][1]
            if all(np.isfinite(d) for d in deltas):
                delta = float(np.mean(deltas))
                value = round(a_lte + delta, 3) if el == 'N' else a_lte + delta
                nlte_txt = (f'{flag} delta {delta:+.4f} applied via pipeline.'
                            f'nlte_corrections (validate-don\'t-tune)')
            else:
                # LOUD, held 1D-LTE, never silently corrected (RYA-409).
                delta, value = float('nan'), a_lte
                nlte_txt = (f'NLTE HELD OUT ({flag}) — grid registered but not '
                            f'interpolable at the solar node; value is 1D-LTE')
                print(f'[two-engine] RYA-695 Kitt Peak {el}: {nlte_txt}')
        else:
            delta, value = float('nan'), a_lte
            nlte_txt = ('1D-LTE — ratified LTE-only-by-design (RYA-458/460); no NLTE '
                        'grid is published for either engine')

        chi2 = [r.get('red_chi2') for r in rows]
        snr = [r.get('cont_snr') for r in rows]
        statuses = {str(r.get('status')) for r in rows}
        if statuses != {'ok'}:
            raise KittPeakEngineBError(
                f"RYA-695: {el} Kitt Peak window(s) {list(keys)} report status "
                f"{sorted(statuses)}, not all 'ok'. A failed flux fit must not be "
                f"averaged into an Engine-B value.")

        out[(el, ion)] = (
            float(value),
            f"{ticket} Kitt Peak {'/'.join(keys)} Turbospectrum flux fit ({nlte_txt})",
            # The KP artifact predates RYA-679 and carries no `reliable` flag: it is a
            # flux fit over an atlas window, so it has neither dEW_dA nor `railed` and
            # the RYA-679 rule is not computable over it. Its own quality statement is
            # recorded instead of asserting a gate that does not exist (the RYA-691 §3A
            # precedent set for the C/O cross-arm artifact).
            f"UNGATED — the RYA-460 Kitt Peak artifact carries no RYA-679 reliability "
            f"flag (an atlas-window flux fit: no dEW_dA / railed). Its own quality "
            f"statement: {len(rows)} window(s) status=ok, red_chi2 {chi2}, "
            f"continuum SNR {snr}; RYA-460 leg VALIDATED by the [O I]6300 / O I 777 "
            f"overlap cross-check vs HARPS/ESPRESSO",
        )
    return out


def assert_matches_phase_c(records, phase_c_path=None):
    """The floor and the verdict must read the SAME Kitt Peak measurement.

    `records` is the mapping returned by `kittpeak_engine_B()`. Raises when a value
    disagrees with the committed phase_c verdict beyond `PHASE_C_TOL`, which would
    mean the two consumers had drifted onto different numbers for one measurement —
    the RYA-669 defect shape, and precisely what RYA-680 forbade for Ba.

    Silent when the verdict is absent (a fresh checkout may legitimately not have
    regenerated it yet); a MISMATCH is never silent.
    """
    path = Path(phase_c_path) if phase_c_path else PHASE_C_JSON
    if not path.exists():
        return {}
    doc = json.loads(path.read_text(encoding='utf-8'))
    pc = {v['element']: v.get('A_measured') for v in doc.get('verdicts', [])}
    checked = {}
    for (el, _ion), (value, _tag, _basis) in records.items():
        ref = pc.get(el)
        if ref is None:
            continue
        if abs(float(value) - float(ref)) > PHASE_C_TOL:
            raise KittPeakEngineBError(
                f"RYA-695: the Kitt Peak Engine-B route computes A({el}) = "
                f"{float(value):.4f} but the phase_c verdict at {path.name} carries "
                f"{float(ref):.4f} (|Δ| = {abs(float(value) - float(ref)):.4f} > "
                f"{PHASE_C_TOL}). Both read {KITTPEAK_JSON.name} through "
                f"pipeline.nlte_corrections, so they cannot legitimately differ: one "
                f"of the two has drifted. The two-engine floor and the verdict must "
                f"report the SAME measurement (RYA-680).")
        checked[el] = (float(value), float(ref))
    return checked
