#!/usr/bin/env python3
"""
[O I] 6300 differential consistency check (eval gate for the measured-O refactor).

Can we report a MEASURED solar O = (our 1D-LTE A(O)) + (Amarsi-2019 3D-NLTE
differential for the 630 nm forbidden line) that is consistent with BOTH
Caffau 2015's independent 3D absolute (8.73) and the Asplund 2021 gate
(8.69 +/- 0.05)?

PASS -> the Ni-blend / atomic-data seam is empirically small; adopt the
        differential (Part B is safe to commit).
STOP -> the differential lands out of gate or disagrees with Caffau beyond the
        cited uncertainty. That divergence is the seam itself; do NOT commit
        Part B, post the numbers to Linear for Ryan.

Validate-don't-tune: nothing here is fitted. The differential is a pure grid
interpolation at solar params; the anchors are comparison targets only.

RYA-447 (refactor [O I] 6300 -> measured O via Amarsi-2019 3D differential).
"""
from __future__ import annotations
import argparse
import numpy as np
from pipeline import nlte_cno
from config.constants import get_star_params

# Solar reference params come from the single canonical source (stars.yaml via
# get_star_params), NOT a re-declared literal, so the diagnostic and the pipeline
# cannot drift. xi (microturbulence) is the pipeline's vmic for the Sun.
_SUN = get_star_params('solar')
SUN = dict(teff=float(_SUN['teff']), logg=float(_SUN['logg']),
           feh=float(_SUN['feh_ref']), vmic=float(_SUN.get('xi', 1.0)))

# Comparison targets (NOT fit targets):
CAFFAU_OI_ABSOLUTE = 8.73       # Caffau et al. 2015 A&A 579 A88 (POSP III), full-3D [O I]+777
CAFFAU_UNC         = 0.05
ASPLUND_OI_GATE    = (8.69, 0.05)  # Asplund 2021 A&A 653 A141 solar O

OI_6300_AIR_A = 6300.30         # [O I] 6300.30 forbidden line (Ni I 6300.34 pinned in synth)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--a-lte', type=float, default=8.80,
                    help='Our measured 1D-LTE A(O) from the [O I] 6300 joint synthesis '
                         '(HARPS run @ 5d3d60c = 8.80).')
    ap.add_argument('--teff', type=float, default=SUN['teff'])
    ap.add_argument('--logg', type=float, default=SUN['logg'])
    ap.add_argument('--feh',  type=float, default=SUN['feh'])
    ap.add_argument('--vmic', type=float, default=SUN['vmic'])
    a = ap.parse_args()

    print(f"[params]  Teff={a.teff:.0f} logg={a.logg:.3f} [Fe/H]={a.feh:+.2f} "
          f"vmic={a.vmic:.2f}  (canonical solar via get_star_params)")

    # 1) resolve the forbidden line to its grid label
    label = nlte_cno.resolve_line('OI', OI_6300_AIR_A)
    leg = nlte_cno.select_leg(a.teff)
    print(f"[resolve] [O I] {OI_6300_AIR_A} A -> grid label = {label!r}; leg = {leg}")
    if label is None:
        print("CRITICAL: [O I] 6300 did not resolve to an O I grid node "
              "(expected the 630 nm node). STOP -- differential unavailable.")
        return

    # 2) interpolate the 3D-NLTE - 1D-LTE differential at solar params
    delta = nlte_cno.cno_nlte_delta('OI', label, a.teff, a.logg, a.feh, a.vmic, a.a_lte)
    print(f"[grid]    Delta(3D-NLTE - 1D-LTE) = {delta:+.4f} dex  (logeps={a.a_lte:.3f})")
    if not np.isfinite(delta):
        print("CRITICAL: solar query is OUTSIDE the 4D grid hull -> NaN. No silent "
              "LTE. STOP -- differential unavailable for the Sun.")
        return

    # 3) our measured, differentially-corrected O
    a_measured = a.a_lte + delta
    print(f"[result]  measured O = A_lte + Delta = {a.a_lte:.3f} {delta:+.4f} = {a_measured:.3f}")

    # 4) consistency vs the two INDEPENDENT references (comparison only, never fit)
    d_caffau = a_measured - CAFFAU_OI_ABSOLUTE
    gate_c, gate_u = ASPLUND_OI_GATE
    d_gate = a_measured - gate_c
    in_gate = abs(d_gate) <= gate_u
    agree_caffau = abs(d_caffau) <= (CAFFAU_UNC + 0.005)  # cited unc + rounding slack

    print("\n-- consistency --------------------------------------------------")
    print(f"  vs Caffau 2015 absolute {CAFFAU_OI_ABSOLUTE:.2f}: d = {d_caffau:+.3f}  "
          f"-> {'AGREE' if agree_caffau else 'DISAGREE'}")
    print(f"  vs Asplund gate {gate_c:.2f}+/-{gate_u:.2f}:  d = {d_gate:+.3f}  "
          f"-> {'IN GATE' if in_gate else 'OUT OF GATE'}")

    print("\n-- verdict ------------------------------------------------------")
    if in_gate and agree_caffau:
        print("  PASS: measured differential O is consistent with BOTH the Asplund "
              "gate and Caffau's independent absolute. Ni-blend/atomic-data seam is "
              "empirically small -> adopt the differential (commit Part B).")
    else:
        print("  STOP: differential diverges (gate and/or Caffau). Real atomic-data / "
              "Ni-blend finding, NOT a number to tune. Do NOT commit Part B. Post the "
              "numbers above to Linear for Ryan.")


if __name__ == '__main__':
    main()
