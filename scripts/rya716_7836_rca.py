#!/usr/bin/env python3
"""RYA-716 step 3 — why is Al I 7836.134 +0.118 dex high? Is it the RYA-704 collision?

The FINISH-Al brief's hypothesis: the GES list carries TWO Al I entries at exactly
7836.134 (loggf -0.494 and -1.795), and that is a concrete instance of RYA-704 (a
per-line map keyed on a rounded wavelength with no collision check, so two lines of one
species inside 0.1 A silently overwrite).

This tests the hypothesis instead of asserting it, on three axes:

  A. WHAT the two entries actually are -- same transition duplicated, or two distinct
     transitions that happen to coincide? (EP / J / reference decide it.)
  B. Does the DUPLICATE move the abundance, and in which direction? Re-invert the same
     measured EW with the second entry removed from the synthesis line list. If the
     duplicate were inflating A, removing it must drop A by ~the observed +0.118.
  C. Is anything ELSE in the window doing the work? List the window's other absorbers
     and report the blend floor, since 7836's window is the busiest of the four.

Read-only w.r.t. the pinned engine and w.r.t. RYA-704's guard: this DISPOSITIONS Al's
line and does not touch the general collision guard (704 owns that).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
ISPEC_SRC = os.environ.get("ISPEC_DIR", "/srv/codex/engines/ispec_src")
if ISPEC_SRC not in sys.path:
    sys.path.insert(0, ISPEC_SRC)

from config.constants import STAR_PARAMS                       # noqa: E402
from pipeline.abundances_derive import (_load_atmosphere, _load_synth_resources,  # noqa: E402
                                        _ISPEC_SOLAR_ABUND_FILE,
                                        _bisect_synth_abundance)
from pipeline.cno_synthesis import _atom_codes                 # noqa: E402

CENTER, EW_OBS, A_RUN = 7836.134, 71.92586415480511, 6.478515625
CLEAN_REF = 6.35546875           # what the other three lines returned


def main() -> None:
    import ispec
    p = STAR_PARAMS["solar"]
    teff, logg, vturb = float(p["teff"]), float(p["logg"]), float(p["xi"])
    feh = float(p.get("feh", p.get("feh_ref", 0.0)))
    linelist, isotopes, chem = _load_synth_resources()
    solar_abund = ispec.read_solar_abundances(_ISPEC_SOLAR_ABUND_FILE)
    atm = _load_atmosphere(teff, logg, feh, vturb)
    atom_code = int(_atom_codes(["Al"], chem, solar_abund)["Al"])

    w_A = np.asarray(linelist["wave_A"], dtype=float)
    elem = np.asarray([str(x) for x in linelist["element"]])

    # ---- A. what ARE the two entries? -------------------------------------------
    dup = np.where(np.abs(w_A - CENTER) < 0.001)[0]
    print("=" * 78)
    print(f"A. GES entries within 0.001 A of {CENTER}: {len(dup)}")
    cols = ["element", "wave_A", "loggf", "lower_state_eV", "upper_state_eV",
            "lower_j", "upper_j", "reference_code", "theoretical_ew", "molecule"]
    for i in dup:
        vals = []
        for c in cols:
            if c in linelist.dtype.names:
                v = linelist[c][i]
                vals.append(f"{c}={v if isinstance(v, (str, np.str_)) else float(v):.6}"
                            if not isinstance(v, (str, np.str_)) else f"{c}={v}")
        print("   " + "  ".join(vals))

    al_in_win = [i for i in np.where(np.abs(w_A - CENTER) <= 0.6)[0]
                 if elem[i].strip().upper().startswith("AL")]
    print(f"\n   Al I entries in the whole +/-0.6 A window: {len(al_in_win)}")

    # ---- B. does removing the weak duplicate move A? -----------------------------
    kw = dict(atmosphere=atm, teff=teff, logg=logg, feh=feh, vturb=vturb,
              isotopes=isotopes, solar_abund=solar_abund,
              element="Al", atom_code=atom_code)
    w_nm = np.linspace(CENTER - 0.6, CENTER + 0.6, 300) / 10.0

    print("\n" + "=" * 78)
    print("B. re-invert the SAME measured EW with the second Al entry removed")
    gf = np.asarray(linelist["loggf"], dtype=float)
    weak = [i for i in dup if gf[i] < -1.0]
    strong = [i for i in dup if gf[i] >= -1.0]
    print(f"   strong entry loggf={gf[strong[0]]:+.3f}   "
          f"weak entry loggf={gf[weak[0]]:+.3f}" if weak and strong else "   (unexpected)")
    if weak and strong:
        eff = np.log10(10 ** gf[strong[0]] + 10 ** gf[weak[0]])
        print(f"   combined effective loggf = {eff:+.4f}  "
              f"(vs {gf[strong[0]]:+.4f} for the strong entry alone)")
        print(f"   => the duplicate makes the line STRONGER by {eff - gf[strong[0]]:+.4f} dex "
              f"in gf,\n      so removing it should RAISE the fitted A by about "
              f"{eff - gf[strong[0]]:+.4f} dex, not lower it.")

    keep = np.ones(len(linelist), dtype=bool)
    keep[weak] = False
    A_nodup, conv, _ = _bisect_synth_abundance(w_nm, EW_OBS, linelist=linelist[keep], **kw)
    print(f"\n   A(with duplicate, = the run) = {A_RUN:.6f}")
    print(f"   A(duplicate removed)         = {A_nodup:.6f}   converged={conv}")
    print(f"   shift from de-duplicating    = {A_nodup - A_RUN:+.6f} dex")
    print(f"   the outlier to explain       = {A_RUN - CLEAN_REF:+.6f} dex")
    explained = (A_RUN - A_nodup) / (A_RUN - CLEAN_REF) * 100 if A_RUN != CLEAN_REF else 0
    print(f"   => de-duplication explains {explained:+.1f}% of the outlier")

    # ---- C. what else is in the window? ------------------------------------------
    print("\n" + "=" * 78)
    print("C. the rest of the window (top absorbers by loggf, non-Al)")
    win = np.where(np.abs(w_A - CENTER) <= 0.6)[0]
    order = win[np.argsort(-gf[win])]
    for i in order[:10]:
        tag = "  <-- Al" if elem[i].strip().upper().startswith("AL") else ""
        print(f"   {w_A[i]:10.4f}  {elem[i]:<9s} loggf={gf[i]:+7.3f}"
              f"  EP={float(linelist['lower_state_eV'][i]):6.3f}{tag}")

    from pipeline.abundances_derive import _synth_ew_at_abund
    floor = _synth_ew_at_abund(w_nm, atmosphere=atm, teff=teff, logg=logg, feh=feh,
                               vturb=vturb, linelist=linelist, isotopes=isotopes,
                               solar_abund=solar_abund, element="Al",
                               atom_code=atom_code, trial_A=4.0)
    print(f"\n   blend floor at A(Al)=4.0: {floor:.3f} mA "
          f"({100 * floor / (floor + EW_OBS):.1f}% of the inversion target)")


if __name__ == "__main__":
    main()
