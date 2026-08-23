#!/usr/bin/env python3
"""
scripts/rya1001_crires_coverage.py — RYA-1001
=============================================
REAL-PIXEL coverage of the Al I NIR target lines in the solar CRIRES+ Vesta IDPs.

Why this exists separately from `pipeline.coverage`: that module cannot see this
holding at all. `coverage.load_registry()` skips any holding whose `manifest_path` is
not a spectrum-location CSV, and the skip is a bare `continue`. Six of the ten solar
holdings are invisible to it that way, including all three CRIRES+ ones — so it reports
ZERO instruments for Al I 13123.42 and 13150.75, the two best-graded Al lines we have.
That is precisely the failure RYA-708 was built to end (RYA-707 published "Al 7835/8772:
NO DATA" off a single hardcoded file), reappearing one layer down.

And header endpoints do not prove coverage (codex-data-audit, Science coverage): a
setting's declared span says nothing about inter-detector gaps or flagged pixels. This
tests actual QUAL==0, finite, non-zero pixels.

Linear issue: RYA-1001
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from config.constants import require_codex_path            # noqa: E402
from pipeline.crires_telluric import load_crires_idp        # noqa: E402
from pipeline.wavelength_util import air_to_vac             # noqa: E402

MANIFEST = ROOT / "data/audit/vesta_crires_plus/vesta_crires_plus_idp_manifest.csv"
SPECTRA = require_codex_path("data.spectra_vesta_crires")   # RYA-810: never a literal
OUT = ROOT / "data" / "results" / "rya1001" / "rya1001_crires_coverage.csv"

#: Al I NIR lines to test, AIR wavelengths as they appear in linelist_solar.csv.
#: The Burheim-graded set plus the two 1.09 um lines the litscan flagged.
TARGETS = {
    10872.973: "Al I 10872.97 — 1.089 um, the line NL17 dropped for tellurics",
    10891.736: "Al I 10891.74",
    11253.189: "Al I 11253.19 — Burheim 5%",
    11254.924: "Al I 11254.92 — Burheim 2%",
    12749.909: "Al I 12749.91 — Burheim 9%",
    12757.275: "Al I 12757.27 — Burheim 11%",
    13123.416: "Al I 13123.42 — Burheim 1.5%",
    13150.753: "Al I 13150.75 — Burheim 3.1%",
}
#: Half-window. A line needs real pixels across its core, not one pixel at the edge.
HALF_A = 1.0


def main() -> int:
    man = pd.read_csv(MANIFEST)
    frames = []
    for f in sorted(set(man.file)):
        p = SPECTRA / f
        if not p.exists():
            print(f"  MISSING {p}")
            continue
        frames.append((f, load_crires_idp(p)))
    if not frames:
        raise SystemExit(f"no CRIRES+ IDP readable under {SPECTRA} — refusing to report "
                         f"'not covered', which is what an unreadable holding would look "
                         f"like (RYA-833: an absence needs a positive control).")
    print(f"loaded {len(frames)} CRIRES+ Vesta IDP frames\n")

    rows = []
    for lam_air, label in sorted(TARGETS.items()):
        lam_vac = float(air_to_vac(np.array([lam_air]))[0])
        best = None
        for fname, fr in frames:
            for s in fr.segments:
                w = s.wave_A
                if not len(w):
                    continue
                m = (w >= lam_vac - HALF_A) & (w <= lam_vac + HALF_A)
                if not m.any():
                    continue
                good = m & np.isfinite(s.flux) & (s.flux != 0) & (s.qual == 0)
                rec = dict(file=fname, wlen_id=fr.wlen_id, order=s.order,
                           detector=s.detector, n_in_window=int(m.sum()),
                           n_good=int(good.sum()),
                           px_lo=float(w[m].min()), px_hi=float(w[m].max()))
                if best is None or rec["n_good"] > best["n_good"]:
                    best = rec
        covered = bool(best and best["n_good"] > 0)
        print(f"{label:56s} lam_vac={lam_vac:10.3f}  "
              + (f"{best['wlen_id']:6s} ord{best['order']}/det{best['detector']} "
                 f"{best['n_good']:3d}/{best['n_in_window']:3d} good px in +-{HALF_A} A"
                 if covered else "NO PIXELS IN ANY SEGMENT"))
        rows.append(dict(line=label, lam_air_A=lam_air, lam_vac_A=round(lam_vac, 4),
                         covered=covered, **(best or {})))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT, index=False)
    print(f"\nwrote {OUT}")
    print("NOTE: coverage is NOT permission. `telluric_policy.gate_holding` REFUSES "
          "solar_vesta_crires_plus_idp (telluric_applied=not-applied, crires_plus "
          "telluric_required=yes). Every line above needs the RYA-424 stage first.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
