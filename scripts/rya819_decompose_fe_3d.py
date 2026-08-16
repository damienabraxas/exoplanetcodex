#!/usr/bin/env python3
"""
RYA-819 — decompose A(Fe) into 1D-LTE / 1D-NLTE / 3D-atmosphere / 3D-NLTE, PER LINE.

THE QUESTION. Gold derives 7.466 as "1D-NLTE 7.516 MINUS the Magic-2013 1D->3D offset
(FE_1D3D_SOLAR_OFFSET = 0.05)". RYA-817 measured the Amarsi-2022 full 3D-NLTE net at
-0.002 dex. With delta_NLTE ~ +0.010 that implies a 3D-atmosphere term near -0.012, not
-0.05. Is the chain right physics, or right-answer-by-compensating-error?

⚠️ EVERY TERM IS PAIRED PER LINE, NEVER AGGREGATE-MINUS-AGGREGATE. The published band
products are computed over DIFFERENT line subsets (1D-LTE n=148, ENGINE-A n=105,
ENGINE-A-3DNLTE n=114), so differencing their medians measures the subsets as much as the
physics — the RYA-785 wrong-referee shape. Here a line contributes to a term only if it
carries BOTH sides of that term, and every n is reported.

⚠️ THE TWO CORRECTIONS ARE DIFFERENT QUANTITIES AND ARE KEPT APART (the ticket's third
"do NOT"). Magic-2013 is a <3D>-minus-1D ATMOSPHERE shift computed in LTE. Amarsi is a full
3D-NLTE-minus-1D-LTE NET, inside which the 3D and NLTE terms partly cancel. The only
like-for-like comparison against Magic is
    3D-atmosphere := A(3D-NLTE) - A(1D-NLTE)
which is what this computes. It is still not exactly Magic's quantity — Magic holds NLTE
out of both sides, this holds it in both — and that residual difference is named in the
verdict rather than glossed.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from config.constants import FE_1D3D_SOLAR_OFFSET, get_star_params      # noqa: E402
from pipeline.nlte_corrections import _mpia_fe_delta                    # noqa: E402

SRC = ROOT / 'data' / 'results' / 'rya817' / 'rya817_3dnlte_per_line.csv'
TRAINING_CSV = (ROOT / 'data' / 'reference' / 'amarsi2022_training'
                / 'amarsi2022_training_lines.csv')
OUT = ROOT / 'data' / 'results' / 'rya819'
#: A line is "the same line" as Amarsi's if the air wavelengths agree this closely. The
#: golden list is a curated set of well-separated lines, so this is a rounding tolerance,
#: not a search radius.
GOLDEN_TOL_A = 0.05


def _stat(v: np.ndarray) -> dict:
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return {'n': 0, 'median': None, 'mad': None, 'sem': None}
    med = float(np.median(v))
    mad = float(1.4826 * np.median(np.abs(v - med)))
    return {'n': int(v.size), 'median': med, 'mad': mad,
            'sem': float(mad / np.sqrt(v.size))}


def main() -> None:
    p = get_star_params('solar')
    teff, logg = float(p['teff']), float(p['logg'])
    feh = float(p.get('feh', p.get('feh_ref', 0.0)))

    d = pd.read_csv(SRC)
    d = d[(d.ion == 'I') & (d.band == 'VIS')].copy()

    # δ_NLTE per line from the SAME MPIA grid Engine-A uses — not a global constant.
    d['delta_nlte'] = [
        _mpia_fe_delta('I', float(w), teff, logg, feh)
        for w in d.wavelength_air_A
    ]
    d['a_1dnlte'] = d.a_1dlte + d.delta_nlte
    # the like-for-like comparator against Magic-2013
    d['term_3d_atmosphere'] = d.a_3dnlte - d.a_1dnlte

    ind = d[d.in_domain == True]                                        # noqa: E712

    # ⚠️ THE SPEC SAYS "THE SAME LINE SET (the Amarsi golden lines)", AND THAT IS NOT THE
    # SAME AS "in-domain". The domain check is a FEATURE-BOX test: a line passes if its
    # (Elo, Eup, log gf) sit inside the trained ranges. Many of our lines pass that box
    # without being lines Amarsi actually trained or published on. Restricting to the
    # recovered golden list is what makes our decomposition comparable to his -0.002,
    # because it removes line selection as an explanation for any difference.
    gold = pd.read_csv(TRAINING_CSV)
    gold_I = gold[gold.ion.astype(str).str.strip() == 'I']
    gw = gold_I.wavelength_air_A.to_numpy(float)
    d['on_golden_list'] = [
        bool(np.min(np.abs(gw - float(w))) <= GOLDEN_TOL_A) for w in d.wavelength_air_A
    ]
    ind = ind.assign(on_golden_list=[
        bool(np.min(np.abs(gw - float(w))) <= GOLDEN_TOL_A)
        for w in ind.wavelength_air_A])

    both = ind[ind.delta_nlte.notna() & ind.aberr.notna()]
    golden = both[both.on_golden_list]

    print(f"Fe I VIS lines in the RYA-817 run : {len(d)}")
    print(f"  in the Amarsi training domain   : {len(ind)}")
    print(f"  ...AND carrying an MPIA delta   : {len(both)}")
    print(f"  lines lost to a missing MPIA node: {int(ind.delta_nlte.isna().sum())}")
    print(f"  ...AND on Amarsi's GOLDEN list  : {len(golden)}   <- the like-for-like set")

    print("\nPER-LINE DECOMPOSITION on the common set (paired, never aggregate-minus-aggregate)")
    print(f"  {'term':<42}{'n':>5}{'median':>10}{'MAD':>9}{'SEM':>8}")
    rows = [
        ('A(1D-LTE)', both.a_1dlte),
        ('A(1D-NLTE) = 1D-LTE + delta_NLTE', both.a_1dnlte),
        ('A(3D-NLTE) = 1D-LTE + aberr', both.a_3dnlte),
        ('delta_NLTE  (1D-NLTE - 1D-LTE)', both.delta_nlte),
        ('NET 3D-NLTE (3D-NLTE - 1D-LTE)', both.aberr),
        ('3D-ATMOSPHERE (3D-NLTE - 1D-NLTE)', both.term_3d_atmosphere),
    ]
    stats = {}
    for name, v in rows:
        st = _stat(v.to_numpy())
        stats[name] = st
        if st['n']:
            print(f"  {name:<42}{st['n']:>5}{st['median']:>+10.4f}"
                  f"{st['mad']:>9.4f}{st['sem']:>8.4f}")

    print("\nSAME DECOMPOSITION, restricted to Amarsi's OWN golden lines")
    print(f"  {'term':<42}{'n':>5}{'median':>10}{'MAD':>9}{'SEM':>8}")
    gstats = {}
    for name, col in (('delta_NLTE  (1D-NLTE - 1D-LTE)', 'delta_nlte'),
                      ('NET 3D-NLTE (3D-NLTE - 1D-LTE)', 'aberr'),
                      ('3D-ATMOSPHERE (3D-NLTE - 1D-NLTE)', 'term_3d_atmosphere')):
        st = _stat(golden[col].to_numpy())
        gstats[name] = st
        if st['n']:
            print(f"  {name:<42}{st['n']:>5}{st['median']:>+10.4f}"
                  f"{st['mad']:>9.4f}{st['sem']:>8.4f}")

    t3d = gstats['3D-ATMOSPHERE (3D-NLTE - 1D-NLTE)'] or stats['3D-ATMOSPHERE (3D-NLTE - 1D-NLTE)']
    if not t3d['n']:
        t3d = stats['3D-ATMOSPHERE (3D-NLTE - 1D-NLTE)']
        print("\n  (no golden-list overlap — falling back to the in-domain set)")
    magic = -FE_1D3D_SOLAR_OFFSET
    disc = t3d['median'] - magic
    print(f"\nTHE COMPARISON")
    print(f"  Magic-2013 as applied by gold      : {magic:+.4f} dex "
          f"(FE_1D3D_SOLAR_OFFSET)")
    print(f"  MEASURED 3D-atmosphere term        : {t3d['median']:+.4f} dex "
          f"(n={t3d['n']}, SEM {t3d['sem']:.4f})")
    print(f"  discrepancy                        : {disc:+.4f} dex")
    n_sem = abs(disc) / t3d['sem'] if t3d['sem'] else float('inf')
    print(f"  significance                       : {n_sem:.0f} x SEM")

    OUT.mkdir(parents=True, exist_ok=True)
    keep = ['wavelength_air_A', 'elo_eV', 'eup_eV', 'loggf', 'ew_mA', 'network',
            'a_1dlte', 'delta_nlte', 'a_1dnlte', 'aberr', 'a_3dnlte',
            'term_3d_atmosphere', 'in_domain']
    d[keep].to_csv(OUT / 'rya819_decomposition_per_line.csv', index=False)
    (OUT / 'rya819_decomposition_summary.json').write_text(json.dumps({
        'ticket': 'RYA-819',
        'line_set': 'RYA-817 in-domain Fe I VIS (Amarsi training domain) with an MPIA node',
        'n_common_in_domain': int(len(both)),
        'n_on_amarsi_golden_list': int(len(golden)),
        'terms_in_domain': stats,
        'terms_golden_lines': gstats,
        'magic2013_as_applied_dex': magic,
        'measured_3d_atmosphere_dex': t3d['median'],
        'discrepancy_dex': disc,
        'significance_sem': n_sem,
        'caveat': ('Magic-2013 is a <3D>-minus-1D shift computed in LTE; this term holds '
                   'NLTE in BOTH sides. They are not identical quantities and the '
                   'residual difference is named in the verdict, not glossed.'),
    }, indent=2, default=float))
    print(f"\n[out] {OUT}/rya819_decomposition_summary.json")


if __name__ == '__main__':
    main()
