#!/usr/bin/env python
"""RYA-1282: O 777-triplet departures of each grid against the Amarsi+2020 GALAH O grid.

Amarsi+2020 (A&A 642, A62; Balder = modified Multi3D, Blue EOS + ~1e6 atomic / 1e8 molecular
background lines, Amarsi+2018 O atom with Barklem H collisions) is the EXTERNAL referee: an
independent code, atom and opacity set. Its abundance axis is [O/Fe] against Grevesse+2007
(A(O)_sun = 8.66), so A(O) = 8.66 + [Fe/H] + abund. Models: standard MARCS, plane-parallel
t01 for log g > 3.5, spherical t02 for log g <= 3.5 (paper Sect. 2.2) - only our nodes on
the SAME MARCS model are compared.

Engines: --gerber (the published deck) and any number of --npz name=DIR
(pipeline.o_departure_grid output). Our grids are interpolated linearly in A(O) (log b), and
Amarsi's tau scale onto the MARCS log tau500 of the node. Output: one CSV row per
(model, Amarsi abund, engine) with the triplet levels' max |dlog b| over log tau -3..0 and
the signed dlog b of 3s 5S and 3p 5P(J=3) at log tau -1.
"""
import argparse
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np                                                          # noqa: E402

A_SUN_G07 = 8.66
TRIPLET = [5, 7, 8, 9]          # 3s 5S2, 3p 5P1,2,3 - same 0-based indices in O41f and Amarsi's atom


def model_convention_ok(mid: str, logg: float) -> bool:
    """Same MARCS model Amarsi used: pp t01 (log g > 3.5) or sph t02 (log g <= 3.5), standard alpha."""
    pp = logg > 3.5
    return (mid.startswith('p') and '_t01_' in mid) if pp else (mid.startswith('s') and '_t02_' in mid)


def gerber_index(aux):
    idx = {}
    for line in open(aux):
        if line.startswith("'"):
            f = line.split()
            idx[(f[0].strip("'"), round(float(f[7]), 2))] = int(f[8])
    return idx


def at_abundance(Avals, B, A):
    """log b interpolated linearly in A(O) between the bracketing nodes (None if outside)."""
    Avals = np.asarray(Avals, float)
    if A < Avals.min() - 1e-6 or A > Avals.max() + 1e-6:
        return None
    j = int(np.clip(np.searchsorted(Avals, A) , 1, len(Avals) - 1))
    w = (A - Avals[j - 1]) / (Avals[j] - Avals[j - 1])
    lb = (1 - w) * np.log10(B[j - 1]) + w * np.log10(B[j])
    return lb


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--amarsi', required=True)
    p.add_argument('--gerber-bin'); p.add_argument('--gerber-aux')
    p.add_argument('--npz', action='append', default=[], help='name=DIR')
    p.add_argument('--models-from', required=True, help='name of the --npz set whose models define the comparison')
    p.add_argument('--abund', default='-0.5,0,0.5')
    p.add_argument('--out', required=True)
    a = p.parse_args()
    from pipeline.nlte_bfactor_synth import PySMEGrid
    from pipeline.o_departure_grid import read_record
    g = PySMEGrid(a.amarsi)
    tau = g.get('tau'); T, G_, M = g.teff, g.logg, g.feh
    sets = dict(s.split('=', 1) for s in a.npz)
    gidx = gerber_index(a.gerber_aux) if a.gerber_aux else None
    rows = []
    for f in sorted(Path(sets[a.models_from]).glob('*.npz')):
        z0 = np.load(f)
        mid = str(z0['model_id']); teff, logg, feh = float(z0['teff']), float(z0['logg']), float(z0['feh'])
        if not model_convention_ok(mid, logg):
            continue
        it, ig, im = (np.where(np.isclose(T, teff))[0], np.where(np.isclose(G_, logg))[0],
                      np.where(np.isclose(M, feh))[0])
        if not (len(it) and len(ig) and len(im)):
            continue
        lt = np.asarray(z0['logtau'], float)
        lta = np.log10(np.asarray(tau[im[0], ig[0], it[0]], float))
        ph = (lt >= -3) & (lt <= 0); k1 = int(np.argmin(abs(lt + 1)))
        for ab in (float(x) for x in a.abund.split(',')):
            key = f't{int(teff)}_g{logg:+.2f}_m{feh:+.2f}_a{ab:+.2f}'
            if key not in g._index:
                key = key.replace('a+0.00', 'a-0.00')
            if key not in g._index:
                continue
            ba = g.get(key)                                       # (nlev, ndep)
            if not np.all(np.isfinite(ba[TRIPLET])) or np.any(ba[TRIPLET] <= 0):
                continue
            lba = np.array([np.interp(lt, lta, np.log10(ba[k])) for k in TRIPLET])
            A = A_SUN_G07 + feh + ab
            cands = {}
            for name, d in sets.items():
                fz = Path(d) / f.name
                if fz.exists():
                    z = np.load(fz)
                    ok = z['ok']
                    lb = at_abundance(z['abundances'][ok], z['b'][ok], A)
                    if lb is not None:
                        cands[name] = lb[TRIPLET]
            if gidx is not None:
                lo = round(np.floor((A - feh) * 10 + 1e-6) / 10 + feh, 2); hi = round(lo + 0.1, 2)
                if (mid, lo) in gidx and (mid, hi) in gidx:
                    Bs = [read_record(a.gerber_bin, gidx[(mid, x)])[2] for x in (lo, hi)]
                    cands['gerber'] = at_abundance([lo, hi], np.array(Bs), A)[TRIPLET]
            for name, lb in cands.items():
                d = lb - lba
                rows.append(dict(model=mid, teff=teff, logg=logg, feh=feh, o_fe=ab, A_O=round(A, 3), engine=name,
                                 max_abs_dlogb=round(float(np.abs(d[:, ph]).max()), 5),
                                 dlogb_3s5S_tau_m1=round(float(d[0, k1]), 5), dlogb_3p5P3_tau_m1=round(float(d[3, k1]), 5),
                                 logb_3s5S_tau_m1_amarsi=round(float(lba[0, k1]), 5)))
    with open(a.out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else ['model'])
        w.writeheader(); w.writerows(rows)
    print(f'{len(rows)} rows -> {a.out}')
    if rows:
        import collections
        agg = collections.defaultdict(list)
        for r in rows:
            agg[(r['engine'], r['feh'])].append(r['max_abs_dlogb'])
        for e in sorted({r['engine'] for r in rows}):
            print(e, ' '.join(f"[Fe/H]{fe:+.2f}: {np.median(v):.3f} (n{len(v)})" for (ee, fe), v in sorted(agg.items()) if ee == e))


if __name__ == '__main__':
    main()
