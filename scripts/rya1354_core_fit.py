"""RYA-1354 v69 (Part B): fit the l-dependent model potential of a one-valence CORE ion exactly as v63 did for Mg II
(rya1354_onevalence_model.fit on NIST J-averaged terms, n <= 9), and write summary.json with potential_fit, which
rya1354_twovalence_model.load_core reads. No line referee at this level (the sister-ion test is on the two-valence ion).
Usage: rya1354_core_fit.py 'Al III'"""
import json
import sys

import rya1354_onevalence_model as M


def main(name):
    sp = M.use(name)
    lev = M.nist_levels()
    M.NMIN = M.nmin_by_l(lev)
    params, rep = M.fit(lev)
    # J-resolved fine structure of the core (used by the ion xi calibration downstream) and fit quality
    terms = M.nist_terms(lev)
    M.OUT.mkdir(parents=True, exist_ok=True)
    s = dict(schema='rya1354.core-fit.v1', species=name, status=M.STATUS, alpha_d=M.ALPHA_D, nmin=M.NMIN,
             n_levels=len(lev), n_terms=len(terms), potential_fit=rep)
    (M.OUT / 'summary.json').write_text(json.dumps(s, indent=1) + '\n')
    return s


if __name__ == '__main__':
    for n in sys.argv[1:]:
        s = main(n)
        print(n, {k: (v['n_terms'], round(v['max_rel_err_start'], 4), round(v['max_rel_err_fit'], 5)) for k, v in s['potential_fit'].items()})
