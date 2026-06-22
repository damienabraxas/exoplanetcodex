"""
scripts/resolve_camn_stops_rya411.py
====================================
RYA-411 — resolve the RYA-410 Ca/Mn NLTE STOPs.

Part A (Mn): the RYA-410 cross-check STOP (+0.018 vs MPIA +0.107) had TWO causes, both
found by probe (not a single "HFS bug"):
  1. LINE-SET mismatch: the harness DIAG_WAVES used the low-excitation, strongly-HFS
     triplet 6013/6016/6021 (EP 3.07); the MPIA grid +0.107 is on high-excitation lines
     4998/6304/6306/6867 (EP 4.4-5.2), of which 6304/6306/6867 are ~single-component.
  2. HFS-collapse: _resolve_atomic gf-SUMS the HFS components of a feature into one line,
     which over-saturates the strong low-EP triplet and suppresses its NLTE delta.
This module fixes BOTH: an HFS-RESOLVED synthesis (each hyperfine component emitted as its
own line, sharing the feature's NLTE level — pysme_nlte._synth_ew RYA-411 path), derived on
the SAME lines MPIA used (apples-to-apples cross-check). validate-don't-tune: reproduce the
MPIA value by fixing the mechanism + line set, never tune to +0.107.

Part B (Ca): same-lines probe + Amarsi-2020-vs-Mashonkina-2017 model-atom adjudication.

    python -m scripts.resolve_camn_stops_rya411 --part mn
    python -m scripts.resolve_camn_stops_rya411 --part mn-demo
    python -m scripts.resolve_camn_stops_rya411 --part ca
"""
from __future__ import annotations
import argparse, sys, warnings
from pathlib import Path
import numpy as np
warnings.filterwarnings('ignore')

_REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO))
from config.constants import ISPEC_DIR, SOLAR_ASPLUND2021          # noqa: E402
from pipeline import pysme_nlte as P                                # noqa: E402

# the lines the MPIA grid uses for its registered +0.107 (the cross-check reference set)
MN_MPIA_LINES = [4998.126, 6304.906, 6306.342, 6867.11]
# the low-EP strongly-HFS triplet the RYA-410 harness picked (the mechanism demo)
MN_HFS_TRIPLET = [6013.51, 6016.67, 6021.80]

_cache = {}


def _ges():
    if 'll' not in _cache:
        sys.path.insert(0, str(ISPEC_DIR)); import ispec
        import pipeline.abundances_derive as ad
        ll = ispec.read_atomic_linelist(ad._SYNTH_LINELIST_FILE)
        _cache['ll'] = ll
        _cache['notes'] = np.array([str(x) for x in ll['element']])
        _cache['w'] = np.asarray(ll['wave_A'], float)
    return _cache['ll'], _cache['notes'], _cache['w']


def resolve_hfs(element, feature_wl, ion='1', dw=0.25, ep_tol=0.03):
    """Dominant-EP-cluster HFS components (NOT gf-summed): returns
    (gf_sum_loggf, EP, Eup, vdW, components=[(wl, gflog), ...])."""
    ll, notes, w = _ges()
    m = np.where((notes == f'{element} {ion}') & (np.abs(w - feature_wl) < dw))[0]
    if len(m) == 0:
        raise ValueError(f"{element} {feature_wl}: not in GES (+/-{dw} A)")
    eps = np.array([float(ll['lower_state_eV'][i]) for i in m])
    gfs = np.array([10 ** float(ll['loggf'][i]) for i in m])
    vws = np.array([float(ll['waals'][i]) for i in m])
    ws = w[m]
    order = np.argsort(eps); groups, cur = [], [order[0]]
    for k in order[1:]:
        if eps[k] - eps[cur[-1]] <= ep_tol:
            cur.append(k)
        else:
            groups.append(cur); cur = [k]
    groups.append(cur)
    g = max(groups, key=lambda idx: gfs[idx].sum())
    comps = sorted([(float(ws[i]), float(np.log10(gfs[i]))) for i in g])
    ep = float(np.median(eps[g])); vw = float(np.median(vws[g]))
    return float(np.log10(gfs[g].sum())), ep, ep + 12398.42 / feature_wl, (vw if vw > 0 else 0.0), comps


def _build_lines(element, feature_waves, hfs=True):
    """Build pysme_nlte line tuples for the features; HFS-resolved (10-tuple with
    component list) or gf-summed (9-tuple)."""
    lines = []
    for fw in feature_waves:
        gf, ep, eup, vw, comps = resolve_hfs(element, fw)
        tl, tu, jl, ju = P.auto_labels(element, ep, eup)
        base = (fw, gf, ep, jl, eup, ju, tl, tu, vw)
        lines.append(base + (comps,) if hfs else base)
    return lines, [len(resolve_hfs(element, fw)[4]) for fw in feature_waves]


def _mpia_ref(element, waves):
    from pipeline.nlte_corrections import _mpia_element_delta
    return {w: _mpia_element_delta(element, float(w), 5772, 4.44, 0.0) for w in waves}


def part_mn(demo=False):
    P._A_SUN.setdefault('Mn', SOLAR_ASPLUND2021['Mn'])
    star = {'teff': 5772, 'logg': 4.44, 'feh': 0.0, 'vmic': 1.0}
    if demo:
        print("\n== Mn MECHANISM DEMO — the RYA-410 triplet 6013/6016/6021 (EP 3.07, strong HFS) ==")
        gfs, ncomp = _build_lines('Mn', MN_HFS_TRIPLET, hfs=False)
        P.NLTE_LINES['Mn'] = gfs
        rs = P.nlte_delta('Mn', star=star)
        print(f"  HFS comps per feature: {dict(zip(MN_HFS_TRIPLET, ncomp))}")
        print(f"  gf-SUMMED (RYA-410 path):   per-line {{wl: round}} = "
              f"{{{', '.join(f'{k}:{v:+.3f}' for k,v in rs['per_line'].items())}}}  median {rs['delta_median']:+.3f}")
        hfsl, _ = _build_lines('Mn', MN_HFS_TRIPLET, hfs=True)
        P.NLTE_LINES['Mn'] = hfsl
        rh = P.nlte_delta('Mn', star=star)
        print(f"  HFS-RESOLVED (RYA-411 path): per-line "
              f"{{{', '.join(f'{k}:{v:+.3f}' for k,v in rh['per_line'].items())}}}  median {rh['delta_median']:+.3f}")
        print(f"  => HFS desaturation moves the strong triplet {rs['delta_median']:+.3f} -> {rh['delta_median']:+.3f}")
        return

    print("\n== Mn CROSS-CHECK — same lines as MPIA (4998/6304/6306/6867), HFS-resolved ==")
    ref = _mpia_ref('Mn', MN_MPIA_LINES)
    lines, ncomp = _build_lines('Mn', MN_MPIA_LINES, hfs=True)
    print(f"  HFS comps per feature: {dict(zip(MN_MPIA_LINES, ncomp))}")
    P.NLTE_LINES['Mn'] = lines
    r = P.nlte_delta('Mn', star=star)
    print(f"  {'line':>10} {'PySME-HFS':>10} {'MPIA':>8} {'diff':>8}")
    diffs = []
    for w in MN_MPIA_LINES:
        d = r['per_line'][w]; m = ref[w]; diffs.append(d - m)
        print(f"  {w:>10} {d:>+10.4f} {m:>+8.4f} {d-m:>+8.4f}")
    med = r['delta_median']; mref = float(np.nanmedian(list(ref.values())))
    ok = abs(med - mref) <= 0.04
    print(f"  MEDIAN PySME-HFS {med:+.4f} vs MPIA {mref:+.4f} (diff {med-mref:+.4f}, tol 0.04) "
          f"-> {'PASS — register Amarsi Mn (HFS-resolved)' if ok else 'STOP'}")
    return dict(median=med, mpia=mref, passed=ok, per_line=r['per_line'], ref=ref)


# Ca: the cleanly grid-matchable lines shared with the MPIA grid (the same-lines probe);
# and the lines the RYA-410 harness actually used (to reproduce its +0.064 and decompose).
CA_MPIA_CLEAN = [5867.57, 6166.44]
CA_RYA410 = [6166.439, 6169.563, 6455.598]


def part_ca():
    P._A_SUN.setdefault('Ca', SOLAR_ASPLUND2021['Ca'])
    star = {'teff': 5772, 'logg': 4.44, 'feh': 0.0, 'vmic': 1.0}
    print("\n== Ca SAME-LINES probe — grid-matchable lines shared with MPIA (5867/6166) ==")
    ref = _mpia_ref('Ca', CA_MPIA_CLEAN)
    lines, _ = _build_lines('Ca', CA_MPIA_CLEAN, hfs=False)
    P.NLTE_LINES['Ca'] = lines
    r = P.nlte_delta('Ca', star=star)
    print(f"  {'line':>9} {'PySME':>9} {'MPIA':>9} {'diff':>8}")
    for w in CA_MPIA_CLEAN:
        print(f"  {w:>9} {r['per_line'][w]:>+9.4f} {ref[w]:>+9.4f} {r['per_line'][w]-ref[w]:>+8.4f}")
    same_med = float(np.median(list(r['per_line'].values())))
    print(f"  same-lines median PySME {same_med:+.4f} vs MPIA {float(np.median(list(ref.values()))):+.4f}")

    print("\n== Ca RYA-410 lines (6166/6169/6455) — reproduce the +0.064 ==")
    lines2, _ = _build_lines('Ca', CA_RYA410, hfs=False)
    P.NLTE_LINES['Ca'] = lines2
    r2 = P.nlte_delta('Ca', star=star)
    for w in CA_RYA410:
        print(f"  {w:>9} {r2['per_line'][w]:>+9.4f}")
    print(f"  RYA-410-lines median {r2['delta_median']:+.4f} (RYA-410 reported +0.064)")
    print(f"\n  DECOMPOSE: MPIA grid +0.012 (8 lines) | RYA-410 +0.064 (6166/6169/6455) | "
          f"same-lines PySME {same_med:+.4f}")
    return dict(same_lines_pysme=same_med, mpia_ref=ref, rya410_repro=r2['delta_median'])


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument('--part', choices=['mn', 'mn-demo', 'ca'], required=True)
    a = p.parse_args(argv)
    if a.part == 'mn-demo':
        part_mn(demo=True)
    elif a.part == 'mn':
        part_mn(demo=False)
    elif a.part == 'ca':
        part_ca()


if __name__ == '__main__':
    main()
