#!/usr/bin/env python3
"""RYA-763: is the GES-vs-atom level disagreement an OFF-BY-ONE, or two different atoms?

A systematic index-origin difference (0-based vs 1-based) would produce exactly the
signature the mapping test reports. Ruling it out is cheap and decisive: scan an offset
range and see whether agreement peaks sharply at some shift. If the best offset is 0 and
agreement stays near half, the atoms genuinely differ.

Also checks agreement vs level index, because two atoms that share a ground-state
ordering and diverge higher up would show a clean decline rather than noise.
"""
from __future__ import annotations

import os
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, os.environ.get("ISPEC_DIR", "/srv/codex/engines/ispec_src"))

from scripts.rya763_level_mapping import read_labels, UNSET  # noqa: E402

EL, ION, LO, HI = "Ti", "I", 6910.0, 9199.9


def core(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9]", "", str(s)).lower()


def main() -> None:
    lab = read_labels(EL).set_index("index")
    from pipeline.abundances_derive import _load_synth_resources
    ll, _, _ = _load_synth_resources()
    w = np.asarray(ll["wave_A"], dtype=float)
    els = np.asarray([str(x).strip() for x in ll["element"]])
    want = f"{EL} {1 if ION.upper() == 'I' else 2}".upper()
    m = np.array([e.upper().startswith(want) for e in els]) & (w >= LO) & (w <= HI)

    pairs = []
    for i in np.where(m)[0]:
        for side in ("low", "up"):
            ri = str(ll[f"nlte_level_{side}"][i]).strip()
            gl = str(ll[f"nlte_label_{side}"][i]).strip()
            if ri in UNSET or gl in UNSET:
                continue
            try:
                pairs.append((int(float(ri)), gl))
            except ValueError:
                pass
    print(f"  {len(pairs)} usable (index, label) endpoints for {EL} {ION}")

    print(f"\n{'offset':>7} {'agree':>7} {'of':>7} {'rate':>8}")
    best = (None, -1.0)
    for off in range(-3, 4):
        ok = tot = 0
        for k, gl in pairs:
            kk = k + off
            if kk not in lab.index:
                continue
            tot += 1
            c, g = core(lab.loc[kk, "term"]), core(gl)
            if g.startswith(c) or c.startswith(g[:len(c)]):
                ok += 1
        rate = ok / tot if tot else 0.0
        flag = "  <-- best" if rate > best[1] else ""
        if rate > best[1]:
            best = (off, rate)
        print(f"{off:>7} {ok:>7} {tot:>7} {rate:>8.3f}{flag}")

    print(f"\n  best offset {best[0]} at {best[1]:.1%} agreement")
    if best[0] == 0 and best[1] < 0.9:
        print("  => NOT an off-by-one. The index origins line up and agreement is still")
        print("     poor, so GES and the Engine-A atom are DIFFERENT LEVEL SETS.")
    elif best[0] != 0 and best[1] > 0.9:
        print(f"  => OFF-BY-{best[0]}. The disagreement is an index-origin convention,")
        print("     not a physics difference.")

    # agreement vs index — do they share the low-lying levels and diverge above?
    print("\n  agreement by index decile (offset 0):")
    ks = np.array([k for k, _ in pairs])
    for q0, q1 in zip(range(0, 100, 20), range(20, 120, 20)):
        lo_k, hi_k = np.percentile(ks, q0), np.percentile(ks, q1)
        sel = [(k, g) for k, g in pairs if lo_k <= k <= hi_k and k in lab.index]
        if not sel:
            continue
        ok = sum(1 for k, g in sel
                 if core(g).startswith(core(lab.loc[k, "term"]))
                 or core(lab.loc[k, "term"]).startswith(core(g)[:len(core(lab.loc[k, 'term']))]))
        print(f"    index {int(lo_k):4d}-{int(hi_k):4d}: {ok:5d}/{len(sel):5d}"
              f"  = {ok/len(sel):.3f}")


if __name__ == "__main__":
    main()
