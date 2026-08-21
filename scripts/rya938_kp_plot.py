#!/usr/bin/env python3
"""RYA-938 figure: the registration defect, and what it did to the lines."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.rya938_kp_crosscheck import (  # noqa: E402
    read_kp1984, read_kurucz2005, read_iag, envelope_normalise)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kp1984", type=Path, required=True)
    ap.add_argument("--kurucz", type=Path, required=True)
    ap.add_argument("--iag", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    panels = [(6694.0, 6699.0, "Al I 6696.185 — RYA-929 measured depth 0.008 here"),
              (7697.0, 7701.0, "K I 7698.964 — RYA-929 measured depth 0.039 here"),
              (6865.0, 6886.0, "O2 B 6867–6884 — the band 1984 does not correct")]
    fig, axes = plt.subplots(len(panels), 1, figsize=(11, 9))
    for ax, (lo, hi, title) in zip(axes, panels):
        wk, fk = read_kp1984(args.kp1984, lo, hi)
        wi, fi = read_iag(args.iag, lo, hi)
        wa, fa = read_kurucz2005(args.kurucz, lo, hi, assume_vacuum=False)
        wv, fv = read_kurucz2005(args.kurucz, lo, hi, assume_vacuum=True)
        ax.plot(wk, fk, lw=0.8, color="#b03030", label="Kitt Peak 1984 (air)")
        ax.plot(wi, fi, lw=0.8, color="#7a7a7a", label="IAG/Baker 2020 telluric-free")
        ax.plot(wa, envelope_normalise(wa, fa), lw=0.9, color="#d08a00", ls="--",
                label="Kurucz 2005 read AS AIR (RYA-929)")
        ax.plot(wv, envelope_normalise(wv, fv), lw=0.9, color="#20609a",
                label="Kurucz 2005 vacuum→air (RYA-938)")
        ax.set_title(title, fontsize=9)
        ax.set_xlim(lo, hi); ax.set_ylim(-0.05, 1.15)
        ax.grid(alpha=0.2); ax.legend(fontsize=7, loc="lower left", ncol=2)
    axes[-1].set_xlabel("air wavelength (Å)")
    fig.suptitle("RYA-938 — Kurucz 2005 is on a VACUUM grid; reading it as air "
                 "moved every line ~1.7 Å", fontsize=11)
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=140)
    print(args.out)


if __name__ == "__main__":
    main()
