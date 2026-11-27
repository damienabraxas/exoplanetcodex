"""RYA-1395 v83 stage 0: native vs GLOBAL reverse orbital order, every symmetry.

Matrix only. Reads the pre-rotation overlap matrices that rbiotransform prints
("Overlap matrix before orbital rotations") and, for every square symmetry
block, reports the no-pivot LU minimum pivot and largest elimination
multiplier in native (n ascending) and reversed (n descending, both axes)
order, plus the order-invariant SVD condition number.

GRASP's KAPDATA accepts only all-normal or all-reversed order across every
kappa, so the reversal is applied to every block at once.

Usage: rya1395_order_stage0_v83.py LABEL=path/to/rbiotransform.stdout ...
"""
import re
import sys
from collections import defaultdict

import numpy as np

ROW = re.compile(r"<\s*(\d+)([a-z]-?)\s*\|\s*(\d+)([a-z]-?)\s*>\s*=\s*(\S+)")


def overlaps(path):
    blocks = defaultdict(dict)
    inside = False
    for line in open(path):
        if "Overlap matrix before orbital rotations" in line:
            inside = True
            continue
        if "Overlap matrix after orbital rotations" in line:
            break
        if inside and (m := ROW.search(line)):
            n1, s1, n2, s2, value = m.groups()
            assert s1 == s2, line
            blocks[s1][int(n1), int(n2)] = float(value.replace("D", "E"))
    assert blocks, f"no pre-rotation overlaps in {path}"
    return blocks


def lu(matrix):
    work = matrix.copy()
    pivots, mult = [], 0.0
    for j in range(len(work)):
        pivot = work[j, j]
        pivots.append(abs(pivot))
        if pivot == 0:
            return 0.0, float("inf")
        factors = work[j + 1:, j] / pivot
        if len(factors):
            mult = max(mult, float(np.max(np.abs(factors))))
        work[j + 1:, j + 1:] -= np.outer(factors, work[j, j + 1:])
    return min(pivots), mult


def main(args):
    print(f"{'layer':6} {'sym':4} {'dim':>3} {'cond':>10} | "
          f"{'native minpiv':>13} {'native mult':>11} | {'rev minpiv':>11} {'rev mult':>11} | verdict")
    for arg in args:
        label, path = arg.split("=", 1)
        for sym, entries in sorted(overlaps(path).items()):
            rows = sorted({a for a, _ in entries})
            cols = sorted({b for _, b in entries})
            if rows != cols or len(entries) != len(rows) ** 2:
                print(f"{label:6} {sym:4} rectangular {len(rows)}x{len(cols)}: skipped")
                continue
            mat = np.array([[entries[a, b] for b in cols] for a in rows])
            nat = lu(mat)
            rev = lu(mat[::-1, ::-1])
            cond = np.linalg.cond(mat)
            worse = rev[0] < nat[0] * 0.1 or rev[1] > max(10 * nat[1], 10)
            verdict = "REVERSE WORSE" if worse else ("reverse better" if rev[1] < nat[1] else "same/ok")
            print(f"{label:6} {sym:4} {len(rows):>3} {cond:>10.3g} | {nat[0]:>13.3g} {nat[1]:>11.3g} | "
                  f"{rev[0]:>11.3g} {rev[1]:>11.3g} | {verdict}")


if __name__ == "__main__":
    main(sys.argv[1:])
