import itertools
import json
from pathlib import Path

import numpy as np


def lu_metrics(matrix):
    work = matrix.copy()
    pivots = []
    multiplier = 0.0
    for j in range(len(work)):
        pivot = work[j, j]
        pivots.append(abs(pivot))
        if pivot == 0:
            return 0.0, float("inf")
        factors = work[j + 1:, j] / pivot
        multiplier = max(multiplier, *abs(factors), 0.0)
        work[j + 1:, j + 1:] -= factors[:, None] * work[j, j + 1:]
    return min(pivots), multiplier


def analyze(entries):
    labels = sorted({left for left, _, _ in entries})
    assert labels == sorted({right for _, right, _ in entries})
    values = {(left, right): value for left, right, value in entries}
    assert len(values) == len(entries) == len(labels) ** 2
    matrix = np.array([[values[left, right] for right in labels] for left in labels])
    scores = {}
    for rows, cols in itertools.product(itertools.permutations(range(len(labels))), repeat=2):
        permuted = matrix[np.ix_(rows, cols)]
        singular = np.linalg.svd(permuted, compute_uv=False)
        np.testing.assert_allclose(singular, np.linalg.svd(matrix, compute_uv=False), rtol=1e-9)
        scores[rows, cols] = (*lu_metrics(permuted), *lu_metrics(np.linalg.inv(permuted)[::-1, ::-1]))
    native = tuple(range(len(labels)))
    reverse = native[::-1]
    max_pivot = max(scores, key=lambda order: (scores[order][0], -scores[order][1]))
    min_multiplier = min(scores, key=lambda order: (scores[order][1], -scores[order][0]))
    inverse_min_multiplier = min(scores, key=lambda order: (scores[order][3], -scores[order][2]))
    print(f"  SVD condition number {np.linalg.cond(matrix):.5g}")
    for name, order in (
        ("native", (native, native)), ("joint reverse", (reverse, reverse)),
        ("max S pivot", max_pivot), ("min S multiplier", min_multiplier),
        ("min reversed inverse multiplier", inverse_min_multiplier),
    ):
        axes = [" ".join(labels[i] for i in axis) for axis in order]
        a, b, c, d = scores[order]
        print(f"  {name:32s} {axes[0]:12s} / {axes[1]:12s} | S {a:.5g} {b:.5g} | reverse inv(S) {c:.5g} {d:.5g}")


if __name__ == "__main__":
    for layer, blocks in json.loads(Path(__file__).with_name("overlaps.json").read_text()).items():
        assert set(blocks) == {"f", "f-"}
        for symmetry, entries in sorted(blocks.items()):
            print(layer, symmetry)
            analyze(entries)
