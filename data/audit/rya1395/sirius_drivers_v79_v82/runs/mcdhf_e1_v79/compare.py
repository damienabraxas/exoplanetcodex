import csv
import math
import statistics
import sys

KEY_WAVELENGTHS = (130.4, 135.6, 615.8, 645.6, 777.4, 844.8, 926.6, 1128.9, 1316.4)


def load(path):
    result = {}
    with open(path, newline="") as source:
        for row in csv.DictReader(source):
            key = (row["lower"], row["upper"])
            if key in result:
                raise ValueError(f"{path}: duplicate transition {key}")
            result[key] = row
    if not result:
        raise ValueError(f"{path}: empty transition table")
    return result


def gauge_gap(row):
    scaled_c = float(row["loggf_C"]) + 2 * math.log10(
        float(row["sigma_model"]) / float(row["sigma_nist"])
    )
    return abs(float(row["loggf_B"]) - scaled_c)


def summary(values):
    finite = [value for value in values if math.isfinite(value)]
    return f"{statistics.median(finite):+.4f} ({len(finite)}/{len(values)})" if finite else "—"


def main(n5_path, n6_path):
    n5, n6 = load(n5_path), load(n6_path)
    common = n5.keys() & n6.keys()
    print("# O I n5 to n6 · HOLD / report only")
    print(f"\nShared transitions: {len(common)}. Only n5: {len(n5.keys() - common)}. Only n6: {len(n6.keys() - common)}.")
    for layer, missing in [("n6", n5.keys() - common), ("n5", n6.keys() - common)]:
        for lower, upper in sorted(missing):
            print(f"- Missing from {layer}: {lower} → {upper}")
    print("\nClass statistics are paired medians in dex with finite/total counts.")
    print("\n| Li class | n5 B−Li | n6 B−Li | n5 scaled gauge gap | n6 scaled gauge gap | n6−n5 B |")
    print("|---|---:|---:|---:|---:|---:|")
    for acc in "ABCDE":
        pairs = [(n5[key], n6[key]) for key in common if n5[key]["li_acc"] == n6[key]["li_acc"] == acc]
        stats = [
            summary([float(a["d_B"]) for a, b in pairs]),
            summary([float(b["d_B"]) for a, b in pairs]),
            summary([gauge_gap(a) for a, b in pairs]),
            summary([gauge_gap(b) for a, b in pairs]),
            summary([float(b["loggf_B"]) - float(a["loggf_B"]) for a, b in pairs]),
        ]
        print("| " + " | ".join([acc, *stats]) + " |")
    changed_classes = [key for key in common if n5[key]["li_acc"] != n6[key]["li_acc"]]
    print(f"\nShared transitions with inconsistent Li class: {len(changed_classes)}.")
    print("\n| Vacuum nm | Lower | Upper | n5 B−Li | n6 B−Li | n6−n5 B | n6 scaled gauge gap |")
    print("|---:|---|---|---:|---:|---:|---:|")
    for key in sorted(common, key=lambda key: (float(n5[key]["lam_vac_nm"]), key)):
        a, b = n5[key], n6[key]
        lam = float(a["lam_vac_nm"])
        if b["li_acc"] and any(abs(lam - center) < 1 for center in KEY_WAVELENGTHS):
            delta = float(b["loggf_B"]) - float(a["loggf_B"])
            print(f"| {lam:.3f} | {key[0]} | {key[1]} | {float(a['d_B']):+.4f} | {float(b['d_B']):+.4f} | {delta:+.4f} | {gauge_gap(b):.4f} |")


if __name__ == "__main__":
    main(*sys.argv[1:])
