"""How a gf source's per-line errors correlate across one product's lines (RYA-1233).

The registry is data/reference/gf_error_model.csv. Ryan, 2026-10-08: uncertainty must be
computed the correct way per element, even where the right way differs. A source listed
`scale_only` contributes its per-line sigma on the diagonal and only its absolute scale
off it -- so the term averages down over lines as the literature says it does (Amarsi &
Asplund 2017 for Garz). An UNLISTED source keeps the conservative fully-correlated rule,
and is named UNREVIEWED in the covariance note so the gap is visible on the product.
"""
from __future__ import annotations

import csv
from pathlib import Path

REGISTRY = Path(__file__).resolve().parents[1] / "data" / "reference" / "gf_error_model.csv"


def load() -> dict:
    with REGISTRY.open(newline="", encoding="utf-8") as fh:
        rows = csv.DictReader(l for l in fh if not l.startswith("#"))
        return {r["lab_source_tag"]: {**r, "scale_sigma_dex": float(r["scale_sigma_dex"] or 0)}
                for r in rows}


#: The registry key for a line graded by its NIST ASD accuracy class and carrying no lab tag.
NIST_KEY = "NIST_ASD"


def source_key(tag: str, src: str = "", nist_grade: str = "") -> str:
    """The registry key of one line's gf source: its lab tag, else NIST_ASD for a NIST-graded
    line (graded by class, no tag), else '' (unclassifiable -> UNREVIEWED)."""
    if tag:
        return tag
    if nist_grade or str(src).startswith("NIST ASD"):
        return NIST_KEY
    return ""


def covariance(sig: list[float], src: list[str], tags: list[str]) -> tuple[list, str, list]:
    """(covariance dex^2, the sentence that says how, unreviewed source names)."""
    model = load()
    n = len(sig)

    tags = [source_key(t, s) for t, s in zip(tags, src)]

    def kind(i):
        m = model.get(tags[i])
        return (m["correlation"], m["group"], m["scale_sigma_dex"]) if m else (
            "fully_correlated", "src:" + src[i], None)

    k = [kind(i) for i in range(n)]
    cov = []
    for a in range(n):
        row = []
        for b in range(n):
            if a == b:
                row.append(sig[a] ** 2)
            elif k[a][1] != k[b][1]:
                row.append(0.0)
            elif k[a][0] == "scale_only":
                s = k[a][2]
                row.append(min(s, sig[a]) * min(s, sig[b]))
            else:
                row.append(sig[a] * sig[b])
        cov.append(row)
    reviewed = sorted({f"{tags[i]}: {model[tags[i]]['correlation']} "
                       f"(scale {model[tags[i]]['scale_sigma_dex']:.4f} dex; "
                       f"{model[tags[i]]['citation']})" for i in range(n) if tags[i] in model})
    unreviewed = sorted({tags[i] or src[i] for i in range(n) if tags[i] not in model})
    note = "gf covariance per data/reference/gf_error_model.csv"
    if reviewed:
        note += " -- " + "; ".join(reviewed)
    if unreviewed:
        note += (" -- UNREVIEWED (default fully correlated, classify in the registry): "
                 + ", ".join(unreviewed))
    return cov, note, unreviewed
