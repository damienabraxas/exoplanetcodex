"""Named laboratory-graded anchor pools — RYA-987.

The anchor is the set of lines whose gf comes from a primary laboratory measurement. It is
what sets the ABSOLUTE scale, and `pipeline.error_budget.zero_point_cap` computes the
accuracy floor from it.

🔴 THIS IS A REGISTRY, NOT A MEASUREMENT. Every value is read back out of the committed
band products; nothing is re-derived. RYA-987 is explicit that the anchor comes from
RYA-984's product and must not be re-measured — a second derivation of the same lines
would be a second home for the anchor (RYA-350/353/954), free to disagree with the
products the numbers were published from.

⚠️ THE ANCHOR IS STILL GROWING, so this is keyed by name and dated rather than written as
one global constant. The HARPS deep leg is owed on both arms, RYA-977 adds ~65 lines below
4200 A, and the UV/IR bands (RYA-847) are not folded in. Add a row; do not edit one.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
_PRODUCTS = ROOT / "data" / "results" / "band_products"
_CANONICAL = ROOT / "data" / "linelists" / "canonical_gf.csv"

#: Same wavelength key the selectors use, so "this line" means one thing project-wide.
_MATCH_TOL_A = 0.005


@dataclass(frozen=True)
class AnchorPool:
    name: str
    parts: tuple[str, ...]          # per-line product filenames that make up the anchor
    species: str
    note: str
    #: 🔴 RYA-1006 — WHAT THE ANCHOR'S SPECTRA WERE, BEFORE ANYONE FITTED THEM.
    #: `'native'` = the holding as it ships. Declared here so `load` can REFUSE a product
    #: conditioned some other way, instead of averaging it in because the filename matched.
    conditioning: str = "native"


ANCHORS: dict[str, AnchorPool] = {
    "rya984_graded_163": AnchorPool(
        name="rya984_graded_163",
        parts=("FeI_4200_6910_kpno_solar_atlas_SYNTH_1D-LTE_lines.csv",
               "FeI_4200_6910_kpno_solar_atlas_SYNTH_DEEPGRADED_1D-LTE_lines.csv"),
        species="Fe I",
        conditioning="native",
        note=("The VIS Kitt-Peak graded anchor as of RYA-984: 55 shallow graded lines "
              "(RYA-967, the EW-comparable set) + 108 deep graded lines (RYA-984, the set "
              "EW can never attempt). Combining them is valid — RYA-984 measured the "
              "deep-vs-shallow offset at p=0.149 with a flat excitation-potential trend, "
              "and the pooled within-set scatter 0.1990 matches the combined 0.1992, so "
              "the combination hides no structure. "
              "🔴 RYA-1006: BOTH PER-LINE PARTS WERE DESTROYED ON 2026-08-23 and must be "
              "regenerated before this anchor loads again. RYA-1000's `--local-renorm` "
              "runs overwrote them under these exact filenames (Kitt Peak 7.417 -> 7.337, "
              "HARPS 7.535 -> 7.339) because the stem carried no conditioning axis. The "
              "committed SUMMARY products were restored from git and read 7.461 (n=55) / "
              "7.417 (n=108) / 7.535 (n=109); the per-line CSVs are gitignored (RYA-469) "
              "and existed nowhere else. The conditioned copies are kept at "
              "data/audit/rya1006_preserved/ under _LOCALRENORM names."),
    ),
}


def _assert_conditioning(d: pd.DataFrame, *, part: str, pool: AnchorPool) -> None:
    """🔴 RYA-1006 — the anchor's own spectra, checked rather than assumed.

    THIS GUARD EXISTS BECAUSE THE FAILURE ALREADY HAPPENED. On 2026-08-23 two RYA-1000
    `--local-renorm` runs wrote `FeI_4200_6910_*_SYNTH_DEEPGRADED_*` — the exact filenames
    this module names — with every fit window divided by its own 95th percentile. The
    anchor's Kitt Peak value moved **7.417 -> 7.337** and NOTHING could tell: the stem did
    not carry the conditioning axis (fixed in `derive_band_products.conditioning_tag`) and
    the provenance file was BYTE-IDENTICAL. A filename matched, so a different measurement
    was loaded as the anchor.

    ⚠️ A MISSING COLUMN IS A REFUSAL, NOT A PASS. A product written before RYA-1006 cannot
    say what was done to its spectrum, and treating silence as `'native'` would re-admit
    exactly the corrupted files this was written for — they are silent too (RYA-833: an
    absence is a hypothesis, never a conclusion). Regenerate the part; the zero-point cap
    the anchor feeds is not worth guessing at.
    """
    if "observed_conditioning" not in d.columns:
        raise SystemExit(
            f"anchor {pool.name!r} part {part} carries no `observed_conditioning` column, "
            f"so it cannot say whether its spectra were conditioned before the fit. It "
            f"predates RYA-1006 and is UNVERIFIABLE — which is the state the RYA-1000 "
            f"`--local-renorm` overwrite left this very anchor in (Kitt Peak 7.417 -> "
            f"7.337 under this filename, with a byte-identical provenance). Regenerate it "
            f"with `scripts/derive_band_products.py` at RYA-1006 or later; refusing to "
            f"compute a laboratory zero point on a pool of unknown provenance.")
    got = sorted(set(d.observed_conditioning.dropna().astype(str)))
    unknown = int(d.observed_conditioning.isna().sum())
    if unknown:
        raise SystemExit(
            f"anchor {pool.name!r} part {part}: {unknown} of {len(d)} rows carry no "
            f"`observed_conditioning`. A blank is UNKNOWN, never 'native' (RYA-833).")
    if got != [pool.conditioning]:
        raise SystemExit(
            f"anchor {pool.name!r} declares conditioning {pool.conditioning!r} but part "
            f"{part} was measured on spectra conditioned {got}. This is a DIFFERENT "
            f"measurement wearing the anchor's filename — the RYA-1006 defect. Point the "
            f"anchor at the product that matches, or regenerate it; do not average them.")


def load(name: str) -> pd.DataFrame:
    """The anchor's per-line abundances, with each line's laboratory source and cited sigma.

    Loud on a missing part: an anchor silently short of its lines would shrink the
    statistical term and report a floor better than the one that was measured.
    """
    if name not in ANCHORS:
        raise SystemExit(f"unknown anchor {name!r}; known: {sorted(ANCHORS)}")
    pool = ANCHORS[name]
    frames = []
    for part in pool.parts:
        p = _PRODUCTS / part
        if not p.exists():
            raise SystemExit(
                f"anchor {name!r} names {part}, which is not in {_PRODUCTS}. The per-line "
                f"products are gitignored (RYA-469); stage them from the run that made "
                f"them rather than computing a cap on a partial anchor.")
        d = pd.read_csv(p)
        _assert_conditioning(d, part=part, pool=pool)
        frames.append(d[d.in_aggregate & d.abundance.notna()].assign(_part=part))
    a = pd.concat(frames, ignore_index=True)

    cg = pd.read_csv(_CANONICAL, low_memory=False)
    lab = (cg[(cg.species == pool.species)
              & cg.gf_tier.astype(str).str.contains("LAB", na=False)]
           .sort_values("wavelength_air_A").reset_index(drop=True))
    W = lab.wavelength_air_A.values.astype(float)
    src, sig = [], []
    for w in a.wavelength_air_A.astype(float):
        i = int(np.argmin(np.abs(W - w)))
        hit = abs(W[i] - w) <= _MATCH_TOL_A
        src.append(str(lab.lab_source_tag.values[i]) if hit else None)
        sig.append(float(lab.gf_sigma_dex.values[i]) if hit else np.nan)
    a["lab_source"] = src
    a["cited_sigma_dex"] = sig

    missing = a.lab_source.isna().sum()
    if missing:
        # Every anchor line is graded BY SELECTION, so an unmatched one means the anchor
        # and canonical_gf disagree about what is graded — a contradiction, not a gap.
        raise SystemExit(
            f"{missing} of {len(a)} anchor lines carry no LAB-tier row in canonical_gf. "
            f"The anchor is defined as the laboratory-graded lines, so this is the pool "
            f"and the line list disagreeing about which lines those are — refusing to "
            f"compute a laboratory zero point on it.")
    return a
