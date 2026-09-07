#!/usr/bin/env python3
"""RYA-1191 — blend/quality check on every graded Fe line, in every band.

    python3 scripts/rya1191_blend_audit.py

Ryan (2026-09-07): "Every graded line, ALL bands: telluric-corrected on a VERIFIED-clean
holding + blend/quality-checked. Drop ONLY genuine blends/artifacts." Telluric is the
other legs; this is the blend half.

🔴 THE SYNTHESIS LINE LIST CANNOT SEE THE CONTAMINANT THAT MATTERS, AND THAT IS THE
FINDING, NOT A LIMITATION OF THIS SCRIPT. The VIS/red-optical synthesis list is
`GESv6_atom_hfs_iso.420_920nm` -- 141,233 rows, 206 species, and **not one of them is
molecular**: no CH, CN, OH, C2, NH or MgH. The stellar catalogue `linelist_solar.csv`
carries all six (CH alone has 4,300 rows), and in the G band 4290-4315 A CH is the LARGEST
species present, 538 rows against Fe's 109.

So Fe II 4303.170 was synthesised against a list in which the dominant absorber at that
wavelength does not exist. The fitter had nowhere to put the CH absorption except Fe, and
A = 10.0 (+2.5 dex) is what that looks like. The line is an artifact of the LINE LIST, not
of the measurement -- which is why it is `curation, not tuning` to drop it, and why the
audit below is run against the CATALOGUE and not against the list the synthesis used.

⚠️ Both are reported per line. The catalogue share is the physical question ("what else
absorbs here?"); the synthesis share is the operational one ("what could the fit have
known about?"), and a line where the two disagree is a line whose synthesis was blind.

🔴 THE CONTROL IS Fe II 4303.170, AND IT COMES FIRST. RYA-515 identified it independently
as an unambiguous artifact -- A = 10.0, +2.5 dex, with 36 of the 65 transitions in its
window belonging to CH (the G band). A blend statistic that cannot rediscover that line
is not measuring blending, so the audit reports where 4303.170 ranks BEFORE any other
line is judged. Rediscovering a known artifact is what licenses the rest of the list
(RYA-161: curation on the data, never a tuned cut).

WHAT IS MEASURED, AND WHY IT IS A SHARE RATHER THAN A COUNT
------------------------------------------------------------
Inside each line's own fit window (the band's `half_width_A`, the same width the
abundance is fitted over), from the SYNTHESIS line list the product actually used:

  * `n_transitions`      -- how many CATALOGUE rows fall in the window at all
  * `depth_share_target` -- the TARGET species' share of the summed catalogue depth
  * `molecular_share`    -- the share belonging to CH/CN/OH/C2/NH/MgH, none of which the
                            synthesis list contains
  * `top_contaminant`    -- the largest non-target species and its share
  * `synth_*`            -- the same two shares computed on the ATOM-ONLY synthesis list,
                            so "what absorbs here" and "what the fit could see" stay apart

⚠️ A COUNT OF NEIGHBOURS IS NOT A BLEND. A window can hold forty transitions that
together absorb nothing, and two that swamp the line. Only the DEPTH share says whether
the flux being fitted is the target's. Counts are reported beside it because that is what
the RYA-515 finding was phrased in, not because they decide anything.

⚠️ AND A LOW SHARE IS A FLAG, NOT A VERDICT. This says what else is in the window; it does
not say the measured abundance is wrong. Nothing is dropped by this script -- it ranks,
and the dropping is a separate, named decision.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

#: The known artifact this audit must rediscover before it is trusted for anything else.
CONTROL_LINE = {"species": "Fe II", "wavelength_air_A": 4303.170, "band": "VIS",
                "ticket": "RYA-515",
                "claim": "A = 10.0 (+2.5 dex); 36 of 65 transitions in window are CH"}

#: Below this share of the window's summed theoretical depth, the flux being fitted is
#: mostly not the target's. NOT tuned to the control: it is the point at which the target
#: stops being the majority absorber, and the control is checked to land well inside it.
SHARE_DOMINATED = 0.50
SHARE_BLENDED = 0.80


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", default="data/results/rya1191")
    a = ap.parse_args(argv)
    out = ROOT / a.out_dir; out.mkdir(parents=True, exist_ok=True)

    from pipeline.nearuv_synth import build_solar_context
    from config.synth_bands import SYNTH_BANDS
    from rya759_nearuv_fe_product import species_token
    import derive_band_products as D

    per_line = pd.read_csv(out / "rya1191_graded_line_telluric.csv")
    fams = (per_line[["ion", "band", "holding", "selector", "synth_band"]]
            .drop_duplicates())

    #: The molecular species the synthesis list does not carry. Named explicitly rather
    #: than detected by heuristic — "is this string a molecule" is not a question a
    #: substring test should answer (RYA-1191: `"O2" in "CO2"`).
    MOLECULES = frozenset({"CH", "CN", "OH", "C2", "NH", "MgH"})

    cat = pd.read_csv(ROOT / "data/linelists/linelist_solar.csv", low_memory=False)
    cat = cat[["element", "ion", "wavelength_air_A", "central_depth"]].dropna(
        subset=["wavelength_air_A"])
    cat = cat.sort_values("wavelength_air_A")
    cw = cat.wavelength_air_A.to_numpy(float)
    cel = cat.element.astype(str).str.strip().to_numpy()
    cion = cat.ion.astype(str).str.strip().to_numpy()
    cd = cat.central_depth.fillna(0.0).to_numpy(float)

    ctx: dict = {}
    def lists_for(bkey):
        if bkey not in ctx:
            b = SYNTH_BANDS[bkey]
            c = build_solar_context("Fe", 500000.0, linelist_file=str(b.linelist),
                                    apply_canonical_gf=True, star="solar")
            ll = c["linelist"]
            names = ll.dtype.names
            w = np.asarray(ll["wave_A"] if "wave_A" in names
                           else ll["wave_nm"] * 10.0, float)
            el = np.asarray([str(x).strip() for x in ll["element"]])
            td = (np.asarray(ll["theoretical_depth"], float)
                  if "theoretical_depth" in names else np.zeros_like(w))
            o = np.argsort(w)
            ctx[bkey] = (w[o], el[o], td[o], float(b.half_width_A))
        return ctx[bkey]

    rows = []
    seen = set()
    for _, f in fams.iterrows():
        bkey = str(f.synth_band)
        if bkey not in SYNTH_BANDS:
            continue
        w, el, td, hw = lists_for(bkey)
        target = species_token("Fe", str(f.ion))
        sub = per_line[(per_line.ion == f.ion) & (per_line.band == f.band)
                       & (per_line.holding == f.holding)
                       & (per_line.selector == f.selector)]
        for lam in sub.wavelength_air_A.unique():
            key = (str(f.ion), bkey, round(float(lam), 4))
            if key in seen:
                continue
            seen.add(key)
            rows.append(_one(float(lam), target, w, el, td, hw, str(f.ion), bkey,
                             f.band, cw, cel, cion, cd, MOLECULES))

    d = pd.DataFrame(rows).sort_values(["ion", "wavelength_air_A"])
    d.to_csv(out / "rya1191_blend_audit.csv", index=False)

    ctl = _control_rank(d)
    doc = {
        "ticket": "RYA-1191 — blend/quality check on every graded Fe line",
        "kind": "VERIFICATION — nothing is dropped by this script",
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "control": ctl,
        "thresholds": {"share_dominated": SHARE_DOMINATED, "share_blended": SHARE_BLENDED,
                       "note": ("a share, not a count — a window can hold forty "
                                "transitions that absorb nothing")},
        "n_lines": int(len(d)),
        "by_verdict": d.verdict.value_counts().to_dict(),
        "synthesis_blind_to_molecular": {
            "n": int(d.synthesis_blind_to_molecular.sum()),
            "meaning": ("🔴 graded lines where molecular species carry >=20% of the "
                        "catalogue absorption in the fit window, on products synthesised "
                        "against an ATOM-ONLY line list. The fit had nowhere to put that "
                        "absorption except the target species."),
            "lines": (d[d.synthesis_blind_to_molecular]
                      .sort_values("molecular_share", ascending=False)
                      [["ion", "band", "wavelength_air_A", "molecular_share",
                        "top_contaminant", "depth_share_target"]]
                      .head(40).to_dict("records")),
        },
        "dominated": (d[d.verdict == "TARGET-NOT-DOMINANT"]
                      .sort_values("depth_share_target")
                      .head(40).to_dict("records")),
    }
    (out / "rya1191_blend_audit.json").write_text(json.dumps(doc, indent=2) + "\n")

    print("=== CONTROL: Fe II 4303.170 (RYA-515: A=10.0, 36/65 transitions CH) ===")
    for k, v in ctl.items():
        print(f"  {k}: {v}")
    print(f"\n=== {len(d)} graded lines ===")
    print(d.verdict.value_counts().to_string())
    print(f"\n=== synthesis blind to molecular absorption: "
          f"{int(d.synthesis_blind_to_molecular.sum())} lines ===")
    bl = d[d.synthesis_blind_to_molecular].sort_values("molecular_share", ascending=False)
    print(bl[["ion", "band", "wavelength_air_A", "molecular_share", "top_contaminant",
              "depth_share_target"]].head(12).to_string(index=False))
    print("\n=== worst 15 by target depth share ===")
    cols = ["ion", "band", "wavelength_air_A", "n_transitions", "depth_share_target",
            "molecular_share", "top_contaminant", "top_contaminant_share", "verdict"]
    print(d.sort_values("depth_share_target")[cols].head(15).to_string(index=False))
    print(f"\nwrote {a.out_dir}/")
    return 0


def _one(lam, target, w, el, td, hw, ion, bkey, band,
         cw, cel, cion, cd, MOLECULES) -> dict:
    rec = {"ion": ion, "band": band, "synth_band": bkey,
           "wavelength_air_A": round(lam, 4), "half_width_A": hw}

    # ── the CATALOGUE, which carries molecules ───────────────────────────────────
    lo, hi = np.searchsorted(cw, lam - hw), np.searchsorted(cw, lam + hw)
    sp = np.array([f"{e} {i}".strip() for e, i in zip(cel[lo:hi], cion[lo:hi])])
    dep = cd[lo:hi]
    rec["n_transitions"] = int(hi - lo)
    # ⚠️ THE TWO LISTS SPELL THE ION DIFFERENTLY AND THE MISMATCH IS SILENT. The synthesis
    # list uses `species_token` form ("Fe 1"/"Fe 2"); the catalogue's `ion` column is ROMAN
    # ("I"/"II"). Comparing "1" against "I" matched nothing, so every target share came
    # back 0.0 and all 353 lines read TARGET-NOT-DOMINANT — a whole-pool verdict that
    # looked like a finding. The catalogue is matched on its own spelling, which is the
    # `ion` value this loop already carries.
    is_tgt = (cel[lo:hi] == "Fe") & (cion[lo:hi] == ion)
    is_mol = np.isin(cel[lo:hi], list(MOLECULES))
    tot = float(dep.sum())
    rec["summed_depth_window"] = round(tot, 4)
    if hi <= lo or tot <= 0:
        rec.update(depth_share_target=None, molecular_share=None, top_contaminant=None,
                   top_contaminant_share=None, verdict="NO-DEPTH-IN-CATALOGUE")
    else:
        share = float(dep[is_tgt].sum()) / tot
        rec["depth_share_target"] = round(share, 4)
        rec["molecular_share"] = round(float(dep[is_mol].sum()) / tot, 4)
        oth = pd.Series(dep[~is_tgt], index=sp[~is_tgt]).groupby(level=0).sum(
            ).sort_values(ascending=False)
        rec["top_contaminant"] = (str(oth.index[0]) if len(oth) else None)
        rec["top_contaminant_share"] = (round(float(oth.iloc[0]) / tot, 4)
                                        if len(oth) else 0.0)
        rec["verdict"] = ("TARGET-NOT-DOMINANT" if share < SHARE_DOMINATED
                          else "BLENDED" if share < SHARE_BLENDED else "CLEAN")

    # ── the SYNTHESIS list, which does not ───────────────────────────────────────
    m = np.abs(w - lam) <= hw
    rec["synth_n_transitions"] = int(m.sum())
    if m.any() and float(td[m].sum()) > 0:
        st = float(td[m].sum())
        rec["synth_depth_share_target"] = round(
            float(td[m][el[m] == target].sum()) / st, 4)
    else:
        rec["synth_depth_share_target"] = None
    # 🔴 the line the synthesis was BLIND to: molecular absorption it cannot represent
    rec["synthesis_blind_to_molecular"] = bool(
        rec.get("molecular_share") is not None and rec["molecular_share"] >= 0.20)
    return rec


def _control_rank(d: pd.DataFrame) -> dict:
    """Where the known artifact lands. If it is not near the bottom, the statistic is
    not measuring what RYA-515 measured and nothing below it should be acted on."""
    m = ((d.ion == "II") & (np.abs(d.wavelength_air_A - CONTROL_LINE["wavelength_air_A"])
                            < 0.01))
    if not m.any():
        return {**CONTROL_LINE, "found": False,
                "verdict": ("🔴 THE CONTROL LINE IS NOT IN THE GRADED POOL — this audit "
                            "cannot be validated against it")}
    r = d[m].iloc[0]
    pool = d[(d.ion == "II") & d.depth_share_target.notna()]
    rank = int((pool.depth_share_target <= r.depth_share_target).sum())
    return {**CONTROL_LINE, "found": True,
            "n_transitions": int(r.n_transitions),
            "depth_share_target": (None if pd.isna(r.depth_share_target)
                                   else float(r.depth_share_target)),
            "top_contaminant": r.top_contaminant,
            "molecular_share": (None if pd.isna(r.molecular_share)
                                else float(r.molecular_share)),
            "synthesis_blind_to_molecular": bool(r.synthesis_blind_to_molecular),
            "top_contaminant_share": (None if pd.isna(r.top_contaminant_share)
                                      else float(r.top_contaminant_share)),
            "verdict_assigned": r.verdict,
            "rank_among_Fe_II_by_share": f"{rank} of {len(pool)} (1 = worst)",
            # 🔴 THE CONTROL PASSES ONLY IF IT REDISCOVERS *WHY*. A verdict of
            # TARGET-NOT-DOMINANT is not enough — 145 of 353 lines carry it. RYA-515's
            # finding is specifically that CH swamps this window, so the control demands
            # the molecular share and the CH identification, not just the label.
            "control_passes": bool(r.verdict == "TARGET-NOT-DOMINANT"
                                   and bool(r.synthesis_blind_to_molecular)
                                   and str(r.top_contaminant).startswith("CH"))}


if __name__ == "__main__":
    raise SystemExit(main())
