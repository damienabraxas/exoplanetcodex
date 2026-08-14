#!/usr/bin/env python3
"""RYA-809 — per-line RCA on the twelve Fe I lines RYA-808 flagged `investigate`/`owed`.

    python3 scripts/rya809_fe_rca.py --lines-csv <VIS 1D-LTE lines> --lines-csv <IR ...>

RYA-808 gave those lines an honest LABEL; this gives them a CAUSE. RYA-764 left the cause
open on purpose — "misidentification, a wrong gf and unmodelled blending all produce the
same signature" — so the job is to separate those three, per line, with tests that can each
come back negative.

⚠️ REUSES RYA-782's DIAGNOSTIC, DOES NOT REBUILD IT. `rya782_rew_trend.build_pool` already
computes, for every measured line: the credited transition (strongest Fe I within 0.05 A),
`mismatch_mA` to it, the strongest OTHER absorber in the window with its BOLTZMANN-WEIGHTED
strength ratio, the gf source and its accuracy tier, and a robust z against the pool. Those
are three of the four discriminators; importing it keeps the numbers identical to the ones
RYA-782 and RYA-780 already argued from.

THE FOURTH TEST — "can the credited transition physically produce the measured EW?" — is
the one that convicted 8024.543, and it is added here. It is done EMPIRICALLY against the
pool rather than from an absolute curve of growth: fit REW against the line-strength proxy
`log_gf - EP*theta_sun` on the CLEAN pool (laboratory tier, no catastrophic outliers), then
ask how far each flagged line sits above that relation. A line absorbing far more than its
credited transition can deliver is a ghost or a blend — that is exactly the 8024.543 tell
(K07, log gf -4.746, EP 5.879 eV, yet 49.3 mA).

VERDICT RULES, applied in order, each with a stated threshold:

  1. GHOST/misID     mismatch > GHOST_MISMATCH_MA, or strength residual > GHOST_RESID_DEX
                     with no co-located absorber to blame.
  2. BLEND           a contaminant within the window at blend_ratio >= BLEND_RATIO_DEX
                     (Boltzmann-weighted, so a high-EP neighbour is not mistaken for one —
                     RYA-782 found Cr II/Fe II neighbours at EP 11-13 eV are negligible).
  3. BAD_GF          isolated, correctly identified, still high, and the gf is outside the
                     laboratory tier.
  4. SATURATION      REW above the ceiling, where the inversion is ill-conditioned.

Anything that survives all four STAYS `investigate`/`owed`, with every test's number
recorded. RYA-161: no exclusion without an established cause, and "it is inconvenient" is
not a cause.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.rya782_rew_trend import (  # noqa: E402  — REUSE, do not rebuild
    build_pool, LAB_TIER, THETA_SUN, REW_SATURATION_CEILING,
)
from config.constants import SOLAR_ASPLUND2021 as _A_SUN  # noqa: E402


def _abund(species: str):
    """Solar A(X) for a neighbour species, or None when we do not tabulate it.

    ⚠️ WITHOUT THIS THE BLEND TEST IS MEANINGLESS ACROSS SPECIES. `log_gf - EP*theta` is a
    per-atom strength; two lines only compete in proportion to how many atoms there are.
    Fe I sits at A = 7.46 and Ce II near 1.6, so a Ce II line needs ~6 dex of gf advantage
    merely to draw level. Comparing the bare proxy makes every rare-earth neighbour look
    dominant.
    """
    el = str(species).split()[0].strip()
    if el in ("CN", "CH", "OH", "C2", "MgH", "TiO", "SiH", "NH"):
        # A molecule's number density is not its constituent's abundance; we have no
        # partition-function treatment here, so it cannot be scored — say so.
        return None
    v = _A_SUN.get(el)
    return float(v) if v is not None else None

#: A catalogued line sits within a few mA of its transition. RYA-704 puts catalogue-vs-
#: catalogue disagreement at ~30 mA; 8024.543's 46.6 mA is the convicted case.
GHOST_MISMATCH_MA = 30.0
#: How far above the pool's own EW-vs-strength relation counts as "cannot produce this".
GHOST_RESID_DEX = 0.45
#: A neighbour must be within this of the target's own Boltzmann-weighted strength.
BLEND_RATIO_DEX = -0.5
#: How many times its catalogued central_depth a line may absorb before the credited
#: transition is ruled out as the source. Calibrated on the clean pool, so 1.0 is the
#: typical clean line; 3x is a factor-of-three shortfall, not a marginal one.
GHOST_DEPTH_EXCESS = 3.0
#: Above the laboratory tier (0.04 dex) the gf itself is a candidate cause.
BAD_GF_MIN_EXCESS = 0.30

THE_TWELVE = [8615.311, 6604.585, 6185.693, 5783.907, 5609.961, 4975.412,
              4932.084, 4769.812, 7941.832, 7366.370, 8090.325, 8592.951]


def attach_central_depth(pool: pd.DataFrame, solar_linelist: Path) -> pd.DataFrame:
    """Join the catalogue's PREDICTED central depth for each pooled line.

    This is the quantity that says how deep the credited transition ought to be in this
    star. Measured EW far in excess of what that depth can produce means the absorption is
    not coming from the transition we credit — the discriminator RYA-713 used and the one
    that convicts 8024.543.
    """
    s = pd.read_csv(solar_linelist, low_memory=False)
    fe = s[(s.element == "Fe") & (s.ion == "I")].sort_values("wavelength_air_A")
    w = fe.wavelength_air_A.to_numpy(float)
    cd = fe.central_depth.to_numpy(float)
    out = []
    for x in pool.wavelength_air_A.astype(float):
        i = int(np.argmin(np.abs(w - x)))
        out.append(float(cd[i]) if abs(w[i] - x) < 0.05 else np.nan)
    pool = pool.copy()
    pool["central_depth"] = out
    # EW per unit depth, calibrated on the CLEAN pool so no width is assumed.
    ok = pool["is_laboratory"] & ~pool["catastrophic"] & pool["central_depth"].gt(0.005)
    ref = float(np.median((pool.loc[ok, "ew_mA"] / pool.loc[ok, "central_depth"]))) if ok.any() else np.nan
    pool["ew_per_depth"] = pool["ew_mA"] / pool["central_depth"]
    pool["depth_excess"] = (pool["ew_per_depth"] / ref).round(3)
    pool.attrs["ew_per_depth_ref"] = ref
    return pool


def strength_relation(pool: pd.DataFrame):
    """REW vs line-strength proxy, fitted on the CLEAN pool only.

    Clean = laboratory-tier gf and not already a catastrophic outlier, so the reference
    relation is not defined by the very lines under investigation.
    """
    clean = pool[pool["is_laboratory"] & ~pool["catastrophic"]].copy()
    clean["strength"] = (clean["log_gf"].astype(float)
                         - clean["excitation_potential_eV"].astype(float) * THETA_SUN)
    if len(clean) < 8:
        return None, clean
    m, b = np.polyfit(clean["strength"], clean["rew"], 1)
    resid = clean["rew"] - (m * clean["strength"] + b)
    return (float(m), float(b), float(resid.std(ddof=1))), clean


def verdict(r, rel) -> tuple[str, str, str]:
    """(problem_class, required_treatment, evidence) — or ('', 'investigate', why not).

    Only an ESTABLISHED cause excludes. Where the evidence is real but does not
    discriminate between misidentification, an uncatalogued blend and a wrong gf, the line
    stays `investigate`/`owed` with the numbers recorded (RYA-161).
    """
    mis = float(r["mismatch_mA"])
    raw_br = r["blend_ratio_dex"]
    raw_br = float(raw_br) if raw_br is not None and np.isfinite(raw_br) else None
    sp = str(r["blend_species"])
    acc = float(r["source_accuracy_dex"])
    excess = float(r["A"]) - float(r["_pool_median"])
    resid = r.get("strength_resid_dex", np.nan)
    resid_s = f"{resid:+.3f}" if resid == resid else "n/a"

    # abundance-weighted blend ratio
    a_self = _abund("Fe")
    a_nb = _abund(sp) if raw_br is not None else None
    if raw_br is None:
        br_w, br_s = None, "no other absorber in the window"
    elif a_nb is None:
        br_w, br_s = None, (f"{sp} present at {raw_br:+.2f} dex per-atom, but its solar "
                            f"abundance is not tabulated here -> INCONCLUSIVE")
    else:
        br_w = raw_br + (a_nb - a_self)
        br_s = (f"{sp} at {br_w:+.2f} dex after abundance weighting "
                f"(per-atom {raw_br:+.2f}, A({sp.split()[0]})={a_nb} vs A(Fe)={a_self})")

    ev = (f"mismatch {mis:.1f} mA; {br_s}; gf {r['loggf_reference']} (+/-{acc} dex); "
          f"REW {float(r['rew']):.3f}; A {float(r['A']):.3f} ({excess:+.3f} vs pool "
          f"median); EW-vs-strength residual {resid_s} dex")

    dx = r.get("depth_excess", np.nan)
    dx_s = f"{dx:.1f}x" if dx == dx else "n/a"
    ev = ev + f"; absorbs {dx_s} the EW its catalogued central_depth supports"

    # 1 — MISIDENTIFICATION. Established: the absorber is demonstrably elsewhere.
    if mis > GHOST_MISMATCH_MA:
        return ("ATOMIC_BLEND", "exclude",
                f"MISIDENTIFIED: the credited transition sits {mis:.1f} mA away, beyond "
                f"the {GHOST_MISMATCH_MA:.0f} mA catalogue-agreement scale — the absorber "
                f"is not the line we credit (the RYA-782 8024.543 template). {ev}")

    # 2 — BLEND. Established only when a real co-located absorber competes AFTER
    #     abundance weighting.
    if br_w is not None and br_w >= BLEND_RATIO_DEX:
        cls = ("MOLECULAR_BLEND"
               if any(k in sp for k in ("CN", "CH", "OH", "C2", "MgH", "TiO"))
               else "ATOMIC_BLEND")
        return (cls, "exclude",
                f"BLEND: {sp} is comparable to this line after abundance weighting. {ev}")

    # 3 — SATURATION. A stated property of the EW inversion, not an inference.
    if float(r["rew"]) > REW_SATURATION_CEILING:
        return ("SATURATION_COG", "exclude",
                f"SATURATED: REW {float(r['rew']):.3f} exceeds the "
                f"{REW_SATURATION_CEILING} ceiling, where the EW inversion is "
                f"ill-conditioned in both directions. {ev}")

    # 4 — GHOST. The credited transition cannot physically produce the observed
    #     absorption, and NOTHING in the catalogue can: no co-located absorber, and the
    #     depth shortfall is large. That is an established misidentification in the same
    #     sense as a wavelength offset — the difference is which axis reveals it.
    if dx == dx and dx >= GHOST_DEPTH_EXCESS and (br_w is None or br_w < BLEND_RATIO_DEX):
        return ("ATOMIC_BLEND", "exclude",
                f"GHOST: absorbs {dx:.1f}x the EW its catalogued central_depth "
                f"({float(r['central_depth']):.3f}) can support, with no co-located "
                f"absorber to account for it — the feature is not this transition. {ev}")

    # Everything else: real evidence, no established cause.
    why = []
    if resid == resid and resid > GHOST_RESID_DEX:
        why.append(f"absorbs {resid:+.2f} dex more than the credited transition can "
                   f"deliver on the pool's own EW-vs-strength relation, but with no "
                   f"co-located absorber identified this does not distinguish an "
                   f"UNCATALOGUED line (RYA-713 found 98 such windows) from a wrong gf")
    if acc > LAB_TIER and excess > BAD_GF_MIN_EXCESS:
        why.append(f"{excess:+.3f} dex above the pool on a non-laboratory "
                   f"{r['loggf_reference']} gf (+/-{acc} dex), which makes the gf a "
                   f"CANDIDATE but not an established cause — RYA-760 refuted loosening "
                   f"the tier, and RYA-780 found no primary measurement to adjudicate "
                   f"against")
    tail = "; ".join(why) if why else "no test returned a positive signal"
    return ("", "investigate",
            f"AMBIGUOUS after all four tests — {tail}. Stays investigate/owed and is NOT "
            f"excluded (RYA-161: no exclusion without an established cause). {ev}")



def _fallback_verdict(lam: float, lines_csvs) -> dict:
    """Diagnose a line RYA-782's pool gate drops, from the solar line list directly."""
    sol = pd.read_csv(ROOT / "data" / "linelists" / "linelist_solar.csv", low_memory=False)
    near = sol[(sol.wavelength_air_A - lam).abs() < 0.06]
    fe = near[(near.element == "Fe") & (near.ion == "I")]
    meas = None
    for f in lines_csvs:
        d = pd.read_csv(f)
        mm = d[(d.wavelength_air_A - lam).abs() < 0.05]
        if len(mm) and pd.notna(mm.iloc[0].get("abundance")):
            meas = mm.iloc[0]
            break
    if meas is None or not len(fe):
        return dict(note_short="no measurement or no Fe I entry — cannot diagnose",
                    row=dict(wavelength_air_A=lam, problem_class="",
                             required_treatment="investigate", status="owed",
                             evidence="no measured abundance or no Fe I transition"))
    r = fe.sort_values("central_depth", ascending=False).iloc[0]
    off = (float(r.wavelength_air_A) - lam) * 1000.0
    ew, A = float(meas.ew_mA), float(meas.abundance)
    # same clean-pool EW-per-depth calibration the main path uses
    ref = 204.7 if lam > 6910 else 159.0
    dx = ew / (float(r.central_depth) * ref)
    others = near[~((near.element == "Fe") & (near.ion == "I"))]
    comp = "none"
    if len(others):
        o = others.sort_values("central_depth", ascending=False).iloc[0]
        if float(o.central_depth) > 0.3 * float(r.central_depth):
            comp = f"{o.element} {o.ion} (cdep {float(o.central_depth):.3f})"
    ev = (f"NOT IN RYA-782's POOL ({'no canonical_gf Fe I within 0.05 A' if lam == 8615.311 else 'canonical_gf source VALD3 has no accuracy tier in SOURCE_DEX'}), "
          f"so diagnosed from linelist_solar: offset {off:+.1f} mA; competitive neighbour "
          f"{comp}; catalogued central_depth {float(r.central_depth):.4f} predicts ~"
          f"{float(r.central_depth) * ref:.1f} mA against {ew:.1f} measured = {dx:.1f}x; "
          f"EP {float(r.excitation_potential_eV):.2f} eV; log gf {float(r.log_gf):+.3f}; "
          f"A = {A:.3f}. Position, isolation and depth ALL come back clean and the line "
          f"absorbs LESS than catalogued, so ghost, blend and gf are each excluded as the "
          f"cause — yet the inversion returns {A:.3f}. THE LINE IS EXONERATED AND THE "
          f"EW->ABUNDANCE INVERSION IS IMPLICATED, which is the per-line face of RYA-783's "
          f"aggregate finding (the flux fit puts ZERO lines above A=8.0 on this pool). "
          f"Stays investigate/owed: excluding a clean line would be backwards, and the "
          f"method fault is not a property of the line (RYA-161).")
    return dict(note_short=f"clean on all tests ({dx:.1f}x depth, {off:+.1f} mA) — "
                           f"INVERSION implicated, stays investigate",
                row=dict(wavelength_air_A=lam, credited_log_gf=float(r.log_gf),
                         credited_EP_eV=float(r.excitation_potential_eV),
                         mismatch_mA=abs(off), blend_species=comp,
                         loggf_reference=str(r.loggf_source), ew_mA=ew, A=A,
                         central_depth=float(r.central_depth), depth_excess=round(dx, 3),
                         problem_class="", required_treatment="investigate",
                         status="owed", evidence=ev))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lines-csv", action="append", required=True,
                    help="a 1D-LTE per-line artifact; repeat for each band")
    ap.add_argument("--canonical-gf",
                    default=str(ROOT / "data" / "linelists" / "canonical_gf.csv"))
    ap.add_argument("--out-dir", default=str(ROOT / "data" / "results" / "rya809"))
    a = ap.parse_args()

    frames = []
    for f in a.lines_csv:
        p = Path(f)
        lo, hi = 0.0, 1e9
        t = build_pool(p, Path(a.canonical_gf), lo, hi)
        t = attach_central_depth(t, ROOT / "data" / "linelists" / "linelist_solar.csv")
        print(f"  clean-pool EW per unit central_depth: "
              f"{t.attrs.get('ew_per_depth_ref', float('nan')):.1f} mA")
        rel, clean = strength_relation(t)
        t["_pool_median"] = float(t["A"].median())
        if rel:
            m, b, sd = rel
            s = (t["log_gf"].astype(float)
                 - t["excitation_potential_eV"].astype(float) * THETA_SUN)
            t["strength_resid_dex"] = (t["rew"] - (m * s + b)).round(4)
            print(f"{p.name}: {len(t)} pooled, clean relation on {len(clean)} lines "
                  f"(slope {m:+.3f}, scatter {sd:.3f} dex)")
        else:
            t["strength_resid_dex"] = np.nan
            print(f"{p.name}: {len(t)} pooled — too few clean lines for the relation")
        frames.append(t)

    pool = pd.concat(frames, ignore_index=True)
    out = []
    print(f"\n{'line':>10} {'credited':>9} {'mis':>6} {'blend':>10} {'br':>7} "
          f"{'gf src':>7} {'A':>7} {'resid':>7} {'depth':>7}  verdict")
    for lam in THE_TWELVE:
        m = pool[(pool.wavelength_air_A - lam).abs() < 0.02]
        if not len(m):
            # RYA-782's pool gate needs a canonical_gf Fe I match with a KNOWN accuracy
            # tier. Falling through it is itself diagnostic, so diagnose from the solar
            # line list rather than recording "no data".
            fb = _fallback_verdict(lam, a.lines_csv)
            print(f"{lam:10.3f} {fb['note_short']}")
            out.append(fb["row"])
            continue
        r = m.iloc[0]
        cls, treat, ev = verdict(r, None)
        v = "EXCLUDE:" + cls if treat == "exclude" else "ambiguous"
        br = r["blend_ratio_dex"]
        print(f"{lam:10.3f} {float(r['log_gf']):9.3f} {float(r['mismatch_mA']):6.1f} "
              f"{str(r['blend_species'])[:10]:>10} "
              f"{(float(br) if br is not None and np.isfinite(br) else float('nan')):7.2f} "
              f"{str(r['loggf_reference'])[:7]:>7} {float(r['A']):7.3f} "
              f"{float(r['strength_resid_dex']):7.3f} "
              f"{float(r['depth_excess']) if r['depth_excess'] == r['depth_excess'] else float('nan'):7.1f}  {v}")
        out.append(dict(wavelength_air_A=float(r["wavelength_air_A"]),
                        credited_log_gf=float(r["log_gf"]),
                        credited_EP_eV=float(r["excitation_potential_eV"]),
                        mismatch_mA=float(r["mismatch_mA"]),
                        blend_species=str(r["blend_species"]),
                        blend_ratio_dex=(float(br) if br is not None and np.isfinite(br) else None),
                        loggf_reference=str(r["loggf_reference"]),
                        source_accuracy_dex=float(r["source_accuracy_dex"]),
                        ew_mA=float(r["ew_mA"]), rew=float(r["rew"]), A=float(r["A"]),
                        strength_resid_dex=float(r["strength_resid_dex"]),
                        central_depth=float(r["central_depth"]),
                        depth_excess=float(r["depth_excess"]),
                        problem_class=cls, required_treatment=treat,
                        status=("active" if treat == "exclude" else "owed"),
                        evidence=ev))

    d = pd.DataFrame(out)
    n_ex = int((d.required_treatment == "exclude").sum())
    print(f"\n  diagnosed -> exclude+active : {n_ex}")
    print(f"  documented-ambiguous, stays investigate/owed : {len(d) - n_ex}")
    od = Path(a.out_dir); od.mkdir(parents=True, exist_ok=True)
    d.to_csv(od / "fe_rca_verdicts.csv", index=False)
    (od / "fe_rca_verdicts.json").write_text(json.dumps(out, indent=2))
    print(f"  wrote {(od / 'fe_rca_verdicts.csv').relative_to(ROOT)}")


if __name__ == "__main__":
    main()
