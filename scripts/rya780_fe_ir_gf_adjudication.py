#!/usr/bin/env python3
"""RYA-780 — adjudicate the disputed Fe I IR gf sources against PRIMARY measurements.

    python3 scripts/rya780_fe_ir_gf_adjudication.py --lines-csv <1D-LTE per-line artifact>

WHAT RYA-760 LEFT OWED
----------------------
RYA-760 established that the Fe I IR lines carrying FMW (and GESB82c) oscillator strengths
return A(Fe) high by +0.294 dex once EP and REW are controlled for, and that the offset is
in the ATOMIC DATA rather than in our analysis. It could not adjudicate the sources,
because it had no independent referee:

  * GES and VALD are not two opinions on FMW -- 96.6% of FMW lines are byte-identical
    between them. Agreement is inheritance.
  * NIST ASD is not independent either: FMW = Fuhr, Martin & Wiese IS a NIST compilation.
  * and ASD has been returning HTTP 500 (lines1.pl line 701) throughout.

Only a PRIMARY laboratory measurement can referee a compilation's scale. This script is
that referee, using two vendored primary sources (data/linelists/primary_gf/):
Den Hartog et al. 2014 (ApJS 215, 23) and Ruffoni et al. 2014 (MNRAS 441, 3127), both
FTS branching fractions normalised on laser-induced-fluorescence lifetimes.

THE THREE THINGS IT ESTABLISHES
-------------------------------
1. TRANSITION IDENTITY. The pool's canonical_gf carries one source per line and no level
   labels, so "the same transition" cannot be checked inside it. The VALD raw extract
   does carry per-line LS terms, J values and its OWN gf source code -- so the identity,
   and VALD's independent source attribution, are read from there.

2. COVERAGE. For each disputed transition: does a primary measurement of it exist? This
   decides RECOVERED-vs-QUARANTINED per RYA-354, and it is a fact about the literature,
   not a judgement call.

3. THE SCALE COMPARISON, for the lines no primary covers. Each disputed line is compared
   against a reference scale fitted to the pool's OWN primary-sourced lines at matched EP
   and REW. That reference is itself validated against Den Hartog (median +0.000 dex on
   the overlap), so the comparison is controlled rather than circular.

THE RYA-161 FIREWALL, EXPLICITLY
--------------------------------
Nothing here compares a line to the optical anchor 7.466, or to any target value. The
reference is the PRIMARY-SOURCE POPULATION. A quarantine reason that cited delta-vs-anchor
would be tuning; a quarantine reason that cites "reads +0.6 dex against laboratory-sourced
lines of the same strength and excitation" is provenance. Only the second is emitted.

MATCHING RULE: wavelength AND lower-level energy. Wavelength alone at 0.06 A pairs our
K07 8876.0059 (EP 4.584 eV) with Den Hartog's 8876.0241 (5.020 eV) -- a different
transition, and a manufactured 2.85 dex 'discrepancy'. The EP filter is load-bearing.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

CANON = ROOT / "data" / "linelists" / "canonical_gf.csv"
VALD = ROOT / "data" / "linelists" / "vald_solar_redopt_6910_9500_hfson_raw.txt"
PRIMARY_DIR = ROOT / "data" / "linelists" / "primary_gf"
OUT_DIR = ROOT / "data" / "results" / "fe_ir_gf"

# cm^-1 -> eV (CODATA via scipy, never a typed-in constant)
import scipy.constants as _sc                                    # noqa: E402
CM_PER_EV = _sc.e / (_sc.h * _sc.c) / 100.0

# Sources that are PRIMARY laboratory measurements. FMW and GESB82c are the disputed
# compilations; K07 is Kurucz semi-empirical (RYA-709's problem, explicitly out of scope).
PRIMARY_PREFIXES = ("BWL", "BK", "BKK", "2014MNRAS", "GESHRL")
DISPUTED_PATTERN = r"FMW|GESB82c"

# Matching tolerances. Wavelength is generous because catalogues disagree at the mA level;
# the EP filter is what actually establishes identity (see the module docstring).
WL_TOL_A = 0.06
EP_TOL_EV = 0.01


# ── inputs ───────────────────────────────────────────────────────────────────

def read_vald_terms(path: Path = VALD) -> pd.DataFrame:
    """Fe I transition identities from the VALD raw extract.

    VALD writes each line as four records: the data row, the lower LS term, the upper LS
    term, and a reference row carrying per-quantity source codes (`gf:FMW`). That gives
    two things canonical_gf cannot: the transition's LEVELS, and VALD's OWN attribution of
    where its gf came from -- which is how we test whether VALD is a second opinion or is
    simply carrying the same compiled number.
    """
    lines = path.read_text(encoding="latin1").splitlines()
    rows, i = [], 0
    while i < len(lines):
        L = lines[i]
        if L.startswith("'") and "'," in L:
            p = [x.strip() for x in L.split(",")]
            try:
                spec = p[0].strip("'").strip()
                wl, gf = float(p[1]), float(p[2])
                elo, jlo, eup, jup = float(p[3]), float(p[4]), float(p[5]), float(p[6])
            except (ValueError, IndexError):
                i += 1
                continue
            def _clean(k):
                s = lines[k].strip().strip("'").strip() if k < len(lines) else ""
                return re.sub(r"^LS\s+", "", s).strip()
            ref = lines[i + 3].strip().strip("'").strip() if i + 3 < len(lines) else ""
            m = re.search(r"gf:(\S+)", ref)
            rows.append(dict(spec=spec, wl=wl, vald_gf=gf, elo=elo, jlo=jlo,
                             eup=eup, jup=jup, lo_term=_clean(i + 1),
                             up_term=_clean(i + 2), vald_src=m.group(1) if m else None))
            i += 4
        else:
            i += 1
    d = pd.DataFrame(rows)
    return d[d.spec == "Fe 1"].reset_index(drop=True)


def read_primary(name: str) -> pd.DataFrame:
    """One vendored VizieR TSV -> (wave_A, ep_eV, loggf, e_loggf, prev, prev_ref)."""
    path = PRIMARY_DIR / name
    raw = [l for l in path.read_text().splitlines() if l and not l.startswith("#")]
    hdr = raw[0].split("\t")
    d = pd.DataFrame([r.split("\t") for r in raw[3:]], columns=hdr)   # rows 1-2 are units/rule
    d = d.rename(columns={"lamAir": "wave_A", "lam.Air": "wave_A",
                          "log(gf)": "loggf", "e_log(gf)": "e_loggf",
                          "log(gf)P": "prev", "r_log(gf)P": "prev_ref",
                          "loggf0": "prev", "e_loggf": "e_loggf"})
    for c in ("wave_A", "loggf", "e_loggf", "prev", "E0", "chi"):
        if c in d.columns:
            d[c] = pd.to_numeric(d[c].astype(str).str.strip(), errors="coerce")
    # Den Hartog gives the lower level in cm^-1 (E0); Ruffoni gives it in eV (chi).
    if "E0" in d.columns:
        d["ep_eV"] = d.E0 / CM_PER_EV
    elif "chi" in d.columns:
        d["ep_eV"] = d.chi
    else:
        d["ep_eV"] = np.nan
    if "prev_ref" in d.columns:
        d["prev_ref"] = d.prev_ref.astype(str).str.strip()
    else:
        d["prev_ref"] = ""
    return d[d.loggf.notna()][["wave_A", "ep_eV", "loggf", "e_loggf", "prev", "prev_ref"]]


def match_primary(prim: pd.DataFrame, wave: float, ep: float):
    """The one primary row for this transition, or None.

    Requires BOTH wavelength and lower-level energy. Wavelength alone is not identity —
    it pairs different transitions that happen to lie together (see the docstring).
    """
    k = prim[(prim.wave_A - wave).abs() <= WL_TOL_A]
    if not len(k):
        return None
    k = k[(k.ep_eV - ep).abs() <= EP_TOL_EV]
    if not len(k):
        return None
    return k.iloc[(k.wave_A - wave).abs().argsort()].iloc[0]


# ── the analysis ─────────────────────────────────────────────────────────────

def build(lines_csv: Path) -> pd.DataFrame:
    """Join the measured Fe IR lines to their catalogue source and transition identity."""
    g = pd.read_csv(CANON, low_memory=False)
    cat = g[(g.species == "Fe I") & g.wavelength_air_A.between(6900, 9210)]
    vald = read_vald_terms()
    meas = pd.read_csv(lines_csv)
    meas = meas[meas.abundance.notna() & meas.in_aggregate]

    rows = []
    for _, r in meas.iterrows():
        c = float(r.wavelength_air_A)
        k = cat[(cat.wavelength_air_A - c).abs() <= 0.05]
        v = vald[(vald.wl - c).abs() <= 0.03]
        if not len(k) or not len(v):
            continue
        k = k.iloc[(k.wavelength_air_A - c).abs().argsort()].iloc[0]
        v = v.iloc[(v.wl - c).abs().argsort()].iloc[0]
        rows.append(dict(wave_A=c, A=float(r.abundance), rew=float(r.rew),
                         ew_mA=float(r.ew_mA), ep_eV=float(k.excitation_potential_eV),
                         ges_src=str(k.loggf_reference), ges_gf=float(k.log_gf),
                         vald_src=str(v.vald_src), vald_gf=float(v.vald_gf),
                         lower_term=v.lo_term, upper_term=v.up_term,
                         J_lo=v.jlo, J_up=v.jup))
    t = pd.DataFrame(rows)
    t["is_primary"] = t.ges_src.str.startswith(PRIMARY_PREFIXES)
    t["is_disputed"] = t.ges_src.str.contains(DISPUTED_PATTERN, na=False)
    t["multiplet"] = t.lower_term + " -> " + t.upper_term
    return t


def primary_scale(t: pd.DataFrame):
    """Fit A(Fe) ~ EP + REW on the pool's PRIMARY-sourced lines only.

    This is the reference the disputed lines are judged against. Fitting it on the primary
    lines alone matters: a fit over everything would absorb the very offset under test, and
    the K07 majority would dominate it.
    """
    p = t[t.is_primary]
    X = np.column_stack([np.ones(len(p)), p.ep_eV, p.rew])
    beta, *_ = np.linalg.lstsq(X, p.A.values, rcond=None)
    rms = float(np.std(p.A.values - X @ beta, ddof=3))
    return beta, rms, len(p)


def validate_scale(t: pd.DataFrame, prims: dict) -> pd.DataFrame:
    """Do the vendored primary measurements corroborate the pool's primary gf values?

    If they do not, the reference scale is not a reference and nothing below is worth
    reading -- so this runs before the adjudication, not as an afterthought.
    """
    out = []
    for _, x in t[t.is_primary].iterrows():
        for name, prim in prims.items():
            k = match_primary(prim, x.wave_A, x.ep_eV)
            if k is not None:
                out.append(dict(wave_A=x.wave_A, ges_src=x.ges_src, ges_gf=x.ges_gf,
                                primary=name, primary_gf=float(k.loggf),
                                delta=float(k.loggf) - x.ges_gf))
    return pd.DataFrame(out)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lines-csv", required=True,
                    help="the Fe IR 1D-LTE per-line artifact (wavelength_air_A, abundance, "
                         "rew, ew_mA, in_aggregate). Passed in rather than hardcoded: it "
                         "is a band product and does not live on main yet (RYA-759/783).")
    ap.add_argument("--out", default=str(OUT_DIR))
    ap.add_argument("--sigma", type=float, default=2.0,
                    help="residual significance, in units of the primary-scale rms, at "
                         "or above which a disputed line is QUARANTINED")
    a = ap.parse_args()

    t = build(Path(a.lines_csv))
    prims = {"DenHartog2014": read_primary("denhartog2014_FeI_6900_9250.tsv"),
             "Ruffoni2014": read_primary("ruffoni2014_FeI_table5.tsv")}
    print(f"measured Fe I IR lines with a transition identity: {len(t)}  "
          f"(primary-sourced {int(t.is_primary.sum())}, disputed {int(t.is_disputed.sum())})")

    # ── 1. can any catalogue referee the disputed set? ───────────────────────
    dis = t[t.is_disputed].sort_values("wave_A").copy()
    same = int((dis.vald_src.str.contains("FMW", na=False)).sum())
    ident = int((np.abs(dis.vald_gf - dis.ges_gf) < 0.005).sum())
    print(f"\n1. CATALOGUE INDEPENDENCE — VALD names FMW for {same}/{len(dis)} disputed "
          f"lines; {ident}/{len(dis)} gf values agree to <0.005 dex.")
    ges82 = dis[dis.ges_src.str.contains("GESB82c")]
    if len(ges82) and ges82.vald_src.str.contains("FMW").all():
        print(f"   NOTE: the {len(ges82)} line(s) GES labels GESB82c "
              f"({', '.join(f'{w:.4f}' for w in ges82.wave_A)}) are FMW in VALD "
              f"=> the disputed population is FMW n={len(dis)}, not 12 + a separate "
              f"GESB82c pair.")

    # ── 2. is the reference scale itself corroborated? ───────────────────────
    val = validate_scale(t, prims)
    print(f"\n2. REFERENCE-SCALE VALIDATION — primary measurements vs the pool's "
          f"primary-sourced gf ({len(val)} overlaps):")
    if len(val):
        for nm, gp in val.groupby("primary"):
            print(f"   {nm}: n={len(gp)} median delta {gp.delta.median():+.3f} dex, "
                  f"max |delta| {gp.delta.abs().max():.3f}")
    beta, rms, npx = primary_scale(t)
    print(f"   primary scale: A = {beta[0]:.3f} {beta[1]:+.4f}*EP {beta[2]:+.4f}*REW "
          f"(n={npx}, rms {rms:.3f})")

    # ── 3. per-line adjudication ─────────────────────────────────────────────
    Xd = np.column_stack([np.ones(len(dis)), dis.ep_eV, dis.rew])
    dis["A_primary_scale"] = Xd @ beta
    dis["residual_dex"] = dis.A - dis.A_primary_scale
    dis["n_sigma"] = dis.residual_dex / rms

    recs = []
    for _, x in dis.iterrows():
        hit = {nm: match_primary(p, x.wave_A, x.ep_eV) for nm, p in prims.items()}
        hit = {k: v for k, v in hit.items() if v is not None}
        sib = t[(t.multiplet == x.multiplet) & t.is_primary]
        # The ticket's RECOVERED/QUARANTINED binary presupposes that a primary measurement
        # of the transition EXISTS. Where none does, neither verdict is available and
        # saying so is the finding -- so coverage is the first axis, and the scale
        # comparison only ranks the lines the literature cannot adjudicate.
        if hit:
            nm, k = next(iter(hit.items()))
            d = float(k.loggf) - x.ges_gf
            disp = "RECOVERED" if abs(d) >= 0.10 else "CONFIRMED-BY-PRIMARY"
            reason = (f"{nm} measures log gf = {k.loggf:+.3f} for this transition vs the "
                      f"pool's {x.ges_gf:+.3f} ({d:+.3f} dex).")
            pg, pe, pn = float(k.loggf), float(k.e_loggf), nm
        else:
            disp = ("QUARANTINED-SCALE-EVIDENCE" if x.n_sigma >= a.sigma
                    else "NO-PRIMARY-NO-EVIDENCE-AGAINST")
            reason = (
                f"No primary laboratory measurement of this transition exists in Den "
                f"Hartog 2014 or Ruffoni 2014, and no catalogue is independent of FMW, so "
                f"the gf cannot be re-sourced. Against the pool's primary-sourced lines at "
                f"matched EP and REW it reads {x.residual_dex:+.3f} dex "
                f"({x.n_sigma:+.2f} sigma)"
                + (f"; its multiplet sibling(s) "
                   f"{', '.join(f'{w:.3f} ({s})' for w, s in zip(sib.wave_A, sib.ges_src))} "
                   f"give median A = {sib.A.median():.3f} vs this line's {x.A:.3f}."
                   if len(sib) else "; no primary-sourced line shares its multiplet.")
                + (" Consistent with the primary scale at this threshold -- the evidence "
                   "does not convict this line individually, which is NOT the same as "
                   "clearing it: the population it belongs to is offset."
                   if disp.startswith("NO-PRIMARY-NO") else ""))
            pg = pe = np.nan
            pn = ""
        recs.append(dict(
            wave_A=round(x.wave_A, 4), ep_eV=x.ep_eV, rew=round(x.rew, 3),
            ew_mA=round(x.ew_mA, 2), multiplet=x.multiplet,
            pool_source=x.ges_src, pool_loggf=x.ges_gf,
            vald_source=x.vald_src, vald_loggf=x.vald_gf,
            primary_source=pn, primary_loggf=pg, primary_e_loggf=pe,
            n_primary_multiplet_siblings=int(len(sib)),
            A_1D_LTE=round(x.A, 4), A_primary_scale=round(x.A_primary_scale, 4),
            residual_dex=round(x.residual_dex, 4), n_sigma=round(x.n_sigma, 2),
            disposition=disp, reason=reason))
    out = pd.DataFrame(recs).sort_values("wave_A").reset_index(drop=True)

    print(f"\n3. PER-LINE DISPOSITION (quarantine at >= {a.sigma:.1f} sigma)")
    print(f"{'line':>10} {'pool':11s} {'A':>7} {'A_pred':>7} {'resid':>7} {'sig':>6}  "
          f"{'sibs':>4}  disposition")
    for _, r in out.iterrows():
        print(f"{r.wave_A:10.4f} {r.pool_source:11s} {r.A_1D_LTE:7.3f} "
              f"{r.A_primary_scale:7.3f} {r.residual_dex:+7.3f} {r.n_sigma:+6.2f}  "
              f"{r.n_primary_multiplet_siblings:4d}  {r.disposition}")
    print("\n  " + "  ".join(f"{k}={v}" for k, v in out.disposition.value_counts().items()))

    # ── 3b. the POPULATION result, which is far stronger than any single line ─
    #
    # Per-line significance is the wrong summary to stop at: with rms 0.16 dex, a single
    # line needs a ~0.33 dex error before it clears 2 sigma, so a real +0.2 dex source
    # offset is individually unprovable on most lines. The SIGN is not: if FMW were on the
    # primary scale, each residual would be positive with probability 1/2.
    n_pos = int((dis.residual_dex > 0).sum())
    n_tot = len(dis)
    from math import comb
    p_sign = sum(comb(n_tot, k) for k in range(n_pos, n_tot + 1)) / 2 ** n_tot
    print(f"\n3b. POPULATION — {n_pos} of {n_tot} disputed lines read HIGH against the "
          f"primary scale (sign test p = {p_sign:.2g}).")
    print(f"    median residual {dis.residual_dex.median():+.3f} dex, "
          f"range {dis.residual_dex.min():+.3f} to {dis.residual_dex.max():+.3f}.")
    print(f"    So the SOURCE is offset with high confidence even where individual lines "
          f"cannot be convicted.")
    print(f"    threshold sensitivity — quarantined at 1.5/2.0/2.5/3.0 sigma: "
          + ", ".join(f"{int((dis.n_sigma >= s).sum())}" for s in (1.5, 2.0, 2.5, 3.0)))

    # ── 4. what it does to the product ───────────────────────────────────────
    quar = out[out.disposition == "QUARANTINED-SCALE-EVIDENCE"].wave_A.tolist()
    before = t.A.median()
    after = t[~t.wave_A.isin(quar)].A.median()
    print(f"\n4. FE IR 1D-LTE PRODUCT, disputed set quarantined")
    print(f"   before : n={len(t):3d}  median A(Fe I) = {before:.4f}")
    print(f"   after  : n={len(t) - len(quar):3d}  median A(Fe I) = {after:.4f}   "
          f"({after - before:+.4f} dex)")
    print(f"   primary-sourced lines only : n={int(t.is_primary.sum())}  "
          f"median A(Fe I) = {t[t.is_primary].A.median():.4f}")

    o = Path(a.out)
    o.mkdir(parents=True, exist_ok=True)
    out.to_csv(o / "FeI_IR_gf_disposition.csv", index=False)
    summary = dict(
        ticket="RYA-780", n_measured=len(t), n_primary=int(t.is_primary.sum()),
        n_disputed=len(out),
        dispositions=out.disposition.value_counts().to_dict(),
        # ROUNDED, deliberately. These are diagnostics quoted to 3-4 dp in the report, and
        # a least-squares solve differs in the last few ULPs between BLAS builds -- the Mac
        # and Sirius disagreed at 1e-14, which is enough to make a committed artifact churn
        # on every cross-machine regeneration. That is the RYA-390 rya390_co_validation.json
        # failure mode and the reason four RYA-761 exact-float tests are red on the Mac
        # today. Six decimals is far more precision than the numbers carry.
        primary_scale=dict(intercept=round(beta[0], 6), ep=round(beta[1], 6),
                           rew=round(beta[2], 6), rms=round(rms, 6), n=npx),
        scale_validation={nm: dict(n=len(gp), median_delta=round(float(gp.delta.median()), 6))
                          for nm, gp in val.groupby("primary")} if len(val) else {},
        median_A_before=round(float(before), 6), median_A_after=round(float(after), 6),
        median_A_primary_only=round(float(t[t.is_primary].A.median()), 6),
        n_with_primary_measurement=int((out.primary_source != "").sum()),
        population_sign_test=dict(n_positive=n_pos, n=n_tot, p=round(p_sign, 12)))
    (o / "FeI_IR_gf_disposition_summary.json").write_text(
        json.dumps(summary, indent=2, default=float) + "\n")
    print(f"\n   wrote {o / 'FeI_IR_gf_disposition.csv'}")
    print(f"   wrote {o / 'FeI_IR_gf_disposition_summary.json'}")


if __name__ == "__main__":
    main()
