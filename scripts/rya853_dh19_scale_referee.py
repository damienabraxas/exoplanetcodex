#!/usr/bin/env python3
"""
RYA-853 scope 3 — is our Fe II gf scale legitimate, or solar-fitted?
====================================================================
RYA-852 found our Fe II pool sits above NIST, and that the current scale is what makes the
solar Fe I - Fe II ionization balance come out right. Two readings, with opposite
consequences for whether the balance — which gates log g — is an INDEPENDENT check or a
CIRCULAR one:

  * MINE (RYA-852): the pool is partly solar-fitted (one member is Meléndez & Barbuy 2009
    outright, and MB09 is built partly by reverse solar analysis), so a scale that makes the
    solar balance work is exactly what you would expect. CIRCULAR.
  * RYA-853 counter-evidence: Den Hartog 2019 reports NIST's Fe II gf sitting ~0.1 dex BELOW
    modern pure-lab values. Then sitting above NIST is correct and NIST is the outlier.
    LEGITIMATE.

THE REFEREE is Den Hartog et al. 2019 (ApJS 243, 33): branching fractions x laser-induced
fluorescence lifetimes. Pure laboratory, NO solar normalisation anywhere in it — which is
the whole point, because a referee that touched the solar spectrum could not settle a
question about solar fitting.

    ours - DH ~ 0        -> our scale IS the pure-lab scale; NIST is low; the balance is
                            an independent check and my solar-fitting hypothesis DIES.
    ours - DH ~ +0.1-0.2 -> we sit above pure lab too; the hypothesis LIVES and the
                            balance may be circular.

and the third leg tests the counter-evidence on its own terms:

    NIST - DH < 0        -> NIST really is low, as Den Hartog 2019 states.

🔴 WHY "OURS" IS A FROZEN SNAPSHOT AND NOT `canonical_gf`
---------------------------------------------------------
[RYA-945] (`b545dc6`, 2026-08-21) ingested DH19 Table 6 INTO `canonical_gf`. Correct for the
line list; fatal for this experiment. Run against live canonical_gf on 2026-08-27 this
script reported

    ours - DH  +0.000   CI [+0.000, +0.000]   sd 0.000   n=10
    VERDICT: LEGITIMATE — our scale IS the pure-lab scale ... REFUTED

because all ten overlap lines had come to cite `PRIMARY LAB DenHartog2019`. The referee was
scoring its own values, and the bootstrap CI added to stop exactly this kind of overreading
did not catch it — a zero-width CI passes every "is it well determined?" test there is.

So "ours" is `fe2_pre945_scale_snapshot.csv`: the Fe II scale as it stood at `b545dc6^`,
which is the scale that underwrote the ionization balance and the only one the referee can
test. Post-945 the overlap lines ARE Den Hartog's by adoption, and asking whether they agree
with Den Hartog is not a question. `--source live` reproduces the degenerate run on purpose,
to show the guard firing.

⚠️ THE OVERLAP DOES NOT CONTAIN THE ARBITER LINES. DH19 stops at 4584 A and the three Fe II
arbiter lines (6147.734 / 6238.386 / 6247.557) are redward of it. This referees the SCALE of
the pool, not those three lines.

⚠️ THE BALANCE THIS FEEDS RESTS ON PRE-CONTINUUM-FIX MEASUREMENTS. See `BALANCE_CAVEAT`.

Usage:
    python3 scripts/rya853_dh19_scale_referee.py                # frozen pre-945 scale
    python3 scripts/rya853_dh19_scale_referee.py --source live  # shows the guard fire
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

CANONICAL_GF = ROOT / "data" / "linelists" / "canonical_gf.csv"
FROZEN = ROOT / "data" / "reference" / "fe_gf_lab" / "fe2_pre945_scale_snapshot.csv"
DH19_TABLE = ROOT / "data" / "reference" / "fe_gf_lab" / "fe2_lab_loggf_dh19.csv"
OUT = ROOT / "data" / "results" / "rya853"

DH19_CITE = "Den Hartog et al. 2019, ApJS 243, 33 (BF x LIF lifetimes; pure lab)"
DH19_VIZIER = "J/ApJS/243/33"

#: EP guard on BOTH sides (RYA-780/846/852). Wavelength alone picks the wrong level.
WAVE_TOL_A = 0.05
EP_TOL_EV = 0.05

#: The referee's own name, in every spelling the reference column uses for it.
REFEREE_MARKERS = ("denhartog2019", "den hartog 2019", "dh19", "2019apjs..243")

#: RYA-852's red Fe II pool. MEASURED here against live NIST rather than carried as the
#: literal +0.106 the previous revision hardcoded: 852/877/945 have all touched the pool
#: since that number was taken, so quoting it beside a freshly measured blue arm and
#: calling the difference band-dependence compares two different pools on two different
#: dates. (These eleven are outside DH19's span, so RYA-945 did not move them — verified.)
RYA852_RED_POOL_A = (5256.932, 5337.722, 5991.371, 6084.102, 6147.734, 6149.246,
                     6238.386, 6247.557, 6369.459, 6432.676, 6456.380)

#: Refuse a verdict below this many surviving lines, and refuse one from an estimator that
#: has collapsed. A spread of exactly zero is not agreement, it is the same number twice.
MIN_LINES_FOR_VERDICT = 5
DEGENERATE_CI_WIDTH = 1e-6

INDEPENDENCE_CAVEAT = (
    "GUARD 1 checks the reference STRING, and it can only catch a value we adopted from the "
    "referee directly. The 12 near-UV lines that carry this verdict are all plain 'VALD3', "
    "'single_source' — and VALD3 is an aggregator, so that label does not establish "
    "independence from Den Hartog, only that we did not ingest Den Hartog for them. The "
    "values are not copies (they differ by up to 0.250 dex), so this is not a laundered "
    "self-match; but a shared upstream ancestor cannot be excluded from what we record. "
    "Closing it needs VALD3's per-line source for these twelve."
)

BALANCE_CAVEAT = (
    "The balance figures this verdict is usually quoted beside — Fe I 7.586 / Fe II 7.568 "
    "(RYA-783, kpno_solar_atlas PROFILEFIT) — are dated 2026-08-16..18 and PREDATE every "
    "continuum fix: RYA-911/913 (08-19), RYA-1000/1006 --local-renorm (08-23), RYA-1026 "
    "prenormalised_guard and RYA-1030 (08-24). RYA-852's test still hardcodes 7.568. "
    "CORRECTED 2026-08-27: post-fix products DO exist — RYA-1045 (Fe II VIS, 08-25) and "
    "RYA-1052 (Fe II near-UV, 08-26) re-ran Fe II on the corrected holdings, and 15 Fe II "
    "products sit in data/products/solar/Fe.json with artifact_mtime 08-25/26. They are "
    "invisible in data/results/band_products because 67 of the 75 feed products carry "
    "provenance.copied_to=None and were never committed (RYA-1080). So the balance needs "
    "RE-DERIVING from the feed, not re-measuring — and the pairing must match tier and "
    "selector, which the old ungraded pair does not share with the current graded ones."
)


def _cites_referee(s: pd.Series) -> pd.Series:
    return s.astype(str).str.lower().str.contains("|".join(REFEREE_MARKERS),
                                                  regex=True, na=False)


def load_dh19() -> pd.DataFrame:
    """The referee, from the vendored machine-readable Table 6 (RYA-945, 131 lines).

    The previous revision carried ten values transcribed from the PDF into a literal tuple
    and said so as a limitation. RYA-945 vendored the whole table with a recorded PDF
    sha256, so the referee is now reproducible from the repo — and 131 lines instead of 10.
    """
    if not DH19_TABLE.exists():
        raise SystemExit(f"missing referee table {DH19_TABLE} — run RYA-945's fetch first")
    d = pd.read_csv(DH19_TABLE)
    d = d[d.species.astype(str).str.strip() == "Fe II"]
    return (d[["wavelength_air_A", "elo_eV", "loggf", "e_loggf_dex"]]
            .rename(columns={"elo_eV": "chi_eV", "loggf": "dh_loggf",
                             "e_loggf_dex": "dh_sigma_dex"})
            .sort_values("wavelength_air_A").reset_index(drop=True))


def ours(dh: pd.DataFrame, source: str) -> pd.DataFrame:
    """Our stored Fe II log gf for the overlap lines, EP-matched, ambiguity refused."""
    if source == "frozen":
        if not FROZEN.exists():
            raise SystemExit(
                f"missing {FROZEN} — run scripts/rya853_freeze_pre945_fe2_scale.py")
        c = pd.read_csv(FROZEN, low_memory=False)
    else:
        c = pd.read_csv(CANONICAL_GF, comment="#", low_memory=False)
    c = c[c.species.astype(str) == "Fe II"]

    rows = []
    for _, d in dh.iterrows():
        rec = {"wavelength_air_A": d.wavelength_air_A, "chi_eV": d.chi_eV,
               "dh_loggf": d.dh_loggf, "dh_sigma_dex": d.dh_sigma_dex}
        cm = c[(np.abs(c.wavelength_air_A - d.wavelength_air_A) <= WAVE_TOL_A)
               & (np.abs(c.excitation_potential_eV - d.chi_eV) <= EP_TOL_EV)]
        rec["n_candidates"] = int(len(cm))
        if len(cm) == 1:
            rec["our_loggf"] = float(cm.iloc[0].log_gf)
            rec["our_ref"] = str(cm.iloc[0].loggf_reference)
        elif len(cm) > 1:
            rec["ambiguous"] = True
        rows.append(rec)
    o = pd.DataFrame(rows)
    return o[o.n_candidates == 1].reset_index(drop=True)


def nist_for(lo: float, hi: float) -> pd.DataFrame:
    """NIST ASD Fe II over a span, in AIR, with Ei kept so the caller can EP-match."""
    from astroquery.nist import Nist
    import astropy.units as u

    t = Nist.query((lo - 5) * u.AA, (hi + 5) * u.AA, linename="Fe II",
                   wavelength_type="vac+air")
    gcol = next((c for c in t.colnames if c.startswith("gi")), None)
    eicol = next((c for c in t.colnames if c.startswith("Ei")), None)
    cand = []
    for r in t:
        def _f(col):
            try:
                return float(str(r[col]).strip())
            except Exception:
                return float("nan")
        w = _f("Observed")
        if not np.isfinite(w):
            w = _f("Ritz")
        fik = _f("fik")
        try:
            gi = float(str(r[gcol]).split("-")[0].strip()) if gcol else np.nan
        except Exception:
            gi = np.nan
        try:
            ei = float(str(r[eicol]).split("-")[0].strip()) if eicol else np.nan
        except Exception:
            ei = np.nan
        if not (np.isfinite(w) and np.isfinite(gi) and np.isfinite(fik) and gi * fik > 0):
            continue
        cand.append({"w": w, "ei": ei, "loggf": float(np.log10(gi * fik)),
                     "acc": str(r["Acc."]).strip()})
    return pd.DataFrame(cand)


def attach_nist(o: pd.DataFrame) -> pd.DataFrame:
    """NIST for the same lines, EP-matched on the NIST side too (RYA-853 scope 4)."""
    cd = nist_for(float(o.wavelength_air_A.min()), float(o.wavelength_air_A.max()))
    out = []
    for _, d in o.iterrows():
        rec = {"wavelength_air_A": d.wavelength_air_A}
        if len(cd):
            m = cd[(np.abs(cd.w - d.wavelength_air_A) <= WAVE_TOL_A)
                   & (np.abs(cd.ei - d.chi_eV) <= EP_TOL_EV)]
            rec["n_nist_candidates"] = int(len(m))
            if len(m) == 1:
                rec["nist_loggf"] = float(m.iloc[0].loggf)
                rec["nist_acc"] = str(m.iloc[0].acc)
        out.append(rec)
    return o.merge(pd.DataFrame(out), on="wavelength_air_A", how="left")


def stat(v: pd.Series, seed: int = 0, n_boot: int = 20000):
    """Median WITH its uncertainty. A median of ten values spanning 0.6 dex is not a number
    you can compare against a 0.05 threshold without saying how well it is determined —
    that is how a test that cannot discriminate gets read as a verdict."""
    v = v.dropna()
    if not len(v):
        return None
    a = v.to_numpy(float)
    rng = np.random.default_rng(seed)
    boot = np.array([np.median(rng.choice(a, len(a))) for _ in range(n_boot)])
    ci = [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))]
    return {"n": int(len(a)), "median": float(np.median(a)), "mean": float(a.mean()),
            "std": float(a.std(ddof=1)) if len(a) > 1 else 0.0,
            "sem": float(a.std(ddof=1) / np.sqrt(len(a))) if len(a) > 1 else 0.0,
            "ci95": ci, "ci_width": float(ci[1] - ci[0]),
            "min": float(a.min()), "max": float(a.max())}


def red_arm(source: str):
    """RYA-852's red pool measured against live NIST, rather than carried as a literal."""
    c = (pd.read_csv(FROZEN, low_memory=False) if source == "frozen"
         else pd.read_csv(CANONICAL_GF, comment="#", low_memory=False))
    c = c[c.species.astype(str) == "Fe II"]
    cd = nist_for(min(RYA852_RED_POOL_A), max(RYA852_RED_POOL_A))
    rows, ambiguous = [], []
    for w in RYA852_RED_POOL_A:
        cm = c[np.abs(c.wavelength_air_A - w) <= 0.01]
        if len(cm) != 1:
            ambiguous.append((w, "ours", int(len(cm))))
            continue
        ep = float(cm.iloc[0].excitation_potential_eV)
        m = cd[(np.abs(cd.w - w) <= WAVE_TOL_A) & (np.abs(cd.ei - ep) <= EP_TOL_EV)]
        if len(m) != 1:
            ambiguous.append((w, "nist", int(len(m))))
            continue
        rows.append({"wavelength_air_A": w, "our_loggf": float(cm.iloc[0].log_gf),
                     "our_ref": str(cm.iloc[0].loggf_reference),
                     "nist_loggf": float(m.iloc[0].loggf),
                     "ours_minus_nist": float(cm.iloc[0].log_gf) - float(m.iloc[0].loggf)})
    return pd.DataFrame(rows), ambiguous


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", choices=("frozen", "live"), default="frozen",
                    help="'frozen' = the pre-945 scale (the real test). 'live' = today's "
                         "canonical_gf, which CONTAINS the referee; kept so the guard can "
                         "be seen firing.")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    dh = load_dh19()
    print(f"[referee] {DH19_CITE}")
    print(f"[referee] {len(dh)} lines from {DH19_TABLE.relative_to(ROOT)}, "
          f"{dh.wavelength_air_A.min():.1f}-{dh.wavelength_air_A.max():.1f} A")
    print(f"[ours]    source = {args.source}"
          + (f"  ({FROZEN.name}, b545dc6^)" if args.source == "frozen"
             else "  (canonical_gf.csv — POST-945, contains the referee)"))

    o = ours(dh, args.source)

    # ── GUARD 1: the referee must not appear on both sides ───────────────────────
    self_ref = _cites_referee(o.our_ref) if "our_ref" in o else pd.Series(dtype=bool)
    n_self = int(self_ref.sum())
    o["self_referential"] = self_ref
    scored = o[~self_ref].copy()
    print(f"\n[guard 1] overlap lines EP-matched: {len(o)}")
    print(f"[guard 1] excluded as SELF-REFERENTIAL (our value cites the referee): {n_self}")
    print(f"[guard 1] scored: {len(scored)}")

    if len(scored):
        scored = attach_nist(scored)
        scored["ours_minus_dh"] = scored.our_loggf - scored.dh_loggf
        scored["nist_minus_dh"] = scored.get("nist_loggf") - scored.dh_loggf
        scored["ours_minus_nist"] = scored.our_loggf - scored.get("nist_loggf")
        scored["band"] = np.where(scored.wavelength_air_A >= 3300, "blue", "near-UV")

        print(f"\n{'wave':>9}{'chi':>7}{'DH':>8}{'ours':>8}{'Δ ours-DH':>11}"
              f"{'NIST':>8}{'Δ NIST-DH':>11}  band     our source")
        for _, r in scored.iterrows():
            print(f"{r.wavelength_air_A:9.3f}{r.chi_eV:7.3f}{r.dh_loggf:8.2f}"
                  f"{r.our_loggf:8.3f}{r.ours_minus_dh:+11.3f}"
                  f"{r.get('nist_loggf', np.nan):8.3f}"
                  f"{r.get('nist_minus_dh', np.nan):+11.3f}  {r.band:<8} "
                  f"{str(r.our_ref)[:20]}")

    s_ours = stat(scored.ours_minus_dh) if len(scored) else None
    s_nist = stat(scored.nist_minus_dh) if len(scored) else None
    s_on = stat(scored.ours_minus_nist) if len(scored) else None
    per_band = {b: stat(g.ours_minus_dh) for b, g in scored.groupby("band")} \
        if len(scored) else {}

    print(f"\n=== the referee (median, 95% bootstrap CI) ===")
    for lab, st in (("ours - DH", s_ours), ("NIST - DH", s_nist),
                    ("ours - NIST (same lines)", s_on)):
        if st:
            print(f"  {lab:<26} {st['median']:+.3f}  CI [{st['ci95'][0]:+.3f}, "
                  f"{st['ci95'][1]:+.3f}]  sd {st['std']:.3f}  n={st['n']}")
    for b, st in per_band.items():
        if st:
            print(f"    ours - DH / {b:<8}   {st['median']:+.3f}  CI "
                  f"[{st['ci95'][0]:+.3f}, {st['ci95'][1]:+.3f}]  sd {st['std']:.3f}  "
                  f"n={st['n']}")

    # ── GUARD 2: refuse a verdict from a collapsed estimator ─────────────────────
    refused = None
    if s_ours is None:
        refused = "no line survived the self-reference guard — there is nothing to score"
    elif s_ours["n"] < MIN_LINES_FOR_VERDICT:
        refused = (f"only {s_ours['n']} line(s) survived the self-reference guard, "
                   f"below the {MIN_LINES_FOR_VERDICT}-line floor")
    elif s_ours["ci_width"] < DEGENERATE_CI_WIDTH or s_ours["std"] == 0.0:
        refused = (f"the estimator has COLLAPSED: sd {s_ours['std']:.4f}, CI width "
                   f"{s_ours['ci_width']:.2e}. A spread of exactly zero is not agreement, "
                   f"it is the same number on both sides")

    verdict, reasoning = "INCONCLUSIVE", ""
    if refused:
        verdict = "UNUSABLE — the comparison is not a test"
        reasoning = refused
        print(f"\n=== 🔴 GUARD 2 FIRED ===\n  {refused}")
    else:
        med = s_ours["median"]
        lo, hi = s_ours["ci95"]
        covers_zero = lo <= 0.0 <= hi
        covers_offset = lo <= 0.13 <= hi
        if covers_zero and covers_offset:
            verdict = "INCONCLUSIVE — the overlap cannot separate the two readings"
            reasoning = (f"ours - DH is {med:+.3f} with a 95% CI of [{lo:+.3f}, {hi:+.3f}], "
                         f"which covers BOTH the pure-lab prediction (~0) and the "
                         f"solar-fitted one (~+0.13). {s_ours['n']} lines spanning "
                         f"{s_ours['min']:+.3f}..{s_ours['max']:+.3f} dex is not enough "
                         f"to choose.")
        elif not covers_offset and abs(med) <= 0.05:
            verdict = "LEGITIMATE — our scale IS the pure-lab scale"
            reasoning = (f"ours sits {med:+.3f} from Den Hartog's pure-lab values with a CI "
                         f"of [{lo:+.3f}, {hi:+.3f}] that EXCLUDES the solar-fitted "
                         f"prediction (~+0.13). Sitting above NIST is then NIST being low, "
                         f"as Den Hartog 2019 reports; the solar Fe I - Fe II balance is an "
                         f"INDEPENDENT check and RYA-852's hypothesis is REFUTED.")
        elif not covers_zero and med >= 0.08:
            verdict = "HYPOTHESIS LIVES — we sit above pure lab too"
            reasoning = (f"ours sits {med:+.3f} ABOVE pure lab, CI [{lo:+.3f}, {hi:+.3f}] "
                         f"excluding zero, so the offset is not NIST being low. A scale "
                         f"above pure lab that happens to make the solar balance work is "
                         f"what RYA-852 suspected, and the balance may be CIRCULAR.")
        else:
            reasoning = (f"ours sits {med:+.3f} from pure lab, CI [{lo:+.3f}, {hi:+.3f}] — "
                         f"between the two predictions, so neither reading is clean.")

    # ── the band-dependence finding, both arms MEASURED on the same source ───────
    band_flip, red = None, pd.DataFrame()
    try:
        red, red_ambig = red_arm(args.source)
        s_red = stat(red.ours_minus_nist) if len(red) else None
        blue_set = scored[scored.band == "blue"] if len(scored) else pd.DataFrame()
        s_blue = stat(blue_set.ours_minus_nist) if len(blue_set) else None
        if s_red and s_blue:
            band_flip = {
                "blue_overlap_ours_minus_nist": s_blue["median"],
                "blue_n": s_blue["n"],
                "red_pool_ours_minus_nist": s_red["median"],
                "red_n": s_red["n"],
                "red_pool_A": [min(RYA852_RED_POOL_A), max(RYA852_RED_POOL_A)],
                "swing_dex": s_red["median"] - s_blue["median"],
                "sign_flips": bool(s_blue["median"] * s_red["median"] < 0),
                "rya852_reported_red_offset": 0.106,
                "ambiguous_lines_dropped": red_ambig,
                "note": ("both arms measured against live NIST on the same source, on this "
                         "run. RYA-852's +0.106 is shown for comparison only — it was taken "
                         "before 852/877/945 touched the pool."),
            }
            print(f"\n=== 🔴 the offset against NIST is BAND-DEPENDENT ===")
            print(f"  blue overlap  ours - NIST = {s_blue['median']:+.3f}  (n={s_blue['n']})")
            print(f"  red pool      ours - NIST = {s_red['median']:+.3f}  (n={s_red['n']}, "
                  f"{min(RYA852_RED_POOL_A):.0f}-{max(RYA852_RED_POOL_A):.0f} A)")
            print(f"  swing {s_red['median'] - s_blue['median']:+.3f} dex"
                  f"{', SIGN FLIPS' if band_flip['sign_flips'] else ''}")
            print(f"  (RYA-852 reported +0.106 for the red arm on the pre-877 pool)")
            if red_ambig:
                print(f"  dropped as ambiguous, not judged: {red_ambig}")
    except Exception as e:
        print(f"\n  band-dependence leg skipped ({type(e).__name__}: {e})")

    print(f"\n=== VERDICT: {verdict} ===")
    print(f"  {reasoning}")
    print(f"\n  ⚠️ This referees the SCALE. The three arbiter lines (6147.734 / 6238.386 /"
          f"\n     6247.557) are REDWARD of DH19's ceiling and are NOT refereed here.")
    print(f"\n  ⚠️ {INDEPENDENCE_CAVEAT}")
    print(f"\n  ⚠️ {BALANCE_CAVEAT}")

    tag = "" if args.source == "frozen" else "_live"
    o.to_csv(OUT / f"rya853_dh19_referee_overlap{tag}.csv", index=False)
    if len(scored):
        scored.to_csv(OUT / f"rya853_dh19_referee{tag}.csv", index=False)
    if len(red):
        red.to_csv(OUT / f"rya853_red_arm_vs_nist{tag}.csv", index=False)
    (OUT / f"rya853_dh19_referee{tag}.json").write_text(json.dumps({
        "ticket": "RYA-853 scope 3",
        "referee": DH19_CITE, "vizier_id": DH19_VIZIER,
        "referee_source": str(DH19_TABLE.relative_to(ROOT)),
        "referee_n_lines": int(len(dh)),
        "ours_source": args.source,
        "ours_source_file": str((FROZEN if args.source == "frozen"
                                 else CANONICAL_GF).relative_to(ROOT)),
        "n_overlap_ep_matched": int(len(o)),
        "n_excluded_self_referential": n_self,
        "n_scored": int(len(scored)),
        "scored_span_A": ([float(scored.wavelength_air_A.min()),
                           float(scored.wavelength_air_A.max())] if len(scored) else None),
        "ours_minus_dh": s_ours, "nist_minus_dh": s_nist,
        "ours_minus_nist_same_lines": s_on,
        "ours_minus_dh_per_band": per_band,
        "verdict": verdict, "reasoning": reasoning,
        "guard_refusal": refused,
        "band_dependence": band_flip,
        "caveat_arbiter_lines": (
            "the overlap does NOT contain the three Fe II arbiter lines, which are redward "
            "of DH19's 4584 A ceiling; this refs the pool SCALE, not those lines"),
        "caveat_independence_is_string_deep_only": INDEPENDENCE_CAVEAT,
        "caveat_balance_is_pre_continuum_fix": BALANCE_CAVEAT,
    }, indent=2))
    print(f"\n[out] {OUT}/rya853_dh19_referee{tag}.json")


if __name__ == "__main__":
    main()
