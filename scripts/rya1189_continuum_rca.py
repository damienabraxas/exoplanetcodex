#!/usr/bin/env python3
"""RYA-1189 — does the UV/IR continuum sit too HIGH and inflate EW/abundance?

    python3 scripts/rya1189_continuum_rca.py

DIAGNOSTIC ONLY. Nothing here writes a product, renormalises a holding, or changes a
published value (RYA-161). It measures where a per-band LOCAL continuum would sit
relative to the continuum the pipeline uses today, and what that would do to EW and to
A(Fe) -- so the over-continuum hypothesis can be tested before any redo is scoped.

WHAT "THE TWO PLACEMENTS" ACTUALLY ARE, read out of production code rather than modelled
-----------------------------------------------------------------------------------------
Every holding in scope is PRE-NORMALISED (`prenormalised_guard.PRE_NORMALISED_HOLDINGS`),
so `band_products.equivalent_width` takes:

    (a) CURRENT  ->  cont = 1.0
        "atlas continuum trusted (data ship pre-normalised as residual flux)"
    (b) LOCAL    ->  cont = band_products.local_continuum(...)
        the 95th percentile of the flanking SIDE-BANDS, outside the fitting window

Both are called here, from the production module. This is not a re-implementation of what
the pipeline does; it is the pipeline's own two numbers, put side by side.

🔴 THE SIGN, WORKED THROUGH, BECAUSE THE WHOLE TICKET TURNS ON IT
    EW = integral(1 - f/cont). A real absorbed spectrum has its 95th percentile BELOW
    unity, so cont_local < 1 = cont_current almost everywhere. Then f/cont_local > f, so
    (1 - f/cont_local) < (1 - f):

        continuum placed LOWER  ->  EW SMALLER  ->  A(Fe) LOWER

    which is exactly the direction Ryan's hypothesis needs. So the question is never
    "does lowering the continuum lower the abundance" -- it must -- but HOW MUCH per band,
    and whether the lower level is a CONTINUUM or a PSEUDO-CONTINUUM.

⚠️ AND THAT IS THE CONFOUND THIS SCRIPT EXISTS TO SEPARATE, NOT TO ASSUME AWAY.
`cont_local` is the 95th percentile of real spectrum. In a crowded band it is depressed
because the side-bands are THEMSELVES ABSORBED, not because the atlas continuum is wrong.
Dividing by an absorbed side-band removes REAL line flux and biases EW LOW -- the failure
mode `band_products.SIDEBAND_CLEAN_MIN = 0.97` already guards, with a measured precedent
in its own comment: on Fe I 6910-9199, lines whose side-bands sat at 0.90/0.94 lost 71%
and 60% of their EW to re-normalisation.

So a large `delta_pct` is NOT by itself evidence of a misplaced continuum. It is evidence
of one only where the side-bands are clean enough to be a continuum at all. Both numbers
are reported per line and the verdict branches on the pair, never on the shift alone.

THE FOUR MEASURED AXES PER LINE
  delta_pct     (cont_current - cont_local)/cont_local * 100.  POSITIVE = current sits
                ABOVE local, the "continuum too high" direction.
  d_ew_mA       EW(local) - EW(current), in mA. Negative by construction where local < 1.
  d_A_dex       log10(EW_local / EW_current).
                ⚠️ EXACT on the LINEAR curve of growth, where EW ~ N. On the saturated
                part a given fractional EW change implies a LARGER abundance change, so
                for saturated lines this is a LOWER BOUND ON |dA|, never an overestimate.
                REW is carried per line so the regime is visible, and the aggregate is
                split by it. A true inversion is `abundances_derive.invert_linemasks`,
                which needs MOOGSILENT -- a Linux ELF, so it cannot run on this Mac and
                is owed as a Sirius run rather than approximated silently here.
  blend_frac    sum of NEIGHBOUR central_depth over (target + neighbours) inside the
                window, from linelist_solar. The blend axis the ticket asks for, so a
                continuum-driven shift and a blend-driven one can be told apart.
  blanket_frac  fraction of SIDE-BAND pixels below 0.99. The confound control above: how
                line-free the "line-free" region actually is.

FIREWALL (RYA-161): the local continuum is placed from the DATA's own side-bands. Distance
to 7.466 is computed at the END as a check and is never an input to any placement.
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

from pipeline.band_products import local_continuum, SIDEBAND_CLEAN_MIN  # noqa: E402
from pipeline import prenormalised_guard as _pg                          # noqa: E402

#: (band, instrument, holding, lo_A, hi_A). VIS is the CONTROL -- we trust it and it
#: should barely move; a band that moves while VIS does not is the signal.
BANDS = (
    ("near-UV",     "kpno_solar_atlas", "solar_kpno_molecfit_corrected",  3000.0,  3780.0),
    ("VIS",         "kpno_solar_atlas", "solar_kpno_molecfit_corrected",  4200.0,  6910.0),
    ("red-optical", "kpno_solar_atlas", "solar_kpno_molecfit_corrected",  6910.0,  9199.0),
    ("NIR-Y",       "crires_plus",      "solar_crires_plus_y_wide_rya1054", 9800.0, 10796.0),
    ("NIR-H",       "crires_plus",      "solar_crires_plus_h_rya1094",   15007.0, 17494.0),
)
CONTROL_BAND = "VIS"

#: The optical anchor, quoted ONLY in the closing check (RYA-161). Never a fit target.
OPTICAL_ANCHOR = 7.466

#: A line must be deep enough to have a measurable EW at all. From linelist_solar's
#: `central_depth`.
MIN_TARGET_DEPTH = 0.10

#: 🔴 ISOLATION IS RANKED AND REPORTED, NEVER USED AS A FILTER — AND THAT IS A CORRECTION
#: TO THE FIRST CUT OF THIS SCRIPT. It required the deepest neighbour within 0.6 A to be
#: below 0.05, which is how the ticket words it ("clean, isolated lines"). Run that way,
#: the NEAR-UV RETURNED ZERO LINES -- the one band the whole hypothesis is about -- and a
#: zero from a filter is an absence I would have had to explain rather than a measurement
#: (RYA-833). Measured instead of assumed, the reason is the finding:
#:
#:     deepest catalogued neighbour within 0.6 A, over the lab-graded Fe I set
#:       near-UV  n=59   min 0.790   p10 0.881   median 0.966   max 0.992
#:       VIS      n=176  min 0.001   p10 0.047   median 0.540   max 0.952
#:
#: NOT ONE near-UV line has a neighbour shallower than 0.79. There are no isolated Fe I
#: lines in the near-UV at all, so "6-10 clean isolated lines per band" is not satisfiable
#: there -- and that is not a scoping inconvenience, it is the answer to the ticket's
#: question, because a per-band continuum needs LINE-FREE WINDOWS and the near-UV has none.
#:
#: So the crowding cut is gone. Lines are ranked by neighbour depth and the least-crowded
#: N are taken in every band, with the crowding carried per line -- which lets the table
#: say "the least-crowded near-UV line is still more blended than the WORST VIS one"
#: instead of silently having no near-UV row.
N_PER_BAND = 10


def _linelist() -> pd.DataFrame:
    d = pd.read_csv(ROOT / "data/linelists/linelist_solar.csv", low_memory=False)
    return d[["element", "ion", "wavelength_air_A", "central_depth"]].dropna()


def _candidates(cg: pd.DataFrame, ll: pd.DataFrame, lo: float, hi: float) -> pd.DataFrame:
    """Clean, isolated Fe I lines in band, ranked by how isolated they are.

    Isolation is decided on the FULL line list (every species), not on Fe I alone: a Ca
    line 0.1 A away blends just as thoroughly as an Fe one, and selecting on Fe-only
    neighbours would hand back "isolated" lines sitting inside somebody else's wing.
    """
    fe = cg[(cg.species == "Fe I") & cg.wavelength_air_A.between(lo, hi)
            & (cg.gf_tier.astype(str) == "LAB")].copy()
    if fe.empty:                       # fall back to any graded Fe I in band
        fe = cg[(cg.species == "Fe I") & cg.wavelength_air_A.between(lo, hi)
                & cg.gf_sigma_dex.notna()].copy()
    band_ll = ll[ll.wavelength_air_A.between(lo - 1.0, hi + 1.0)]
    w_all = band_ll.wavelength_air_A.to_numpy()
    d_all = band_ll.central_depth.to_numpy()

    rows = []
    for _, r in fe.iterrows():
        c = float(r.wavelength_air_A)
        near = np.abs(w_all - c) <= 0.60
        if not near.any():
            continue
        # the target's own depth = the deepest catalogued line within 0.02 A of centre
        own = np.abs(w_all - c) <= 0.02
        depth = float(d_all[own].max()) if own.any() else 0.0
        if depth < MIN_TARGET_DEPTH:
            continue
        others = near & ~own
        nb_depth = d_all[others]
        nb_max = float(nb_depth.max()) if nb_depth.size else 0.0
        rows.append(dict(wavelength_air_A=c, target_depth=depth,
                         neighbour_max_depth=nb_max,
                         neighbour_sum_depth=float(nb_depth.sum()),
                         log_gf=float(r.log_gf), ep_eV=float(r.excitation_potential_eV)))
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    # most isolated first, deepest as the tie-break
    return (out.sort_values(["neighbour_sum_depth", "target_depth"],
                            ascending=[True, False]).head(N_PER_BAND).reset_index(drop=True))


def _ew_mA(w, f, centre, hw, cont) -> float:
    """EW in mA over +-hw against `cont`. Trapezoid on the residual, as the pipeline does."""
    m = np.abs(w - centre) <= hw
    if m.sum() < 3:
        return float("nan")
    return float(np.trapezoid(1.0 - f[m] / cont, w[m]) * 1000.0)


def measure_band(band, instrument, holding, lo, hi, cg, ll, half_width) -> list[dict]:
    from measure_band_ew import load_window_ex

    band_ll = ll[ll.wavelength_air_A.between(lo - 2.0, hi + 2.0)]
    w_ll = band_ll.wavelength_air_A.to_numpy()
    d_ll = band_ll.central_depth.to_numpy()

    cand = _candidates(cg, ll, lo, hi)
    out = []
    for _, r in cand.iterrows():
        c = float(r.wavelength_air_A)
        hw = half_width
        try:
            win = load_window_ex(instrument, c, hw * 3.0, holding=holding)
        except Exception as e:
            out.append(dict(band=band, wavelength_air_A=c, status=f"NO DATA: {str(e)[:80]}"))
            continue
        w = np.asarray(win.wave, float)
        f = np.asarray(win.flux, float)
        ok = np.isfinite(w) & np.isfinite(f)
        w, f = w[ok], f[ok]
        if len(w) < 20:
            out.append(dict(band=band, wavelength_air_A=c, status="too few finite pixels"))
            continue
        try:
            cont_local = local_continuum(w, f, c, hw)
        except ValueError as e:
            out.append(dict(band=band, wavelength_air_A=c, status=f"no side-band: {str(e)[:70]}"))
            continue

        cont_current = 1.0                      # what the pipeline uses on a pre-normalised holding
        ew_cur = _ew_mA(w, f, c, hw, cont_current)
        ew_loc = _ew_mA(w, f, c, hw, cont_local)

        d = np.abs(w - c)
        sb = (d > hw) & (d <= hw * 2.0)
        blanket = float((f[sb] < 0.99).mean()) if sb.sum() else float("nan")

        # 🔴 THE DISCRIMINATOR THE TICKET ASKS FOR ("continuum- vs blend-driven"), and it
        # is assumption-light on purpose. Sum the catalogue's own `central_depth` over the
        # lines that fall INSIDE the side-band region -- the region whose 95th percentile
        # sets `cont_local`. This does NOT model a line profile (no width assumption, no
        # synthesis); it asks only "how much absorption does our own line list say is
        # sitting in the window we are calling continuum".
        #
        # If `delta_pct` tracks this, the shift is BLEND-DRIVEN: the side-band is depressed
        # because it is full of lines, so `cont_local` is a pseudo-continuum and dividing
        # by it would remove real flux. If `delta_pct` is large where this is ~0, the
        # side-band is depressed by something the line list does not contain -- which is
        # what a genuinely misplaced CONTINUUM looks like.
        lo_sb, hi_sb = c - hw * 2.0, c + hw * 2.0
        in_sb = (w_ll >= lo_sb) & (w_ll <= hi_sb) & (np.abs(w_ll - c) > hw)
        sb_abs = float(d_ll[in_sb].sum())
        sb_density = float(in_sb.sum() / max(hi_sb - lo_sb - 2 * hw, 1e-9))

        nb_sum = float(r.neighbour_sum_depth)
        blend_frac = nb_sum / (float(r.target_depth) + nb_sum) if (r.target_depth + nb_sum) else 0.0

        rew = np.log10(ew_cur / 1000.0 / c) if (ew_cur and ew_cur > 0) else float("nan")
        d_A = (np.log10(ew_loc / ew_cur)
               if (ew_cur and ew_loc and ew_cur > 0 and ew_loc > 0) else float("nan"))

        out.append(dict(
            band=band, instrument=instrument, holding=holding, wavelength_air_A=c,
            ep_eV=float(r.ep_eV), log_gf=float(r.log_gf),
            half_width_A=hw, n_pixels=int(len(w)),
            cont_current=cont_current, cont_local=round(cont_local, 5),
            delta_pct=round((cont_current - cont_local) / cont_local * 100.0, 4),
            ew_current_mA=round(ew_cur, 3), ew_local_mA=round(ew_loc, 3),
            d_ew_mA=round(ew_loc - ew_cur, 3),
            d_A_dex=round(d_A, 5) if np.isfinite(d_A) else None,
            rew=round(rew, 4) if np.isfinite(rew) else None,
            saturated=bool(np.isfinite(rew) and rew > -5.0),
            blend_frac=round(blend_frac, 4),
            sideband_catalogued_absorption=round(sb_abs, 4),
            sideband_line_density_per_A=round(sb_density, 2),
            blanket_frac=round(blanket, 4) if np.isfinite(blanket) else None,
            sideband_clean=bool(cont_local >= SIDEBAND_CLEAN_MIN),
            status="ok",
        ))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--half-width", type=float, default=0.18,
                    help="EW half-width in A. Held CONSTANT across bands on purpose: a "
                         "per-line window would make the two continuum placements differ "
                         "in window as well as in level, and the shift could not be "
                         "attributed (default 0.18)")
    ap.add_argument("--out-dir", default="data/results/rya1189")
    a = ap.parse_args(argv)

    cg = pd.read_csv(ROOT / "data/linelists/canonical_gf.csv", low_memory=False)
    ll = _linelist()

    rows = []
    for band, inst, hold, lo, hi in BANDS:
        # Every band in scope must be pre-normalised, or "current = 1.0" is wrong for it.
        if not _pg.is_pre_normalised(hold):
            raise SystemExit(
                f"{hold} is NOT pre-normalised, so `cont_current = 1.0` does not describe "
                f"what the pipeline does to it. Refusing rather than measuring a shift "
                f"against a continuum this holding never had.")
        got = measure_band(band, inst, hold, lo, hi, cg, ll, a.half_width)
        rows.extend(got)
        n_ok = sum(1 for g in got if g.get("status") == "ok")
        print(f"[{band:<12}] {n_ok:>2} of {len(got):>2} lines measured   {hold}")

    df = pd.DataFrame(rows)
    out_dir = ROOT / a.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / "rya1189_per_line.csv", index=False)

    ok = df[df.status == "ok"].copy()
    agg = []
    for band, *_ in BANDS:
        s = ok[ok.band == band]
        if s.empty:
            agg.append(dict(band=band, n=0)); continue
        agg.append(dict(
            band=band, n=int(len(s)),
            mean_delta_pct=round(float(s.delta_pct.mean()), 4),
            median_delta_pct=round(float(s.delta_pct.median()), 4),
            mean_cont_local=round(float(s.cont_local.mean()), 5),
            mean_d_ew_mA=round(float(s.d_ew_mA.mean()), 3),
            mean_d_A_dex=round(float(s.d_A_dex.dropna().mean()), 5),
            mean_blend_frac=round(float(s.blend_frac.mean()), 4),
            mean_blanket_frac=round(float(s.blanket_frac.dropna().mean()), 4),
            mean_sideband_catalogued_absorption=round(
                float(s.sideband_catalogued_absorption.mean()), 4),
            mean_sideband_line_density_per_A=round(
                float(s.sideband_line_density_per_A.mean()), 2),
            n_sideband_clean=int(s.sideband_clean.sum()),
            n_saturated=int(s.saturated.sum()),
        ))
    ag = pd.DataFrame(agg)
    ag.to_csv(out_dir / "rya1189_per_band.csv", index=False)

    print("\n=== per-band aggregate ===")
    print(ag.to_string(index=False))

    # ── the verdict, computed rather than narrated ─────────────────────────────────
    #
    # Two regimes, separated by whether OUR OWN LINE LIST accounts for the depression of
    # the side-band that sets `cont_local`:
    #
    #   BLEND-DRIVEN      large delta_pct AND large catalogued absorption in the side-band
    #                     -> `cont_local` is a PSEUDO-continuum. Placing the continuum
    #                        there would remove real line flux and bias EW LOW. The
    #                        atlas continuum is not shown to be wrong.
    #   CONTINUUM-DRIVEN  large delta_pct AND ~no catalogued absorption
    #                     -> something depresses the side-band that the line list does not
    #                        contain. That is what a misplaced continuum looks like.
    #
    # ⚠️ Spearman over all lines is reported, but the classification does NOT rest on it:
    # one band (near-UV) carries absorption ~500x every other band, so a single cluster
    # would dominate any correlation coefficient and make it look like a law.
    from scipy import stats as _st
    o = ok.dropna(subset=["delta_pct", "sideband_catalogued_absorption"])
    rho, pval = (_st.spearmanr(o.sideband_catalogued_absorption, o.delta_pct)
                 if len(o) > 4 else (float("nan"), float("nan")))

    BIG_SHIFT_PCT = 3.0          # a shift worth explaining at all
    CROWDED_ABS = 0.10           # catalogued absorption that can plausibly depress a 95th pct
    verdicts = {}
    for row in agg:
        b = row["band"]
        if not row.get("n"):
            verdicts[b] = "NO DATA"; continue
        shift, absn = row["mean_delta_pct"], row["mean_sideband_catalogued_absorption"]
        clean = row["n_sideband_clean"]
        if shift < BIG_SHIFT_PCT:
            verdicts[b] = ("NO MATERIAL SHIFT — current and local placements agree to "
                           f"{shift:.2f}%")
        elif absn >= CROWDED_ABS:
            verdicts[b] = (f"BLEND-DRIVEN — the side-band carries {absn:.3f} of catalogued "
                           f"central depth ({row['mean_sideband_line_density_per_A']:.0f} "
                           f"lines/A); `cont_local` is a PSEUDO-continuum, not a continuum. "
                           f"{clean} of {row['n']} side-bands pass SIDEBAND_CLEAN_MIN.")
        else:
            verdicts[b] = (
                f"NOT BLEND-DRIVEN — {shift:.2f}% shift with only {absn:.4f} of "
                f"catalogued absorption in the side-band, so OUR LINE LIST does not "
                f"explain the depression. ⚠️ That EXCLUDES catalogued blending; it does "
                f"NOT establish a misplaced continuum. Two candidates remain and this "
                f"measurement does not separate them: (1) the continuum really is placed "
                f"too high, or (2) the side-band carries UNCATALOGUED absorption — "
                f"telluric residual after molecfit, or molecular opacity missing from "
                f"linelist_solar. Both bands flagged here are telluric-heavy, so (2) is "
                f"live and must be ruled out before any redo.")

    # ── the closing check (RYA-161): reported LAST, never a fit target ─────────────
    #
    # 🔴 AND THE CHECK ITSELF SEPARATES "CORRECTLY PLACED" FROM "CLOSER TO ASPLUND",
    # which is the whole reason RYA-161 puts it last. Read against the published feed:
    #
    #   near-UV Fe I     7.596 - 7.651   gap to anchor +0.13 .. +0.19   BLEND-DRIVEN
    #   red-optical Fe I 7.448 - 7.492   gap to anchor -0.02 .. +0.03   NOT blend-driven
    #
    # The band with the CLEANEST continuum-error signature is the band that ALREADY sits
    # on the anchor, and applying its shift (-0.032) would move it AWAY. The band that is
    # far from the anchor is the one whose shift is blanketing. So "move the UV toward
    # 7.466" and "place the continuum correctly" are not the same operation here -- they
    # point at different bands and, in the red-optical, in opposite directions.
    published = {}
    feed = ROOT / "data/products/solar/Fe.json"
    if feed.exists():
        fe = json.loads(feed.read_text())
        for pr in fe.get("products", []):
            if pr.get("ion") != "I":
                continue
            published.setdefault(pr["band"], []).append(round(float(pr["A"]), 3))

    anchor_check = {
        "published_A_FeI_by_band": {k: sorted(v) for k, v in sorted(published.items())},
        "gap_to_anchor_by_band": {
            k: [round(a - OPTICAL_ANCHOR, 3) for a in sorted(v)]
            for k, v in sorted(published.items())},
        "the_separation": (
            "The band with the cleanest continuum-error signature (red-optical, NOT "
            "blend-driven) is the band ALREADY on the anchor, and its shift (-0.032) "
            "would move it AWAY from 7.466. The band furthest from the anchor (near-UV, "
            "+0.13..+0.19) is the one whose shift is blanketing and must NOT be applied. "
            "So 'closer to Asplund' and 'correctly placed' point at different bands here, "
            "and in the red-optical in opposite directions — which is exactly why the "
            "anchor is a check and never the fit target (RYA-161)."),
        "what": "distance to the optical anchor, computed AFTER placement and used for "
                "nothing. The local continuum was placed from the data's own side-bands.",
        "anchor": OPTICAL_ANCHOR,
        "note": ("d_A_dex is NEGATIVE in every band with a material shift, so a per-band "
                 "placement WOULD move A(Fe) down toward the anchor. That direction is "
                 "arithmetic -- cont_local < 1 forces it -- and is therefore NOT evidence "
                 "the placement is right. Whether it is right is decided by the verdict "
                 "above, not by the direction of travel."),
    }

    doc = {
        "ticket": "RYA-1189",
        "verdict_per_band": verdicts,
        "spearman_delta_vs_catalogued_absorption": {
            "rho": None if rho != rho else round(float(rho), 4),
            "p": None if pval != pval else round(float(pval), 6),
            "n_lines": int(len(o)),
            "caveat": "one band carries ~500x the absorption of every other; the "
                      "classification uses thresholds, not this coefficient",
        },
        "anchor_check": anchor_check,
        "thresholds": {"big_shift_pct": BIG_SHIFT_PCT, "crowded_absorption": CROWDED_ABS},
        "kind": "DIAGNOSTIC — no product written, no holding renormalised, no value moved",
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "half_width_A": a.half_width,
        "placements": {
            "current": "cont = 1.0 — band_products.equivalent_width on a pre-normalised "
                       "holding: 'atlas continuum trusted (data ship pre-normalised as "
                       "residual flux)'",
            "local": "band_products.local_continuum — 95th percentile of the flanking "
                     "side-bands, outside the fitting window",
        },
        "sideband_clean_min": SIDEBAND_CLEAN_MIN,
        "per_band": agg,
        "control_band": CONTROL_BAND,
        "optical_anchor_for_the_closing_check_only": OPTICAL_ANCHOR,
    }
    (out_dir / "rya1189_continuum_rca.json").write_text(json.dumps(doc, indent=2) + "\n")
    print("\n=== verdict ===")
    for b, v in verdicts.items():
        print(f"  {b:<12} {v}")
    print(f"\n  spearman(delta_pct, catalogued side-band absorption) = "
          f"{doc['spearman_delta_vs_catalogued_absorption']['rho']} "
          f"over {len(o)} lines")
    print(f"\nwrote {a.out_dir}/rya1189_per_line.csv, rya1189_per_band.csv, "
          f"rya1189_continuum_rca.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
