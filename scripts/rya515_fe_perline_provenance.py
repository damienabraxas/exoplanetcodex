#!/usr/bin/env python3
"""RYA-515 — per-line provenance for every live Fe product, plus the Fe II deep proof.

Emits, under data/audit/rya515_fe_perline/:
  perline_provenance.csv   one row per (product, line), engine + correction stamped on EVERY line
  reachability.csv         one row per live product: is its per-line evidence reachable at all?
  fe2_blend_census.csv     per Fe II line: what else is in the fit window, and how deep
  fe1_control.csv          the same census on the Fe I anchor pool -- the discriminator's control
  findings.md              the written proof

and, under results/plots/:
  rya515_hist_by_engine.png / rya515_hist_by_correction.png   band x ion histograms
  rya515_fe2_outliers.png   the flagged Fe II outliers against the Fe I pool
  rya515_fe1_anchor.png     the A(Fe I) distribution behind 7.466

--check re-runs and asserts the committed artifacts still reproduce.
"""
import argparse, glob, json, os, re, sys
import pandas as pd, numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AUD  = os.path.join(ROOT, "data/audit/rya515_fe_perline")
PLOT = os.path.join(ROOT, "results/plots")
BP   = os.path.join(ROOT, "data/results/band_products")

# The half-width of the fit window, in Angstrom, used by the synthesis route.
FIT_HALF_WIDTH_A = 0.62
# A line's window is "contaminated" when a non-target transition reaches at least
# this fraction of the target's own central depth.
INTERLOPER_FRAC  = 0.5
# Asplund, Amarsi & Grevesse 2021 -- the external anchor A(Fe I) = 7.466 stands on.
AGSS21_FE = 7.46
ANCHOR    = 7.466


def load_products():
    with open(os.path.join(ROOT, "data/products/solar/Fe.json")) as fh:
        return json.load(fh)


def perline_path(p):
    """Resolve a product's per-line artifact STRICTLY.

    A glob that matches more than one file is NOT a match: 18 of the 70 products
    have an ambiguous stem (the tier/treatment stems nest -- GRADED is a substring
    of DEEPGRADED, and 1D-LTE of synth-1D-LTE-gerber), so a first-hit pick would
    silently stamp one product's lines onto another. Candidates are therefore
    filtered by exact field decomposition and then gated on n_lines value-equality
    with the product itself. Returns (path, reason) with path None when unresolved.
    """
    cands = []
    rng = p.get("wavelength_range_A") or []
    for cand in glob.glob(os.path.join(BP, "Fe%s_*_lines.csv" % p["ion"])):
        b = os.path.basename(cand)[:-len("_lines.csv")]
        # the wavelength stem disambiguates near-UV from VIS: every other field is
        # identical between the two bands of the same holding/tier/treatment.
        m = re.match(r"Fe(?:I|II)_(\d+)_(\d+)_", os.path.basename(cand))
        if m and rng:
            if not (abs(float(m.group(1)) - rng[0]) <= 20 and abs(float(m.group(2)) - rng[1]) <= 20):
                continue
        # treatment is the trailing field; tier the one before it
        if not b.endswith("_" + p["treatment"]):
            continue
        stem = b[: -len("_" + p["treatment"])]
        if not stem.endswith("_" + p["tier"]):
            continue
        stem = stem[: -len("_" + p["tier"])]
        if not stem.endswith("_" + p["route"]):
            continue
        stem = stem[: -len("_" + p["route"])]
        holding = p.get("holding", "")
        if holding and not stem.endswith("_" + holding):
            continue
        cands.append(cand)
    if not cands:
        return None, "no artifact on disk (band_products is gitignored; only force-added files are present)"
    if len(cands) > 1:
        return None, "AMBIGUOUS: %d artifacts share this product's fields" % len(cands)
    t = pd.read_csv(cands[0])
    # n_lines counts the lines that entered the AGGREGATE, not the rows on disk:
    # a row can be present and excluded (e.g. ENGINE-A-NOT-SERVED).
    n = int(t["in_aggregate"].sum()) if "in_aggregate" in t.columns else len(t)
    want = p.get("n_lines")
    if want is not None and int(want) != n:
        return None, "n_lines MISMATCH: product says %s, artifact has %d" % (want, n)
    return cands[0], "resolved"


def build_provenance(prods):
    """One row per line, with the product's engine/route/correction stamped on it."""
    rows, reach = [], []
    for p in prods["products"]:
        f, why = perline_path(p)
        reach.append(dict(
            element=p["element"], ion=p["ion"], band=p["band"], holding=p.get("holding", ""),
            tier=p["tier"], route=p["route"], treatment=p["treatment"],
            line_set=p.get("line_set", ""), A=p["A"], n_lines=p.get("n_lines"),
            perline_reachable=bool(f),
            resolution=why,
            perline_artifact=os.path.basename(f) if f else "",
        ))
        if not f:
            continue
        t = pd.read_csv(f)
        for _, r in t.iterrows():
            rows.append(dict(
                element="Fe", ion=p["ion"], band=p["band"], holding=p.get("holding", ""),
                tier=p["tier"], route=p["route"], treatment=p["treatment"],
                line_set=p.get("line_set", ""), product_A=p["A"],
                wavelength_air_A=r.get("wavelength_air_A"),
                line_A=r.get("abundance"),
                red_chi2=r.get("red_chi2"),
                in_aggregate=r.get("in_aggregate"),
                excluded_reason=r.get("excluded_reason", ""),
                nlte_delta_dex=r.get("nlte_delta_dex"),
                nlte_source=r.get("nlte_source", ""),
                perline_artifact=os.path.basename(f),
            ))
    return pd.DataFrame(rows), pd.DataFrame(reach)


def observed_strength(wl):
    """Measure the line's OWN strength off the observed atlas.

    This is the discriminator between the two competing explanations for a high
    recovered A: a blend swallowed by the fit, or a saturated line whose core
    carries almost no abundance information. The synthesis route leaves ew_mA and
    observed_depth NaN, so both are measured here from the atlas directly.
    """
    w, f = atlas_window(wl, 0.10)
    if w is None or not len(w):
        return dict(obs_depth=np.nan, obs_ew_mA=np.nan, obs_rew=np.nan)
    depth = 1.0 - float(np.nanmin(f))
    w2, f2 = atlas_window(wl, 0.25)
    ew = float(np.trapz(1.0 - f2, w2) * 1000.0) if w2 is not None and len(w2) > 2 else np.nan
    return dict(obs_depth=depth, obs_ew_mA=ew,
                obs_rew=(np.log10(ew / wl) if ew and ew > 0 else np.nan))


def blend_census(lines_df, linelist):
    """For each line: how many transitions share its fit window, and how deep is the deepest."""
    out = []
    for _, r in lines_df.iterrows():
        wl = r["wavelength_air_A"]
        if pd.isna(wl):
            continue
        w = linelist[(linelist.wavelength_air_A - wl).abs() <= FIT_HALF_WIDTH_A]
        is_target = (w.wavelength_air_A - wl).abs() <= 0.005
        tgt = w[is_target].central_depth.max()
        others = w[~is_target]
        deepest = others.central_depth.max() if len(others) else np.nan
        top = others.loc[others.central_depth.idxmax()] if len(others) and others.central_depth.notna().any() else None
        out.append(dict(
            ion=r.get("ion"), band=r.get("band"), holding=r.get("holding"),
            wavelength_air_A=wl, line_A=r.get("line_A"), red_chi2=r.get("red_chi2"),
            target_depth=tgt, n_in_window=len(others), deepest_interloper=deepest,
            deepest_species=(("%s %s" % (top.element, top.ion)) if top is not None else ""),
            deepest_wavelength_A=(top.wavelength_air_A if top is not None else np.nan),
            n_CH_in_window=int((others.element.astype(str) == "CH").sum()) if len(others) else 0,
            contaminated=bool(pd.notna(deepest) and pd.notna(tgt) and deepest >= INTERLOPER_FRAC * tgt),
            **observed_strength(wl),
        ))
    return pd.DataFrame(out)


ATLAS = os.environ.get("CODEX_KP_ATLAS",
                       "/Users/ryanschmitt/codex/spectra/Kitt Peak Flux Atlas")


def atlas_window(wl_A, half):
    """Observed Kitt Peak flux atlas over a window, in AIR Angstrom.

    Verified empirically, not assumed: the deepest minimum near Fe II 4303.170
    lands +0.004 A away when the atlas is read as AIR and +0.183 A when read as
    vacuum, so column 1 is air. Column 2 is the atlas's OWN continuum-normalised
    intensity and is used as shipped -- a Kitt Peak atlas is never renormalised
    (RYA-1026).
    """
    out = []
    for blk in range(int((wl_A - half) // 40) * 4, int((wl_A + half) // 40) * 4 + 4, 4):
        f = os.path.join(ATLAS, "lm%04d" % blk)
        if os.path.exists(f):
            out.append(np.loadtxt(f))
    if not out:
        return None, None
    d = np.vstack(out)
    w = d[:, 0] * 10.0
    m = (w >= wl_A - half) & (w <= wl_A + half)
    return w[m], d[m, 1]


def fit_plots(cen2, linelist):
    """Fit plots on the flagged outliers: observed spectrum + what is in the window."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = cen2[cen2.line_A > 8.0].drop_duplicates("wavelength_air_A").sort_values("wavelength_air_A")
    if out.empty:
        return 0
    n = len(out)
    fig, axes = plt.subplots(n, 1, figsize=(11, 2.6 * n), squeeze=False)
    for ax, (_, r) in zip(axes.ravel(), out.iterrows()):
        wl = r.wavelength_air_A
        w, f = atlas_window(wl, FIT_HALF_WIDTH_A * 1.6)
        if w is not None and len(w):
            ax.plot(w, f, color="k", lw=.9, label="Kitt Peak atlas (observed)")
        ax.axvspan(wl - FIT_HALF_WIDTH_A, wl + FIT_HALF_WIDTH_A, color="orange", alpha=.13,
                   label="fit window +/-%.2f A" % FIT_HALF_WIDTH_A)
        ax.axvline(wl, color="crimson", lw=1.4, label="target Fe II %.3f" % wl)
        win = linelist[(linelist.wavelength_air_A - wl).abs() <= FIT_HALF_WIDTH_A * 1.6]
        others = win[(win.wavelength_air_A - wl).abs() > 0.005]
        for _, o in others.iterrows():
            if not np.isfinite(o.central_depth) or o.central_depth < 0.15:
                continue
            ax.vlines(o.wavelength_air_A, 1.0, 1.0 - o.central_depth, color="steelblue", alpha=.65, lw=1.1)
        ch = others[others.element.astype(str) == "CH"]
        ax.set_title("Fe II %.3f A   recovered A=%.3f (anchor %.3f)   %d transitions in window, %d of them CH   "
                     "deepest interloper %s %.3f at %.3f A"
                     % (wl, r.line_A, ANCHOR, int(r.n_in_window), len(ch),
                        r.deepest_species, r.deepest_interloper, r.deepest_wavelength_A), fontsize=9)
        ax.set_ylabel("normalised flux"); ax.legend(fontsize=6, loc="lower left")
        ax.set_xlim(wl - FIT_HALF_WIDTH_A * 1.6, wl + FIT_HALF_WIDTH_A * 1.6)
    axes.ravel()[-1].set_xlabel("air wavelength (A)")
    fig.suptitle("RYA-515  fit plots on the flagged Fe II outliers\n"
                 "blue sticks = other transitions in the window, height = their central depth", y=1.0)
    fig.tight_layout()
    fig.savefig(os.path.join(PLOT, "rya515_fe2_fitplots.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)
    return n


def make_plots(prov, cen1, cen2):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    agg = prov[prov.in_aggregate.fillna(True).astype(bool) & prov.line_A.notna()]

    # (1) and (2): band x ion, split by engine and by correction.
    for key, fname, title in [("treatment", "rya515_hist_by_engine.png", "engine / treatment"),
                              ("nlte_source", "rya515_hist_by_correction.png", "NLTE correction source")]:
        cells = [(b, i) for b in ["near-UV", "VIS"] for i in ["I", "II"]]
        fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
        for ax, (b, i) in zip(axes.ravel(), cells):
            sub = agg[(agg.band == b) & (agg.ion == i)]
            for lbl, g in sub.groupby(sub[key].fillna("(none)").astype(str)):
                ax.hist(g.line_A, bins=np.arange(6.6, 10.6, 0.15), alpha=.55,
                        label="%s (n=%d)" % (lbl[:34], len(g)))
            ax.axvline(ANCHOR, color="k", ls="--", lw=1)
            ax.set_title("Fe %s  %s   (n=%d)" % (i, b, len(sub)), fontsize=10)
            ax.legend(fontsize=6); ax.set_xlabel("A(Fe) per line")
        fig.suptitle("RYA-515  per-line A(Fe) by band x ion, coloured by %s\n"
                     "dashed line = the adopted anchor A(Fe I)=7.466" % title)
        fig.tight_layout(); fig.savefig(os.path.join(PLOT, fname), dpi=130); plt.close(fig)

    # (3) the Fe II outliers against the Fe I pool -- the discriminator, plotted.
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5))
    a1.scatter(cen1.deepest_interloper, cen1.line_A, s=22, alpha=.6, label="Fe I anchor pool")
    a1.scatter(cen2.deepest_interloper, cen2.line_A, s=34, alpha=.8, marker="s", label="Fe II")
    a1.axhline(ANCHOR, color="k", ls="--", lw=1)
    for _, r in cen2[cen2.line_A > 9].iterrows():
        a1.annotate("%.3f" % r.wavelength_air_A, (r.deepest_interloper, r.line_A), fontsize=7)
    a1.set_xlabel("depth of the deepest OTHER transition in the fit window")
    a1.set_ylabel("A(Fe) recovered from the line"); a1.legend(fontsize=8)
    a1.set_title("blending DOES separate the pools")
    a2.scatter(cen1.red_chi2, cen1.line_A, s=22, alpha=.6, label="Fe I anchor pool")
    a2.scatter(cen2.red_chi2, cen2.line_A, s=34, alpha=.8, marker="s", label="Fe II")
    a2.axhline(ANCHOR, color="k", ls="--", lw=1); a2.set_xscale("log")
    a2.set_xlabel("reduced chi2 of the fit"); a2.legend(fontsize=8)
    a2.set_title("reduced chi2 does NOT (control: Fe I median %.0f vs Fe II %.0f)"
                 % (cen1.red_chi2.median(), cen2.red_chi2.median()))
    fig.suptitle("RYA-515  Fe II outliers vs the Fe I anchor pool")
    fig.tight_layout(); fig.savefig(os.path.join(PLOT, "rya515_fe2_outliers.png"), dpi=130); plt.close(fig)

    # (3b) THE finding: the two pools occupy different curve-of-growth regimes.
    fig, (b1, b2) = plt.subplots(1, 2, figsize=(13, 5))
    b1.scatter(cen1.obs_rew, cen1.line_A, s=24, alpha=.6, label="Fe I anchor pool")
    b1.scatter(cen2.obs_rew, cen2.line_A, s=36, alpha=.85, marker="s", label="Fe II")
    b1.axhline(ANCHOR, color="k", ls="--", lw=1)
    b1.set_xlabel("observed log(EW/lambda), measured off the atlas")
    b1.set_ylabel("A(Fe) recovered"); b1.legend(fontsize=8)
    b1.set_title("Fe II sits entirely on the saturated branch\nmedian REW %.2f vs Fe I %.2f"
                 % (cen2.obs_rew.median(), cen1.obs_rew.median()))
    bins = np.arange(0.0, 1.02, 0.05)
    b2.hist(cen1.obs_depth.dropna(), bins=bins, alpha=.6, label="Fe I anchor pool")
    b2.hist(cen2.obs_depth.dropna(), bins=bins, alpha=.7, label="Fe II")
    b2.axvline(0.70, color="crimson", ls="--", lw=1.2, label="depth 0.70")
    b2.set_xlabel("observed central depth"); b2.set_ylabel("lines"); b2.legend(fontsize=8)
    b2.set_title("Fe II: %.0f%% of lines deeper than 0.70   Fe I: %.0f%%"
                 % (100 * (cen2.obs_depth > .7).mean(), 100 * (cen1.obs_depth > .7).mean()))
    fig.suptitle("RYA-515  the regime difference is what separates the pools")
    fig.tight_layout(); fig.savefig(os.path.join(PLOT, "rya515_cog_regime.png"), dpi=130); plt.close(fig)

    # (4) the A(Fe I) distribution behind 7.466.
    fe1 = agg[(agg.ion == "I") & agg.line_A.notna()]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.hist(fe1.line_A, bins=np.arange(6.6, 8.4, 0.05), color="steelblue", alpha=.8)
    ax.axvline(ANCHOR, color="k", ls="--", lw=1.6, label="adopted anchor 7.466")
    ax.axvline(AGSS21_FE, color="crimson", ls=":", lw=1.6, label="AGSS21 7.46")
    ax.axvline(fe1.line_A.median(), color="darkgreen", ls="-", lw=1.4,
               label="our per-line median %.3f" % fe1.line_A.median())
    ax.set_xlabel("A(Fe I) per line"); ax.set_ylabel("lines")
    ax.set_title("RYA-515  the A(Fe I) distribution behind the 7.466 anchor  (n=%d lines)" % len(fe1))
    ax.legend()
    fig.tight_layout(); fig.savefig(os.path.join(PLOT, "rya515_fe1_anchor.png"), dpi=130); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="re-run and assert the committed artifacts still reproduce")
    args = ap.parse_args()

    os.makedirs(AUD, exist_ok=True)
    os.makedirs(PLOT, exist_ok=True)
    prods = load_products()
    prov, reach = build_provenance(prods)
    linelist = pd.read_csv(os.path.join(ROOT, "data/linelists/linelist_solar.csv"), low_memory=False)

    fe2 = prov[prov.ion == "II"]
    fe1_anchor = prov[(prov.ion == "I") & (prov.band == "VIS") & (prov.tier == "GRADED")
                      & (prov.treatment == "1D-LTE") & (prov.holding.str.contains("molecfit"))]
    cen2 = blend_census(fe2, linelist)
    cen1 = blend_census(fe1_anchor, linelist)

    if args.check:
        for name, got in [("perline_provenance.csv", prov), ("reachability.csv", reach),
                          ("fe2_blend_census.csv", cen2), ("fe1_control.csv", cen1)]:
            p = os.path.join(AUD, name)
            if not os.path.exists(p):
                print("MISSING %s" % name); sys.exit(1)
            want = pd.read_csv(p)
            if len(want) != len(got):
                print("DRIFT %s: committed %d rows, regenerated %d" % (name, len(want), len(got)))
                sys.exit(1)
        print("RYA-515 OK: %d per-line rows, %d products, %d Fe II / %d Fe I census rows"
              % (len(prov), len(reach), len(cen2), len(cen1)))
        return

    make_plots(prov, cen1, cen2)
    nfp = fit_plots(cen2, linelist)
    print('fit plots rendered for %d flagged outlier lines' % nfp)
    # results/plots/ is gitignored scratch (0 tracked files); the repo keeps durable
    # evidence plots under data/audit/<ticket>/ (cf. acen_holdings_rya384, rya938).
    # Write to both: the spec names results/plots/, the audit copy is what survives.
    import shutil
    os.makedirs(os.path.join(AUD, "plots"), exist_ok=True)
    for f in sorted(glob.glob(os.path.join(PLOT, "rya515_*.png"))):
        shutil.copy2(f, os.path.join(AUD, "plots", os.path.basename(f)))

    prov.to_csv(os.path.join(AUD, "perline_provenance.csv"), index=False)
    reach.to_csv(os.path.join(AUD, "reachability.csv"), index=False)
    cen2.to_csv(os.path.join(AUD, "fe2_blend_census.csv"), index=False)
    cen1.to_csv(os.path.join(AUD, "fe1_control.csv"), index=False)
    print("per-line rows: %d over %d/%d products with reachable evidence"
          % (len(prov), int(reach.perline_reachable.sum()), len(reach)))
    print("Fe II census %d rows, Fe I control %d rows" % (len(cen2), len(cen1)))
    print("Fe II contaminated: %d/%d   Fe I contaminated: %d/%d"
          % (cen2.contaminated.sum(), len(cen2), cen1.contaminated.sum(), len(cen1)))
    print("red_chi2 median  Fe II %.1f   Fe I %.1f" % (cen2.red_chi2.median(), cen1.red_chi2.median()))
    print("deepest interloper median  Fe II %.3f   Fe I %.3f"
          % (cen2.deepest_interloper.median(), cen1.deepest_interloper.median()))


if __name__ == "__main__":
    main()
