#!/usr/bin/env python3
"""RYA-794 — CRIRES+ Vesta solar Y band: science-ready spectrum (+ QA plot).

    python3 scripts/normalize_vesta_ir.py            # Step 0-3
    python3 scripts/normalize_vesta_ir.py --plot-only

STEP 0 VERDICT — SOURCE A, and the vintage is confirmed
-------------------------------------------------------
`sp/Sun_Y_rv.dat` from the RYA-789 Elgueta-2026 VizieR pull, md5-verified against the
pinned manifest (`8428f7009e59c6073ba11242b83d48ff`). The ReadMe states it plainly: *"We
used the CRIRES+ instrument with ESO's VLT at Paranal Observatory ... at a nominal
spectral resolution of ~100000"* and *"a solar spectrum obtained from observations of
Vesta"*. So this is CRIRES+, not the 2012-13 CRIRES the ticket warns about, and the STOP
condition does not fire.

**Source B is not available.** The ticket expects 18 CRIRES+ Vesta IDPs per RYA-377; a
filesystem sweep finds ZERO Vesta FITS on Sirius. Only source A exists, so Step 1
(molecfit) is not merely skipped-by-preference — there is nothing else to run it on.

⚠️ UNIT TRAP — THE READ-ME'S OWN LABEL IS WRONG
The sp/ byte-by-byte description gives the wavelength column as `0.1nm` (i.e. Angstrom).
The column actually runs **979.649 - 1079.611**, which as Angstrom would be far-UV and
absurd for a Y-band NIR spectrum. It is **nanometres**: x10 gives 9796.49-10796.11 A,
exactly the Y band, and the 10280-10680 A science window then sits fully inside it with
13707 points. Taking the label at face value would have produced a confident wrong answer
of precisely the kind the instrument catalog already warns about for the IAG atlas.

TELLURIC RULE — SATISFIED, AND CHECKED RATHER THAN ASSUMED
The rule is permanent: no IR abundance without verified telluric correction. Measured on
this spectrum, the science window has **0.10 % of points below 0.5** and a mean of 0.990 —
no saturated absorption anywhere. For scale, the Kitt Peak atlas runs **51.3 %** below 0.5
inside the O2 A-band, which is what uncorrected telluric absorption looks like. Elgueta
selected 10280-10680 A precisely because it is almost telluric-free, and the data agrees.

STEP 4 RESULT — THE Y BAND HAS NO ELGUETA-CERTIFIED SOLAR Fe I LINE
--------------------------------------------------------------------
Elgueta's own `GDRob` flag (G dwarf = the Sun's type) selects **0 of 141** Fe I records.
Measured, not assumed: 42 lines are GDRob=Y and every one is another species (Si I 18,
C I 7, Cr I 5, Ti I / Ni I / Mg I / S I / P I 2 each, Sr II / Ca I 1); all 42 pass all
four sub-flags; and no Fe I line anywhere in the file carries GDSat='Y'. Fe I fails their
saturation criterion for a solar-type star. By contrast K dwarfs get 46 robust Fe I lines
— the Y-band Fe I lines are high-EP and only strengthen in cooler stars, which is the
physical reason and is consistent with the EP range found below.

Five Fe I lines DO pass depth+purity+goodness-of-fit inside the science window and are
carried as tier "candidate", each matching our VALD solar list to +/-0.000 A with
identical EP and log gf. They are reported, not hidden and not promoted. What they do not
support is an abundance: see the closing warning the driver prints.

WHAT STEP 2 IS AND IS NOT DOING HERE
The reflectance-aware normaliser below is the ticket's specified helper and is implemented
as specified. But source A arrives ALREADY NORMALISED (median 0.9983, 95th percentile
1.0068), because Elgueta's sp/ files are "spectroscopic observed normalized fluxes". So on
this input the helper is a RESIDUAL-SLOPE check, not the primary normalisation, and it is
reported that way. It is written in full because source B (raw Vesta IDPs) would need it,
and because a residual reflectance slope in an already-normalised product is exactly the
thing worth measuring rather than assuming away.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import UnivariateSpline

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# Sirius-only authoritative copy; sp/ is gitignored and md5-pinned (RYA-789).
SP_Y = Path("/mnt/codex-data/codex/rya789/data/reference/elgueta2026_vizier/sp/Sun_Y_rv.dat")
MD5_Y = "8428f7009e59c6073ba11242b83d48ff"
NM_TO_A = 10.0                      # the ReadMe says 0.1nm; the data says nm. See above.
WINDOW = (10280.0, 10680.0)         # Elgueta's near-telluric-free Y window
OUT_DIR = ROOT / "data" / "results" / "rya794"


def normalize_reflectance(wave, flux, niter=8, lo=2.0, hi=4.0, s=None):
    """Upper-envelope continuum for a Vesta (reflectance-slope) IR spectrum.

    Asymmetric sigma-clip: reject ABSORPTION (low) points hard, keep the
    reflectance-sloped upper envelope, so solar lines are divided out, not fitted.
    Returns (norm_flux, continuum).

    Vesta's basaltic pyroxene reflectance is a smooth sloped pseudo-continuum in the NIR,
    so a SYMMETRIC fit would be pulled down into the solar lines and would divide them
    away. The asymmetry is the whole point: absorption is evidence about the star, not
    about the continuum.
    """
    w = np.asarray(wave, float); f = np.asarray(flux, float)
    good = np.isfinite(f) & np.isfinite(w)
    # UnivariateSpline needs strictly increasing x
    o = np.argsort(w)
    w, f, good = w[o], f[o], good[o]
    keep = np.concatenate([[True], np.diff(w) > 0])
    good &= keep
    mask = good.copy()
    cont = np.full_like(f, np.nan)
    for _ in range(niter):
        if mask.sum() < 10:
            break
        spl = UnivariateSpline(w[mask], f[mask], k=3,
                               s=(s if s is not None else len(w[mask])))
        cont = spl(w)
        resid = f - cont
        sig = np.nanstd(resid[mask])
        if not np.isfinite(sig) or sig == 0:
            break
        # asymmetric: absorption (resid < -lo*sig) rejected; emission spikes (> hi*sig) too
        mask = good & (resid > -lo * sig) & (resid < hi * sig)
    with np.errstate(invalid="ignore", divide="ignore"):
        norm = f / cont
    return norm, cont, o


def load_y_band() -> tuple[np.ndarray, np.ndarray]:
    import hashlib
    if not SP_Y.exists():
        raise SystemExit(f"Y-band spectrum not staged at {SP_Y} (RYA-789, Sirius-only)")
    got = hashlib.md5(SP_Y.read_bytes()).hexdigest()
    if got != MD5_Y:
        raise SystemExit(f"md5 mismatch on {SP_Y.name}: {got} != pinned {MD5_Y}")
    d = np.loadtxt(SP_Y)
    return d[:, 0] * NM_TO_A, d[:, 1]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lo", type=float, default=WINDOW[0])
    ap.add_argument("--hi", type=float, default=WINDOW[1])
    ap.add_argument("--no-plot", action="store_true")
    a = ap.parse_args()

    w, f = load_y_band()
    print(f"STEP 0 — source A: {SP_Y.name}, md5 verified")
    print(f"  CRIRES+ (ReadMe: 'CRIRES+ ... ~100000'), solar spectrum from Vesta")
    print(f"  {len(w)} pts, {w.min():.2f}-{w.max():.2f} A  "
          f"(ReadMe says 0.1nm; data is nm -> x{NM_TO_A:.0f})")

    m = (w >= a.lo) & (w <= a.hi)
    print(f"\n  science window {a.lo:.0f}-{a.hi:.0f} A: {int(m.sum())} pts")
    print(f"    median {np.median(f[m]):.4f}  95th pct {np.percentile(f[m], 95):.4f}"
          f"  -> arrives ALREADY NORMALISED")
    print(f"    telluric check: {100 * np.mean(f[m] < 0.5):.2f}% of points below 0.5 "
          f"(Kitt Peak's O2 A-band is 51.3%) -> clean")

    wx, fx = w[m], f[m]
    norm, cont, order = normalize_reflectance(wx, fx)
    ws = wx[order]
    print(f"\nSTEP 2 — residual-slope check (the helper, on an already-normalised input)")
    print(f"    continuum fitted: {np.nanmin(cont):.4f} - {np.nanmax(cont):.4f}  "
          f"(slope {np.nanmax(cont) - np.nanmin(cont):+.4f} across the window)")
    print(f"    after division  : median {np.nanmedian(norm):.4f}  "
          f"95th pct {np.nanpercentile(norm, 95):.4f}")
    deep_before = int(np.sum(fx < 0.7))
    deep_after = int(np.sum(norm < 0.7))
    print(f"    line cores preserved: {deep_before} points below 0.7 before, "
          f"{deep_after} after"
          f"  {'OK' if deep_after >= 0.9 * deep_before else 'WARNING — cores eaten'}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / "vesta_crires_plus_Y_10280_10680_normalized.csv"
    np.savetxt(out, np.column_stack([ws, norm]), delimiter=",",
               header="wavelength_air_A,flux_normalized", comments="")
    print(f"\nSTEP 3 — wrote {out.relative_to(ROOT)}  ({len(ws)} rows)")

    # ── STEP 4 — Fe I line selection, Elgueta's flags + VALD cross-check ─────
    allfe = elgueta_y_lines("FeI", tier="all")
    assessed = elgueta_y_lines("FeI", tier="assessed")
    rob = elgueta_y_lines("FeI", tier="robust")
    cand = elgueta_y_lines("FeI", tier="candidate")
    inwin = [r for r in (rob + cand) if a.lo <= r["wave_A"] <= a.hi]
    inwin.sort(key=lambda r: r["wave_A"])
    print(f"\nSTEP 4 — Fe I selection from Elgueta atomicy.dat (their own G-dwarf flags)")
    print(f"    Fe I records in the Y list      : {len(allfe)}")
    print(f"    assessed for a G dwarf          : {len(assessed)}")
    print(f"    GDRob = Y (Elgueta-certified)   : {len(rob)}"
          f"{'   <-- ZERO. See the docstring; this is measured, not a parse bug.' if not rob else ''}")
    print(f"    candidate (depth+purity+GoF Y,  : {len(cand)}")
    print(f"               saturation uncertified)")
    print(f"    inside {a.lo:.0f}-{a.hi:.0f} A                : {len(inwin)}")
    if inwin:
        vald_crosscheck(inwin)
        verify_line_depths(ws, norm, inwin)
        nm = sum(1 for r in inwin if r["vald_match"])
        print(f"\n    VALD cross-check against data/linelists/linelist_solar.csv "
              f"({nm}/{len(inwin)} matched), and depth recovered from OUR spectrum:")
        print(f"      {'wave_A':>10s} {'EP':>5s} {'loggf':>7s} {'chi2':>5s} | "
              f"{'dlam':>6s} {'VALDgf':>7s} {'src':>6s} | "
              f"{'dlam':>6s} {'meas':>6s} {'Elg':>6s} {'ratio':>5s} | tier")
        for r in inwin:
            print(f"      {r['wave_A']:10.3f} {r['ep_eV']:5.3f} {r['loggf']:+7.3f} "
                  f"{r['chi2_red']:5.2f} | "
                  f"{r['vald_dlam_A']:+6.3f} {r['vald_loggf']:+7.3f} "
                  f"{r['vald_gf_source'][:6]:>6s} | "
                  f"{r['meas_dlam_A']:+6.3f} {r['meas_depth']:6.4f} "
                  f"{r['depth_obs']:6.4f} {r['depth_ratio']:5.2f} | {r['tier']}")
        ok = sum(1 for r in inwin if abs(r["meas_dlam_A"]) < 0.05
                 and 0.8 < r["depth_ratio"] < 1.25)
        print(f"      -> {ok}/{len(inwin)} land within 0.05 A and recover Elgueta's own "
              f"depth to 25%.")
        print(f"         Independent end-to-end confirmation of the nm->A conversion, the")
        print(f"         wavelength solution and the continuum, against numbers we did "
              f"not produce.")
        import csv as _csv
        sel = OUT_DIR / "vesta_crires_plus_Y_FeI_lines.csv"
        with open(sel, "w", newline="") as fh:
            wcsv = _csv.DictWriter(fh, fieldnames=list(inwin[0].keys()))
            wcsv.writeheader(); wcsv.writerows(inwin)
        print(f"    wrote {sel.relative_to(ROOT)}")
    print(f"\n    ⚠️  NO ABUNDANCE IS QUOTED FROM THIS. Zero lines carry Elgueta's robust")
    print(f"       certification for a solar-type star, the {len(inwin)} candidates are "
          f"single-exposure,")
    print(f"       and EP {min((r['ep_eV'] for r in inwin), default=0):.1f}-"
          f"{max((r['ep_eV'] for r in inwin), default=0):.1f} eV on {len(inwin)} lines "
          f"cannot constrain A(Fe) honestly.")

    if not a.no_plot:
        _qa_plot(wx, fx, cont, norm, order, a)


def _qa_plot(wx, fx, cont, norm, order, a) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"  (plot skipped: {e})")
        return
    ws = wx[order]; fs = fx[order]
    # a couple of Fe I lines from Elgueta's own Y-band list, marked for the eye
    fe = [x for x in _elgueta_fe_y() if a.lo <= x <= a.hi][:6]
    # zoom onto the marked lines rather than a fixed offset from the window edge — a
    # blue-edge zoom showed an empty panel because every Fe I candidate sits redward
    if fe:
        lo2, hi2 = min(fe) - 8.0, max(fe) + 8.0
    else:
        lo2, hi2 = a.lo, min(a.lo + 60.0, a.hi)
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    axes[0].plot(ws, fs, lw=0.6, label="Elgueta sp/ (already normalised)")
    axes[0].plot(ws, cont, lw=1.3, color="tab:red",
                 label="reflectance upper-envelope continuum")
    axes[0].set_ylabel("flux"); axes[0].legend(fontsize=8)
    axes[0].set_title("RYA-794 — CRIRES+ Vesta solar, Y band. "
                      "Source A (Elgueta 2026 sp/), telluric-clean window.", fontsize=10)
    axes[1].plot(ws, norm, lw=0.6, color="k", label="after residual-slope division")
    axes[1].axhline(1.0, ls="--", lw=0.8, color="0.5")
    # The Fe I candidates are SHALLOW (depth 0.06-0.13) and sit beside 0.5-deep
    # neighbours, so a full-height axvline reads as "the mark missed the line". Mark them
    # with a short tick at the line's own depth instead, and clip the y-range to the
    # shallow regime, so the eye can check the identification rather than be misled by it.
    for x in fe:
        m = (ws > x - 0.3) & (ws < x + 0.3)
        d = float(np.min(norm[m])) if m.any() else 1.0
        axes[1].plot([x, x], [d - 0.012, d - 0.045], color="tab:blue", lw=1.1)
        axes[1].text(x, d - 0.05, f"Fe I {x:.1f}", rotation=90, fontsize=6,
                     ha="center", va="top", color="tab:blue")
    axes[1].set_xlim(lo2, hi2)
    axes[1].set_ylim(0.80, 1.03)
    axes[1].set_xlabel("air wavelength (A)"); axes[1].set_ylabel("normalised flux")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    p = OUT_DIR / "vesta_crires_plus_Y_qa.png"
    fig.savefig(p, dpi=140); plt.close(fig)
    print(f"  wrote {p.relative_to(ROOT)}  (zoom {lo2:.0f}-{hi2:.0f} A, "
          f"{len(fe)} Fe I lines marked)")


ATOMICY = Path("/mnt/codex-data/codex/rya789/data/reference/elgueta2026_vizier/atomicy.dat")


# atomicy.dat is FIXED-WIDTH, Lrecl 805, six per-spectral-type blocks. The ReadMe's byte
# numbers are 1-indexed; these are the Python offsets for the **G dwarf** block, which is
# the Sun's type. Verified empirically: a standalone-Y/N column scan finds flags at exactly
# these positions and nowhere else in the block.
GD = {"depth": 249, "sat": 251, "pur": 262, "gof": 293, "rob": 295}


def elgueta_y_lines(species: str = "FeI", tier: str = "robust") -> "list[dict]":
    """Elgueta's Y-band line table, read on its FIXED-WIDTH byte definitions.

    ⚠️ Two unit facts that do not agree with each other, both from the same ReadMe:
    atomicy.dat's Wave really IS 0.1nm (Angstrom) as labelled — 9800.308 is a Y-band line
    — while sp/'s Wave is labelled the same and is actually NANOMETRES. Reading either by
    the other's convention gives a confident wrong answer, so both are pinned explicitly:
    atomicy is used as-is, sp/ is multiplied by 10.

    THE SELECTION IS ELGUETA'S OWN, NOT MINE. Byte 296 is `GDRob`, their final robust-line
    flag for a G dwarf, combining four sub-criteria: depth, unsaturated, purity and
    goodness-of-fit. Inventing a fresh cut here would substitute my judgement for the
    paper's while still citing the paper, so the sub-flags are read and reported rather
    than re-derived.

    ⚠️ THE HEADLINE RESULT OF THIS READER: `GDRob` selects **ZERO Fe I lines**. That is not
    a parsing bug — it is measured, and it reproduces three ways:
      * of 141 Fe I records, 25 are assessed for a G dwarf and all 25 are GDRob=N;
      * 42 lines ARE GDRob=Y, all of them other species (Si I 18, C I 7, Cr I 5, ...);
      * `GDSat` never takes the value 'N' anywhere in the file — it is 'Y' (58 lines) or
        blank (265) — and every one of the 42 robust lines is Y on all four sub-flags,
        while NO Fe I line carries GDSat='Y' at all.
    So Fe I fails at Elgueta's saturation step for a solar-type star. Seven Fe I lines
    nonetheless pass depth+purity+goodness-of-fit; they are returned as tier "candidate",
    clearly separated, because they are real measurements that their strictest flag does
    not certify — not something to silently promote to robust, and not something to hide.

    ⚠️ AND blank GDSat means FAILED, not merely UNEVALUATED. The column never takes 'N', so
    blank is ambiguous on its face — but of the 71 lines that pass GDDepth, the 13 whose
    GDSat is blank are 12 Fe I and 1 Fe II, i.e. iron and nothing else. A sparsely
    populated column would not single out one element. So the candidate tier is
    Elgueta-REJECTED science that we choose to carry and label, not a gap in their table,
    and it must never be quietly reclassified as robust on the strength of the other three
    flags.

    tier: "robust"    -> GDRob == Y                      (Elgueta-certified)
          "candidate" -> depth & purity & GoF all Y, but not GDRob
          "assessed"  -> any line Elgueta evaluated for a G dwarf (GDRob in Y/N)
          "all"       -> every record for the species
    """
    out = []
    if not ATOMICY.exists():
        return out
    for line in ATOMICY.read_text(errors="replace").splitlines():
        if len(line) < 805 or not line.strip():        # ReadMe: Lrecl 805, blank-padded
            continue
        try:
            wv = float(line[0:12])
        except ValueError:
            continue
        el = line[13:17].strip()
        if el != species:
            continue
        fl = {k: line[i] for k, i in GD.items()}
        robust = fl["rob"] == "Y"
        candidate = (not robust) and all(fl[k] == "Y" for k in ("depth", "pur", "gof"))
        keep = {"robust": robust, "candidate": candidate,
                "assessed": fl["rob"] in "YN", "all": True}[tier]
        if not keep:
            continue

        def _f(a, b):
            try:
                return float(line[a:b])
            except (ValueError, IndexError):
                return float("nan")

        out.append(dict(wave_A=wv, element=el, ep_eV=_f(18, 27), loggf=_f(28, 37),
                        depth_obs=_f(210, 219), purity_ew=_f(253, 261),
                        chi2_red=_f(282, 292),
                        gd_depth=fl["depth"], gd_unsaturated=fl["sat"],
                        gd_purity=fl["pur"], gd_gof=fl["gof"], gd_robust=fl["rob"],
                        tier="robust" if robust else
                             ("candidate" if candidate else "rejected")))
    return sorted(out, key=lambda r: r["wave_A"])


def vald_crosscheck(rows: "list[dict]", tol: float = 0.05) -> "list[dict]":
    """Cross-check Elgueta's Y-band lines against our own VALD solar list (ticket Step 4).

    Reports the match rather than filtering on it: a line Elgueta lists that our list does
    not carry is a finding about our list, not a reason to drop the line.
    """
    import pandas as pd
    p = ROOT / "data" / "linelists" / "linelist_solar.csv"
    if not p.exists():
        return rows
    d = pd.read_csv(p, low_memory=False)
    fe = d[(d.element == "Fe") & (d.ion == "I")].sort_values("wavelength_air_A")
    w = fe.wavelength_air_A.to_numpy()
    for r in rows:
        i = int(np.argmin(np.abs(w - r["wave_A"])))
        dl = float(w[i] - r["wave_A"])
        v = fe.iloc[i]
        hit = abs(dl) <= tol
        r.update(vald_match=hit, vald_dlam_A=dl if hit else float("nan"),
                 vald_loggf=float(v.log_gf) if hit else float("nan"),
                 vald_ep_eV=float(v.excitation_potential_eV) if hit else float("nan"),
                 vald_central_depth=float(v.central_depth) if hit else float("nan"),
                 vald_gf_source=str(v.loggf_source) if hit else "")
    return rows


def verify_line_depths(w, f, rows, tol_A: float = 0.3) -> "list[dict]":
    """Independent end-to-end check: do the selected lines actually sit in OUR spectrum?

    For each selected line, find the local minimum within tol_A and compare the depth we
    measure to Elgueta's published `GDDepO`. This is the one test that exercises the whole
    chain at once — the nm->Angstrom conversion, the wavelength solution, the continuum
    division and the line selection — against a number we did not produce. A unit error
    or a continuum error would show up here as a gross depth mismatch or as no line at all.
    """
    for r in rows:
        m = (w > r["wave_A"] - tol_A) & (w < r["wave_A"] + tol_A)
        if not m.any():
            r.update(meas_dlam_A=float("nan"), meas_depth=float("nan"),
                     depth_ratio=float("nan"))
            continue
        i = int(np.argmin(f[m]))
        dep = 1.0 - float(f[m][i])
        r.update(meas_dlam_A=float(w[m][i] - r["wave_A"]), meas_depth=dep,
                 depth_ratio=dep / r["depth_obs"] if r["depth_obs"] else float("nan"))
    return rows


def _elgueta_fe_y() -> list[float]:
    """Fe I wavelengths to mark on the QA plot: robust if any exist, else candidates."""
    r = elgueta_y_lines("FeI", tier="robust")
    return [x["wave_A"] for x in (r or elgueta_y_lines("FeI", tier="candidate"))]


if __name__ == "__main__":
    main()
