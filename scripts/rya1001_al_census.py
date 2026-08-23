#!/usr/bin/env python3
"""
scripts/rya1001_al_census.py — RYA-1001 Al PHASE 0: full-band line census
=========================================================================
CENSUS ONLY. No synthesis, no abundance, no gf adoption. This inventories what
gradeable Al I/II/III lines exist per band, what the best gf source for each is,
whether HFS is carried, which instrument each line has a home on, and whether a
telluric correction stands between us and it.

WHY A SCRIPT AND NOT A TABLE
----------------------------
Same reason as RYA-709: the answer changes when a holding is registered, a line
list is re-pulled, or a lab paper lands. A number that was true once and is quoted
forever is how RYA-707 published "Al 7835/8772: NO DATA".

THE FIREWALL (RYA-161)
----------------------
Nothing here adopts a gf. That a source moves A(Al) toward 6.43 is NOT a reason to
prefer it. Every ranking below is on PROVENANCE ONLY: experiment > independent
theory > compilation grade > semi-empirical. The census reports the ladder; the
substitution, if any, is a pool rebuild owned by a later ticket.

THE gf LADDER, declared in advance
----------------------------------
  1 EXP-BURHEIM23   Burheim, Hartman & Nilsson 2023, A&A 672 A197 (arXiv:2309.06273)
                    — experimental branching fractions (FTS + hollow cathode lamp)
                    x radiative lifetimes. sigma = the paper's own per-line Unc_gf.
  2 THEORY-P19      Papoulia, Ekman & Jonsson 2019, A&A 621 A16 — multiconfiguration
                    (MCDHF/RCI). INDEPENDENT of NIST. Quoted where Burheim tabulates it.
  3 NIST-<grade>    NIST ASD compilation grade (Kelleher & Podobedova 2008 for Al).
                    sigma from the ASD accuracy ladder.
  4 THEORY-OP95     Mendoza, Eissner, Le Dourneuf & Zeippen 1995, J.Phys.B 28 3485 —
                    Opacity Project close-coupling. This is what the GES/canonical
                    tag `1995JPhB..` actually is (see OQ2). THEORY, not a lab value.
  5 SEMIEMP-K75     Kurucz 1975 semi-empirical (GES tag `K75`). RYA-161 territory,
                    0.1-0.3 dex.
  6 UNGRADED

Nordlander & Lind 2017 is deliberately ABSENT from the ladder: its Al uncertainties
trace to NIST/Kelleher-Podobedova 2008, so counting it as an independent source would
double-count the NIST rung (RYA-835 established this; re-stated here so the ladder
cannot be mis-read).

Linear issue: RYA-1001
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.wavelength_util import air_to_vac, vac_to_air      # noqa: E402  SSOT
from pipeline import band_policy                                  # noqa: E402
from pipeline.coverage import coverage_at, CoverageError          # noqa: E402
from pipeline.telluric_policy import gate_holding                 # noqa: E402

OUT = ROOT / "data" / "results" / "rya1001"

# ─────────────────────────────────────────────────────────────────────────────
# BURHEIM 2023 — Table 3 (`tab:loggf_comp`), the 12 lines with a DERIVED log gf.
#
# TRANSCRIBED FROM THE arXiv LaTeX SOURCE (2309.06273, file 45394corr.tex), not from
# a rendered PDF and not from the abstract's stated 670-4200 nm RANGE. A range is not
# a line list.
#
# ⚠️ THE COLUMN TRAP THAT DEFEATED RYA-835. Tables 1-3 carry TWO numeric columns of
# similar magnitude: `sigma [cm^-1]` (wavenumber) and `lambda_vac [A]`. RYA-835's
# BURHEIM_LINES tuple is the UNION of both columns — 35 numbers, wavenumbers and
# wavelengths mixed — so every "nearest Burheim line" it computed was meaningless.
# That is what produced `burheim_nearest_A = 7836.521` (a WAVENUMBER) 0.387 A from
# our 7836.134 A line, and `6697.864` (a VACUUM WAVELENGTH) "1.679 A" from 6696.185.
# Here only lambda_vac is used, and it is converted to air before any match.
#
# Unc_gf (%) is from Table 1 (`tab:BFtable`), which is the table that carries them.
# ─────────────────────────────────────────────────────────────────────────────
BURHEIM_2023 = pd.DataFrame([
    # lam_vac_A, loggf_exp, unc_pct, loggf_P19, loggf_K95, loggf_M00, transition
    (6697.864,  -1.46,    9.0,  -1.499, -1.347, -1.569, "4s 2S1/2 - 5p 2P3/2"),
    (6700.522,  -1.76,    8.0,  -1.808, -1.647, -1.870, "4s 2S1/2 - 5p 2P1/2"),
    (11256.270,  0.167,   5.0,   0.206,  0.276,  0.170, "3d 2D3/2 - 4f 2F5/2"),
    (11258.008,  0.327,   2.0,   0.362,  0.431,  0.325, "3d 2D5/2 - 4f 2F7/2"),
    (12753.397, -2.29,    9.0,  -2.348, -2.217, -2.257, "3d 2D5/2 - 5p 2P3/2"),
    (12760.765, -2.62,   11.0,  -2.588, -2.472, -2.513, "3d 2D3/2 - 5p 2P1/2"),
    (13127.005,  0.232,   1.5,   0.215,  0.270,  0.219, "4s 2S1/2 - 4p 2P3/2"),
    (13154.350, -0.0980,  3.1,  -0.0867, -0.030, -0.083, "4s 2S1/2 - 4p 2P1/2"),
    (38632.83,   0.259,   7.0,   0.361,  0.358,  0.365, "5s 2S1/2 - 5p 2P3/2"),
    (38721.19,  -0.0396,  6.0,   0.0596, 0.056,  0.064, "5s 2S1/2 - 5p 2P1/2"),
    (41841.42,   0.607,  10.0,   0.513,  0.578,  0.500, "nd 2D3/2 - 4f 2F5/2"),
    (41920.7,    0.737,  11.0,   0.667,  0.732,  0.654, "nd 2D5/2 - 4f 2F7/2"),
], columns=["lam_vac_A", "loggf_exp", "unc_pct", "loggf_P19", "loggf_K95",
            "loggf_M00", "transition"])

#: Burheim Table 2 (`tab:BRtable`) tabulates BRANCHING FRACTIONS ONLY for these
#: lam_vac. They carry NO derived log gf, so Burheim does NOT grade them. Kept
#: explicitly so "Burheim covers 10875.953 / 16723.541 / ..." (the RYA-835 docstring
#: claim) can be shown to be a Table-2 hit, i.e. not a gf.
BURHEIM_BF_ONLY_VAC = (21169.58, 21098.84, 16767.948, 16723.541, 53575.35, 53406.34,
                       10894.716, 10875.953, 51128.64, 50974.74, 10771.313,
                       51022.78, 10785.000, 25029.69, 24992.75, 8843.705)

#: Match tolerance, line-centre. Generous enough to absorb the air<->vac step and
#: NIST-vs-VALD centre differences (the largest real offset seen is 0.036 A), tight
#: enough to exclude the 6696.023-vs-6696.185 pair (0.162 A) which are DIFFERENT
#: transitions. Cross-checked by the EP test below, which is the real discriminant.
MATCH_TOL_A = 0.06
#: Two lines are the same transition only if their lower level agrees to this.
#: RYA-835 matched 6696.185 (EP 4.0215) to NIST 6696.015 (EP 3.1427) on wavelength
#: alone: a 0.88 eV mismatch, i.e. a different lower level entirely.
MATCH_TOL_EP_EV = 0.02

#: NIST ASD accuracy ladder -> worst-case % on A_ki (enumerated in FULL including the
#: '+' tiers; omitting them once put B+ (<=7%) below B (<=10%), RYA-592).
NIST_ACC_PCT = {"AAA": 0.3, "AA": 1.0, "A+": 2.0, "A": 3.0, "B+": 7.0, "B": 10.0,
                "C+": 18.0, "C": 25.0, "D+": 40.0, "D": 50.0, "E": 100.0}

#: RYA-709's declared usable-depth triage window. Re-used, not re-invented.
DEPTH_LO, DEPTH_HI = 0.05, 0.60
#: RYA-709's HFS/feature grouping distance, kept as the DEFAULT separation for
#: unrelated lines.
GROUP_A = 0.05
#: Maximum wavelength span of ONE transition's hyperfine multiplet. 27Al (I=5/2, 100%
#: abundant) splits the 4s/3d lines over a few hundredths of an Angstrom; the widest
#: in our own list is 16750.4549-16750.6514 = 0.197 A. A distance-only rule at
#: GROUP_A=0.05 tears that into three "lines" and tears 13123.4 (internal gap 0.0577 A)
#: into two — which is what RYA-709's per_line.csv actually does. The multiplet is
#: instead held together by its PHYSICS key (same ion, same lower level, same predicted
#: central depth, which is a per-parent-line quantity), with this span as the outer bound.
#: MEASURED, not chosen. Hyperfine splitting is bounded in ENERGY, so the separation
#: test must be in wavenumber; a fixed Angstrom span is 54 km/s at the Al II 1670 A
#: lines and 5 km/s at 16750 A, and no single value can work at both ends (a 0.30 A
#: rule merged 438 FUV features down to 165, and still split 16750 in three).
#: In wavenumber the two populations separate with a factor >4 of clear air on each
#: side, measured on our own Al list:
#:     largest gap INSIDE a known multiplet   16750.45-16750.65  0.070 cm-1
#:                                            13123.37-13123.45  0.045 cm-1
#:                                            6696.00 -6696.02   0.042 cm-1
#:     smallest gap BETWEEN distinct lines    3092.710-3092.839  1.349 cm-1
#:                                            6696.185-6696.788  1.345 cm-1
#: 0.30 cm-1 sits in that empty band. It is a measured separation, not a tuned knob.
HFS_SPAN_CM1 = 0.30
#: Lower-level match tolerance for "same transition". The linelist rounds EP to 4 dp,
#: so 3d 2D5/2 appears as both 4.0216 and 4.0217; 2e-4 absorbs that without reaching
#: 3d 2D3/2 (4.0215 vs 4.0216 is separated by the 0.30 A span rule instead).
HFS_EP_TOL_EV = 2.0e-4

#: RYA-161 semi-empirical Kurucz systematic, carried when nothing better ties.
K75_SYSTEMATIC_DEX = 0.20
#: Opacity-Project close-coupling: no per-line published sigma. Carried as the same
#: order as the semi-empirical bar rather than blank — a hidden bar is a defect.
OP95_SYSTEMATIC_DEX = 0.20


def pct_to_dex(pct: float) -> float:
    """Fractional accuracy on gf -> dex. sigma_dex = log10(1 + p)."""
    return float(np.log10(1.0 + pct / 100.0))


# ─────────────────────────────────────────────────────────────────────────────
# 1. FEATURES — collapse HFS components into physical lines
# ─────────────────────────────────────────────────────────────────────────────
def features(ll: pd.DataFrame, element: str = "Al") -> pd.DataFrame:
    """One row per physical transition, HFS components SUMMED.

    ⚠️ RYA-709's `features()` reports `gf = max` over the cluster — the STRONGEST
    HFS COMPONENT, not the line. For Al that is not a rounding difference: the
    6696.02 multiplet's strongest component is -1.886 while the six components sum
    to -1.460, a 0.43 dex under-report on a line whose whole point is its gf.
    Al I is 100% 27Al with I=5/2, so every 4s/3d line is genuinely split; summing
    is the only correct collapse.

    Grouping: consecutive-in-wavelength rows within HFS_SPAN_A that share ion and
    lower level. Two things that look like they belong in the key and do NOT:

    * `GROUP_A = 0.05` (RYA-709's distance) is too tight, and ANY fixed-Angstrom
      distance is the wrong unit — see HFS_SPAN_CM1. At 0.05 A the 13123.4 multiplet
      (widest internal gap 0.0577 A) splits into a 3+3 pair and 16750 (span 0.197 A)
      into three. RYA-709's per_line.csv does exactly this today.
    * `central_depth` equality looks like a safe extra key and is WRONG. VALD gives
      true hyperfine components the parent's depth (so they match) but computes
      fine-structure partners separately (so they do not). Keying on it splits
      8773.896/8773.899 and the two 7836.134 entries — pairs 0.003 A and 0.000 A
      apart that no spectrograph can resolve and that NIST itself lists as one
      transition's components.

    POSITIVE CONTROL: `verify_against_canonical()` checks the summed gf and component
    count against `canonical_gf.csv`, which carries its own independent
    `hfs_n_components` and an already-summed `log_gf`. The rule is not accepted on
    the strength of the reasoning above; it is accepted because it reproduces that
    column.
    """
    a = ll[ll.element == element].sort_values("wavelength_air_A").copy()
    if a.empty:
        return a
    dsigma = (1e8 * a.wavelength_air_A.diff().abs()
              / (a.wavelength_air_A * a.wavelength_air_A.shift())).fillna(9e9)
    same = ((dsigma <= HFS_SPAN_CM1)
            & (a.excitation_potential_eV.diff().abs().fillna(9e9) <= HFS_EP_TOL_EV)
            & (a.ion == a.ion.shift()))
    a["_k"] = (~same).cumsum()
    g = a.groupby("_k")
    out = pd.DataFrame({
        "ion": g.ion.first(),
        # gf-weighted centroid: the wavelength the blended feature actually presents
        "wave_air_A": g.apply(lambda d: float(np.average(d.wavelength_air_A,
                                                         weights=10 ** d.log_gf)),
                              include_groups=False),
        "log_gf_sum": g.log_gf.apply(lambda s: float(np.log10(np.sum(10.0 ** s)))),
        "log_gf_max_component": g.log_gf.max(),
        "ep_eV": g.excitation_potential_eV.min(),
        "central_depth": g.central_depth.max(),
        "hfs_n_components": g.size(),
        "blend_flag": g.blend_flag.any(),
        "vald_proximity_min": g.vald_proximity_flag.min(),
        "loggf_source": g.loggf_source.first(),
    }).reset_index(drop=True).sort_values("wave_air_A").reset_index(drop=True)
    return out


def verify_against_canonical(feats: pd.DataFrame, can: pd.DataFrame) -> dict:
    """POSITIVE CONTROL for the HFS collapse. `canonical_gf.csv` independently carries
    `hfs_n_components` and an already-summed `log_gf` for the Al lines it covers, from
    the GES/VALD seed rather than from this grouping. If the collapse is right, the two
    agree; if it over- or under-merges, they diverge and this says so."""
    hits, gf_ok, n_ok, bad = 0, 0, 0, []
    for _, c in can.iterrows():
        d = (feats.wave_air_A - float(c.wavelength_air_A)).abs()
        if not len(d) or d.min() > 0.05:
            continue
        f = feats.loc[d.idxmin()]
        if abs(float(f.ep_eV) - float(c.excitation_potential_eV)) > HFS_EP_TOL_EV:
            continue
        hits += 1
        dgf = abs(float(f.log_gf_sum) - float(c.log_gf))
        same_n = (pd.isna(c.hfs_n_components)
                  or int(c.hfs_n_components) == int(f.hfs_n_components))
        gf_ok += dgf <= 0.005
        n_ok += bool(same_n)
        if dgf > 0.005 or not same_n:
            bad.append(dict(wave=round(float(f.wave_air_A), 4),
                            feat_gf=round(float(f.log_gf_sum), 4),
                            canonical_gf=round(float(c.log_gf), 4),
                            feat_n=int(f.hfs_n_components),
                            canonical_n=(None if pd.isna(c.hfs_n_components)
                                         else int(c.hfs_n_components))))
    return dict(n_compared=hits, n_gf_agree=gf_ok, n_ncomp_agree=n_ok,
                disagreements=bad)


# ─────────────────────────────────────────────────────────────────────────────
# 2. BAND
# ─────────────────────────────────────────────────────────────────────────────
def band_of(w: float) -> tuple[str, str]:
    """(band_name, permitted_methods) from the declared band policy, or the honest
    'no declared policy' answer. 492 of our 979 Al I linelist rows sit below 3000 A,
    where `band_policy` has NO entry — that is a real census result, not an error to
    swallow: a regime with no declared policy has no validated method."""
    try:
        p = band_policy.resolve(w)
        return p.name, "|".join(p.permitted_methods)
    except band_policy.BandPolicyError:
        return ("FUV/<3000A (NO DECLARED POLICY)" if w < 3000.0
                else f">{band_policy.POLICIES[-1].hi_A:.0f}A (NO DECLARED POLICY)"), ""


# ─────────────────────────────────────────────────────────────────────────────
# 3. gf SOURCES
# ─────────────────────────────────────────────────────────────────────────────
def _nearest(df, w, ep, wcol, epcol, tol=MATCH_TOL_A, eptol=MATCH_TOL_EP_EV):
    """Nearest row within tolerance in BOTH wavelength and lower-level energy."""
    if df is None or not len(df):
        return None
    d = (df[wcol] - w).abs()
    ok = d <= tol
    if epcol is not None and epcol in df.columns:
        ok = ok & ((df[epcol] - ep).abs() <= eptol)
    cand = df[ok]
    if not len(cand):
        return None
    return cand.loc[(cand[wcol] - w).abs().idxmin()]


def load_sources() -> dict:
    can = pd.read_csv(ROOT / "data/linelists/canonical_gf.csv", low_memory=False)
    can = can[can.species.astype(str).str.startswith("Al")].copy()
    nist_p = ROOT / "data/linelists/primary_gf/nist_asd_AlI_6600_8800.tsv"
    nist = pd.read_csv(nist_p, sep="\t") if nist_p.exists() else None
    if nist is not None:
        # NIST lists the fine-structure/HFS components of one transition on separate
        # rows, so they must be summed the way the feature collapse sums the linelist.
        #
        # ⚠️ GROUP ON THE RITZ WAVELENGTH, NOT THE OBSERVED ONE. At 6906 the two
        # gi=6 components share ritz 6906.279 but one carries wavelength_obs 6906.4
        # — 0.121 A away. Keyed on the observed value they never group, and a
        # nearest-match then returns the E-graded 4% component alone (log gf -2.481)
        # for a line whose total is -1.159. Same transition, wrong number, wrong grade.
        nist = nist[nist.log_gf.notna()].copy()
        nist["_ritz"] = nist.wavelength_ritz_A.fillna(nist.wavelength_A)
        nist = nist.sort_values(["ei_eV", "_ritz"])
        brk = ((nist._ritz.diff().abs().fillna(9e9) > 0.01)
               | (nist.ei_eV.diff().abs().fillna(9e9) > 1e-4))
        nist["_k"] = brk.cumsum()

        def _agg(d):
            gf = 10.0 ** d.log_gf.values
            tot = float(gf.sum())
            frac = gf / tot
            # sigma of the SUM: each component's fractional accuracy, gf-weighted.
            # Quoting the worst component's grade for the total would call 6906 an
            # "E" line when 96% of its gf is C+; quoting the strongest would hide the
            # rest. Weighting is the honest middle and is stated, not inferred.
            pct = np.array([NIST_ACC_PCT.get(str(g), 100.0) for g in d.nist_grade])
            eff = float((frac * pct).sum())
            grades = [str(g) for g in d.nist_grade if str(g) in NIST_ACC_PCT]
            worst = max(grades, key=lambda g: NIST_ACC_PCT[g]) if grades else "E"
            dom = str(d.nist_grade.iloc[int(np.argmax(gf))])
            return pd.Series(dict(
                wavelength_A=float(d._ritz.iloc[int(np.argmax(gf))]),
                ei_eV=float(d.ei_eV.iloc[0]),
                log_gf=float(np.log10(tot)),
                nist_grade=dom, nist_grade_worst=worst,
                nist_eff_pct=eff, ref_line=str(d.ref_line.iloc[0]),
                n_components=int(len(d))))
        nist = nist.groupby("_k").apply(_agg, include_groups=False).reset_index(drop=True)
    b = BURHEIM_2023.copy()
    b["lam_air_A"] = vac_to_air(b.lam_vac_A.values)
    b["ep_eV"] = np.nan     # Burheim tabulates term labels, not EP; matched on lambda
    return dict(canonical=can, nist=nist, burheim=b)


def assign_unique(feats: pd.DataFrame, src: pd.DataFrame, wcol: str,
                  tol: float) -> dict[int, pd.Series]:
    """Assign each SOURCE line to at most ONE feature, nearest first.

    ⚠️ A per-feature `nearest within tolerance` scan is not a matching — it lets one
    source line be claimed by several features. Seen twice in the first run of this
    census: Burheim's 11258.008 (air 11254.926) was awarded BOTH to 11254.926 and to
    the separate, three-times-shallower 11254.891 line 0.035 A away, and NIST's
    8773.896 to both 8773.896 and 8773.899. Each time the weaker line inherited a
    log gf ~1.2 dex too strong and a 2% error bar it has no claim to. Greedy nearest-
    pair assignment with each side consumed once removes the whole class.
    """
    pairs = []
    for si, srow in src.iterrows():
        d = (feats.wave_air_A - float(srow[wcol])).abs()
        for fi in d[d <= tol].index:
            pairs.append((float(d[fi]), fi, si))
    pairs.sort()
    used_f, used_s, out = set(), set(), {}
    for dist, fi, si in pairs:
        if fi in used_f or si in used_s:
            continue
        used_f.add(fi); used_s.add(si)
        out[fi] = src.loc[si]
    return out


def gf_ladder(row, src, bm=None, nm=None) -> dict:
    """Best gf source for one line, on PROVENANCE. Returns the whole ladder so a
    reader can see what was passed over and why. `bm`/`nm` are the UNIQUE Burheim /
    NIST assignments from `assign_unique` — never a per-row nearest scan."""
    w, ep = float(row.wave_air_A), float(row.ep_eV)
    r = {}

    if bm is not None:
        r["burheim_lam_vac_A"] = float(bm.lam_vac_A)
        r["burheim_lam_air_A"] = float(bm.lam_air_A)
        r["burheim_gap_A"] = float(bm.lam_air_A - w)
        r["burheim_log_gf"] = float(bm.loggf_exp)
        r["burheim_unc_pct"] = float(bm.unc_pct)
        r["burheim_sigma_dex"] = pct_to_dex(float(bm.unc_pct))
        r["burheim_transition"] = str(bm.transition)
        r["papoulia19_log_gf"] = float(bm.loggf_P19)
        r["kurucz95_log_gf"] = float(bm.loggf_K95)
        r["topbase_m00_log_gf"] = float(bm.loggf_M00)

    if nm is not None:
        r["nist_lam_A"] = float(nm.wavelength_A)
        r["nist_gap_A"] = float(nm.wavelength_A - w)
        r["nist_log_gf"] = float(nm.log_gf)
        r["nist_grade"] = str(nm.nist_grade)
        r["nist_grade_worst"] = str(nm.nist_grade_worst)
        r["nist_sigma_dex"] = pct_to_dex(float(nm.nist_eff_pct))
        r["nist_n_components"] = int(nm.n_components)
        r["nist_ref_line"] = str(nm.ref_line)

    cm = _nearest(src["canonical"], w, ep, "wavelength_air_A", "excitation_potential_eV")
    if cm is not None:
        r["canonical_log_gf"] = float(cm.log_gf)
        r["canonical_reference"] = str(cm.loggf_reference)
        r["canonical_gf_tier"] = str(cm.gf_tier)
        r["canonical_hfs_n"] = (int(cm.hfs_n_components)
                                if pd.notna(cm.hfs_n_components) else None)
        r["canonical_lam_A"] = float(cm.wavelength_air_A)

    # ── the ladder ───────────────────────────────────────────────────────────
    ref = str(r.get("canonical_reference", row.loggf_source or ""))
    if "burheim_log_gf" in r:
        r["best_gf_source"] = "EXP-BURHEIM23"
        r["best_log_gf"] = r["burheim_log_gf"]
        r["best_sigma_dex"] = r["burheim_sigma_dex"]
    elif "nist_grade" in r and r["nist_grade"] in ("AAA", "AA", "A+", "A", "B+", "B",
                                                   "C+", "C"):
        r["best_gf_source"] = f"NIST-{r['nist_grade']}"
        r["best_log_gf"] = r["nist_log_gf"]
        r["best_sigma_dex"] = r["nist_sigma_dex"]
    elif ref.startswith("1995JPhB"):
        r["best_gf_source"] = "THEORY-OP95"
        r["best_log_gf"] = r.get("canonical_log_gf")
        r["best_sigma_dex"] = OP95_SYSTEMATIC_DEX
    elif "nist_grade" in r:
        r["best_gf_source"] = f"NIST-{r['nist_grade']}"
        r["best_log_gf"] = r["nist_log_gf"]
        r["best_sigma_dex"] = r["nist_sigma_dex"]
    elif ref.upper().startswith("K") and ref.upper() != "KP":
        r["best_gf_source"] = "SEMIEMP-K75"
        r["best_log_gf"] = r.get("canonical_log_gf")
        r["best_sigma_dex"] = K75_SYSTEMATIC_DEX
    else:
        r["best_gf_source"] = "UNGRADED"
        r["best_log_gf"] = float(row.log_gf_sum)
        r["best_sigma_dex"] = K75_SYSTEMATIC_DEX

    # SCALE-MISMATCH (the RYA-799/825 class): a better gf exists in-repo and the
    # linelist the pool would be measured on does NOT carry it.
    bl = r.get("best_log_gf")
    if bl is not None and np.isfinite(bl):
        r["delta_best_minus_linelist"] = float(bl) - float(row.log_gf_sum)
        r["scale_mismatch"] = bool(abs(r["delta_best_minus_linelist"]) > 0.02)
    return r


# ─────────────────────────────────────────────────────────────────────────────
# 4. COVERAGE + TIER
# ─────────────────────────────────────────────────────────────────────────────
#: Holdings that `pipeline.coverage` structurally cannot see (their manifest is not a
#: spectrum-location CSV, so `load_registry` `continue`s past them — SILENTLY). Their
#: real-pixel coverage is measured separately by scripts/rya1001_crires_coverage.py.
#: Listed here so the census never reports "no instrument" for a line we in fact hold.
COVERAGE_BLIND_SPOT = {
    "solar_vesta_crires_plus_idp": (9479.3, 24855.0, "crires_plus"),
    "solar_crires_plus_y_rya794": (10280.0, 10680.0, "crires_plus"),
    "solar_harps_molecfit_corrected": (3782.6, 6910.0, "harps"),
    "solar_kpno_molecfit_corrected": (2960.0, 13000.0, "kpno_solar_atlas"),
    "solar_kpno_kurucz2005_corrected": (2960.0, 13000.0, "kpno_solar_atlas"),
}


def tier_line(row, inst_ids: tuple[str, ...]) -> str:
    """Per-line disposition, declared in advance, on OBSERVABLE properties only.

      EXCLUDED-NO-HOME     no registered holding covers it
      EXCLUDED-SHALLOW     predicted central depth < DEPTH_LO — no signal to carry
      EXCLUDED-SATURATED   predicted central depth > DEPTH_HI — off the linear COG
      EXCLUDED-NO-POLICY   band has no declared analytical method
      CANDIDATE-BLENDED    in-window, in-depth, but blend_flag / VALD proximity
      GRADEABLE            in-window, in-depth, unblended, and has a home
    """
    if not inst_ids:
        return "EXCLUDED-NO-HOME"
    d = float(row.central_depth)
    if d < DEPTH_LO:
        return "EXCLUDED-SHALLOW"
    if d > DEPTH_HI:
        return "EXCLUDED-SATURATED"
    if not row.band_methods:
        return "EXCLUDED-NO-POLICY"
    if bool(row.blend_flag) or float(row.vald_proximity_min) < 0.10:
        return "CANDIDATE-BLENDED"
    return "GRADEABLE"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="RYA-1001 Al Phase-0 census")
    ap.add_argument("--star", default="solar")
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args(argv)

    ll = pd.read_csv(ROOT / "data/linelists/linelist_solar.csv", low_memory=False)
    feats = features(ll, "Al")
    src = load_sources()

    cov_cache: dict[float, tuple[str, ...]] = {}

    def reach(w: float) -> tuple[tuple[str, ...], tuple[str, ...]]:
        """(all covering instrument_ids, those whose data is actually on THIS host).
        The split matters: the Kitt Peak atlas is registered host=mac, so a line whose
        only cover is KP cannot be measured on Sirius today."""
        k = round(float(w), 1)
        if k not in cov_cache:
            try:
                c = coverage_at(w, args.star).covering
                cov_cache[k] = (tuple(i.instrument_id for i in c),
                                tuple(i.instrument_id for i in c if i.host == "sirius"))
            except CoverageError:
                cov_cache[k] = ((), ())
        return cov_cache[k]

    ctrl = verify_against_canonical(feats, src["canonical"])
    print(f"[positive control] HFS collapse vs canonical_gf: {ctrl['n_gf_agree']}"
          f"/{ctrl['n_compared']} summed log gf agree to 0.005 dex, "
          f"{ctrl['n_ncomp_agree']}/{ctrl['n_compared']} component counts agree")
    for b in ctrl["disagreements"][:15]:
        print(f"    DISAGREE {b}")

    burheim_of = assign_unique(feats, src["burheim"], "lam_air_A", tol=0.03)
    nist_of = (assign_unique(feats, src["nist"], "wavelength_A", tol=MATCH_TOL_A)
               if src["nist"] is not None else {})

    rows = []
    for fi, f in feats.iterrows():
        band, methods = band_of(float(f.wave_air_A))
        seen, here = reach(float(f.wave_air_A))
        blind = tuple(sorted({v[2] for h, v in COVERAGE_BLIND_SPOT.items()
                              if v[0] <= f.wave_air_A <= v[1]} - set(seen)))
        rec = dict(
            element="Al", ion=f.ion, wave_air_A=round(float(f.wave_air_A), 4),
            wave_vac_A=round(float(air_to_vac(np.array([f.wave_air_A]))[0]), 4),
            band=band, band_methods=methods, ep_eV=round(float(f.ep_eV), 4),
            log_gf_linelist_sum=round(float(f.log_gf_sum), 4),
            log_gf_strongest_component=round(float(f.log_gf_max_component), 4),
            hfs_n_components=int(f.hfs_n_components),
            hfs_carried=bool(f.hfs_n_components > 1),
            central_depth=round(float(f.central_depth), 4),
            blend_flag=bool(f.blend_flag),
            vald_proximity_min=round(float(f.vald_proximity_min), 4),
            instruments_coverage_module="|".join(seen),
            instruments_on_sirius="|".join(here),
            instruments_coverage_blind_spot="|".join(blind),
        )
        rec.update(gf_ladder(f, src, burheim_of.get(fi), nist_of.get(fi)))
        rec["tier"] = tier_line(pd.Series(rec | dict(central_depth=f.central_depth,
                                                     blend_flag=f.blend_flag,
                                                     vald_proximity_min=f.vald_proximity_min,
                                                     band_methods=methods)),
                                seen + blind)
        # telluric: the BAND policy decides whether a correction stage is owed at all
        try:
            rec["telluric_required_band"] = bool(band_policy.resolve(f.wave_air_A).telluric_required)
        except band_policy.BandPolicyError:
            rec["telluric_required_band"] = None
        rows.append(rec)

    df = pd.DataFrame(rows)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "rya1001_al_line_census.csv", index=False)

    # ── summary ──────────────────────────────────────────────────────────────
    bands = df.groupby("band", sort=False)
    summary = bands.apply(lambda g: pd.Series({
        "n_features": len(g),
        "n_AlI": int((g.ion == "I").sum()),
        "n_AlII": int((g.ion == "II").sum()),
        "n_AlIII": int((g.ion == "III").sum()),
        "n_hfs_carried": int(g.hfs_carried.sum()),
        "n_with_home": int((g.instruments_coverage_module.astype(bool)
                            | g.instruments_coverage_blind_spot.astype(bool)).sum()),
        "n_home_on_sirius": int((g.instruments_on_sirius.astype(bool)
                                 | g.instruments_coverage_blind_spot.astype(bool)).sum()),
        "GRADEABLE": int((g.tier == "GRADEABLE").sum()),
        "CANDIDATE-BLENDED": int((g.tier == "CANDIDATE-BLENDED").sum()),
        "EXCLUDED-SHALLOW": int((g.tier == "EXCLUDED-SHALLOW").sum()),
        "EXCLUDED-SATURATED": int((g.tier == "EXCLUDED-SATURATED").sum()),
        "EXCLUDED-NO-HOME": int((g.tier == "EXCLUDED-NO-HOME").sum()),
        "EXCLUDED-NO-POLICY": int((g.tier == "EXCLUDED-NO-POLICY").sum()),
        "n_burheim_graded": int(g.get("burheim_log_gf", pd.Series(dtype=float)).notna().sum())
        if "burheim_log_gf" in g else 0,
    }), include_groups=False).reset_index()
    summary.to_csv(out / "rya1001_band_summary.csv", index=False)

    print(summary.to_string(index=False))
    print()
    grad = df[df.tier.isin(["GRADEABLE", "CANDIDATE-BLENDED"])]
    cols = ["ion", "wave_air_A", "band", "ep_eV", "log_gf_linelist_sum",
            "hfs_n_components", "central_depth", "best_gf_source", "best_log_gf",
            "best_sigma_dex", "tier", "instruments_coverage_module",
            "instruments_on_sirius", "instruments_coverage_blind_spot"]
    print("=== GRADEABLE + CANDIDATE pool ===")
    print(grad[[c for c in cols if c in grad.columns]].to_string(index=False))

    meta = dict(
        ticket="RYA-1001", phase="0 (census only — no synthesis, no gf adoption)",
        gf_ladder=["EXP-BURHEIM23", "THEORY-P19", "NIST-<grade>", "THEORY-OP95",
                   "SEMIEMP-K75", "UNGRADED"],
        depth_window=[DEPTH_LO, DEPTH_HI], group_A=GROUP_A,
        match_tol_A=MATCH_TOL_A, match_tol_ep_eV=MATCH_TOL_EP_EV,
        n_linelist_rows=int((ll.element == "Al").sum()), n_features=len(df),
        hfs_collapse_control=ctrl,
        firewall="RYA-161 — nothing here adopts a gf; the ranking is provenance-only.",
    )
    (out / "rya1001_census_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"\nwrote {out}/rya1001_al_line_census.csv, rya1001_band_summary.csv, "
          f"rya1001_census_meta.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
