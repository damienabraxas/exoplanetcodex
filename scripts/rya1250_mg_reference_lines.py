"""Solar Mg Reference Grade lines, repaired at the level the synthesis reads them (RYA-1250 fix).

The synthesis does not read canonical_gf per published component. `gf_resolver.apply_to_synth_array`
clusters each synthesis list into physical lines (same lower level, gaps < HFS_GAP), then scales
every cluster to the ONE canonical row nearest its centroid. A published component set is
therefore only honoured when the canonical row of each cluster carries the published TOTAL of the
components inside that cluster. PR #585 wrote component gf onto whole clusters and added
component rows that no cluster ever resolves to:

  12083.28  cluster {.278,.346} carried only .278 (-1.347) instead of their sum (-1.116)
  10965.45  cluster {.414,.450} carried only .450
  K 19723/19733, near-UV 3092.98/3096.89: a NIST grade of the DOMINANT component (whose value
            VALD already carries exactly) was written as the line total: -0.125/-0.076/-0.125/-0.076
  10811/11828/12083.66/15765.84/15886: duplicate rows the resolver never picks

and left the AGSS21 Mg I lines ungraded for want of a per-line sigma. This script:

  * sets each feature's canonical row to the published total of its components (plus the VALD
    gf of cluster components a PARTIAL source does not list), asserting the row is the one the resolver picks for that cluster
  * sigma = gf-weighted quadrature of the listed components (blank if an unlisted-sigma component
    carries > 10% of the weight)
  * tags Pehlivan Rhodin+2017 Table 2 EXPERIMENTAL values LAB (the Reference Grade lab pool);
    Bergemann+2017 / AGSS21 values that rest on PR17 theoretical gf stay PUBLISHED-SET. Sigma 0.05
    for PR17 theoretical gf is Bergemann+2017's own Table 2 convention (note a).
  * deletes the duplicate component rows, splits Mg II 10092 into its two components
    (0.122 A apart = two synthesis clusters) scaled to the AGSS21 total
  * writes the two line-set source tables at feature level

Run from the repo root; idempotent. The literature is in docs/audit/rya1250_mg_fix.md.
"""
from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from pipeline import gf_resolver as gfr  # noqa: E402

CANON = ROOT / "data/linelists/canonical_gf.csv"
LISTS = {"ges": Path(os.environ.get("ISPEC_DIR", Path.home() / "codex/ispec"))
         / "input/linelists/transitions/GESv6_atom_hfs_iso.420_920nm/atomic_lines.tsv",
         "ir": ROOT / "data/linelists/ispec_ir_9200_13000/atomic_lines.tsv",
         "h": ROOT / "data/linelists/ispec_h_15007_17494/atomic_lines.tsv",
         "k": ROOT / "data/linelists/ispec_k_19452_24846/atomic_lines.tsv"}

PR17 = ("Pehlivan Rhodin et al. 2017, A&A 598 A102", "10.1051/0004-6361/201629849")
B17 = ("Bergemann et al. 2017, ApJ 847, 15, Table 2", "10.3847/1538-4357/aa88cb")
AGSS = ("Asplund, Amarsi & Grevesse 2021, A&A 653 A141, Table A.1", "10.1051/0004-6361/202140445")

# (row, species, published components [(wl, loggf, sigma)], tier, tag, (ref, doi), gf note)
# row: an existing canonical line_id.  Components are what the SOURCE lists; VALD components of the
# same cluster that the source does not list are kept at their VALD gf.
FEATURES = [
    # -- AGSS21 Mg I: Table A.1 gf (PR17 theoretical; sigma per Bergemann+2017 note a) --------
    ("gf_087258", "Mg I", [(6318.717, -2.020, 0.05)], "PUBLISHED-SET", "MG_AGSS21_PR17calc", AGSS, "PR17 theo"),
    ("gf_087259", "Mg I", [(6319.237, -2.242, 0.05)], "PUBLISHED-SET", "MG_AGSS21_PR17calc", AGSS, "PR17 theo"),
    ("gf_087349", "Mg I", [(8712.6857, -1.152, 0.05)], "PUBLISHED-SET", "MG_AGSS21_PR17calc", AGSS, "PR17 theo, feature total"),
    ("gf_087350", "Mg I", [(8717.825, -0.930, 0.05)], "PUBLISHED-SET", "MG_AGSS21_PR17calc", AGSS, "PR17 theo, feature total"),
    ("gf_087357", "Mg I", [(8923.569, -1.679, 0.05)], "PUBLISHED-SET", "MG_AGSS21_PR17calc", AGSS, "PR17 theo"),
    ("gf_179261", "Mg I", [(9983.188, -2.177, 0.05)], "PUBLISHED-SET", "MG_AGSS21_PR17calc", AGSS, "PR17 theo"),
    ("gf_167377", "Mg I", [(10312.531, -1.718, 0.05)], "PUBLISHED-SET", "MG_AGSS21_PR17calc", AGSS, "PR17 theo"),
    ("gf_173504", "Mg I", [(15879.558, -1.226, 0.06)], "PUBLISHED-SET", "MG_AGSS21_PR17calc", AGSS,
     "AGSS21 feature total; sigma = PR17 exp error of the dominant 15879.567 component"),
    # -- Bergemann+2017 Table 2 on PR17 theoretical / Butler gf ------------------------------
    ("gf_087240", "Mg I", [(4571.096, -5.397, 0.05)], "PUBLISHED-SET", "MG_BERGEMANN2017", B17, "PR17 theo"),
    ("gf_087243", "Mg I", [(4702.991, -0.456, 0.05)], "PUBLISHED-SET", "MG_BERGEMANN2017", B17, "PR17 theo"),
    ("gf_087252", "Mg I", [(5711.088, -1.742, 0.05)], "PUBLISHED-SET", "MG_BERGEMANN2017", B17,
     "PR17 theo, as Bergemann+2017 and Magg+2022 adopt; PR17 exp is -1.842+/-0.05"),
    ("gf_167460", "Mg I", [(10811.053, 0.052, 0.04), (10811.076, -0.137, 0.043), (10811.097, -1.038, 0.079),
                           (10811.122, -1.036, 0.079), (10811.143, -2.587, 0.30), (10811.158, -0.321, 0.04)],
     "PUBLISHED-SET", "MG_BERGEMANN2017", B17, "PR17 exp (.053,.158) + Butler 1993 (10/18/50 % -> dex)"),
    ("gf_167494", "Mg I", [(10965.386, -2.184, 0.05), (10965.414, -1.008, 0.05), (10965.450, -0.260, 0.05)],
     "PUBLISHED-SET", "MG_BERGEMANN2017", B17, "PR17 theo"),
    ("gf_167651", "Mg I", [(12083.278, -1.347, 0.05), (12083.346, -1.500, None)],
     "PUBLISHED-SET", "MG_BERGEMANN2017", B17,
     "3d1D2-4f3F intercombination pair: PR17 exp .278 + Bergemann+2015 .346 (no sigma). VALD's +0.45 at .278 is refuted by PR17"),
    ("gf_173378", "Mg I", [(15748.988, 0.129, 0.05)], "PUBLISHED-SET", "MG_BERGEMANN2017", B17,
     "PR17 theo; Bergemann excludes it from the solar mean (blends ~half the EW)"),
    ("gf_173393", "Mg I", [(15765.645, -1.524, 0.05)], "PUBLISHED-SET", "MG_BERGEMANN2017", B17, "PR17 theo"),
    ("gf_173394", "Mg I", [(15765.747, -0.348, 0.05), (15765.842, 0.400, 0.05)],
     "PUBLISHED-SET", "MG_BERGEMANN2017", B17, "PR17 theo"),
    ("gf_173512", "Mg I", [(15886.183, -1.465, 0.05), (15886.261, -1.942, 0.05)],
     "PUBLISHED-SET", "MG_BERGEMANN2017", B17,
     "PR17 theo as Bergemann adopts; PR17 exp -1.896+/-0.10 / -1.760+/-0.14 (total -1.52, 0.18 weaker)"),
    # -- PR17 Table 2 EXPERIMENTAL gf: the Reference Grade lab pool --------------------------
    ("gf_087245", "Mg I", [(5167.322, -0.854, 0.05)], "LAB", "PR2017_exp", PR17, "exp"),
    ("gf_087246", "Mg I", [(5172.684, -0.363, 0.04)], "LAB", "PR2017_exp", PR17, "exp"),
    ("gf_087247", "Mg I", [(5183.604, -0.168, 0.04)], "LAB", "PR2017_exp", PR17,
     "exp; replaces +0.18 (NIST -0.180 with its sign dropped; every line list has -0.24)"),
    ("gf_087251", "Mg I", [(5528.405, -0.547, 0.02)], "LAB", "PR2017_exp", PR17, "exp"),
    ("gf_087288", "Mg I", [(7657.603, -1.292, 0.05)], "LAB", "PR2017_exp", PR17, "exp"),
    ("gf_087289", "Mg I", [(7659.152, -1.425, 0.05)], "LAB", "PR2017_exp", PR17, "exp"),
    ("gf_087290", "Mg I", [(7659.901, -1.948, 0.05)], "LAB", "PR2017_exp", PR17, "exp"),
    ("gf_087355", "Mg I", [(8806.757, -0.144, 0.03)], "LAB", "PR2017_exp", PR17, "exp"),
    ("gf_167184", "Mg I", [(9255.778, -0.187, 0.04)], "LAB", "PR2017_exp", PR17,
     "exp 1F component; VALD 9255.722 3F intercombination (5% of the cluster) kept"),
    ("gf_167615", "Mg I", [(11828.171, -0.350, 0.03)], "LAB", "PR2017_exp", PR17, "exp (PR17 11828.185)"),
    ("gf_167652", "Mg I", [(12083.649, 0.377, 0.04)], "LAB", "PR2017_exp", PR17, "exp (PR17 12083.662)"),
    ("gf_172816", "Mg I", [(15024.997, 0.334, 0.03)], "LAB", "PR2017_exp", PR17, "exp"),
    ("gf_172829", "Mg I", [(15040.246, 0.113, 0.04)], "LAB", "PR2017_exp", PR17, "exp"),
    ("gf_172836", "Mg I", [(15047.714, -0.364, 0.04)], "LAB", "PR2017_exp", PR17, "exp (PR17 15047.705)"),
    ("gf_173518", "Mg I", [(15889.521, -1.876, 0.08)], "LAB", "PR2017_exp", PR17, "exp (PR17 15889.485)"),
    # K: 4d1D2-5f1F3, PR17 24572.92 vac. NIST ASD (T5539) grades the same transition A at
    # -0.225, 0.17 dex stronger. VALD's 24565.585 (+0.110, 4d1D2-5f3F3, Kurucz) would be the
    # strongest line of the feature; the CRIRES+ K spectrum shows only a 9% shoulder there
    # against 22% at 24566.37, the 12083 pattern. Kept at VALD; read the fit's residual.
    ("gf_179143", "Mg I", [(24566.296, -0.391, 0.06)], "LAB", "PR2017_exp", PR17,
     "exp (PR17 24572.92 vac); NIST A -0.225"),
]
#: Features whose source lists only SOME of the cluster's transitions; every other feature is a
#: complete published list (AGSS21 feature totals, Bergemann's full multiplets), so VALD-only
#: components of the cluster are not added on top (VALD 10811.198/.219 are not in B17 or PR17).
PARTIAL = {"gf_167184"}
#: PR #585 rows that duplicate a cluster above (no synthesis cluster resolves to them).
DUPLICATES = [f"gf_1792{n}" for n in range(63, 77)]
#: NIST graded the dominant component only; VALD carries exactly that value, so the line total
#: is the VALD total (verified on the raw extracts 2026-10-23). The grade stays.
RESTORE_TOTAL = {"gf_148749": -0.447, "gf_148896": -0.225, "gf_176848": -0.739, "gf_176858": -0.517}
#: Mg II 10092: VALD components (= NIST ritz 1009.209 / 1009.222 nm), scaled to the AGSS21 total.
MG2_10092 = [(10092.095, 0.91), (10092.217, 1.02)]
MG2_10092_TOTAL, MG2_10092_SIGMA = 1.283, 0.01


def _synth(species: str, wl: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lst = "ges" if wl < 9200 else "ir" if wl < 13000 else "h" if wl < 19000 else "k"
    d = pd.read_csv(LISTS[lst], sep="\t", usecols=["element", "wave_A", "loggf", "lower_state_eV"],
                    low_memory=False)
    tag = f"{species.split()[0]} {1 if species.endswith(' I') else 2}"
    d = d[d.element == tag]
    return d.wave_A.to_numpy(float), d.lower_state_eV.to_numpy(float), d.loggf.to_numpy(float)


def feature_total(species: str, comps: list, partial: bool) -> tuple[float, float | None, float, float]:
    """(total log gf, sigma, cluster centroid, cluster EP) for the synthesis cluster holding comps."""
    wl, ep, gf = _synth(species, comps[0][0])
    key = [(0, 0)] * len(wl)
    cl = next((c for c in gfr.cluster_physical_lines(key, wl, ep)
               if any(abs(wl[i] - comps[0][0]) < 0.04 for i in c)), None)
    if cl is None:
        raise SystemExit(f"{species} {comps[0][0]}: no synthesis component")
    covered = {i for i in cl for w, _, _ in comps if abs(wl[i] - w) < 0.04}
    lin = [10 ** g for _, g, _ in comps] + ([10 ** gf[i] for i in cl if i not in covered] if partial else [])
    tot = float(np.sum(lin))
    if any(s is None and 10 ** g / tot > 0.10 for (_, g, s) in comps):
        sigma = None
    else:
        sigma = float(np.sqrt(sum((10 ** g / tot * (s or 0.0)) ** 2 for (_, g, s) in comps)))
    w = 10 ** gf[cl]
    return float(np.log10(tot)), sigma, float(np.sum(wl[cl] * w) / w.sum()), float(np.mean(ep[cl]))


def main() -> int:
    with CANON.open(newline="") as fh:
        rows = list(csv.reader(fh))
    head, body = rows[0], rows[1:]
    col = {k: i for i, k in enumerate(head)}
    by_id = {r[0]: r for r in body}

    def put(r, **kv):
        for k, v in kv.items():
            r[col[k]] = "" if v is None else str(v)

    report = []
    for lid, sp, comps, tier, tag, (ref, doi), note in FEATURES:
        r = by_id[lid]
        total, sigma, cen, ep = feature_total(sp, comps, lid in PARTIAL)
        # the row must be the one the resolver picks for this cluster
        near = [x for x in body if x[col["species"]] == sp
                and abs(float(x[col["wavelength_air_A"]]) - cen) <= gfr.WTOL
                and abs(float(x[col["excitation_potential_eV"]]) - ep) <= gfr.EPTOL
                and x[0] not in DUPLICATES]
        pick = min(near, key=lambda x: abs(float(x[col["wavelength_air_A"]]) - cen))[0] if near else None
        if pick != lid:
            raise SystemExit(f"{sp} {cen:.3f}: resolver picks {pick}, not {lid}")
        # nist_grade is a fact about NIST's value and other tiers read it: keep it, and keep the
        # superseded reference so the cross-check source stays on record (RYA-592's rule)
        prior = r[col["loggf_reference"]].split(" | prior: ")[-1]
        put(r, log_gf=round(total, 3), gf_sigma_dex=None if sigma is None else round(sigma, 4),
            gf_tier=tier, lab_source_tag=tag, gf_source_doi=doi,
            loggf_reference=f"{ref} ({note})" + (f" | prior: {prior}" if prior else ""),
            adjudication_status="rya1250_fix_feature_total")
        report.append((lid, sp, round(cen, 3), round(total, 3), sigma, tier))
    for lid, total in RESTORE_TOTAL.items():
        r = by_id[lid]
        put(r, log_gf=total, loggf_reference=r[col["loggf_reference"]].split(" (")[0]
            + " (grade of the dominant component; line total = VALD components, RYA-1250 fix)")
    # Mg II 10092: two components, two clusters
    r0 = by_id["gf_179262"]
    shift = MG2_10092_TOTAL - float(np.log10(sum(10 ** g for _, g in MG2_10092)))
    nxt = 1 + max(int(x[0].split("_")[1]) for x in body)
    new = []
    from pipeline.physical_line_id import physical_id
    for j, (w, g) in enumerate(MG2_10092):
        # A re-run finds the component row it added before (the docstring's idempotence:
        # appending each time stacked duplicate 10092.217 rows, RYA-1250 2026-10-24).
        prev = [x for x in body if x[0] != r0[0] and x[col["species"]] == "Mg II"
                and abs(float(x[col["wavelength_air_A"]]) - w) < 1e-4]
        if j and prev:
            r = prev[0]
            body = [x for x in body if x not in prev[1:]]
        elif j:
            r = list(r0)
            r[0] = f"gf_{nxt:06d}"
            new.append(r)
        else:
            r = r0
        put(r, wavelength_air_A=w, log_gf=round(g + shift, 3), hfs_n_components=1,
            physical_id=physical_id("Mg II", w, 11.6297), gf_sigma_dex=MG2_10092_SIGMA,
            gf_linelist_vald=g, adjudication_status="rya1250_fix_component_of_agss21_feature",
            loggf_reference=f"{AGSS[0]} (10092.20 feature {MG2_10092_TOTAL} split over the two "
                            f"NIST/VALD components in VALD ratio; sigma = AGSS21 Sect. 2.4 Mg II bound)")
    body = [x for x in body if x[0] not in DUPLICATES] + new
    with CANON.open("w", newline="") as fh:
        csv.writer(fh, lineterminator="\n").writerows([head] + body)
    for x in report:
        print(*x)
    print(f"deleted {len(DUPLICATES)} duplicate rows; restored {len(RESTORE_TOTAL)} line totals; "
          f"Mg II 10092 -> 2 components (shift {shift:+.3f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
