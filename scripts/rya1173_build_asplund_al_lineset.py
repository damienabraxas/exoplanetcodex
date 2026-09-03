#!/usr/bin/env python3
"""RYA-1173 - build the AGSS21-lineage solar Al reference line set from its primary sources.

    python3 scripts/rya1173_build_asplund_al_lineset.py [--check]

RYA-1109 built `data/reference/asplund2021_fe/` by transcribing AGSS21's OWN Table A.2. That route
does not exist for Al: **AGSS21 publishes no Al line list.** Its Aluminium section adopts Nordlander
& Lind (2017), who "adopted the same lines and line data as in Scott et al. (2015b), except that
they excluded the 1089.1 nm line due to telluric contamination". So the reference set is
reconstructed from two primaries, and this script is the reconstruction.

    AGSS21 (A&A 653, A141)          adopts the VALUE, publishes no lines
      -> Nordlander & Lind 2017     the analysis: six lines, one published exclusion
         (A&A 607, A75)             Table A.1 = the line DATA (levels, EP, log gf, sigma)
           -> Scott et al. 2015b    the SELECTION: seven lines, EWs, weights, per-model abundances
              (A&A 573, A25)        Table 2 + Table 3 (level J identity + HFS)

🔴 THE SET IS SIX USED + ONE EXCLUDED, AND THE EXCLUDED ROW SHIPS. RYA-946: "Preserve explicit
negative selections and scope statements." 10891.732 A is IN this file with
`selection_status = EXCLUDED_BY_SOURCE_ANALYSIS` and the published reason. Dropping it would erase
a published decision and make the set indistinguishable from a set that never had the line.

🔴 AGSS21 SAYS "FIVE" AND THAT NUMBER REPRODUCES FROM NOTHING. Scott retains seven; NL2017 removes
one and Fig. 8 names the remaining six individually (6696 6698 7835 8912 10768 10872). This script
carries SIX on the authority of the primary and records the conflict in
`published_line_count_conflict` -- it does not pick a number quietly. See `raw/lineage_quotations.md`.

🔴 THE TWO PRIMARIES ARE JOINED ON THE LEVEL, NOT THE WAVELENGTH. Scott prints lambda in nm to 4 dp
and NL2017 in A to 3 dp, and they DISAGREE -- 6696.023 vs 6696.015 A, 8 mA apart, which is 16x the
`line_match` default window. Joining them on wavelength would have to invent a tolerance. Both tables
publish the lower and upper level term AND its J (Scott's Table 3 exists precisely for that), and the
level identity is exact and unambiguous, so that is the key. RYA-1151's lesson, applied at the source
rather than downstream: THE FIX IS THE LEVEL.

EXTRACTION CONTROLS (see `controls` in the .prov.json). A transcription must be checked against
claims the source makes about columns the transcription does NOT go on to use, and against an
identity the source never tabulates -- a silent column shift survives every check built from the
headline number (RYA-1110).

  C1  Scott's Table 1 publishes logeps(Al) for FIVE model atmospheres plus two differences.
      The weight-weighted means of the five transcribed LTE columns must reproduce all five, and
      3D-HM and 3D-<3D> must reproduce too. Seven checks, and the weight column and the model
      columns are used by nothing downstream. This also fixes the DIRECTION of Scott's weight
      scale: on w=Wt the 3D+NLTE mean is 6.4298 -> 6.43 as published; unweighted it is 6.4187,
      and on inverted weights 6.4076. Only "larger weight = better line" reproduces the paper.
  C2  Per line, a_nlte_3d - a_lte_3d == delta_nlte_3d (a column relation Scott never states).
  C3  Cross-paper: on the LEVEL key, NL2017 Table A.1 and Scott Table 2 must agree EXACTLY on
      log gf and on E_low for all seven lines. Two independently typeset tables in two papers.
  C4  An identity NEITHER paper tabulates: E_up = E_low + hc/lambda_vac must be single-valued
      per upper level, and E_low single-valued per lower level, across all 55 rows of Table A.1.
      It is simultaneously the AIR/VACUUM control -- the levels say the column is air, 4.5x
      tighter than the vacuum reading -- and it reports where its own sensitivity runs out.
  C5  The six-line set derived here must be exactly the six NL2017's Fig. 8 axis names.
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

REF = ROOT / "data" / "reference" / "asplund2021_al"
RAW = REF / "raw"
OUT_LINES = REF / "asplund2021_al_lines.csv"
OUT_ANALYSIS = REF / "nordlander_lind_2017_analysis_lines.csv"
OUT_PROV = REF / "asplund2021_al_lines.prov.json"

#: 🔴 RYA-1084: artifacts must not carry raw float repr. Every derived float here comes out of a
#: division or a weighted mean, and one ULP of numpy drift becomes a visible byte change.
FLOAT_DECIMALS = 10

#: hc in eV.A -- the same value RYA-1109 used for the Fe set's energy-residual control.
HC_EV_A = 12398.419843320026


def eup_eV(elo_eV, wavelength_air_A):
    """E_up = E_low + hc/lambda_VACUUM. Never hc/lambda_air.

    🔴 E IS A VACUUM QUANTITY AND THE PRINTED WAVELENGTHS ARE AIR. Skipping the conversion is a
    quiet 0.0003-0.0017 eV error that looks like nothing and is 4.5x the residual C4 tolerates --
    C4 measures exactly that gap, so an air-based E_up here would contradict the control that
    validated the column. Conversion is `pipeline.wavelength_util`, the declared SSOT.
    """
    from pipeline.wavelength_util import air_to_vac
    return np.asarray(elo_eV, float) + HC_EV_A / air_to_vac(np.asarray(wavelength_air_A, float))

#: Scott et al. 2015b's abundance summary table, the Al i row -- the numbers C1 must reproduce.
#:
#: ⚠️ THE SAME TABLE HAS TWO NUMBERS IN THE TWO COPIES. It is **Table 1** in the published A&A
#: article (aanda.org .../aa24109-14/T1.html) and **Table 5** in arXiv:1405.0279v2. Citing "Table 1"
#: at someone holding the preprint sends them to the model stratification.
#:
#: 🔴 EVERY COLUMN HERE IS AN NLTE RESULT, NOT AN LTE ONE. Its caption, verbatim: "A summary of the
#: NLTE results obtained in this analysis with our 3D model, and with the four different 1D models
#: we used." So the `HM` entry is not the mean of Table 2's HM column -- it is the mean of
#: (HM + dNLTE), the single 3D-derived NLTE correction applied to every model. Reading these as LTE
#: columns puts all five exactly 0.01 dex low while the two DIFFERENCE columns still reproduce,
#: which is precisely the shape of a wrong reading that looks nearly right. Measured, not assumed:
#: under the LTE reading 0 of 5 reproduce; under the NLTE reading 5 of 5 do, exactly.
#:
#: The caption also warns about its own rounding: "because all means were computed using abundances
#: accurate to three decimal places, entries in columns 8 and 9 differ in some cases from the
#: differences between the entries in columns 3-5." For Al they do not differ.
#: 3D  <3D>  HM  MARCS  MISS   3D-HM  3D-<3D>   recommended    meteoritic
SCOTT2015B_SUMMARY_AL = {
    "a_lte_3d": 6.43, "a_lte_mean3d": 6.45, "a_lte_hm": 6.48,
    "a_lte_marcs": 6.40, "a_lte_miss": 6.47,
    "3D_minus_HM": -0.05, "3D_minus_mean3D": -0.02,
    "recommended": 6.43, "recommended_sigma": 0.04,
    "meteoritic": 6.43, "meteoritic_sigma": 0.01,
}

#: Nordlander & Lind 2017 Fig. 8's x-axis, verbatim: the six lines of the solar abundance set.
NL2017_FIG8_LINES = (6696, 6698, 7835, 8912, 10768, 10872)

#: The gf reference legend of NL2017 Table A.1, verbatim from its Notes.
NL_GF_REFS = {
    "1": "Wiese, Smith & Miles 1969, NSRDS-NBS 22 (Sodium through Calcium)",
    "2": ("TOPbase: Mendoza et al. 1995, and C. Mendoza, W. Eissner, M. Le Dourneuf & "
          "C. J. Zeippen (unpublished) -- Opacity Project, THEORY"),
    "3": "Davidson et al. 1990",
    "4": "Vujnovic et al. 2002, A&A 388, 704",
    "5": "Tachiev & Froese Fischer 2002, A&A 385, 716",
    "6": "Kurucz 2012, online data, http://kurucz.harvard.edu/atoms/1300/",
}

#: 🔴 ROLES ARE ASSIGNED ONLY WHERE THE PAPER SAYS SO, PER LINE. NL2017 Table A.1's own note is
#: "Only lines used in the spectrum analyses are listed" -- a statement about the TABLE, not about
#: any row. These are the rows the running text names individually, with the sentence that names
#: them. Every other row gets ANALYSIS_LINE_ROLE_NOT_STATED_PER_LINE. Inferring a role from a
#: wavelength's familiarity is exactly what RYA-946 forbids.
NL_NAMED_ROLES = {
    3944.006: ("SOLAR_DIAGNOSTIC_NAMED_IN_TEXT",
               "Sect. 3.3: 'We also include the resonance line at 3944 A, which is heavily blended "
               "with strong lines of CH in carbon-rich stars and therefore not normally recommended "
               "for use in abundance analyses.'"),
    3961.520: ("SOLAR_DIAGNOSTIC_NAMED_IN_TEXT",
               "Sect. 3.1.1 is titled 'The solar 3961 A resonance line'; Fig. 3 compares the KPNO "
               "flux atlas to NLTE synthesis for it."),
    7835.309: ("SOLAR_DIAGNOSTIC_NAMED_IN_TEXT",
               "Sect. 3.1.2 'Center-to-limb variations in the 7835 A line'; also one of the two "
               "lines cited as constraining electron collisional rates. ALSO a member of the "
               "solar abundance set."),
    13123.416: ("SOLAR_DIAGNOSTIC_NAMED_IN_TEXT",
                "Abstract/Sect. 3.1.4: cited as the line indicating hydrogen collisional rates, and "
                "Fig. 7 shows it as one of 'two near-IR lines sensitive to hyperfine splitting'."),
    16750.519: ("SOLAR_DIAGNOSTIC_NAMED_IN_TEXT",
                "Fig. 7, the second of the 'two near-IR lines sensitive to hyperfine splitting': "
                "'assuming A(Al) = 6.43 at 13123 A and 6.33 at 16750 A for the Sun'."),
    123349.6: ("SOLAR_DIAGNOSTIC_NAMED_IN_TEXT",
               "Sect. 3.1.3 'The 12.33 micron emission line'; the abstract names it alongside 7835 A "
               "as what the 3D NLTE modelling reproduces."),
}


class BuildError(RuntimeError):
    """A control failed. The artifact is NOT written."""


def _read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="#", dtype=str).apply(
        lambda s: s.str.strip() if s.dtype == object else s)


def _f(s: pd.Series) -> np.ndarray:
    return pd.to_numeric(s, errors="coerce").to_numpy(float)


def _gf_ref(code: str) -> str:
    """'2' or '2,3' -> the resolved citation(s). NL2017 cites two sources on some rows and the
    combination is part of the provenance, so it is joined, never reduced to the first."""
    parts = [c.strip() for c in str(code).split(",") if c.strip()]
    unknown = [c for c in parts if c not in NL_GF_REFS]
    if unknown:
        raise BuildError(f"NL2017 Table A.1 gf reference code(s) {unknown} are not in its own "
                         f"legend {sorted(NL_GF_REFS)}")
    return " + ".join(NL_GF_REFS[c] for c in parts)


def _level_key(term: str, j: str) -> str:
    """'4p 2Po' + '3/2' -> '4p 2Po J=3/2'. The join key between the two primaries."""
    return f"{' '.join(str(term).split())} J={str(j).strip()}"


def load_sources() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    nl = _read_tsv(RAW / "nordlander_lind_2017_tableA1.tsv")
    sc = _read_tsv(RAW / "scott2015b_table2_al.tsv")
    hf = _read_tsv(RAW / "scott2015b_table3_al_hfs.tsv")
    if len(nl) != 55:
        raise BuildError(f"NL2017 Table A.1 has {len(nl)} rows, expected 55")
    if len(sc) != 7 or len(hf) != 7:
        raise BuildError(f"Scott Al block has {len(sc)}/{len(hf)} rows, expected 7/7")
    return nl, sc, hf


# ── controls ──────────────────────────────────────────────────────────────────────────

def control_c1_scott_table1(sc: pd.DataFrame) -> dict:
    """Scott's five model columns + two differences + the recommended value, reproduced.

    Each entry is the Wt-weighted mean of (that model's LTE abundance + dNLTE_3D) -- see the
    SCOTT2015B_SUMMARY_AL note for why the dNLTE term is there and what it costs to omit it.
    """
    w = _f(sc["weight"])
    dn = _f(sc["delta_nlte_3d"])
    got = {}
    for col in ("a_lte_3d", "a_lte_mean3d", "a_lte_hm", "a_lte_marcs", "a_lte_miss"):
        got[col] = float(np.average(_f(sc[col]) + dn, weights=w))
    got["3D_minus_HM"] = got["a_lte_3d"] - got["a_lte_hm"]
    got["3D_minus_mean3D"] = got["a_lte_3d"] - got["a_lte_mean3d"]
    #: Table 2's own 3D+NLTE column, averaged the same way. C2 checks it equals a_lte_3d + dNLTE
    #: per line, so this must equal got["a_lte_3d"] -- and it is the published `Recommended`.
    got["recommended"] = float(np.average(_f(sc["a_nlte_3d"]), weights=w))

    bad = []
    for k, pub in SCOTT2015B_SUMMARY_AL.items():
        if k not in got:
            continue
        if abs(round(got[k], 2) - pub) > 1e-9:
            bad.append(f"{k}: transcription gives {got[k]:.4f} -> {round(got[k], 2)}, "
                       f"Scott's summary table publishes {pub}")
    if bad:
        raise BuildError("C1 FAILED (Scott's abundance summary does not reproduce):\n  "
                         + "\n  ".join(bad))

    #: ⚠️ THE LTE READING IS CARRIED AS THE MEASURED NEGATIVE. Recording what the wrong reading
    #: gives is what makes "5 of 5 reproduce" mean something: the two readings differ, and only
    #: one of them is right. Without this the control could be passing for a trivial reason.
    lte_only = {c: float(np.average(_f(sc[c]), weights=w))
                for c in ("a_lte_3d", "a_lte_mean3d", "a_lte_hm", "a_lte_marcs", "a_lte_miss")}
    n_lte_ok = sum(abs(round(v, 2) - SCOTT2015B_SUMMARY_AL[c]) < 1e-9 for c, v in lte_only.items())
    if n_lte_ok != 0:
        raise BuildError(f"C1 is not discriminating: the LTE reading also reproduces "
                         f"{n_lte_ok} of 5 columns")

    #: 🔴 The direction of the weight scale is MEASURED here, not assumed. Only one of the three
    #: readings reproduces the published 6.43, so the transcribed weights cannot be silently
    #: inverted (which would have left every other check in this script green).
    inv = 4.0 - w                            # 1<->3, 2<->2 : "smaller is better"
    alts = {
        "weights_as_transcribed_larger_is_better": got["recommended"],
        "unweighted": float(np.mean(_f(sc["a_nlte_3d"]))),
        "weights_inverted_smaller_is_better": float(np.average(_f(sc["a_nlte_3d"]), weights=inv)),
    }
    if sum(abs(round(v, 2) - 6.43) < 1e-9 for v in alts.values()) != 1:
        raise BuildError(f"C1 weight-direction control is NOT discriminating: {alts}")
    return {"weighted_means_nlte": {k: round(v, 4) for k, v in got.items()},
            "published": SCOTT2015B_SUMMARY_AL,
            "n_checks": 8,
            "lte_reading_measured_negative": {
                "means": {k: round(v, 4) for k, v in lte_only.items()},
                "n_of_5_reproducing": n_lte_ok,
                "note": ("Every LTE-only mean lands exactly 0.01 dex below the published value "
                         "while 3D-HM and 3D-<3D> still reproduce -- a wrong reading that looks "
                         "nearly right. This is why the caption, not the arithmetic, decides.")},
            "weight_direction_discriminant": {k: round(v, 4) for k, v in alts.items()},
            "status": "PASS"}


def control_c2_nlte_column(sc: pd.DataFrame) -> dict:
    """a_nlte_3d - a_lte_3d == delta_nlte_3d, per line. A relation Scott never writes down."""
    d = _f(sc["a_nlte_3d"]) - _f(sc["a_lte_3d"])
    pub = _f(sc["delta_nlte_3d"])
    resid = np.abs(d - pub)
    if resid.max() > 5e-4:      # the columns are printed to 3 dp / 2 dp
        raise BuildError(f"C2 FAILED: max |(3D+NLTE - 3D_LTE) - dNLTE| = {resid.max():.5f} at "
                         f"{sc['lambda_nm'].iloc[int(resid.argmax())]} nm")
    return {"max_residual_dex": round(float(resid.max()), 6), "n_lines": int(len(sc)),
            "status": "PASS"}


def control_c3_cross_paper(nl: pd.DataFrame, sc: pd.DataFrame, hf: pd.DataFrame) -> dict:
    """The two primaries agree on log gf and E_low for all seven lines, joined on the LEVEL."""
    hfi = hf.set_index("lambda_nm")
    nl_key = {}
    for _, r in nl.iterrows():
        nl_key[(_level_key(r.lower_level, r.j_lower), _level_key(r.upper_level, r.j_upper))] = r

    rows, bad = [], []
    for _, r in sc.iterrows():
        h = hfi.loc[r.lambda_nm]
        # Scott Table 2 gives no level TERM, only Table 3's J. The term comes from NL2017, so the
        # key is built the other way round: find the NL row whose (J_low, J_upp) and E_low match.
        cands = [x for _, x in nl.iterrows()
                 if str(x.j_lower) == str(h.j_lower) and str(x.j_upper) == str(h.j_upper)
                 and abs(float(x.elo_eV) - float(r.elo_eV)) < 1e-9]
        # more than one level pair can share (J,J,E_low); the log gf is what separates them, and
        # requiring an EXACT match is the control, so an ambiguous case must be reported not picked.
        exact = [x for x in cands if abs(float(x.loggf) - float(r.loggf)) < 1e-9]
        if len(exact) != 1:
            bad.append(f"{r.lambda_nm} nm: {len(cands)} candidates on (J_low={h.j_lower}, "
                       f"J_upp={h.j_upper}, E_low={r.elo_eV}), {len(exact)} of them agree on "
                       f"log gf {r.loggf}")
            continue
        x = exact[0]
        rows.append({
            "scott_lambda_nm": r.lambda_nm,
            "scott_lambda_A": round(float(r.lambda_nm) * 10.0, 4),
            "nl_wavelength_A": float(x.wavelength_A),
            "delta_A": round(float(x.wavelength_A) - float(r.lambda_nm) * 10.0, 4),
            "level": f"{_level_key(x.lower_level, x.j_lower)} -> {_level_key(x.upper_level, x.j_upper)}",
            "loggf": float(r.loggf), "elo_eV": float(r.elo_eV),
        })
    if bad:
        raise BuildError("C3 FAILED (the two primaries do not resolve 1:1 on the level key):\n  "
                         + "\n  ".join(bad))
    dmax = max(abs(x["delta_A"]) for x in rows)
    return {"n_matched": len(rows), "n_expected": 7,
            "max_wavelength_disagreement_A": dmax,
            "wavelength_disagreement_note":
                ("MEASURED, and it is why the join is on the LEVEL: the two primaries print the same "
                 "transition at wavelengths up to %.3f A apart, %.0fx the line_match 0.005 A default."
                 % (dmax, dmax / 0.005)),
            "pairs": rows, "status": "PASS"}


def control_c4_level_energy_identity(nl: pd.DataFrame) -> dict:
    """E_up = E_low + hc/lambda_VAC is single-valued per upper level; E_low per lower level.

    Neither paper tabulates E_up, and no published number depends on it, so this cannot be
    satisfied by copying a headline value. A shifted wavelength column, a shifted E_low column,
    or a row that slipped against its level designation all break it.

    🔴 IT IS ALSO THE AIR/VACUUM CONTROL, AND THAT IS NOT A SIDE EFFECT. NL2017's Table A.1 never
    says which medium its wavelengths are on, and the table runs from 2103 A to 12.33 micron --
    right across the region where different sources switch convention. E_up is a VACUUM quantity
    (E = hc/lambda_vac), so reading the column as air and reading it as vacuum give different
    answers, and the levels arbitrate:

        as printed, treated as VACUUM   max E_up spread  0.001672 eV
        as printed, treated as AIR      max E_up spread  0.000369 eV     <- 4.5x tighter

    So the table is AIR throughout, INCLUDING its 2103-2660 A ultraviolet rows, and this is
    measured rather than assumed from a convention. The residual that survives is what E_low's
    3-decimal printing allows and no more.

    ⚠️ WHERE THIS CONTROL IS BLIND, IT SAYS SO. Its sensitivity is dE/dlambda = hc/lambda^2, which
    falls as lambda^-2: at 21208 A a wavelength would have to be wrong by tens of angstroms before
    the identity noticed. `detection_floor_A` reports that per upper level rather than letting a
    single global PASS imply uniform coverage.
    """
    from pipeline.wavelength_util import air_to_vac   # the SSOT converter, never a local copy

    d = nl.copy()
    d["elo"] = _f(d["elo_eV"])
    d["lam_air"] = _f(d["wavelength_A"])
    d["lam_vac"] = air_to_vac(d["lam_air"].to_numpy())
    d["eup"] = d["elo"] + HC_EV_A / d["lam_vac"]
    d["lo_key"] = [_level_key(t, j) for t, j in zip(d.lower_level, d.j_lower)]
    d["up_key"] = [_level_key(t, j) for t, j in zip(d.upper_level, d.j_upper)]

    #: DERIVED, not chosen. E_low is printed to 3 decimals, so each value carries +/-0.0005 eV and
    #: two rows' E_up may legitimately differ by up to 0.001 eV on that account alone. Lambda is
    #: printed to 3 decimals too, worth at most 1.4e-6 eV at the bluest row -- not the binding term.
    tol = 0.001

    multi = {k: g for k, g in d.groupby("up_key") if len(g) > 1}
    if not multi:
        raise BuildError("C4 has no multi-row upper level; it would be vacuous")

    offenders, floors = [], {}
    worst_up = 0.0
    for key, g in multi.items():
        spread = float(np.ptp(g["eup"]))
        worst_up = max(worst_up, spread)
        if spread > tol:
            offenders.append(f"{key}: E_up spread {spread:.6f} eV over "
                             f"{list(np.round(g['lam_air'].to_numpy(), 3))} A")
        # how wrong a wavelength in this group would have to be before the identity noticed
        sens = HC_EV_A / float(g["lam_air"].max()) ** 2      # eV per A, at the least sensitive row
        floors[key] = round((tol - spread) / sens, 4)
    if offenders:
        raise BuildError("C4 FAILED (E_up is not single-valued per upper level):\n  "
                         + "\n  ".join(offenders))

    worst_lo = max((float(np.ptp(g["elo"])) for _, g in d.groupby("lo_key") if len(g) > 1),
                   default=0.0)
    if worst_lo > 1e-9:
        raise BuildError(f"C4 FAILED: E_low is not single-valued per lower level "
                         f"(spread {worst_lo:.6f} eV)")

    #: ⚠️ A CONTROL THAT CANNOT FAIL IS NOT A CONTROL (RYA-853/1080). Two ways of showing this one
    #: can: the AIR-vs-VACUUM discriminant above, and a mutation sized to the group's OWN
    #: sensitivity -- a fixed nudge is not a probe, because 0.5 A is 18x the detection floor at
    #: 3082 A and 1/45th of it at 21208 A.
    as_vac = float(max(np.ptp(g["elo"] + HC_EV_A / g["lam_air"]) for g in multi.values()))
    if as_vac <= tol:
        raise BuildError(f"C4 air/vacuum discriminant is NOT discriminating: reading the column "
                         f"as vacuum also passes ({as_vac:.6f} eV <= {tol})")

    sharp = min(floors, key=floors.get)
    probe = d.copy()
    i = int(probe.index[probe.up_key == sharp][0])
    displacement = 2.0 * floors[sharp]
    probe.loc[i, "lam_vac"] = air_to_vac(np.array([probe.loc[i, "lam_air"] + displacement]))[0]
    probe["eup"] = probe["elo"] + HC_EV_A / probe["lam_vac"]
    mutated = float(np.ptp(probe.loc[probe.up_key == sharp, "eup"]))
    if mutated <= tol:
        raise BuildError(f"C4 MUTATION CONTROL FAILED: displacing {sharp} by "
                         f"{displacement} A did not break the identity")

    blind = max(floors, key=floors.get)
    return {"n_rows": int(len(d)),
            "n_upper_levels_with_2plus_rows": len(multi),
            "medium": "AIR, measured not assumed (see air_reading_measured_negative)",
            "max_Eup_spread_eV": round(worst_up, 8),
            "max_Elow_spread_eV": round(worst_lo, 12),
            "tolerance_eV": tol,
            "tolerance_basis": ("E_low is printed to 3 decimals: +/-0.0005 eV per value, so two "
                                "rows' E_up may differ by up to 0.001 eV. Lambda's 3-decimal "
                                "printing contributes at most 1.4e-6 eV and is not binding."),
            "air_reading_measured_negative": {
                "max_Eup_spread_eV_if_column_read_as_vacuum": round(as_vac, 8),
                "ratio_to_air_reading": round(as_vac / worst_up, 2),
                "note": ("Treating the printed wavelengths as VACUUM inflates the worst spread "
                         "past the tolerance. The two readings are separable and the table is "
                         "AIR -- across the 2103-2660 A UV rows as well.")},
            "detection_floor_A": floors,
            "sharpest_level": {"level": sharp, "floor_A": floors[sharp]},
            "⚠️ blindest_level": {
                "level": blind, "floor_A": floors[blind],
                "note": ("Sensitivity is hc/lambda^2. This control cannot see a wavelength error "
                         "smaller than this in this group, and PASS must not be read as uniform "
                         "coverage across the table.")},
            "mutation_probe": {"level": sharp, "displacement_A": displacement,
                               "spread_after_eV": round(mutated, 6),
                               "note": "sized to the group's own detection floor, not a fixed nudge"},
            "status": "PASS"}


def control_c5_fig8_set(used: pd.DataFrame, nl: pd.DataFrame) -> dict:
    """The derived six ARE the six NL2017's Fig. 8 axis names, and those names are unambiguous.

    Fig. 8 labels each line by its wavelength truncated to the integer angstrom (6696 6698 7835
    8912 10768 10872). Two things have to hold for that axis to be evidence of a SET:

      a) the six derived here truncate to exactly those six labels; and
      b) no OTHER row of NL2017's own Table A.1 truncates to any of them -- otherwise a label
         would not identify a line and the axis could not settle the count. 8912 vs 8923 and
         10872 vs 10891 are the close calls, and they separate.
    """
    got = sorted(int(w) for w in _f(used["wavelength_air_A"]))   # truncate, as the axis does
    want = sorted(NL2017_FIG8_LINES)
    if got != want:
        raise BuildError(f"C5 FAILED: derived set {got} does not match NL2017 Fig. 8 {want}")

    collisions = {}
    used_lam = set(np.round(_f(used["wavelength_air_A"]), 3))
    for w in _f(nl["wavelength_A"]):
        if int(w) in want and round(float(w), 3) not in used_lam:
            collisions.setdefault(int(w), []).append(round(float(w), 3))
    if collisions:
        raise BuildError(f"C5 FAILED: Fig. 8 labels are ambiguous within Table A.1 -- {collisions}")

    return {"derived": got, "fig8_axis": list(want),
            "label_rule": "wavelength truncated to the integer angstrom, as the axis prints it",
            "labels_unique_within_tableA1": True,
            "nearest_non_member_A": {
                "8912 vs": 8923.555, "10872 vs": 10891.732, "7835 vs": 7836.134,
                "note": "the nearest Table A.1 rows that do NOT carry a Fig. 8 label"},
            "status": "PASS"}


# ── build ─────────────────────────────────────────────────────────────────────────────

def build() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    nl, sc, hf = load_sources()

    c1 = control_c1_scott_table1(sc)
    c2 = control_c2_nlte_column(sc)
    c3 = control_c3_cross_paper(nl, sc, hf)
    c4 = control_c4_level_energy_identity(nl)

    by_level = {p["level"]: p for p in c3["pairs"]}
    hfi = hf.set_index("lambda_nm")
    nl_by_lam = {float(r.wavelength_A): r for _, r in nl.iterrows()}

    rows = []
    for _, r in sc.iterrows():
        pair = next(p for p in c3["pairs"] if p["scott_lambda_nm"] == r.lambda_nm)
        x = nl_by_lam[pair["nl_wavelength_A"]]
        h = hfi.loc[r.lambda_nm]
        excluded = abs(float(x.wavelength_A) - 10891.732) < 0.05
        rows.append({
            "line_set": "asplund_agss21_al",
            "species": "Al I",
            "wavelength_air_A": float(x.wavelength_A),
            "wavelength_air_A_scott2015b": pair["scott_lambda_A"],
            "lambda_air_nm_scott2015b": float(r.lambda_nm),
            "wavelength_disagreement_A": pair["delta_A"],
            "elo_eV": float(r.elo_eV),
            "eup_eV": round(float(eup_eV(float(r.elo_eV), float(x.wavelength_A))), 6),
            "lower_level": _level_key(x.lower_level, x.j_lower),
            "upper_level": _level_key(x.upper_level, x.j_upper),
            "j_lower": str(h.j_lower), "j_upper": str(h.j_upper),
            "loggf": float(r.loggf),
            "loggf_sigma_dex": (np.nan if pd.isna(x.loggf_sigma_dex)
                                else float(x.loggf_sigma_dex)),
            "loggf_sigma_basis": ("Kelleher & Podobedova 2008, JPCRD 37, 709 -- 90% uncertainties, "
                                 "per NL2017 Table A.1 note (a)"),
            "gf_source_per_line": _gf_ref(x.gf_ref),
            "gf_source_collective": (
                "Scott et al. 2015b Sect. 5.3: 'The data for our adopted lines come from theoretical "
                "calculations by the OP (Mendoza et al. 1995), under the assumption of LS-coupling.' "
                "NL2017 Table A.1 gives the SAME reference per line (code 2 = TOPbase/Mendoza). "
                "🔴 THIS IS THEORY, NOT LABORATORY: an AGSS21 abundance is not a gf grade (RYA-946), "
                "and the whole AGSS21-lineage Al set rests on Opacity Project LS-coupling values."),
            "ew_pm_scott2015b": float(r.ew_pm),
            "ew_mA_scott2015b": round(float(r.ew_pm) * 10.0, 4),
            "weight_scott2015b": int(float(r.weight)),
            "weight_convention": "1-3, LARGER IS BETTER (Scott Sect. 2 + Sect. 6; verified by C1)",
            "a_lte_3d": float(r.a_lte_3d), "a_lte_mean3d": float(r.a_lte_mean3d),
            "a_lte_hm": float(r.a_lte_hm), "a_lte_marcs": float(r.a_lte_marcs),
            "a_lte_miss": float(r.a_lte_miss),
            "delta_nlte_3d": float(r.delta_nlte_3d),
            "a_nlte_3d": float(r.a_nlte_3d),
            "hfs_a_lower_MHz": float(h.a_lower_MHz), "hfs_b_lower_MHz": float(h.b_lower_MHz),
            "hfs_a_upper_MHz": (np.nan if str(h.a_upper_MHz).strip() in ("", "nan")
                                else float(h.a_upper_MHz)),
            "hfs_b_upper_MHz": (np.nan if str(h.b_upper_MHz).strip() in ("", "nan")
                                else float(h.b_upper_MHz)),
            "hfs_status": ("SOURCE_PUBLISHES_HFS_CONSTANTS; Al is 100% 27Al, I = 5/2 "
                           "(Scott Table 3 header)"),
            "selection_status": ("EXCLUDED_BY_SOURCE_ANALYSIS" if excluded
                                 else "USED_BY_SOURCE_ANALYSIS"),
            "selection_reason": (
                "Nordlander & Lind 2017 Sect. 3.1.5, verbatim: 'Adopting a line selection and "
                "weights from Scott et al. (2015), but disregarding the line at 10 891 A due to "
                "TELLURIC CONTAMINATION'. Published NEGATIVE selection, preserved not dropped "
                "(RYA-946)." if excluded else
                "Retained by Scott et al. 2015b Table 2 and carried into Nordlander & Lind 2017's "
                "solar abundance; named individually on the Fig. 8 axis."),
            "adopted_solar_value_lineage": (
                "AGSS21 (A&A 653, A141) Sect. Aluminium adopts Nordlander & Lind 2017 (A&A 607, "
                "A75) A(Al) = 6.43 +/- 0.03, whose line selection and weights are Scott et al. "
                "2015b (A&A 573, A25) Table 2 minus 10891 A"),
            "source": ("Nordlander & Lind 2017, A&A 607, A75, Table A.1 (line data) + "
                       "Scott et al. 2015b, A&A 573, A25, Tables 2 and 3 (selection, EW, weight, "
                       "level J)"),
            "source_key": "nordlander_lind2017;scott2015b;asplund2021",
            "source_band": "",   # filled below
        })

    lines = pd.DataFrame(rows).sort_values("wavelength_air_A").reset_index(drop=True)

    from scripts.build_al_intake_rya1132 import band  # ONE band vocabulary, not a second one
    lines["source_band"] = [band(w) for w in lines["wavelength_air_A"]]

    used = lines[lines.selection_status.eq("USED_BY_SOURCE_ANALYSIS")]
    c5 = control_c5_fig8_set(used, nl)

    # ── the wider Table A.1, with roles ONLY where the paper states them ───────────────
    an = pd.DataFrame({
        "line_set": "nordlander_lind2017_analysis_lines",
        "species": "Al I",
        "wavelength_air_A": _f(nl["wavelength_A"]),
        "lower_level": [_level_key(t, j) for t, j in zip(nl.lower_level, nl.j_lower)],
        "upper_level": [_level_key(t, j) for t, j in zip(nl.upper_level, nl.j_upper)],
        "elo_eV": _f(nl["elo_eV"]),
        "loggf": _f(nl["loggf"]),
        "loggf_sigma_dex": _f(nl["loggf_sigma_dex"]),
        "gf_source_per_line": [_gf_ref(k) for k in nl["gf_ref"]],
        "gamma6_sigma_au": _f(nl["gamma6_sigma_au"]),
        "gamma6_alpha": _f(nl["gamma6_alpha"]),
    })
    an["eup_eV"] = np.round(eup_eV(an["elo_eV"], an["wavelength_air_A"]), 6)
    an["source_band"] = [band(w) for w in an["wavelength_air_A"]]

    used_lam = set(np.round(used["wavelength_air_A"].to_numpy(float), 3))
    excl_lam = set(np.round(
        lines.loc[lines.selection_status.eq("EXCLUDED_BY_SOURCE_ANALYSIS"),
                  "wavelength_air_A"].to_numpy(float), 3))
    roles, notes = [], []
    for w in an["wavelength_air_A"]:
        k = round(float(w), 3)
        if k in used_lam:
            roles.append("SOLAR_ABUNDANCE_USED")
            notes.append("In the six-line set NL2017 Fig. 8 names.")
        elif k in excl_lam:
            roles.append("SOLAR_ABUNDANCE_EXCLUDED_TELLURIC")
            notes.append("Scott retained it; NL2017 disregarded it for telluric contamination.")
        elif float(w) in NL_NAMED_ROLES:
            r, n = NL_NAMED_ROLES[float(w)]
            roles.append(r)
            notes.append(n)
        else:
            roles.append("ANALYSIS_LINE_ROLE_NOT_STATED_PER_LINE")
            notes.append("Table A.1's note says only that the table lists 'lines used in the "
                         "spectrum analyses'. The paper states no role for THIS row, and RYA-946 "
                         "forbids inferring one from the wavelength. NOT part of the solar "
                         "abundance set.")
    an["role"] = roles
    an["role_basis"] = notes
    an["source"] = "Nordlander & Lind 2017, A&A 607, A75, Table A.1"
    an["source_key"] = "nordlander_lind2017"
    an = an.sort_values(["wavelength_air_A", "upper_level"]).reset_index(drop=True)

    prov = {
        "ticket": "RYA-1173",
        "artifact": "asplund2021_al_lines.csv",
        "what": ("The AGSS21-lineage SOLAR Al reference line set: the seven Al I lines Scott et al. "
                 "2015b retained, six of which carry the A(Al) = 6.43 that AGSS21 adopts, plus the "
                 "seventh as a PRESERVED published exclusion."),
        "generator": "scripts/rya1173_build_asplund_al_lineset.py",
        "🔴 agss21_publishes_no_al_line_list": (
            "Unlike Fe (AGSS21 Table A.2, RYA-1109), AGSS21 publishes no Al table at all. This set "
            "is RECONSTRUCTED from the two papers AGSS21 cites. `SOURCE_LINE_LIST_NOT_PUBLISHED` "
            "does NOT apply: the primaries publish the rows, and they are here."),
        "lineage": [
            {"step": 1, "paper": "Asplund, Amarsi & Grevesse 2021, A&A 653, A141",
             "doi": "10.1051/0004-6361/202140445",
             "contributes": "the adopted solar value A(Al) = 6.43 +/- 0.03, and NOTHING per line"},
            {"step": 2, "paper": "Nordlander & Lind 2017, A&A 607, A75",
             "doi": "10.1051/0004-6361/201730427",
             "contributes": ("the analysis AGSS21 adopts; Table A.1 = per-line level identity, "
                             "E_low, log gf, log gf sigma, VdW; Sect. 3.1.5 = the telluric "
                             "exclusion; Fig. 8 = the six used lines, named")},
            {"step": 3, "paper": "Scott, Grevesse, Asplund et al. 2015b, A&A 573, A25",
             "doi": "10.1051/0004-6361/201424109",
             "contributes": ("the SELECTION (seven lines), the solar EWs, the line weights, the "
                             "five-model LTE abundances and the 3D+NLTE result, and Table 3's "
                             "level J identity")},
        ],
        "source_copies_read": {
            "note": ("Two copies of one paper are not interchangeable (RYA-1110): both were read "
                     "for both primaries, and here they agree row for row."),
            "nordlander_lind2017": [
                "PUBLISHED: https://www.aanda.org/articles/aa/full_html/2017/11/aa30427-17/ "
                "(article + T1 + T2 = Table A.1 + T3 = Table A.2), retrieved 2026-09-03",
                "arXiv:1708.01949v2 PDF text layer",
            ],
            "scott2015b": [
                "PUBLISHED: https://www.aanda.org/articles/aa/full_html/2015/01/aa24109-14/ "
                "(T1 = Table 1, T4 = Table 2, T6 = Table 3), retrieved 2026-09-03",
                "arXiv:1405.0279v2 PDF, local at 'Reference documents/1405.0279v2.pdf'",
            ],
            "asplund2021": ["local PDF 'Reference documents/aa40445-21.pdf'"],
        },
        "🔴 published_line_count_conflict": {
            "scott2015b": 7, "nordlander_lind2017": 6, "agss21_prose": 5,
            "adopted_here": 6,
            "authority": ("the PRIMARY. NL2017 Sect. 3.1.5 removes exactly one line from Scott's "
                          "seven, and Fig. 8's axis names the remaining six individually: "
                          "6696 6698 7835 8912 10768 10872."),
            "why_not_five": ("AGSS21's 'these five Al i lines' reproduces from neither source it "
                             "cites. It is recorded, NOT resolved: see raw/lineage_quotations.md "
                             "for the two candidate explanations, neither of which any of the "
                             "three papers states."),
            "consequence": ("A downstream reader must not inherit 'five'. A replication measuring "
                            "this set measures SIX lines."),
        },
        "🔴 the_gf_are_theory": (
            "Every log gf in this set is Opacity Project / TOPbase (Mendoza et al. 1995) under the "
            "LS-coupling assumption -- theory, not laboratory. Scott says so in Sect. 5.3 and "
            "NL2017's per-line reference code agrees on all seven. So the AGSS21-lineage Al set is "
            "NOT a lab-gf set, and matching a Codex line to it is NOT evidence of a gf grade "
            "(RYA-946: 'an AGSS21 abundance value is not a gf grade')."),
        "selection_note": (
            "10891.732 A is IN this file, flagged EXCLUDED_BY_SOURCE_ANALYSIS with NL2017's "
            "published reason (telluric contamination). Scott's Table 2 EW, weight and abundances "
            "for it are carried too, so the exclusion can be audited rather than believed."),
        "band_scope": (
            "The six used lines span VIS (6696, 6698), red-optical (7835, 8912) and NIR "
            "(10768, 10872) on this repo's band vocabulary; the excluded line is NIR. There is NO "
            "FUV/NUV/near-UV/J/H/K line in the solar abundance set -- and that is a published "
            "selection, not a coverage gap: NL2017 analyses the 3961 A resonance line, the "
            "13123/16750 A IR lines and the 12.33 micron emission line as DIAGNOSTICS, and none of "
            "them enters the abundance. `nordlander_lind_2017_analysis_lines.csv` carries all 55 "
            "Table A.1 rows with that distinction made per row."),
        "not_the_same_as": (
            "data/reference/asplund2021_fe/ is AGSS21's OWN published Fe table. This set is a "
            "RECONSTRUCTION from AGSS21's cited primaries, and its line_set axis value is "
            "`asplund-al`, never `asplund` -- the two are different provenance chains and must "
            "never share a product-key value."),
        "controls": {"C1_scott_table1": c1, "C2_nlte_column": c2, "C3_cross_paper_level_join": c3,
                     "C4_level_energy_identity": c4, "C5_fig8_set": c5},
        "counts": {
            "n_rows": int(len(lines)),
            "n_used": int(len(used)),
            "n_excluded_by_source": int(len(lines) - len(used)),
            "n_analysis_lines_tableA1": int(len(an)),
            "analysis_roles": {k: int(v) for k, v in an["role"].value_counts().items()},
        },
    }
    return lines.round(FLOAT_DECIMALS), an.round(FLOAT_DECIMALS), prov


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="rebuild and compare to the committed artifacts; write nothing")
    a = ap.parse_args()

    lines, an, prov = build()

    if a.check:
        bad = []
        for path, got in ((OUT_LINES, lines), (OUT_ANALYSIS, an)):
            if not path.exists():
                bad.append(f"{path.name}: MISSING")
                continue
            want = pd.read_csv(path)
            if not got.reset_index(drop=True).equals(
                    want.astype(got.dtypes.to_dict(), errors="ignore").reset_index(drop=True)):
                # compare through the CSV round trip, which is what "committed" means
                import io
                buf = io.StringIO()
                got.to_csv(buf, index=False)
                if buf.getvalue() != path.read_text():
                    bad.append(f"{path.name}: DRIFTED from its generator")
        if OUT_PROV.exists():
            if json.loads(OUT_PROV.read_text()) != prov:
                bad.append(f"{OUT_PROV.name}: DRIFTED from its generator")
        else:
            bad.append(f"{OUT_PROV.name}: MISSING")
        if bad:
            print("DRIFT:\n  " + "\n  ".join(bad))
            return 1
        print("committed artifacts reproduce from their generator")
        return 0

    REF.mkdir(parents=True, exist_ok=True)
    lines.to_csv(OUT_LINES, index=False)
    an.to_csv(OUT_ANALYSIS, index=False)
    OUT_PROV.write_text(json.dumps(prov, indent=2, ensure_ascii=False) + "\n")

    c1 = prov['controls']['C1_scott_table1']
    print(f"C1 Scott summary table        PASS  8/8; weighted 3D+NLTE "
          f"{c1['weighted_means_nlte']['recommended']} -> 6.43 "
          f"(LTE misreading reproduces {c1['lte_reading_measured_negative']['n_of_5_reproducing']}/5)")
    print(f"C2 dNLTE column               PASS  max residual "
          f"{prov['controls']['C2_nlte_column']['max_residual_dex']} dex")
    print(f"C3 cross-paper level join     PASS  7/7, max lambda disagreement "
          f"{prov['controls']['C3_cross_paper_level_join']['max_wavelength_disagreement_A']} A")
    c4 = prov['controls']['C4_level_energy_identity']
    print(f"C4 E_up per-level identity    PASS  max spread {c4['max_Eup_spread_eV']} eV vs tol "
          f"{c4['tolerance_eV']}; vacuum misreading "
          f"{c4['air_reading_measured_negative']['max_Eup_spread_eV_if_column_read_as_vacuum']} "
          f"({c4['air_reading_measured_negative']['ratio_to_air_reading']}x) => the column is AIR")
    print(f"C5 NL2017 Fig. 8 set          PASS  {prov['controls']['C5_fig8_set']['derived']}")
    print()
    print(f"wrote {OUT_LINES.relative_to(ROOT)}      {len(lines)} rows "
          f"({prov['counts']['n_used']} used + {prov['counts']['n_excluded_by_source']} excluded)")
    print(f"wrote {OUT_ANALYSIS.relative_to(ROOT)}  {len(an)} rows")
    print(f"wrote {OUT_PROV.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
