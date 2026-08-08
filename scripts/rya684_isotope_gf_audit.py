#!/usr/bin/env python3
"""RYA-684 — isotope-fraction double-application audit.

WHAT THIS DECIDES
=================
Turbospectrum multiplies the number density of every isotope-coded species by
``isotopfrac(Z, A)`` before it forms the line opacity::

    bsyn.f:1350   ntot(k) = ntot(k) * isotopfrac(atom(i), isotope(i))
    bsyn.f:1351   ntt(k)  = ntt(k)  * isotopfrac(atom(i), isotope(i))

``makeabund.f`` states the convention that goes with that, verbatim::

    * isotopfrac(x,0)==1.0; this is used in the case of no isotopes wanted
    * in calculation. It will ensure that the total abundance of the species
    * will be used.(if for example isotopic factors were included in gf-values)

So there are exactly two self-consistent ways to ship an isotope-split line:

  (A) ISOTOPE-CODED, gf carries the FULL (fraction-free) oscillator strength.
      The engine applies the fraction.  Correct.
  (B) NOT isotope-coded (species written ``Z.000``), gf carries the fraction.
      ``isotopfrac(Z,0) == 1`` so the engine applies nothing.  Also correct.

Shipping a list that is isotope-coded AND has the fraction folded into log gf
applies the fraction TWICE.  The effective total gf of the feature is then

    sum_i f_i^2   instead of   sum_i f_i == 1

and a fit that inverts an observed EW against that weakened line has to raise
the abundance by

    dA = -log10( sum_i f_i^2 )        (linear part of the curve of growth)

which for Eu II (f151=0.478, f153=0.522) is +0.3002 dex — the number RYA-565
measured on its two-leg VALD-vs-GES comparison (+0.300).

WHAT THIS SCRIPT DOES
=====================
It does not assume the fold; it MEASURES it, per species, on the shipped files.

For every species that has at least one isotope-coded block, group the lines
into *features* (same species, ion and transition designation).  Within a
feature, sum gf per isotope.  Then reconstruct, per isotope,

    unfolded_i = log10( sum gf_i ) - log10( f_i )

If the shipped gf already carries the fraction, ``unfolded_i`` is the same
physical number for every isotope and the spread collapses.  If the shipped gf
is fraction-free, it is ``log10(sum gf_i)`` that collapses instead.  The two
hypotheses are distinguishable whenever the feature is seen in two isotopes
whose fractions differ, and the verdict is read off the smaller spread.

Both consumer surfaces are audited:

  * ``vald``  — the TSFitPy ``linelist_vald`` "for-grid" lists that every
    ``scripts/rya*_synth_sirius.py`` bsyn harness copies verbatim.  Those
    harnesses emit ``'ISOTOPES : ' '0'``, i.e. no override, so the engine uses
    the ``makeabund.f`` defaults this script parses.
  * ``ges``   — ``GESv6_atom_hfs_iso.420_920nm/atomic_lines.tsv``, the list
    ``pipeline/abundances_derive`` feeds to iSpec.  iSpec writes an explicit
    ``ISOTOPES :`` block from ``input/isotopes/SPECTRUM.lst``
    (``ispec/synth/turbospectrum.py:260``), so that surface's fractions come
    from SPECTRUM.lst, not from makeabund.f.
  * ``ts-nlte-ges`` — ``COM/linelists/nlte_ges_linelist_jmg17feb2022_I_II``,
    the deck the RYA-533/534 Engine-B gate (``scripts/ts_gerber_gate.py``)
    synthesises on.  It also emits ``'ISOTOPES : ' '0'``, so makeabund.f
    defaults apply there too.

Nothing is hardcoded that can drift: the isotope fractions are parsed out of
the engine source and out of SPECTRUM.lst, and the script loud-fails if either
is missing rather than falling back to a stale copy.

Usage (Sirius)::

    /mnt/codex-data/venv_ci/bin/python scripts/rya684_isotope_gf_audit.py \
        --out data/audit/rya684_isotope_gf_audit.json
"""

from __future__ import annotations

import argparse
import collections
import glob
import json
import math
import os
import re
import sys

# ── Engine + line-list locations (Sirius) ────────────────────────────────────
TS_SOURCE_DIR = "/mnt/codex-data/engines/Turbospectrum_NLTE/source"
MAKEABUND     = os.path.join(TS_SOURCE_DIR, "makeabund.f")
BSYN          = os.path.join(TS_SOURCE_DIR, "bsyn.f")
VALD_DIR      = "/mnt/codex-data/engines/TSFitPy/input_files/linelists/linelist_vald"
ISPEC_DIR     = "/mnt/codex-data/engines/ispec_src"
GES_TSV       = os.path.join(ISPEC_DIR, "input", "linelists", "transitions",
                             "GESv6_atom_hfs_iso.420_920nm", "atomic_lines.tsv")
SPECTRUM_LST  = os.path.join(ISPEC_DIR, "input", "isotopes", "SPECTRUM.lst")
TS_NLTE_GES   = os.path.join("/mnt/codex-data/engines/Turbospectrum_NLTE", "COM",
                             "linelists", "nlte_ges_linelist_jmg17feb2022_I_II")

# A species is "molecular" in the Turbospectrum encoding when the species code
# is >= 100 (2-atom codes like 606 = C2, 822 = TiO) — audited separately from
# the atomic sweep because their "isotope" field is an isotopologue pair.
ATOMIC_Z_MAX = 92

# Spread (dex) under which a set of per-isotope reconstructions is called
# consistent.  0.02 dex is well below the smallest thing this project acts on
# and well above VALD/GES 3-decimal log gf rounding.
CONSISTENT_DEX = 0.02

# Per-isotope components closer than this blend into a single observed feature
# at solar resolution, so the whole feature shares one offset.  Atomic isotope
# shifts are a few mA; molecular isotopologue lines are whole angstroms apart.
COLOCATED_MA = 50.0

# Consecutive lines of one designation further apart than this are different
# fine-structure lines, not isotope/HFS components of one feature.  Atomic
# isotope shifts top out around 0.2 A (Li 6707 is the extreme at 0.16 A) and
# HFS patterns span ~0.3 A, so 0.5 A separates the two cases cleanly.
CLUSTER_GAP_A = 0.5

# 5040 / 5777 K — the Boltzmann factor exponent used to weight a line's strength
# by its excitation potential when scoring how much of a fit window is missing.
SOLAR_THETA = 5040.0 / 5777.0

VALD_HDR = re.compile(r"^'\s*(\d+)\.(\d+)\s*'\s+(\d+)\s+(\d+)")


# ── Engine-side isotope fractions ────────────────────────────────────────────

def parse_makeabund_isotopfrac(path: str = MAKEABUND) -> dict:
    """Parse the ``isotopfrac(Z,A)=f`` defaults straight out of makeabund.f."""
    if not os.path.exists(path):
        raise SystemExit(f"RYA-684: Turbospectrum source not found at {path} — the "
                         f"engine-side isotope fractions cannot be sourced. STOP.")
    pat = re.compile(r"^\s*isotopfrac\((\d+),(\d+)\)\s*=\s*([0-9.eEdD+-]+)")
    out: dict = {}
    for ln in open(path):
        if ln.lstrip().startswith(("*", "c", "C", "!")):
            continue
        m = pat.match(ln)
        if not m:
            continue
        Z, A = int(m.group(1)), int(m.group(2))
        f = float(m.group(3).replace("d", "e").replace("D", "e"))
        if A == 0:
            continue          # the "no isotope wanted" sentinel, always 1.0
        out[(Z, A)] = f
    if not out:
        raise SystemExit(f"RYA-684: parsed zero isotopfrac assignments from {path} — "
                         f"the parser and the engine source have drifted. STOP.")
    return out


def parse_spectrum_lst(path: str = SPECTRUM_LST) -> dict:
    """Parse iSpec's SPECTRUM.lst isotope table (Z, A -> solar-system fraction)."""
    if not os.path.exists(path):
        raise SystemExit(f"RYA-684: iSpec isotope table not found at {path}. STOP.")
    out: dict = {}
    for ln in open(path):
        parts = ln.split()
        if len(parts) < 4:
            continue
        try:
            Z = int(float(parts[0])); A = int(parts[1]); f = float(parts[3])
        except ValueError:
            continue
        out[(Z, A)] = f
    if not out:
        raise SystemExit(f"RYA-684: parsed zero isotopes from {path}. STOP.")
    return out


def verify_bsyn_application_site(path: str = BSYN) -> dict:
    """Assert that bsyn.f still multiplies the populations by isotopfrac.

    The whole finding rests on this one statement existing; if a future engine
    bump removes or renames it the audit must fail loudly, not quietly pass.
    """
    if not os.path.exists(path):
        raise SystemExit(f"RYA-684: bsyn.f not found at {path}. STOP.")
    hits = []
    for i, ln in enumerate(open(path), start=1):
        s = ln.strip()
        if s.startswith(("*", "c", "C", "!")):
            continue
        if "isotopfrac(atom(i),isotope(i))" in s.replace(" ", "") and "=" in s:
            hits.append((i, s))
    if not hits:
        raise SystemExit("RYA-684: bsyn.f no longer contains the "
                         "ntot*isotopfrac(atom,isotope) application — the premise of "
                         "this audit has changed. STOP.")
    return {"file": path, "sites": [{"line": i, "code": s} for i, s in hits]}


# ── Line-list readers ────────────────────────────────────────────────────────

def read_ts_lists(paths: list) -> list:
    """Read Turbospectrum-format line lists into flat line records.

    Every consumer surface except the iSpec one is in this format: a
    ``'  Z.AAA  '  ion  nlines`` header, a species-name line, then the lines.
    """
    if not paths:
        raise SystemExit("RYA-684: no Turbospectrum-format line lists to read. STOP.")
    rows = []
    for p in paths:
        base = os.path.basename(p)
        Z = iso = ion = None
        with open(p, errors="replace") as fh:
            for ln in fh:
                m = VALD_HDR.match(ln)
                if m:
                    Z, iso, ion = int(m.group(1)), int(m.group(2)), int(m.group(3))
                    iso_raw = m.group(2)
                    continue
                if ln.startswith("'"):
                    continue          # the species-name line
                if Z is None:
                    continue
                parts = ln.split()
                if len(parts) < 3:
                    continue
                try:
                    wave = float(parts[0]); ep = float(parts[1]); loggf = float(parts[2])
                except ValueError:
                    continue
                q = ln.rfind("'")
                p0 = ln.rfind("'", 0, q)
                label = ln[p0 + 1:q] if (q > 0 and p0 >= 0) else ""
                rows.append(dict(Z=Z, iso=iso, iso_raw=iso_raw, ion=ion, wave=wave,
                                 ep=ep, loggf=loggf, label=label, source=base))
    return rows


def split_molecular_code(Z: int, iso_raw: str) -> list:
    """Decompose a Turbospectrum molecular species code into (atom, isotope) pairs.

    ``getinfospecies.f`` reads the part before the dot as consecutive 2-digit
    atomic numbers and the part after it as consecutive 3-digit mass numbers, so
    ``822.000048`` is 48Ti + all-isotope O and ``606.012012`` is 12C12C.
    """
    zs = str(Z)
    if len(zs) % 2:
        zs = "0" + zs
    atoms = [int(zs[i:i + 2]) for i in range(0, len(zs), 2)]
    raw = iso_raw.rjust(3 * len(atoms), "0")
    isos = [int(raw[i:i + 3]) for i in range(0, 3 * len(atoms), 3)]
    return list(zip(atoms, isos))


def vald_for_grid_paths(vald_dir: str = VALD_DIR) -> list:
    files = sorted(glob.glob(os.path.join(vald_dir, "vald-*-for-grid.list")))
    if not files:
        raise SystemExit(f"RYA-684: no VALD for-grid lists under {vald_dir}. STOP.")
    return files


def read_ges_tsv(path: str = GES_TSV) -> list:
    """Read the GES v6 HFS/ISO list that iSpec is fed."""
    if not os.path.exists(path):
        raise SystemExit(f"RYA-684: GES linelist not found at {path}. STOP.")
    rows = []
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        c = {name: i for i, name in enumerate(head)}
        need = ("turbospectrum_species", "spectrum_synthe_isotope", "ion",
                "wave_A", "lower_state_eV", "loggf", "element", "molecule")
        missing = [n for n in need if n not in c]
        if missing:
            raise SystemExit(f"RYA-684: GES list is missing columns {missing} — "
                             f"schema drift. STOP.")
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < len(head):
                continue
            try:
                Z = int(float(f[c["turbospectrum_species"]]))
                iso = int(f[c["spectrum_synthe_isotope"]])
                ion = int(f[c["ion"]])
                wave = float(f[c["wave_A"]])
                ep = float(f[c["lower_state_eV"]])
                loggf = float(f[c["loggf"]])
            except ValueError:
                continue
            rows.append(dict(Z=Z, iso=iso, iso_raw=str(iso).zfill(3), ion=ion,
                             wave=wave, ep=ep, loggf=loggf,
                             label=f[c["element"]] + "|" + f[c["lower_state_eV"]] +
                                   "|" + f[c["upper_state_eV"]] if "upper_state_eV" in c
                                   else f[c["element"]],
                             molecule=f[c["molecule"]] == "T",
                             source=os.path.basename(path)))
    return rows


# ── The measurement ──────────────────────────────────────────────────────────

def group_features(lines: list) -> dict:
    """Group lines that are the same physical FEATURE across isotopes.

    The transition designation plus excitation potential is not enough on its
    own: several fine-structure lines of one term transition share both (Ca II
    3933 and 3968 carry the same label and the same lower state), and merging
    them makes the per-isotope components look tens of angstroms apart.  So
    within a designation the lines are additionally clustered in wavelength,
    breaking a cluster wherever consecutive lines are further apart than any
    isotope shift or HFS spread can account for.

    Returns {feature_id: [line, ...]}.
    """
    buckets = collections.defaultdict(list)
    for r in lines:
        buckets[(r["Z"], r["ion"], r["label"], round(r["ep"], 3))].append(r)

    out = {}
    for key, rows in buckets.items():
        rows.sort(key=lambda r: r["wave"])
        cluster, prev = [], None
        for r in rows:
            if prev is not None and (r["wave"] - prev) > CLUSTER_GAP_A:
                out[key + (round(cluster[0]["wave"], 3),)] = cluster
                cluster = []
            cluster.append(r)
            prev = r["wave"]
        if cluster:
            out[key + (round(cluster[0]["wave"], 3),)] = cluster
    return out


def audit_species(rows: list, fracs: dict, surface: str) -> list:
    """Per-species verdict: is the shipped log gf fraction-folded or not?"""
    by_species = collections.defaultdict(list)
    for r in rows:
        by_species[(r["Z"], r["ion"])].append(r)

    results = []
    for (Z, ion), lines in sorted(by_species.items()):
        isos_present = sorted({r["iso"] for r in lines if r["iso"] != 0})
        if not isos_present:
            continue                      # no isotope coding => no exposure
        kind = "molecular" if Z > ATOMIC_Z_MAX else "atomic"

        # Fractions the engine will apply to each coded block on this surface.
        # For a molecule the engine multiplies the population once per constituent
        # atom, so the applied factor is the product over the code's (atom, A) pairs
        # with isotopfrac(atom, 0) == 1 for an uncoded constituent.
        raw_for = {r["iso"]: r["iso_raw"] for r in lines if r["iso"] != 0}
        applied: dict = {}
        for A in isos_present:
            if kind == "atomic":
                applied[A] = fracs.get((Z, A))
            else:
                prod = 1.0
                for at, a in split_molecular_code(Z, raw_for[A]):
                    if a == 0:
                        continue
                    f = fracs.get((at, a))
                    if f is None:
                        prod = None
                        break
                    prod *= f
                applied[A] = prod
        unknown = [A for A, f in applied.items() if f is None]

        # Group into features and keep those resolved in >= 2 isotopes whose
        # fractions actually differ (otherwise the two hypotheses coincide).
        feats = collections.defaultdict(lambda: collections.defaultdict(float))
        feat_iso_wave = collections.defaultdict(lambda: collections.defaultdict(list))
        feat_wave, feat_label = {}, {}
        for k, rows in group_features([r for r in lines if r["iso"] != 0]).items():
            for r in rows:
                feats[k][r["iso"]] += 10.0 ** r["loggf"]
                feat_iso_wave[k][r["iso"]].append(r["wave"])
            feat_wave[k] = (min(r["wave"] for r in rows), max(r["wave"] for r in rows))
            feat_label[k] = rows[0]["label"]

        # How far apart are the per-isotope components of one feature?  This is
        # what decides which offset arithmetic applies, so it is measured, not
        # assumed: the per-isotope mean wavelength within a feature, spread.
        seps = []
        for k, per_iso_w in feat_iso_wave.items():
            if len(per_iso_w) < 2:
                continue
            means = [sum(ws) / len(ws) for ws in per_iso_w.values()]
            seps.append((max(means) - min(means)) * 1000.0)
        max_sep_mA = round(max(seps), 2) if seps else None

        folded_spreads, raw_spreads, usable = [], [], []
        for k, per_iso in feats.items():
            if len(per_iso) < 2:
                continue
            fs = [applied.get(A) for A in per_iso]
            if any(f is None or f <= 0 for f in fs):
                continue
            if max(fs) / min(fs) < 1.05:
                continue                  # fractions too close to discriminate
            raw = [math.log10(g) for g in per_iso.values()]
            unf = [math.log10(g) - math.log10(applied[A])
                   for A, g in per_iso.items()]
            raw_spreads.append(max(raw) - min(raw))
            folded_spreads.append(max(unf) - min(unf))
            usable.append(dict(feature=str(feat_label[k])[:70], ep=k[3],
                               wave_lo=round(feat_wave[k][0], 3),
                               wave_hi=round(feat_wave[k][1], 3),
                               isotopes={str(A): round(math.log10(g), 4)
                                         for A, g in sorted(per_iso.items())},
                               spread_if_folded=round(max(unf) - min(unf), 4),
                               spread_if_fraction_free=round(max(raw) - min(raw), 4)))

        if not usable:
            verdict = "UNDECIDABLE"
            detail = ("no feature resolved in two isotopes with distinguishable "
                      "fractions on this surface")
        else:
            med_f = sorted(folded_spreads)[len(folded_spreads) // 2]
            med_r = sorted(raw_spreads)[len(raw_spreads) // 2]
            if med_f <= CONSISTENT_DEX and med_r > CONSISTENT_DEX:
                verdict = "FOLDED"
                detail = ("per-isotope gf reconstructs to one physical gf only after "
                          "dividing by the isotope fraction => the fraction is already "
                          "in log gf, and the engine applies it a second time")
            elif med_r <= CONSISTENT_DEX and med_f > CONSISTENT_DEX:
                verdict = "FRACTION_FREE"
                detail = ("per-isotope gf is equal across isotopes => the fraction is "
                          "not in log gf; the engine applying it is correct")
            else:
                verdict = "AMBIGUOUS"
                detail = (f"median spread folded={med_f:.4f} vs fraction-free="
                          f"{med_r:.4f}; neither hypothesis collapses")

        # Linear-CoG offset the double-application WOULD cost if the shipped gf is
        # folded.  Reported for every species, but it only bites where the verdict
        # is FOLDED — the field name says so, so a FRACTION_FREE row cannot be
        # misread as a correction that is owed.
        #
        # TWO DIFFERENT ARITHMETICS, and using the wrong one misstates the answer
        # by dex:
        #   * CO-LOCATED (atomic isotope/HFS structure — the components of one
        #     feature sit within a few mA and blend into one observed line).  The
        #     feature should carry sum_i f_i * gf == gf and instead carries
        #     sum_i f_i^2 * gf, so the offset is -log10(sum f_i^2), ONE number for
        #     the whole feature.
        #   * SEPARATED (molecular isotopologues — 12CH and 13CH are different
        #     lines at different wavelengths).  Each line should carry f_i * gf_i
        #     and instead carries f_i^2 * gf_i, so the offset is -log10(f_i),
        #     PER ISOTOPOLOGUE, and for a trace isotopologue it is enormous while
        #     the co-located sum would have called it negligible.
        known = [f for f in applied.values() if f]
        sum_f2 = sum(f * f for f in known)
        sum_f = sum(known)
        colocated = (max_sep_mA is not None and max_sep_mA <= COLOCATED_MA)
        pred = round(-math.log10(sum_f2), 4) if (colocated and sum_f2 > 0) else None
        per_iso_offset = {str(A): (round(-math.log10(f), 4) if f else None)
                          for A, f in sorted(applied.items())}

        results.append(dict(
            surface=surface, Z=Z, ion=ion, kind=kind,
            isotopes_coded=isos_present,
            fractions_applied={str(A): applied[A] for A in isos_present},
            fractions_unknown_to_engine=unknown,
            coded_line_count=sum(1 for r in lines if r["iso"] != 0),
            uncoded_line_count=sum(1 for r in lines if r["iso"] == 0),
            coverage_sum_f=round(sum_f, 6),
            sum_f_squared=round(sum_f2, 6),
            components_max_separation_mA=max_sep_mA,
            components_colocated=colocated,
            offset_if_folded_dex=pred,
            offset_if_folded_dex_per_isotope=per_iso_offset,
            verdict=verdict, detail=detail,
            features_tested=len(usable),
            evidence=sorted(usable, key=lambda d: d["wave_lo"])[:6],
        ))
    return results


# ── Blast radius: which fit windows actually see an exposed line ─────────────
#
# Every bsyn harness copies the in-window VALD blocks VERBATIM as its blend
# model, so an exposed block does not have to belong to the target species to
# bias it — a blend line synthesised ~0.24-0.30 dex too weak leaves absorption
# for the target to absorb.  These are the fit windows (centre +- fit_hw) of
# every committed harness, read off their LINES tables.
HARNESS_WINDOWS = [
    ("RYA-551", "Sr II", 4077.709, 1.6), ("RYA-551", "Sr II", 4161.792, 0.9),
    ("RYA-551", "Sr II", 4215.519, 1.6), ("RYA-551", "Sr II", 4305.443, 0.9),
    ("RYA-560/585", "Zr II", 4208.980, 1.6), ("RYA-560/585", "Zr II", 4258.041, 1.2),
    ("RYA-560/585", "Zr II", 4442.992, 1.4), ("RYA-560/585", "Zr II", 4629.079, 0.9),
    ("RYA-560/585", "Zr II", 5350.089, 0.9), ("RYA-560/585", "Zr II", 5372.466, 0.9),
    ("RYA-564", "Co I", 4813.476, 0.9), ("RYA-564", "Co I", 5212.688, 0.9),
    ("RYA-564", "Co I", 5352.040, 0.9), ("RYA-564", "Co I", 5647.234, 0.8),
    ("RYA-564", "Co I", 6454.994, 0.9), ("RYA-564", "Co I", 6632.439, 0.8),
    ("RYA-564", "Co I", 6771.034, 0.9), ("RYA-564", "Co I", 5301.041, 0.9),
    ("RYA-564", "Co I", 5331.453, 0.9), ("RYA-564", "Co I", 5483.354, 1.0),
    ("RYA-564", "Co I", 5915.550, 0.8), ("RYA-564", "Co I", 6188.997, 0.8),
    ("RYA-592", "Mg I", 5528.405, 1.4), ("RYA-592", "Mg I", 5711.088, 0.6),
    ("RYA-581", "Ba II", 5853.668, 0.45),
    ("RYA-565", "Eu II", 6645.064, 0.45), ("RYA-565", "Eu II", 6437.640, 0.45),
    ("RYA-565", "Eu II", 6303.410, 0.45), ("RYA-565", "Eu II", 6173.030, 0.45),
    ("RYA-565", "Eu II", 5818.760, 0.45),
]


def audit_windows(rows: list, exposed: dict) -> list:
    """List every harness fit window that contains an exposed isotope-coded line.

    Broken out per isotope, because the aggregate hides the thing that decides
    materiality: a 12CH blend is 0.005 dex too weak (nothing) while a 13CH blend
    is 1.96 dex too weak (deleted), and both are "CH lines in the window".
    """
    out = []
    for ticket, target, centre, hw in HARNESS_WINDOWS:
        hits = [r for r in rows
                if r["iso"] != 0 and (r["Z"], r["ion"]) in exposed
                and abs(r["wave"] - centre) <= hw]
        if not hits:
            continue
        per = collections.defaultdict(list)
        for r in hits:
            per[(r["Z"], r["ion"], r["iso"])].append(r)
        # How much of the window's absorption is actually missing?  A line whose
        # gf is double-folded contributes f^2*gf where it should contribute f*gf,
        # so it is short by a factor (1/f - 1) of what it currently contributes.
        # Weighting by the LTE line-strength proxy gf * 10^(-theta * EP) keeps a
        # deep-lying trace line from being scored as if it were a strong one.
        def strength(r):
            return 10.0 ** (r["loggf"] - SOLAR_THETA * r["ep"])

        in_window = [r for r in rows if abs(r["wave"] - centre) <= hw]
        total = sum(strength(r) for r in in_window)
        missing = 0.0
        for r in hits:
            rec = exposed[(r["Z"], r["ion"])]
            f = rec["fractions_applied"].get(str(r["iso"]))
            if f:
                missing += strength(r) * (1.0 / f - 1.0)

        contam = []
        for (Z, ion, A), v in sorted(per.items()):
            name = f"{ELEM.get(Z, Z)} {ROMAN.get(ion, ion)}"
            rec = exposed[(Z, ion)]
            off = (rec["offset_if_folded_dex"] if rec["components_colocated"]
                   else rec["offset_if_folded_dex_per_isotope"].get(str(A)))
            contam.append(dict(
                species=name, isotope=A,
                is_the_target_species=(name == target),
                n_lines=len(v),
                wave_lo=round(min(x["wave"] for x in v), 3),
                wave_hi=round(max(x["wave"] for x in v), 3),
                line_too_weak_by_dex=off))
        foreign = sum(strength(r) * (1.0 / exposed[(r["Z"], r["ion"])]
                                     ["fractions_applied"][str(r["iso"])] - 1.0)
                      for r in hits
                      if f"{ELEM.get(r['Z'], r['Z'])} {ROMAN.get(r['ion'], r['ion'])}" != target
                      and exposed[(r["Z"], r["ion"])]["fractions_applied"].get(str(r["iso"])))
        out.append(dict(ticket=ticket, target=target, centre=centre,
                        fit_halfwidth=hw, contaminants=contam,
                        window_strength_missing_frac=round(missing / total, 6) if total else None,
                        foreign_blend_strength_missing_frac=(
                            round(foreign / total, 6) if total else None)))
    return out


ELEM = {3: "Li", 20: "Ca", 25: "Mn", 27: "Co", 29: "Cu", 38: "Sr", 40: "Zr",
        56: "Ba", 57: "La", 59: "Pr", 60: "Nd", 62: "Sm", 63: "Eu",
        106: "CH", 108: "OH", 112: "MgH", 114: "SiH",
        606: "C2", 607: "CN", 608: "CO", 822: "TiO"}
ROMAN = {1: "I", 2: "II", 3: "III"}


def render_table(results: list) -> str:
    lines = ["| species | surface | isotopes coded | shipped log gf | exposed | dex if folded |",
             "|---|---|---|---|---|---|"]
    for r in sorted(results, key=lambda d: (d["kind"], d["Z"], d["ion"], d["surface"])):
        name = f"{ELEM.get(r['Z'], 'Z' + str(r['Z']))} {ROMAN.get(r['ion'], r['ion'])}"
        exposed = {"FOLDED": "YES", "FRACTION_FREE": "no"}.get(r["verdict"], "undecided")
        shipped = {"FOLDED": "fraction folded in",
                   "FRACTION_FREE": "fraction-free (full gf)"}.get(r["verdict"], r["verdict"])
        if r["components_colocated"]:
            pred = (f"+{r['offset_if_folded_dex']:.4f} (blended feature)"
                    if r["offset_if_folded_dex"] is not None else "—")
        else:
            per = r["offset_if_folded_dex_per_isotope"]
            vals = [abs(v) for v in per.values() if v is not None]
            pred = (f"+{min(vals):.4f}..+{max(vals):.4f} per isotopologue"
                    if vals else "—")
        lines.append(f"| {name} | {r['surface']} | "
                     f"{','.join(str(a) for a in r['isotopes_coded'])} | {shipped} | "
                     f"{exposed} | {pred} |")
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="RYA-684 isotope-fraction double-application audit")
    ap.add_argument("--out", default=None, help="write the full JSON record here")
    ap.add_argument("--vald-dir", default=VALD_DIR)
    ap.add_argument("--ges", default=GES_TSV)
    ap.add_argument("--ts-nlte-ges", default=TS_NLTE_GES)
    args = ap.parse_args(argv)

    site = verify_bsyn_application_site()
    makeabund = parse_makeabund_isotopfrac()
    spectrum = parse_spectrum_lst()
    print(f"[engine] isotopfrac applied at {site['file']}:"
          f"{','.join(str(s['line']) for s in site['sites'])}")
    print(f"[engine] makeabund.f defaults: {len(makeabund)} isotopes; "
          f"iSpec SPECTRUM.lst: {len(spectrum)} isotopes")

    results = []
    print(f"\n[vald] reading {args.vald_dir}")
    vald_rows = read_ts_lists(vald_for_grid_paths(args.vald_dir))
    results += audit_species(vald_rows, makeabund, "vald(bsyn harnesses)")
    print(f"[tsng] reading {args.ts_nlte_ges}")
    results += audit_species(read_ts_lists([args.ts_nlte_ges]), makeabund,
                             "ts-nlte-ges(Engine-B)")
    print(f"[ges]  reading {args.ges}")
    results += audit_species(read_ges_tsv(args.ges), spectrum, "ges(iSpec path)")

    print()
    print(render_table(results))
    print()
    for r in results:
        if r["verdict"] != "FOLDED":
            continue
        if r["components_colocated"]:
            mag = f"+{r['offset_if_folded_dex']:.4f} dex on the blended feature"
        else:
            per = r["offset_if_folded_dex_per_isotope"]
            mag = ("per isotopologue " +
                   ", ".join(f"{A}:+{v:.4f}" for A, v in per.items() if v is not None))
        print(f"  ! {ELEM.get(r['Z'], r['Z'])} {ROMAN.get(r['ion'], r['ion'])} "
              f"[{r['surface']}] DOUBLE-APPLIED — {mag} "
              f"({r['features_tested']} features tested, components "
              f"{r['components_max_separation_mA']} mA apart)")

    # Which harness fit windows actually see one of those exposed blocks?
    exposed = {(r["Z"], r["ion"]): r for r in results
               if r["verdict"] == "FOLDED" and r["surface"] == "vald(bsyn harnesses)"}
    windows = audit_windows(vald_rows, exposed)
    print("\n[windows] harness fit windows containing an exposed isotope-coded line:")
    if not windows:
        print("  (none)")
    for w in windows:
        own_hit = any(c["is_the_target_species"] for c in w["contaminants"])
        sp = ", ".join(f"{c['species']}({c['isotope']})x{c['n_lines']}"
                       for c in w["contaminants"])
        print(f"  {w['ticket']:12s} {w['target']:6s} {w['centre']:9.3f} "
              f"+-{w['fit_halfwidth']:.2f}  "
              f"{'TARGET EXPOSED' if own_hit else 'blends only  '}  "
              f"window absorption missing: {w['window_strength_missing_frac']:.5%} "
              f"(foreign blends alone {w['foreign_blend_strength_missing_frac']:.5%})")
        print(f"                 {sp}")

    record = dict(
        ticket="RYA-684",
        engine_application_site=site,
        makeabund_isotopfrac_count=len(makeabund),
        spectrum_lst_isotope_count=len(spectrum),
        consistent_dex=CONSISTENT_DEX,
        results=results,
        harness_windows_hit=windows,
    )
    if args.out:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        with open(args.out, "w") as fh:
            json.dump(record, fh, indent=2, sort_keys=True)
            fh.write("\n")
        print(f"\n[out] {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
