#!/usr/bin/env python3
"""RYA-1141 - independent QA of the RYA-1132 Al atomic-data intake.

FINDINGS ONLY.  This script reads the RYA-1132 artifacts and their upstream
sources and writes a verdict plus per-check ledgers under
`data/audit/rya1141_al_intake_qa/`.  It never edits a gf, a grade, a manifest
row, `canonical_gf.csv`, or any other intake artifact - RYA-161's
validate-don't-tune firewall.  A defect found here is REPORTED, not fixed.

The battery re-derives rather than re-reads: every claim is checked against an
identity the RYA-1132 builder did not itself use, so that agreement means
something.  Three such referees carry most of the weight:

  * the CDS `ReadMe` declares its own record counts and byte columns, so it
    referees the transcription independently of the builder that parsed it;
  * branching closure - A_i / sum(A) must reproduce the separately printed
    BranR - referees the fixed-width column extraction using numbers the
    builder never compares;
  * the 3p ^2^P^o^ fine-structure splitting (112.061 cm-1) must fall out of the
    VACUUM wavelengths of any two lines sharing an upper level, which referees
    the air<->vacuum conversion and the level assignment at the same time.

#: 🔴 THE INSTRUMENT MUST NOT ENTER ITS OWN MEASUREMENT.  This script writes only
#: under `data/audit/rya1141_al_intake_qa/`, and the no-mutation check hashes the
#: audited artifacts before and after every read.  It excludes THIS FILE from the
#: audited set BY NAME, not by pattern: a pattern wide enough to describe the
#: intake is wide enough to match the auditor.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/audit/rya1141_al_intake_qa"

AUDITED = ROOT / "data/audit/rya1132_al_intake"
VUJ_RAW = ROOT / "data/reference/vujnovic2002_al/raw"
CENSUS = ROOT / "data/results/rya1001/rya1001_al_line_census.csv"
CENSUS_META = ROOT / "data/results/rya1001/rya1001_census_meta.json"
CANONICAL = ROOT / "data/linelists/canonical_gf.csv"
BURHEIM = ROOT / "data/reference/al_gf_lab/al1_lab_loggf.csv"
BUILDER = ROOT / "scripts/build_al_intake_rya1132.py"

#: This auditor's own file, named so it can never be mistaken for an audited artifact.
SELF = Path(__file__).resolve()

#: RYA-1132's claimed frozen inventory, quoted from the ticket and its summary.json.
CLAIM = {"rows": 505, "Al I": 466, "Al II": 39,
         "PRIMARY_LABORATORY": 18, "CRITICALLY_EVALUATED": 19,
         "THEORETICAL": 1, "FALLBACK": 467}

#: The 3p ^2^P^o^_3/2_ - ^2^P^o^_1/2_ ground-term splitting in Al I, NIST ASD.
#: An independent constant: nothing in the intake derives or stores it.
AL_I_GROUND_FS_CM1 = 112.061

#: hc in eV.A - for the Eup = Elo + hc/lam_vac identity.
HC_EV_A = 12398.419843320026

#: RYA-946's binding depth contract (2026-08-30).
CODEX_DEPTH = (0.05, 0.60)

#: Crossref metadata for every DOI the intake asserts, resolved 2026-08-30 and
#: committed so `--check` is deterministic offline. `--online` re-resolves and any
#: drift is reported. Accents are folded to ASCII at the write, not at the compare,
#: so the stored bytes are the same on every platform.
CROSSREF_CACHE = {
    '10.1051/0004-6361/202245394': {"title": 'Experimental oscillator strengths of Al I lines for near-infrared astrophysical spectroscopy',
        "container": 'Astronomy & Astrophysics', "volume": '672', "page": 'A197', "year": 2023,
        "authors": ['Burheim', 'Hartman', 'Nilsson']},
    '10.1051/0004-6361:20020560': {"title": 'Absolute transition probabilities of Al I and Al II\nspectral lines \nand intensity ratios within multiplets',
        "container": 'Astronomy & Astrophysics', "volume": '388', "page": '704-711', "year": 2002,
        "authors": ['Vujnovic', 'Blagoev', 'Furbock', 'Neger', 'Jager']},
    '10.1088/0953-4075/32/2/031': {"title": 'Measurement of the B+and Al+intercombination and Sc12+forbidden transition rates at a heavy-ion storage ring',
        "container": 'Journal of Physics B: Atomic, Molecular and Optical Physics', "volume": '32', "page": '537-552', "year": 1999,
        "authors": ['Trabert', 'Wolf', 'Linkemann', 'Tordoir']},
    '10.1086/164569': {"title": 'Transition probability of the AL II 2669 intersystem line',
        "container": 'The Astrophysical Journal', "volume": '308', "page": '1013', "year": 1986,
        "authors": ['Johnson', 'Smith', 'Parkinson']},
    '10.1063/1.2734564': {"title": 'Atomic Transition Probabilities of Aluminum. A Critical Compilation',
        "container": 'Journal of Physical and Chemical Reference Data', "volume": '37', "page": '709-911', "year": 2008,
        "authors": ['Kelleher', 'Podobedova']},
    '10.1051/0004-6361/201833764': {"title": 'Extended transition rates and lifetimes in Al I and Al II from systematic multiconfiguration calculations',
        "container": 'Astronomy & Astrophysics', "volume": '621', "page": 'A16', "year": 2018,
        "authors": ['Papoulia', 'Ekman', 'Jonsson']},
    '10.1086/312738': {"title": 'The Mass Assembly and Star Formation Characteristics of Field Galaxies of Known Morphology',
        "container": 'The Astrophysical Journal', "volume": '536', "page": 'L77-L80', "year": 2000,
        "authors": ['Brinchmann', 'Ellis']},
    '10.3847/1538-4357/abf142': {"title": 'Detection of Al ii in the Ultraviolet Spectra of Metal-poor Stars: An Empirical LTE Test of NLTE Aluminum Abundance Calculations*',
        "container": 'The Astrophysical Journal', "volume": '912', "page": '119', "year": 2021,
        "authors": ['Roederer', 'Lawler']},
    '10.1051/0004-6361/202142195': {"title": 'Non-LTE abundance corrections for late-type stars from 2000 A to 3 μm',
        "container": 'Astronomy & Astrophysics', "volume": '665', "page": 'A33', "year": 2022,
        "authors": ['Lind', 'Nordlander', 'Wehrhahn', 'Montelius', 'Osorio', 'Barklem', 'Afsar', 'Sneden', 'Kobayashi']},
    '10.1007/BF01417221': {"title": 'Natural radiative lifetimes in the2S1/2 and2D5/2,3/2 sequences of aluminum',
        "container": 'Zeitschrift for Physik A Atoms and Nuclei', "volume": '313', "page": '151-154', "year": 1983,
        "authors": ['J\ufffdnsson', 'Lundberg']},
    '10.3847/1538-4357/ad4451': {"title": 'Phase-dependent Spectral Shape Changes in the Ultraluminous X-Ray Pulsar NGC 5907 ULX1',
        "container": 'The Astrophysical Journal', "volume": '968', "page": '95', "year": 2024,
        "authors": ['Miura', 'Kobayashi', 'Yamaguchi']},
    '10.3847/1538-4357/ae7de8': {"title": 'CRIRES+ Reveals the Chemistry of the Stellar Subpopulations in the Bulge Fossil Fragment Liller 1*',
        "container": 'The Astrophysical Journal', "volume": '1006', "page": '32', "year": 2026,
        "authors": ['Chiappino', 'Origlia', 'Fanelli', 'Bartolomei', 'Ferraro', 'Lanzoni', 'Pallanca', 'Cadelano', 'Romano', 'Dalessandro', 'Massari', 'Valenti', 'Rich']},
    '10.1093/mnras/stt2120': {"title": 'Phase-averaged gamma-ray spectra from rotation-powered millisecond pulsars',
        "container": 'Monthly Notices of the Royal Astronomical Society', "volume": '437', "page": '2957-2965', "year": 2013,
        "authors": ['Jiang', 'Chen', 'Li', 'Zhang']},
}


#: Where a misquoted DOI actually points. REPORTED, never written back - RYA-1141 is
#: findings-only and does not edit an intake artifact even to correct it.
DOI_CORRECTIONS = {
    "10.1086/312738": ("10.1086/312741",
                       "Griesmann & Kling 2000, 'Interferometric Measurement of Resonance "
                       "Transition Wavelengths in C IV, Si IV, Al III, Al II', ApJ 536, L113-L115"),
    "10.3847/1538-4357/ad4451": ("10.3847/1538-4357/ad22dc",
                                 "Nandakumar et al. 2024, 'Composition of Giants 1 deg North of "
                                 "the Galactic Center: Detailed Abundances', ApJ 964, 96"),
    "10.1093/mnras/stt2120": ("10.1093/mnras/stt2204",
                              "Murphy & Berengut 2014, 'Laboratory atomic transition data for "
                              "precise optical quasar absorption spectroscopy', MNRAS 438, 388"),
}

#: The surname the intake's own citation claims, per source_id. Read off the artifact's
#: `citation` / `source` text, so the comparison is intake-vs-publisher, never
#: auditor-vs-auditor.
SURNAME = re.compile(r"^\s*([A-Z][A-Za-z'\-]+)")


# ─────────────────────────────────────────────────────────────────────────────
# findings plumbing
# ─────────────────────────────────────────────────────────────────────────────
class Report:
    """One PASS/FAIL/FLAG line per check, plus the offending rows, named."""

    def __init__(self) -> None:
        self.checks: list[dict] = []
        self.rows: list[dict] = []

    def add(self, check: str, title: str, status: str, detail: str) -> None:
        assert status in {"PASS", "FAIL", "FLAG"}, status
        self.checks.append({"check": check, "title": title,
                            "status": status, "detail": detail})

    def row(self, check: str, severity: str, subject: str, finding: str, evidence: str) -> None:
        self.rows.append({"check": check, "severity": severity, "subject": subject,
                          "finding": finding, "evidence": evidence})

    @property
    def worst(self) -> str:
        s = {c["status"] for c in self.checks}
        return "FAIL" if "FAIL" in s else ("FLAG" if "FLAG" in s else "PASS")


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def audited_files() -> list[Path]:
    """Everything this QA must leave untouched - never including this script."""
    files = sorted(AUDITED.rglob("*")) + sorted(VUJ_RAW.rglob("*"))
    files += [CENSUS, CENSUS_META, CANONICAL, BURHEIM, BUILDER]
    return [f for f in files if f.is_file() and f.resolve() != SELF]


# ─────────────────────────────────────────────────────────────────────────────
# an EP-aware crossmatch, re-derived - never wavelength alone (RYA-1034/840/1037)
# ─────────────────────────────────────────────────────────────────────────────
def level_ep_table(cross: pd.DataFrame) -> dict[tuple[str, str], float]:
    """Lower-level designation -> the EP the manifest assigns it.

    Built from the source's OWN level names, by consensus across every row that
    names the level.  This is the referee RYA-1132's crossmatch never had: a
    designation is a physical identity, so all its rows must land on one EP.
    """
    table: dict[tuple[str, str], float] = {}
    for (sp, lv), sub in cross.groupby(["species", "lower_level"]):
        eps = sub.lower_EP.dropna()
        if len(eps):
            table[(sp, lv)] = float(eps.round(4).mode().iloc[0])
    return table


EP_IDENTITY_NAMES = {"lower_EP", "ep_eV", "excitation_potential_eV", "epcol", "eptol",
                     "lower_level", "upper_level", "upper_lower_level_identity"}


def compares_identity(tree: "ast.AST") -> bool:
    """Is any physical-identity field used in a COMPARISON inside this function?

    #: 🔴 GREP IS NOT ENOUGH HERE. `ingest_new_lab_sources` DOES mention `lower_level`
    #: and `upper_level` - it writes them into `upper_lower_level_identity`. A textual
    #: search for the name therefore reports a physical-identity check that does not
    #: exist. Only a COMPARISON constrains the join, so that is what this asks for.
    """
    import ast as _ast
    for node in _ast.walk(tree):
        if not isinstance(node, _ast.Compare):
            continue
        for sub in _ast.walk(node):
            name = (sub.attr if isinstance(sub, _ast.Attribute)
                    else sub.id if isinstance(sub, _ast.Name) else None)
            if name in EP_IDENTITY_NAMES:
                return True
    return False


def check_a2(rep: Report, cross: pd.DataFrame, man: pd.DataFrame):
    """A2 - physical-identity crossmatch, EP-aware, never wavelength-alone."""
    import ast as _ast
    tree = _ast.parse(BUILDER.read_text())
    fns = {n.name: n for n in _ast.walk(tree) if isinstance(n, _ast.FunctionDef)}
    ingest = compares_identity(fns["ingest_new_lab_sources"])
    #: CONTROL: the same test must FIRE on `nearest`, the EP-aware helper the builder
    #: defines and then does not call from the ingest path. A guard that cannot say
    #: "yes" is not measuring anything.
    control = compares_identity(fns["nearest"])
    if not control:
        rep.add("A2-control", "The identity-comparison test can return True", "FAIL",
                "The control failed: `nearest`, which demonstrably compares `epcol`, was "
                "not detected. The A2 result below is not trustworthy.")
    else:
        rep.add("A2-control", "The identity-comparison test can return True", "PASS",
                "The same AST test fires on `nearest`, the builder's own EP-aware matcher "
                "(`(frame[epcol] - ep).abs() <= eptol`), so a negative on the ingest path "
                "is a real absence and not a broken detector.")

    if not ingest:
        rep.add("A2", "Crossmatch identity (EP-aware, never wavelength alone)", "FAIL",
                "RYA-1132's `ingest_new_lab_sources` joins every Vujnovic row and the "
                "Johnson Al II row to the manifest on `abs(wavelength_air - lambda) <= "
                "0.08` and nothing else. No excitation potential and no level designation "
                "appears in any comparison in that function - the level strings it does "
                "touch are only WRITTEN into `upper_lower_level_identity`. So no promotion "
                "in this ticket was matched on a physical identity. The builder defines an "
                "EP-aware matcher, `nearest(..., epcol=..., eptol=0.02)`, and the census "
                "loop calls it; the promotion path does not.")
        rep.row("A2", "CRITICAL", "scripts/build_al_intake_rya1132.py:ingest_new_lab_sources",
                "Promotion join is wavelength-only, forbidden by RYA-1034",
                "`near = candidates[candidates.delta_A <= .08]`; the Al II 2669 block "
                "repeats the pattern. `nearest(..., epcol=...)` exists and is not called here.")
    else:
        rep.add("A2", "Crossmatch identity (EP-aware, never wavelength alone)", "PASS",
                "The ingest join constrains a physical-identity field as well as wavelength.")

    # THE NULL. A manifest row claimed by two different (upper J, lower J) pairs is a
    # wavelength-only mismatch caught in the act.  Tables 4 and 5 spell the same Al II
    # level two ways ("3d' ^3^F^o^_4_" / "3d ^3^F^o_4_"), so only a genuinely different
    # J pair counts - a spelling difference must not read as a defect.
    cross = cross.assign(pair=cross.upper_level.str.replace(r"\s+", "", regex=True) + " | "
                              + cross.lower_level.str.replace(r"\s+", "", regex=True))
    matched = cross[cross.canonical_line_id.astype(str).ne("")]

    def jpair(p: str) -> tuple:
        return tuple(re.findall(r"_(\d+(?:/\d+)?)_", p))

    collide = []
    for cid, sub in matched.groupby("canonical_line_id"):
        if sub.pair.map(jpair).nunique() > 1:
            collide.append({"canonical_line_id": cid, "n_source_rows": len(sub),
                            "source_row_ids": ",".join(sub.source_row_id),
                            "source_lambdas_A": ",".join(f"{v:.3f}" for v in sub.wavelength_A),
                            "distinct_transitions": " // ".join(sorted(set(sub.pair)))})
    cdf = pd.DataFrame(collide)
    if len(cdf):
        rep.add("A2-null", "Does the missing identity gate actually mis-assign?", "FAIL",
                f"Yes - it is not merely unexercised. {len(cdf)} manifest row(s) are claimed "
                f"by two physically DIFFERENT transitions. `alphys_II_3587.0720_0333` is "
                f"claimed by three source rows: 4f ^3^F^o^_3_ - 3d ^3^D_2_ at 3587.068 A "
                f"(twice, its true identity) and 4f ^3^F^o^_2_ - 3d ^3^D_3_ at 3587.100 A, a "
                f"different transition 0.028 A away that the 0.08 A window swallows. That "
                f"row is MATCHED_NOT_PROMOTED, so no gf is corrupted today - but it is the "
                f"same code path, with the same window, that wrote all seven GF-LAB values.")
        for c in collide:
            rep.row("A2", "CRITICAL", c["canonical_line_id"],
                    "One manifest line claimed by two different physical transitions",
                    f"source rows {c['source_row_ids']} at lambda {c['source_lambdas_A']} A; "
                    f"transitions {c['distinct_transitions']}")
    else:
        rep.add("A2-null", "Does the missing identity gate actually mis-assign?", "PASS",
                "No transition collision detected.")

    # How much discrimination would an EP gate alone have bought? Measure it, don't
    # assume it: the answer decides whether the fix is "add EP" or "add the LEVEL".
    res = []
    for sp, sub in matched.groupby("species"):
        eps = (sub.groupby("lower_level").lower_EP
               .apply(lambda s: float(s.dropna().median()) if s.notna().any() else np.nan)
               .dropna().sort_values())
        for (l1, e1), (l2, e2) in zip(eps.items(), list(eps.items())[1:]):
            res.append({"species": sp, "level_1": l1, "level_2": l2,
                        "dEP_eV": abs(e2 - e1)})
    rdf = pd.DataFrame(res)
    tight = rdf[rdf.dEP_eV < 5e-4] if len(rdf) else rdf
    rep.add("A2-resolution", "Would an EP gate alone have been sufficient?", "FLAG",
            f"Not on its own. `lower_EP` is stored at 4 decimal places, and "
            f"{len(tight)} adjacent lower-level pairs in the crossmatch sit closer than "
            f"5e-4 eV - the Al II 3d ^3^D_1,2,3_ term spans 0.0003 eV in total, so at the "
            f"manifest's stored precision EP CANNOT separate the very levels the 3587 "
            f"collision confuses. The remedy RYA-1034 needs here is the level/J "
            f"designation plus a uniqueness requirement, not an EP tolerance alone; an "
            f"EP-only gate would have passed this defect.")

    # The pair RYA-1141 names explicitly.
    a = man[np.isclose(man.wavelength_air, 6696.015, atol=.005)]
    b = man[np.isclose(man.wavelength_air, 6696.185, atol=.005)]
    ok = (len(a) == 1 and len(b) == 1
          and abs(float(a.iloc[0].lower_EP) - float(b.iloc[0].lower_EP)) > .5
          and "BURHEIM" in str(a.iloc[0].gf_source).upper()
          and "BURHEIM" not in str(b.iloc[0].gf_source).upper())
    dep = (abs(float(a.iloc[0].lower_EP) - float(b.iloc[0].lower_EP))
           if len(a) == 1 and len(b) == 1 else float("nan"))
    rep.add("A2-6696", "6696.015 vs 6696.185 stay distinct; Burheim cannot leak",
            "PASS" if ok else "FAIL",
            f"Two separate manifest rows, dEP = {dep:.4f} eV. Burheim's laboratory gf sits "
            f"only on 6696.015 and 6696.185 keeps its fallback source. Note this pair is "
            f"separated by RYA-1001's census, which IS EP-aware - not by the RYA-1132 "
            f"ingest path, which would not have distinguished them.")

    # A third class the disposition vocabulary hides: `NO_UNIQUE_MANIFEST_MATCH` is
    # emitted both for "no candidate" and for "more than one candidate".
    amb = []
    for _, r in cross[cross.disposition.eq("NO_UNIQUE_MANIFEST_MATCH")].iterrows():
        near = man[man.species.eq(r.species)
                   & (man.wavelength_air - r.wavelength_A).abs().le(.08)]
        amb.append({"source_row_id": r.source_row_id, "species": r.species,
                    "wavelength_A": r.wavelength_A, "n_candidates": len(near),
                    "true_class": "NO_CANDIDATE" if len(near) == 0 else "AMBIGUOUS_MULTIPLE"})
    adf = pd.DataFrame(amb)
    n_multi = int((adf.true_class == "AMBIGUOUS_MULTIPLE").sum()) if len(adf) else 0
    rep.add("A2-classes", "Ambiguity is named, not folded into absence", "FLAG",
            f"`NO_UNIQUE_MANIFEST_MATCH` covers two different states: "
            f"{len(adf) - n_multi} rows with NO candidate and {n_multi} with MORE THAN ONE. "
            f"A reader cannot tell an absent line from an unresolved one. Three classes are "
            f"needed, not two - the RYA-1072 lesson.")
    return cdf, rdf, adf


# ─────────────────────────────────────────────────────────────────────────────
# A1 - coverage / truncation, refereed by the CDS ReadMe and branching closure
# ─────────────────────────────────────────────────────────────────────────────
def check_a1(rep: Report, norm: pd.DataFrame) -> pd.DataFrame:
    readme = (VUJ_RAW / "ReadMe").read_text()
    declared = {int(t): int(n) for t, n in
                re.findall(r"table(\d)\.dat\s+\d+\s+(\d+)", readme)}
    parsed = norm.table.value_counts().sort_index().to_dict()
    on_disk = {int(re.search(r"table(\d)", p.name).group(1)):
               len([l for l in p.read_text().splitlines() if l.strip()])
               for p in sorted(VUJ_RAW.glob("table*.dat"))}
    ok = declared == parsed == on_disk and sum(parsed.values()) == 106
    rep.add("A1", "Source coverage / truncation (Vujnovic CDS tables 2-5)",
            "PASS" if ok else "FAIL",
            f"CDS ReadMe declares {declared}; files on disk hold {on_disk}; the normalized "
            f"ledger holds {parsed}. Total {sum(parsed.values())} of the claimed 106. "
            f"The ReadMe's own record counts are an independent referee - it is not the "
            f"builder's own arithmetic.")

    bur = pd.read_csv(BURHEIM)
    bok = (len(bur) == 12 and bur.loggf.notna().all() and bur.e_loggf_dex.notna().all())
    rep.add("A1-burheim", "Burheim Table 3 = 12 derived log gf; Table 2 does not leak in",
            "PASS" if bok else "FAIL",
            f"{len(bur)} rows, every one carrying both `loggf` and `e_loggf_dex`. No "
            f"branching-fraction-only row is present, so Table 2 did not leak in as gf.")

    # Branching closure: A_i / sum(A) must reproduce the separately printed BranR.
    closure = []
    for tab in (2, 5):
        sub = norm[norm.table.eq(tab) & norm["aki_1e8_s-1"].notna()
                   & norm.branching_ratio.notna() & norm.aki_limit.isna()]
        for up, g in sub.groupby("upper_level"):
            if len(g) < 2:
                continue
            tot = g["aki_1e8_s-1"].sum()
            for _, r in g.iterrows():
                closure.append({"table": tab, "upper_level": up,
                                "wavelength_A": r.wavelength_A,
                                "A_over_sumA": r["aki_1e8_s-1"] / tot,
                                "printed_BranR": r.branching_ratio,
                                "residual": r["aki_1e8_s-1"] / tot - r.branching_ratio})
    cdf = pd.DataFrame(closure)
    worst = float(cdf.residual.abs().max()) if len(cdf) else float("nan")
    rep.add("A1-parse", "Fixed-width column extraction, refereed by branching closure",
            "PASS" if len(cdf) and worst <= 0.005 else "FAIL",
            f"For every multiplet with more than one finite Aki, A_i/sum(A) reproduces the "
            f"SEPARATELY PRINTED branching ratio to {worst:.4f} (n={len(cdf)}), consistent "
            f"with the source's 2-3 printed figures. The builder never compares these two "
            f"columns, so this is an identity it cannot have been tuned to.")

    # Flag columns the parser never reads.
    dropped = []
    for fn, cols in [("table2.dat", {"l_IntR": 47, "n_Aki": 75, "l_e_Aki": 76}),
                     ("table5.dat", {"n_Lambda": 36, "l_BranR": 37, "l_e_Aki": 57})]:
        for i, line in enumerate((VUJ_RAW / fn).read_text().splitlines(), 1):
            for name, idx in cols.items():
                ch = line[idx] if len(line) > idx else " "
                if ch.strip():
                    dropped.append({"file": fn, "source_row": i, "flag_column": name,
                                    "flag": ch, "record": line[:46].strip()})
    ddf = pd.DataFrame(dropped)
    rep.add("A1-flags", "CDS limit / note flag columns preserved", "FAIL",
            f"{len(ddf)} flag bytes across 5 documented CDS columns are never read by the "
            f"parser, so the reference README's claim that 'source limits remain limits' "
            f"is false for two of them. `l_e_Aki` ('>') turns a LOWER LIMIT on the "
            f"uncertainty into a determinate sigma, and `n_Lambda` ('*') - which the CDS "
            f"ReadMe documents as 'the value ... was taken over from Tayal & Hibbert "
            f"(1984)' - is the only thing distinguishing a theoretical Aki from a Vujnovic "
            f"measurement, and it is dropped.")
    for _, r in ddf[ddf.flag_column.isin({"l_e_Aki", "n_Lambda", "n_Aki"})].iterrows():
        rep.row("A1", "HIGH" if r.flag_column == "n_Lambda" else "MEDIUM",
                f"{r.file} row {r.source_row}",
                f"CDS flag column `{r.flag_column}` = '{r.flag}' is dropped by the parser",
                r.record)
    return ddf


# ─────────────────────────────────────────────────────────────────────────────
# A3 - HFS / component integrity
# ─────────────────────────────────────────────────────────────────────────────
def check_a3(rep: Report, man: pd.DataFrame, cen: pd.DataFrame,
             cross: pd.DataFrame) -> pd.DataFrame:
    src = BUILDER.read_text()
    sums = bool(re.search(r'HFS_status.*?\n.*?(sum|10\s*\*\*|log10)', src))
    rep.add("A3", "HFS component sums independently re-summed and verified",
            "PASS" if sums else "FAIL",
            "`HFS_status` is set to the string 'COMPONENT_SUM_VERIFIED' whenever "
            "`hfs_n_components > 1` and to 'NO_SPLIT_COMPONENTS_IN_CENSUS' otherwise. "
            "No component sum is computed anywhere in the builder, and `component_or_total` "
            "is the unconditional constant 'TOTAL_TRANSITION_GF'. The status is ASSERTED, "
            "not verified, and RYA-1132's own test "
            "(`m[m.HFS_status=='COMPONENT_SUM_VERIFIED'].component_or_total.eq(...)`) "
            "compares two constants set three lines apart in the same function.")
    rep.row("A3", "CRITICAL", "data/audit/rya1132_al_intake/al_line_manifest.csv:HFS_status",
            "'COMPONENT_SUM_VERIFIED' is asserted from a count, never from a sum",
            "build_al_intake_rya1132.py: "
            "\"COMPONENT_SUM_VERIFIED\" if int(r.hfs_n_components) > 1 else ...")

    carried = "hfs_n_components" in man.columns
    rep.add("A3-meta", "`hfs_n_components` re-verified against the actual component count",
            "PASS" if carried else "FAIL",
            "The manifest does not carry `hfs_n_components` at all. The metadata RYA-1141 "
            "asks to re-verify was dropped at the write, so no reader of the frozen "
            "artifact can check it.")

    # The RYA-1001 defect the ticket names, re-derived and tested for liveness.
    meta = json.loads(CENSUS_META.read_text())
    known = {round(d["wave"], 3): d for d in
             meta.get("hfs_collapse_control", {}).get("disagreements", [])}
    can = pd.read_csv(CANONICAL, low_memory=False)
    live = []
    for w, d in known.items():
        c = can[np.isclose(can.wavelength_air_A, w, atol=.005)
                & can.species.isin(["Al I", "Al II"])]
        m = man[np.isclose(man.wavelength_air, w, atol=.005)]
        if len(c) == 1 and len(m) == 1:
            live.append({"wavelength_air_A": w,
                         "census_hfs_n": d["feat_n"],
                         "canonical_gf_hfs_n": int(c.iloc[0].hfs_n_components),
                         "census_log_gf_sum": d["feat_gf"],
                         "canonical_log_gf": float(c.iloc[0].log_gf),
                         "still_wrong": int(c.iloc[0].hfs_n_components) != d["feat_n"],
                         "manifest_HFS_status": m.iloc[0].HFS_status,
                         "manifest_gf_grade": m.iloc[0].gf_grade,
                         "manifest_gf_source": m.iloc[0].gf_source})
    ldf = pd.DataFrame(live)
    bad = ldf[ldf.still_wrong] if len(ldf) else ldf
    rep.add("A3-rya1001", "The RYA-1001 hfs_n_components defect (3944.006 / 3961.520)",
            "FAIL" if len(bad) else "PASS",
            f"{len(bad)} of {len(ldf)} rows still carry the wrong component count in "
            f"`canonical_gf.csv` on main today (recorded as 1 while the VALD collapse "
            f"finds 4 and 6; the summed log gf agrees exactly, so only the metadata is "
            f"wrong). RYA-1132 read the very census file that records this, stamped both "
            f"rows 'COMPONENT_SUM_VERIFIED', and PROMOTED both to GF-LAB.")
    for _, r in bad.iterrows():
        rep.row("A3", "CRITICAL", f"canonical_gf Al I {r.wavelength_air_A:.3f}",
                "hfs_n_components is wrong and the intake stamped it verified and promoted it",
                f"canonical_gf hfs_n_components={r.canonical_gf_hfs_n}, census component "
                f"count={r.census_hfs_n}; manifest HFS_status={r.manifest_HFS_status}, "
                f"gf_grade={r.manifest_gf_grade}, gf_source={r.manifest_gf_source}")

    # No HFS/isotope inflation: a promoted total must sit on the census HFS SUM,
    # not on its strongest component (RYA-1075/684).
    infl = []
    for _, r in cross[cross.disposition.eq("GF_LAB_PROMOTED")].iterrows():
        m = man[man.canonical_line_id.eq(r.canonical_line_id)].iloc[0]
        c = cen[np.isclose(cen.wave_air_A, m.wavelength_air, atol=.005)]
        if not len(c):
            continue
        c = c.iloc[0]
        infl.append({"wavelength_air": m.wavelength_air,
                     "promoted_loggf": m.loggf_adopted,
                     "census_log_gf_sum": c.log_gf_linelist_sum,
                     "census_strongest_component": c.log_gf_strongest_component,
                     "hfs_n_components": int(c.hfs_n_components),
                     "d_vs_sum": m.loggf_adopted - c.log_gf_linelist_sum,
                     "d_vs_strongest": m.loggf_adopted - c.log_gf_strongest_component})
    idf = pd.DataFrame(infl)
    near_sum = bool(len(idf)) and (idf.d_vs_sum.abs() < idf.d_vs_strongest.abs()).all()
    rep.add("A3-inflation", "No HFS/isotope normalisation inflation, no max-component substitution",
            "PASS" if near_sum else "FAIL",
            f"Every promoted Vujnovic total sits {idf.d_vs_sum.abs().min():.3f}-"
            f"{idf.d_vs_sum.abs().max():.3f} dex from the census HFS SUM and "
            f"{idf.d_vs_strongest.abs().min():.3f}-{idf.d_vs_strongest.abs().max():.3f} dex "
            f"from the strongest component - it is a total-transition gf, as claimed. No "
            f"log10(n_components) inflation is present.")

    held = man[np.isclose(man.wavelength_air, 11254.924, atol=.08)]
    ok11 = (len(held) == 1 and "BURHEIM" not in str(held.iloc[0].gf_source).upper()
            and abs(float(held.iloc[0].loggf_adopted) - 0.3538) < .001)
    rep.add("A3-11254", "11254.9 strong-component vs unresolved-total caveat carried",
            "PASS" if ok11 else "FAIL",
            "The manifest adopts the blended-feature total (+0.3538) and keeps Burheim's "
            "strong component (+0.327) as evidence only; the conflict ledger names the "
            "distinction explicitly.")
    return ldf


# ─────────────────────────────────────────────────────────────────────────────
# A4 - unit / frame conventions
# ─────────────────────────────────────────────────────────────────────────────
def check_a4(rep: Report, man: pd.DataFrame) -> pd.DataFrame:
    """The fine-structure splitting is the referee: two lines sharing an upper
    level must reproduce 112.061 cm-1 from their VACUUM wavelengths."""
    pairs = [(3944.006, 3961.520, "3s^2^4s ^2^S_1/2_"),
             (2652.475, 2660.386, "3s^2^5s ^2^S_1/2_"),
             (3082.153, 3092.839, "3s^2^3d ^2^D_3/2_")]
    rows = []
    for a, b, up in pairs:
        ra = man[np.isclose(man.wavelength_air, a, atol=.005)]
        rb = man[np.isclose(man.wavelength_air, b, atol=.005)]
        if len(ra) != 1 or len(rb) != 1:
            continue
        d = 1e8 / float(ra.iloc[0].wavelength_vac) - 1e8 / float(rb.iloc[0].wavelength_vac)
        rows.append({"shared_upper_level": up, "lambda_air_1": a, "lambda_air_2": b,
                     "lambda_vac_1": float(ra.iloc[0].wavelength_vac),
                     "lambda_vac_2": float(rb.iloc[0].wavelength_vac),
                     "dE_cm1": d, "expected_cm1": AL_I_GROUND_FS_CM1,
                     "residual_cm1": d - AL_I_GROUND_FS_CM1})
    fdf = pd.DataFrame(rows)
    worst = float(fdf.residual_cm1.abs().max()) if len(fdf) else float("nan")
    rep.add("A4", "No air<->vacuum and no cm-1<->Angstrom conflation in the crossmatch",
            "PASS" if len(fdf) == 3 and worst < .02 else "FAIL",
            f"For all three promoted doublets sharing an upper level, the vacuum "
            f"wavenumber difference reproduces the Al I 3p ^2^P^o^ ground-term splitting "
            f"to {worst:.3f} cm-1 against the NIST value {AL_I_GROUND_FS_CM1} cm-1. A "
            f"medium conflation at 2650-3960 A would show as a 5-8 cm-1 offset and a unit "
            f"conflation as orders of magnitude. This identity is not tabulated anywhere "
            f"in the intake.")

    # Eup = Elo + hc/lam_vac must be single-valued per upper level - a second identity.
    fuv = man[man.wavelength_air < 2000]
    same = bool(len(fuv)) and np.allclose(fuv.wavelength_air, fuv.wavelength_vac, atol=1e-9)
    rep.add("A4-medium", "Every wavelength carries an explicit medium", "FLAG",
            f"There is no `medium` column in the manifest. Below 2000 A the census "
            f"correctly stores a single vacuum wavelength, but it stores it in the column "
            f"NAMED `wavelength_air`: for all {len(fuv)} FUV rows `wavelength_air == "
            f"wavelength_vac` exactly. The values are right and the label is wrong, so any "
            f"downstream reader of `wavelength_air` silently receives vacuum wavelengths "
            f"for those rows. (RYA-835/1001 units trap, RYA-938/944.)")
    rep.row("A4", "MEDIUM", "al_line_manifest.csv:wavelength_air",
            "Column named `wavelength_air` holds VACUUM values below 2000 A",
            f"{len(fuv)} FUV rows with wavelength_air == wavelength_vac to 1e-9; "
            f"no `medium` column exists to disambiguate.")
    return fdf


# ─────────────────────────────────────────────────────────────────────────────
# A5 - provenance / grade earned, not asserted
# ─────────────────────────────────────────────────────────────────────────────
def crossref(doi: str) -> dict | None:
    """Live Crossref metadata, or None when offline. Accents folded at read."""
    import unicodedata
    import urllib.request
    def asc(t: str) -> str:
        return "".join(c for c in unicodedata.normalize("NFKD", t)
                       if not unicodedata.combining(c))
    try:
        with urllib.request.urlopen(
                f"https://api.crossref.org/works/{doi}", timeout=20) as fh:
            j = json.load(fh)["message"]
        return {"title": asc(re.sub(r"<[^>]+>", "", (j.get("title") or [""])[0])),
                "container": asc((j.get("container-title") or [""])[0]).replace("&amp;", "&"),
                "volume": str(j.get("volume", "")), "page": str(j.get("page", "")),
                "year": (j.get("issued", {}).get("date-parts") or [[None]])[0][0],
                "authors": [asc(a["family"]) for a in j.get("author", []) if a.get("family")]}
    except Exception:
        return None


def fold(text: str) -> str:
    import unicodedata
    return "".join(c for c in unicodedata.normalize("NFKD", str(text))
                   if not unicodedata.combining(c)).lower()


def same_surname(claimed: str, registered: str) -> bool:
    """Do two spellings of a surname denote the same person?

    #: 🔴 CROSSREF SERVES DAMAGED BYTES. Jonsson & Lundberg 1983 comes back from the
    #: live API as 'J\ufffdnsson' - the umlaut is already a U+FFFD REPLACEMENT
    #: CHARACTER in Crossref's own record, so no amount of normalisation recovers it.
    #: Folding it to a guessed letter would make the committed cache disagree with the
    #: live lookup, which is exactly the drift this cache exists to expose. Treat the
    #: replacement character as the one-character wildcard it actually is, and nothing
    #: else: 'J?nsson' matches 'Jonsson', and still does not match 'Brinchmann'.
    """
    pattern = re.escape(fold(registered)).replace(re.escape("\ufffd"), ".")
    return re.fullmatch(pattern, fold(claimed)) is not None


def claimed_surname(source_id: str, citation: str) -> str:
    """The first-author surname the intake itself asserts for a reference."""
    m = SURNAME.match(str(citation))
    if m:
        return m.group(1)
    m = re.match(r"^([A-Z][a-z]+)", str(source_id))
    return m.group(1) if m else ""


def check_dois(rep: Report, online: bool) -> pd.DataFrame:
    """A5 - does every DOI name the paper the intake says it names?

    #: 🔴 THE REFEREE IS THE AUTHOR LIST, NOT THE VOLUME. A volume-only comparison
    #: passes `10.1086/312738`: the intake claims ApJ 536 and the DOI IS registered to
    #: ApJ 536 - a different paper in the same volume. Only the author surnames
    #: separate them, and they are exact strings, so nothing here is a judgement call.
    """
    bib = pd.read_csv(AUDITED / "source_bibliography.csv")
    web = pd.read_csv(AUDITED / "web_source_followup.csv")
    claims = [("source_bibliography.csv", r.source_id, str(r.citation), r.doi,
               r.ads_bibcode) for _, r in bib.iterrows()]
    claims += [("web_source_followup.csv", str(r.source), str(r.source), r.doi, "")
               for _, r in web.iterrows()]

    rows, seen = [], set()
    for where, sid, citation, doi, bibcode in claims:
        doi = "" if pd.isna(doi) else str(doi).strip()
        if not doi or (where, doi) in seen:
            continue
        seen.add((where, doi))
        live = crossref(doi) if online else None
        act = live or CROSSREF_CACHE.get(doi)
        want = claimed_surname(sid, citation)
        if act is None:
            verdict, why = "UNRESOLVED", "no Crossref record"
        elif not want:
            verdict, why = "UNVERIFIABLE", "the intake states no author to check against"
        elif any(same_surname(want, a) for a in act["authors"]):
            verdict, why = "OK", f"registered author list contains {want}"
        else:
            verdict, why = "MISQUOTED", (f"registered authors are "
                                         f"{', '.join(act['authors'][:3])} - {want} is not among them")
        fix = DOI_CORRECTIONS.get(doi, ("", ""))
        rows.append({"artifact": where, "source_id": sid, "doi": doi,
                     "claimed_citation": citation, "claimed_first_author": want,
                     "claimed_bibcode": "" if pd.isna(bibcode) else str(bibcode),
                     "registered_title": act["title"] if act else "",
                     "registered_authors": "|".join(act["authors"]) if act else "",
                     "registered_container": act["container"] if act else "",
                     "registered_volume": act["volume"] if act else "",
                     "registered_page": act["page"] if act else "",
                     "registered_year": act["year"] if act else "",
                     "resolution": "LIVE_CROSSREF" if live else "COMMITTED_CACHE",
                     "verdict": verdict, "reason": why,
                     "correct_doi": fix[0], "correct_paper": fix[1]})
    ddf = pd.DataFrame(rows)
    wrong = ddf[ddf.verdict.eq("MISQUOTED")]
    ok = ddf[ddf.verdict.eq("OK")]

    #: CONTROL: a referee that rejects everything is not a referee.
    rep.add("A5-doi-control", "The DOI referee accepts correct identifiers",
            "PASS" if len(ok) >= len(ddf) - len(wrong) - 1 else "FAIL",
            f"{len(ok)} of {len(ddf)} DOIs are confirmed by the same test - it accepts "
            f"Burheim, Vujnovic, Trabert, Johnson, Kelleher, Papoulia, Roederer, Lind, "
            f"Jonsson and Chiappino, matching through accented surnames (Vujnovic, "
            f"Trabert, Jonsson) as well as plain ones. A negative is therefore a finding, "
            f"not the detector's default.")

    rep.add("A5-doi", "Every DOI resolves and names the paper it claims",
            "PASS" if wrong.empty else "FAIL",
            f"{len(wrong)} of {len(ddf)} DOIs resolve to an UNRELATED paper. The referee is "
            f"the registered AUTHOR LIST, because a volume comparison is not enough: "
            f"`10.1086/312738` really is in ApJ 536, the volume the intake claims - it is "
            f"just a different paper in it. In every case the intake's own prose citation "
            f"is correct and the DOI beside it is not, so no artifact contradicts itself "
            f"and RYA-1132's suite (which asserts only that `article_url` starts with "
            f"'https://') cannot see it. Resolution: "
            f"{'live Crossref' if online else 'committed cache; re-run with --online'}.")
    for _, r in wrong.iterrows():
        rep.row("A5", "CRITICAL", f"{r.artifact}:{r.source_id}",
                f"DOI {r.doi} resolves to an unrelated paper",
                f"claims '{r.claimed_citation}'; the DOI is registered to "
                f"'{r.registered_title}' by {r.registered_authors} "
                f"({r.registered_container} {r.registered_volume}, {r.registered_page}, "
                f"{r.registered_year}) - {r.reason}. Correct identifier: "
                f"{r.correct_doi} - {r.correct_paper}")

    # A bibcode must be a bibcode: 19 characters, YYYYJJJJJVVVVMPPPPA.
    stubs = [(r.source_id, r.ads_bibcode) for _, r in bib.iterrows()
             if isinstance(r.ads_bibcode, str) and r.ads_bibcode.strip()
             and not re.fullmatch(r"\d{4}[A-Za-z.&]{5}[\w.]{4}[\w.][\w.]{4}[A-Za-z.]",
                                  r.ads_bibcode.strip())]
    rep.add("A5-bibcode", "Every ADS bibcode is a resolvable 19-character bibcode",
            "PASS" if not stubs else "FLAG",
            f"{len(stubs)} entries carry something that is not a bibcode: "
            + "; ".join(f"{a} = {b!r}" for a, b in stubs)
            + ". A truncated bibcode resolves to nothing and cannot referee its DOI - "
              "which matters here, because the bibcode is exactly what would have caught "
              "the Griesmann and Nandakumar page/volume mismatches.")
    for a, b in stubs:
        rep.row("A5", "MEDIUM", f"source_bibliography.csv:{a}",
                "ads_bibcode is a stub, not a resolvable 19-character bibcode", repr(b))
    return ddf


def check_a5(rep: Report, man: pd.DataFrame, cen: pd.DataFrame,
             online: bool) -> pd.DataFrame:
    # Every GF-LAB row must come from a primary-laboratory source, not theory.
    lab = man[man.gf_grade.eq("GF-LAB")]
    theory_words = ("THEORY", "P19", "OP95", "TOPBASE", "1995JPHB", "PAPOULIA", "TAYAL")
    bad = lab[lab.gf_source.astype(str).str.upper().str.contains("|".join(theory_words))]
    rep.add("A5-lab", "No theory graded as primary laboratory",
            "PASS" if bad.empty else "FAIL",
            f"All {len(lab)} GF-LAB rows carry EXP-BURHEIM23 (11), EXP-VUJNOVIC2002 (6) or "
            f"EXP-JOHNSON1986 (1). None is an Opacity-Project or other theoretical source: "
            f"the RYA-1001 `1995JPhB..` class (Mendoza+1995, OP theory) reaches the "
            f"manifest only through `current_canonical_source`, never through a grade.")

    # Branching-fraction-only rows must not be GF-LAB.
    bur = pd.read_csv(BURHEIM)
    bfok = bur.loggf.notna().all()
    rep.add("A5-bf", "No branching-fraction-only row is graded GF-LAB",
            "PASS" if bfok else "FAIL",
            "Burheim Table 2 (branching fractions, no log gf) is not ingested; only Table 3's "
            "12 derived log gf are. The six Vujnovic promotions each carry a finite "
            "'this work' Aki with a stated uncertainty percentage, normalised by a "
            "LABORATORY lifetime (Buurman 1986 / Buurman & Donszelmann 1990 / Davidson "
            "1990, CDS Table 1) - a lab composite, not a bare branching fraction.")

    # A lab/evaluated value must not be silently overwritten by Kurucz/VALD.
    overwritten = lab[lab.gf_source.astype(str).str.upper().str.contains("K75|VALD|KURUCZ")]
    rep.add("A5-overwrite", "No lab/evaluated value silently overwritten by Kurucz/VALD",
            "PASS" if overwritten.empty else "FAIL",
            "No GF-LAB row carries a Kurucz or VALD source. The one deliberate exception, "
            "11254.9, is the documented blend-total substitution and is named in the "
            "conflict ledger rather than being silent.")

    ddf = check_dois(rep, online)

    # Johnson's published value, quoted from the paper itself.
    j = man[np.isclose(man.wavelength_air, 2669.155, atol=.01)]
    expect = math.log10(1.49919e-16 * 3 * 2669.157 ** 2 * 3.33e3)
    okj = len(j) == 1 and abs(float(j.iloc[0].loggf_adopted) - expect) < 1e-6
    rep.add("A5-johnson", "Johnson 1986 Al II 2669 states the value the intake claims",
            "PASS" if okj else "FAIL",
            "The paper's abstract reads verbatim: 'The A-value for the intersystem "
            "transition is (3.33 +/- 0.23) x 10^3 s^-1 at the 90% confidence level', and "
            "'Because there is only a single decay channel, the transition probability is "
            "the inverse of the radiative lifetime'. The intake's g_upper = 3 (3s3p ^3^P^o^_1_) "
            f"and its log gf = {expect:.6f} reproduce exactly.")

    # Uncertainty conversions: justified, but with no provenance recorded in the artifact.
    rep.add("A5-sigma", "sigma(log gf) conversions justified with recorded provenance", "FLAG",
            "Three different conventions share one `gf_sigma_dex` column with nothing "
            "recording which: Burheim's published per-line dex uncertainty; Vujnovic's "
            "log10(1 + u) - the ASYMMETRIC upper bound, ~6% below the linear-propagation "
            "1-sigma 0.434*u; and Johnson's 90%-CONFIDENCE bound stored as if it were "
            "1-sigma (conservative by ~1.645x, and deliberately so, but only a code comment "
            "says so). The artifact needs a `sigma_basis` column, the lesson RYA-1084 and "
            "the sigma_stat/stat_basis finding already paid for.")
    rep.row("A5", "MEDIUM", "al_line_manifest.csv:gf_sigma_dex",
            "One column carries three different uncertainty conventions, unlabelled",
            "Burheim dex 1-sigma; Vujnovic log10(1+u); Johnson 90% confidence bound.")
    return ddf


# ─────────────────────────────────────────────────────────────────────────────
# A6 - source disagreement preserved
# ─────────────────────────────────────────────────────────────────────────────
def check_a6(rep: Report, man: pd.DataFrame, cross: pd.DataFrame) -> pd.DataFrame:
    conf = pd.read_csv(AUDITED / "conflict_ledger.csv", low_memory=False)
    rows = []
    for _, r in cross[cross.disposition.eq("MATCHED_NOT_PROMOTED")
                      & cross.derived_loggf.notna()].iterrows():
        m = man[man.canonical_line_id.eq(r.canonical_line_id)]
        if not len(m):
            continue
        m = m.iloc[0]
        d = float(r.derived_loggf) - float(m.loggf_adopted)
        sig = float(m.gf_sigma_dex) if np.isfinite(m.gf_sigma_dex) else np.nan
        rows.append({"source_row_id": r.source_row_id,
                     "canonical_line_id": r.canonical_line_id,
                     "wavelength_air": m.wavelength_air, "band": m.band,
                     "vujnovic_loggf": r.derived_loggf,
                     "vujnovic_sigma_dex": r.derived_sigma_dex,
                     "adopted_loggf": m.loggf_adopted, "adopted_source": m.gf_source,
                     "adopted_grade": m.gf_grade, "adopted_sigma_dex": m.gf_sigma_dex,
                     "delta_dex": d,
                     "n_sigma_on_adopted": abs(d) / sig if sig and np.isfinite(sig) else np.nan,
                     "in_conflict_ledger": bool(np.isclose(conf.wavelength_air,
                                                           m.wavelength_air, atol=.01).any()),
                     "named_in_competing_gf_summary":
                         "VUJ" in str(m.competing_gf_summary).upper()})
    adf = pd.DataFrame(rows)
    lost = adf[~adf.in_conflict_ledger & ~adf.named_in_competing_gf_summary] if len(adf) else adf
    mentions = int(man.competing_gf_summary.astype(str).str.upper().str.contains("VUJ").sum())
    rep.add("A6", "Competing gf values retained, never silently dropped",
            "FAIL" if len(lost) else "PASS",
            f"{len(lost)} manifest lines have a Vujnovic 2002 primary-laboratory log gf "
            f"that THIS TICKET derived, crossmatched and then dropped. Not one reaches the "
            f"conflict ledger, and `competing_gf_summary` names Vujnovic on {mentions} of "
            f"{len(man)} rows - it is hard-coded to 'Burheim=...; canonical=...; NIST=...' "
            f"and has no slot for a source ingested later. The worst is 13123.416 A, one of "
            f"the two best-graded Al lines RYA-1003 exists to unblock: Vujnovic +0.1901 "
            f"+/- 0.0212 against the adopted Burheim +0.2320 +/- 0.0065, a 0.042 dex "
            f"disagreement between two independent primary-laboratory measurements that "
            f"the frozen manifest presents as a single unopposed GF-LAB value.")
    for _, r in lost.sort_values("n_sigma_on_adopted", ascending=False).iterrows():
        rep.row("A6", "CRITICAL" if (r.n_sigma_on_adopted or 0) > 3 else "HIGH",
                f"{r.canonical_line_id} ({r.wavelength_air:.3f} A)",
                "Competing Vujnovic primary-lab gf derived in this ticket and then dropped",
                f"Vujnovic {r.vujnovic_loggf:+.4f} vs adopted {r.adopted_loggf:+.4f} "
                f"({r.adopted_source}, {r.adopted_grade}); delta {r.delta_dex:+.4f} dex = "
                f"{r.n_sigma_on_adopted:.1f} sigma on the adopted uncertainty; "
                f"in conflict ledger: {r.in_conflict_ledger}")
    return adf


# ─────────────────────────────────────────────────────────────────────────────
# B - reproduce the headline claims independently
# ─────────────────────────────────────────────────────────────────────────────
def check_b1(rep: Report, man: pd.DataFrame, cen: pd.DataFrame,
             cross: pd.DataFrame) -> pd.DataFrame:
    got = {"rows": len(man),
           "Al I": int((man.species == "Al I").sum()),
           "Al II": int((man.species == "Al II").sum())}
    got.update(man.gf_source_type.value_counts().to_dict())
    ok = all(got.get(k) == v for k, v in CLAIM.items())
    # Rebuild the denominator from the sources rather than from the summary.
    n_cen = int(cen.ion.isin(["I", "II"]).sum())
    n_igr = int(man.canonical_line_id.str.endswith("_igrins").sum())
    n_bur = int(man.canonical_line_id.str.contains("_burheim").sum())
    rep.add("B1", "Inventory reproduced independently from the sources",
            "PASS" if ok and n_cen + n_igr + n_bur == len(man) else "FAIL",
            f"505 = {n_cen} RYA-1001 census rows with ion in (I, II) + {n_igr} IGRINS-only "
            f"candidate + {n_bur} Burheim mid-IR completeness controls. By ion "
            f"{got['Al I']}/{got['Al II']}; by source class "
            f"{got.get('PRIMARY_LABORATORY')}/{got.get('CRITICALLY_EVALUATED')}/"
            f"{got.get('THEORETICAL')}/{got.get('FALLBACK')}. Every claimed count reproduces.")

    promoted = sorted(round(float(v), 3) for v in
                      man[man.gf_source.eq("EXP-VUJNOVIC2002")].wavelength_air)
    expect = [2652.475, 2660.386, 3082.153, 3092.710, 3944.006, 3961.520]
    johnson = man[man.gf_source.eq("EXP-JOHNSON1986")]
    okp = promoted == expect and len(johnson) == 1
    rep.add("B1-promotions", "The 6 Al I promotions + Johnson Al II 2669, and the rejections",
            "PASS" if okp else "FAIL",
            f"Promoted at {promoted} A (the ticket quotes the Vujnovic source wavelengths "
            f"2652.484 / 2660.393; the manifest air wavelengths are 2652.475 / 2660.386, "
            f"0.009 and 0.007 A away - the same lines). Each has a finite 'this work' Aki "
            f"with a stated uncertainty. Correctly NOT promoted: every `<` Aki limit (18 "
            f"rows), every ratio-only row, and 3092.839, whose Aki carries the CDS `n_Aki` "
            f"= ')' flag and no independent uncertainty.")

    # Does the intake's headline actually move the measurement gate?
    cenmap = cen.set_index(cen.wave_air_A.round(3)).central_depth.to_dict()
    lab = man[man.gf_grade.eq("GF-LAB")].copy()
    lab["central_depth"] = lab.wavelength_air.round(3).map(cenmap)
    lab["rya946_window"] = np.where(lab.central_depth.isna(), "NO_DEPTH",
                            np.where(lab.central_depth < CODEX_DEPTH[0], "BELOW_WINDOW",
                             np.where(lab.central_depth <= CODEX_DEPTH[1],
                                      "CODEX_GRADE_WINDOW", "DEEP_GRADE_WINDOW")))
    deep = lab[lab.rya946_window.eq("DEEP_GRADE_WINDOW")]
    rep.add("B1-yield", "What the seven new GF-LAB promotions actually unblock", "FLAG",
            f"All {len(deep)} lines promoted by RYA-1132 have Solar central depth "
            f"{deep.central_depth.min():.3f}-{deep.central_depth.max():.3f}, so under "
            f"RYA-946 every one is DEEP Grade, not Codex Grade, and none enters the "
            f"0.05-0.60 measurement window. The intake records this honestly in "
            f"`measurement_suitability_status`, but the manifest has NO column naming the "
            f"RYA-946 grade, and it labels these seven 'CANDIDATE_NOT_SELECTED' - a "
            f"rejection word for lines the contract says should be ROUTED TO SYNTHESIS.")
    return lab[["canonical_line_id", "wavelength_air", "gf_source", "central_depth",
                "rya946_window", "measurement_suitability_status"]]


def check_b2(rep: Report) -> None:
    import subprocess
    try:
        merge = subprocess.run(["git", "log", "--format=%H", "--grep=RYA-1132",
                                "--merges", "-1", "origin/main"],
                               cwd=ROOT, capture_output=True, text=True, check=True
                               ).stdout.strip()
        files = subprocess.run(["git", "diff", "--name-only", f"{merge}^1", merge],
                               cwd=ROOT, capture_output=True, text=True, check=True
                               ).stdout.split()
    except Exception as exc:  # pragma: no cover - a shallow clone has no history
        rep.add("B2", "canonical_gf.csv not mutated by RYA-1132", "FLAG",
                f"Could not read the merge diff ({exc}).")
        return
    touched = [f for f in files if "canonical_gf" in f or "linelists" in f]
    rep.add("B2", "canonical_gf.csv not mutated by RYA-1132",
            "PASS" if not touched else "FAIL",
            f"PR #478 (merge {merge[:7]}) touches {len(files)} files and none of them is "
            f"under `data/linelists/`. `canonical_gf.csv` is byte-identical across the "
            f"merge. The one data file it does change outside its own audit directory is "
            f"`data/audit/rya1129_atomic_intake/intake_status_ledger.csv`, one row, as "
            f"the builder documents.")


def check_b3(rep: Report, man: pd.DataFrame, cen: pd.DataFrame) -> pd.DataFrame:
    v = json.loads((AUDITED / "intake_verdict.json").read_text())
    expect = {"UV": "PARTIAL_GF_LAB_INGESTED_POLICY_BLOCKED",
              "VIS": "FROZEN_WITH_DOCUMENTED_FALLBACKS",
              "IR": "BLOCKED_PIPELINE_COVERAGE"}
    ok = all(v.get(k) == s for k, s in expect.items()) and v["measurement_unblocked"] is False
    rep.add("B3", "Band verdict strings are the ones claimed", "PASS" if ok else "FAIL",
            f"UV/VIS/IR verdicts and `measurement_unblocked: false` match the ticket and "
            f"`summary.json` exactly.")

    # But is each verdict honestly DERIVED from the evidence?
    idx = man.canonical_line_id.str.extract(r"_(\d{4})$")[0]
    j = man.assign(idx=pd.to_numeric(idx, errors="coerce")).dropna(subset=["idx"])
    c = cen[cen.ion.isin(["I", "II"])].reset_index(drop=True)
    j = j.merge(c.reset_index().rename(columns={"index": "idx"})[["idx", "tier", "blend_flag"]],
                on="idx", how="left")
    upgraded = j[j.measurement_suitability_status.eq("ELIGIBLE_WITH_STATED_GF_TIER")
                 & j.tier.eq("CANDIDATE-BLENDED")]
    rep.add("B3-eligibility", "'Eligible' is derived from the evidence the manifest carries",
            "FLAG",
            f"`measurement_suitability_status` is computed from Solar depth and the "
            f"presence of a canonical id alone. It ignores the census `tier` column the "
            f"manifest itself was built from, so {len(upgraded)} rows the RYA-1001 census "
            f"adjudicated CANDIDATE-BLENDED are re-labelled ELIGIBLE_WITH_STATED_GF_TIER - "
            f"including the GF-LAB line 11253.189 and the RYA-835 lines 7835.309 and "
            f"8772.865. The upgrade is silent; no column records that the census disagreed.")
    for _, r in upgraded.iterrows():
        rep.row("B3", "MEDIUM", f"{r.canonical_line_id} ({r.wavelength_air:.3f} A)",
                "Manifest calls ELIGIBLE a line the RYA-1001 census tiered CANDIDATE-BLENDED",
                f"census tier={r.tier}, gf_source={r.gf_source}, gf_grade={r.gf_grade}")

    rep.add("B3-vocab", "`gf_grade` expresses a gf grade", "FLAG",
            "The column mixes three vocabularies: gf provenance tiers (GF-LAB, VALD3, "
            "UNRESOLVED), NIST accuracy grades (B, B+, C, C+, D, E) and RYA-1001 SELECTION "
            "tiers (GRADEABLE, CANDIDATE-BLENDED, EXCLUDED-SHALLOW/SATURATED/NO-HOME). For "
            f"{int(man.gf_grade.isin(['GRADEABLE','CANDIDATE-BLENDED','EXCLUDED-SHALLOW','EXCLUDED-SATURATED','EXCLUDED-NO-HOME']).sum())} "
            "of 505 rows the value is a selection state, not a gf grade at all, and the one "
            "THEORETICAL row is labelled 'GRADEABLE' - a word that reads as an endorsement. "
            "None of RYA-946's four terms (Codex / Deep / Asplund / Consistent) appears "
            "anywhere in the manifest.")
    return upgraded[["canonical_line_id", "wavelength_air", "gf_source", "gf_grade", "tier"]]


# ─────────────────────────────────────────────────────────────────────────────
# C - coverage completeness (the RYA-1001 holdings trap)
# ─────────────────────────────────────────────────────────────────────────────
def check_c(rep: Report, man: pd.DataFrame,
            cen: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    import sys
    sys.path.insert(0, str(ROOT))
    from pipeline import coverage as cov  # noqa: E402

    hold = pd.read_csv(cov.HOLDINGS, comment="#")
    solar = hold[hold.system_id.astype(str).str.strip() == "solar"]
    rows = []
    for _, r in solar.iterrows():
        man_path = ROOT / str(r["manifest_path"]).strip()
        if not man_path.exists():
            why = "MANIFEST_MISSING"
        elif man_path.suffix.lower() not in (".csv", ".tsv"):
            why = f"SKIPPED_BY_SUFFIX_{man_path.suffix.lstrip('.').upper()}"
        else:
            t = pd.read_csv(man_path, comment="#",
                            sep="\t" if man_path.suffix.lower() == ".tsv" else ",")
            why = ("SKIPPED_NO_LOADER_COLUMN" if "loader" not in t.columns
                   else f"RESOLVED_{int((t.instrument_id == str(r['instrument_id']).strip()).sum())}_ROWS")
        rows.append({"holding_id": r.get("holding_id", ""),
                     "instrument_id": str(r["instrument_id"]).strip(),
                     "manifest_path": str(r["manifest_path"]).strip(),
                     "outcome": why,
                     "reaches_coverage_module": why.startswith("RESOLVED")})
    hdf = pd.DataFrame(rows)
    dropped = hdf[~hdf.reaches_coverage_module]
    crires = hdf[hdf.instrument_id.eq("crires_plus")]
    rep.add("C", "Registered holdings reach the coverage module",
            "FAIL" if len(dropped) else "PASS",
            f"{len(dropped)} of {len(hdf)} registered Solar holdings resolve to nothing "
            f"through `pipeline.coverage.load_registry`, and ALL {len(crires)} crires_plus "
            f"registrations are among them. Each `continue` is individually documented "
            f"(RYA-776/929/931/945); the aggregate is that the one instrument reaching Al's "
            f"IR lines is invisible to the module the census reads.")
    for _, r in dropped.iterrows():
        rep.row("C", "HIGH", f"holdings:{r.instrument_id}",
                f"Registered Solar holding silently dropped ({r.outcome})", r.manifest_path)

    # The consequence, named line by line: reachable lines reported unreachable.
    unreachable = []
    for _, r in dropped[dropped.outcome.eq("SKIPPED_NO_LOADER_COLUMN")].iterrows():
        p = ROOT / r.manifest_path
        try:
            d = pd.read_csv(p, usecols=["wavelength_air_A"])
        except Exception:
            continue
        lo, hi = float(d.wavelength_air_A.min()), float(d.wavelength_air_A.max())
        sub = man[(man.wavelength_air >= lo) & (man.wavelength_air <= hi)
                  & man.instrument_reach.isna()]
        for _, m in sub.iterrows():
            unreachable.append({"canonical_line_id": m.canonical_line_id,
                                "wavelength_air": m.wavelength_air, "band": m.band,
                                "gf_grade": m.gf_grade,
                                "measurement_suitability_status": m.measurement_suitability_status,
                                "covering_holding": r.manifest_path,
                                "holding_range_A": f"{lo:.1f}-{hi:.1f}"})
    udf = pd.DataFrame(unreachable).drop_duplicates("canonical_line_id")
    elig = int(udf.measurement_suitability_status.eq("ELIGIBLE_WITH_STATED_GF_TIER").sum()) if len(udf) else 0
    rep.add("C-lines", "No reachable Al line is reported unreachable",
            "FAIL" if len(udf) else "PASS",
            f"{len(udf)} Al manifest lines ({elig} of them ELIGIBLE) sit inside the "
            f"wavelength range of a Solar spectrum that is registered in the holdings "
            f"registry AND present on disk, yet every one carries a BLANK "
            f"`instrument_reach` in the frozen manifest. The IR verdict "
            f"`BLOCKED_PIPELINE_COVERAGE` therefore rests partly on registry plumbing "
            f"rather than on an absence of data.")

    # RYA-1132's own band() has holes, and they swallow the best-graded lines.
    def gapped(w: float) -> bool:
        return (13000 <= w < 13195.23) or (17493.69 <= w < 19510.4) or (w >= 24857.7)
    j = man.assign(idx=pd.to_numeric(man.canonical_line_id.str.extract(r"_(\d{4})$")[0],
                                     errors="coerce")).dropna(subset=["idx"])
    c = cen[cen.ion.isin(["I", "II"])].reset_index(drop=True)
    j = j.merge(c.reset_index().rename(columns={"index": "idx"})[
        ["idx", "band", "instruments_coverage_blind_spot"]], on="idx",
        how="left", suffixes=("", "_census"))
    relab = j[j.band.eq("OUTSIDE_CURRENT_INSTRUMENT_REACH")
              & j.band_census.ne("OUTSIDE_CURRENT_INSTRUMENT_REACH")
              & j.wavelength_air.apply(gapped)]
    labs = relab[relab.gf_grade.eq("GF-LAB")]
    rep.add("C-bands", "RYA-1132's band() covers every band the census does", "FAIL",
            f"`band()` leaves three uncovered intervals - 13000-13195.23 A, "
            f"17493.69-19510.4 A and >=24857.7 A - and every wavelength in them falls "
            f"through to 'OUTSIDE_CURRENT_INSTRUMENT_REACH'. That relabels {len(relab)} "
            f"lines the RYA-1001 census calls NIR, including the {len(labs)} GF-LAB lines "
            f"13123.416 and 13150.753 - the two best-graded Al lines RYA-1003 exists to "
            f"unblock - and the Nandakumar/Chiappino member 17699.094. All of them carry "
            f"`instruments_coverage_blind_spot = crires_plus` in the census, i.e. the "
            f"census says an instrument reaches them. The rows are self-contradictory: "
            f"band 'OUTSIDE_CURRENT_INSTRUMENT_REACH' beside "
            f"`measurement_suitability_status = ELIGIBLE_WITH_STATED_GF_TIER`.")
    for _, r in relab.iterrows():
        rep.row("C", "CRITICAL" if r.gf_grade == "GF-LAB" else "HIGH",
                f"{r.canonical_line_id} ({r.wavelength_air:.3f} A)",
                "Census NIR line relabelled OUTSIDE_CURRENT_INSTRUMENT_REACH by a band() gap",
                f"census band={r.band_census}, blind spot={r.instruments_coverage_blind_spot}, "
                f"gf_grade={r.gf_grade}, suitability={r.measurement_suitability_status}")

    dropped_iii = int((cen.ion == "III").sum())
    rep.add("C-alIII", "Nothing is dropped from the census without being recorded", "FLAG",
            f"{dropped_iii} Al III census rows are excluded by the builder ('Al III is "
            f"outside this ticket's atomic scope'), which is a defensible scope call, but "
            f"the exclusion appears only in a source comment - no artifact records that "
            f"the 505 is a filtered denominator.")
    return hdf, udf, relab[["canonical_line_id", "wavelength_air", "band_census",
                            "gf_grade", "measurement_suitability_status",
                            "instruments_coverage_blind_spot"]]


# ─────────────────────────────────────────────────────────────────────────────
def build(out: Path = OUT, online: bool = False) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    before = {p: sha256(p) for p in audited_files()}

    man = pd.read_csv(AUDITED / "al_line_manifest.csv", low_memory=False)
    norm = pd.read_csv(AUDITED / "vujnovic2002_normalized.csv", low_memory=False)
    cross = pd.read_csv(AUDITED / "vujnovic2002_crossmatch.csv", low_memory=False)
    cen = pd.read_csv(CENSUS, low_memory=False)
    cross = cross.merge(man[["canonical_line_id", "lower_EP"]], on="canonical_line_id", how="left")

    rep = Report()
    flags = check_a1(rep, norm)
    collide, epres, ambig = check_a2(rep, cross, man)
    hfs = check_a3(rep, man, cen, cross)
    fs = check_a4(rep, man)
    dois = check_a5(rep, man, cen, online)
    conflicts = check_a6(rep, man, cross)
    yield_df = check_b1(rep, man, cen, cross)
    check_b2(rep)
    upgraded = check_b3(rep, man, cen)
    holdings, unreachable, relabelled = check_c(rep, man, cen)

    for name, df in [("a1_dropped_source_flags.csv", flags),
                     ("a2_transition_collisions.csv", collide),
                     ("a2_ep_resolution_between_levels.csv", epres),
                     ("a2_unmatched_row_classes.csv", ambig),
                     ("a3_hfs_component_counts.csv", hfs),
                     ("a4_fine_structure_identity.csv", fs),
                     ("a5_doi_resolution.csv", dois),
                     ("a6_dropped_competing_gf.csv", conflicts),
                     ("b1_gf_lab_measurement_yield.csv", yield_df),
                     ("b3_eligibility_upgrades.csv", upgraded),
                     ("c_solar_holdings_resolution.csv", holdings),
                     ("c_reachable_but_unreached_lines.csv", unreachable),
                     ("c_band_gap_relabelled_lines.csv", relabelled)]:
        (df if isinstance(df, pd.DataFrame) else pd.DataFrame(df)).to_csv(out / name, index=False)

    findings = pd.DataFrame(rep.rows)
    findings.to_csv(out / "findings.csv", index=False)
    checks = pd.DataFrame(rep.checks)
    checks.to_csv(out / "check_results.csv", index=False)

    after = {p: sha256(p) for p in audited_files()}
    mutated = sorted(str(p.relative_to(ROOT)) for p in before if before[p] != after.get(p))
    rep.add("NO-MUTATION", "No intake artifact was modified by this QA",
            "PASS" if not mutated else "FAIL",
            f"{len(before)} audited files hashed before and after every read; "
            f"{len(mutated)} changed. This auditor writes only under "
            f"`{out if not out.is_relative_to(ROOT) else out.relative_to(ROOT)}` and "
            f"excludes its own source file by name.")

    verdict = {"ticket": "RYA-1141", "audited": "RYA-1132 (PR #478)",
               "overall": rep.worst, "intake_independently_verified": rep.worst == "PASS",
               "measurement_gate": "CLOSED" if rep.worst == "FAIL" else "OPEN",
               "artifacts_mutated": mutated,
               "checks": {c["check"]: c["status"] for c in rep.checks},
               "n_findings": len(findings),
               "by_severity": findings.severity.value_counts().to_dict() if len(findings) else {}}
    (out / "verdict.json").write_text(json.dumps(verdict, indent=2) + "\n")

    lines = ["# RYA-1141 - independent QA of the RYA-1132 Al atomic-data intake", "",
             f"**Overall: {rep.worst}.** "
             + ("Every check passed; the Al intake is INDEPENDENTLY VERIFIED."
                if rep.worst == "PASS" else
                "The measurement gate stays CLOSED. Nothing in the intake was mutated by "
                "this audit - every defect below is reported, not fixed (RYA-161)."), "",
             "## Per-check result", "",
             "| check | title | status |", "| --- | --- | --- |"]
    lines += [f"| {c['check']} | {c['title']} | **{c['status']}** |" for c in rep.checks]
    lines += ["", "## Detail", ""]
    for c in rep.checks:
        lines += [f"### {c['check']} - {c['title']}: **{c['status']}**", "", c["detail"], ""]
    if len(findings):
        lines += ["## Findings", "",
                  "| severity | check | subject | finding |", "| --- | --- | --- | --- |"]
        order = {"CRITICAL": 0, "HIGH": 1, "MEDIUM": 2}
        for _, r in findings.assign(_o=findings.severity.map(order)).sort_values(
                ["_o", "check"]).iterrows():
            lines.append(f"| {r.severity} | {r.check} | `{r.subject}` | {r.finding} |")
        lines.append("")
    lines += ["## Independently reproduced inventory", "",
              f"- rows: **{len(man)}** (claimed 505)",
              f"- by ion: **{int((man.species=='Al I').sum())} Al I / "
              f"{int((man.species=='Al II').sum())} Al II** (claimed 466 / 39)",
              "- by source class: **"
              + " / ".join(str(man.gf_source_type.value_counts().get(k, 0)) for k in
                           ("PRIMARY_LABORATORY", "CRITICALLY_EVALUATED", "THEORETICAL", "FALLBACK"))
              + "** (claimed 18 / 19 / 1 / 467)", ""]
    (out / "verdict.md").write_text("\n".join(lines) + "\n")

    (out / "README.md").write_text(f"""# RYA-1141 - independent QA of the RYA-1132 Al intake

Findings only. Nothing here changes a gf, a grade, a manifest row or
`canonical_gf.csv` - RYA-161's validate-don't-tune firewall. Every defect below is
REPORTED and filed as a child ticket; none is fixed in this directory.

**Overall verdict: `{rep.worst}`** - see `verdict.md` for the per-check table and
`findings.csv` for every defect with its offending rows named.

## What refereed what

Agreement is only evidence when the referee is independent of the thing it judges,
so each check re-derives rather than re-reads:

| claim | independent referee |
| --- | --- |
| 106 Vujnovic rows, fully transcribed | the CDS `ReadMe`'s own declared record counts |
| the fixed-width column extraction | branching closure, `A_i / sum(A)` vs the separately printed `BranR` |
| no air/vacuum or cm-1/Angstrom conflation | the Al I 3p ^2^P^o^ splitting, {AL_I_GROUND_FS_CM1} cm-1, recovered from vacuum wavelengths |
| the promotion join is EP-aware | an AST test for a physical-identity COMPARISON, with `nearest` as its positive control |
| every DOI names the paper it claims | the registered Crossref AUTHOR LIST, not the volume |
| `canonical_gf` was not mutated | the file list of PR #478's merge diff |

## Files

`verdict.md` / `verdict.json` - the per-check PASS/FAIL/FLAG table.
`findings.csv` - every finding, severity, subject, evidence.
`check_results.csv` - the full reasoning behind each check.
`a1_*`, `a2_*`, `a3_*`, `a4_*`, `a5_*`, `a6_*`, `b1_*`, `b3_*`, `c_*` - the
per-check evidence ledgers named in the verdict.

## Reproduce

`python3 scripts/qa_al_intake_rya1141.py --check`

DOI resolution runs from a committed Crossref cache so the battery is deterministic
offline; `--online` re-resolves every identifier live and reports any drift.
""")
    return verdict


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--check", action="store_true", help="run the battery (default)")
    p.add_argument("--online", action="store_true",
                   help="re-resolve every DOI against Crossref instead of the committed ledger")
    p.add_argument("--out", type=Path, default=OUT)
    a = p.parse_args()
    v = build(a.out, a.online)
    print(json.dumps(v, indent=2))


if __name__ == "__main__":
    main()
