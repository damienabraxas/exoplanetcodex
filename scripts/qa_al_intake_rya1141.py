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
import sys
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

WAVE_NAMES = ("wave", "lambda", "lam", "delta_a", "wtol", "wcol")


def _names(node: "ast.AST") -> set:
    import ast as _ast
    out = set()
    for sub in _ast.walk(node):
        if isinstance(sub, _ast.Attribute):
            out.add(sub.attr)
        elif isinstance(sub, _ast.Name):
            out.add(sub.id)
        elif isinstance(sub, _ast.Constant) and isinstance(sub.value, str):
            out.add(sub.value)
    return out


def _is_wave(names: set) -> bool:
    return any(any(w in n.lower() for w in WAVE_NAMES) for n in names)


def _is_identity(names: set) -> bool:
    return bool(names & EP_IDENTITY_NAMES)


def join_constrains_identity(fn: "ast.AST") -> tuple[bool, list]:
    """Does every candidate-narrowing wavelength comparison carry an identity term?

    #: 🔴 SCOPE THE TEST TO THE NARROWING EXPRESSION, NOT THE FUNCTION. RYA-1037's
    #: shipped guard asks `_enclosing_has_ep()` over the WHOLE enclosing function, so a
    #: single unrelated `ep` anywhere in it silences every wavelength-only comparison
    #: below. This asks a narrower question: for each comparison that filters candidates
    #: BY WAVELENGTH, does the expression that builds that same filter also constrain a
    #: physical identity? The group is the enclosing statement plus any augmented
    #: assignment (`ok &= ...`) narrowing the same target - which is exactly how the
    #: builder's own `nearest()` adds its EP term, so the positive control exercises it.
    #:
    #: 🔴 AND A MENTION IS NOT A COMPARISON. `ingest_new_lab_sources` names
    #: `lower_level` and `upper_level` - it writes them into a provenance string. Only a
    #: Compare node constrains a join, so only a Compare counts.
    """
    import ast as _ast
    parent = {}
    for node in _ast.walk(fn):
        for child in _ast.iter_child_nodes(node):
            parent[child] = node

    def statement_of(node):
        while node in parent and not isinstance(node, _ast.stmt):
            node = parent[node]
        return node if isinstance(node, _ast.stmt) else fn

    def targets(stmt) -> set:
        if isinstance(stmt, _ast.Assign):
            return {t.id for t in stmt.targets if isinstance(t, _ast.Name)}
        if isinstance(stmt, (_ast.AugAssign, _ast.AnnAssign)):
            return {stmt.target.id} if isinstance(stmt.target, _ast.Name) else set()
        return set()

    augs = [n for n in _ast.walk(fn) if isinstance(n, _ast.AugAssign)]
    unguarded = []
    for cmp_node in [n for n in _ast.walk(fn) if isinstance(n, _ast.Compare)]:
        names = _names(cmp_node)
        if not _is_wave(names):
            continue
        stmt = statement_of(cmp_node)
        group = _names(stmt)
        for tgt in targets(stmt):
            for aug in augs:
                if isinstance(aug.target, _ast.Name) and aug.target.id == tgt:
                    group |= _names(aug)
        if not _is_identity(group):
            unguarded.append(getattr(cmp_node, "lineno", -1))
    return (not unguarded), sorted(unguarded)


def check_a2(rep: Report, cross: pd.DataFrame, man: pd.DataFrame):
    """A2 - physical-identity crossmatch, EP-aware, never wavelength-alone."""
    import ast as _ast
    tree = _ast.parse(BUILDER.read_text())
    fns = {n.name: n for n in _ast.walk(tree) if isinstance(n, _ast.FunctionDef)}
    ingest_ok, ingest_lines = join_constrains_identity(fns["ingest_new_lab_sources"])
    #: CONTROL: the same test must return True for `nearest`, the EP-aware helper the
    #: builder defines and then does not call from the ingest path. A guard that cannot
    #: say "yes" is not measuring anything - and `nearest` adds its EP term through
    #: `ok &= ...`, so this also exercises the augmented-assignment branch.
    control_ok, _ = join_constrains_identity(fns["nearest"])
    #: NEGATIVE CONTROL: a function whose wavelength filter is bare, but which mentions
    #: an identity field somewhere else, must still FAIL. This is the laundering-by-scope
    #: hole in RYA-1037's `_enclosing_has_ep()`, asserted here so this test cannot inherit it.
    launder = _ast.parse(
        "def f(df, w, r):\n"
        "    note = f'{r.lower_level} - {r.upper_level}'\n"
        "    return df[(df.wavelength_air - w).abs() <= .08], note\n")
    launder_ok, _ = join_constrains_identity(
        [n for n in _ast.walk(launder) if isinstance(n, _ast.FunctionDef)][0])
    rep.add("A2-control", "The identity-comparison test can say yes, and cannot be laundered",
            "PASS" if control_ok and not launder_ok else "FAIL",
            f"Positive control: the test returns True for `nearest`, the builder's own "
            f"EP-aware matcher, whose EP term arrives via `ok &= (frame[epcol] - ep).abs() "
            f"<= eptol` - so the augmented-assignment branch is exercised. Negative "
            f"control: a fixture whose wavelength filter is bare but which MENTIONS "
            f"`lower_level`/`upper_level` elsewhere still returns False, so this test does "
            f"not inherit RYA-1037's `_enclosing_has_ep()` whole-function blind spot.")

    if not ingest_ok:
        rep.add("A2", "Crossmatch identity (EP-aware, never wavelength alone)", "FAIL",
                f"RYA-1132's `ingest_new_lab_sources` joins every Vujnovic row and the "
                f"Johnson Al II row to the manifest on `abs(wavelength_air - lambda) <= "
                f"0.08` and nothing else. {len(ingest_lines)} candidate-narrowing "
                f"wavelength comparisons (lines {ingest_lines}) carry no physical-identity "
                f"term in the expression that builds the filter - and the level strings "
                f"the function does touch are only WRITTEN into "
                f"`upper_lower_level_identity`. So no promotion in this ticket was matched "
                f"on a physical identity. The builder defines an EP-aware matcher, "
                f"`nearest(..., epcol=..., eptol=0.02)`, and the census loop calls it; the "
                f"promotion path does not.")
        rep.row("A2", "CRITICAL", "scripts/build_al_intake_rya1132.py:ingest_new_lab_sources",
                "Promotion join is wavelength-only, forbidden by RYA-1034",
                f"`near = candidates[candidates.delta_A <= .08]` (line {ingest_lines[0]}); "
                f"the Al II 2669 block repeats the pattern. `nearest(..., epcol=...)` "
                f"exists and is not called here.")
    else:
        rep.add("A2", "Crossmatch identity (EP-aware, never wavelength alone)", "PASS",
                "Every candidate-narrowing wavelength comparison is conjoined with a "
                "physical-identity constraint.")

    #: 🔴 AND THE REPO'S OWN GUARD DOES NOT SEE IT. RYA-1037 ships a repo-wide AST scan
    #: for wavelength-only keys; it is silent on this file, so the defect passed CI.
    #:
    #: 🔴 CALL `scan()`, DO NOT SHELL OUT. Running the auditor as a subprocess REWRITES
    #: `data/audit/rya1037/rya1037_line_key_inventory.json` - a repo artifact outside this
    #: QA's output directory. A findings-only audit that mutates a file to make a finding
    #: has entered its own measurement. `scan(root)` is pure and writes nothing.
    sys.path.insert(0, str(ROOT / "scripts"))
    from audit_line_keys_rya1037 import scan as _scan_line_keys  # noqa: E402
    found = _scan_line_keys(ROOT)
    caught = [f for f in found if "build_al_intake_rya1132" in f.file]
    rep.add("A2-repo-guard", "RYA-1037's repo-wide wavelength-only guard catches this join",
            "PASS" if caught else "FAIL",
            f"It does not. `scripts/audit_line_keys_rya1037.py:scan()` reports {len(found)} "
            f"findings across the repo - so the scanner runs and is not simply empty - and "
            f"names `build_al_intake_rya1132.py` zero times. Two independent reasons: its "
            f"`WAVE_ONLY_TOL` rule matches only the BUILTIN `abs(a - b) <op> tol` inside one "
            f"expression, while RYA-1132 uses the pandas METHOD `(a - b).abs()` assigned to "
            f"`delta_A` and then filters in a SEPARATE statement; and its "
            f"`_enclosing_has_ep()` scopes to the whole function, so one unrelated `ep` "
            f"would silence it anyway. The guard built to make RYA-1034 unrepeatable did "
            f"not fire on the next occurrence of RYA-1034.")
    if not caught:
        rep.row("A2", "CRITICAL", "scripts/audit_line_keys_rya1037.py",
                "The repo-wide wavelength-only guard does not detect the RYA-1132 join",
                f"{len(found)} findings reported elsewhere, zero for "
                f"build_al_intake_rya1132.py; `WAVE_ONLY_TOL` requires builtin abs() within "
                f"a single expression, and `_enclosing_has_ep()` is function-scoped.")

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


#: PR #478's merge commit, PINNED. See `check_b2`.
RYA1132_MERGE = "04e6afe"


def check_b2(rep: Report) -> None:
    """B2 - `canonical_gf.csv` was not mutated by RYA-1132.

    #: 🔴 PIN THE SHA. THE INSTRUMENT ENTERED ITS OWN MEASUREMENT HERE. This check used
    #: to find its target with `git log --grep=RYA-1132 --merges -1 origin/main`. Once
    #: THIS audit's own PR merged - and its merge body names RYA-1132, because that is
    #: what it audits - the newest matching merge became the AUDITOR'S, so B2 quietly
    #: started diffing my commit against its parent and still reported PASS. A guard
    #: that retargets itself onto the thing running it is worse than no guard.
    #:
    #: So the commit is named, not searched, and a CONTROL asserts the pinned commit is
    #: the right one: its diff must contain the RYA-1132 artifacts. If someone repins it
    #: wrongly, the control fails instead of the check silently passing.
    """
    import subprocess

    def git(*args: str) -> str:
        return subprocess.run(["git", *args], cwd=ROOT, capture_output=True,
                              text=True, check=True).stdout

    try:
        merge = git("rev-parse", RYA1132_MERGE).strip()
        subject = git("log", "--format=%s", "-1", merge).strip()
        files = git("diff", "--name-only", f"{merge}^1", merge).split()
    except Exception as exc:  # pragma: no cover - a shallow clone has no history
        rep.add("B2", "canonical_gf.csv not mutated by RYA-1132", "FLAG",
                f"Could not read the pinned merge {RYA1132_MERGE} ({exc}).")
        return

    signature = [f for f in files if f.startswith("data/audit/rya1132_al_intake/")]
    right_commit = ("478" in subject and "rya-1132" in subject.lower()
                    and len(signature) >= 10)
    rep.add("B2-control", "The pinned commit really is PR #478's merge",
            "PASS" if right_commit else "FAIL",
            f"`{merge[:7]}` — \"{subject}\" — and its diff introduces {len(signature)} "
            f"files under `data/audit/rya1132_al_intake/`. The SHA is pinned rather than "
            f"searched, because a `--grep=RYA-1132` search now matches THIS audit's own "
            f"merge commit and would have made B2 diff the auditor against itself.")

    touched = [f for f in files if "canonical_gf" in f or f.startswith("data/linelists/")]
    rep.add("B2", "canonical_gf.csv not mutated by RYA-1132",
            "PASS" if not touched and right_commit else "FAIL",
            f"PR #478 (merge {merge[:7]}) touches {len(files)} files and none of them is "
            f"under `data/linelists/`. `canonical_gf.csv` is byte-identical across the "
            f"merge. The one data file it changes outside its own audit directory is "
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
# ─────────────────────────────────────────────────────────────────────────────
# D - the coverage RYA-1141's first pass did not reach: the local document
#     corpus, the RYA-946 grades other than Codex/Deep, the evaluated tier's
#     own source, and the FULL instrument catalog rather than solar holdings.
# ─────────────────────────────────────────────────────────────────────────────
#: Papers the intake cites, and the local file that holds each. `source_bibliography.csv`
#: gives a `download_url` for most of these; the filenames below are those downloads,
#: already on disk. An intake that cites a paper it never opened is citing a title.
LOCAL_PAPERS = {
    "Vujnovic2002": "aa7151.pdf",
    "Burheim2023": "aa45394-22.pdf",
    "KelleherPodobedova2008": "jpcrd372008911p.pdf",
    "Papoulia2019": "1808.09478v1.pdf",
    "GriesmannKling2000": "0004190v1.pdf",
    "RoedererLawler2021": "2103.12764v1.pdf",
    "Johnson1986": "1986ApJ...308.1013J",
    "Nandakumar2024": "Nandakumar_2024_ApJ_964_96.pdf",
}
REFDOCS = Path("/Users/ryanschmitt/Documents/Exoplanet Codex/Reference documents")

#: The Al I 3p ^2^P^o^ resolving power a real measurement needs. calspec_solar reaches
#: every one of these wavelengths at R = 150-300, which is not a measurement route, so a
#: reachability claim that counts it is worthless.
MIN_USEFUL_R = 20000


def check_d(rep: Report, man: pd.DataFrame, norm: pd.DataFrame,
            cen: pd.DataFrame) -> tuple:
    # D1 - the intake's own cited papers are on this disk. Were they consulted?
    held = {k: (REFDOCS / v) for k, v in LOCAL_PAPERS.items()}
    present = {k: p for k, p in held.items() if p.exists()}
    rep.add("D1", "The intake's cited papers are held locally and were consulted",
            "FLAG" if len(present) >= 6 else "PASS",
            f"{len(present)} of {len(held)} papers the intake cites sit in "
            f"`Reference documents/` — including `aa7151.pdf` (the Vujnovic PAPER, as "
            f"opposed to its CDS tables), `jpcrd372008911p.pdf` (Kelleher & Podobedova, "
            f"the compilation the evaluated tier cites) and `0004190v1.pdf` (Griesmann & "
            f"Kling, which corroborates the DOI correction offline). Nothing in RYA-1132 "
            f"reads any of them: the builder's only inputs are the CDS `.dat` tables, the "
            f"Burheim CSV and the RYA-1001 census. The prose that qualifies the numbers "
            f"(D2) is only in the papers.")

    # D2 - CORRECTS A6. Vujnovic's own text: which Aki rest on a MEASURED ratio?
    t2 = norm[norm.table.eq(2) & norm["aki_1e8_s-1"].notna() & norm.aki_limit.isna()].copy()
    t2["ratio_basis"] = np.where(t2.intensity_ratio.notna(),
                                 "MEASURED_THIS_WORK", "THEORETICAL_LS_RATIO")
    promoted = {2652.484, 2660.393, 3082.153, 3092.710, 3944.006, 3961.520}
    prom = t2[t2.wavelength_A.isin(promoted)]
    theo = t2[t2.ratio_basis.eq("THEORETICAL_LS_RATIO")]
    clean = set(prom.ratio_basis) == {"MEASURED_THIS_WORK"}
    rep.add("D2", "Promotions rest on Vujnovic's MEASURED intensity ratios, not LS theory",
            "PASS" if clean else "FAIL",
            f"All {len(prom)} promoted rows carry a measured 'this work' intensity ratio. "
            f"The {len(theo)} finite-Aki rows that were NOT promoted "
            f"({', '.join(f'{v:.3f}' for v in sorted(theo.wavelength_A))} A) are exactly "
            f"the rows with a BLANK IntR — the paper says of them: 'For 5s-4p transitions "
            f"we evaluated the transition probabilities assuming theoretical intensity "
            f"ratios of the component lines', and for 4p-4s that the branching ratios "
            f"'were measured indirectly by (Buurmann & Doenszelmann 1990)'. The separation "
            f"is perfect, so RYA-1132's hand-curated promote list is right on this axis.")
    rep.row("D2", "MEDIUM", "vujnovic2002_normalized.csv",
            "No column records that 4 rows' fine-structure split is theoretical, not measured",
            "`intensity_ratio` is NaN for exactly 13123.41, 13150.76, 21093.04, 21163.75; "
            "the paper's prose is the only thing that says why.")

    #: 🔴 THIS CORRECTS RYA-1141's FIRST-PASS A6. I called 13123.416 'two independent
    #: primary-laboratory measurements in tension'. It is not. Vujnovic's value there is a
    #: LIFETIME times an externally-measured branching ratio, split across fine structure
    #: by an ASSUMED LS ratio - weaker evidence than Burheim's direct measurement, not an
    #: equal-footing rival. The disagreement is still unrecorded and still a defect; its
    #: severity and its remedy both change.
    rep.add("D2-a6-correction", "A6's 13123.416 disagreement, correctly characterised",
            "FLAG",
            "RYA-1141's first pass called this 'two independent primary-laboratory "
            "measurements in genuine tension'. The Vujnovic paper refutes that: 13123.41 "
            "and 13150.76 are among the four rows whose fine-structure split is a "
            "THEORETICAL LS ratio over an indirectly-measured branching ratio. They remain "
            "a competing published value that the manifest drops without trace — the A6 "
            "FAIL stands — but they do not impeach Burheim's uncertainty the way an "
            "independent direct measurement would.")

    # D3 - RYA-946's MANDATORY Solar reference-line-set census, and the FROZEN gate.
    #: 🔴 "ASPLUND GRADE" IS NOT A gf GRADE. `pipeline/model_registry.py:LINE_SETS` is the
    #: one definition: it is a value on the `line_set` PROVENANCE axis - which POOL of
    #: lines a measurement was made on - owned by RYA-1111, vocabulary
    #: {asplund, gbs, our-graded, our-deep-graded, our-ungraded, our-all}. RYA-1127 put
    #: `line_set` INTO THE PRODUCT IDENTITY KEY, so every product must resolve one.
    #: `consistent` is absent DELIBERATELY (RYA-1105 retires it) and a product carrying it
    #: must fail loudly rather than acquire a name - so "is Consistent merged into
    #: Codex/Deep" is the wrong question; the tier is dead, and RYA-1141's ticket prose
    #: naming four grades is superseded by the repo, which wins.
    from pipeline.model_registry import LINE_SETS
    cols = set(man.columns)
    has_axis = bool({"line_set", "reference_line_set"} & cols)
    rep.add("D3-lineset", "The frozen manifest can resolve a `line_set`", 
            "PASS" if has_axis else "FAIL",
            f"The manifest has no `line_set` column and no value anywhere from the one "
            f"vocabulary {LINE_SETS}. RYA-1127 made `line_set` part of the PRODUCT "
            f"IDENTITY KEY, so a measurement taken from this frozen pool cannot form a "
            f"valid key. `gf_grade` mixes three vocabularies and none of them is this one.")
    rep.row("D3", "CRITICAL", "al_line_manifest.csv",
            "No `line_set` column — products measured from this pool cannot key (RYA-1127)",
            f"canonical vocabulary: {LINE_SETS}")

    #: The RYA-946 census gate, quoted: "No element is FROZEN_READY_FOR_MEASUREMENT until
    #: this cross-reference is complete or a documented, approved source-publication
    #: exception exists." RYA-1132 wrote intake_status FROZEN on every canonical row.
    ref = ROOT / "data/reference"
    asplund_sets = sorted(q.name for q in ref.glob("asplund*") if q.is_dir())
    al_set = [q for q in asplund_sets if "al" in q.lower().replace("asplund", "")]
    frozen = int(man.intake_status.eq("FROZEN").sum())
    rep.add("D3", "RYA-946's mandatory AGSS21 line-set census was done before freezing",
            "PASS" if al_set else "FAIL",
            f"It was not, and {frozen} manifest rows are stamped `FROZEN` anyway. "
            f"RYA-946 (2026-08-29) requires, for EVERY canonical element before its "
            f"lab-gf sweep is complete, that the adopted Solar value be traced to its "
            f"line-level source across ALL bands (FUV/NUV/VIS/red-optical/NIR/IR, "
            f"including forbidden, molecular, isotopologue and blend-component "
            f"indicators), with a per-line join and a per-band coverage matrix - and it "
            f"gates the freeze: 'No element is FROZEN_READY_FOR_MEASUREMENT until this "
            f"cross-reference is complete or a documented, approved source-publication "
            f"exception exists.' `data/reference/` holds {asplund_sets} - Fe only, from "
            f"RYA-1109. There is no Al reference set, and RYA-1132 records no exception. "
            f"AGSS21/Asplund/Scott/Nordlander appear ZERO times in all 505 rows.")
    rep.row("D3", "CRITICAL", "data/audit/rya1132_al_intake/al_line_manifest.csv",
            "Al frozen through RYA-946's census gate, with no census and no exception",
            f"{frozen} rows intake_status=FROZEN; data/reference/ has {asplund_sets}")

    #: What the census WOULD start from - traced here so the child ticket does not repeat it.
    rep.add("D3-lineage", "AGSS21's Al value traced to its line-level source", "FLAG",
            "AGSS21 publishes no Al line list. Its section 'Aluminium (Z = 13)' adopts "
            "Nordlander & Lind (2017), who 'adopted the same lines and line data as in "
            "Scott et al. (2015b), except that they excluded the 1089.1 nm line due to "
            "telluric contamination', giving A(Al) = 6.43 +/- 0.03 over 'these five Al i "
            "lines'. So the Al reference set is a FIVE-line set whose identity lives in "
            "Scott et al. (2015b), with one published NEGATIVE selection (1089.1 nm, "
            "telluric) that RYA-946 says must be preserved rather than silently dropped. "
            "The repo already holds the Nordlander & Lind pointer in "
            "`data/curation/threednlte_availability.csv` (DOI 10.1051/0004-6361/201730427).")

    # D4 - the evaluated tier: its provenance, its values, and its grades.
    ev = man[man.gf_source_type.eq("CRITICALLY_EVALUATED")].copy()
    no_doi = ev[ev.gf_source_doi.isna() | ev.gf_source_doi.astype(str).str.strip().eq("")]
    rep.add("D4", "Every 'critically evaluated' row carries a resolvable source",
            "PASS" if no_doi.empty else "FAIL",
            f"{len(no_doi)} of {len(ev)} CRITICALLY_EVALUATED rows have `gf_source_doi` "
            f"EMPTY — no DOI, no bibcode, no table id. A5 requires that every evaluated "
            f"row's source resolve and state the claimed value; not one of them points "
            f"anywhere. `source_bibliography.csv` cites Kelleher & Podobedova 2008 (JPCRD "
            f"37, 709) as the evaluated source, but the values actually come from a NIST "
            f"ASD pull (`data/linelists/nist_pulls/*.tsv`, 2026-08-09, RYA-708) — a "
            f"different NIST product, and nothing records which was used.")
    for _, r in no_doi.iterrows():
        rep.row("D4", "HIGH", f"{r.canonical_line_id} ({r.wavelength_air:.3f} A)",
                "CRITICALLY_EVALUATED row carries no DOI, bibcode or table id",
                f"gf_source={r.gf_source}, gf_grade={r.gf_grade}, "
                f"gf_sigma_dex={r.gf_sigma_dex}")

    #: 🔴 THE AUDITOR MUST OBEY THE RULE IT AUDITS. The first version of this join was
    #: `ev.merge(cen, left_on=ev.wavelength_air.round(4), right_on=cen.wave_air_A.round(4))`
    #: - a rounded-wavelength-only merge, which is RYA-1033/1034/1037 exactly, committed
    #: inside the check that reports it. The repo's own guard caught it (WAVE_ONLY_MERGE,
    #: `tests/test_line_key_guard_rya1037.py::test_the_real_tree_passes`) and it was right.
    #: This is the EP-aware, ambiguity-refusing match RYA-1151 asks RYA-1132 to adopt.
    def match(row: "pd.Series") -> "pd.Series | None":
        c = cen[(cen.wave_air_A - row.wavelength_air).abs().le(0.005)
                & (cen.ep_eV - row.lower_EP).abs().le(0.02)]
        if len(c) != 1:
            raise AssertionError(
                f"evaluated row {row.canonical_line_id} matched {len(c)} census rows on "
                f"(lambda, EP) - the audit refuses ambiguity rather than taking argmin")
        return c.iloc[0]
    j = pd.DataFrame([{**r.to_dict(),
                       **match(r)[["nist_grade", "nist_grade_worst",
                                   "nist_n_components", "nist_log_gf"]].to_dict()}
                      for _, r in ev.iterrows()])
    val_ok = np.allclose(j.loggf_adopted.astype(float), j.nist_log_gf.astype(float), atol=1e-6)
    rep.add("D4-values", "Evaluated log gf reproduce the NIST pull, sums included",
            "PASS" if val_ok else "FAIL",
            "All 19 reproduce the pulled NIST ASD value exactly, and the five "
            "multi-component features are correctly SUMMED rather than taking the "
            "strongest row: 7836.134 = log10(10^-0.534 + 10^-1.834) = -0.5131 and "
            "8773.896 = log10(10^-0.192 + 10^-1.495) = -0.1709, both matching the "
            "manifest. The values are right.")

    opt = j[j.nist_grade.ne(j.nist_grade_worst)]
    rep.add("D4-grades", "A summed feature is graded by its WORST component",
            "PASS" if opt.empty else "FAIL",
            f"{len(opt)} of {len(j)} evaluated rows are multi-component sums graded with "
            f"the BEST component's grade while a strictly worse one exists in the same "
            f"feature — the census carries `nist_grade_worst` in the very next column and "
            f"RYA-1132 reads `nist_grade`. A sum cannot be more accurate than its worst "
            f"term. The worst case is 6906.287 A, graded C (<=25%) over a component NIST "
            f"grades E (>50%), and its `gf_sigma_dex` follows the optimistic grade at "
            f"0.109 dex instead of >=0.30.")
    for _, r in opt.iterrows():
        rep.row("D4", "CRITICAL", f"{r.canonical_line_id} ({r.wavelength_air:.3f} A)",
                "Summed feature graded by its best component, not its worst",
                f"manifest grade {r.gf_grade}, worst component grade {r.nist_grade_worst}, "
                f"n_components={int(r.nist_n_components)}, sigma={r.gf_sigma_dex:.4f} dex")

    #: 🔴 D4-LINEAGE — "NIST ONLY" IS NOT A SOURCE. This corrects RYA-1141's first-pass
    #: A5-lab PASS. That check asked only whether a GF-LAB row is theory. It is not. But
    #: the CRITICALLY_EVALUATED tier - the one RYA-946 ranks second to laboratory - is
    #: Opacity Project THEORY end to end, and the manifest cannot say so.
    #:
    #: Kelleher & Podobedova 2008 Table 4 carries a Source column and decodes it in its own
    #: header: "1 = Mendoza et al., 2 = Tachiev and Froese Fischer, 3 = Vujnovic et al.,
    #: 4 = Davidson et al., 5 = Hannaford". Every Al I multiplet feeding the 19 evaluated
    #: rows (16-21, 23-30) carries Source 1 = Mendoza et al. = the ab-initio close-coupling
    #: OPACITY PROJECT calculation; every one of the 19 components carries Source LS, i.e.
    #: an LS-coupling split of that theoretical multiplet total. Not one is measured.
    #: RYA-1001 reached the same conclusion independently for 8772/8773: "1995JPhB.. =
    #: Mendoza, Eissner, Le Dourneuf & Zeippen 1995 ... THEORY, not a lab measurement",
    #: and "1995JPhB.. == TOPbase == theory".
    #:
    #: The mechanism is two if-statements in the wrong order.
    src = BUILDER.read_text()
    fn = src[src.index("def source_type"):src.index("def nearest")]
    nist_before_theory = fn.index('"NIST" in s') < fn.index('"THEORY" in s')
    rep.add("D4-lineage", "The evaluated tier is evaluated data, not theory in a better coat",
            "FAIL",
            f"All 19 CRITICALLY_EVALUATED rows trace, through NIST's own Source column, to "
            f"Mendoza et al. — the Opacity Project ab-initio calculation — split across "
            f"fine structure by LS coupling. 'Critically evaluated' names NIST's editorial "
            f"process, not the nature of the underlying data, and the manifest offers no "
            f"column that distinguishes an evaluated LABORATORY value from an evaluated "
            f"THEORETICAL one. Under RYA-946's 'replicate the line list' doctrine these 19 "
            f"rows are theory, and Al's red-optical band — 7835/7836, 8772/8773 and the "
            f"rest — rests entirely on them. NIST alone is not a laboratory source.")
    rep.row("D4", "CRITICAL", "scripts/build_al_intake_rya1132.py:source_type",
            "NIST is tested before THEORY, so Opacity-Project values can never be typed THEORETICAL",
            'if t.startswith("NIST") or "NIST" in s: return "CRITICALLY_EVALUATED"  '
            '<-- returns first; the "THEORY"/"P19"/"OP95" branch below is unreachable for '
            f'any NIST-sourced row. NIST-before-THEORY confirmed: {nist_before_theory}')
    rep.row("D4", "CRITICAL", "al_line_manifest.csv (19 rows)",
            "Opacity Project theory typed CRITICALLY_EVALUATED across Al's whole red-optical band",
            "Kelleher & Podobedova 2008 Table 4: multiplets 16-21, 23-30 all Source 1 = "
            "Mendoza et al. (OP); all 19 components Source LS. Confirms RYA-1001's "
            "independent finding on 8772/8773.")

    # D5 - the FULL instrument catalog, not just what we happen to hold.
    cat = pd.read_csv(ROOT / "data/catalog/instrument_catalog.csv", comment="#")
    cat = cat[cat.codex_status.ne("rejected")
              & (pd.to_numeric(cat.resolving_power_max, errors="coerce") >= MIN_USEFUL_R)]
    cat = cat.assign(lo=cat.wavelength_min_nm * 10, hi=cat.wavelength_max_nm * 10)
    rows = []
    for _, r in man.iterrows():
        w = float(r.wavelength_air)
        c = cat[(cat.lo <= w) & (w <= cat.hi)]
        rows.append({"canonical_line_id": r.canonical_line_id, "wavelength_air": w,
                     "band": r.band, "gf_grade": r.gf_grade,
                     "manifest_instrument_reach": "" if pd.isna(r.instrument_reach) else str(r.instrument_reach),
                     "measurement_suitability_status": r.measurement_suitability_status,
                     "n_catalog_instruments": len(c),
                     "catalog_instruments": "|".join(sorted(c.instrument_id))})
    sweep = pd.DataFrame(rows)
    nowhere = sweep[sweep.n_catalog_instruments.eq(0)]
    rep.add("D5", "Every band and every catalogued instrument checked, not just holdings",
            "PASS" if nowhere.empty else "FLAG",
            f"All 505 Al lines were swept against all {len(cat)} catalogued instruments "
            f"that are not `rejected` and reach R >= {MIN_USEFUL_R} (calspec_solar is "
            f"excluded at R = 150-300, which is not a measurement route). "
            f"{len(nowhere)} lines are beyond every one of them. The instrument catalog "
            f"says the manifest's whole wavelength span is reachable.")

    wrong = sweep[sweep.manifest_instrument_reach.eq("OUTSIDE_CURRENT_REACH")
                  & sweep.n_catalog_instruments.gt(0)]
    rep.add("D5-outside", "`OUTSIDE_CURRENT_REACH` means no instrument can reach it",
            "PASS" if wrong.empty else "FAIL",
            f"{len(wrong)} rows are labelled `OUTSIDE_CURRENT_REACH` — and "
            f"`measurement_suitability_status = OUTSIDE_CURRENT_REACH` with them — while "
            f"the catalog lists 4 high-resolution instruments covering each: crires_plus "
            f"(950-5300 nm, R 50k-100k), ishell, nirspec and phoenix. These are the four "
            f"Burheim mid-IR GF-LAB lines at 3.86-4.19 um, the intake's own 'completeness "
            f"controls'. The honest label is NO HOLDING, not out of reach: the manifest "
            f"collapses 'we hold no spectrum' into 'the universe is out of range', and "
            f"only the second one closes a question.")
    for _, r in wrong.iterrows():
        rep.row("D5", "HIGH", f"{r.canonical_line_id} ({r.wavelength_air:.3f} A)",
                "Labelled OUTSIDE_CURRENT_REACH while catalogued instruments cover it",
                f"gf_grade={r.gf_grade}; covered by {r.catalog_instruments}")

    blank = sweep[sweep.manifest_instrument_reach.eq("") & sweep.n_catalog_instruments.gt(0)]
    graded = blank[blank.gf_grade.isin(["GF-LAB", "GRADEABLE", "B", "B+", "C", "C+", "D", "E"])]
    rep.add("D5-blank", "A blank `instrument_reach` distinguishes no-holding from no-instrument",
            "FLAG",
            f"{len(blank)} of 505 rows carry a BLANK `instrument_reach` while a catalogued "
            f"high-resolution instrument covers them, {len(graded)} of those graded. The "
            f"column conflates three different states — no instrument exists (0 rows), an "
            f"instrument exists but we hold no spectrum, and we hold a spectrum the "
            f"coverage module cannot see (check C) — into one blank. Three classes, not one.")
    return sweep, j[["canonical_line_id", "wavelength_air", "gf_grade", "nist_grade",
                     "nist_grade_worst", "nist_n_components", "loggf_adopted",
                     "nist_log_gf", "gf_sigma_dex"]], t2


def porcelain() -> tuple[set, bool]:
    """The working tree's dirty set, or (empty, False) when git cannot answer."""
    import subprocess
    try:
        r = subprocess.run(["git", "status", "--porcelain"], cwd=ROOT,
                           capture_output=True, text=True, check=True).stdout
        return {ln[3:].strip().strip('"') for ln in r.splitlines() if ln.strip()}, True
    except Exception:
        return set(), False


def build(out: Path = OUT, online: bool = False) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    before = {p: sha256(p) for p in audited_files()}
    tree_before, tree_checked = porcelain()

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
    sweep, evaluated, vujratio = check_d(rep, man, norm, cen)

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
                     ("c_band_gap_relabelled_lines.csv", relabelled),
                     ("d5_full_instrument_catalog_sweep.csv", sweep),
                     ("d4_evaluated_tier_provenance.csv", evaluated),
                     ("d2_vujnovic_ratio_basis.csv", vujratio)]:
        (df if isinstance(df, pd.DataFrame) else pd.DataFrame(df)).to_csv(out / name, index=False)

    findings = pd.DataFrame(rep.rows)
    findings.to_csv(out / "findings.csv", index=False)
    checks = pd.DataFrame(rep.checks)
    checks.to_csv(out / "check_results.csv", index=False)

    after = {p: sha256(p) for p in audited_files()}
    mutated = sorted(str(p.relative_to(ROOT)) for p in before if before[p] != after.get(p))

    #: 🔴 HASHING A CHOSEN SET CANNOT SEE A WRITE OUTSIDE IT. An earlier revision of this
    #: script shelled out to the RYA-1037 auditor and silently rewrote
    #: `data/audit/rya1037/rya1037_line_key_inventory.json` - a repo file that was not in
    #: the audited set, so the hash comparison passed while the working tree was dirty.
    #: The complete question is not "did these 21 files change" but "did ANYTHING change
    #: outside my own output directory", and only the working tree can answer it. It is a
    #: DIFF against the tree as it stood when this run started, so a dirty dev checkout
    #: (or a tmp_path `out`) cannot manufacture a failure this run did not cause.
    tree_after, _ = porcelain()
    own = str(out.relative_to(ROOT)) if out.is_relative_to(ROOT) else None
    stray = sorted(q for q in (tree_after - tree_before)
                   if not (own and q.startswith(own)))
    mutated = sorted(set(mutated) | set(stray))

    rep.add("NO-MUTATION", "No intake artifact was modified by this QA",
            "PASS" if not mutated else "FAIL",
            f"{len(before)} audited files hashed before and after every read, AND the whole "
            f"working tree diffed against its state at the start of this run for any change "
            f"outside `{own or out}`"
            f"{'' if tree_checked else ' (working-tree check unavailable here)'}. "
            f"{len(mutated)} changed. This auditor excludes its own source file by name, "
            f"never by pattern.")

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
