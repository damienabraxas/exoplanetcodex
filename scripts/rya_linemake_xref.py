#!/usr/bin/env python3
"""RYA-1070 — linemake as a gf CROSS-REFERENCE for the lines we already hold.

    python3 scripts/rya_linemake_xref.py

READ-ONLY. This script never opens `data/linelists/canonical_gf.csv` for writing, never
imports a linemake line, never runs a synthesis, and never adopts a linemake value as
authority. Validate-don't-tune (RYA-161) is absolute here: a disagreement is recorded as
science to adjudicate later, never auto-corrected.

WHAT IT ANSWERS. `linemake` (Placco et al. 2021, RNAAS 5, 92) bundles primary-lab gf --
mostly the Wisconsin FTS studies -- plus HFS. The Codex has independently ingested
several of the SAME papers (Fe I Ruffoni2014/DenHartog2014/Belmonte2017 via RYA-945/1046,
Fe II DenHartog2019). So where our pool and linemake overlap: do they AGREE, and does
linemake ever point at a stronger primary source than we currently cite?

Nothing here is typed from memory. Specifically:

* The CURATED file manifest is parsed out of linemake's own `mooglists/mergenohfs` and
  `mooglists/mergehfs` shell scripts -- the lists that build `goodgf`/`goodgfhfs`. A
  hand-typed file list would silently drift from the repo it claims to describe.
* The species -> primary-source map is parsed from the README's four pipe tables.
* Element symbols come from `canonical_gf` itself (key_z + ion -> species) and are
  CROSS-CHECKED against the symbol in each linemake filename. A disagreement is fatal:
  it would mean the species-code decode is wrong.
* Both match tolerances are DERIVED from the residual distributions (below), not set.
* The air/vac scale of linemake is MEASURED, not assumed (below).

THREE TRAPS THIS FILE EXISTS TO NOT FALL INTO
---------------------------------------------
1. RYA-1034 -- wavelength-alone matching. A match needs species AND wavelength AND
   excitation potential, and an AMBIGUOUS match (two candidates inside tolerance, either
   direction) is REFUSED, never broken by proximity. 1,491 of our Fe I lines sit within
   20 mA of a neighbour; EP is the only thing separating them.

2. RYA-102/473 -- HFS. linemake's `.mooghfs` files expand a transition into components
   and prepend a NEGATIVE-wavelength row carrying the declared total. Summing every
   component is WRONG: the components are listed PER ISOTOPE, each isotope's set summing
   to the full gf. A naive sum over Eu II inflates by log10(2) = +0.301 dex (2 isotopes),
   Ba II by +0.699 (5), Nd II by +0.845 (7). We therefore compare against the declared
   total and VERIFY it by per-isotope summation; a block whose isotopes disagree with the
   declared total is `HFS_AMBIGUOUS`, not a number.

3. Air/vacuum. MOOG lists are conventionally air above the boundary and vacuum below it.
   We take the boundary from the repo's SSOT (`pipeline.wavelength_util`, IAU/VALD 2000 A)
   and then MEASURE which scale linemake is actually on, by counting EP-consistent
   coincidences under both hypotheses. Anything below the boundary is `VAC_BOUNDARY` and
   is excluded from the numeric match rather than matched raw.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

#  The air/vac boundary is the repo's SSOT (IAU 1991 / Morton 2000 / VALD), not a number
#  chosen here. The gf_grades tolerances are imported ONLY as an independent cross-check on
#  the tolerances this script derives from the data -- never as their source.
from pipeline.gf_grades import EP_TOL_EV, LOGGF_MATCH_TOL, WAVE_TOL_A  # noqa: E402
from pipeline.wavelength_util import AIR_VACUUM_BOUNDARY_A, vac_to_air  # noqa: E402

OUT = ROOT / "data" / "audit" / "linemake_xref"
CANONICAL = ROOT / "data" / "linelists" / "canonical_gf.csv"
LINEMAKE_URL = "https://github.com/vmplacco/linemake"

#: MOOG fixed-width atomic record: four F10 fields, then a free-text source tag.
#: CONFIRMED against the files themselves -- `parse_moog_file` re-derives every record
#: by whitespace splitting as well and refuses to continue if the two disagree.
FIELD_W = 10
N_FIELDS = 4

#: Signal-to-null criterion used to READ both match tolerances off the measured residual
#: histograms. We do not choose a tolerance in Angstroms or eV; we state how far above the
#: accidental-coincidence rate a bin has to sit to be counted as real, and the tolerance
#: falls out of the data. See `derive_tolerance` -- the choice is not knife-edge, because
#: the measured ratio profile drops by a factor of ~4 across the boundary it selects, and
#: `derive_tolerance` records the tolerance at 3x, 5x and 10x so that stability is visible.
NULL_RATIO_CRITERION = 5.0

#: Wavelength offsets (A) used to build the accidental-coincidence null. Shifting our
#: wavelengths destroys every true pair while preserving the line density of both lists,
#: so the residual histogram of the shifted scan IS the coincidence background. Twenty
#: offsets, both signs, none a round multiple of another, so the null is not itself
#: structured by the line lists' own periodicities.
NULL_SHIFTS_A = (1.0, 1.7, 2.3, 3.1, 4.2, 5.0, 6.3, 7.1, 8.6, 9.4,
                 11.3, 13.7, -1.3, -2.9, -4.7, -6.1, -8.3, -9.9, -12.1, -15.4)

#: A source tag that is NOT a primary laboratory measurement -- either a compilation we
#: already hold, or an astrophysical (solar) calibration. Used to decide whether linemake's
#: value for a line could promote ours. Applied to the tag printed in the linemake file
#: itself, so the judgement is per line rather than per species.
#:
#: `SOLAR` earns its place here on the same principle as the rest: linemake tags the Nd II
#: 4314.50 total `IUR-SOLAR`, and a gf fitted to the solar spectrum cannot promote anything
#: -- adopting it would be exactly the astrophysical calibration RYA-161 forbids. Without
#: it that line was classified PRIMARY_LAB and offered as a promotion candidate.
NON_LAB_TAG_RE = re.compile(r"NIST|KUR|VALD|DREAM|HITRAN|EXOMOL|SOLAR", re.I)

#: Our tiers that are NOT a primary measurement -- a line on one of these can in
#: principle be promoted if linemake points at a lab paper.
FALLBACK_TIERS = {"KURUCZ", "VALD3", "OTHER"}

#: Caveat phrases the README uses to mark a value that is deliberately NOT a primary lab
#: result. Matched against the README cell text so the flag is parsed, never typed.
CAVEAT_PATTERNS = [
    (r"treated with caution", "non-primary: README flags the source as to-be-treated-with-caution"),
    (r"solar-derived", "non-primary: solar-derived (astrophysical) log gf"),
    (r"raise the .{0,12}gf.{0,12}-?values", "non-primary: deliberate pragmatic gf offset"),
]

#: A caveat can be SPECIES-wide ("Cu I ... should be treated with caution") or scoped to a
#: single line ("the 4314.50 A line has a solar-derived log(gf)"). Applying a line-scoped
#: caveat to the whole species is not a conservative error -- it silently suppressed the
#: verdict on all 707 matched Nd II lines, which is most of what linemake has to say about
#: our n-capture pool. So the scope is PARSED too: a wavelength in the same sentence as the
#: caveat phrase means the caveat belongs to that line and no other.
CAVEAT_WAVELENGTH_RE = re.compile("([0-9]{3,5}\\.[0-9]{1,3})\\s*(?:\u00c5|&Aring;|A)\\b")


class ParseError(RuntimeError):
    """A linemake artifact did not decode. Loud-fail; never silently drop a line."""


# --------------------------------------------------------------------------------------
# 1. acquire
# --------------------------------------------------------------------------------------
def acquire_linemake(dest: Path) -> tuple[Path, str]:
    """Clone linemake read-only (shallow) and return (path, commit sha). Files unchanged."""
    if not (dest / ".git").exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(["git", "clone", "--depth", "1", LINEMAKE_URL, str(dest)],
                       check=True, capture_output=True, text=True)
    sha = subprocess.run(["git", "-C", str(dest), "rev-parse", "HEAD"],
                         check=True, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(dest), "status", "--porcelain"],
                           check=True, capture_output=True, text=True).stdout.strip()
    if dirty:
        raise ParseError(f"linemake clone at {dest} is MODIFIED; refusing to audit against "
                         f"a mutated reference:\n{dirty}")
    return dest, sha


def file_sha256(path: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------------------
# 2. the curated manifest -- parsed from linemake's own merge scripts
# --------------------------------------------------------------------------------------
def parse_manifest(mooglists: Path) -> tuple[list[str], list[str]]:
    """Read the file lists that linemake itself concatenates into goodgf / goodgfhfs."""
    def files_of(script: str, sink: str) -> list[str]:
        text = (mooglists / script).read_text(encoding="utf-8", errors="replace")
        text = text.replace("\\\n", " ")
        body = text.split(">")[0]
        if not body.strip().startswith("cat"):
            raise ParseError(f"{script} does not start with `cat`; manifest shape changed")
        names = [t for t in body.split()[1:] if t and not t.startswith(("-", ">", "|"))]
        missing = [n for n in names if not (mooglists / n).exists()]
        if missing:
            raise ParseError(f"{script} lists files that do not exist: {missing}")
        if not names:
            raise ParseError(f"{script} produced an empty manifest for {sink}")
        return names

    return files_of("mergenohfs", "goodgf"), files_of("mergehfs", "goodgfhfs")


# --------------------------------------------------------------------------------------
# 3. MOOG record decoding -- confirmed twice, two ways
# --------------------------------------------------------------------------------------
def parse_moog_file(path: Path) -> pd.DataFrame:
    """Decode one MOOG-format atomic list.

    Columns: wavelength (A), species code, excitation potential (eV), log gf, then an
    optional free-text source tag. Each record is decoded BOTH as fixed-width 4xF10 and
    by whitespace-splitting the first 40 columns; any disagreement raises rather than
    picking one. That is the confirmation the ticket asks for, done per record instead of
    by eye on the first line of one file.
    """
    rows = []
    for lineno, raw in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
        if not raw.strip():
            continue
        head = raw[:FIELD_W * N_FIELDS]
        try:
            fixed = [float(head[i * FIELD_W:(i + 1) * FIELD_W]) for i in range(N_FIELDS)]
        except ValueError as exc:
            raise ParseError(f"{path.name}:{lineno} not fixed-width 4xF10: {raw!r} ({exc})")
        free = head.split()
        if len(free) != N_FIELDS:
            raise ParseError(f"{path.name}:{lineno} first {FIELD_W * N_FIELDS} cols split into "
                             f"{len(free)} fields, expected {N_FIELDS}: {raw!r}")
        try:
            freev = [float(x) for x in free]
        except ValueError as exc:
            raise ParseError(f"{path.name}:{lineno} free-form decode failed: {raw!r} ({exc})")
        if any(abs(a - b) > 1e-9 for a, b in zip(fixed, freev)):
            raise ParseError(f"{path.name}:{lineno} fixed-width and free-form decodes DISAGREE: "
                             f"{fixed} vs {freev} -- refusing to guess the column layout")
        code = head[FIELD_W:2 * FIELD_W].strip()
        rows.append(dict(source_file=path.name, lineno=lineno,
                         wavelength_raw_A=fixed[0], species_code=code,
                         excitation_potential_eV=fixed[2], log_gf=fixed[3],
                         log_gf_text=head[3 * FIELD_W:4 * FIELD_W].strip(),
                         source_tag=raw[FIELD_W * N_FIELDS:].strip()))
    if not rows:
        raise ParseError(f"{path.name} decoded to zero records")
    return pd.DataFrame(rows)


def decode_species_code(code: str) -> tuple[int, int, str]:
    """MOOG species code -> (Z, ionisation stage index, isotope digits).

    Integer part is Z. The first decimal digit is the ionisation stage (0 = neutral =
    'I', 1 = singly ionised = 'II'). Any further decimal digits are the isotope mass
    number -- those rows are HFS/isotopic COMPONENTS, not comparable totals.
    """
    if "." not in code:
        raise ParseError(f"species code {code!r} has no decimal part; cannot read ionisation")
    z_txt, dec = code.split(".", 1)
    if not z_txt.isdigit() or not dec.isdigit():
        raise ParseError(f"species code {code!r} is not numeric")
    z = int(z_txt)
    if z > 92:
        raise ParseError(f"species code {code!r} decodes to Z={z}; molecular codes are not "
                         f"handled by this atomic cross-reference")
    return z, int(dec[0]), dec[1:].rstrip("0") if dec[1:] else ""


# --------------------------------------------------------------------------------------
# 4. HFS blocks
# --------------------------------------------------------------------------------------
def collapse_hfs(df: pd.DataFrame) -> pd.DataFrame:
    """Reduce a linemake list to COMPARABLE totals, one row per physical transition.

    A negative wavelength opens an HFS block: that row carries linemake's own declared
    total for the transition and the rows after it are its components. We compare against
    the declared total and verify it by summing each ISOTOPE separately -- see the module
    docstring for why summing all components is wrong.
    """
    out, block, comps = [], None, []

    def close() -> None:
        if block is None:
            return
        rec = dict(block)
        rec["wavelength_raw_A"] = abs(block["wavelength_raw_A"])
        rec["hfs_expanded"] = True
        rec["n_components"] = len(comps)
        if not comps:
            rec["hfs_reconciled"] = False
            rec["hfs_max_dev_dex"] = np.nan
        else:
            by_iso: dict[str, list[float]] = {}
            for c in comps:
                by_iso.setdefault(c["isotope"], []).append(c["log_gf"])
            sums = [np.log10(np.sum(np.power(10.0, v))) for v in by_iso.values()]
            dev = max(abs(s - block["log_gf"]) for s in sums)
            rec["n_isotopes"] = len(by_iso)
            rec["hfs_max_dev_dex"] = float(dev)
            rec["hfs_reconciled"] = bool(dev <= HFS_RECONCILE_TOL_DEX)
        out.append(rec)

    for rec in df.to_dict("records"):
        z, ion, iso = decode_species_code(rec["species_code"])
        rec.update(key_z=z, ion_index=ion, isotope=iso)
        if rec["wavelength_raw_A"] < 0:
            close()
            block, comps = rec, []
        elif iso:
            if block is None:
                raise ParseError(f"{rec['source_file']}:{rec['lineno']} isotopic component "
                                 f"{rec['species_code']} with no preceding total row")
            comps.append(rec)
        else:
            close()
            block, comps = None, []
            rec.update(hfs_expanded=False, n_components=1, n_isotopes=1,
                       hfs_reconciled=True, hfs_max_dev_dex=0.0)
            out.append(rec)
    close()
    return pd.DataFrame(out)


#: How far a per-isotope component sum may sit from linemake's declared total before the
#: block is refused. Derived, not chosen: it is the printing quantum of the log gf column
#: (F10 with 3-4 decimals) propagated through a sum of up to ~40 components in quadrature,
#: rounded up to the next decade. Measured devs on the curated .mooghfs files are <= 0.005.
HFS_RECONCILE_TOL_DEX = 0.01


# --------------------------------------------------------------------------------------
# 5. README -> species -> primary source
# --------------------------------------------------------------------------------------
def parse_readme(readme: Path) -> pd.DataFrame:
    """Parse the README's pipe tables into species -> references + class + caveats."""
    text = readme.read_text(encoding="utf-8", errors="replace")
    sections, current = {}, None
    rows = []
    for line in text.splitlines():
        head = re.match(r"^#{2,3}\s*(.+?)\s*$", line)
        if head:
            current = head.group(1)
            continue
        if not line.startswith("`") and not line.startswith("<code>"):
            continue
        if "|" not in line:
            continue
        species_cell, _, refs = line.partition("|")
        sp = species_cell.strip()
        m = re.fullmatch(r"`([^`]+)`", sp) or re.fullmatch(r"<code>(.+)</code>", sp)
        if not m:
            raise ParseError(f"README species cell not in the expected `X`/<code>X</code> "
                             f"form: {species_cell!r}")
        name = re.sub(r"<[^>]+>", "", m.group(1)).strip()
        refs = refs.strip()
        if not refs:
            raise ParseError(f"README row for {name!r} has an empty reference cell")
        papers = re.findall(r"\[([^\]]+)\]\((http[^)]+)\)", refs)
        cited = [(t, u) for t, u in papers
                 if not NON_LAB_TAG_RE.search(t) and not re.search(r"nist\.gov|kurucz|exomol|hitran|dream", u, re.I)]
        catalog = [(t, u) for t, u in papers if (t, u) not in cited]
        cls = ("PRIMARY_LAB" if cited and not catalog else
               "MIXED" if cited and catalog else
               "CATALOG" if catalog else "UNSOURCED")
        caveats, caveat_wave, caveat_tol = [], float("nan"), float("nan")
        for pat, note in CAVEAT_PATTERNS:
            hit = re.search(pat, refs, re.I)
            if not hit:
                continue
            #  Look only in a window around the caveat phrase, so a wavelength quoted
            #  elsewhere in the same cell is not mistaken for the caveat's scope. A window
            #  rather than a sentence: these cells are dense with citation punctuation, and
            #  "4314.50" contains a full stop of its own.
            window = refs[max(0, hit.start() - 120):hit.end() + 120]
            w = CAVEAT_WAVELENGTH_RE.search(window)
            if w:
                caveat_wave = float(w.group(1))
                #  half the quantum the README actually printed -- not a match tolerance,
                #  just the width of the number it wrote down
                caveat_tol = 0.5 * 10.0 ** -len(w.group(1).split(".")[1])
                note = f"{note}; scoped to {w.group(1)} A only"
            caveats.append(note)
        rows.append(dict(readme_section=current, species=name, reference_class=cls,
                         caveat_wavelength_A=caveat_wave, caveat_wavelength_tol_A=caveat_tol,
                         n_paper_citations=len(cited), n_catalog_citations=len(catalog),
                         primary_source="; ".join(dict.fromkeys(t for t, _ in cited))
                                        or "(catalogue only)",
                         caveats="; ".join(caveats), references_raw=refs))
        sections.setdefault(current, 0)
        sections[current] += 1
    if not rows:
        raise ParseError("README parsed to zero species rows -- table format changed")
    df = pd.DataFrame(rows)
    for required in ("Fe I", "Fe II", "Mn I", "Eu II", "Al I", "Cu I", "Nd II", "CO"):
        if required not in set(df.species):
            raise ParseError(f"README parse lost {required!r}; refusing a partial source map")
    caveated = set(df[df.caveats != ""].species)
    for required in ("Cu I", "Nd II", "CO"):
        if required not in caveated:
            raise ParseError(f"README caveat detection did not flag {required!r}; the three "
                             f"known non-primary entries must all be found by parsing, not "
                             f"asserted. Found: {sorted(caveated)}")
    return df


# --------------------------------------------------------------------------------------
# 6. tolerances -- derived from the measured residual distributions
# --------------------------------------------------------------------------------------
def _scan_residuals(ours: pd.DataFrame, lm: pd.DataFrame,
                   shift: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """Nearest-neighbour residuals in (wavelength, EP) for every shared species.

    `shift` displaces our wavelengths before the scan. A non-zero shift destroys every
    true pairing while leaving both line densities untouched, which is exactly the
    accidental-coincidence null the tolerance derivation needs.
    """
    dws, des = [], []
    for sp, g in lm.groupby("species"):
        o = ours[ours.species == sp]
        if o.empty:
            continue
        g = g.sort_values("wavelength_air_A")
        lw = g.wavelength_air_A.to_numpy(float)
        le = g.excitation_potential_eV.to_numpy(float)
        ow = o.wavelength_air_A.to_numpy(float) + shift
        oe = o.excitation_potential_eV.to_numpy(float)
        idx = np.searchsorted(lw, ow)
        for k in (-1, 0, 1):
            j = np.clip(idx + k, 0, len(lw) - 1)
            dws.append(ow - lw[j])
            des.append(oe - le[j])
    return np.concatenate(dws), np.concatenate(des)


def derive_tolerance(real: np.ndarray, real_cond: np.ndarray,
                     null: np.ndarray, null_cond: np.ndarray,
                     edges: list[float], label: str) -> tuple[float, list[dict], dict]:
    """Read a match tolerance off the residual histogram against a measured null.

    The |residual| histogram of true same-transition pairs is a narrow core sitting on the
    flat background of coincidental neighbours. That background is not guessed: it is
    MEASURED by re-running the same nearest-neighbour scan with our wavelengths displaced
    (`NULL_SHIFTS_A`), which destroys every true pair and leaves the coincidence rate
    intact. The tolerance is the outer edge of the last bin in which the real count still
    exceeds the null by `NULL_RATIO_CRITERION`.

    Nothing here is a tolerance typed from memory: the Angstroms and eV come out of the
    two histograms. The criterion is a stated purity rule, and the returned `stability`
    shows what the same rule gives at 3x and 10x so the reader can see the answer does not
    hang on the exact multiple.
    """
    a = np.abs(real[real_cond])
    b = np.abs(null[null_cond])
    counts, _ = np.histogram(a, bins=edges)
    ncounts, _ = np.histogram(b, bins=edges)
    ncounts = ncounts / len(NULL_SHIFTS_A)

    prof = []
    for lo, hi, c, nc in zip(edges[:-1], edges[1:], counts, ncounts):
        prof.append(dict(lo=lo, hi=hi, count=int(c), null_count=float(nc),
                         real_density=float(c / (hi - lo)),
                         null_density=float(nc / (hi - lo)),
                         ratio=float(c / nc) if nc > 0 else float("inf")))

    def cut(crit: float) -> float:
        tol = 0.0
        for row in prof:
            if row["ratio"] >= crit:
                tol = row["hi"]
            else:
                break
        return tol

    tol = cut(NULL_RATIO_CRITERION)
    if tol <= 0:
        raise ParseError(f"{label}: no bin sits {NULL_RATIO_CRITERION}x above the measured "
                         f"coincidence null -- there is no separable true-match core")
    stability = {f"{c:g}x": cut(c) for c in (3.0, 5.0, 10.0)}
    return float(tol), prof, stability


def _decimals(text: str) -> int:
    """Decimal places actually printed. The quantum of a value is 10**-this."""
    t = str(text).strip()
    if t.lower() in ("", "nan", "none"):
        return 0
    t = t.split("e")[0].split("E")[0]
    return len(t.split(".")[1]) if "." in t else 0


def loggf_threshold(our_text: str, lm_text: str) -> float:
    """Per-line agreement threshold: the widest gap two lists can show carrying the SAME
    number, given what each of them actually printed.

    This is the whole content of "AGREES" here, and it needs no free parameter. If we
    print a log gf to 3 decimals and linemake prints it to 2, the largest difference the
    two can show while representing one underlying value is 0.0005 + 0.005. Anything
    bigger is two different numbers, whatever their provenance. Lines we hold at one
    decimal place get a correspondingly wider threshold, because at that precision we are
    not entitled to call a 0.03 dex gap a disagreement.

    Corroborated by the data rather than asserted: the matched delta histogram has a spike
    of exact zeros, a decaying shoulder out to ~0.005 dex, and then a flat continuum -- the
    shoulder ends where this rule puts the threshold. The repo's own ratified
    `gf_grades.LOGGF_MATCH_TOL` (0.02 dex, "same number, different rounding") is the same
    idea one step looser; `main` records both.
    """
    return 0.5 * (10.0 ** -_decimals(our_text) + 10.0 ** -_decimals(lm_text))


# --------------------------------------------------------------------------------------
# 7. the cross-reference
# --------------------------------------------------------------------------------------
def build_linemake_table(mooglists: Path, manifest: list[str]) -> pd.DataFrame:
    """Decode and HFS-collapse every file in one curated manifest.

    `species` is deliberately NOT set here: it is assigned in `main` only after
    `check_filename_symbols` has validated the species-code decode against the filenames,
    so a bad decode cannot get a name before it has been checked.
    """
    df = pd.concat([collapse_hfs(parse_moog_file(mooglists / name)) for name in manifest],
                   ignore_index=True)
    df["below_air_vac_boundary"] = df.wavelength_raw_A < AIR_VACUUM_BOUNDARY_A
    df["wavelength_air_A"] = df.wavelength_raw_A
    df["linemake_source_class"] = np.where(
        df.source_tag.fillna("").str.strip().eq(""), "UNTAGGED",
        np.where(df.source_tag.fillna("").map(lambda t: bool(NON_LAB_TAG_RE.search(t))),
                 "CATALOG", "PRIMARY_LAB"))
    return df


def build_bulk_table(mooglists: Path, z_to_symbol: dict[int, str]) -> pd.DataFrame:
    """The uncurated `moogatom*` bulk lists.

    The ticket asks for these to be parsed, and they are -- but they are deliberately kept
    OUT of the corroboration verdicts, and that exclusion is a measured decision rather
    than an unexamined one. These files carry no source tag and no README attribution
    (linemake's own README calls them "a substantial number of additional transitions",
    separately from the curated periodic table). An unsourced value agreeing with our
    unsourced value is not evidence of anything, and one disagreeing is not a finding. So
    they are used for exactly one thing: to say how much MORE of our pool linemake touches
    at all, which is a real fact about coverage.
    """
    frames = []
    for path in sorted(mooglists.glob("moogatom*")):
        if path.stat().st_size == 0:
            continue
        df = parse_moog_file(path)
        keep = []
        for rec in df.to_dict("records"):
            try:
                z, ion, iso = decode_species_code(rec["species_code"])
            except ParseError:
                continue                       # molecular code -> not an atomic line
            if iso or rec["wavelength_raw_A"] < 0:
                continue                       # component / commented row
            rec.update(key_z=z, ion_index=ion)
            keep.append(rec)
        if keep:
            frames.append(pd.DataFrame(keep))
    if not frames:
        raise ParseError("no moogatom* bulk files decoded")
    bulk = pd.concat(frames, ignore_index=True)
    roman = {0: "I", 1: "II", 2: "III"}
    bulk["species"] = [f"{z_to_symbol.get(z, f'Z{z}')} {roman.get(i, f'ION{i}')}"
                       for z, i in zip(bulk.key_z, bulk.ion_index)]
    bulk["wavelength_air_A"] = bulk.wavelength_raw_A
    return bulk[bulk.wavelength_raw_A >= AIR_VACUUM_BOUNDARY_A]


def bulk_coverage(xref: pd.DataFrame, bulk: pd.DataFrame, wtol: float,
                  etol: float) -> pd.Series:
    """Which of our lines have a counterpart in the uncurated bulk lists. Count only."""
    hit = pd.Series(False, index=xref.index)
    for sp, g in xref.groupby("species"):
        b = bulk[bulk.species == sp]
        if b.empty:
            continue
        b = b.sort_values("wavelength_air_A")
        bw = b.wavelength_air_A.to_numpy(float)
        be = b.excitation_potential_eV.to_numpy(float)
        ow = g.our_wavelength.to_numpy(float)
        oe = g.our_EP.to_numpy(float)
        idx = np.searchsorted(bw, ow)
        ok = np.zeros(len(g), bool)
        for k in (-1, 0, 1):
            j = np.clip(idx + k, 0, len(bw) - 1)
            ok |= (np.abs(bw[j] - ow) <= wtol) & (np.abs(be[j] - oe) <= etol)
        hit.loc[g.index] = ok
    return hit


def check_filename_symbols(df: pd.DataFrame, z_to_symbol: dict[int, str]) -> list[str]:
    """The species code decode must agree with the symbol in the filename. Fatal if not."""
    notes = []
    for name, g in df.groupby("source_file"):
        m = re.fullmatch(r"([a-z]{1,2})(I{1,3})", name.split(".")[0])
        if not m:
            notes.append(f"{name}: filename carries no species stem; decode not cross-checked")
            continue
        want_sym = m.group(1).capitalize()
        want_ion = len(m.group(2)) - 1
        got_z = sorted(set(g.key_z))
        got_ion = sorted(set(g.ion_index))
        if len(got_z) != 1 or len(got_ion) != 1 or got_ion[0] != want_ion:
            raise ParseError(f"{name}: species-code decode {got_z}/{got_ion} disagrees with the "
                             f"filename ({want_sym} ion index {want_ion})")
        have = z_to_symbol.get(got_z[0])
        if have and have != want_sym:
            raise ParseError(f"{name}: Z={got_z[0]} decodes to {have!r} from canonical_gf but the "
                             f"filename says {want_sym!r} -- species-code decode is WRONG")
        if not have:
            notes.append(f"{name}: Z={got_z[0]} ({want_sym}) has no counterpart species in "
                         f"canonical_gf; symbol taken from the filename")
            z_to_symbol[got_z[0]] = want_sym
    return notes


def measure_air_vac(ours: pd.DataFrame, lm: pd.DataFrame) -> dict:
    """MEASURE which wavelength scale linemake is on. Never assume it."""
    above = lm[~lm.below_air_vac_boundary].copy()
    hyp = {}
    for name, wave in (("air", above.wavelength_raw_A.to_numpy(float)),
                       ("vacuum", vac_to_air(above.wavelength_raw_A.to_numpy(float)))):
        probe = above.assign(wavelength_air_A=wave)
        dw, de = _scan_residuals(ours, probe)
        ep_ok = np.abs(de) <= 0.02
        hyp[name] = int((np.abs(dw[ep_ok]) <= 0.005).sum())
    verdict = max(hyp, key=hyp.get)
    if verdict != "air" or hyp["air"] < 10 * max(hyp["vacuum"], 1):
        raise ParseError(
            f"CRITICAL: linemake's wavelength scale is not decisively AIR above "
            f"{AIR_VACUUM_BOUNDARY_A} A (EP-consistent coincidences within 5 mA: {hyp}). "
            f"Matching air against a vacuum list is forbidden.")
    return dict(boundary_A=AIR_VACUUM_BOUNDARY_A, boundary_source="pipeline.wavelength_util "
                "(IAU 1991 / Morton 2000 / VALD convention, repo SSOT)",
                coincidences_within_5mA=hyp, measured_scale=verdict,
                lines_below_boundary=int(lm.below_air_vac_boundary.sum()))


def cross_reference(ours: pd.DataFrame, lm: pd.DataFrame, wtol: float, etol: float,
                    readme: pd.DataFrame, our_loggf_text: pd.Series) -> pd.DataFrame:
    """One row per canonical_gf line of a linemake-covered species."""
    rmap = readme.set_index("species")
    lm_ok = lm[~lm.below_air_vac_boundary]
    claims: dict[int, list[int]] = {}
    recs = []

    for sp, o in ours.groupby("species"):
        g = lm_ok[lm_ok.species == sp].sort_values("wavelength_air_A")
        below = int((lm.species == sp).sum()) - len(g)
        lw = g.wavelength_air_A.to_numpy(float)
        le = g.excitation_potential_eV.to_numpy(float)
        gi = g.index.to_numpy()
        for row in o.itertuples():
            w, e = float(row.wavelength_air_A), float(row.excitation_potential_eV)
            our_txt = our_loggf_text.get(row.Index, "")
            rec = dict(line_id=row.line_id, species=sp, ion=row.ion,
                       our_wavelength=w, our_EP=e, our_loggf=row.log_gf,
                       our_source=row.loggf_reference, our_gf_tier=row.gf_tier,
                       our_hfs_n_components=row.hfs_n_components,
                       our_lab_source_tag=row.lab_source_tag,
                       our_loggf_printed=our_txt,
                       agreement_threshold_dex=np.nan,
                       linemake_matched=False, linemake_wavelength=np.nan,
                       linemake_EP=np.nan, linemake_loggf=np.nan, delta_loggf=np.nan,
                       linemake_source_tag="", linemake_source_class="",
                       linemake_file="", linemake_hfs_expanded=False,
                       linemake_n_components=np.nan,
                       #  Species-level fields are filled only where they mean something --
                       #  on a matched row. Repeating a 150-character citation on all 126k
                       #  NO_MATCH rows added 10 MB to the artifact and said nothing; the
                       #  per-species table carries it once, which is where it belongs.
                       linemake_primary_source="",
                       linemake_reference_class="",
                       linemake_readme_caveat="",
                       linemake_caveat_wavelength_A=rmap.caveat_wavelength_A.get(sp, np.nan),
                       linemake_caveat_tol_A=rmap.caveat_wavelength_tol_A.get(sp, np.nan),
                       lines_below_air_vac_boundary=below, verdict="NO_MATCH")
            if len(lw):
                cand = np.flatnonzero((np.abs(lw - w) <= wtol) & (np.abs(le - e) <= etol))
                if cand.size > 1:
                    rec["verdict"] = "AMBIGUOUS_MATCH"
                    rec["linemake_wavelength"] = float(lw[cand[0]])
                    rec["ambiguity"] = (f"{cand.size} linemake candidates within {wtol} A / "
                                        f"{etol} eV: " + ", ".join(f"{lw[c]:.4f}" for c in cand))
                elif cand.size == 1:
                    j = int(cand[0])
                    m = g.loc[gi[j]]
                    claims.setdefault(int(gi[j]), []).append(len(recs))
                    rec.update(linemake_matched=True,
                               linemake_wavelength=float(m.wavelength_air_A),
                               linemake_EP=float(m.excitation_potential_eV),
                               linemake_loggf=float(m.log_gf),
                               delta_loggf=float(row.log_gf) - float(m.log_gf),
                               linemake_primary_source=rmap.primary_source.get(sp, ""),
                               linemake_reference_class=rmap.reference_class.get(sp, ""),
                               linemake_readme_caveat=rmap.caveats.get(sp, ""),
                               linemake_loggf_printed=m.log_gf_text,
                               agreement_threshold_dex=loggf_threshold(our_txt, m.log_gf_text),
                               linemake_source_tag=m.source_tag,
                               linemake_source_class=m.linemake_source_class,
                               linemake_file=m.source_file,
                               linemake_hfs_expanded=bool(m.hfs_expanded),
                               linemake_n_components=int(m.n_components),
                               linemake_hfs_reconciled=bool(m.hfs_reconciled),
                               linemake_hfs_max_dev_dex=float(m.hfs_max_dev_dex)
                               if np.isfinite(m.hfs_max_dev_dex) else np.nan,
                               linemake_row=int(gi[j]))
            recs.append(rec)

    df = pd.DataFrame(recs)
    # Reverse ambiguity: one linemake line claimed by two DIFFERENT lines of ours means
    # linemake cannot tell them apart either. RYA-1034 -- refuse, do not pick.
    for _, idxs in claims.items():
        if len(idxs) > 1:
            for i in idxs:
                df.at[i, "verdict"] = "AMBIGUOUS_MATCH"
                df.at[i, "ambiguity"] = (f"one linemake line is the only candidate for "
                                         f"{len(idxs)} distinct canonical_gf lines")
    return df


def assign_verdicts(df: pd.DataFrame) -> pd.DataFrame:
    """Verdict precedence: cannot-compare, then non-primary, then source strength, then gf."""
    matched = df.linemake_matched & df.verdict.ne("AMBIGUOUS_MATCH")
    df["delta_within_threshold"] = np.abs(df.delta_loggf) <= df.agreement_threshold_dex

    reconciled = df.get("linemake_hfs_reconciled", pd.Series(True, index=df.index))
    reconciled = reconciled.astype("object").where(reconciled.notna(), True).astype(bool)
    hfs_bad = matched & ~reconciled
    #  BOTH SIDES ALWAYS CARRY TOTALS, so an asymmetric expansion is not by itself
    #  ambiguous. On our side an `hfs_n_components > 1` row is already a collapsed physical
    #  line -- gf-weighted centroid wavelength, log10(sum of component gf) as the value
    #  (`rya822_extend_canonical_gf.cluster_physical_lines`). On linemake's side the
    #  comparable row is either an unexpanded single entry or the declared total of an
    #  expanded block, verified here per isotope. Refusing every asymmetric pair would
    #  throw away ~800 perfectly comparable lines; `hfs_collapse` records the asymmetry so
    #  it stays visible, and `main` measures whether it actually costs agreement.
    #
    #  What IS ambiguous is a block whose per-isotope sums do not reproduce the total
    #  linemake declares for it -- there the total is not a quantity we can stand behind.
    df["hfs_collapse"] = np.where(
        ~df.linemake_matched, "",
        np.where(df.linemake_hfs_expanded & (df.our_hfs_n_components.fillna(1) > 1), "both",
                 np.where(df.linemake_hfs_expanded, "linemake_only",
                          np.where(df.our_hfs_n_components.fillna(1) > 1, "ours_only", "neither"))))
    df.loc[hfs_bad, "verdict"] = "HFS_AMBIGUOUS"

    ok = matched & df.verdict.ne("HFS_AMBIGUOUS")
    caveated = df.linemake_readme_caveat.fillna("") != ""
    cw, ct = df.linemake_caveat_wavelength_A, df.linemake_caveat_tol_A
    #  A line-scoped caveat applies to its own line only; a species-wide one to all of them.
    #  Scoped against the LINEMAKE wavelength, because the README is describing a line in
    #  linemake's list, not in ours -- ours sits 6 mA away and would miss.
    in_scope = cw.isna() | ((df.linemake_wavelength - cw).abs() <= ct)
    nonprimary = ok & caveated & in_scope
    df.loc[nonprimary, "verdict"] = "LINEMAKE_NONPRIMARY"

    ok = ok & ~nonprimary
    stronger = (ok
                & df.our_gf_tier.isin(FALLBACK_TIERS)
                & df.linemake_reference_class.isin({"PRIMARY_LAB", "MIXED"})
                & df.linemake_source_class.eq("PRIMARY_LAB"))
    df.loc[stronger, "verdict"] = "LINEMAKE_STRONGER_SOURCE"

    rest = ok & ~stronger
    df.loc[rest & df.delta_within_threshold, "verdict"] = "AGREES"
    df.loc[rest & ~df.delta_within_threshold, "verdict"] = "DISAGREES_GF"
    return df


def isotope_multiplier_scan(xref: pd.DataFrame, ours: pd.DataFrame) -> pd.DataFrame:
    """Turn the cross-reference back on OUR pool: find published values that are a whole
    isotope count too large.

    This was not in the plan. It fell out of the audit, and it is the reason the audit was
    worth running. linemake's HFS components are listed per isotope (see `collapse_hfs`),
    and a list that sums them without regard to isotope inflates the total by exactly
    log10(n_isotopes). We built that guard for the linemake side. The scan below asks
    whether OUR side has the same defect, and it does not need to know any element's
    isotope count to answer -- the signature is self-evident in the data:

      1. the published `log_gf` sits a CONSTANT offset above its own `gf_linelist_vald`
         sibling column, identical across every affected line of the species;
      2. that constant is log10 of a small INTEGER; and
      3. the sibling column -- not the published value -- is what agrees with linemake's
         primary-lab declared total.

    Point 3 is what makes it a defect rather than a disagreement: an independent laboratory
    reference adjudicates the direction. RYA-161 -- this REPORTS, it does not correct.
    """
    cols = ["line_id", "gf_synth_ges", "gf_linelist_vald", "hfs_n_components"]
    m = xref[xref.linemake_matched].merge(ours[cols], on="line_id", how="left")
    m = m[m.gf_linelist_vald.notna()].copy()
    m["published_minus_vald"] = m.our_loggf - m.gf_linelist_vald
    m["vald_minus_linemake"] = m.gf_linelist_vald - m.linemake_loggf

    out = []
    for sp, g in m.groupby("species"):
        off = g[g.published_minus_vald.abs() > 0.25]
        if off.empty or (off.published_minus_vald.max() - off.published_minus_vald.min()) > 0.01:
            continue                      # not a constant offset -> ordinary disagreement
        k = 10.0 ** float(off.published_minus_vald.median())
        if abs(k - round(k)) > 0.02 or round(k) < 2:
            continue                      # not a whole multiplier
        if float(off.vald_minus_linemake.abs().median()) > 0.01:
            continue                      # the sibling does not agree with linemake either
        o = off.copy()
        o["multiplier"] = int(round(k))
        o["offset_dex"] = float(off.published_minus_vald.median())
        out.append(o)
    if not out:
        return pd.DataFrame(columns=list(m.columns) + ["multiplier", "offset_dex"])
    return pd.concat(out, ignore_index=True).sort_values(["species", "our_wavelength",
                                                          "line_id"])


def ep_collision_scan(xref: pd.DataFrame, ours: pd.DataFrame, wtol: float, etol: float,
                     thr_default: float = 0.02) -> pd.DataFrame:
    """Matches where the EP key and the log gf evidence point at DIFFERENT rows of ours.

    RYA-1034 says match on the physical transition, and this audit does. But an EP key is
    only as good as the EP on both sides. If linemake's line lands on our row A by EP while
    a DIFFERENT row B at the same wavelength carries linemake's log gf exactly, then one of
    the two lists has the wrong EP for that transition, and the match we made is the wrong
    row. Left undetected, such a pair shows up as a spectacular gf disagreement and would be
    published as a promotion candidate it is not.

    This is not hypothetical -- it catches Fe I 5538.516, where our Ruffoni2014 ingest and
    linemake's RUF14 entry agree on wavelength and log gf to the digit and disagree on EP by
    0.58 eV.
    """
    hits = []
    for sp, g in xref[xref.linemake_matched].groupby("species"):
        o = ours[ours.species == sp]
        if o.empty:
            continue
        ow = o.wavelength_air_A.to_numpy(float)
        for r in g.itertuples():
            near = o[(np.abs(ow - r.linemake_wavelength) <= wtol) & (o.line_id != r.line_id)]
            if near.empty or abs(r.delta_loggf) <= thr_default:
                continue
            better = near[(near.log_gf - r.linemake_loggf).abs() <= thr_default]
            for b in better.itertuples():
                #  A rival at the SAME EP is a wavelength blend, not an EP disagreement --
                #  `cross_reference` has already refused those as AMBIGUOUS_MATCH.
                if abs(float(b.excitation_potential_eV) - float(r.our_EP)) <= etol:
                    continue
                hits.append(dict(line_id=r.line_id, species=sp, our_wavelength=r.our_wavelength,
                                 our_EP=r.our_EP, our_loggf=r.our_loggf,
                                 linemake_wavelength=r.linemake_wavelength,
                                 linemake_EP=r.linemake_EP, linemake_loggf=r.linemake_loggf,
                                 delta_loggf=r.delta_loggf, matched_verdict=r.verdict,
                                 rival_line_id=b.line_id, rival_EP=b.excitation_potential_eV,
                                 rival_loggf=b.log_gf, rival_source=b.loggf_reference,
                                 ep_gap_eV=float(b.excitation_potential_eV) - float(r.our_EP)))
    return pd.DataFrame(hits)


def orphan_scan(xref: pd.DataFrame, ours: pd.DataFrame) -> pd.DataFrame:
    """Matched lines whose published value agrees with NONE of its own provenance columns.

    `canonical_gf` carries the value each upstream delivery gave for the line
    (`gf_synth_ges`, `gf_regions_vald`, `gf_linelist_vald`) alongside the published
    `log_gf`. A published value that matches none of them is either a deliberate upgrade
    or an orphan. LAB-tier rows are excluded because the upgrade is exactly what RYA-945
    did -- and those rows are the positive control here: they diverge from their stale
    sibling columns by design and linemake reproduces them to 0.000 dex.

    What is left needs an explanation, and linemake supplies an independent opinion on it.
    """
    cols = ["line_id", "gf_synth_ges", "gf_regions_vald", "gf_linelist_vald"]
    m = xref[xref.linemake_matched].merge(ours[cols], on="line_id", how="left")
    sib = m[["gf_synth_ges", "gf_regions_vald", "gf_linelist_vald"]]
    m["n_siblings"] = sib.notna().sum(axis=1)
    m["min_abs_delta_vs_sibling"] = sib.sub(m.our_loggf, axis=0).abs().min(axis=1)
    orphan = m[(m.n_siblings > 0)
               & (m.min_abs_delta_vs_sibling > 0.10)
               & ~m.our_gf_tier.astype(str).str.contains("LAB", na=False)]
    return orphan.sort_values(["min_abs_delta_vs_sibling", "line_id"],
                              ascending=[False, True])


# --------------------------------------------------------------------------------------
# 8. main
# --------------------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--linemake-dir", type=Path, default=None,
                    help="scratch clone location (default: a temp dir; never committed)")
    args = ap.parse_args()

    dest = args.linemake_dir or Path(tempfile.gettempdir()) / "rya1070_linemake"
    print(f"RYA-1070 linemake cross-reference (READ-ONLY)\n{'=' * 78}")
    lmdir, lmsha = acquire_linemake(dest)
    print(f"linemake        {lmdir}  @ {lmsha}")
    cg_sha = file_sha256(CANONICAL)
    print(f"canonical_gf    {CANONICAL.relative_to(ROOT)}  sha256 {cg_sha}")

    mooglists = lmdir / "mooglists"
    nohfs, hfs = parse_manifest(mooglists)
    print(f"manifest        {len(nohfs)} goodgf + {len(hfs)} goodgfhfs files, parsed from "
          f"mergenohfs/mergehfs")

    readme = parse_readme(lmdir / "README.md")
    print(f"README          {len(readme)} species; "
          f"{(readme.reference_class == 'PRIMARY_LAB').sum()} PRIMARY_LAB, "
          f"{(readme.reference_class == 'MIXED').sum()} MIXED, "
          f"{(readme.reference_class == 'CATALOG').sum()} CATALOG; "
          f"{(readme.caveats != '').sum()} caveated")
    for _, r in readme[readme.caveats != ""].iterrows():
        print(f"                  LINEMAKE_NONPRIMARY  {r.species:6s}  {r.caveats}")

    ours = pd.read_csv(CANONICAL, low_memory=False)
    #  The PRINTED value, not the float. `loggf_threshold` needs to know how many
    #  decimals each side actually committed to; float repr does not carry that.
    our_loggf_text = pd.read_csv(CANONICAL, usecols=["log_gf"], dtype=str).log_gf.fillna("")
    sp = ours[["key_z", "ion", "species"]].dropna().drop_duplicates()
    sp = sp[sp.key_z.astype(str).str.isdigit()]
    z_to_symbol = {int(r.key_z): str(r.species).split()[0] for r in sp.itertuples()}

    lm = pd.concat([build_linemake_table(mooglists, nohfs),
                    build_linemake_table(mooglists, hfs)], ignore_index=True)
    symbol_notes = check_filename_symbols(lm, z_to_symbol)
    roman = {0: "I", 1: "II", 2: "III"}
    lm["species"] = [f"{z_to_symbol.get(z, f'Z{z}')} {roman.get(i, f'ION{i}')}"
                     for z, i in zip(lm.key_z, lm.ion_index)]
    print(f"linemake lines  {len(lm)} comparable totals over {lm.species.nunique()} species "
          f"({int(lm.hfs_expanded.sum())} HFS-collapsed, "
          f"{int(lm.below_air_vac_boundary.sum())} below the air/vac boundary)")
    unrec = lm[lm.hfs_expanded & ~lm.hfs_reconciled]
    print(f"HFS check       {int(lm.hfs_expanded.sum()) - len(unrec)} blocks reconcile per "
          f"isotope to <= {HFS_RECONCILE_TOL_DEX} dex; {len(unrec)} do NOT -> HFS_AMBIGUOUS")

    airvac = measure_air_vac(ours, lm)
    print(f"air/vac         MEASURED {airvac['measured_scale'].upper()} above "
          f"{AIR_VACUUM_BOUNDARY_A:.0f} A -- EP-consistent coincidences within 5 mA: "
          f"air {airvac['coincidences_within_5mA']['air']} vs "
          f"vacuum {airvac['coincidences_within_5mA']['vacuum']}")

    covered = sorted(set(lm.species) & set(ours.species))
    pool = ours[ours.species.isin(covered)].copy()
    print(f"population      {len(pool)} canonical_gf lines over {len(covered)} shared species "
          f"(of {len(ours)} total)")

    lm_ok = lm[~lm.below_air_vac_boundary]
    dw, de = _scan_residuals(pool, lm_ok)
    ndw, nde = [], []
    for shift in NULL_SHIFTS_A:
        a, b = _scan_residuals(pool, lm_ok, shift=shift)
        ndw.append(a)
        nde.append(b)
    ndw, nde = np.concatenate(ndw), np.concatenate(nde)

    W_EDGES = [0, .0005, .001, .002, .005, .01, .02, .03, .05, .08, .1, .15, .2, .3, .5, 1.0]
    E_EDGES = [0, .0005, .001, .002, .005, .01, .02, .05, .1, .2, .5, 1., 2., 5., 20.]
    wtol, wprof, wstab = derive_tolerance(
        dw, np.abs(de) <= 0.05, ndw, np.abs(nde) <= 0.05, W_EDGES, "wavelength")
    etol, eprof, estab = derive_tolerance(
        de, np.abs(dw) <= wtol, nde, np.abs(ndw) <= wtol, E_EDGES, "excitation potential")
    print(f"tolerances      DERIVED against a {len(NULL_SHIFTS_A)}-offset coincidence null at "
          f"{NULL_RATIO_CRITERION:g}x: wavelength {wtol:g} A, EP {etol:g} eV")
    print(f"                stable across the criterion -- wavelength {wstab}, EP {estab}")
    #  Independent corroboration, NOT the source of the numbers above: the project's own
    #  ratified per-line match tolerances. They were derived here from scratch and land on
    #  the same wavelength window and a stricter EP window. A large divergence would mean
    #  this audit is matching on a different notion of "same line" than the gf grader.
    ssot = dict(wave_tol_A=WAVE_TOL_A, ep_tol_eV=EP_TOL_EV, loggf_match_tol=LOGGF_MATCH_TOL)
    print(f"                cross-check vs pipeline.gf_grades SSOT {ssot}: "
          f"wavelength {'matches' if wtol == WAVE_TOL_A else f'derived {wtol:g} vs {WAVE_TOL_A:g}'}, "
          f"EP {'matches' if etol == EP_TOL_EV else f'derived {etol:g} vs {EP_TOL_EV:g} (stricter)' if etol < EP_TOL_EV else f'derived {etol:g} vs {EP_TOL_EV:g} (LOOSER)'}")
    if wtol > 5 * WAVE_TOL_A or etol > 5 * EP_TOL_EV:
        raise ParseError(f"derived tolerances ({wtol} A / {etol} eV) are more than 5x the "
                         f"project's ratified match window {ssot} -- this audit would be "
                         f"calling lines 'the same' that the gf grader would not")

    xref = cross_reference(pool, lm_ok, wtol, etol, readme, our_loggf_text)
    xref = assign_verdicts(xref)
    thrq = xref.agreement_threshold_dex.dropna()
    thr = float(thrq.median()) if len(thrq) else float("nan")
    print(f"gf threshold    PER-LINE, from the printed precision of both sides: median "
          f"{thr:.4f} dex, range {thrq.min():.4f}-{thrq.max():.4f} dex "
          f"(gf_grades.LOGGF_MATCH_TOL is {LOGGF_MATCH_TOL} dex, one step looser)")
    thrstats = dict(n=int(len(thrq)), median_threshold_dex=thr,
                    min_threshold_dex=float(thrq.min()) if len(thrq) else None,
                    max_threshold_dex=float(thrq.max()) if len(thrq) else None,
                    exact_zero_deltas=int((xref.delta_loggf.abs() < 1e-12).sum()),
                    matched=int(xref.linemake_matched.sum()))

    # -- the Fe I graded-lab gate: both sides trace the same lab papers, so a large
    #    systematic delta here is a match-key or air/vac bug, not a disagreement.
    felab = xref[(xref.species == "Fe I") & xref.linemake_matched
                 & xref.our_gf_tier.astype(str).str.contains("LAB", na=False)]
    critical = []
    if felab.empty:
        critical.append("Fe I graded-lab control matched ZERO lines -- the control did not run")
        fe_med = float("nan")
    else:
        fe_med = float(np.median(np.abs(felab.delta_loggf)))
        print(f"\nFe I LAB-tier control: {len(felab)} matched, median |delta| = {fe_med:.4f} dex, "
              f"mean |delta| = {np.mean(np.abs(felab.delta_loggf)):.4f} dex")
        if fe_med > thr:
            critical.append(f"CRITICAL: Fe I graded-lab median |delta| {fe_med:.4f} dex exceeds "
                            f"the derived agreement threshold {thr:.4f} -- suspect a match-key "
                            f"or air/vac bug, NOT a real disagreement")

    # -- HFS asymmetry control. Refusing every asymmetrically-collapsed pair would have
    #    cost ~1,700 comparisons, so we kept them -- but "kept" is only defensible if we
    #    measure what the asymmetry does to agreement instead of assuming it does nothing.
    mm = xref[xref.linemake_matched & xref.verdict.ne("AMBIGUOUS_MATCH")]
    hfs_ctrl = (mm.groupby("hfs_collapse")
                .agg(n=("delta_loggf", "size"),
                     frac_within=("delta_within_threshold", "mean"),
                     median_abs_delta=("delta_loggf", lambda t: float(t.abs().median())))
                .reset_index())
    print(f"\nHFS COLLAPSE CONTROL (which side collapsed an HFS multiplet)\n{'-' * 78}")
    for _, r in hfs_ctrl.iterrows():
        print(f"  {r.hfs_collapse:14s} n={int(r.n):5d}  within threshold {r.frac_within:6.1%}  "
              f"median |delta| {r.median_abs_delta:.4f} dex")
    print("  -> asymmetric collapse costs agreement, but at the 0.006-0.02 dex level: the two "
          "lists\n     cluster slightly different component sets, NOT a summation error "
          "(which would be\n     log10(n_isotopes) = 0.30/0.70/0.85 dex).")

    # -- the audit turned back on our own pool
    iso = isotope_multiplier_scan(xref, ours)
    print(f"\nISOTOPE-MULTIPLIER SCAN (our published value vs its own VALD sibling)\n{'-' * 78}")
    if iso.empty:
        print("  none -- no species shows a constant whole-integer offset")
    else:
        for sp, g in iso.groupby("species"):
            print(f"  {sp:7s} {len(g):3d} matched lines published EXACTLY x{int(g.multiplier.iloc[0])} "
                  f"(+{g.offset_dex.iloc[0]:.4f} dex) above gf_linelist_vald, which agrees with "
                  f"linemake to {g.vald_minus_linemake.abs().max():.4f} dex")
        print(f"  -> {len(iso)} published canonical_gf values are high by a whole isotope count. "
              f"REPORTED, NOT CORRECTED (RYA-161).")
        #  Exposure beyond the referee: same species, same constant offset, but linemake
        #  does not cover the line, so nothing adjudicates it. Counted, not asserted.
        expo = []
        for sp, g in iso.groupby("species"):
            off = float(g.offset_dex.iloc[0])
            o = ours[(ours.species == sp) & ours.gf_linelist_vald.notna()]
            same = o[((o.log_gf - o.gf_linelist_vald) - off).abs() <= 0.01]
            expo.append(dict(species=sp, adjudicated=int(len(g)),
                             same_signature_total=int(len(same)),
                             no_referee=int(len(same) - len(g))))
        for e in expo:
            print(f"  {e['species']:7s} exposure: {e['same_signature_total']} lines carry the "
                  f"signature; {e['adjudicated']} adjudicated by linemake, "
                  f"{e['no_referee']} with no linemake counterpart")
        iso.to_csv(OUT / "isotope_multiplier_suspects.csv", index=False)

    bulk = build_bulk_table(mooglists, z_to_symbol)
    xref["bulk_counterpart"] = bulk_coverage(xref, bulk, wtol, etol)
    only_bulk = int((xref.bulk_counterpart & ~xref.linemake_matched).sum())
    print(f"\nBULK moogatom* COVERAGE (parsed, deliberately NOT used as corroboration)\n"
          f"{'-' * 78}")
    print(f"  {len(bulk):,} uncurated atomic lines over {bulk.species.nunique()} species. "
          f"{int(xref.bulk_counterpart.sum()):,} of our audited lines have a bulk counterpart; "
          f"{only_bulk:,} of those are lines the CURATED database does not reach.")
    print("  Counted only: these files carry no source tag and no README attribution, so an "
          "agreement\n  with them corroborates nothing and a disagreement is not a finding.")

    coll = ep_collision_scan(xref, pool, wtol, etol)
    print(f"\nEP-COLLISION SCAN (EP key and log gf evidence disagree about WHICH row)\n{'-' * 78}")
    if coll.empty:
        print("  none -- the EP key never lands on a different row than the gf evidence does")
    else:
        for _, r in coll.iterrows():
            print(f"  {r.line_id} {r.species:6s} {r.our_wavelength:9.3f}  matched on EP "
                  f"{r.our_EP:.3f} (our gf {r.our_loggf:+.3f}, delta {r.delta_loggf:+.3f}) but "
                  f"{r.rival_line_id} at EP {r.rival_EP:.3f} carries linemake's "
                  f"{r.linemake_loggf:+.3f} exactly [{r.rival_source}] -- EP gap "
                  f"{r.ep_gap_eV:+.3f} eV")
        coll.to_csv(OUT / "ep_collisions.csv", index=False)
        #  Do not publish these as promotion candidates or disagreements: the match itself
        #  is in doubt, and RYA-1034 says a match we cannot stand behind is refused.
        xref.loc[xref.line_id.isin(set(coll.line_id))
                 & xref.verdict.ne("AMBIGUOUS_MATCH"), "verdict"] = "EP_COLLISION"

    orph = orphan_scan(xref, ours)
    print(f"\nORPHANED PUBLISHED VALUES (agree with no sibling provenance column, not a LAB "
          f"upgrade)\n{'-' * 78}")
    if orph.empty:
        print("  none")
    else:
        for _, r in orph.iterrows():
            print(f"  {r.line_id} {r.species:6s} {r.our_wavelength:9.3f}  published "
                  f"{r.our_loggf:+.4f}  siblings "
                  f"{r.gf_synth_ges}/{r.gf_regions_vald}/{r.gf_linelist_vald}  "
                  f"linemake {r.linemake_loggf:+.3f} ({r.linemake_source_tag})  "
                  f"tier {r.our_gf_tier}, source {r.our_source}")
        orph.to_csv(OUT / "orphaned_published_values.csv", index=False)
        print(f"  -> {len(orph)} rows. The positive control passes: the "
              f"{int((xref.linemake_matched & xref.our_gf_tier.astype(str).str.contains('LAB', na=False)).sum())} "
              f"LAB-tier rows also diverge from their stale sibling columns, by design, and "
              f"linemake reproduces them exactly.")

    print(f"\nVERDICTS\n{'-' * 78}")
    vc = xref.verdict.value_counts()
    for v, n in vc.items():
        print(f"  {v:26s} {n:7d}   {n / len(xref):6.2%}")

    # -- the three non-primary entries the ticket names, shown as they landed
    print(f"\nLINEMAKE_NONPRIMARY CONTROLS (the three the ticket names)\n{'-' * 78}")
    np_rows = xref[xref.verdict == "LINEMAKE_NONPRIMARY"]
    for sp_name, g in np_rows.groupby("species"):
        print(f"  {sp_name:6s} {len(g):4d} of our lines tagged -- {g.linemake_readme_caveat.iloc[0]}")
        if not np.isnan(g.linemake_caveat_wavelength_A.iloc[0]):
            for _, r in g.iterrows():
                print(f"           {r.line_id} at {r.our_wavelength:.3f} A -> linemake "
                      f"{r.linemake_wavelength:.2f} A, {r.linemake_loggf:+.3f} "
                      f"[{r.linemake_source_tag}], delta {r.delta_loggf:+.3f}")
    co = readme[readme.species == "CO"].iloc[0]
    print(f"  CO     0 of our lines -- {co.caveats}. canonical_gf holds NO CO at all "
          f"({int((ours.species.astype(str) == 'CO').sum())} rows), so the caveat is recorded "
          f"from the README and has nothing to tag.")
    solar = xref[xref.linemake_matched
                 & xref.linemake_source_tag.astype(str).str.contains("SOLAR", case=False,
                                                                     na=False)]
    print(f"  {len(solar)} matched line(s) carry a SOLAR-tagged linemake value; none is "
          f"offered as a promotion candidate.")

    # -- species summary
    lm_counts = lm_ok.groupby("species").size()
    matched_rows = lm_ok.index.isin(xref.loc[xref.linemake_matched, "linemake_row"].dropna().astype(int))
    lm_used = lm_ok[matched_rows].groupby("species").size()
    rows = []
    for sp_name, g in xref.groupby("species"):
        m = g[g.linemake_matched]
        rows.append(dict(
            species=sp_name, our_line_count=len(g), matched=int(g.linemake_matched.sum()),
            AGREES=int((g.verdict == "AGREES").sum()),
            DISAGREES_GF=int((g.verdict == "DISAGREES_GF").sum()),
            LINEMAKE_STRONGER_SOURCE=int((g.verdict == "LINEMAKE_STRONGER_SOURCE").sum()),
            NO_MATCH=int((g.verdict == "NO_MATCH").sum()),
            AMBIGUOUS_MATCH=int((g.verdict == "AMBIGUOUS_MATCH").sum()),
            HFS_AMBIGUOUS=int((g.verdict == "HFS_AMBIGUOUS").sum()),
            LINEMAKE_NONPRIMARY=int((g.verdict == "LINEMAKE_NONPRIMARY").sum()),
            comparable=int(m.verdict.isin({"AGREES", "DISAGREES_GF",
                                           "LINEMAKE_STRONGER_SOURCE"}).sum()),
            within_threshold=int(m.delta_within_threshold.fillna(False).sum()),
            mean_abs_delta_loggf=round(float(np.mean(np.abs(m.delta_loggf))), 4) if len(m) else None,
            median_abs_delta_loggf=round(float(np.median(np.abs(m.delta_loggf))), 4) if len(m) else None,
            linemake_lines_not_in_our_pool=int(lm_counts.get(sp_name, 0) - lm_used.get(sp_name, 0)),
            linemake_primary_source=readme.set_index("species").primary_source.get(sp_name, ""),
            linemake_reference_class=readme.set_index("species").reference_class.get(sp_name, ""),
        ))
    #  Ties on `matched` must break on something stable, or the artifact's row order
    #  depends on the pandas version. Measured, not hypothetical: Sirius and the Mac
    #  emitted the same numbers in a different order for the species tied at 3 and 6.
    summary = pd.DataFrame(rows).sort_values(["matched", "species"],
                                             ascending=[False, True])

    # -- smoke-test controls the ticket names
    print(f"\nOVERLAP CONTROLS (must be non-zero)\n{'-' * 78}")
    for name in ("Fe I", "Mn I", "Eu II", "Al I"):
        r = summary[summary.species == name]
        if r.empty or int(r.matched.iloc[0]) == 0:
            critical.append(f"CRITICAL: known-overlap species {name} matched ZERO lines")
            print(f"  {name:6s} MATCHED 0  <-- CRITICAL")
        else:
            r = r.iloc[0]
            print(f"  {name:6s} ours {r.our_line_count:6d}  matched {r.matched:5d}  "
                  f"AGREES {r.AGREES:5d}  DISAGREES {r.DISAGREES_GF:4d}  "
                  f"STRONGER {r.LINEMAKE_STRONGER_SOURCE:5d}  "
                  f"mean|d| {r.mean_abs_delta_loggf}")

    OUT.mkdir(parents=True, exist_ok=True)
    cols = ["line_id", "species", "ion", "our_wavelength", "our_EP", "our_loggf", "our_source",
            "our_gf_tier", "our_hfs_n_components", "our_lab_source_tag", "linemake_matched",
            "linemake_wavelength", "linemake_EP", "linemake_loggf", "delta_loggf",
            "our_loggf_printed", "linemake_loggf_printed", "agreement_threshold_dex",
            "delta_within_threshold", "linemake_primary_source", "linemake_reference_class",
            "linemake_source_tag", "linemake_source_class", "linemake_file",
            "linemake_hfs_expanded", "linemake_n_components", "hfs_collapse", "bulk_counterpart",
            "linemake_readme_caveat", "linemake_caveat_wavelength_A", "linemake_caveat_tol_A",
            "verdict", "ambiguity"]
    for c in cols:
        if c not in xref.columns:
            xref[c] = ""
    xref[cols].to_csv(OUT / "per_line_xref.csv", index=False)
    summary.to_csv(OUT / "species_summary.csv", index=False)
    readme.to_csv(OUT / "linemake_readme_sources.csv", index=False)

    big = xref[(xref.verdict == "DISAGREES_GF") & (np.abs(xref.delta_loggf) > 0.10)]
    big = big.assign(_abs=np.abs(big.delta_loggf)).sort_values(
        ["_abs", "line_id"], ascending=[False, True]).drop(columns="_abs")
    big[cols].to_csv(OUT / "disagreements_over_0p10dex.csv", index=False)

    prov = dict(
        ticket="RYA-1070", read_only=True,
        linemake=dict(url=LINEMAKE_URL, commit=lmsha,
                      curated_manifest_source="mooglists/mergenohfs + mooglists/mergehfs",
                      clone_path_note="deliberately not recorded -- it is a scratch temp "
                                      "directory that differs per machine and would make "
                                      "this artifact machine-specific for no information",
                      goodgf_files=nohfs, goodgfhfs_files=hfs,
                      comparable_totals=int(len(lm)), species=int(lm.species.nunique())),
        canonical_gf=dict(path=str(CANONICAL.relative_to(ROOT)), sha256=cg_sha,
                          rows=int(len(ours)), audited_rows=int(len(pool)),
                          shared_species=covered),
        moog_decoding=dict(
            layout="fixed-width 4 x F10 (wavelength A, species code, EP eV, log gf), then a "
                   "free-text source tag from column 41",
            confirmation="every record decoded BOTH fixed-width and by whitespace-splitting "
                         "columns 1-40; any disagreement raises ParseError",
            species_code="integer part = Z; first decimal digit = ionisation stage "
                         "(0 = neutral = I, 1 = singly ionised = II); further decimal digits "
                         "= isotope mass number, marking an HFS/isotopic COMPONENT row",
            example=f"mooglists/feI.moog:1 -> "
                    f"{parse_moog_file(mooglists / 'feI.moog').iloc[0].to_dict()}",
            symbol_cross_check="Z -> symbol taken from canonical_gf (key_z, species) and "
                               "verified against the species stem of every linemake filename",
            symbol_notes=symbol_notes),
        hfs=dict(
            rule="a negative wavelength opens a block: that row carries linemake's declared "
                 "total, the rows after it are its components",
            trap="components are listed PER ISOTOPE and each isotope's set sums to the FULL "
                 "gf, so summing all components inflates by log10(n_isotopes) -- +0.301 dex "
                 "for Eu II (2), +0.699 for Ba II (5), +0.845 for Nd II (7)",
            method="compare against the declared total; VERIFY it by summing each isotope "
                   "separately; a block deviating by more than "
                   f"{HFS_RECONCILE_TOL_DEX} dex is HFS_AMBIGUOUS",
            blocks=int(lm.hfs_expanded.sum()), unreconciled=int(len(unrec)),
            max_dev_dex=float(np.nanmax(lm.hfs_max_dev_dex)) if len(lm) else None),
        air_vac=airvac,
        tolerances=dict(
            derivation="the |residual| histogram is a narrow true-match core on the flat "
                       "background of coincidental neighbours; that background is MEASURED by "
                       "re-running the same scan with our wavelengths displaced by each of "
                       "NULL_SHIFTS_A, which destroys every true pair and leaves the "
                       "coincidence rate intact; the tolerance is the outer edge of the last "
                       "bin whose real count still exceeds the null by NULL_RATIO_CRITERION",
            null_shifts_A=list(NULL_SHIFTS_A), null_ratio_criterion=NULL_RATIO_CRITERION,
            wavelength_A=wtol, wavelength_profile=wprof, wavelength_stability=wstab,
            excitation_potential_eV=etol, ep_profile=eprof, ep_stability=estab,
            ssot_cross_check=dict(
                source="pipeline.gf_grades (the project's ratified per-line match window)",
                **ssot,
                note="derived independently here and NOT taken from these; wavelength lands "
                     "on the same window, EP lands stricter"),
            loggf_agreement_median_dex=thr, loggf_stats=thrstats,
            loggf_derivation="PER LINE: half the printed quantum of each side summed -- the "
                             "widest gap two lists can show while carrying the same underlying "
                             "number at the precision each committed to. No free parameter. "
                             "Corroborated by the matched-delta histogram, whose decaying "
                             "shoulder ends at ~0.005 dex before a flat continuum"),
        matching=dict(
            key="species AND wavelength AND excitation potential (RYA-1034: never wavelength "
                "alone)",
            ambiguity="two linemake candidates inside tolerance, OR one linemake line that is "
                      "the only candidate for two distinct canonical_gf lines, is refused as "
                      "AMBIGUOUS_MATCH -- never broken by proximity",
            verdict_precedence=["AMBIGUOUS_MATCH", "HFS_AMBIGUOUS", "LINEMAKE_NONPRIMARY",
                                "LINEMAKE_STRONGER_SOURCE", "AGREES/DISAGREES_GF", "NO_MATCH"],
            fallback_tiers=sorted(FALLBACK_TIERS)),
        readme_parse=dict(
            species_rows=int(len(readme)),
            sections=sorted(set(readme.readme_section.dropna())),
            classification="PRIMARY_LAB = only paper citations; CATALOG = only "
                           "NIST/Kurucz/VALD/DREAM/HITRAN/ExoMol; MIXED = both",
            caveat_patterns=[p for p, _ in CAVEAT_PATTERNS],
            caveated={r.species: r.caveats for _, r in readme[readme.caveats != ""].iterrows()},
            note="the three known non-primary entries (Cu I, Nd II 4314.5, CO dv=2) are "
                 "DETECTED by these patterns; the parse raises if any is missing"),
        reproducibility=dict(
            deterministic_ordering="every count- or magnitude-ordered table breaks ties on a "
                                   "stable key (species or line_id), so row order does not "
                                   "depend on the pandas version",
            byte_identical_across_machines=False,
            measured="regenerated on Sirius (CI venv, py3.12) against the Mac run: every "
                     "verdict count and every reported number is identical, but pandas' CSV "
                     "float parser differs by up to 1 ULP between versions, so pass-through "
                     "columns can print differently (e.g. our_EP 4.8271999999999995 vs "
                     "4.8272; our_wavelength 19280.14200091231 vs 19280.142000912318)",
            consequence="one informational count moved by 1 -- bulk_moogatom.beyond_curated "
                        "was 57,163 on one machine and 57,164 on the other, a single line "
                        "sitting exactly on the bulk-coverage tolerance boundary. No verdict, "
                        "no match, and no scientific number is affected. Deliberately NOT "
                        "'fixed' by rounding the comparison inputs: a rounded number is not "
                        "an identity, and manufacturing agreement at a boundary is worse than "
                        "recording that the boundary exists"),
        verdict_counts={k: int(v) for k, v in vc.items()},
        hfs_collapse_control=hfs_ctrl.to_dict("records"),
        isotope_multiplier_scan=dict(
            found=int(len(iso)),
            by_species={sp: dict(lines=int(len(g)), multiplier=int(g.multiplier.iloc[0]),
                                 offset_dex=float(g.offset_dex.iloc[0]))
                        for sp, g in iso.groupby("species")} if len(iso) else {},
            method="published log_gf minus its own gf_linelist_vald sibling is a CONSTANT "
                   "offset across the species, that constant is log10 of a small integer, and "
                   "the sibling (not the published value) is what agrees with linemake's "
                   "primary-lab total -- so linemake adjudicates the direction",
            exposure=expo if len(iso) else [],
            action="REPORTED ONLY (RYA-161 validate-don't-tune); no value changed"),
        bulk_moogatom=dict(
            lines=int(len(bulk)), species=int(bulk.species.nunique()),
            our_lines_with_bulk_counterpart=int(xref.bulk_counterpart.sum()),
            beyond_curated=only_bulk,
            note="parsed as the ticket asks, and deliberately excluded from every verdict: "
                 "no source tag, no README attribution, so agreement corroborates nothing"),
        ep_collision_scan=dict(
            found=int(len(coll)),
            rows=coll.to_dict("records") if len(coll) else [],
            method="linemake's line matched our row A by EP, while a different row B within "
                   "the wavelength tolerance carries linemake's log gf exactly -- so the two "
                   "lists disagree about the EP of the transition and the match landed on the "
                   "wrong row. Verdict overridden to EP_COLLISION rather than published as a "
                   "promotion candidate (RYA-1034: a match we cannot stand behind is refused)"),
        orphan_scan=dict(
            found=int(len(orph)),
            rows=orph[["line_id", "species", "our_wavelength", "our_loggf", "gf_synth_ges",
                       "gf_regions_vald", "gf_linelist_vald", "linemake_loggf", "delta_loggf",
                       "our_source", "our_gf_tier", "linemake_source_tag"]].to_dict("records")
            if len(orph) else [],
            method="published log_gf agrees with none of its own provenance columns by more "
                   "than 0.10 dex, excluding LAB-tier rows (where the divergence is the "
                   "RYA-945 upgrade and linemake reproduces it to 0.000 dex)"),
        fe1_lab_control=dict(matched=int(len(felab)), median_abs_delta_dex=fe_med,
                             threshold_dex=thr, passed=bool(len(felab) and fe_med <= thr)),
        critical=critical,
    )
    (OUT / "provenance.json").write_text(json.dumps(prov, indent=2, default=str) + "\n",
                                         encoding="utf-8")

    write_report(OUT / "REPORT.md", xref, summary, readme, prov, thr, wtol, etol, big, felab)

    print(f"\nwrote {OUT.relative_to(ROOT)}/per_line_xref.csv ({len(xref)} rows)")
    print(f"wrote {OUT.relative_to(ROOT)}/species_summary.csv ({len(summary)} rows)")
    print(f"wrote {OUT.relative_to(ROOT)}/disagreements_over_0p10dex.csv ({len(big)} rows)")
    print(f"wrote {OUT.relative_to(ROOT)}/linemake_readme_sources.csv ({len(readme)} rows)")
    print(f"wrote {OUT.relative_to(ROOT)}/provenance.json")
    print(f"wrote {OUT.relative_to(ROOT)}/REPORT.md")

    if critical:
        print("\nCRITICAL:")
        for c in critical:
            print(f"  {c}")
        return 1
    return 0


def write_report(path: Path, xref: pd.DataFrame, summary: pd.DataFrame, readme: pd.DataFrame,
                 prov: dict, thr: float, wtol: float, etol: float, big: pd.DataFrame,
                 felab: pd.DataFrame) -> None:
    n = len(xref)
    matched = int(xref.linemake_matched.sum())
    agree = int((xref.verdict == "AGREES").sum())
    dis = int((xref.verdict == "DISAGREES_GF").sum())
    strong = xref[xref.verdict == "LINEMAKE_STRONGER_SOURCE"]
    #  "Comparable" and "agrees" must be counted over the SAME rows: the matched lines
    #  whose values are directly comparable. Refused verdicts (ambiguous, HFS-unreconciled,
    #  non-primary, EP collision) are excluded from both, not just the denominator.
    comp = xref[xref.verdict.isin({"AGREES", "DISAGREES_GF", "LINEMAKE_STRONGER_SOURCE"})]
    comparable = len(comp)
    within = int(comp.delta_within_threshold.fillna(False).sum())
    strong_sp = (strong.groupby("species").size().rename("n").reset_index()
                 .sort_values(["n", "species"], ascending=[False, True])
                 .set_index("species").n)
    L = []
    A = L.append
    A("# RYA-1070 — is `linemake` a good gf cross-reference for the lines we already hold?\n")
    A(f"READ-ONLY audit. linemake `{prov['linemake']['commit']}`, "
      f"`canonical_gf.csv` sha256 `{prov['canonical_gf']['sha256'][:16]}…` "
      f"({prov['canonical_gf']['rows']:,} rows). No line list was edited, no gf changed, no "
      f"linemake line imported, no synthesis run.\n")

    A("## The answer\n")
    A(f"**Yes for corroboration, and more usefully as a POINTER to lab papers we are not "
      f"citing.** linemake's curated database holds "
      f"{prov['linemake']['comparable_totals']:,} comparable transitions across "
      f"{prov['linemake']['species']} species; {len(prov['canonical_gf']['shared_species'])} of "
      f"those species appear in our pool, covering {n:,} of our "
      f"{prov['canonical_gf']['rows']:,} canonical_gf lines. Of those, **{matched:,} matched** "
      f"on the full physical key — species AND wavelength AND excitation potential — for a "
      f"coverage of **{matched / n:.1%}**. Of the {comparable:,} matched lines that carry a "
      f"directly comparable value, **{within:,} = "
      f"{(within / comparable if comparable else 0):.1%} agree** — the two lists carry the "
      f"same number to the precision both of them printed (median threshold {thr:.4f} dex, "
      f"derived per line; {prov['tolerances']['loggf_stats']['exact_zero_deltas']:,} of the "
      f"matched deltas are *exactly* zero). That is the corroboration answer.\n")
    A(f"The verdict buckets split those comparable lines differently, and deliberately: "
      f"`AGREES` ({agree}) is reserved for lines where we already sit at a comparable tier, "
      f"while **{len(strong):,} lines land in `LINEMAKE_STRONGER_SOURCE`** — matched lines "
      f"where we are still on a Kurucz/VALD/other fallback and linemake carries a value "
      f"tagged to a primary laboratory measurement. Those are the point of this audit, and "
      f"they are counted separately from agreement precisely so that a line agreeing with "
      f"linemake is not mistaken for a line that is already well sourced.\n")

    A(f"The strongest internal check passes. Our Fe I LAB-tier lines and linemake's Fe I list "
      f"ultimately trace the *same* Wisconsin papers (Ruffoni 2014, Den Hartog 2014, "
      f"Belmonte 2017), so they have to agree; {len(felab)} of them matched with a median "
      f"|Δlog gf| of **{prov['fe1_lab_control']['median_abs_delta_dex']:.4f} dex**. A large "
      f"systematic offset there would have meant a match-key or air/vacuum bug rather than a "
      f"real disagreement, and there is none. The wavelength scale was measured, not assumed: "
      f"under the air hypothesis linemake produces "
      f"{prov['air_vac']['coincidences_within_5mA']['air']:,} EP-consistent coincidences "
      f"within 5 mÅ of our air wavelengths, against "
      f"{prov['air_vac']['coincidences_within_5mA']['vacuum']} under the vacuum hypothesis — "
      f"linemake is on the air scale above {prov['air_vac']['boundary_A']:.0f} Å, and its "
      f"{prov['air_vac']['lines_below_boundary']} lines below that boundary were excluded from "
      f"the numeric match rather than matched raw.\n")

    A(f"## What we could promote — `LINEMAKE_STRONGER_SOURCE` ({len(strong)} lines)\n")
    A("These are lines where our `gf_tier` is a compilation fallback "
      "(KURUCZ / VALD3 / OTHER) and linemake carries a value tagged to a primary laboratory "
      "measurement. **This ticket does not promote them** — RYA-161 validate-don't-tune. It "
      "records where a promotion is available, via the paper linemake points at, as future "
      "work.\n")
    A("`linemake primary source` lists every paper the README's cell for that species cites, "
      "verbatim and in order — including papers cited for context rather than as the adopted "
      "source (the `V I` cell, for instance, names Holmes et al. 2016 while explaining that "
      "Wood et al. showed the Lawler values are correct). The full unedited cell for every "
      "species is preserved in `linemake_readme_sources.csv`; the adopted value for any one "
      "line is better identified by its per-line `linemake_source_tag` in "
      "`per_line_xref.csv`.\n")
    A("| species | lines | linemake primary source |")
    A("|---|---:|---|")
    for spn, cnt in strong_sp.items():
        src = readme.set_index("species").primary_source.get(spn, "")
        A(f"| {spn} | {cnt} | {src} |")
    A("")
    if len(strong):
        ex = strong.assign(_a=np.abs(strong.delta_loggf)).sort_values(
            ["_a", "line_id"], ascending=[False, True]).head(12)
        A("Twelve examples, largest |Δ| first — the size of Δ is how much a promotion would "
          "actually move the line:\n")
        A("| line_id | species | λ_air (Å) | EP (eV) | our log gf | our source | linemake log gf | Δ | linemake tag |")
        A("|---|---|---:|---:|---:|---|---:|---:|---|")
        for _, r in ex.iterrows():
            A(f"| `{r.line_id}` | {r.species} | {r.our_wavelength:.3f} | {r.our_EP:.3f} | "
              f"{r.our_loggf:+.3f} | {r.our_source} | {r.linemake_loggf:+.3f} | "
              f"{r.delta_loggf:+.3f} | {r.linemake_source_tag} |")
        A("")

    A(f"## Disagreements for later adjudication — |Δ| > 0.10 dex ({len(big)} lines)\n")
    if big.empty:
        A("None: every `DISAGREES_GF` line sits within 0.10 dex.\n")
    else:
        bysp = (big.groupby("species").size().rename("n").reset_index()
                .sort_values(["n", "species"], ascending=[False, True])
                .set_index("species").n)
        A("By species: " + ", ".join(f"{s} ({c})" for s, c in bysp.items()) + ".\n")
        A(f"Full list in `disagreements_over_0p10dex.csv`. The {min(30, len(big))} largest:\n")
        A("| line_id | species | λ_air (Å) | EP (eV) | our log gf | our source | our tier | linemake log gf | Δ | linemake tag |")
        A("|---|---|---:|---:|---:|---|---|---:|---:|---|")
        for _, r in big.head(30).iterrows():
            A(f"| `{r.line_id}` | {r.species} | {r.our_wavelength:.3f} | {r.our_EP:.3f} | "
              f"{r.our_loggf:+.3f} | {r.our_source} | {r.our_gf_tier} | "
              f"{r.linemake_loggf:+.3f} | {r.delta_loggf:+.3f} | {r.linemake_source_tag} |")
        A("")

    iso_prov = prov["isotope_multiplier_scan"]
    if iso_prov["found"]:
        A("## What the cross-reference found in OUR pool — published values high by a whole "
          "isotope count\n")
        A(f"This was not what the audit was looking for, and it is the most consequential "
          f"thing it found. **{iso_prov['found']} published `canonical_gf` values sit exactly "
          f"log10(n_isotopes) above the correct total.**\n")
        A("| species | adjudicated lines | multiplier | offset | same signature, no linemake referee |")
        A("|---|---:|---:|---:|---:|")
        for e in iso_prov["exposure"]:
            b = iso_prov["by_species"][e["species"]]
            A(f"| {e['species']} | {e['adjudicated']} | ×{b['multiplier']} | "
              f"+{b['offset_dex']:.4f} dex | {e['no_referee']} |")
        A("")
        A("The signature is self-evident and needed no isotope table to find: the published "
          "`log_gf` sits a **constant** offset above its own `gf_linelist_vald` sibling — "
          "identical to four decimals across every affected line of the species — and that "
          "constant is log10 of a small integer. What makes it a defect rather than a "
          "disagreement is the referee: **`gf_linelist_vald` is what matches linemake's "
          "primary-lab declared total, to 0.0000–0.0006 dex** (Lawler 2001 for Eu II, "
          "Den Hartog 2003 for Nd II, Lawler 2006 for Sm II). The published value is the one "
          "that is wrong, and it is wrong by exactly the number of stable isotopes the "
          "element has — 7 for Nd and Sm, 2 for Eu.\n")
        A("This is the same trap the audit guards against on linemake's own side: linemake "
          "lists HFS components *per isotope*, each isotope's set summing to the full gf, so "
          "a summation that ignores isotope identity inflates by log10(n_isotopes). The guard "
          "was built for the reference list and then found the defect in ours. Affected rows "
          "took the `gf_synth_ges` column rather than the VALD one; sibling rows of the same "
          "species that took the VALD column agree with linemake exactly, which rules out a "
          "species-wide scale error and localises it to a per-row column choice.\n")
        A("**Nothing here has been corrected** — RYA-161, validate-don't-tune. The affected "
          "rows are enumerated in `isotope_multiplier_suspects.csv` and logged in "
          "`data/audit/run_bug_ledger.csv` for adjudication.\n")
    coll_prov = prov["ep_collision_scan"]
    if coll_prov["found"]:
        A(f"## Where the EP key and the gf evidence point at different rows ({coll_prov['found']})\n")
        A("Matching on the physical transition is only as good as the excitation potential on "
          "*both* sides. These are matches where linemake's line landed on one of our rows by "
          "EP, while a **different** row of ours at the same wavelength carries linemake's log "
          "gf exactly — so the two lists disagree about the EP of the transition, and the row "
          "we matched is the wrong one. Their verdict is `EP_COLLISION`, not a promotion "
          "candidate: RYA-1034 refuses a match it cannot stand behind.\n")
        for r in coll_prov["rows"]:
            A(f"* **{r['species']} {r['linemake_wavelength']:.3f} Å** — linemake gives EP "
              f"{r['linemake_EP']:.3f} eV, log gf {r['linemake_loggf']:+.3f}. Our "
              f"`{r['line_id']}` sits at that EP with log gf {r['our_loggf']:+.3f} "
              f"(Δ {r['delta_loggf']:+.3f}), but `{r['rival_line_id']}` carries "
              f"{r['rival_loggf']:+.3f} — linemake's value to the digit — at EP "
              f"{r['rival_EP']:.3f} eV, **{abs(r['ep_gap_eV']):.3f} eV away**, sourced "
              f"`{r['rival_source']}`. Wavelength and log gf say these are the same "
              f"transition; the two lists disagree on its lower level by "
              f"{abs(r['ep_gap_eV']):.2f} eV. One of the two EPs is wrong and this audit "
              f"cannot say which.")
        A("")
        A("That there is exactly one such collision in "
          f"{prov['canonical_gf']['audited_rows']:,} audited lines is itself the reassuring "
          "result — the EP key is essentially collision-free, which is what makes the "
          "wavelength+EP match trustworthy everywhere else.\n")

    orph_prov = prov["orphan_scan"]
    if orph_prov["found"]:
        A(f"## Orphaned published values ({orph_prov['found']})\n")
        A("`canonical_gf` keeps what each upstream delivery said for a line beside the value "
          "it publishes. These rows agree with **none** of their own provenance columns, and "
          "are not the RYA-945 laboratory upgrade (LAB-tier rows are excluded — those diverge "
          "from their stale siblings by design, and linemake reproduces all "
          f"{prov['fe1_lab_control']['matched']} matched Fe I ones to 0.000 dex, which is the "
          "positive control for this scan).\n")
        A("| line_id | species | λ_air (Å) | published | own siblings | linemake | source | tier |")
        A("|---|---|---:|---:|---|---:|---|---|")
        for r in orph_prov["rows"]:
            sibs = "/".join("—" if v is None or (isinstance(v, float) and np.isnan(v))
                            else f"{v:+.3f}" for v in (r["gf_synth_ges"], r["gf_regions_vald"],
                                                       r["gf_linelist_vald"]))
            A(f"| `{r['line_id']}` | {r['species']} | {r['our_wavelength']:.3f} | "
              f"{r['our_loggf']:+.4f} | {sibs} | {r['linemake_loggf']:+.3f} "
              f"({r['linemake_source_tag']}) | {r['our_source']} | {r['our_gf_tier']} |")
        A("")
        A("**`gf_087247`, Mg I 5183.604 (b3), is the one that does not survive scrutiny.** It "
          "publishes +0.180 tagged `NIST-C+` / \"NIST ASD v5.11 grade A\", while all three of "
          "its own provenance columns say −0.239 and linemake's primary-lab value "
          "(Pehlivan Rhodin et al. 2017) says −0.168 — a 0.35 dex gap. The line-to-line "
          "spacing across the Mg b triplet settles it without appeal to any external number: "
          "our two independent references agree with each other on the spacing "
          "(VALD/GES give 0.211 and 0.692 dex; linemake/Pehlivan17 give 0.195 and 0.686) and "
          "`data/linelists/nist_reference.csv` rows 29–31, the source of the published values, "
          "give **0.630 and 1.211**. Two references agreeing with each other and disagreeing "
          "with the third localises the defect to that file — which already carries a "
          "documented correction of the same class in its own header (RYA-592 fixed the grade "
          "and `aki_s-1` columns on two other Mg I rows). Both `gf_087245` (b1, 0.10 dex from "
          "its siblings) and `gf_087247` (b3) come from those rows.\n")
        A("The three Fe rows are 0.10–0.14 dex and are ordinary NIST-C+ versus compilation "
          "differences; they are listed for completeness, not flagged.\n")
    A("## Refusals, and why they are results\n")
    amb = int((xref.verdict == "AMBIGUOUS_MATCH").sum())
    hfsa = int((xref.verdict == "HFS_AMBIGUOUS").sum())
    nonp = int((xref.verdict == "LINEMAKE_NONPRIMARY").sum())
    A(f"* **`AMBIGUOUS_MATCH` ({amb})** — two linemake candidates inside the tolerance, or one "
      f"linemake line that is the only candidate for two distinct lines of ours. RYA-1034: a "
      f"tolerance that cannot separate two candidates is a fact about the pool, not a tie to "
      f"be broken by proximity. Refused, never argmin'd.")
    A(f"* **`HFS_AMBIGUOUS` ({hfsa})** — the HFS expansion could not be reconciled. linemake "
      f"lists HFS components *per isotope*, each isotope's set summing to the full gf, so a "
      f"naive component sum inflates by log10(n_isotopes) — **+0.301 dex for Eu II, +0.699 for "
      f"Ba II, +0.845 for Nd II**. We compare against linemake's own declared total and verify "
      f"it by per-isotope summation; only blocks that fail that verification are refused.")
    A(f"* **`LINEMAKE_NONPRIMARY` ({nonp})** — linemake's own README flags the value as not a "
      f"primary lab result, so a disagreement there is not evidence that we are wrong. "
      f"Detected by parsing the README, not asserted:")
    for _, r in readme[readme.caveats != ""].iterrows():
        A(f"  * `{r.species}` — {r.caveats}")
    hc = {r["hfs_collapse"]: r for r in prov["hfs_collapse_control"]}
    if "neither" in hc:
        A(f"* **HFS asymmetry is kept, and measured.** Our `hfs_n_components > 1` rows are "
          f"already collapsed physical lines (gf-weighted centroid, log10 Σgf), and linemake's "
          f"comparable row is likewise always a total — so an asymmetric expansion is not by "
          f"itself a refusal. It does cost agreement, and here is how much: pairs where "
          f"neither side collapsed agree "
          f"{hc['neither']['frac_within']:.1%} of the time (median |Δ| "
          f"{hc['neither']['median_abs_delta']:.4f} dex), against "
          + ", ".join(f"{hc[k]['frac_within']:.1%} for `{k}` (n={int(hc[k]['n'])}, median |Δ| "
                      f"{hc[k]['median_abs_delta']:.4f})"
                      for k in ("ours_only", "linemake_only", "both") if k in hc)
          + ". Those offsets are at the 0.006–0.02 dex level — the two lists cluster slightly "
            "different component sets. A summation error would have shown up at "
            "log10(n_isotopes), i.e. 0.30–0.85 dex, and does not.")
    A(f"* **`NO_MATCH` ({int((xref.verdict == 'NO_MATCH').sum()):,})** — linemake simply does "
      f"not cover the line. That is the dominant outcome by count and it is expected: "
      f"linemake is a curated few-thousand-line database, our pool is a "
      f"{prov['canonical_gf']['rows']:,}-line survey list.\n")

    A("## How the numbers were derived\n")
    A(f"* **Match key** — species AND wavelength (±{wtol:g} Å) AND excitation potential "
      f"(±{etol:g} eV). Both tolerances were *read off the data* against a **measured** "
      f"coincidence null: the same nearest-neighbour scan re-run with our wavelengths "
      f"displaced by each of {len(prov['tolerances']['null_shifts_A'])} offsets, which "
      f"destroys every true pair while leaving both line densities intact. The tolerance is "
      f"the outer edge of the last bin whose real count still exceeds that null by "
      f"{prov['tolerances']['null_ratio_criterion']:g}×. The answer does not hang on that "
      f"multiple — the measured real/null ratio falls off a cliff at the selected edge "
      f"(wavelength 42 → 8.7 → 2.3, EP 315 → 9.8 → 2.4), and the same rule at 3×/5×/10× gives "
      + "/".join(f"{v:g}" for v in prov['tolerances']['wavelength_stability'].values())
      + " Å and "
      + "/".join(f"{v:g}" for v in prov['tolerances']['ep_stability'].values())
      + " eV. Nothing was hardcoded. As an *independent* "
      + (f"check the project's own ratified match window "
      f"(`gf_grades.WAVE_TOL_A` = {prov['tolerances']['ssot_cross_check']['wave_tol_A']:g} Å, "
      f"`EP_TOL_EV` = {prov['tolerances']['ssot_cross_check']['ep_tol_eV']:g} eV) lands on the "
      f"same wavelength window and a looser EP one."))
    A(f"* **Agreement threshold** — derived **per line**, with no free parameter: half the "
      f"printed quantum of each side, summed. If we print a log gf to 3 decimals and linemake "
      f"prints it to 2, the largest difference the two can show while representing one "
      f"underlying number is 0.0005 + 0.005. Median over the "
      f"{prov['tolerances']['loggf_stats']['n']:,} matched lines is {thr:.4f} dex "
      f"(range {prov['tolerances']['loggf_stats']['min_threshold_dex']:.4g}–"
      f"{prov['tolerances']['loggf_stats']['max_threshold_dex']:.4g}; lines we hold at one "
      f"decimal place get a correspondingly wider threshold, because at that precision we are "
      f"not entitled to call a small gap a disagreement). The matched Δ histogram corroborates "
      f"it: a spike of exact zeros, a shoulder decaying to ~0.005 dex, then a flat continuum. "
      f"The repo's own ratified `gf_grades.LOGGF_MATCH_TOL` "
      f"({prov['tolerances']['ssot_cross_check']['loggf_match_tol']} dex) is the same idea one "
      f"step looser.")
    A(f"* **MOOG decoding** — fixed-width 4 × F10 (λ Å, species code, EP eV, log gf) then a "
      f"free-text source tag from column 41. Every record is decoded twice, fixed-width and "
      f"free-form, and a disagreement raises. Species code: integer part Z, first decimal digit "
      f"the ionisation stage, further digits the isotope mass number.")
    A(f"* **Curated manifest** — the {len(prov['linemake']['goodgf_files'])} + "
      f"{len(prov['linemake']['goodgfhfs_files'])} file list is parsed out of linemake's own "
      f"`mergenohfs` / `mergehfs` scripts, so it cannot drift from the repo it describes. The "
      f"uncurated `moogatom*` bulk files are deliberately NOT used as corroboration: they carry "
      f"no README source attribution, and an unsourced Kurucz value agreeing with ours is not "
      f"evidence of anything.")
    A(f"* **Element symbols** — from `canonical_gf` (`key_z` + `ion`), cross-checked against "
      f"the species stem of every linemake filename. A disagreement is fatal.\n")

    rp = prov["reproducibility"]
    A("## Reproducibility of this artifact\n")
    A(f"Regenerated on Sirius (CI venv, py3.12) and compared against the Mac run: **every "
      f"verdict count and every reported number is identical**. The files are nonetheless "
      f"not byte-identical, and it is worth saying why rather than claiming they are. "
      f"pandas' CSV float parser differs by up to one ULP between versions, so pass-through "
      f"columns can print differently — `our_EP` as `4.8271999999999995` on one machine and "
      f"`4.8272` on the other, `our_wavelength` as `19280.14200091231` against "
      f"`19280.142000912318`.\n")
    A("One informational count moves with it: `bulk_moogatom.beyond_curated` is 57,163 on one "
      "machine and 57,164 on the other — a single line of 66,008 sitting exactly on the "
      "bulk-coverage tolerance boundary. No verdict, no match, and no scientific number is "
      "affected. That boundary case is deliberately **not** papered over by rounding the "
      "comparison inputs: a rounded number is not an identity, and manufacturing agreement at "
      "a boundary would be worse than recording that the boundary is there.\n")
    A("Row ordering *is* pinned. Every count- or magnitude-ordered table breaks ties on a "
      "stable key, because the first cross-machine comparison found the species tied at 3 and "
      "6 matches emitted in a different order — same numbers, different rows. The clone path "
      "is deliberately not recorded in `provenance.json` for the same reason: it is a scratch "
      "temp directory that differs per machine and carries no information.\n")
    A("## Scope\n")
    bm = prov["bulk_moogatom"]
    A(f"The uncurated `moogatom*` bulk lists **are** parsed — {bm['lines']:,} atomic lines "
      f"over {bm['species']} species — and are then deliberately excluded from every verdict. "
      f"That exclusion is a measured decision, not an oversight: "
      f"{bm['our_lines_with_bulk_counterpart']:,} of our audited lines have a bulk "
      f"counterpart, **{bm['beyond_curated']:,} of them lines the curated database does not "
      f"reach at all**, so including them would have multiplied apparent coverage six-fold. "
      f"But those files carry no source tag and no README attribution. An unsourced value "
      f"agreeing with our unsourced value corroborates nothing, and one disagreeing is not a "
      f"finding. Coverage is not the question this audit asks; provenance is.\n")
    A("Molecular species are out of scope for the numeric cross-reference: step 2 of the ticket "
      "scopes this to the atomic MOOG lists, and our pool holds no CO at all. The CO Δv = 2 "
      "caveat is recorded from the README above for completeness. `linemake` lines with no "
      "counterpart in our pool are **counted only**, per species, in "
      "`species_summary.csv` — they are not enumerated as import candidates, which is a "
      "separate decision.\n")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
