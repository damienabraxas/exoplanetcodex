"""RYA-1015 — the authoritative element x model-type availability matrix.

Ryan's call: STOP re-deriving grid/model availability in conversation. This module
builds ONE reconciled matrix from the sources that already exist, so "what atoms have
what models" is READ OFF A GRID, never re-litigated.

It PROMOTES rather than rebuilds. Four independent sources are reconciled:

  1. CSV claim   -- data/curation/nlte_grid_availability.csv        (RYA-462)
  2. 3D claim    -- data/curation/threednlte_availability.csv       (RYA-817)
  3. CODE truth  -- config.constants NLTE/THREED_CORRECTION_ELEMENTS (what actually runs)
  4. DISK truth  -- a Sirius `find -L` snapshot                      (RYA-1015 scan)

**Every disagreement becomes a loud PROBLEM cell carrying the triggering fact.** A
matrix built from any ONE source reproduces the loop this ticket exists to end — the
RYA-597 Ti drift (CSV said Bergemann2011, code ran Mallinson2024) is exactly what a
single-source matrix cannot see.

THE DISK HALF IS A SNAPSHOT, NOT A LIVE READ. Only Mr. Code can reach Sirius, and CI
cannot, so the scan is committed as a dated artifact. The snapshot records the
`find -L` control result; a snapshot whose control FAILED is refused (see
`load_disk_snapshot`) because a blind scan reports absences that are pure artifact.
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

#: The canonical 28 = 27 (incl. Fe II counted separately, RYA-109) + Zn (RYA-757).
#: See pipeline/element_freeze.py — the record key is (element, ion), not element.
CANONICAL_28: tuple[str, ...] = (
    "Li", "C", "N", "O", "Na", "Mg", "Al", "Si", "P", "S", "K", "Ca", "Sc", "Ti",
    "V", "Cr", "Mn", "Fe", "Fe II", "Co", "Ni", "Cu", "Zn", "Sr", "Y", "Zr", "Ba", "Eu",
)

#: Model-type axis (RYA-400 "The Beast": 1D-vs-3D x LTE-vs-NLTE).
MODEL_TYPES: tuple[str, ...] = ("1D_LTE", "1D_NLTE", "MEAN3D_NLTE", "FULL_3D_NLTE")

#: WHAT A CELL MEANS -- fixed here so it is never re-argued.
#:
#: A cell states the availability of a RUNNABLE OR APPLICABLE GRID/DECK of that model
#: type for that element. It is a CAPABILITY statement, not a literature statement.
#:
#: The distinction is load-bearing and it is the thing that keeps getting re-litigated:
#: a VENDORED POST-HOC CORRECTION TABLE (e.g. data/nlte_grids/amarsi2019_cno/ for C/O,
#: vendor/1L-3NErrors/ for Fe) lets us APPLY somebody's 3D result. It does NOT give us
#: a full-3D model we can run. Per RYA-1008 no public full-3D NLTE stellar RT code
#: exists at all, so FULL_3D_NLTE is NONE for every element -- while the vendored
#: correction that IS applied in production is recorded in the cell's `facts` so the
#: capability gap and the production reality are both visible.
#:
#: The RYA-817 CSV (threednlte_availability.csv) encodes LITERATURE availability, which
#: is a different question and therefore reports different values for the same element.
#: Both are kept; neither is silently resolved into the other.

#: Cell states.
HAVE = "HAVE"                  # CSV says + disk confirms + code uses
CODE_USES = "CODE_USES"        # the pipeline applies it now (code is ground truth)
CSV_ONLY = "CSV_ONLY"          # claimed but not on disk -> PROBLEM
DISK_ONLY = "DISK_ONLY"        # on disk but unregistered/unwired -> PROBLEM
REQUEST_ONLY = "REQUEST_ONLY"  # exists in literature, no public download
NONE = "NONE"                  # genuinely absent -> acquisition task
PROBLEM = "PROBLEM"            # sources disagree; `facts` carries why

#: RYA-1015 disk-parse traps. `CO` in a Gerber filename is the CO MOLECULE, not
#: cobalt, and `MN` is Mn shouted. Mapping CO -> Co would invent a cobalt NLTE grid
#: we do not have; this is a real defect the naive parse produces.
_FILENAME_ELEMENT_FIXUPS = {"MN": "Mn"}
_NOT_ELEMENTS = {"CO"}  # molecular decks, never an atomic element cell


@dataclass
class Cell:
    element: str
    model_type: str
    state: str
    csv_claim: str | None = None
    code_grid: str | None = None
    disk_paths: list[str] = field(default_factory=list)
    facts: list[str] = field(default_factory=list)

    @property
    def is_problem(self) -> bool:
        return self.state in (PROBLEM, CSV_ONLY, DISK_ONLY)

    def as_dict(self) -> dict:
        return {
            "element": self.element,
            "model_type": self.model_type,
            "state": self.state,
            "csv_claim": self.csv_claim,
            "code_grid": self.code_grid,
            "disk_paths": self.disk_paths,
            "facts": self.facts,
            "problem": self.is_problem,
        }


class DiskSnapshotError(RuntimeError):
    """The Sirius snapshot is missing, or its find -L control did not pass."""


def load_disk_snapshot(path: Path) -> dict[tuple[str, str], list[str]]:
    """Parse the committed Sirius `find -L` snapshot into {(element, model_type): paths}.

    LOUD-FAILS if the snapshot's positive control did not pass. A `find` that does not
    follow symlinks returns NOTHING for /srv/codex/grids (it is a symlink to the ntfs3
    drive), so an uncontrolled scan manufactures absences -- the RYA-1013 trap. We
    refuse to build a matrix on a blind scan rather than silently report NONE.
    """
    if not path.exists():
        raise DiskSnapshotError(
            f"Sirius disk snapshot missing: {path}. Regenerate with the RYA-1015 "
            f"scan (find -L + positive control) -- do NOT build the matrix without it."
        )
    text = path.read_text()
    if "CONTROL=PASS" not in text:
        raise DiskSnapshotError(
            f"{path}: positive control did not PASS. The scan was blind (find without "
            f"-L returns nothing through the /srv/codex/grids symlink), so every "
            f"absence in it is an artifact. Refusing to build the matrix."
        )
    out: dict[tuple[str, str], list[str]] = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("|")
        if len(parts) < 2:
            continue
        directory, fname = parts[0], parts[1]
        parsed = classify_disk_file(fname)
        if parsed is None:
            continue
        element, model_type = parsed
        out.setdefault((element, model_type), []).append(f"{directory}/{fname}")
    return out


def classify_disk_file(fname: str) -> tuple[str, str] | None:
    """Map a grid filename to (element, model_type), or None if it is not an atomic grid.

    Model atoms (atom.*) are supporting inputs, not availability cells, so they map to
    None -- an atom without a departure grid does not make an element NLTE-capable.
    """
    # Gerber TS-native <3D> deck: NLTEgrid[4TS]_<El>_STAGGERmean3D_<date>.bin
    m = re.match(r"NLTEgrid(?:4TS)?_([A-Za-z]{1,2})_STAGGERmean3D_", fname)
    if m:
        el = _normalise_element(m.group(1))
        return (el, "MEAN3D_NLTE") if el else None
    # Gerber TS-native 1D deck: NLTEgrid4TS_<El>_MARCS_<date>.bin
    m = re.match(r"NLTEgrid(?:4TS)?_([A-Za-z]{1,2})_MARCS_", fname)
    if m:
        el = _normalise_element(m.group(1))
        return (el, "1D_NLTE") if el else None
    # Amarsi GALAH PySME departure grid: nlte_<El>_*_pysme.grd
    m = re.match(r"nlte_([A-Za-z]{1,2})_.*pysme\.grd$", fname)
    if m:
        el = _normalise_element(m.group(1))
        return (el, "1D_NLTE") if el else None
    return None


def _normalise_element(token: str) -> str | None:
    if token in _NOT_ELEMENTS:
        return None
    if token in _FILENAME_ELEMENT_FIXUPS:
        return _FILENAME_ELEMENT_FIXUPS[token]
    return token[0].upper() + token[1:].lower() if len(token) == 2 else token.upper()


def load_csv_claims(path: Path) -> dict[str, list[dict]]:
    """Rows of the RYA-462 availability CSV, grouped by element."""
    out: dict[str, list[dict]] = {}
    with path.open() as fh:
        for row in csv.DictReader(fh):
            out.setdefault(row["element"].strip(), []).append(row)
    return out


def load_threed_claims(path: Path) -> dict[str, dict]:
    """Rows of the RYA-817 3D availability CSV, keyed by element."""
    out: dict[str, dict] = {}
    with path.open() as fh:
        for row in csv.DictReader(fh):
            out[row["element"].strip()] = row
    return out


def _code_truth() -> tuple[dict[str, dict], dict[str, dict]]:
    """What the pipeline ACTUALLY applies right now -- the ground truth."""
    import sys
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))
    from config import constants as C
    return dict(C.NLTE_CORRECTION_ELEMENTS), dict(C.THREED_CORRECTION_ELEMENTS)


def build_matrix(snapshot_path: Path | None = None) -> dict:
    """Reconcile CSV claim vs disk reality vs code usage into the full matrix."""
    curation = ROOT / "data" / "curation"
    snapshot_path = snapshot_path or (ROOT / "data" / "audit" / "rya1015"
                                      / "sirius_scan_raw.txt")
    disk = load_disk_snapshot(snapshot_path)
    csv_claims = load_csv_claims(curation / "nlte_grid_availability.csv")
    threed = load_threed_claims(curation / "threednlte_availability.csv")
    nlte_code, threed_code = _code_truth()

    cells: list[Cell] = []
    for element in CANONICAL_28:
        base = element.split()[0] if element.startswith("Fe") else element
        for mt in MODEL_TYPES:
            cells.append(_reconcile(element, base, mt, disk, csv_claims,
                                    threed, nlte_code, threed_code))

    problems = [c for c in cells if c.is_problem]
    return {
        "generated": date.today().isoformat(),
        "generator": "pipeline/model_availability_matrix.py (RYA-1015)",
        "snapshot": str(snapshot_path.relative_to(ROOT)),
        "sources": {
            "csv_claim": "data/curation/nlte_grid_availability.csv (RYA-462)",
            "threed_claim": "data/curation/threednlte_availability.csv (RYA-817)",
            "code_truth": "config.constants NLTE/THREED_CORRECTION_ELEMENTS",
            "disk_truth": f"{snapshot_path.name} (Sirius find -L, control PASS)",
        },
        "elements": list(CANONICAL_28),
        "model_types": list(MODEL_TYPES),
        "cells": [c.as_dict() for c in cells],
        "problem_count": len(problems),
        "problems": [c.as_dict() for c in problems],
    }


def _reconcile(element: str, base: str, mt: str, disk, csv_claims, threed,
               nlte_code, threed_code) -> Cell:
    cell = Cell(element=element, model_type=mt, state=NONE)
    disk_paths = disk.get((base, mt), [])
    cell.disk_paths = disk_paths

    if mt == "1D_LTE":
        # Every element is reachable by 1D-LTE synthesis; this is the baseline route,
        # not a grid. Stated explicitly so the column is not mistaken for a gap.
        cell.state = HAVE
        cell.facts.append("1D-LTE synthesis is the universal baseline route (no grid).")
        return cell

    if mt == "1D_NLTE":
        code_entry = nlte_code.get(base)
        # cno-3dnlte is a WIRED PRODUCTION ROUTE for C/N/O, not a registry entry --
        # omitting it made O read DISK_ONLY when O is in fact live in production.
        rows = [r for r in csv_claims.get(base, [])
                if r["subsystem"] in ("registry-nlte", "fe-nlte", "cno-3dnlte")]
        claim = rows[0]["grid_file"] if rows else None
        # Not every wired element goes through NLTE_CORRECTION_ELEMENTS: C/N/O run via
        # the cno-3dnlte subsystem and Fe via fe-nlte. Treating the registry as the only
        # wiring route reports live production elements as unwired DISK_ONLY.
        wired_elsewhere = [r for r in rows
                           if r["subsystem"] in ("cno-3dnlte", "fe-nlte")
                           and r["wired"].strip().lower() == "true"]
        if wired_elsewhere and not code_entry:
            cell.csv_claim = claim
            cell.state = HAVE
            cell.facts.append(
                f"Wired via the {wired_elsewhere[0]['subsystem']} subsystem "
                f"({claim}), not NLTE_CORRECTION_ELEMENTS. "
                f"{len(disk_paths)} departure grid(s) on Sirius.")
            return cell
        cell.csv_claim = claim
        cell.code_grid = code_entry.get("grid") if code_entry else None

        if code_entry and claim and claim != code_entry["grid"]:
            cell.state = PROBLEM
            cell.facts.append(
                f"DRIFT: CSV claims '{claim}' but code applies "
                f"'{code_entry['grid']}'. Code is ground truth.")
        elif code_entry:
            cell.state = HAVE if (claim or disk_paths) else CODE_USES
            if not disk_paths:
                cell.facts.append(
                    "Applied from the vendored CSV in data/nlte_grids/; the "
                    "departure grid itself lives on Sirius only.")
        elif disk_paths:
            cell.state = DISK_ONLY
            cell.facts.append(
                f"{len(disk_paths)} departure grid(s) on Sirius but NO code entry in "
                f"NLTE_CORRECTION_ELEMENTS -- unwired.")
        elif claim:
            cell.state = CSV_ONLY
            cell.facts.append(f"CSV claims '{claim}' but neither disk nor code has it.")
        else:
            cell.state = NONE
        return cell

    # --- 3D axis, from the RYA-817 CSV + disk ---
    row = threed.get(base, {})
    if mt == "MEAN3D_NLTE":
        offsolar = (row.get("offsolar_3d_nlte") or "").strip()
        solar = (row.get("solar_3d_nlte") or "").strip()
        if disk_paths:
            cell.state = DISK_ONLY
            cell.facts.append(
                f"<3D> STAGGERmean3D deck present on Sirius ({len(disk_paths)}) but no "
                f"code path consumes a <3D> departure deck -- unwired capability.")
        elif solar in ("FULL_3D_NLTE", "MEAN3D_NLTE") or offsolar == "GRID_MEAN3D":
            # A published solar 3D/<3D> treatment exists for this element, but we hold
            # no <3D> deck -> the corrections are obtainable only by request/from a
            # paper table (e.g. <3D> O = Amarsi 2016 Table 5, not a public download).
            cell.state = REQUEST_ONLY
            cell.facts.append(
                f"RYA-817 literature: solar={solar or '-'}, offsolar={offsolar or '-'}. "
                f"Published, but no <3D> deck on disk -> paper-table / request only.")
        else:
            cell.state = NONE
            if solar or offsolar:
                cell.facts.append(
                    f"RYA-817 literature: solar={solar or '-'}, offsolar={offsolar or '-'}.")
        return cell

    # FULL_3D_NLTE -- capability, not literature. See the MODEL_TYPES note above.
    solar = (row.get("solar_3d_nlte") or "").strip()
    offsolar = (row.get("offsolar_3d_nlte") or "").strip()
    holding = (row.get("our_holding") or "").strip()
    code_entry = threed_code.get(base)
    cell.code_grid = code_entry.get("grid") if code_entry else None

    # No full-3D deck exists on disk for any element, and RYA-1008 established that no
    # public full-3D NLTE stellar RT code exists to produce one. So this is NONE
    # everywhere -- but never SILENTLY: record what IS applied in production.
    cell.state = NONE
    if holding and holding.lower() != "none":
        cell.facts.append(
            f"NO runnable full-3D deck. A vendored POST-HOC CORRECTION is applied in "
            f"production from {holding} -- that applies someone else's 3D result, it "
            f"is not a full-3D capability (RYA-1008: no public full-3D NLTE RT code).")
    if solar or offsolar:
        cell.facts.append(
            f"RYA-817 literature: solar={solar or '-'}, offsolar={offsolar or '-'}.")
    return cell


def render_markdown(matrix: dict) -> str:
    """Compact element x model-type grid for humans."""
    by = {(c["element"], c["model_type"]): c for c in matrix["cells"]}
    mts = matrix["model_types"]
    out = [f"# Element x model availability (RYA-1015)", "",
           f"Generated {matrix['generated']} by `{matrix['generator']}`.",
           f"Disk half: `{matrix['snapshot']}` (Sirius `find -L`, control PASS).", "",
           "| element | " + " | ".join(mts) + " |",
           "|---|" + "---|" * len(mts)]
    for el in matrix["elements"]:
        row = [el]
        for mt in mts:
            c = by[(el, mt)]
            mark = "**" + c["state"] + "**" if c["problem"] else c["state"]
            row.append(mark)
        out.append("| " + " | ".join(row) + " |")
    out += ["", f"**PROBLEM cells: {matrix['problem_count']}**", ""]
    for c in matrix["problems"]:
        out.append(f"- `{c['element']}` / `{c['model_type']}` -> **{c['state']}**: "
                   + " ".join(c["facts"]))
    return "\n".join(out) + "\n"
