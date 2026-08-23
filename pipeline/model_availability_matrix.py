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



#: REPO-SIDE 3D HOLDINGS. The Sirius scan covers the departure decks only; the actual
#: 3D capability lives IN THE REPO and was invisible to a Sirius-only scan. This is
#: what we CAN RUN, keyed by the element it can correct.
#:
#: A matrix that reported FULL_3D_NLTE = NONE everywhere was WRONG: we hold a 3D-NLTE
#: Fe engine, 3D-NLTE C/O tables, and a 3D metals increment.
THREED_HOLDINGS: dict[str, dict] = {
    "Fe": {
        "path": "vendor/1L-3NErrors/",
        "kind": "FULL_3D_NLTE",
        "engine": "ENGINE-A-3DNLTE",
        "what": "Amarsi, Liljegren & Nissen 2022 (A&A 668 A68) 3D-NLTE Fe MLP "
                "(fe1_model_gt02.p / fe1_model_lt02.p / fe2_model.p)",
        "blocked_by": "RYA-923",
        "blocker": "URGENT/OPEN: the MLP returns NaN for EVERY in-domain line on main "
                   "(114 in-domain -> n=0). 1D-LTE legs still PASS, so only the "
                   "correction path regressed. Committed cells carry values from when "
                   "it worked (Fe I 7.604 n=114, Fe II 7.642 n=7) -- so the capability "
                   "is REAL but currently UNRUNNABLE.",
    },
    "C": {"path": "data/nlte_grids/amarsi2019_cno/", "kind": "FULL_3D_NLTE",
          "engine": "cno-3dnlte",
          "what": "Amarsi, Nissen & Skuladottir 2019 (A&A 630 A104) line-by-line "
                  "3D-NLTE / 1D-NLTE tables; 3D leg below Teff 6500 K"},
    "O": {"path": "data/nlte_grids/amarsi2019_cno/", "kind": "FULL_3D_NLTE",
          "engine": "cno-3dnlte",
          "what": "Amarsi 2019 3D-NLTE O I 777; [O I] 6300 is forbidden-LTE by "
                  "construction (RYA-447)"},
    "N": {"path": "data/nlte_grids/amarsi2019_cno/", "kind": "FULL_3D_NLTE",
          "engine": "cno-3dnlte",
          "what": "Amarsi 2019 CNO synthesis leg (N atomic departures are the separate "
                  "1D registry grid)"},
    "Si": {"path": "data/threed_grids/solar3d_metals_rya399.csv", "kind": "MEAN3D_NLTE",
           "engine": "THREED_CORRECTION_ELEMENTS",
           "what": "Amarsi & Asplund 2017 (MNRAS 464, 264) solar 3D increment (RYA-399)"},
    "Ti": {"path": "data/threed_grids/solar3d_metals_rya399.csv", "kind": "MEAN3D_NLTE",
           "engine": "THREED_CORRECTION_ELEMENTS",
           "what": "Scott et al. 2015 Paper II (A&A 573, A26) solar 3D Ti (RYA-399)"},
    "Cr": {"path": "data/threed_grids/solar3d_metals_rya399.csv", "kind": "MEAN3D_NLTE",
           "engine": "THREED_CORRECTION_ELEMENTS",
           "what": "Scott et al. 2015 Paper II (A&A 573, A26) solar 3D Cr (RYA-399)"},
}

#: The <3D> STAGGER solar atmosphere we hold -- the model any <3D> route needs.
STAGGER_MEAN3D_ATMOSPHERE = "data/atmospheres/stagger_avg3d_rya442/sun_avg3d_stagger.mod"


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
        if (h := THREED_HOLDINGS.get(base)) and h["kind"] == "MEAN3D_NLTE":
            cell.state = HAVE
            cell.code_grid = h["path"]
            cell.facts.append(f"CAN RUN via {h['engine']}: {h['what']} ({h['path']}).")
            if disk_paths:
                cell.facts.append(
                    f"ALSO holds {len(disk_paths)} unwired <3D> STAGGERmean3D deck(s) "
                    f"on Sirius -- a second, richer route nothing consumes yet.")
        elif disk_paths:
            cell.state = DISK_ONLY
            cell.facts.append(
                f"<3D> STAGGERmean3D deck present on Sirius ({len(disk_paths)}) but no "
                f"code path consumes a <3D> departure deck -- unwired capability.")
        elif False:
            cell.state = HAVE
            cell.code_grid = h["path"]
            cell.facts.append(f"CAN RUN via {h['engine']}: {h['what']} ({h['path']}).")
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

    # What we CAN RUN, from the repo-side holdings -- NOT a blanket NONE. We hold a
    # 3D-NLTE Fe MLP and 3D-NLTE C/N/O tables; a Sirius-only scan could not see them.
    h = THREED_HOLDINGS.get(base)
    if h and h["kind"] == "FULL_3D_NLTE":
        cell.code_grid = h["path"]
        if h.get("blocked_by"):
            cell.state = PROBLEM
            cell.facts.append(
                f"CAPABILITY EXISTS BUT IS BROKEN ({h['blocked_by']}): {h['blocker']}")
            cell.facts.append(f"Engine {h['engine']}: {h['what']} ({h['path']}).")
        else:
            cell.state = HAVE
            cell.facts.append(f"CAN RUN via {h['engine']}: {h['what']} ({h['path']}).")
    else:
        cell.state = NONE
        cell.facts.append(
            "No 3D-NLTE correction or engine we can run for this element. "
            "(We cannot COMPUTE full 3D from scratch for anything -- RYA-1008: no "
            "public full-3D NLTE RT code -- so every 3D capability here is a "
            "published grid/model we apply.)")
    if solar or offsolar:
        cell.facts.append(
            f"RYA-817 literature: solar={solar or '-'}, offsolar={offsolar or '-'}.")
    return cell


#: The two engines (pipeline/engine_selection.py). There is no "primary" — both are
#: products that get presented, and higher reach is BROADER, not better.
ENGINE_A = "Engine-A (EW + grid delta)"
ENGINE_B = "Engine-B (synthesis)"


def build_engine_matrix(matrix: dict) -> list[dict]:
    """Per element: what each ENGINE can actually run, and where its grid lives.

    Engine-A runs 1D-NLTE from a VENDORED departure CSV in data/nlte_grids/ (in-repo).
    Engine-B runs synthesis, and goes NLTE only when a TS-native Gerber deck is on
    Sirius. Either engine falls back to LTE, so "no grid" never means "no engine" --
    it means that engine is LTE-only for that element, which is the distinction that
    keeps getting lost.
    """
    import sys
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))
    from config import constants as C
    nlte_code = dict(C.NLTE_CORRECTION_ELEMENTS)

    by = {(c["element"], c["model_type"]): c for c in matrix["cells"]}
    grids_dir = ROOT / "data" / "nlte_grids"
    rows = []
    for element in CANONICAL_28:
        base = element.split()[0] if element.startswith("Fe") else element
        code_entry = nlte_code.get(base)
        onedee = by[(element, "1D_NLTE")]
        mean3d = by[(element, "MEAN3D_NLTE")]

        # --- Engine-A: needs a vendored CSV it can read from the repo ---
        vendored = code_entry.get("grid") if code_entry else onedee.get("csv_claim")
        a_present = bool(vendored) and (grids_dir / vendored).exists() \
            if vendored and vendored.endswith(".csv") else bool(vendored)
        if code_entry:
            a_mode, a_where = "1D-NLTE", f"data/nlte_grids/{code_entry['grid']}"
        elif onedee["state"] == HAVE and vendored:
            a_mode, a_where = "1D-NLTE", f"wired via subsystem ({vendored})"
        else:
            a_mode, a_where = "LTE only", "no departure grid"

        # --- Engine-B: NLTE only with a TS-native deck on Sirius ---
        ts_1d = [p for p in onedee["disk_paths"] if "gerber_ts" in p]
        ts_3d = [p for p in mean3d["disk_paths"] if "gerber_ts" in p]
        if ts_3d:
            b_mode = "<3D>-NLTE deck on disk (UNWIRED)"
            b_where = "; ".join(Path(p).name for p in ts_3d)
        elif ts_1d:
            b_mode = "1D-NLTE (TS-native)"
            b_where = "; ".join(Path(p).name for p in ts_1d)
        else:
            b_mode, b_where = "LTE only", "no TS-native deck"

        # Engine-A-3DNLTE / the 3D route: a published 3D model we APPLY.
        h3 = THREED_HOLDINGS.get(base)
        if h3 and h3.get("blocked_by"):
            c_mode = f"3D-NLTE BROKEN ({h3['blocked_by']})"
            c_where = h3["path"]
        elif h3:
            c_mode = ("3D-NLTE" if h3["kind"] == "FULL_3D_NLTE" else "<3D> increment")
            c_where = h3["path"]
        else:
            c_mode, c_where = "none", "no 3D model held"

        rows.append({
            "engine_c_mode": c_mode,
            "engine_c_where": c_where,
            "element": element,
            "engine_a_mode": a_mode,
            "engine_a_where": a_where,
            "engine_a_present": a_present,
            "engine_b_mode": b_mode,
            "engine_b_where": b_where,
            "two_engine": a_mode != "LTE only" and not b_mode.startswith("LTE"),
        })
    return rows


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


#: Molecules matter for the M-dwarf tier and for C/N/O, and they live in a DIFFERENT
#: place from the atomic grids, which is exactly why they get forgotten: the LTE
#: linelists are vendored in-repo while the only molecular NLTE deck is on Sirius.
MOLECULES_EXPECTED: tuple[str, ...] = (
    # carried today (C/N/O coupling, RYA-360)
    "C2", "CH", "CN", "CO", "NH", "OH",
    # NOT carried -- the M-dwarf / cool-star gap, listed so the hole is visible
    "TiO", "VO", "ZrO", "MgH", "FeH", "SiH", "CaH", "H2O",
)


def build_molecule_matrix(snapshot_path: Path | None = None) -> list[dict]:
    """What we hold per MOLECULE: vendored LTE linelist, and any NLTE deck.

    Absent molecules are listed explicitly rather than omitted -- an inventory that
    only shows what you have cannot show you a gap.
    """
    snapshot_path = snapshot_path or (ROOT / "data" / "audit" / "rya1015"
                                      / "sirius_scan_raw.txt")
    lists_dir = ROOT / "data" / "linelists" / "molecular" / "turbospectrum"
    have_lists = {p.name for p in lists_dir.iterdir() if p.is_dir()} \
        if lists_dir.exists() else set()

    decks: dict[str, list[str]] = {}
    for line in snapshot_path.read_text().splitlines():
        if line.startswith("#") or "|" not in line:
            continue
        _, fname = line.split("|")[0], line.split("|")[1]
        m = re.match(r"NLTEgrid(?:4TS)?_([A-Za-z0-9]{2,4})_MARCS_", fname)
        if m and m.group(1).upper() in {x.upper() for x in MOLECULES_EXPECTED}:
            decks.setdefault(m.group(1).upper(), []).append(fname)

    rows = []
    for mol in MOLECULES_EXPECTED:
        lte = mol in have_lists
        nlte = decks.get(mol.upper(), [])
        rows.append({
            "molecule": mol,
            "lte_linelist": ("data/linelists/molecular/turbospectrum/" + mol)
                            if lte else None,
            "nlte_deck": nlte[0] if nlte else None,
            "state": "HAVE" if lte else "NONE",
        })
    return rows


def _molecule_table(molecules: list[dict] | None) -> str:
    if not molecules:
        return ""
    rows = []
    for m in molecules:
        lte = (f'<span class="st" style="color:#5fd38d">HAVE</span>'
               f'<div class=g>{m["lte_linelist"]}</div>') if m["lte_linelist"] \
              else '<span class="st" style="color:#6a7690">NONE</span>'
        nlte = (f'<span class="st" style="color:#5fd38d">HAVE</span>'
                f'<div class=g>{m["nlte_deck"]}</div>') if m["nlte_deck"] \
               else '<span class="st" style="color:#6a7690">NONE</span>'
        cls = "" if m["lte_linelist"] else ' class="prob"'
        rows.append(f'<tr><th>{m["molecule"]}</th><td{cls}>{lte}</td><td>{nlte}</td></tr>')
    return ('<h2 style="font-size:1.15rem;margin:2rem 0 .5rem">Molecules</h2>'
            '<p class="sub">Molecular data lives apart from the atomic grids, which is '
            'why it goes missing. LTE linelists are vendored in-repo; the only molecular '
            'NLTE deck is on Sirius. Absent molecules are listed, not omitted &mdash; an '
            'inventory that shows only what you have cannot show a gap.</p>'
            '<table><thead><tr><th>molecule</th><th>LTE linelist</th>'
            '<th>NLTE deck</th></tr></thead><tbody>'
            + "".join(rows) + '</tbody></table>')


def render_html(matrix: dict, engine_rows: list[dict],
                molecules: list[dict] | None = None) -> str:
    """Self-contained page: every canonical element x model type, with REAL grid names.

    Static by design -- the live site is GitHub Pages, so a fetch()-driven page would
    show nothing. What you can see is the point of a tracker.
    """
    by = {(c["element"], c["model_type"]): c for c in matrix["cells"]}
    eng = {r["element"]: r for r in engine_rows}

    def td(el, mt):
        c = by[(el, mt)]
        cls = {"HAVE": "have", "REQUEST_ONLY": "req", "NONE": "none",
               "DISK_ONLY": "prob", "CSV_ONLY": "prob", "PROBLEM": "prob",
               "CODE_USES": "have"}.get(c["state"], "none")
        names = "<br>".join(Path(p).name for p in c["disk_paths"]) or ""
        claim = c["code_grid"] or c["csv_claim"] or ""
        detail = "<br>".join(x for x in (claim, names) if x)
        return (f'<td class="{cls}"><span class="st">{c["state"]}</span>'
                f'{"<div class=g>" + detail + "</div>" if detail else ""}</td>')

    rows = []
    for el in matrix["elements"]:
        e = eng[el]
        rows.append(
            f'<tr><th>{el}</th>'
            + "".join(td(el, mt) for mt in matrix["model_types"])
            + f'<td class="eng">A: {e["engine_a_mode"]}<br>B: {e["engine_b_mode"]}</td></tr>')

    heads = "".join(f"<th>{mt.replace('_', '-')}</th>" for mt in matrix["model_types"])
    return f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Element x model availability - The Exoplanet Codex</title>
<style>
:root{{color-scheme:dark}}
body{{background:#0a0e1a;color:#d8e0f0;font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;margin:0;padding:2rem}}
h1{{font-size:1.5rem;margin:0 0 .25rem}} .sub{{color:#8fa0c0;margin:0 0 1.5rem}}
table{{border-collapse:collapse;width:100%;font-size:13px}}
th,td{{border:1px solid #223;padding:.5rem .6rem;vertical-align:top;text-align:left}}
thead th{{background:#141c30;position:sticky;top:0}}
tbody th{{background:#101830;font-weight:700;width:4rem}}
.st{{font-weight:700;font-size:11px;letter-spacing:.04em}}
.g{{color:#8fa0c0;font-size:11px;font-family:ui-monospace,Menlo,monospace;margin-top:.3rem;word-break:break-all}}
.have .st{{color:#5fd38d}} .req .st{{color:#e8c060}} .none .st{{color:#6a7690}}
.prob{{background:#2a1420}} .prob .st{{color:#ff7b8a}}
.eng{{color:#a8b6d0;font-size:11px}}
.legend{{margin:1.25rem 0;color:#8fa0c0;font-size:12px}}
.legend b{{color:#d8e0f0}}
</style></head><body>
<h1>Element &times; model availability</h1>
<p class="sub">All {len(matrix['elements'])} canonical species &times; {len(matrix['model_types'])} model types, with the actual grid on disk.
Generated {matrix['generated']} by <code>{matrix['generator']}</code> &mdash; reconciled across
the RYA-462 CSV, the RYA-817 3D CSV, the code, and a Sirius <code>find -L</code> scan (control PASS).</p>
<table><thead><tr><th>species</th>{heads}<th>engines</th></tr></thead>
<tbody>{''.join(rows)}</tbody></table>
{_molecule_table(molecules)}
<p class="legend">
<b>HAVE</b> wired and usable &nbsp;|&nbsp; <b>REQUEST_ONLY</b> published, no public deck &nbsp;|&nbsp;
<b>NONE</b> genuinely absent &nbsp;|&nbsp; <b class="pr" style="color:#ff7b8a">DISK_ONLY</b> on Sirius but nothing consumes it
({matrix['problem_count']} such cells).<br>
<b>Engine-A</b> = EW + grid delta (vendored CSV, in repo). <b>Engine-B</b> = synthesis
(Turbospectrum; TS-native NLTE from a Gerber deck on Sirius). Neither is primary.<br>
<b>1D-LTE</b> is available for every species by synthesis and needs no grid.
<b>full-3D-NLTE</b> is NONE everywhere: no public full-3D NLTE stellar RT code exists (RYA-1008).
</p></body></html>
"""


def write_findings_csv(matrix: dict, engine_rows: list[dict],
                       molecules: list[dict], out: Path) -> Path:
    """One flat CSV of EVERY finding: species x model type, engines, and molecules.

    Long format, one row per (subject, model_type), so it sorts and filters in a
    spreadsheet without unpacking anything.
    """
    eng = {r["element"]: r for r in engine_rows}
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["kind", "subject", "model_type", "state", "can_run",
                    "engine", "grid_or_path", "disk_decks", "blocker", "notes"])
        for c in matrix["cells"]:
            e = eng.get(c["element"], {})
            if c["model_type"] == "1D_LTE":
                engine, can = "Engine-B (synthesis)", "yes"
            elif c["model_type"] == "1D_NLTE":
                engine = f'Engine-A: {e.get("engine_a_mode","?")} | ' \
                         f'Engine-B: {e.get("engine_b_mode","?")}'
                can = "yes" if c["state"] in ("HAVE", "CODE_USES") else "no"
            else:
                engine = e.get("engine_c_mode", "none")
                can = "yes" if c["state"] in ("HAVE", "CODE_USES") else "no"
            blocker = next((f for f in c["facts"] if "BROKEN" in f or "RYA-923" in f), "")
            w.writerow([
                "element", c["element"], c["model_type"], c["state"], can, engine,
                c["code_grid"] or c["csv_claim"] or "",
                "; ".join(Path(p).name for p in c["disk_paths"]),
                blocker, " ".join(c["facts"]),
            ])
        for m in molecules:
            w.writerow(["molecule", m["molecule"], "LTE_linelist",
                        "HAVE" if m["lte_linelist"] else "NONE",
                        "yes" if m["lte_linelist"] else "no",
                        "Engine-B (synthesis)", m["lte_linelist"] or "", "", "", ""])
            w.writerow(["molecule", m["molecule"], "NLTE_deck",
                        "HAVE" if m["nlte_deck"] else "NONE",
                        "yes" if m["nlte_deck"] else "no",
                        "Engine-B (TS-native)", "", m["nlte_deck"] or "", "", ""])
    return out
