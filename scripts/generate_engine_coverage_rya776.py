#!/usr/bin/env python3
"""RYA-776 — generate `data/catalog/engine_coverage.csv`, the engine x wavelength reach
reference, so "do we have Engine A on Fe in the IR?" is a LOOKUP and not a re-derivation.

    python3 scripts/generate_engine_coverage_rya776.py            # write the table
    python3 scripts/generate_engine_coverage_rya776.py --check    # regenerate + diff only

MUST RUN ON SIRIUS. The departure grids and model atoms live at /srv/codex/grids/nlte/
and are never on the Mac. Without them this refuses to write rather than emitting a
table of REACH-UNKNOWN, because a table that says "unknown" everywhere because you ran
it in the wrong place is indistinguishable from a real finding.

WHAT IT ANSWERS, AND WHY THE TRACKER CANNOT
--------------------------------------------
`element_status_tracker.csv` carries WHICH grid and WHAT STATE per engine (Fe:
`Fe_Bergemann_MPIA.csv [A] PRODUCTION`). It carries no WAVELENGTH REACH, so an
(element x engine x wavelength) question is unanswerable from it and has been re-derived
by hand every time it came up. This emits the missing axis as a SIBLING table — the
tracker references it and keeps its own row shape and its own verdict logic untouched.

THE READER IS RYA-763'S, NOT A SECOND ONE
------------------------------------------
Reach comes from `scripts/rya763_level_mapping.py` — `read_levels` (which deck holds
this element's levels) and `resolve_by_label` (identify a level by (J, energy), never by
the index, because the index is native to one deck and silently addresses a different
level in the other). That logic is merged on main and is imported here, not copied.
Forking it would give the project two answers to one question, which is the failure mode
this reference exists to end. This module adds only the SWEEP: every element x ion x
engine x band, plus the extract side, plus the state machine.

THE STATE MACHINE
-----------------
Per (element, ion, engine, grid, band), with `n_served` from the per-line extract and
`n_reach` from the model atom's levels against the in-band catalogued lines:

    n_served > 0                         -> SERVED
    n_served = 0, n_reach > 0            -> REACHABLE-NOT-EXTRACTED
    n_served = 0, n_reach = 0, decidable -> UNCOVERED
    reach not locally decidable          -> REACH-UNKNOWN

"Decidable" is the load-bearing word and it is why there are four states and not three.
Two cases are NOT decidable from any local file, and calling either one UNCOVERED would
manufacture exactly the false "no coverage" this ticket exists to kill:

  * NO LOCAL LEVEL ASSET. Fe's Engine A is the Bergemann/MPIA web service; there is no
    local Fe departure grid to interrogate. The committed extract stops at 6843.7 A, but
    RYA-763 measured the LIVE service still answering 46.7% of probes in 6910-9199 A.
    "The extract is not the model" — writing UNCOVERED at 8000 A would be a lie about a
    service that answers there.
  * NO CATALOGUED LINE IN BAND. Where the GES linelist carries no line, there is nothing
    to resolve, so a zero reach measures OUR LINELIST's span, not the grid's. The
    9199.9 A wall is a linelist limit; `atom.fe607a` reaches 20000 A.
  * THE KEY DOES NOT ADDRESS THE TABLE. The Gerber atoms pack their higher ionisation
    stages into SUPER-LEVELS, whose statistical weight is the sum of the merged levels'
    -- so the reader's J = (G-1)/2 is not a J (Fe II comes out at J = 14.5, 27.5) and
    NOTHING resolves. Measured, not assumed: this generator's first run was about to
    write UNCOVERED against 8870 catalogued Fe II optical lines, a species we measure in
    production as the ionisation arbiter. `classify` separates it from a real absence by
    whether ANY endpoint resolved at all.

DETERMINISM (RYA-768 discipline)
--------------------------------
Rows are sorted on the FULL key (element, ion, engine, grid_id, band_lo_A) — grid_id is
in the key because an element can hold two extracts for one engine (Mg and Si each have
both a Bergemann/MPIA and an Amarsi/PySME CSV), so leaving it out would leave those pairs
ordered by whatever the filesystem returned. That is precisely the RYA-768 defect: a sort
key that is a total order only by accident. Every enumeration below is a sorted() of a
glob or an explicit tuple; no set or dict iteration order reaches the output.
"""
from __future__ import annotations

import argparse
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline import band_policy                                      # noqa: E402
from pipeline.coverage import (                                       # noqa: E402
    ENGINE_COVERAGE, REACH_UNKNOWN, REACHABLE_NOT_EXTRACTED, SERVED, UNCOVERED,
)
from pipeline.nlte_corrections import detect_placeholder_zero_lines   # noqa: E402
# RYA-763's reach-reader, imported and not reimplemented (SSOT — see the docstring).
from scripts.rya763_level_mapping import (                            # noqa: E402
    GERBER_DIR, GRID_DIR, read_gerber_atom, read_labels, resolve_by_label,
)

EXTRACT_DIR = ROOT / "data" / "nlte_grids"
AVAILABILITY = ROOT / "data" / "curation" / "nlte_grid_availability.csv"

COLUMNS = ("element", "ion", "engine", "grid_id", "band", "band_lo_A", "band_hi_A",
           "state", "n_lines_served", "n_lines_reachable", "n_lines_catalogued",
           "level_asset", "grid_asset", "note")

# The (J, energy) tolerance is RYA-763's MEASURED optimum, not a pick: `--scan-tol` on
# Ti I put the resolution rate at 76.6% at 0.001 eV and showed it failing BOTH ways
# (0.0005 -> ABSENT climbs; 0.01 -> AMBIGUOUS explodes and the yield halves).
TOL_EV = 0.001

# Roman is the project's ion convention (`Fe I`). The committed extracts disagree with
# themselves about it -- Al_Amarsi2020_PySME.csv writes `1` where every other extract
# writes `I` -- so it is normalised on read. Left alone it would split Al into two
# species that never join, which is a silent half-empty answer rather than a loud one.
_ROMAN = {"1": "I", "2": "II", "3": "III", "4": "IV", "5": "V",
          "I": "I", "II": "II", "III": "III", "IV": "IV", "V": "V"}
_ARABIC = {v: i + 1 for i, v in enumerate(("I", "II", "III", "IV", "V"))}


def _norm_ion(raw) -> str:
    k = str(raw).strip().upper()
    if k.endswith(".0"):
        k = k[:-2]
    if k not in _ROMAN:
        raise SystemExit(f"unrecognised ion {raw!r} — extend _ROMAN deliberately")
    return _ROMAN[k]


def ges_ion_code(ion: str) -> int:
    """Roman ion -> the GES linelist's arabic stage.

    Written as a lookup because the obvious inline shortcut (`1 if ion == 'I' else 2`)
    silently maps EVERY higher stage onto Fe 2 — it made the Fe III rows duplicates of
    the Fe II ones, with Fe II's catalogued count attached to a species that does not
    have it.
    """
    return _ARABIC[_norm_ion(ion)]


#: Committed paths this table is derived FROM. The provenance line pins these, not HEAD.
_INPUT_PATHS = ("data/nlte_grids", "data/curation/nlte_grid_availability.csv")


def _inputs_commit() -> str:
    """The commit that last changed the INPUTS, not the current checkout.

    Stamping `git rev-parse HEAD` here was self-defeating: committing the artifact moves
    HEAD, so the very next `--check` reported DIFFERS on a byte-identical table and the
    determinism contract could never be satisfied. The tracker generator already models
    the right answer -- it records the provenance of its SOURCE artifact rather than of
    the working tree -- so this pins the extract set instead. It moves when the thing
    being described moves, which is the only time this table should change.
    """
    try:
        return subprocess.run(
            ["git", "log", "-1", "--format=%h", "--", *_INPUT_PATHS],
            cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip() or "none"
    except Exception:                                   # noqa: BLE001
        return "unknown"


# ── the asset side ───────────────────────────────────────────────────────────

def engine_a_extracts() -> list[tuple[str, str, Path]]:
    """(element, grid_id, path) for every committed Engine-A per-line extract.

    Enumerated from the ASSETS, not from `nlte_grid_availability.csv`, because the
    ticket asks for a table generated from what we hold. The availability registry is
    read afterwards only to ANNOTATE the role — where the two disagree, the file on disk
    wins the coverage question and the disagreement shows up in the note.
    """
    out = []
    for p in sorted(EXTRACT_DIR.glob("*.csv")):
        if p.name.endswith(".prov.json") or "_" not in p.name:
            continue
        el = p.name.split("_", 1)[0]
        if not el.isalpha():
            continue
        out.append((el, p.stem, p))
    return out


def engine_b_atoms() -> list[tuple[str, Path, Path | None]]:
    """(element, atom_path, grid_bin_or_None) for the Gerber TS-native deck.

    Engine B reads the model atom and the departure grid DIRECTLY at synthesis time --
    there is no per-line extract in between, which is why its SERVED test is "atom
    resolves the line AND the .bin is provisioned" rather than "a CSV lists it".
    """
    out = []
    for atom in sorted(GERBER_DIR.glob("atom.*")):
        # atom.fe607a -> Fe, atom.na_qmh -> Na. The element is the LEADING alphabetic
        # run of the suffix, not every letter in it: `"".join(c for c in ... isalpha())`
        # would read atom.fe607a as 'Fea'.
        stem = atom.name.split(".", 1)[1]
        el = ""
        for c in stem:
            if not c.isalpha():
                break
            el += c
        if not el:
            continue
        el = el.capitalize()
        hits = sorted(GERBER_DIR.glob(f"NLTEgrid4TS_{el}_*.bin")) or \
            sorted(GERBER_DIR.glob(f"NLTEgrid4TS_{el.upper()}_*.bin"))
        out.append((el, atom, hits[0] if hits else None))
    return out


def _registered_ion(element: str, fname: str) -> int:
    """The ion a per-line extract applies to, for the extracts that do not say.

    Read from `config.constants.NLTE_CORRECTION_ELEMENTS`, which is already the project's
    single source for which ion each element's departure grid corrects (Ba 2, Sr 2, the
    rest 1). Not re-declared here -- a second copy of this mapping is exactly the kind of
    duplicate that drifts. Where an element registers a NEWER grid than the file in hand
    (Ti registers Mallinson-2024, and Ti_Bergemann2011_MPIA.csv is the superseded
    vintage), the ION is still the element's and is what we want: same species, older
    model.
    """
    from config.constants import NLTE_CORRECTION_ELEMENTS
    rec = NLTE_CORRECTION_ELEMENTS.get(element)
    if not rec or rec.get("ion") is None:
        raise SystemExit(
            f"{fname} carries no `ion` column and {element} is not in "
            f"NLTE_CORRECTION_ELEMENTS, so the species this extract corrects cannot be "
            f"resolved. Register it there rather than letting the row vanish — a coverage "
            f"table that silently omits a grid is worse than one that fails loudly.")
    return int(rec["ion"])


def availability_roles() -> dict[tuple[str, str], str]:
    """(element, grid_stem) -> role, from the curation registry. Annotation only."""
    if not AVAILABILITY.exists():
        return {}
    roles = {}
    for r in pd.read_csv(AVAILABILITY, comment="#").to_dict("records"):
        gf = str(r.get("grid_file") or "").strip()
        if not gf.endswith(".csv"):
            continue
        roles[(str(r["element"]).strip(), Path(gf).stem)] = str(r.get("role") or "").strip()
    return roles


# ── the reach side (RYA-763's reader, swept) ─────────────────────────────────

def levels_for(element: str, engine: str) -> tuple[pd.DataFrame | None, str]:
    """Level table + the asset that supplied it, or (None, why-not).

    Engine A's model lives in the Amarsi/PySME deck (`label_{El}.txt`); Engine B's in the
    Gerber TS-native deck (`atom.*`). They are NOT interchangeable: RYA-763 measured the
    GES level index at 0% mismatch against Gerber and 48.7% against Amarsi, so mixing
    decks is the actual hazard. This keeps each engine on its own deck rather than
    calling `read_levels(deck='auto')`, whose fallback would silently hand Engine A a
    Gerber atom and report a reach the MPIA route does not have.
    """
    if engine == "A":
        lab = GRID_DIR / f"label_{element}.txt"
        if lab.exists():
            return read_labels(element), lab.name
        return None, (f"no local Engine-A level table (no {lab.name}) — this element's "
                      f"Engine A is a web-service supplier, so its reach is not locally "
                      f"decidable")
    hits = sorted(GERBER_DIR.glob(f"atom.{element.lower()}*"))
    if hits:
        return read_gerber_atom(hits[0]), hits[0].name
    return None, f"no Gerber model atom for {element}"


#: Sanity bound on a per-stage energy zero point, in eV. NOT a tuned parameter: the
#: first ionisation potentials of every element in these decks fall inside it (K 4.34 at
#: the bottom, O 13.62 at the top), so an offset outside it means the stage's ground
#: state is missing from the atom and the rebase below would be inventing a coordinate.
_IP_PLAUSIBLE_EV = (3.0, 30.0)


def _ion_filter(lab: pd.DataFrame, ion: str) -> tuple[pd.DataFrame, float]:
    """Restrict a level table to one ion, on that ion's OWN energy zero point.

    Returns (levels, offset_applied_eV).

    Both decks carry the ionisation stage, in different columns: the Gerber atom reader
    emits an integer `ion`, `label_{El}.txt` an `Fe 1`-style `species` string. Filtering
    matters -- resolve_by_label keys on (J, energy) alone, so an unfiltered table lets a
    neutral line pick up a singly-ionised level of coincidentally similar energy and
    report a reach that does not exist.

    THE ZERO POINT IS NOT SHARED, AND THAT COST A FALSE ABSENCE. A Gerber model atom
    numbers every stage's energies CUMULATIVELY from the NEUTRAL ground state -- the
    Mn II levels start at 7.434 eV, which is Mn I's ionisation potential exactly -- while
    the GES linelist measures each ion's levels from THAT ION's ground state. Compared
    raw, not one of 3386 Mn II endpoints resolved, and this generator wrote UNCOVERED
    over 1693 catalogued optical lines. Rebased, 171 lines resolve both endpoints. The
    gap was a coordinate mismatch wearing the costume of absent physics.

    The rebase reads the atom's own convention rather than fitting anything: a model
    atom's lowest level in a stage IS that stage's ground state. The one way that can be
    false is a stage whose ground state was omitted, which would silently shift every
    energy and manufacture matches -- so the offset is bounds-checked against the range
    real ionisation potentials occupy, and a stage that fails the check is left on the
    raw scale, where it will read as REACH-UNKNOWN rather than as a confident answer.
    """
    want = _ARABIC[ion]
    if "ion" in lab.columns:
        sub = lab[lab["ion"].astype(int) == want]
    elif "species" in lab.columns:
        tail = lab["species"].astype(str).str.strip().str.split().str[-1]
        sub = lab[tail == str(want)]
    else:
        return lab, 0.0
    if want == 1 or sub.empty:
        return sub, 0.0                 # the neutral stage already starts at zero
    off = float(sub["energy_eV"].min())
    if not (_IP_PLAUSIBLE_EV[0] <= off <= _IP_PLAUSIBLE_EV[1]):
        return sub, 0.0
    return sub.assign(energy_eV=sub["energy_eV"] - off), off


def reach_in_band(ll, lab: pd.DataFrame, element: str, ion: str,
                  lo: float, hi: float) -> tuple[int, int, int]:
    """(n_catalogued, n_both_endpoints, n_either_endpoint) for one species in one band.

    n_both is the reach proper: a departure coefficient needs the lower AND the upper
    level, so half a mapping is not a partial answer, it is no answer. That is the
    both-ends test RYA-763 used and this reproduces.

    n_either exists to tell two very different zeroes apart, and it is load-bearing --
    see `classify`. A level table where SOME endpoints resolve and some do not is
    answering the question and reporting real absences (RYA-763's Fe I IR result was
    exactly this: 4189 upper levels genuinely missing from the 607-level atom). A table
    where NOTHING resolves, over thousands of attempts, is not reporting absence -- it is
    telling you the key does not address it at all.
    """
    w = np.asarray(ll["wave_A"], dtype=float)
    els = np.asarray([str(x).strip() for x in ll["element"]])
    want = f"{element} {ges_ion_code(ion)}".upper()
    m = np.array([e.upper() == want for e in els]) & (w >= lo) & (w < hi)
    n_cat = int(m.sum())
    if not n_cat or lab is None or lab.empty:
        return n_cat, 0, 0
    n_both = n_either = 0
    for i in np.where(m)[0]:
        vlo, _, _ = resolve_by_label(lab, float(ll["lower_state_eV"][i]),
                                     float(ll["lower_j"][i]), TOL_EV)
        vup, _, _ = resolve_by_label(lab, float(ll["upper_state_eV"][i]),
                                     float(ll["upper_j"][i]), TOL_EV)
        n_either += int(vlo == "UNIQUE" or vup == "UNIQUE")
        n_both += int(vlo == "UNIQUE" and vup == "UNIQUE")
    return n_cat, n_both, n_either


# ── the state machine ────────────────────────────────────────────────────────

def classify(n_served: int, n_reach: int, n_cat: int,
             level_asset: str, why_no_levels: str, n_either: int = 0) -> tuple[str, str]:
    """(state, note). The four outcomes, and never a fifth dressed as UNCOVERED.

    The last guard is the one that had to be MEASURED into existence. The Gerber model
    atoms pack their higher ionisation stages into SUPER-LEVELS -- composites whose
    statistical weight is the sum of the merged levels', so the reader's J = (G-1)/2 is
    not a J at all (Fe II comes out with J = 14.5, 27.5). The (J, energy) key therefore
    does not address those rows, and the first run of this generator was about to write
    UNCOVERED against 8870 catalogued Fe II optical lines -- a species we measure in
    production, as the ionisation arbiter. That is precisely the false "no coverage" this
    table exists to end, so a zero that comes with ZERO partial matches is treated as the
    key failing, not as the physics being absent.

    The distinction is measured, not thresholded: `n_either > 0` means the table answered
    for some endpoints and really lacks the others (RYA-763's Fe I IR result, 4189 upper
    levels genuinely missing); `n_either == 0` over thousands of attempts means nothing
    in this table is addressable by the key at all.
    """
    if n_served > 0:
        return SERVED, ""
    if not level_asset:
        return REACH_UNKNOWN, why_no_levels
    if n_cat == 0:
        return REACH_UNKNOWN, ("no catalogued line in band — a zero here measures the "
                               "LINELIST's span, not the grid's, so absence cannot be "
                               "asserted from it")
    if n_reach > 0:
        return REACHABLE_NOT_EXTRACTED, ("levels present in the model; no per-line "
                                         "extract exposes this band")
    if n_either == 0:
        return REACH_UNKNOWN, (
            f"not one endpoint of {n_cat} catalogued line(s) resolves in {level_asset} — "
            f"the (J, energy) key does not address this table (super-levels carry a "
            f"summed statistical weight, so G is not 2J+1), so absence cannot be "
            f"asserted from it")
    return UNCOVERED, (f"{n_cat} catalogued line(s) in band; {n_either} resolve one "
                       f"endpoint, none resolve both in {level_asset}")


# ── the sweep ────────────────────────────────────────────────────────────────

def build_rows(ll) -> list[dict]:
    roles = availability_roles()
    bands = sorted(band_policy.POLICIES, key=lambda p: p.lo_A)
    rows: list[dict] = []

    # ---- Engine A: the per-line departure extracts --------------------------
    for element, grid_id, path in engine_a_extracts():
        ext = pd.read_csv(path)
        if "wave_A" not in ext.columns:
            raise SystemExit(f"{path.name} has no wave_A column — it is not a per-line "
                             f"extract; exclude it deliberately rather than by accident")
        # SIX of the committed extracts carry no `ion` column at all (the Bergemann/MPIA
        # vintage: Ca, Cr, Mg, Mn, Si, Ti). An earlier pass `continue`d past them, which
        # silently dropped Ca and Cr -- both in PRODUCTION -- out of the table entirely,
        # and the tracker then read "(no engine rows)" for them. A silent drop is the one
        # outcome a coverage reference may never produce, so the ion is resolved from the
        # project's own registry instead, and an element missing from BOTH is loud.
        if "ion" in ext.columns:
            ions = [_norm_ion(v) for v in ext["ion"]]
        else:
            ions = [_norm_ion(_registered_ion(element, path.name))] * len(ext)
        ext = ext.assign(_ion=ions,
                         _w=pd.to_numeric(ext["wave_A"], errors="coerce"))
        # RYA-413/417: a line that is identically 0.000 across EVERY node of the extract
        # is a placeholder, not a correction -- MPIA offers it in its dropdown and returns
        # zero for it. Counting it as SERVED would be this table's own failure mode, an
        # extract entry mistaken for a modelled line, so it is subtracted and SAID.
        zeros = set(detect_placeholder_zero_lines(ext)) if "delta_nlte" in ext else set()
        lab_all, level_asset_or_why = levels_for(element, "A")
        has_levels = lab_all is not None
        role = roles.get((element, grid_id), "")
        for ion in sorted(set(ext["_ion"])):
            sub = ext[ext["_ion"] == ion]
            lab, _off = _ion_filter(lab_all, ion) if has_levels else (None, 0.0)
            for b in bands:
                served = sub[(sub["_w"] >= b.lo_A) & (sub["_w"] < b.hi_A)]
                waves = {round(float(x), 3) for x in served["_w"].dropna().unique()}
                n_zero = len(waves & zeros)
                n_served = len(waves) - n_zero
                n_cat, n_reach, n_either = reach_in_band(ll, lab, element, ion,
                                                         b.lo_A, b.hi_A)
                state, note = classify(n_served, n_reach, n_cat,
                                       level_asset_or_why if has_levels else "",
                                       level_asset_or_why, n_either)
                if n_zero:
                    note = (f"{n_zero} placeholder-zero line(s) in band excluded from "
                            f"served (RYA-413)" + (f"; {note}" if note else ""))
                if role and note:
                    note = f"{note}; registry role={role}"
                elif role:
                    note = f"registry role={role}"
                rows.append(dict(
                    element=element, ion=ion, engine="A", grid_id=grid_id,
                    band=b.name, band_lo_A=b.lo_A, band_hi_A=b.hi_A, state=state,
                    n_lines_served=n_served, n_lines_reachable=n_reach,
                    n_lines_catalogued=n_cat,
                    level_asset=level_asset_or_why if has_levels else "",
                    grid_asset=path.name, note=note))

    # ---- Engine B: the Gerber TS-native deck --------------------------------
    for element, atom, grid_bin in engine_b_atoms():
        lab_all = read_gerber_atom(atom)
        if lab_all.empty:
            continue
        for ion in sorted({_norm_ion(v) for v in lab_all["ion"]}):
            lab, off = _ion_filter(lab_all, ion)
            # A stage carrying fewer than two levels is the atom's IONISATION CONTINUUM
            # reservoir, not a modelled ion -- every one of these atoms ends in a single
            # top-stage level (atom.fe607a: 548 Fe I + 58 Fe II + 1 Fe III). A line needs
            # two levels of the same stage, so such a stage can never serve one BY
            # CONSTRUCTION, and emitting a measured verdict against it is a category
            # error: the first run wrote 8 spurious UNCOVERED rows for Fe III, Ba III,
            # Ti III and friends, each of which reads as a real modelling gap.
            if len(lab) < 2:
                continue
            for b in bands:
                n_cat, n_reach, n_either = reach_in_band(ll, lab, element, ion,
                                                         b.lo_A, b.hi_A)
                # Engine B has no extract layer: a line is SERVED when the atom resolves
                # it AND the departure grid is actually provisioned. An atom without its
                # .bin is the cheap-to-unlock class, not a gap -- which is exactly Fe's
                # Engine-B history (grid gettable, not yet pulled).
                n_served = n_reach if grid_bin is not None else 0
                state, note = classify(n_served, n_reach, n_cat, atom.name, "", n_either)
                if grid_bin is None and state == REACHABLE_NOT_EXTRACTED:
                    note = ("model atom carries the levels but the departure grid .bin "
                            "is NOT provisioned — a pull away, not a data gap")
                if off:
                    # Said out loud, per row: an energy shift applied silently would be
                    # unauditable, and this one decides SERVED vs UNCOVERED for Mn II.
                    note = (f"levels rebased to this ion's ground state (-{off:.3f} eV; "
                            f"the atom counts from the NEUTRAL ground state)"
                            + (f"; {note}" if note else ""))
                rows.append(dict(
                    element=element, ion=ion, engine="B", grid_id=atom.name,
                    band=b.name, band_lo_A=b.lo_A, band_hi_A=b.hi_A, state=state,
                    n_lines_served=n_served, n_lines_reachable=n_reach,
                    n_lines_catalogued=n_cat, level_asset=atom.name,
                    grid_asset=grid_bin.name if grid_bin else "", note=note))

    # RYA-768: the FULL key, grid_id included. Two extracts for one (element, engine)
    # is real (Mg, Na, Si, Ti), so a key without grid_id is a partial order and the
    # artifact never byte-diffs clean.
    rows.sort(key=lambda r: (r["element"], r["ion"], r["engine"], r["grid_id"],
                             r["band_lo_A"]))
    return rows


def render(rows: list[dict]) -> str:
    head = (
        "# ENGINE x WAVELENGTH COVERAGE — **GENERATED by "
        "scripts/generate_engine_coverage_rya776.py — do not hand-edit**\n"
        "# =====================================================================\n"
        "# RYA-776. The sibling of instrument_coverage: that says which INSTRUMENT sees a\n"
        "# wavelength, this says which ENGINE/GRID can model it. Read it through\n"
        "# pipeline.coverage.engine_reach(), never by eye.\n"
        "#\n"
        "# state:\n"
        "#   SERVED                   the per-line extract resolves it. Usable now.\n"
        "#   REACHABLE-NOT-EXTRACTED  the model atom's levels cover it, no extract does.\n"
        "#                            A derivation away, NOT a data gap.\n"
        "#   UNCOVERED                catalogued lines exist here and the model carries\n"
        "#                            none of them. The only state that means 'we cannot\n"
        "#                            model this'.\n"
        "#   REACH-UNKNOWN            reach is not locally decidable — a service-only\n"
        "#                            supplier (Fe/MPIA) or no catalogued line in band.\n"
        "#                            NOT a claim of absence.\n"
        "#\n"
        "# n_lines_catalogued is the GES in-band denominator; a 0 there is OUR linelist's\n"
        "# limit, not the grid's (the 9199.9 A wall is GES, atom.fe607a reaches 20000 A).\n"
        f"# bands: pipeline.band_policy · reach reader: scripts/rya763_level_mapping.py\n"
        f"# extracts @ {_inputs_commit()} · decks on {platform.node()}: "
        f"{GRID_DIR.name} + {GERBER_DIR.name}\n")
    body = pd.DataFrame(rows, columns=list(COLUMNS)).to_csv(index=False,
                                                            lineterminator="\n")
    return head + body


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true",
                    help="regenerate and diff against the committed file; write nothing")
    ap.add_argument("--out", default=str(ENGINE_COVERAGE))
    a = ap.parse_args(argv)

    if not GRID_DIR.exists() and not GERBER_DIR.exists():
        raise SystemExit(
            f"neither {GRID_DIR} nor {GERBER_DIR} exists. This must run on SIRIUS — the "
            f"grids and model atoms are never on the Mac. REFUSING to emit a table of "
            f"REACH-UNKNOWN that would be indistinguishable from a real finding.")

    from pipeline.abundances_derive import _load_synth_resources
    ll, _, _ = _load_synth_resources()

    rows = build_rows(ll)
    text = render(rows)
    out = Path(a.out)

    n = len(rows)
    by_state: dict[str, int] = {}
    for r in rows:
        by_state[r["state"]] = by_state.get(r["state"], 0) + 1
    print(f"\n  {n} rows over "
          f"{len({(r['element'], r['ion']) for r in rows})} species x "
          f"{len({r['grid_id'] for r in rows})} grids x "
          f"{len(band_policy.POLICIES)} bands")
    for k in sorted(by_state, key=lambda s: -by_state[s]):
        print(f"    {k:<26} {by_state[k]:5d}")

    if a.check:
        if not out.exists():
            print(f"\n  {out} does not exist yet — nothing to diff")
            return 1
        same = out.read_text(encoding="utf-8") == text
        print(f"\n  {'BYTE-IDENTICAL' if same else 'DIFFERS'} vs {out}")
        return 0 if same else 1

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text, encoding="utf-8")
    print(f"\n  wrote {out.relative_to(ROOT) if out.is_relative_to(ROOT) else out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
