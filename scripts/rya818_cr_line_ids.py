#!/usr/bin/env python3
"""RYA-818 — identify Cr NLTE levels against the term-resolved atom.cr374.

    python3 scripts/rya818_cr_line_ids.py --report          # reach, no writes
    python3 scripts/rya818_cr_line_ids.py --emit --apply    # + write the linelist

RUNS ON SIRIUS. The model atom and the GES linelist live on the grid/engine
volumes, which are Sirius-only (never the Mac). Paths resolve through the RYA-810
register — no literals.

WHAT THIS PRODUCES
------------------
1. A reach report, per species and PER FLAG, into data/audit/rya818/. A line whose
   endpoints did not both resolve is flagged `x` and runs in LTE; it is reported as
   LTE, never folded into coverage (RYA-764).
2. Optionally, the NLTE-tagged linelist block, written next to the source on
   Sirius and md5-pinned. NOT committed: it is a large engine input, and grid/
   engine payloads stay on the data volume.

WHY THE DEFAULT IS --report
---------------------------
The identification rests on an approximation (term-shared departure coefficients)
that changes what the synthesis does. Emitting is opt-in so the reach can be
reviewed before anything downstream can consume it.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os as _os_boot
import sys as _sys_boot
from datetime import datetime, timezone
from pathlib import Path

# Standalone-script bootstrap (RYA-313): repo root on sys.path BEFORE importing
# config/pipeline, so this runs from any cwd. Derived from __file__, never cwd.
_sys_boot.path.insert(0, _os_boot.path.dirname(_os_boot.path.dirname(
    _os_boot.path.abspath(__file__))))

from config.constants import codex_path, require_codex_path  # noqa: E402
from pipeline.model_atom import (atom_resolution,  # noqa: E402
                                 ion_stage_histogram, read_gerber_atom)
from pipeline.nlte_line_identification import (  # noqa: E402
    energy_route_agreement, identification_provenance, identify_lines,
    reach_report, read_species_lines, render_identification_fields)

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "data" / "audit" / "rya818"

ATOM_NAME = "atom.cr374"
ELEMENT = "Cr"
IONS = (1, 2)
ROMAN = {1: "I", 2: "II"}

#: The stage histogram atom.cr374 must present. Asserted, not assumed: reading a
#: differently-named stage attribute returns None for every level and silently
#: pools Cr I with Cr II, which is the failure this pins.
EXPECTED_STAGES = {1: 148, 2: 225, 3: 1}


def _atom_path() -> Path:
    return require_codex_path('grids.gerber_ts') / ATOM_NAME


def _linelist_path() -> Path:
    return require_codex_path('engines.ges_nlte_linelist')


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run(emit: bool, apply: bool, energy_tol_cm: float,
        energy_fallback: bool) -> int:
    atom, linelist = _atom_path(), _linelist_path()
    levels = read_gerber_atom(atom)
    res = atom_resolution(levels)
    stages = ion_stage_histogram(levels)

    print(f"RYA-818 — Cr NLTE line identification")
    print(f"  atom      {atom}")
    print(f"  linelist  {linelist}")
    print(f"  {res.describe()}")
    print(f"  ion stages {stages}")

    if stages != EXPECTED_STAGES:
        print(f"\n  REFUSED: expected stage histogram {EXPECTED_STAGES}, got {stages}.")
        print("  A changed stage split means the atom is not the one this was measured")
        print("  against, or the stage column was misread — either way the per-stage")
        print("  join would silently mix ionisation stages.")
        return 1
    if not res.is_term_resolved:
        print(f"\n  REFUSED: {ATOM_NAME} classifies as {res.verdict}, not term-resolved.")
        print("  This path exists FOR term-resolved atoms. A fine-structure atom should")
        print("  use the upstream converter's label+J route instead.")
        return 1

    summary, frames, controls = {}, {}, {}
    for ion in IONS:
        lines = read_species_lines(linelist, ELEMENT, ion)
        if lines.empty:
            print(f"\n  {ELEMENT} {ROMAN[ion]}: no lines in the linelist")
            continue
        ident = identify_lines(lines, levels, ELEMENT, ion,
                               energy_tol_cm=energy_tol_cm,
                               energy_fallback=energy_fallback)
        rep = reach_report(ident)
        summary[f"{ELEMENT} {ROMAN[ion]}"] = rep
        frames[ion] = ident

        # POSITIVE CONTROL: the energy fallback only fires where the label route
        # failed, so nothing it produces is checkable in production. Run it where
        # the label DID answer and compare. This is what demoted the energy route
        # from "second pass" to "off by default" -- it disagrees with the known
        # answer far too often to run silently.
        ctl = energy_route_agreement(lines, levels, ELEMENT, ion,
                                     energy_tol_cm=energy_tol_cm)
        controls[f"{ELEMENT} {ROMAN[ion]}"] = ctl

        print(f"\n  {ELEMENT} {ROMAN[ion]}: {rep['n_lines']} lines")
        print(f"    NLTE (both endpoints identified) : {rep['n_nlte']} "
              f"({rep['reach_pct']}%)")
        print(f"    by label alone                   : {rep['n_label_both']}")
        print(f"    energy-assisted                  : {rep['n_energy_assisted']}")
        print(f"    RUNS IN LTE despite the NLTE tag : "
              f"{rep['n_lte_despite_nlte_block']}")
        print(f"    lower flags {rep['lower_by_flag']}")
        print(f"    upper flags {rep['upper_by_flag']}")
        print(f"    CONTROL — energy route vs label route where both answer:")
        print(f"      agree {ctl['agree']}  disagree {ctl['disagree']}  "
              f"silent {ctl['energy_route_silent']}  -> {ctl['agreement_pct']}%")
        for ex in ctl["disagreement_examples"][:3]:
            print(f"        {ex['wave_A']:.3f} {ex['end']:<3s} {ex['term']:<8s} "
                  f"label->{ex['label_level']} ({ex['label_E']:.1f} cm-1)  "
                  f"energy->{ex['energy_level']} ({ex['energy_E']:.1f})  "
                  f"computed {ex['computed_E']:.1f}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    prov = identification_provenance(ELEMENT, ATOM_NAME, res, energy_tol_cm,
                                     energy_fallback=energy_fallback)
    prov.update({
        "ticket": "RYA-818",
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "atom_md5": _md5(atom),
        "linelist_md5": _md5(linelist),
        "ion_stage_histogram": stages,
        "reach": summary,
        "energy_route_control": controls,
    })
    (OUT_DIR / "cr_line_identification_rya818.json").write_text(
        json.dumps(prov, indent=2) + "\n", encoding="utf-8")
    print(f"\n  wrote {OUT_DIR / 'cr_line_identification_rya818.json'}")

    for ion, ident in frames.items():
        cols = ["wave_A", "ep_eV", "loggf", "term_low", "term_up", "level_low",
                "level_up", "label_low", "label_up", "flag_low", "flag_up",
                "reason_low", "reason_up", "nlte"]
        p = OUT_DIR / f"cr{ROMAN[ion].lower()}_identifications_rya818.csv"
        ident[cols].to_csv(p, index=False)
        print(f"  wrote {p}")

    if emit:
        if not apply:
            print("\n  [dry-run] --emit needs --apply to write the linelist block.")
            return 0
        out = linelist.parent / f"nlte_ges_cr_identified_rya818.list"
        with open(out, "w", encoding="utf-8") as fh:
            for ion, ident in frames.items():
                fh.write(f"'{ELEMENT} {ROMAN[ion]:<4s} NLTE'\n")
                for rec in ident.to_dict("records"):
                    fh.write(rec["raw"].rstrip() +
                             render_identification_fields(rec) + "\n")
        print(f"\n  wrote {out}  md5 {_md5(out)}")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", action="store_true", help="reach report only (default)")
    ap.add_argument("--emit", action="store_true", help="also build the linelist block")
    ap.add_argument("--apply", action="store_true", help="actually write with --emit")
    ap.add_argument("--energy-tol-cm", type=float, default=50.0,
                    help="energy fallback tolerance, cm-1 (default 50)")
    ap.add_argument("--energy-fallback", action="store_true",
                    help="ALSO identify by energy where the label fails. OFF by "
                         "default: controlled against the label route it names the "
                         "wrong level 18%% (Cr I) / 35%% (Cr II) of the time.")
    a = ap.parse_args()
    raise SystemExit(run(a.emit, a.apply, a.energy_tol_cm, a.energy_fallback))
