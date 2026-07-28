"""
pipeline/molecular_lists.py
===========================
RYA-360 — the securing layer for the C/N/O molecular line lists.

RYA-236 acquired/verified the Turbospectrum molecular `.bsyn` lists (CH/¹³CH, CN
isotopologues, C2, OH, NH) plus the converted CO IR list (CO_IR_Li2015.dat), but
they live ONLY in the iSpec install tree
(`ispec/input/linelists/turbospectrum/molecules/`) — the project repo tracked
**zero** molecular artifacts. An iSpec reinstall/rebuild would silently reset them to
the stock bundle and wipe the RYA-236 CO addition, with nothing to catch it (the same
failure class as the gf / blend_flag / STAR_PARAMS arc).

This module is the single source of the vendored-copy paths + the `.bsyn` line
counter, shared by the vendoring script (`scripts/vendor_molecular_lists_rya360.py`)
and the stewardship guard (`scripts/check_stewardship.py` `[molecular]` invariant), so
neither re-implements the file-format reader or hard-codes the layout.

Turbospectrum `.bsyn` format: two header lines (line 1 `'  <species_code>' 1 <nlines>`,
line 2 `'<Name Source>'`) then one data row per line, first whitespace-delimited token
= wavelength in Å. `CO_IR_Li2015.dat` uses the same layout.
"""
from __future__ import annotations

import json
from pathlib import Path

from config.constants import ISPEC_DIR

_REPO = Path(__file__).resolve().parents[1]

# The vendored, git-tracked secure record of truth (RYA-360 Step 1).
VENDORED_DIR = _REPO / 'data' / 'linelists' / 'molecular' / 'turbospectrum'
MANIFEST_PATH = VENDORED_DIR / 'MOLECULAR_MANIFEST.json'

# The iSpec install tree the pipeline (Turbospectrum, use_molecules=True) reads from.
ISPEC_MOLECULES_DIR = Path(str(ISPEC_DIR)) / 'input' / 'linelists' / 'turbospectrum' / 'molecules'


def count_bsyn_lines(path) -> int:
    """Count data rows in a Turbospectrum `.bsyn` / `.dat` molecular list — every line
    except the two header lines (those begin with a quote after left-strip). Robust to
    the exact header count; returns 0 for an empty/truncated file."""
    n = 0
    with open(path, 'r', errors='ignore') as fh:
        for ln in fh:
            s = ln.lstrip()
            if not s or s[0] == "'":
                continue
            # a data row starts with the wavelength (a float)
            tok = s.split(None, 1)[0]
            try:
                float(tok)
            except ValueError:
                continue
            n += 1
    return n


def bsyn_wavelengths(path):
    """Yield the wavelength (Å, first column) of every data row — for span/window
    measurement. Uses the same header rule as ``count_bsyn_lines``."""
    with open(path, 'r', errors='ignore') as fh:
        for ln in fh:
            s = ln.lstrip()
            if not s or s[0] == "'":
                continue
            tok = s.split(None, 1)[0]
            try:
                yield float(tok)
            except ValueError:
                continue


def source_label(path) -> str:
    """The `.bsyn` line-2 provenance label (e.g. 'CH PGopher', 'ExoMol Li2015'),
    read from the file itself (machine provenance)."""
    with open(path, 'r', errors='ignore') as fh:
        fh.readline()                       # species-code header
        second = fh.readline().strip()
    return second.strip("'").strip()


def load_manifest() -> dict:
    """Load the vendored provenance manifest (the guard registry). Raises if absent —
    the manifest IS the securing record; a missing one is a loud failure, not a skip."""
    if not MANIFEST_PATH.exists():
        raise FileNotFoundError(
            f"molecular manifest absent at {MANIFEST_PATH} — run "
            f"scripts/vendor_molecular_lists_rya360.py to vendor + record the C/N/O "
            f"molecular lists (RYA-360).")
    with open(MANIFEST_PATH) as fh:
        return json.load(fh)
