"""RYA-682 — the two-engine driver's inputs: one place, checked before any compute.

WHY THIS EXISTS
===============
`scripts/rya527_two_engine_run.py` is the Beta gate's driver. It consumes:

  * ONE **generated** input — the Engine-B synthesis-v2 per-line table. It lives
    under `data/outputs/`, which is **gitignored**, so it is regenerable by
    design and never committed. Regenerate with `REGEN_CMD` (below).
  * SEVEN **committed** inputs — the dedicated Engine-B synthesis-harness
    measurements for the synthesis-required elements (CNO, Mn, Cu/V, Sr II,
    Zr II ×2, Mg 5528). These are tracked in git; a missing one means a broken
    checkout, not a missing run.

Three failure modes were reachable before this module, all of them quiet in the
way RYA-518 exists to forbid:

1. **Late failure.** The driver checked the generated artifact only *after* the
   whole Engine-A leg (GES linelist load, EW triage, MOOG baseline — minutes of
   compute). The message was correct; you just paid for it first.
2. **Silent partial record set.** The dedicated committed inputs were read
   behind bare `if path.exists()`, so a missing one dropped that element's
   Engine-B value and the run still emitted a record set — smaller, unlabelled.
3. **Silently EMPTY generated artifact.** This is the one that bites hardest,
   because the file is present and parses fine. See below.

THE EMPTY-ARTIFACT TRAP (the RYA-682 root cause)
================================================
`ispec/abundances.py:132` assigns a size-1 array into a scalar recarray slot::

    free_abundances['Abund'][i] = solar_abundances['Abund'][solar_abundances['code'] == int(specie)]

NumPy deprecated that in 1.25 and made it a hard error in 2.3. On numpy 2.2.x it
warns and works; on numpy >= 2.3 it raises ``ValueError: setting an array element
with a sequence``. `_run_synthesis_v2_mode` catches that per element, prints
``WARNING: element 'X' not in chem table``, and every line then takes the
``element not in atom_codes`` branch — which appends a row with
``status='failed'``. The frame is NOT empty, so the RYA-342 empty-set guard does
not fire, and the run writes a full-length per-line CSV in which **no row is
usable**, then exits 0.

The two-engine driver reads that file, filters ``status == 'ok'``, gets nothing,
and reports Engine-A-only values for every element that is not synthesis-required
— a one-engine result wearing the two-engine floor's name.

**docs/SCIENCE_STANDARDS.md §RYA-517 ratifies the reference stack as Python 3.12
+ numpy 2.2.x**, "the newest stack with a proven iSpec C-extension build". A venv
that has drifted past numpy 2.2.x cannot generate this artifact. That is an
environment defect, not a code defect, and `assert_synthesis_stack()` names it
instead of letting it turn into a silent empty table.
"""

from __future__ import annotations

import os
from pathlib import Path

from pipeline import data_namespace as ns

ROOT = Path(__file__).resolve().parent.parent

# The bare product suffix; `data_namespace.output_path` adds the {star}_ prefix
# and the per-star directory, so the RYA-469 layout is never hand-built.
ENGINE_B_PRODUCT = 'per_line_synth_v2.csv'

REGEN_CMD = 'python -m pipeline.abundances_derive {star} ATLAS9.Castelli synthesis-v2 --pin'

# RYA-517 / docs/SCIENCE_STANDARDS.md. iSpec's abundance-structure call breaks on
# numpy >= 2.3 (see module docstring), so this is not a stylistic pin.
REFERENCE_NUMPY = '2.2'


def engine_b_per_line_path(star: str = 'solar', *, create: bool = False) -> Path:
    """Canonical (RYA-469 namespaced) path of the generated Engine-B per-line table."""
    return ns.output_path(star, ENGINE_B_PRODUCT, create=create)


def regen_command(star: str = 'solar') -> str:
    return REGEN_CMD.format(star=star)


def _usable_row_count(path: Path) -> int:
    """Rows the two-engine driver would actually consume (status == 'ok')."""
    import pandas as pd
    try:
        df = pd.read_csv(path)
    except Exception as exc:                       # noqa: BLE001 - reported, not swallowed
        raise SystemExit(
            f"RYA-682: Engine-B per-line table at {path} could not be parsed: "
            f"{type(exc).__name__}: {exc}. Regenerate it with:\n    {regen_command()}")
    if 'status' not in df.columns:
        return 0
    return int((df['status'].astype(str) == 'ok').sum())


def assert_engine_b_artifact(star: str = 'solar') -> dict:
    """The generated Engine-B input exists AND carries usable rows.

    Raises SystemExit naming the artifact, the regeneration command, and — when
    the file is present but unusable — the environment cause, which is the case
    that otherwise reads as a successful run.
    """
    path = engine_b_per_line_path(star)
    if not path.exists():
        raise SystemExit(
            f"RYA-682 MISSING GENERATED INPUT: the Engine-B synthesis-v2 per-line table\n"
            f"    {path}\n"
            f"is absent. It lives under data/outputs/, which is gitignored — it is a "
            f"GENERATED input, never committed, so a clean checkout will not have it.\n"
            f"Regenerate it on Sirius (RYA-567: computation is Sirius-only) with:\n"
            f"    {regen_command(star)}\n"
            f"and run it on the RYA-517 reference stack (python 3.12 + numpy "
            f"{REFERENCE_NUMPY}.x); see docs/CONVENTIONS.md.")

    n_ok = _usable_row_count(path)
    if n_ok == 0:
        raise SystemExit(
            f"RYA-682 EMPTY GENERATED INPUT: the Engine-B synthesis-v2 per-line table\n"
            f"    {path}\n"
            f"exists and parses, but has ZERO rows with status == 'ok' — every line "
            f"failed. Running on it would emit an Engine-A-only record set under the "
            f"two-engine floor's name.\n"
            f"The known cause is a numpy-version break: ispec/abundances.py:132 assigns "
            f"a size-1 array into a scalar slot, which numpy >= 2.3 rejects, so every "
            f"element loses its atom code and every line is written status='failed'.\n"
            f"Check the generating environment is the RYA-517 reference stack "
            f"(python 3.12 + numpy {REFERENCE_NUMPY}.x), then regenerate:\n"
            f"    {regen_command(star)}")
    return dict(path=str(path), usable_rows=n_ok)


def assert_committed_inputs(paths: dict) -> dict:
    """Every declared COMMITTED Engine-B input is present.

    `paths` maps a human label -> Path. These are tracked in git, so a missing one
    is a broken checkout; failing loudly beats emitting a quietly smaller record
    set (RYA-518: never a partial result that looks complete).
    """
    missing = {k: v for k, v in paths.items() if not Path(v).exists()}
    if missing:
        lines = '\n'.join(f"    {k}: {v}" for k, v in sorted(missing.items()))
        raise SystemExit(
            f"RYA-682 MISSING COMMITTED INPUT(S): the dedicated Engine-B synthesis "
            f"measurements below are tracked in git but absent from this working tree:\n"
            f"{lines}\n"
            f"These are committed artifacts, not generated ones — restore them "
            f"(`git checkout -- <path>`) rather than re-running a harness. Proceeding "
            f"would drop those elements' Engine-B values and emit a partial record set.")
    return dict(checked=sorted(paths))


def assert_synthesis_stack() -> dict:
    """The running interpreter can actually build iSpec atom codes.

    Exercises the exact call that breaks off the RYA-517 reference stack, so the
    failure surfaces here — named, with a remedy — instead of as 129 warnings and
    an artifact full of status='failed'.
    """
    import numpy as np
    try:
        import pipeline.abundances_derive as ad
        ispec = ad.ispec
        chem = ispec.read_chemical_elements(ad._SYNTH_CHEM_FILE)
        sab = ispec.read_solar_abundances(
            str(Path(ad.ISPEC_DIR) / 'input' / 'abundances' / 'Asplund.2009' / 'stdatom.dat'))
        ispec.create_free_abundances_structure(['Fe'], chem, sab)
    except SystemExit:
        raise
    except Exception as exc:                       # noqa: BLE001 - converted to a loud stop
        raise SystemExit(
            f"RYA-682 SYNTHESIS STACK UNUSABLE: ispec.create_free_abundances_structure "
            f"failed under numpy {np.__version__} —\n"
            f"    {type(exc).__name__}: {exc}\n"
            f"ispec/abundances.py:132 assigns a size-1 array into a scalar recarray slot. "
            f"numpy deprecated that in 1.25 and made it an ERROR in 2.3, so every element "
            f"loses its atom code and the synthesis-v2 run writes a per-line table in "
            f"which every row is status='failed' — while exiting 0.\n"
            f"docs/SCIENCE_STANDARDS.md (RYA-517) ratifies the reference stack as "
            f"python 3.12 + numpy {REFERENCE_NUMPY}.x, 'the newest stack with a proven "
            f"iSpec C-extension build'. This interpreter is off that stack. Use a "
            f"numpy {REFERENCE_NUMPY}.x environment to generate synthesis products.")
    return dict(numpy=np.__version__, atom_codes='ok')


def env_summary() -> str:
    import numpy as np
    return (f"python-numpy {np.__version__} (RYA-517 reference {REFERENCE_NUMPY}.x), "
            f"ISPEC_DIR={os.environ.get('ISPEC_DIR', '<default>')}")
