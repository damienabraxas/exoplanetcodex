"""
tests/test_runtime_fork_blas_rya514.py
======================================
RYA-514 — the centralized runtime init forces fork + single-thread BLAS at import, on
every pipeline entry point (covers macOS spawn AND Linux-3.14 forkserver).
"""
import importlib
import multiprocessing as mp
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from config.constants import ISPEC_DIR as _ISPEC_DIR  # noqa: E402
BLAS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
        'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')

# Entry points that must carry the runtime init (touch iSpec child / multiprocessing).
ENTRY_POINTS = ['pipeline._runtime', 'pipeline.lines_fit', 'pipeline.spectra_normalize',
                'pipeline.abundances_derive', 'pipeline.cno_synthesis',
                'pipeline.pysme_nlte', 'pipeline.nlte_bfactor_synth']


def test_runtime_sets_fork_and_pins_in_process():
    sys.path.insert(0, str(ROOT))
    importlib.import_module('pipeline._runtime')
    assert mp.get_start_method(allow_none=True) == 'fork'
    for v in BLAS:
        assert os.environ.get(v) == '1', f'{v} not pinned to 1'


@pytest.mark.parametrize('mod', ENTRY_POINTS)
def test_each_entry_point_yields_fork_in_fresh_process(mod):
    """A FRESH process importing only this entry point must end up on fork (start method is
    process-global, so each is checked in its own interpreter)."""
    code = (f"import sys; sys.path.insert(0, {str(ROOT)!r}); "
            f"sys.path.insert(0, {str(_ISPEC_DIR)!r}); "   # RYA-1140
            "import warnings; warnings.filterwarnings('ignore'); "
            f"import importlib; importlib.import_module({mod!r}); "
            "import multiprocessing as mp, os; "
            "assert mp.get_start_method(allow_none=True)=='fork', mp.get_start_method(); "
            "assert all(os.environ.get(v)=='1' for v in "
            f"{BLAS!r}); print('OK')")
    r = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, timeout=120)
    assert r.returncode == 0 and 'OK' in r.stdout, f"{mod}: {r.stderr[-400:]}"


def test_assert_helper():
    importlib.import_module('pipeline._runtime')
    from pipeline import _runtime
    info = _runtime.assert_fork_and_pins()
    assert info['start_method'] == 'fork'
