"""RYA-1395 v74: v73 with method 1 (node counting) also for odd 5d. Declared on RYA-1395 before any run.
v73 FAILED_RETAINED at 5d: method-3 stagnation (self-consistency 3.6e-7/5.2e-7 vs ACCY 6e-8) on a physically correct
5d (2 nodes, eps -1% vs NIST). Everything else is v73."""
from pathlib import Path
import shutil
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import rya1395_mcdhf_s0_v73 as V73  # noqa: E402

V73.METHOD1 = V73.METHOD1 | {'5d'}
V73.SCHEMA = 'rya1395.mcdhf-s0-v74'
_run_parity = V73.V71.run_parity


def _run(a):
    shutil.copyfile(__file__, a.out / Path(__file__).name)
    return _run_parity(a)


V73.V71.run_parity = _run

if __name__ == '__main__':
    V73.main()
