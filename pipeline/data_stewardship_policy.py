"""
What may enter git, and what may not -- RYA-1228 standing policy.

🔴 DELIBERATELY DEPENDENCY-FREE. This module is imported by the pre-commit hook, which
must run on whatever interpreter the developer happens to have. `check_stewardship.py`
cannot serve that role: it imports `pipeline.abundances_derive`, which uses PEP 604
annotations and therefore fails to import at all on Python 3.9 -- so a hook that reached
the policy through it would crash on this very machine and be silently disabled.

One definition, two consumers: the CI-side invariant check and the commit-time gate.
"""

from __future__ import annotations

import fnmatch
import os

#: Raw external inputs. Not Codex measurements -- they live under
#: ~/Documents/Exoplanet Codex (RYA-461) and/or Zenodo (RYA-1124).
RAW_DATA_PATTERNS = (
    'data/linelists/vald_*_raw.txt',
    'data/linelists/vald_*.vald',
    '*.fits', '*.fits.gz', '*.fit', '*.fts',
    'data/model_atmospheres/*.mod',
    'data/model_atmospheres/*.dat',
)

#: A ceiling catches the shape the patterns do not anticipate.
SIZE_CEILING_BYTES = 5 * 1024 * 1024

#: ⚠️ The Codex's own products are legitimately large -- linelist_solar.csv is 25 MB.
#: Exempting tabular and code suffixes is what stops the ceiling from refusing our own
#: measurements, which would train people to pass --no-verify.
CEILING_EXEMPT_SUFFIXES = frozenset({
    '.csv', '.tsv', '.py', '.md', '.json', '.yaml', '.yml', '.txt', '.cfg', '.toml', '.sh',
})


def matched_pattern(path: str) -> str | None:
    """The raw-data pattern this path violates, or None."""
    base = os.path.basename(path)
    for pat in RAW_DATA_PATTERNS:
        if fnmatch.fnmatch(path, pat) or fnmatch.fnmatch(base, pat):
            return pat
    return None


def over_ceiling(path: str, size: int) -> bool:
    """True when a non-product file is too large to be anything but raw input."""
    return (os.path.splitext(path)[1].lower() not in CEILING_EXEMPT_SUFFIXES
            and size > SIZE_CEILING_BYTES)
