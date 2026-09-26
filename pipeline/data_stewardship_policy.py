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
import pathlib

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


# --- where raw external data actually lives, post-purge -------------------------------
#
# RYA-1228 removed data/linelists/vald_*_raw.txt from every commit. The content is not
# lost: 54 objects are preserved and hash-verified under the archive below. But code that
# still resolved those deliveries against the repo root started raising FileNotFoundError
# the moment the purge landed, and eight al_grade tests have been erroring in CI ever
# since -- read as "known red" rather than as the purge's own fallout.
#
# 🔴 THE POLICY IS NOT "THE FILES ARE GONE", IT IS "THEY LIVE OUTSIDE GIT". Anything that
# needs a raw delivery resolves it HERE, so there is one answer to where it is, and a
# caller that cannot find it gets a message naming the archive rather than a bare path.

#: Where a preserved delivery may be found, in order. More than one entry because the
#: repository is worked from more than one machine and they do NOT share a filesystem:
#: the curated archive lives on the Mac, while the synthesis host and the CI runner are
#: the same Linux box and need their own staged copy.
#:
#: ⚠️ A MAC-ONLY PATH IS A CI FAILURE. The first version of this resolver named only the
#: Mac archive, so eight al_grade tests kept erroring on the runner while passing locally
#: -- a fix that is only a fix on the machine you tested it on. The staged Linux copy is
#: hash-verified against the archive manifest.
VALD_SEARCH_ROOTS = (
    pathlib.Path("/mnt/codex-data/linelists/vald_raw"),
    pathlib.Path.home() / "Documents" / "Exoplanet Codex" / "vald_raw_archive" / "current",
    pathlib.Path.home() / "Documents" / "Exoplanet Codex" / "vald_raw_archive" / "superseded",
)

#: Kept for the message text and for callers that want to name the curated archive.
VALD_ARCHIVE = (pathlib.Path.home() / "Documents" / "Exoplanet Codex" / "vald_raw_archive")


class RawDeliveryMissing(FileNotFoundError):
    """A preserved raw delivery is not reachable on this machine."""


def vald_delivery_path(name: str, root=None):
    """Resolve a raw VALD delivery by filename, or None if it is not on this machine.

    Looks in the repository first -- a checkout predating the purge still has them, and a
    caller should get the file it is standing on -- then in the preservation archive.
    Returns None rather than raising so a test can skip with a reason; use
    `require_vald_delivery` where absence should be an error.
    """
    if root is not None:
        in_repo = pathlib.Path(root) / "data" / "linelists" / name
        if in_repo.exists():
            return in_repo
    for base in VALD_SEARCH_ROOTS:
        candidate = base / name
        if candidate.exists():
            return candidate
    return None


def require_vald_delivery(name: str, root=None):
    path = vald_delivery_path(name, root)
    if path is None:
        raise RawDeliveryMissing(
            f"{name}: not in the repository (RYA-1228 removed raw VALD deliveries from "
            f"git) and not in any preserved location. Searched: "
            f"{', '.join(str(r) for r in VALD_SEARCH_ROOTS)}. The content IS preserved "
            f"and hash-verified in the archive at {VALD_ARCHIVE} -- stage a copy on this "
            f"machine, do not re-commit the file.")
    return path
