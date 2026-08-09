#!/usr/bin/env python3
"""
pipeline/data_namespace.py — per-star output namespacing + the frozen, versioned
gold-standard solar reference (RYA-469).

THE LAW (see docs/design/adr_data_namespacing_and_gold_reference.md):

  1. NO BARE PER-STAR FILENAMES. Every per-star product carries the star in its PATH:
        data/outputs/{star}/{star}_abundances.csv
        data/outputs/{star}/{star}_per_line.csv
        data/outputs/{star}/{star}_ew_integrity.csv
        data/outputs/{star}/{star}_verdict.json
        data/outputs/{star}/diagnostics/...
     Because the star is in the path, two stars CANNOT write the same file — Sirius
     float-drift runs and what-if branches each get their own namespaced artifact and
     diff *named* files instead of fighting over one generic name.
     `data/outputs/` is gitignored and regenerable.

  2. THE SUN IS THE GOLD-STANDARD DIFFERENTIAL DENOMINATOR — it is FROZEN and VERSIONED:
        data/reference/solar/solar_abundances_v{N}.csv   (write-once, immutable, committed)
        data/reference/solar/CURRENT                      (names the active version)
        data/reference/solar/hash_manifest.json           (the immutability guard)
     Each version embeds a provenance header (commit, date, frozen verdict, changelog).
     Re-baselining the Sun BUMPS the version (v{N+1}); it NEVER overwrites an existing
     vN. Promotion is deliberate (scripts/promote_solar_reference.py + reviewed PR).
     A working solar run writes to the namespaced working path (outputs/solar/...), NOT
     the reference. Targets pin a version and record it in their own provenance, so
     re-baselining the Sun never silently changes a target's already-derived numbers.

This module is the single accessor. Don't open these paths by hand elsewhere.
"""
from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

import config.constants as const

ROOT = Path(str(const.ROOT))
OUTPUTS_ROOT = Path(str(const.PATHS['outputs_root']))
SOLAR_REFERENCE_DIR = Path(str(const.PATHS['solar_reference_dir']))
CURRENT_POINTER = SOLAR_REFERENCE_DIR / 'CURRENT'
HASH_MANIFEST = SOLAR_REFERENCE_DIR / 'hash_manifest.json'

# Provenance is embedded as leading '#'-comment lines so each version is a single,
# self-describing, hashable artifact; readers go through read_solar_reference (comment-aware).
_COMMENT = '#'
_VERSION_RE = re.compile(r'^v(\d+)$')
_REF_STEM = 'solar_abundances_'


class ImmutableReferenceError(RuntimeError):
    """A frozen gold reference version was modified (or its hash is missing)."""


class ReferenceVersionExists(RuntimeError):
    """Refused to overwrite an existing immutable reference version."""


# ── star slug ────────────────────────────────────────────────────────────────
def star_slug(star: str) -> str:
    """Canonical filesystem-safe star token. 'Sun'/'sol' -> 'solar'; lowercased, spaces->_."""
    s = str(star).strip().lower().replace(' ', '_')
    if s in ('sun', 'sol'):
        return 'solar'
    if not s:
        raise ValueError("star id is empty — refusing to build a namespaced path")
    return s


# ── per-star working outputs (Deliverable B) ─────────────────────────────────
def outputs_dir(star: str, *, create: bool = True) -> Path:
    d = OUTPUTS_ROOT / star_slug(star)
    if create:
        d.mkdir(parents=True, exist_ok=True)
    return d


def output_path(star: str, name: str, *, create: bool = True) -> Path:
    """Namespaced per-star product path: data/outputs/{star}/{star}_{name}.

    `name` is the bare product suffix, e.g. 'abundances.csv', 'per_line.csv',
    'ew_integrity.csv', 'verdict.json'. The {star}_ prefix is added here so callers
    cannot accidentally emit a bare generic filename.
    """
    slug = star_slug(star)
    fn = name if name.startswith(f'{slug}_') else f'{slug}_{name}'
    return outputs_dir(star, create=create) / fn


def diagnostics_dir(star: str, *, create: bool = True) -> Path:
    d = outputs_dir(star, create=create) / 'diagnostics'
    if create:
        d.mkdir(parents=True, exist_ok=True)
    return d


# ── gold reference: versions (Deliverable C) ─────────────────────────────────
def reference_path(version: str) -> Path:
    """Path to a specific frozen version, e.g. version='v1'."""
    if not _VERSION_RE.match(str(version)):
        raise ValueError(f"bad version token {version!r} (expected 'v<N>', e.g. 'v1')")
    return SOLAR_REFERENCE_DIR / f'{_REF_STEM}{version}.csv'


def list_versions() -> list[str]:
    """Existing frozen versions, ascending by N (['v1', 'v2', ...])."""
    if not SOLAR_REFERENCE_DIR.exists():
        return []
    vs = []
    for p in SOLAR_REFERENCE_DIR.glob(f'{_REF_STEM}v*.csv'):
        m = _VERSION_RE.match(p.stem[len(_REF_STEM):])
        if m:
            vs.append(int(m.group(1)))
    return [f'v{n}' for n in sorted(vs)]


def next_version() -> str:
    vs = list_versions()
    return 'v1' if not vs else f'v{int(vs[-1][1:]) + 1}'


def current_version() -> str:
    """The active gold version named by the CURRENT pointer (fails loud if unset)."""
    if not CURRENT_POINTER.exists():
        raise ImmutableReferenceError(
            f"no CURRENT pointer at {CURRENT_POINTER} — the gold solar reference is not "
            f"initialised (run scripts/promote_solar_reference.py)")
    v = CURRENT_POINTER.read_text().strip()
    if not _VERSION_RE.match(v):
        raise ImmutableReferenceError(f"CURRENT names {v!r}, not a valid 'v<N>' token")
    return v


def resolve_version(version: str = 'CURRENT') -> str:
    return current_version() if version == 'CURRENT' else version


# ── reading the differential denominator (Deliverable D) ─────────────────────
#: RYA-695 — the prefix that reads a CANDIDATE reference instead of a frozen version.
#:
#: WHY A CANDIDATE NEEDS A DOOR OF ITS OWN. RYA-681 made the scale-provenance guard a
#: hard refusal, and gold v3 trips it: it carries A(Fe) 7.466 (the POST 1D->3D value)
#: under method_scale '1D-NLTE (Fe I)', so `resolve_gold_scale` sees a row that
#: contradicts itself and declines to load. That is the guard working — the fix is a
#: re-freeze, never a code exemption — but it also means NO re-emit can run against
#: CURRENT until v4 is frozen, and the re-emit is the evidence the freeze is decided
#: on. RYA-653 already produced the corrected table
#: (`solar_abundances_corrected_candidate_rya653.csv`: same values, `scale_state`
#: 3D-NLTE, `corrections_applied` populated); it simply had no way to be read.
#:
#: Deliberately NOT a version. `reference_path()` still accepts only 'v<N>', the
#: CURRENT pointer is untouched, and the returned provenance token keeps the
#: `candidate:` prefix — so anything that stamps a version records, permanently and
#: visibly, that it was differenced against a candidate and not against frozen gold.
#: Gold stays write-once (RYA-469).
CANDIDATE_PREFIX = 'candidate:'


def candidate_path(name: str) -> Path:
    """Path to a CANDIDATE (unfrozen) solar reference table under the reference dir.

    Refuses anything that is not plainly a candidate: a bare filename, living in the
    reference directory, whose stem contains 'candidate'. Without that last check this
    function is a second way to open a FROZEN version — one that bypasses the version
    resolver and the CURRENT pointer, which is exactly the ambiguity the RYA-469
    write-once rule exists to prevent.
    """
    p = Path(name)
    if p.name != str(name) or p.is_absolute():
        raise ValueError(
            f"candidate reference must be a bare filename inside {SOLAR_REFERENCE_DIR}, "
            f"got {name!r}")
    if 'candidate' not in p.stem:
        raise ValueError(
            f"{name!r} is not a candidate table (its name has no 'candidate'). Frozen "
            f"versions are read through read_solar_reference('v<N>'), never this path — "
            f"one file must not have two doors (RYA-469 write-once).")
    return SOLAR_REFERENCE_DIR / p.name


def read_solar_reference(version: str = 'CURRENT') -> tuple[pd.DataFrame, str]:
    """Read a frozen gold solar reference. Returns (DataFrame, resolved_version).

    Default reads CURRENT. The resolved version string is returned so callers can
    STAMP it into a target's provenance — re-baselining the Sun later never silently
    changes a target's already-derived numbers.

    RYA-695: `version='candidate:<file>.csv'` reads an UNFROZEN candidate table from
    the reference directory and returns the same `candidate:` token as its version, so
    the distinction survives into every artifact downstream. See CANDIDATE_PREFIX.
    """
    if str(version).startswith(CANDIDATE_PREFIX):
        name = str(version)[len(CANDIDATE_PREFIX):]
        path = candidate_path(name)
        if not path.exists():
            raise ImmutableReferenceError(
                f"candidate solar reference missing at {path}")
        return pd.read_csv(path, comment=_COMMENT), f'{CANDIDATE_PREFIX}{path.name}'
    v = resolve_version(version)
    path = reference_path(v)
    if not path.exists():
        raise ImmutableReferenceError(f"gold solar reference {v} missing at {path}")
    df = pd.read_csv(path, comment=_COMMENT)
    return df, v


def differential_denominator(version: str = 'CURRENT') -> tuple[pd.DataFrame, str]:
    """Semantic alias of read_solar_reference for target [X/H]-vs-our-Sun code."""
    return read_solar_reference(version)


def stamp_solar_ref_version(df: pd.DataFrame, version: str = 'CURRENT') -> pd.DataFrame:
    """Record which gold solar version a target was differenced against, in its output."""
    out = df.copy()
    out['solar_ref_version'] = resolve_version(version)
    return out


def read_provenance(version: str = 'CURRENT') -> dict:
    """Parse the embedded '# key: value' provenance header of a version."""
    v = resolve_version(version)
    path = reference_path(v)
    prov: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if not line.startswith(_COMMENT):
            break
        body = line[len(_COMMENT):].strip()
        if ':' in body:
            k, _, val = body.partition(':')
            prov[k.strip()] = val.strip()
    return prov


# ── immutability guard (Deliverable E) ───────────────────────────────────────
def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(Path(path).read_bytes())
    return h.hexdigest()


def load_manifest() -> dict:
    if not HASH_MANIFEST.exists():
        return {}
    return json.loads(HASH_MANIFEST.read_text())


def _write_manifest(manifest: dict) -> None:
    HASH_MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')


def verify_frozen_references() -> list[tuple[str, bool, str]]:
    """Re-hash every committed version against the manifest.

    Returns [(filename, ok, detail)]. ok=False on hash mismatch OR a version present
    on disk but absent from the manifest (an un-recorded write).
    """
    manifest = load_manifest()
    results = []
    for v in list_versions():
        fn = f'{_REF_STEM}{v}.csv'
        path = SOLAR_REFERENCE_DIR / fn
        actual = _sha256(path)
        recorded = manifest.get(fn)
        if recorded is None:
            results.append((fn, False, 'present on disk but ABSENT from hash_manifest.json'))
        elif recorded != actual:
            results.append((fn, False, f'HASH MISMATCH (frozen file edited): '
                                       f'manifest {recorded[:12]} != actual {actual[:12]}'))
        else:
            results.append((fn, True, 'ok'))
    return results


def assert_frozen_references() -> None:
    """Fail LOUD if any frozen gold version was modified. The CI immutability guard."""
    bad = [(fn, detail) for fn, ok, detail in verify_frozen_references() if not ok]
    if bad:
        lines = '\n'.join(f"  {fn}: {detail}" for fn, detail in bad)
        raise ImmutableReferenceError(
            "RYA-469 immutability guard TRIPPED — a frozen gold solar reference changed.\n"
            "Frozen versions are WRITE-ONCE; re-baselining BUMPS the version, never edits "
            "an existing one. Revert the file, or promote a new version via "
            "scripts/promote_solar_reference.py.\n" + lines)


# ── promotion (used by scripts/promote_solar_reference.py) ────────────────────
def _provenance_header(version: str, provenance: dict) -> str:
    lines = [f"{_COMMENT} solar gold-standard abundance reference — RYA-469 (write-once, immutable)",
             f"{_COMMENT} version: {version}"]
    for k, val in provenance.items():
        lines.append(f"{_COMMENT} {k}: {val}")
    lines.append(f"{_COMMENT} immutable: editing this file trips the RYA-469 immutability guard")
    return '\n'.join(lines) + '\n'


def write_reference_version(version: str, df: pd.DataFrame, provenance: dict,
                            *, overwrite: bool = False) -> Path:
    """Write a frozen version (provenance header + CSV) and record its hash. Refuses to
    overwrite an existing version unless overwrite=True (never used in production)."""
    path = reference_path(version)
    if path.exists() and not overwrite:
        raise ReferenceVersionExists(
            f"{path.name} already exists — refusing to overwrite an immutable version. "
            f"Bump to {next_version()} instead.")
    SOLAR_REFERENCE_DIR.mkdir(parents=True, exist_ok=True)
    prov = {'frozen_utc': datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'), **provenance}
    body = df.to_csv(index=False)
    path.write_text(_provenance_header(version, prov) + body)
    manifest = load_manifest()
    manifest[path.name] = _sha256(path)
    _write_manifest(manifest)
    return path


def set_current(version: str) -> None:
    if not reference_path(version).exists():
        raise ImmutableReferenceError(f"cannot point CURRENT at missing version {version}")
    SOLAR_REFERENCE_DIR.mkdir(parents=True, exist_ok=True)
    CURRENT_POINTER.write_text(version + '\n')
