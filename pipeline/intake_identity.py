"""
pipeline/intake_identity.py
===========================
RYA-973 — THE STANDING INTAKE IDENTITY PROCEDURE. Run this on every frame we download,
every time. It is not an audit step to remember; it is the definition of "we know what
star this is."

**Two routes, always both, and disagreement is the signal.**

* **ASTROMETRY is the ARBITER** — the pointing (`RA`/`DEC`) against SIMBAD positions
  propagated by proper motion to the exposure's epoch (`pipeline.audit_crires`, RYA-952).
* **The LABEL is CORROBORATION** — `resolve_star` on `OBJECT` and `ESO OBS TARG NAME`
  (`pipeline.star_id`, RYA-964).

Never the other way round, and never one alone. RYA-964 proved why in its own field test:
two tau Ceti files carry `OBJECT='HD18884'`, and HD 18884 is **α Ceti — a different star,
20° away**. That is not a misspelling an alias table can catch; it is the CORRECT name of
the WRONG star, and if α Cet were ever added as a system the label would resolve
*confidently and wrongly*. Only astrometry caught it. Conversely RYA-952 found tau Ceti
hiding under `OBJECT='STD'`, where the label says nothing at all and only astrometry
identifies the star. Each route fails in a way the other catches.

⚠️ **Proper motion is load-bearing, not a nicety.** tau Ceti moves 1.92″/yr, so a J2000
position is 42″ stale by 2022 — larger than any usable match radius. Skipping the
propagation does not blur the answer, it changes it.

🔴 **BOTH IDS ARE NORMALISED THROUGH `resolve_star` BEFORE THEY ARE COMPARED.** Comparing
raw strings is how this check destroys itself: until RYA-973, `audit_crires` returned
`tau_cet` (from `data/reference/crires_target_astrometry.csv`) while `system_catalog.csv`,
`stars.yaml` and the holdings registry all said `tau_ceti`. Compared raw, the two routes
DISAGREED on every tau Ceti frame ever taken — a check that cries wolf every time is worse
than no check, because it trains the reader to ignore the one real alarm. The id was
fixed; the normalisation stays, because it is what makes the comparison meaningful rather
than a string coincidence. `assert_star_id_namespace()` below is the guard that stops the
namespaces drifting apart again.

Verdicts:
  `<system_id>`     both routes agree, or astrometry confirms and the labels are
                    unresolved (a role label like `STD`, or a coordinate name)
  `CONTRADICTION`   a label resolves to a DIFFERENT star than the pointing — the HD18884
                    shape. Quarantine; do not guess.
  `INDETERMINATE`   astrometry could not confirm any catalogued target
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

#: Ids that are legitimately outside `system_catalog.star_params_key`, each for a stated
#: reason. Anything else drifting out of the canonical namespace is a defect, not an
#: entry here — this list is deliberately short and every line carries its why.
NAMESPACE_EXEMPT = {
    'tau_boo': 'has data but no star_params_key yet (RYA-957); resolves UNRESOLVED on '
               'purpose until cited parameters are adopted',
    'vesta': 'solar-system body, not a stellar system with parameters',
}


@dataclass
class IntakeIdentity:
    path: str
    verdict: str
    evidence: str
    astrometry_star: str = ''
    astrometry_star_raw: str = ''
    astrometry_sep_arcsec: float = float('nan')
    astrometry_status: str = ''
    label_object: str = ''
    label_object_resolved: str = ''
    label_targ_name: str = ''
    label_targ_resolved: str = ''
    id_namespace_mismatch: bool = False

    @property
    def confirmed(self) -> bool:
        return self.verdict not in ('CONTRADICTION', 'INDETERMINATE')


def identify_at_intake(fits_path, repo_root=None, radius: float | None = None
                       ) -> IntakeIdentity:
    """Run BOTH identity routes on one frame and adjudicate. See the module docstring."""
    from astropy.io import fits
    from config.constants import codex_root
    from pipeline.audit_crires import (MATCH_RADIUS_ARCSEC, identify, load_astrometry,
                                       read_frame)
    from pipeline.star_id import resolve_star

    root = Path(repo_root) if repo_root else Path(codex_root('repo'))
    rad = MATCH_RADIUS_ARCSEC if radius is None else radius
    fr = identify(read_frame(Path(fits_path)), load_astrometry(root), radius=rad)

    astro_raw = fr.star_id or ''
    astro = resolve_star(astro_raw) if astro_raw else 'UNRESOLVED'
    if astro == 'UNRESOLVED' and astro_raw in NAMESPACE_EXEMPT:
        astro = astro_raw                      # exempt ids stand for themselves

    hdr = fits.getheader(str(fits_path))
    obj = str(hdr.get('OBJECT', '') or '').strip()
    targ = next((str(hdr[k]).strip() for k in hdr.keys()
                 if 'OBS TARG NAME' in str(k)), '')
    lab_obj, lab_targ = resolve_star(obj), resolve_star(targ)
    named = [(n, r) for n, r in (('OBJECT', lab_obj), ('OBS TARG NAME', lab_targ))
             if r != 'UNRESOLVED']

    common = dict(path=str(fits_path), astrometry_star=astro,
                  astrometry_star_raw=astro_raw,
                  astrometry_sep_arcsec=float(fr.id_sep_arcsec),
                  astrometry_status=fr.id_status,
                  label_object=obj, label_object_resolved=lab_obj,
                  label_targ_name=targ, label_targ_resolved=lab_targ,
                  id_namespace_mismatch=bool(astro_raw and astro_raw != astro
                                             and astro_raw not in NAMESPACE_EXEMPT))

    if fr.id_status != 'confirmed' or astro in ('', 'UNRESOLVED'):
        return IntakeIdentity(verdict='INDETERMINATE', evidence=(
            f"astrometry confirmed no catalogued target within {rad:g} arcsec "
            f"(status={fr.id_status!r}, nearest={astro_raw or 'none'}). The label is "
            f"corroboration only and cannot stand in for a pointing."), **common)

    disagree = [f"{n}={obj if n == 'OBJECT' else targ!r}->{r}"
                for n, r in named if r != astro]
    if disagree:
        return IntakeIdentity(verdict='CONTRADICTION', evidence=(
            f"astrometry says {astro} at {fr.id_sep_arcsec:.1f} arcsec, but "
            f"{'; '.join(disagree)}. A label can be the CORRECT name of a DIFFERENT star "
            f"(RYA-964: HD18884 on tau Ceti files is alpha Ceti, 20 deg away) — "
            f"astrometry is the arbiter. QUARANTINE; do not guess."), **common)

    silent = [n for n, r in (('OBJECT', lab_obj), ('OBS TARG NAME', lab_targ))
              if r == 'UNRESOLVED']
    return IntakeIdentity(verdict=astro, evidence=(
        f"astrometry: {fr.id_sep_arcsec:.1f} arcsec from {astro}; "
        + (f"labels agree ({', '.join(n for n, _ in named)})" if named
           else "no label names this star")
        + (f"; {', '.join(silent)} carries a role or coordinate name, which names no "
           f"star and is not evidence either way" if silent else "")), **common)


def assert_star_id_namespace(repo_root=None) -> dict:
    """Guard: every star id any intake path can emit must be IDENTICAL to its
    `resolve_star` form, or be a declared exemption.

    This is the tripwire for the defect RYA-973 found. Two registries drifted to
    `tau_cet` and `tau_ceti` for one star, and nothing caught it because nothing joined
    them — until the first consumer ran both identity routes and compared. The drift is
    invisible right up to the moment it makes the identity check useless, so it is
    checked directly rather than waited for."""
    import csv
    from config.constants import codex_root
    from pipeline.star_id import resolve_star

    root = Path(repo_root) if repo_root else Path(codex_root('repo'))
    sources = {
        'data/reference/crires_target_astrometry.csv': 'star_id',
        'data/catalog/holdings_manifest_registry.csv': 'system_id',
    }
    bad = []
    seen = {}
    for rel, col in sources.items():
        path = root / rel
        if not path.exists():
            continue
        for row in csv.DictReader(open(path, newline='')):
            sid = (row.get(col) or '').strip()
            if not sid:
                continue
            seen.setdefault(rel, set()).add(sid)
            if sid in NAMESPACE_EXEMPT:
                continue
            resolved = resolve_star(sid)
            if resolved != sid:
                bad.append(f"{rel}:{col}={sid!r} resolves to {resolved!r} — two names "
                           f"for one star split the identity check")
    if bad:
        raise AssertionError(
            "star-id namespace drift (RYA-973):\n  " + "\n  ".join(bad)
            + "\nEvery intake id must equal its resolve_star() form, or be declared in "
              "NAMESPACE_EXEMPT with a reason.")
    return {'sources': {k: sorted(v) for k, v in seen.items()},
            'exempt': dict(NAMESPACE_EXEMPT)}
