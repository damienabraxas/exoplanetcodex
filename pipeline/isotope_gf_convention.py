"""RYA-684 — the isotope-fraction gf convention, and the guard that enforces it.

THE CONVENTION
==============
Turbospectrum multiplies an isotope-coded species' population by
``isotopfrac(Z, A)`` before forming the line opacity (``bsyn.f:1350``).
``makeabund.f`` sets ``isotopfrac(Z, 0) == 1.0`` and says why, verbatim:

    isotopfrac(x,0)==1.0; this is used in the case of no isotopes wanted
    in calculation. It will ensure that the total abundance of the species
    will be used.(if for example isotopic factors were included in gf-values)

So there are exactly TWO self-consistent ways to hand an isotope-split feature
to the engine, and a line list must commit to one of them:

  (A) **ISOTOPE-CODED + FRACTION-FREE gf** — the block is written ``Z.AAA`` and
      each isotope's components carry the FULL oscillator strength.  The engine
      applies the fraction.  This is what the GES v6 HFS/ISO list and the
      Gerber NLTE deck do.

  (B) **UNCODED + FOLDED gf** — the block is written ``Z.000`` and the gf values
      already carry the isotopic split.  ``isotopfrac(Z, 0) == 1`` so the engine
      applies nothing.  This is what ``scripts/rya581_ba2_deblend_sirius.py``
      does when it writes its own Ba II 5853 HFS block.

**ISOTOPE-CODED + FOLDED gf is the error.**  The fraction lands twice, the
feature is synthesised ``sum_i f_i^2`` instead of ``sum_i f_i == 1`` too weak,
and a fit that inverts an observed EW against it has to raise the abundance by
``-log10(sum_i f_i^2)`` to compensate — +0.3002 dex for Eu II, +0.2694 for
Ba II, +0.2415 for Cu I.

The shipped TSFitPy ``linelist_vald`` "for-grid" lists are form (B) gf values
written with form (A) headers.  That is the defect RYA-684 measured, and it is
a property of those files, not of any harness.

WHAT THIS MODULE ENFORCES
=========================
A harness may legitimately carry an exposed block in its BLEND model — the bias
is second-order there and RYA-684 measured it as <0.01 % of window absorption in
every window that feeds a live value.  What is never acceptable is fitting a
TARGET species against a double-folded block, because then the whole offset
lands directly on the reported abundance.

``assert_target_convention()`` is the preflight: call it with the line list
about to go to ``bsyn`` and the species being fitted, and it loud-fails if that
species appears isotope-coded on a surface known to ship folded gf.

The list of folded species is NOT hardcoded here — it is read from the committed
RYA-684 audit record, which ``scripts/rya684_isotope_gf_audit.py`` regenerates
by measuring the shipped files.  Re-vendor a line list, re-run the audit, and
this guard tracks it.
"""

from __future__ import annotations

import json
import math
import os
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
AUDIT_RECORD = ROOT / 'data' / 'audit' / 'rya684_isotope_gf_audit.json'

# The surface whose files ship folded gf.  Named so a future surface added to
# the audit cannot be silently swept into the same bucket.
VALD_FOR_GRID_SURFACE = 'vald(bsyn harnesses)'

_HDR = re.compile(r"^'\s*(\d+)\.(\d+)\s*'\s+(\d+)\s+(\d+)")


def load_audit_record(path: os.PathLike | str = AUDIT_RECORD) -> dict:
    """Read the committed RYA-684 audit record (loud-fail, never a default)."""
    p = Path(path)
    if not p.exists():
        raise SystemExit(
            f"RYA-684: isotope-gf audit record missing at {p}. Regenerate it with "
            f"scripts/rya684_isotope_gf_audit.py on Sirius — the convention guard "
            f"must not fall back to a hardcoded species list.")
    with open(p) as fh:
        rec = json.load(fh)
    if not rec.get('results'):
        raise SystemExit(f"RYA-684: audit record at {p} carries no results. STOP.")
    return rec


def folded_species(surface: str = VALD_FOR_GRID_SURFACE,
                   record: dict | None = None) -> dict:
    """{(Z, ion): result} for every species shipping folded gf on `surface`."""
    rec = record if record is not None else load_audit_record()
    return {(r['Z'], r['ion']): r for r in rec['results']
            if r['surface'] == surface and r['verdict'] == 'FOLDED'}


def isotope_coded_blocks(linelist_path: os.PathLike | str) -> list:
    """Every isotope-coded ``'Z.AAA' ion nlines`` block header in a TS line list."""
    out = []
    with open(linelist_path, errors='replace') as fh:
        for i, ln in enumerate(fh, start=1):
            m = _HDR.match(ln)
            if m and int(m.group(2)) != 0:
                out.append(dict(line_number=i, Z=int(m.group(1)),
                                isotope=int(m.group(2)), ion=int(m.group(3)),
                                n_lines=int(m.group(4))))
    return out


def double_application_offset(Z: int, ion: int, surface: str = VALD_FOR_GRID_SURFACE,
                              record: dict | None = None):
    """Dex the double-application costs this species, or None if it is not exposed.

    Returns the blended-feature offset for co-located atomic isotope structure
    and the per-isotope mapping otherwise — the two are different arithmetics
    and conflating them misstates a trace isotopologue by dex.
    """
    hit = folded_species(surface, record).get((Z, ion))
    if hit is None:
        return None
    if hit['components_colocated']:
        return hit['offset_if_folded_dex']
    return hit['offset_if_folded_dex_per_isotope']


def assert_target_convention(linelist_path: os.PathLike | str, target_z: int,
                             target_ion: int, *, surface: str = VALD_FOR_GRID_SURFACE,
                             record: dict | None = None) -> dict:
    """Preflight a bsyn line list: the FITTED species must not be double-folded.

    Call this immediately before handing `linelist_path` to bsyn.  Returns a
    record describing what was checked (so a harness can bank it as provenance);
    raises SystemExit if the target species is isotope-coded on a folded surface.
    """
    exposed = folded_species(surface, record)
    coded = isotope_coded_blocks(linelist_path)
    offending = [b for b in coded
                 if (b['Z'], b['ion']) == (target_z, target_ion)
                 and (target_z, target_ion) in exposed]
    if offending:
        off = double_application_offset(target_z, target_ion, surface, record)
        raise SystemExit(
            f"RYA-684 ISOTOPE CONVENTION VIOLATION: {linelist_path} fits Z={target_z} "
            f"ion={target_ion} against isotope-coded block(s) "
            f"{[b['isotope'] for b in offending]} sourced from '{surface}', whose log gf "
            f"ALREADY carries the isotope fraction. Turbospectrum will apply it a second "
            f"time (bsyn.f:1350) and the fitted abundance will come out {off} dex high. "
            f"Write the target as Z.000 with the fraction folded in (the RYA-581 Ba "
            f"pattern), or divide the fractions back out of the components.")

    blend_exposure = [dict(Z=b['Z'], ion=b['ion'], isotope=b['isotope'],
                           n_lines=b['n_lines'])
                      for b in coded if (b['Z'], b['ion']) in exposed]
    return dict(
        linelist=str(linelist_path), surface=surface,
        target=dict(Z=target_z, ion=target_ion),
        target_isotope_coded=False,
        blend_blocks_exposed=blend_exposure,
        note=("blend-model blocks on a folded surface are second-order and recorded, "
              "not fatal; RYA-684 measured them at <0.01% of window absorption in "
              "every window feeding a live value"),
    )


def colocated_offset(fractions) -> float:
    """-log10(sum f_i^2): the dex a blended isotope feature is short by."""
    fs = [f for f in fractions if f]
    if not fs:
        raise ValueError("no isotope fractions given")
    return -math.log10(sum(f * f for f in fs))


# ══════════════════════════════════════════════════════════════════════════════
# RYA-1075 — the OTHER half of the class: consumer-side aggregation
# ══════════════════════════════════════════════════════════════════════════════
#
# RYA-684 (above) closed the ENGINE side: a folded gf handed to Turbospectrum on an
# isotope-coded header, so ``isotopfrac`` lands twice.  Offset ``-log10(sum f_i^2)``.
#
# It did NOT close the consumer side, and that is how 54 ``canonical_gf`` rows came to be
# published too high.  A form-(A) source is correct as delivered — isotope-coded, every
# isotope carrying the FULL gf — but anything that aggregates it to "one physical line"
# with ``log10(sum 10**gf)`` over an isotope-BLIND cluster gets ``n_isotopes x`` the gf.
# ``gf_resolver.cluster_physical_lines`` is isotope-blind by construction: it groups on
# species + EP + wavelength gap, and two isotopes of one transition share the lower level
# and sit milli-Angstroms apart.
#
# The offset here is ``log10(n)`` — a COUNT.  It does not depend on the abundances, and it
# is NOT RYA-684's term.  La II settles it: La is 99.911% La-139, so RYA-684's term is
# +0.0008 while a count gives +0.3010.  Nd separates them too (+0.7258 vs +0.8451).  The
# measured offsets were +0.3010 and +0.8451.  Reusing RYA-684's term here would have left
# the defect essentially untouched — which is why this needs its own detector rather than
# an extension of ``double_application_offset``.
#
# ⚠️ n is the number of isotopes IN THE COMPONENT SET, not the element's catalogued
# isotope count.  Ba II 4934 codes 5 isotopes in the GES v6 delivery while ``makeabund.f``
# lists 7, and the measured offset is log10(5).

#: Per-isotope totals must agree within this for the source to be form (A). See
#: ``gf_resolver.ISOTOPE_FORM_A_SPREAD_DEX`` — the same threshold, same rationale.
FORM_A_SPREAD_DEX = 0.02

#: A stored total counts as "built by the naive sum" within this. canonical_gf stores 4
#: decimals; this is that quantum with room.
NAIVE_MATCH_DEX = 0.001


def isotope_inflated_rows(canonical, ges_tsv_path) -> list[dict]:
    """Rows of a canonical gf table that carry ``n_isotopes x`` the physical gf.

    Re-detects the RYA-1075 defect from the SOURCE every time, rather than pinning the 54
    line_ids that were corrected — a pinned list would pass forever while a new ingest
    reintroduced the defect on different lines, which is exactly how RYA-684 came to be
    "Done" with 54 live instances still in the table.

    Returns one dict per offending row. Empty list == clean.
    """
    import numpy as np                                             # noqa: PLC0415
    import pandas as pd                                            # noqa: PLC0415
    from pipeline.gf_resolver import (cluster_physical_lines,      # noqa: PLC0415
                                      species_key)

    src = pd.read_csv(ges_tsv_path, sep='\t', low_memory=False)
    need = {'element', 'ion', 'molecule', 'wave_A', 'lower_state_eV', 'loggf',
            'spectrum_synthe_isotope'}
    missing = need - set(src.columns)
    if missing:
        raise ValueError(f"{ges_tsv_path}: missing columns {sorted(missing)} — the isotope "
                         f"inflation guard cannot run and must not silently pass")

    keys = [species_key(src.element.values[i], src.ion.values[i], src.molecule.values[i])
            for i in range(len(src))]
    wls = src.wave_A.to_numpy(float)
    eps = src.lower_state_eV.to_numpy(float)
    gf = src.loggf.to_numpy(float)
    iso = src.spectrum_synthe_isotope.to_numpy(int)

    canon = canonical[canonical.key_z.astype(str).str.isdigit()].copy()
    canon['_Z'] = canon.key_z.astype(int)

    out: list[dict] = []
    for cl in cluster_physical_lines(keys, wls, eps):
        present = sorted({int(iso[i]) for i in cl if int(iso[i]) != 0})
        if len(present) < 2:
            continue
        per = [float(np.log10(np.sum(10.0 ** gf[[i for i in cl if int(iso[i]) == a]])))
               for a in present]
        if max(per) - min(per) > FORM_A_SPREAD_DEX:
            continue                       # not form (A) — the count rule is undefined
        w = 10.0 ** gf[cl]
        naive = float(np.log10(w.sum()))
        physical = float(np.mean(per))
        centroid = float((wls[cl] * w).sum() / w.sum())
        ep_mean = float(np.mean(eps[cl]))
        z, ion = keys[cl[0]]
        m = canon[(canon._Z == z) & (canon.ion == ion)
                  & ((canon.wavelength_air_A - centroid).abs() <= 0.02)
                  & ((canon.excitation_potential_eV - ep_mean).abs() <= 0.02)]
        if len(m) != 1:
            continue                       # no unique counterpart; nothing to assert about
        r = m.iloc[0]
        if abs(float(r.log_gf) - naive) <= NAIVE_MATCH_DEX and \
                abs(naive - physical) > NAIVE_MATCH_DEX:
            out.append(dict(
                line_id=str(r.line_id), species=str(r.species),
                wavelength_air_A=float(r.wavelength_air_A),
                stored_log_gf=float(r.log_gf), physical_log_gf=round(physical, 4),
                n_isotopes=len(present), isotopes=present,
                inflation_dex=round(naive - physical, 4),
                expected_log10_n=round(math.log10(len(present)), 4)))
    return out
