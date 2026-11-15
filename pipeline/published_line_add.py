"""Add a published solar-literature line that neither canonical_gf nor the synthesis list has
(RYA-1233, governing step 6 "lines secured").

Al's AGSS21 set (Nordlander & Lind 2017) uses Al I 10768.363; the repo had it nowhere, so the
prepare report could only say "run the VALD extraction by hand". That is a step the
orchestrator can take itself, the same way RYA-1214 added AGSS21's missing CNO lines to the
NIR list -- made generic here instead of another per-element script:

  1  the transition, from our own SOLAR VALD extracts (hfs-on: every hyperfine component),
     keyed on wavelength AND lower-level energy (RYA-1037), never wavelength alone
  2  the synthesis list for that wavelength (config: element_prepare.SYNTH_LISTS) gains the
     components, through iSpec's own writer (`nearuv_linelist.to_ispec_array` / `write`,
     round-trip checked)
  3  canonical_gf gains ONE physical row: the authors' total log gf and sigma, tier
     PUBLISHED-SET -- the synthesis scales the components to that total at load
     (`gf_resolver.apply_to_synth_array`)

GUARD (RYA-1214's, generalised): the list is shared by every element's products in that band.
A committed pool line of ANOTHER element within the band's fit half-width of an added
component would move that product, so the add refuses and names it.

The raw VALD extracts left git in RYA-1228; they are read through the registered
`data.vald_raw` path (``CODEX_VALD_RAW`` on Sirius), with the repository staging
and legacy archive as compatibility fallbacks.
"""
from __future__ import annotations

import glob
import io
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CANONICAL = ROOT / "data" / "linelists" / "canonical_gf.csv"
VALD_ARCHIVE = Path(os.environ.get(
    "CODEX_VALD_ARCHIVE", Path.home() / "Documents/Exoplanet Codex/vald_raw_archive/current"))


def _vald_roots() -> tuple[Path, ...]:
    """Return configured VALD locations, in authoritative-to-compatibility order."""
    roots: list[Path] = []
    try:
        from config.constants import codex_path
        roots.append(codex_path("data.vald_raw"))
    except (ImportError, KeyError, OSError):
        pass
    # Keep explicit overrides useful for older host profiles and local mirrors.
    for raw in (os.environ.get("CODEX_VALD_RAW"), os.environ.get("CODEX_VALD_ARCHIVE")):
        if raw:
            roots.append(Path(raw).expanduser())
    roots.extend((ROOT / "data" / "linelists", VALD_ARCHIVE))
    return tuple(dict.fromkeys(roots))
#: Our solar VALD3 Extract Stellar deliveries with hyperfine structure ON, by span (A, air).
SOLAR_EXTRACTS = ((2000.0, 3780.0, "vald_solar_nearuv_2000_3780_hfson_raw.txt"),
                  (6910.0, 9500.0, "vald_solar_redopt_6910_9500_hfson_raw.txt"),
                  (9500.0, 17000.0, "vald_solar_ir_9500_17000_hfson_raw.txt"),
                  (17000.0, 25000.0, "vald_solar_ir_17000_25000_hfson_raw.txt"))
#: Hyperfine components of one transition spread < 0.1 A; the published wavelength is the
#: centroid, printed to 0.001-0.01 A.
COMPONENT_WINDOW_A = 0.10
EP_TOL_EV = 0.01
TIER = "PUBLISHED-SET"
SCOPED_ADDITIONS = "element_scoped_additions.tsv"
UNSOLD_FDAMP = {"Na": 2.0}           # default 2.5 (RYA-1232, GESv6)


class LineAddError(RuntimeError):
    pass


def extract_for(w: float) -> Path:
    for lo, hi, name in SOLAR_EXTRACTS:
        if lo <= w < hi:
            searched = _vald_roots()
            for root in searched:
                for p in (root / name, root.parent / name):
                    if p.exists():
                        return p
            raise LineAddError(
                f"{name} (covers {w:.3f} A) not found; searched VALD roots: "
                + ", ".join(str(p) for p in searched))
    raise LineAddError(f"no solar hfs-on VALD extract covers {w:.3f} A")


def components(species: str, w: float, ep: float) -> tuple[list[dict], Path]:
    """Every VALD component of the transition at (w, ep) -- one lower level, one window."""
    from pipeline import nearuv_linelist as nl
    path = extract_for(w)
    band, _ = nl.read_band(path, w - 1.0, w + 1.0)
    el, ion = species.split()
    hit = [r for r in band if str(r["element"]).strip() == el and str(r["ion"]).strip() == ion
           and abs(float(r["wavelength"]) - w) <= COMPONENT_WINDOW_A
           and abs(float(r["e_low_eV"]) - ep) <= EP_TOL_EV]
    if not hit:
        raise LineAddError(f"{species} {w:.3f} (EP {ep}): no component in {path.name}")
    if len({round(float(r["e_up_eV"]), 3) for r in hit}) != 1:
        raise LineAddError(f"{species} {w:.3f}: components reach {len({round(float(r['e_up_eV']), 3) for r in hit})} "
                           f"upper levels -- two transitions in the window, refusing to pick")
    return hit, path


def _band_of(w: float):
    from config.synth_bands import SYNTH_BANDS
    return next((b for b in SYNTH_BANDS.values() if b.lo_A <= w < b.hi_A), None)


def clashes(element: str, waves: list[float]) -> list[str]:
    """Committed pool lines of OTHER elements within the band half-width of an added line."""
    out = []
    for f in glob.glob(str(ROOT / "data/results/band_products/*_lines.csv")):
        stem = Path(f).name
        sp = stem.split("_")[0]                        # e.g. FeI, SiII, AlI
        if sp.rstrip("I") == element:
            continue
        d = pd.read_csv(f, usecols=lambda c: c in ("wavelength_air_A", "in_aggregate"))
        if "wavelength_air_A" not in d:
            continue
        if "in_aggregate" in d:
            d = d[d.in_aggregate.astype(str).str.lower() == "true"]
        for x in d.wavelength_air_A.astype(float):
            b = _band_of(x)
            if b and any(abs(x - w) <= b.half_width_A for w in waves):
                out.append(f"{sp} {x:.3f} ({stem})")
    return sorted(set(out))


def plan(row: dict, *, allow_ungraded: bool = False) -> dict:
    """What adding one missing published line would write (no files touched)."""
    sigma = row.get('authors_sigma_dex')
    if not allow_ungraded and (sigma is None or not np.isfinite(float(sigma))):
        raise LineAddError("Published line has no resolved gf uncertainty; cannot admit as graded")
    from pipeline.element_prepare import SYNTH_LISTS
    sp, w, ep = row["species"], float(row["wavelength_A"]), float(row["ep_eV"])
    comps, src = components(sp, w, ep)
    gf = np.array([float(r["log_gf"]) for r in comps])
    wl = np.array([float(r["wavelength"]) for r in comps])
    # RYA-1250: one canonical row is ONE synthesis cluster (gf_resolver.HFS_GAP). Mg II 10092.20
    # is two components 0.122 A apart: one row would scale one of them to the feature total and
    # leave the other unresolved -- ~0.25 dex too much opacity. Register the components instead.
    from pipeline.gf_resolver import HFS_GAP
    if len(wl) > 1 and np.max(np.diff(np.sort(wl))) > HFS_GAP:
        raise LineAddError(f"{sp} {w:.3f}: components span more than one synthesis cluster "
                           f"(gap > {HFS_GAP} A) -- register each component, not the feature")
    centroid = float(np.sum(wl * 10 ** gf) / np.sum(10 ** gf))
    lst = None
    for lo, f in SYNTH_LISTS:
        if w >= lo:
            lst = f
    if lst is None:
        raise LineAddError(f"{sp} {w:.3f}: its synthesis list is iSpec's own GES list, outside "
                           f"the repo -- not writable here")
    # A published transition can already be present in the shared synthesis list while
    # missing only from canonical_gf (common for IR VALD lists).  Keep the canonical
    # admission path generic, but do not duplicate the physical components in that list.
    existing = row.get("in_synthesis_list")
    return {"species": sp, "wavelength_A": w, "ep_eV": ep, "set": row["set"],
            "authors_loggf": row["authors_loggf"], "authors_sigma_dex": row["authors_sigma_dex"],
            "authors_gf_source": row.get("authors_gf_source"),
            "extract": src.name, "n_components": len(comps), "centroid_A": round(centroid, 4),
            "vald_total_loggf": round(float(np.log10(np.sum(10 ** gf))), 3),
            "synthesis_list": existing or str(lst.relative_to(ROOT)),
            "synthesis_existing": bool(existing),
            "clashes": clashes(sp.split()[0], list(wl)), "_components": comps}


def scoped_additions_path(synthesis_list: Path) -> Path:
    return Path(synthesis_list).parent / SCOPED_ADDITIONS


def read_scoped_additions(synthesis_list: Path) -> pd.DataFrame:
    """Rows added to a shared list for one element (RYA-728). `clashes` has already shown they
    cannot move another element's product, so they do not count in other elements' cell identity."""
    p = scoped_additions_path(synthesis_list)
    if not p.exists():
        return pd.DataFrame(columns=["element", "ion_element", "wave_A"])
    return pd.read_csv(p, sep="\t", dtype={"element": str, "ion_element": str, "wave_A": float})


def _append_canonical_rows(rows: list[dict], columns: list[str]) -> None:
    """Append text only. A full to_csv reformats floats in rows nobody touched (RYA-728: 157 Ti I
    and Ca I lines changed text and no value), and those bytes key finished elements' cells."""
    buf = io.StringIO(newline="")
    pd.DataFrame(rows).reindex(columns=columns).to_csv(buf, index=False, header=False)
    with CANONICAL.open("a", encoding="utf-8", newline="") as fh:
        fh.write(buf.getvalue())


def apply(plans: list[dict], *, citation: dict, tier: str | None = None, status: str | None = None,
          element_scoped: bool = False, at_centroid: bool = False,
          accept_clash: tuple = ()) -> list[dict]:
    """Write each plan (synthesis list + canonical row). Refuses on any clash.

    `tier` and `status` name a lab line's own grade rather than the published-set defaults.
    `element_scoped` records the list rows, so other elements' cells keep their identity.
    `at_centroid` places the canonical row on VALD's wavelength: a lab table's 2-decimal value
    can sit 0.09 A from the transition (Sobeck 9667.2 is VALD 9667.286), and the fit centres there.
    `accept_clash` names lines a person decided to add next to another element's pool line (Ryan
    2026-11-08: a line leaves only when it cannot be resolved). The caller records the clash."""
    from pipeline import nearuv_linelist as nl
    from pipeline.element_prepare import adopted_reference
    from pipeline.physical_line_id import physical_id
    bad = [p for p in plans if p["clashes"] and p["wavelength_A"] not in accept_clash]
    if bad:
        raise LineAddError("; ".join(f"{p['species']} {p['wavelength_A']}: would move "
                                     f"{', '.join(p['clashes'][:5])}" for p in bad))
    nl._ensure_ispec_on_path()
    import ispec
    by_list: dict = {}
    for p in plans:
        if not p.get("synthesis_existing"):
            by_list.setdefault(p["synthesis_list"], []).extend(p["_components"])
    for rel, recs in by_list.items():
        path = ROOT / rel
        base = ispec.read_atomic_linelist(str(path))
        add = nl.to_ispec_array(recs)
        # fdamp 0.0 is NOT a default: our Turbospectrum DROPS the line (RYA-1232, measured).
        # VALD gives no vdW for these -> Unsöld x 2.5, GESv6's own factor (Na I 2.0), the
        # value RYA-1232's nearuv_linelist.unsold_fdamp writes.
        z = add["turbospectrum_fdamp"] == 0.0
        add["turbospectrum_fdamp"][z] = [UNSOLD_FDAMP.get(str(e).split()[0], 2.5)
                                         for e in add["element"][z]]
        merged = np.concatenate([base, add.astype(base.dtype)])
        nl.write(merged[np.argsort(merged["wave_A"], kind="stable")], path)
        if element_scoped:
            keep = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False,
                               usecols=["element", "wave_A"])
            keep["w"] = keep.wave_A.astype(float)
            added = np.round(np.asarray(add["wave_A"], dtype=float), 4)
            names = [str(e).strip() for e in add["element"]]
            mine = keep[[any(n == e and abs(a - w) < 1e-3 for n, a in zip(names, added))
                         for e, w in zip(keep.element, keep.w)]]
            row = pd.DataFrame({"element": mine.element.str.split().str[0], "ion_element": mine.element,
                                "wave_A": mine.wave_A})
            sp = scoped_additions_path(path)
            old = read_scoped_additions(path).astype(str) if sp.exists() else None
            pd.concat([old, row] if old is not None else [row]).to_csv(sp, sep="\t", index=False)
    canon = pd.read_csv(CANONICAL, low_memory=False)
    start = 1 + max(int(x.split("_")[1]) for x in canon.line_id if str(x).startswith("gf_"))
    rows = []
    for i, p in enumerate(plans):
        el, ion = p["species"].split()
        same = canon[canon.species == p["species"]]
        z = int(same.key_z.iloc[0]) if len(same) else None
        wave = round(p["centroid_A"] if at_centroid else p["wavelength_A"], 3)
        cit = citation.get(p["set"], {})
        rows.append({
            "line_id": f"gf_{start + i:06d}",
            "physical_id": physical_id(p["species"], wave, round(p["ep_eV"], 4)),
            "key_z": z, "ion": 1.0 if ion == "I" else 2.0, "species": p["species"],
            "wavelength_air_A": wave,
            "excitation_potential_eV": round(p["ep_eV"], 4),
            "hfs_n_components": p["n_components"],
            "log_gf": p["authors_loggf"], "loggf_reference": adopted_reference(p, citation),
            "seed_source": f"{p['extract']} (orchestrator prepare --apply, RYA-1233)",
            "adjudication_status": status or f"added_rya1233_{p['set'].lower()}",
            "gf_linelist_vald": p["vald_total_loggf"], "in_synth": True, "in_linelist": True,
            "lab_source_tag": p["set"], "gf_sigma_dex": p["authors_sigma_dex"],
            "gf_source_doi": cit.get("doi"), "is_diagnostic": False, "gf_tier": tier or TIER})
    _append_canonical_rows(rows, list(canon.columns))
    return [{k: v for k, v in p.items() if not k.startswith("_")} for p in plans]
