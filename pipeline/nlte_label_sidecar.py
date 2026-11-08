"""NLTE level labels for the VALD-built synthesis lists, kept BESIDE the list (RYA-1250).

WHY A SIDECAR
-------------
`nearuv_linelist` writes `nlte_label_* = 'none'` on every row it builds, so Gerber NLTE can
reach none of Mg / Al in the IR, H and K lists and no Fe line in K (`assert_linelist_
supports_nlte` refuses; 11 Mg cells FAILED on that). Writing the labels INTO the list would
change `run_matrix._rows_fingerprint` for every cell of every element in those bands -- Fe,
Si, Al and CNO LTE products alike -- although an LTE synthesis never reads a label. So the
labels live in `nlte_labels.tsv` next to `atomic_lines.tsv`, are merged when the list is
loaded, and are fingerprinted only for the NLTE cells of the element they label.

TWO ROUTES, UNIONED PER LEVEL
-----------------------------
* (energy, J) -- `model_atom.resolve_level`, RYA-1094's route. Right for a fine-structure
  atom whose energies agree with VALD's (Al: 173/175; Fe K: 433).
* TERM LABEL -- VALD's own LS designation for each level, from the raw extract. Needed for
  `atom.mg86d`: it is term-resolved above 3p (3d3D g=15, 4f3F g=21) and its energies above
  4d sit ~0.02 eV below VALD's (7p3P* 58322.8 vs 58477 cm-1), so (energy, J) answers ABSENT
  for 95% of Mg I lines while the term is plainly in the atom.

A level is labelled only when the routes name ONE atom level; two different levels is a
contradiction and stays unlabelled. A line is labelled only when BOTH endpoints are (bsyn
needs the pair; RYA-1094). A label match is further required to agree with VALD's energy
within LABEL_DE_EV, so a homonym in a differently-labelled atom cannot slip through.

STAGE ENERGIES ARE MATCHED FROM THE STAGE'S OWN GROUND (RYA-1311)
------------------------------------------------------------------
VALD gives a Ca II level's energy above the Ca II ground; `atom.ca105b` stores it above the
Ca I ground (Ca II 4s 2S at 6.113 eV). Compared raw, no Ca II line could ever match, and every
Ca II list row came back unlabelled although the atom carries 224 Ca II bound-bound
transitions. `_stage_levels` re-zeroes each stage on its lowest level, which is a no-op for an
atom that already stores a stage that way.

A LEVEL'S TERM BELONGS TO THE LEVEL, NOT THE LINE (RYA-728)
-----------------------------------------------------------
A depth-cut VALD extract omits weak lines, so a line can be absent while both its levels appear
on other transitions (Cr II 4848.235: no VALD record, but a4F and z4F* are on dozens). Hence
`read_vald_levels`: (species, energy, J) -> term, from every transition of every extract, kept
only when one term answers. This is a term join keyed on the level; it never picks the nearest
atom energy, which RYA-818 measured naming a wrong level 18% (Cr I) / 35% (Cr II) of the time.

THE GES LIST CAN NEED A SIDECAR TOO (RYA-728)
---------------------------------------------
GES v6 labels most elements well, but its Cr labels were energy-assigned: 0 of the 27 graded
Cr I lines and 1 of 8 Cr II carry both levels, and that one disagrees with VALD's terms. The GES
list lives in the engine install, so its sidecar is committed at `GES_SIDECAR`. It is
AUTHORITATIVE for an element it names (that element's GES labels are cleared first), and it is
applied only when the caller names the cell's own element, so every other element's list, and
therefore its synthesis-cache digest, stays byte-identical.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

SIDECAR_NAME = "nlte_labels.tsv"
ROOT = Path(__file__).resolve().parents[1]
GES_SIDECAR = ROOT / "data" / "linelists" / "ges_v6_nlte_labels" / SIDECAR_NAME
GES_LIST_DIRS = {"GESv6_atom_hfs_iso.420_920nm", "ges_v6"}
#: A term must be this close (eV) to be the same level across VALD extracts (rounding only).
LEVEL_KEY_EV = 0.0015
#: Below this share of the element's graded lines labelled, GES cannot serve its NLTE cells.
GES_MIN_GRADED_SHARE = 0.5
KEY = ["element", "wave_A", "lower_state_eV", "lower_j", "upper_j"]
FIELDS = ["nlte", "nlte_level_low", "nlte_level_up", "nlte_label_low", "nlte_label_up"]
#: A label match must sit this close to VALD's energy. mg86d's largest term offset is
#: 0.020 eV (n=7 triplets); 0.05 still separates every same-label level of every atom here.
LABEL_DE_EV = 0.05
_REC = re.compile(r"^'(\w+ \d)',\s*([-\d.]+),\s*([-\d.]+),\s*([-\d.]+),\s*([-\d.]+),"
                  r"\s*([-\d.]+),\s*([-\d.]+),")
_TERM = re.compile(r"^([a-zA-Z]?)(\d+)([SPDFGHIKLMNOQRT])(\*?)$")


def is_ges_list(list_path) -> bool:
    return Path(list_path).parent.name in GES_LIST_DIRS


def sidecar_path(list_path) -> Path:
    return GES_SIDECAR if is_ges_list(list_path) else Path(list_path).parent / SIDECAR_NAME


def read_vald_terms(raw_path, species: set[str], lo_A: float, hi_A: float) -> pd.DataFrame:
    """(species, wave, E_low, J_low, E_up, J_up, LS designation low/up) from a long-format
    VALD extract. A JK- or jj-coupled level gets '' (no LS term to join on)."""
    out = []
    with open(raw_path, encoding="utf-8", errors="replace") as fh:
        lines = iter(fh)
        for ln in lines:
            m = _REC.match(ln)
            if not m or m.group(1) not in species:
                continue
            w = float(m.group(2))
            lo, up = next(lines, ""), next(lines, "")
            if not lo_A <= w <= hi_A:
                continue
            des = [s.strip().strip("'").strip() for s in (lo, up)]
            terms = [d[2:].strip() if d.startswith("LS") else "" for d in des]
            out.append(dict(species=m.group(1), wave=w, e_low=float(m.group(4)),
                            j_low=float(m.group(5)), e_up=float(m.group(6)),
                            j_up=float(m.group(7)), term_low=terms[0], term_up=terms[1]))
    return pd.DataFrame(out)


def read_vald_levels(raw_paths, species: set[str]) -> dict:
    """{(species, 2J): (energies, terms, designations)} from every transition of every extract.
    A level seen with two different terms within LEVEL_KEY_EV answers nothing (`_level_term`)."""
    t = pd.concat([read_vald_terms(p, species, 0.0, 1e9) for p in raw_paths], ignore_index=True)
    if t.empty:
        return {}
    ends = pd.concat([
        t[["species", "e_low", "j_low", "term_low"]].set_axis(["sp", "e", "j", "des"], axis=1),
        t[["species", "e_up", "j_up", "term_up"]].set_axis(["sp", "e", "j", "des"], axis=1)])
    ends = ends[ends.des != ""].drop_duplicates()
    ends["des"] = ends.des.str.split().str.join(" ")
    out = {}
    for (sp, j), g in ends.groupby(["sp", "j"]):
        out[(sp, int(round(2 * j)))] = (g.e.to_numpy(float), g.des.str.split().str[-1].to_numpy(),
                                        g.des.to_numpy())
    return out


def _level_term(table: dict, species: str, e: float, j: float) -> str:
    """The VALD designation of this level, '' unless exactly one term answers."""
    hit = table.get((species, int(round(2 * j))))
    if hit is None:
        return ""
    m = np.abs(hit[0] - e) < LEVEL_KEY_EV
    return str(hit[2][m][0]) if m.any() and len(set(hit[1][m])) == 1 else ""


def ges_term_join(element: str, levels: pd.DataFrame, nlte_list) -> dict:
    """RYA-818's label join on the GES NLTE list: {(wave mA, ep cdex, ion): (low, up)} for lines
    with both endpoints identified. A key answered two ways is dropped. {} if not applicable."""
    from pipeline import nlte_line_identification as nli
    out, bad = {}, set()
    for ion in (1, 2):
        try:
            d = nli.identify_lines(nli.read_species_lines(nlte_list, element, ion), levels, element, ion)
        except (KeyError, ValueError):
            continue
        for r in d[d.nlte].itertuples():
            k = (int(round(r.wave_A * 1000)), int(round(r.ep_eV * 100)), ion)
            v = (int(r.level_low), int(r.level_up))
            if out.get(k, v) != v:
                bad.add(k)
            out[k] = v
    return {k: v for k, v in out.items() if k not in bad}


def label_candidates(designation: str, j: float) -> set[str]:
    """Atom-label spellings of one VALD LS designation ('3s.7p 3P*', '3d7.(4F).4d f3D',
    '3s2.4s 2S'): Moore letter + term, last orbital + term, full configuration + term,
    each with and without J (an atom may resolve the level or only its term)."""
    parts = designation.split()
    if len(parts) < 2:
        return set()
    conf, term = " ".join(parts[:-1]), parts[-1]
    m = _TERM.match(term)
    if not m:
        return set()
    letter, mult, L, star = m.groups()
    jj = f"{int(j)}" if float(j).is_integer() else f"{j:g}"
    base = f"{mult}{L}"
    out = set()
    orb = re.split(r"[.]", conf)[-1]
    orb = re.sub(r"\d+$", "", orb) if re.search(r"[spdfghik]\d+$", orb) else orb
    for pre in ([letter] if letter else []) + [orb]:
        out |= {f"{pre}{base}{star}", f"{pre}{base}{jj}{star}"}
    full = re.sub(r"\s+", "", f"{conf}{letter}{base}{star}")
    out |= {full, re.sub(r"\s+", "", f"{conf} {letter}{base}{star}")}
    return out


def _norm(s: str) -> str:
    return re.sub(r"\s+", "", str(s))


def _staged_candidates(designation: str, j: float, ion: int) -> set[str]:
    """`atom.ca105b` spells a level stage + last orbital + term (+ J, or 2J when J is half
    integer): '14S21S' Ca I 4s2 1S, '14P3P1' Ca I 4p 3P1, '23D2D3' Ca II 3d 2D3/2, and,
    above its lowest levels, the fine-structure-averaged term alone ('25P2P' Ca II 5p 2P)."""
    parts = designation.split()
    if len(parts) < 2:
        return set()
    m = _TERM.match(parts[-1])
    if not m:
        return set()
    _, mult, L, _ = m.groups()
    orb = re.split(r"[.]", " ".join(parts[:-1]))[-1].upper()
    jj = f"{int(j)}" if float(j).is_integer() else f"{int(round(2 * j))}"
    return {f"{ion}{orb}{mult}{L}", f"{ion}{orb}{mult}{L}{jj}"}


def _stage_levels(levels: pd.DataFrame, ion: int) -> pd.DataFrame:
    """One stage's levels with energies measured from that stage's own ground, as VALD gives them."""
    sub = levels[levels.ion == ion].copy()
    if len(sub):
        sub["energy_eV"] = sub.energy_eV - sub.energy_eV.min()
    return sub


def _by_term(levels: pd.DataFrame, ion: int, designation: str, j: float, e: float):
    if not designation:
        return set()
    cand = {_norm(c) for c in label_candidates(designation, j) | _staged_candidates(designation, j, ion)}
    sub = _stage_levels(levels, ion)
    c = sub[sub.term.map(_norm).isin(cand) & ((sub.energy_eV - e).abs() <= LABEL_DE_EV)]
    return set(int(x) for x in c["index"])


def label_list(list_path, element: str, levels: pd.DataFrame, raw_paths,
               *, tol_eV: float = 0.001, level_paths=None,
               term_join: dict | None = None) -> tuple[pd.DataFrame, dict]:
    """Sidecar rows for `element` in one synthesis list, and a per-ion verdict count.

    `level_paths` adds the level table for lines with no VALD record; `term_join` (from
    `ges_term_join`) adds RYA-818's route. Routes that name different levels leave the line
    unlabelled and count it as `conflict`."""
    from pipeline.model_atom import resolve_level
    d = pd.read_csv(list_path, sep="\t", low_memory=False)
    sp = {f"{element} {i}" for i in (1, 2)}
    d = d[d.element.isin(sp)]
    if d.empty:
        return pd.DataFrame(columns=KEY + FIELDS + ["route"]), {}
    lo, hi = float(d.wave_A.min()) - 1, float(d.wave_A.max()) + 1
    terms = pd.concat([read_vald_terms(p, sp, lo, hi) for p in raw_paths], ignore_index=True)
    level_table = read_vald_levels(level_paths, sp) if level_paths else {}
    rows, stats = [], {}
    for _, r in d.iterrows():
        ion = int(r.element.split()[1])
        t = terms[(terms.species == r.element) & ((terms.wave - r.wave_A).abs() < 0.0015)
                  & ((terms.e_low - r.lower_state_eV).abs() < 0.0006)
                  & ((terms.j_low - r.lower_j).abs() < .01) & ((terms.j_up - r.upper_j).abs() < .01)]
        t = t.iloc[0] if len(t) else None
        ends, routes = [], set()
        for e, j, des in ((r.lower_state_eV, r.lower_j, t.term_low if t is not None else ""),
                          (r.upper_state_eV, r.upper_j, t.term_up if t is not None else "")):
            if not des and level_table:
                des = _level_term(level_table, r.element, float(e), float(j))
                routes |= {"level-table"} if des else set()
            sub = _stage_levels(levels, ion)
            v, idx, _ = resolve_level(sub, energy_eV=float(e), j=float(j), tol_eV=tol_eV)
            hit = {int(idx)} if v == "UNIQUE" else set()
            by_t = _by_term(levels, ion, des, j, float(e))
            # A label can name several J levels of one term (al_qmh '3s2.3p 2P*' x2): the
            # (energy, J) answer picks among them; outside them it is a contradiction.
            if hit and (not by_t or hit <= by_t):
                ends.append(next(iter(hit)))
                routes |= {"energy"} | ({"term"} if by_t else set())
            elif not hit and len(by_t) == 1:
                ends.append(next(iter(by_t)))
                routes.add("term")
            else:
                ends.append(None)
        key = f"{r.element}"
        st = stats.setdefault(key, {"lines": 0, "labelled": 0, "conflict": 0})
        st["lines"] += 1
        joined = (term_join or {}).get((int(round(r.wave_A * 1000)),
                                        int(round(r.lower_state_eV * 100)), ion))
        if joined is not None:
            if any(a is not None and a != b for a, b in zip(ends, joined)):
                st["conflict"] += 1
                continue
            ends = list(joined)
            routes.add("ges-term-join")
        if None in ends:
            continue
        st["labelled"] += 1
        lab = lambda i: str(levels.loc[levels["index"] == i, "term"].iloc[0]).replace(" ", "")
        rows.append({**{k: r[k] for k in KEY}, "nlte": "T", "nlte_level_low": ends[0],
                     "nlte_level_up": ends[1], "nlte_label_low": lab(ends[0]),
                     "nlte_label_up": lab(ends[1]), "route": "+".join(sorted(routes))})
    # A list may carry one transition twice (Al I 10872.973 in IR); apply_sidecar labels
    # every matching row, so one sidecar row per key.
    out = pd.DataFrame(rows, columns=KEY + FIELDS + ["route"]).drop_duplicates()
    return out, stats


def write_sidecar(list_path, frames: list[pd.DataFrame]) -> Path:
    """Write the elements' rows, replacing those elements' previous rows and keeping others'."""
    p = sidecar_path(list_path)
    new = pd.concat([f for f in frames if len(f)], ignore_index=True) if frames else None
    old = pd.read_csv(p, sep="\t") if p.exists() else pd.DataFrame(columns=KEY + FIELDS + ["route"])
    if new is not None and len(new):
        els = set(new.element.str.split().str[0])
        old = old[~old.element.str.split().str[0].isin(els)]
        old = pd.concat([old, new], ignore_index=True)
    old.sort_values(["wave_A", "element"], kind="stable").to_csv(p, sep="\t", index=False)
    return p


def apply_sidecar(linelist, list_path, element: str | None = None):
    """Merge `nlte_labels.tsv` into a loaded iSpec array. A row the list already labels is
    left alone (Fe IR/H from RYA-1094), and a sidecar row matching no list row raises:
    the sidecar must describe THIS list.

    On the GES list (RYA-728) only `element`'s rows apply, after that element's own GES labels
    are cleared; with no element, or an element the sidecar does not name, nothing changes."""
    p = sidecar_path(list_path)
    if not p.exists():
        return linelist, 0
    s = pd.read_csv(p, sep="\t")
    el = np.char.strip(linelist["element"].astype(str))
    if is_ges_list(list_path):
        s = s[s.element.str.split().str[0] == element] if element else s.iloc[0:0]
        if s.empty:
            return linelist, 0
        clear = np.char.startswith(el, f"{element} ")
        linelist["nlte_label_low"][clear] = "none"
        linelist["nlte_label_up"][clear] = "none"
        linelist["nlte_level_low"][clear] = 0
        linelist["nlte_level_up"][clear] = 0
    wl = linelist["wave_A"].astype(float)
    n = 0
    for _, r in s.iterrows():
        hit = np.where((el == r.element) & (np.abs(wl - r.wave_A) < 5e-4)
                       & (np.abs(linelist["lower_state_eV"] - r.lower_state_eV) < 5e-4)
                       & (np.abs(linelist["lower_j"] - r.lower_j) < .01)
                       & (np.abs(linelist["upper_j"] - r.upper_j) < .01))[0]
        if len(hit) == 0:
            raise ValueError(f"{p}: {r.element} {r.wave_A} matches no row of {Path(list_path).name}")
        for i in hit:
            if str(linelist["nlte_label_low"][i]).strip() not in ("", "none"):
                continue
            for f in FIELDS:
                linelist[f][i] = r[f]
            n += 1
    return linelist, n


def ges_graded_share(element: str, ges_list, graded_csvs) -> tuple[int, int]:
    """(graded lines labelled at both ends in the GES list itself, graded lines it holds)."""
    d = pd.read_csv(ges_list, sep="\t", low_memory=False,
                    usecols=["element", "wave_A", "nlte_label_low", "nlte_label_up"])
    d = d[d.element.astype(str).str.split().str[0] == element]
    n = lab = 0
    for c in graded_csvs:
        for w in pd.read_csv(c).wavelength_air_A:
            hit = d[(d.wave_A - float(w)).abs() < 0.02]
            if hit.empty:
                continue
            n += 1
            lab += bool(((hit.nlte_label_low.astype(str) != "none")
                         & (hit.nlte_label_up.astype(str) != "none")).any())
    return lab, n


def fingerprint_rows(list_path, element: str, lo: float, hi: float) -> str | None:
    """sha256 of THIS element's sidecar rows in [lo, hi] -- the NLTE cells' identity only.
    None when it has none there, so a cell the sidecar does not touch keeps its hash."""
    import hashlib
    p = sidecar_path(list_path)
    if not p.exists():
        return None
    s = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
    w = pd.to_numeric(s.wave_A, errors="coerce")
    sub = s[(s.element.str.split().str[0] == element) & (w >= lo) & (w <= hi)]
    if sub.empty:
        return None
    body = "\n".join(sorted(sub.to_csv(index=False, header=False).splitlines()))
    return hashlib.sha256(body.encode()).hexdigest()
