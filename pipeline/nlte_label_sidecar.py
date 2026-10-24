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
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

SIDECAR_NAME = "nlte_labels.tsv"
KEY = ["element", "wave_A", "lower_state_eV", "lower_j", "upper_j"]
FIELDS = ["nlte", "nlte_level_low", "nlte_level_up", "nlte_label_low", "nlte_label_up"]
#: A label match must sit this close to VALD's energy. mg86d's largest term offset is
#: 0.020 eV (n=7 triplets); 0.05 still separates every same-label level of every atom here.
LABEL_DE_EV = 0.05
_REC = re.compile(r"^'(\w+ \d)',\s*([-\d.]+),\s*([-\d.]+),\s*([-\d.]+),\s*([-\d.]+),"
                  r"\s*([-\d.]+),\s*([-\d.]+),")
_TERM = re.compile(r"^([a-zA-Z]?)(\d+)([SPDFGHIKLMNOQRT])(\*?)$")


def sidecar_path(list_path) -> Path:
    return Path(list_path).parent / SIDECAR_NAME


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


def _by_term(levels: pd.DataFrame, ion: int, designation: str, j: float, e: float):
    if not designation:
        return set()
    cand = {_norm(c) for c in label_candidates(designation, j)}
    c = levels[(levels.ion == ion) & levels.term.map(_norm).isin(cand)
               & ((levels.energy_eV - e).abs() <= LABEL_DE_EV)]
    return set(int(x) for x in c["index"])


def label_list(list_path, element: str, levels: pd.DataFrame, raw_paths,
               *, tol_eV: float = 0.001) -> tuple[pd.DataFrame, dict]:
    """Sidecar rows for `element` in one synthesis list, and a per-ion verdict count."""
    from pipeline.model_atom import resolve_level
    d = pd.read_csv(list_path, sep="\t", low_memory=False)
    sp = {f"{element} {i}" for i in (1, 2)}
    d = d[d.element.isin(sp)]
    if d.empty:
        return pd.DataFrame(columns=KEY + FIELDS + ["route"]), {}
    lo, hi = float(d.wave_A.min()) - 1, float(d.wave_A.max()) + 1
    terms = pd.concat([read_vald_terms(p, sp, lo, hi) for p in raw_paths], ignore_index=True)
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
            sub = levels[levels.ion == ion]
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
        st = stats.setdefault(key, {"lines": 0, "labelled": 0})
        st["lines"] += 1
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


def apply_sidecar(linelist, list_path):
    """Merge `nlte_labels.tsv` into a loaded iSpec array. A row the list already labels is
    left alone (Fe IR/H from RYA-1094), and a sidecar row matching no list row raises:
    the sidecar must describe THIS list."""
    p = sidecar_path(list_path)
    if not p.exists():
        return linelist, 0
    s = pd.read_csv(p, sep="\t")
    el = np.char.strip(linelist["element"].astype(str))
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
