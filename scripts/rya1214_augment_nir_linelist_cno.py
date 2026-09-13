#!/usr/bin/env python3
"""RYA-1214 — add AGSS21's missing CNO indicator lines to the NIR synthesis list.

`data/linelists/ispec_ir_9200_13000/atomic_lines.tsv` (RYA-762) was built from a depth-cut
VALD extract and holds 34 C I lines, ZERO N I and ZERO O I. So an AGSS21-set run on the NIR
band refused with "no 'O 1' rows in the synthesis list" — not because the lines are absent
from our holdings, but because the band's list was never built for CNO.

WHAT IS ADDED — exactly the AGSS21 reference-set lines the list lacks, nothing else:
N I 10108.892 and the nine O I 926.1 / 926.3 / 926.6 nm components. They are weak
(catalogued depth 0.002-0.038), which is why the depth cut dropped them.

WHERE THE ROWS COME FROM — our own solar VALD hfs-on extracts, which carry J and upper-level
energies (`vald_solar_redopt_6910_9500_hfson_raw.txt` for O I 926,
`vald_solar_ir_9500_17000_hfson_raw.txt` for N I 10108.9), converted by the SAME
`pipeline.nearuv_linelist.to_ispec_array` / `write` the near-UV list uses (iSpec's own
writer, round-trip checked). gf is not decided here: the synthesis route applies canonical
gf at load, so these lines take their NIST-adjudicated values (RYA-1214 Step 1).

🔴 GATE: IT CANNOT MOVE AN FE PRODUCT. The list is shared with every Fe NIR product. Before
writing, every committed Fe I NIR pool line is checked against the added lines at the NIR
fit half-width (config.synth_bands NIR, 1.40 A); any overlap refuses. Measured: zero.
Must run where iSpec is importable (Sirius).
"""
from __future__ import annotations

import glob
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline import nearuv_linelist as nl  # noqa: E402

LIST = ROOT / "data/linelists/ispec_ir_9200_13000/atomic_lines.tsv"
SETS = ROOT / "data/linelists/reference_sets"
VALD = {"O": ROOT / "data/linelists/vald_solar_redopt_6910_9500_hfson_raw.txt",
        "N": ROOT / "data/linelists/vald_solar_ir_9500_17000_hfson_raw.txt"}
AUDIT = ROOT / "data/audit/rya1214_cno_products"
LO_A, HI_A = 9199.0, 13000.0
MATCH_TOL_A = 0.05          # the reference sets' own tolerance (half AGSS21's printed precision)
#: Identifying ONE VALD component for a set row. The O I 926.1 nm feature is three components
#: 0.04 A apart (9260.80/.84/.93), so 0.05 A matches two of them and the unique-match guard
#: refused. The bound that CANNOT match two components is half the closest spacing, 0.02 A;
#: 0.015 A leaves margin for the set's 0.01 A printing (set 9262.66 vs VALD 9262.670 sits
#: exactly on 0.010, so a 0.01 A bound failed on float rounding). Uniqueness still required.
COMPONENT_TOL_A = 0.015


def _half_width() -> float:
    from config.synth_bands import SYNTH_BANDS
    return float(SYNTH_BANDS["NIR"].half_width_A)


def main() -> int:
    nl._ensure_ispec_on_path()
    import ispec
    base = ispec.read_atomic_linelist(str(LIST))
    have = pd.DataFrame({"element": base["element"], "wave_A": base["wave_A"]})

    want = []
    manifest = json.loads((SETS / "agss21_cno_lineset_rya1214.prov.json").read_text())
    for f in manifest["files"]:
        s = pd.read_csv(ROOT / f["file"])
        s = s[s.wavelength_air_A.between(LO_A, HI_A)]
        for _, r in s.iterrows():
            sp = f"{r.element} 1"
            if not ((have.element == sp) & ((have.wave_A - r.wavelength_air_A).abs()
                                             <= MATCH_TOL_A)).any():
                want.append((r.element, float(r.wavelength_air_A), r.line_label))
    if not want:
        print("every AGSS21 set line in the NIR band is already in the list — nothing to do")
        return 0

    recs, reports = [], {}
    for el in sorted({w[0] for w in want}):
        band, rep = nl.read_band(VALD[el], LO_A, HI_A)
        reports[el] = {k: rep[k] for k in ("n_parsed", "n_in_band", "source") if k in rep}
        for e, w, lab in [x for x in want if x[0] == el]:
            # parse_vald_long splits the species: element 'N', ion 'I' ('species' is 'N 1').
            hit = [r for r in band if str(r["element"]).strip() == e
                   and str(r["ion"]).strip() == "I"
                   and abs(float(r["wavelength"]) - w) <= COMPONENT_TOL_A]
            if len(hit) != 1:
                raise SystemExit(f"{e} {lab} {w}: {len(hit)} VALD rows within {COMPONENT_TOL_A} "
                                 f"A in {VALD[el].name} — refusing to pick one")
            if any(hit[0] is r for r in recs):
                raise SystemExit(f"{e} {lab} {w}: VALD row already taken by another set row")
            recs.append(hit[0])

    hw = _half_width()
    fe = []
    for f in glob.glob(str(ROOT / "data/results/band_products/FeI_*_lines.csv")):
        d = pd.read_csv(f)
        if "wavelength_air_A" in d:
            fe += [float(x) for x in d.wavelength_air_A if LO_A <= float(x) <= HI_A]
    clash = sorted({round(w, 3) for w in fe for r in recs
                    if abs(w - float(r["wavelength"])) <= hw})
    if clash:
        raise SystemExit(f"Fe I NIR pool lines within {hw} A of an added line: {clash} — "
                         f"adding opacity there would move a committed Fe product. Refusing.")

    add = nl.to_ispec_array(recs)
    merged = np.concatenate([base, add.astype(base.dtype)])
    merged = merged[np.argsort(merged["wave_A"], kind="stable")]
    nl.write(merged, LIST)
    print(f"added {len(add)} line(s) to {LIST.relative_to(ROOT)}: {len(base)} -> {len(merged)}")
    for r in recs:
        print(f"  {r['species']} {float(r['wavelength']):.3f}")

    AUDIT.mkdir(parents=True, exist_ok=True)
    (AUDIT / "nir_linelist_cno_augment.prov.json").write_text(json.dumps({
        "ticket": "RYA-1214",
        "list": str(LIST.relative_to(ROOT)),
        "why": "RYA-762's depth-cut VALD extract left the NIR list with 0 N I / 0 O I; the "
               "AGSS21 set runs refused for lack of synthesis rows",
        "added": [{"species": str(r["species"]).strip(), "wavelength_air_A":
                   float(r["wavelength"])} for r in recs],
        "n_before": int(len(base)), "n_after": int(len(merged)),
        "source_extracts": reports,
        "row_converter": "pipeline.nearuv_linelist.to_ispec_array + write (iSpec writer, "
                         "round-trip checked)",
        "gf": "canonical gf applied by the synthesis route at load (NIST-adjudicated)",
        "fe_guard": {"half_width_A": hw, "fe_nir_pool_lines_checked": len(fe),
                     "within_half_width": 0},
    }, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
