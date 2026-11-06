#!/usr/bin/env python3
"""RYA-1250: write `nlte_labels.tsv` beside every VALD-built synthesis list, for every
element with a Gerber deck (`gerber_nlte.DECKS`). See `pipeline/nlte_label_sidecar.py`.

    python scripts/rya1250_nlte_label_sidecars.py            # report only
    python scripts/rya1250_nlte_label_sidecars.py --apply    # write the sidecars
    python scripts/rya1250_nlte_label_sidecars.py --apply --element Ti   # one element's rows

RUNS ON SIRIUS (model atoms and raw VALD extracts are on the data volumes).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from config.constants import codex_path                       # noqa: E402
from pipeline import nlte_label_sidecar as nls                # noqa: E402
from pipeline.gerber_nlte import DECKS                        # noqa: E402
from pipeline.model_atom import read_gerber_atom              # noqa: E402

#: list dir -> the raw extracts it was built from (nearuv_linelist / RYA-762 / RYA-1094).
LISTS = {"ispec_nearuv_3000_3780": ["vald_solar_nearuv_2000_3780_hfson_raw.txt"],
         "ispec_ir_9200_13000": ["vald_solar_redopt_6910_9500_hfson_raw.txt",
                                 "vald_solar_ir_9500_17000_hfson_raw.txt"],
         "ispec_h_15007_17494": ["vald_solar_ir_9500_17000_hfson_raw.txt"],
         "ispec_k_19452_24846": ["vald_solar_ir_17000_25000_hfson_raw.txt"]}
#: Labelled IN the list by RYA-1094 and feeding published Fe products. The routes here
#: reproduce those labels (IR 191/191, H 498/502, 0 disagreements) and would add 12 IR and
#: 296 H lines -- a change to finished Fe work, so not made here (Ryan's call).
SKIP = {("ispec_ir_9200_13000", "Fe"), ("ispec_h_15007_17494", "Fe")}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--element", action="append",
                    help="label only these elements (others' rows are kept as they are)")
    a = ap.parse_args()
    raw = Path(str(codex_path("data.vald_raw")))
    atoms = Path(str(codex_path("grids.gerber_ts")))
    els = sorted({k.split("@")[0]: v["atom"] for k, v in DECKS.items()}.items())
    if a.element:
        els = [(el, atom) for el, atom in els if el in a.element]
    report = {}
    for name, extracts in LISTS.items():
        lst = ROOT / "data" / "linelists" / name / "atomic_lines.tsv"
        if not lst.exists():
            print(f"{name}: list not built here -- skipped")
            continue
        frames = []
        for el, atom in els:
            if (name, el) in SKIP:
                continue
            lv = read_gerber_atom(atoms / atom)
            rows, st = nls.label_list(lst, el, lv, [raw / e for e in extracts])
            frames.append(rows)
            report[f"{name}/{el}"] = {**{k: v for k, v in st.items()},
                                      "routes": rows.route.value_counts().to_dict()}
            print(name, el, st, rows.route.value_counts().to_dict())
        if a.apply:
            print("  wrote", nls.write_sidecar(lst, frames).relative_to(ROOT))
    out = ROOT / "data" / "audit" / "rya1250" / (
        "nlte_label_sidecars.json" if not a.element else f"nlte_label_sidecars_{'_'.join(a.element)}.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=1, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
