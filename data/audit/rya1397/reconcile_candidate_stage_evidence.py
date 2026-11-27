import csv
from collections import defaultdict
from pathlib import Path
from math import isfinite

audit = Path("data/audit/rya1397")
with (audit / "solar_line_review.csv").open() as source:
    rows = list(csv.DictReader(source))
ids = defaultdict(list)
for row in rows:
    key = row["canonical_candidate_ids"]
    if len(key.split()) == 1 and key.startswith("pk_") and row["decision"] != "PRESERVE_BASELINE":
        ids[key].append(row)
index = defaultdict(list)
with Path("data/linelists/canonical_gf.csv").open() as source:
    for row in csv.DictReader(source):
        try:
            w = float(row["wavelength_air_A"])
            ep = float(row["excitation_potential_eV"])
        except ValueError:
            continue
        if isfinite(w) and isfinite(ep) and row["physical_id"]:
            index[row["species"], round(w * 100)].append((row["physical_id"], w, ep))
hits = defaultdict(lambda: defaultdict(set))
root = Path("data/results/orchestrator/budget")
for path in sorted(root.glob("solar_*/stage/*_lines.csv")):
    with path.open() as source:
        for line in csv.DictReader(source):
            try:
                w = float(line["wavelength_air_A"])
                ep = float(line["ep_eV"])
            except (KeyError, ValueError, TypeError):
                continue
            if not (isfinite(w) and isfinite(ep)):
                continue
            species = f'{line["element"]} {line["ion"]}'
            centers = range(round(w * 100) - 1, round(w * 100) + 2)
            matches = {key for center in centers for key, cw, ce in index[species, center]
                       if abs(w - cw) <= .005 and abs(ep - ce) <= .02}
            if len(matches) != 1:
                continue
            key, = matches
            flag = line["in_aggregate"].lower()
            if key not in ids or flag not in {"true", "false"}:
                continue
            state = "accepted" if flag == "true" else "excluded"
            hits[key][state].add(str(path))
fields = ("target", "species", "canonical_id", "decisions_in_audit", "accepted_stage_files",
          "excluded_stage_files", "accepted_example", "excluded_example")
with (audit / "candidate_stage_evidence.csv").open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fields, lineterminator="\n")
    writer.writeheader()
    for key, evidence in sorted(hits.items(), key=lambda item: (ids[item[0]][0]["target"], item[0])):
        a, e = evidence["accepted"], evidence["excluded"]
        writer.writerow(dict(target=ids[key][0]["target"], species=ids[key][0]["species"],
                             canonical_id=key,
                             decisions_in_audit=";".join(sorted({row["decision"] for row in ids[key]})),
                             accepted_stage_files=len(a), excluded_stage_files=len(e),
                             accepted_example=min(a, default=""), excluded_example=min(e, default="")))
