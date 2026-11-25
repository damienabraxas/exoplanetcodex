"""Join held Papoulia Al I/II source levels and compare to canonical transitions."""
import argparse
from collections import defaultdict
import csv
import hashlib
from pathlib import Path
import re

from pipeline.wavelength_util import vac_to_air

AUDIT = Path("data/audit/rya1397")
EV_PER_WAVENUMBER = 1 / 8065.544


def level_key(label):
    parity = "*" if "^o^" in label else ""
    return re.sub(r"[\s^_]", "", label.replace("^o^", "")) + parity


def read_tables(path):
    tables = defaultdict(list)
    name = None
    columns = None
    for line in path.read_text().splitlines():
        if line.startswith("#Table"):
            name = line.split("_")[-1].rstrip(":")
            columns = None
        elif line.startswith("recno\t"):
            columns = line.split("\t")
        elif name and columns and re.match(r"^\s*\d+\t", line):
            tables[name].append(dict(zip(columns, (x.strip() for x in line.split("\t")))))
    return tables


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    args = parser.parse_args()
    tables = read_tables(args.source)
    levels = defaultdict(list)
    for species, table in (("AlI", "table1"), ("AlII", "table8")):
        configuration = ""
        for level in tables[table]:
            configuration = level["Conf"] or configuration
            levels[species, level_key(configuration + level["LSJ"])].append(level)

    canonical = defaultdict(list)
    with Path("data/linelists/canonical_gf.csv").open() as stream:
        for row in csv.DictReader(stream):
            if row["species"] in {"Al I", "Al II"}:
                canonical[row["species"].replace(" ", "")].append(row)

    fields = ["species", "source_row", "source_sha256", "lower_level_pos", "upper_level_pos",
              "observed_ritz_air_A", "published_calculated_A", "candidate_canonical_ids",
              "state"]
    source_hash = hashlib.sha256(args.source.read_bytes()).hexdigest()
    with (AUDIT / "al_papoulia_level_join.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields, lineterminator="\n")
        writer.writeheader()
        for transition in tables["trans"]:
            species = transition["El"]
            low = levels[species, level_key(transition["Lower"])]
            up = levels[species, level_key(transition["Upper"])]
            row = dict(species=species, source_row=transition["recno"], source_sha256=source_hash,
                       published_calculated_A=transition["lambda"], state="HOLD_LEVEL_IDENTITY")
            if len(low) == len(up) == 1:
                lo, = low
                hi, = up
                row.update(lower_level_pos=lo["Pos"], upper_level_pos=hi["Pos"],
                           state="HOLD_OBSERVED_LEVEL_ENERGY")
                if lo["Eobs"].isdigit() and hi["Eobs"].isdigit():
                    delta = int(hi["Eobs"]) - int(lo["Eobs"])
                    if delta > 0:
                        vac = 1e8 / delta
                        if 3000 <= vac <= 25010:
                            wave = float(vac_to_air(vac))
                            ep = int(lo["Eobs"]) * EV_PER_WAVENUMBER
                            # Eobs is rounded to 1 cm^-1. Allow the resulting Ritz wavelength
                            # uncertainty, then require the excitation to agree independently.
                            tol = max(.015, 1e8 / delta**2)
                            ids = sorted({c["physical_id"] for c in canonical[species]
                                          if abs(float(c["wavelength_air_A"]) - wave) <= tol
                                          and abs(float(c["excitation_potential_eV"]) - ep) <= .02
                                          and c["physical_id"]})
                            row.update(observed_ritz_air_A=f"{wave:.4f}",
                                       candidate_canonical_ids=";".join(ids),
                                       state=("THEORY_COMPARISON_REQUIRED" if len(ids) == 1 else
                                              "HOLD_CANONICAL_TUPLE"))
                        else:
                            row["state"] = "OUTSIDE_SCOPED_WAVELENGTHS"
            writer.writerow(row)


if __name__ == "__main__":
    main()
