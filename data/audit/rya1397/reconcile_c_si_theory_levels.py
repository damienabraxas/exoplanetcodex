"""Join the held carbon and silicon theory candidates to published level tables."""
import argparse
from collections import defaultdict
import csv
from itertools import islice
from pathlib import Path
import re

from pipeline.wavelength_util import vac_to_air

AUDIT = Path("data/audit/rya1397")
CARBON = Path("data/reference/cno_atomic_primary/c_li2021/CI-IV_supplementary_tables")


def carbon_key(label):
    label = re.sub(r"\([^)]*\)|~\^\{2\}P(?=\\,)", "", label)
    return re.sub(r"[\s\$\\,{}~^_]", "", label.replace(r"\circ", "o"))


def silicon_key(label):
    return re.sub(r"\s+", "", label)


def read_silicon_levels(path):
    levels = defaultdict(list)
    with path.open() as stream:
        columns = next(line for line in stream if line.startswith("recno\tLevel"))
        for line in islice(stream, 2, None):
            if not line.strip() or line.startswith("#"):
                break
            row = dict(zip(columns.strip().split("\t"), (s.strip() for s in line.split("\t"))))
            levels[silicon_key(row["Level"])].append(row)
    return levels


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("silicon_source_directory", type=Path)
    args = parser.parse_args()
    carbon = defaultdict(list)
    with (CARBON / "tableA1_level.tex").open() as stream:
        for line in stream:
            row = [s.strip() for s in line.split("&")]
            if row[0] == "C~I":
                carbon[carbon_key(row[2])].append((row[1], row[5].replace("~", "").strip()))
    silicon = {ion: read_silicon_levels(args.silicon_source_directory / table)
               for ion, table in (("Si I", "tableb1.tsv"), ("Si II", "tableb5.tsv"))}
    canonical = defaultdict(list)
    with Path("data/linelists/canonical_gf.csv").open() as stream:
        for row in csv.DictReader(stream):
            species = row["species"]
            if species in {"C I", "Si I", "Si II"} and row["physical_id"]:
                canonical[species].append((float(row["wavelength_air_A"]),
                                           float(row["excitation_potential_eV"]), row["physical_id"]))
    fields = ["target", "species", "source_path", "source_row", "source_kind", "lower_level_id",
              "upper_level_id", "observed_lower_cm-1", "observed_upper_cm-1", "ritz_air_A", "canonical_candidates", "state"]
    with (AUDIT / "solar_line_review.csv").open() as source, \
            (AUDIT / "c_si_source_level_join.csv").open("w", newline="") as dest:
        writer = csv.DictWriter(dest, fields, lineterminator="\n")
        writer.writeheader()
        for candidate in csv.DictReader(source):
            species = candidate["species"]
            theory = candidate["source_kind"] == "THEORY_CANDIDATE"
            si_lab = candidate["source_kind"] == "ACQUIRED_LAB_TABLE"
            if species == "C I" and theory:
                lo = carbon.get(carbon_key(candidate["lower_level"]), [])
                hi = carbon.get(carbon_key(candidate["upper_level"]), [])
            elif species in silicon and (theory or si_lab):
                table = silicon[species]
                column = "Energy5" if species == "Si I" else "Energy4"
                lo = [(r["recno"], r[column]) for r in table.get(silicon_key(candidate["lower_level"]), [])]
                hi = [(r["recno"], r[column]) for r in table.get(silicon_key(candidate["upper_level"]), [])]
            else:
                continue
            row = {k: candidate[k] for k in fields if k in candidate}
            row["state"] = "HOLD_SOURCE_LEVEL_JOIN"
            if len(lo) == len(hi) == 1:
                (lo_id, lo_energy), = lo
                (hi_id, hi_energy), = hi
                row.update(lower_level_id=lo_id, upper_level_id=hi_id,
                           state="HOLD_OBSERVED_ENERGY")
                try:
                    gap = float(hi_energy) - float(lo_energy)
                except ValueError:
                    gap = 0
                if gap > 0:
                    vac = 1e8 / gap
                    row.update(state="OUTSIDE_SCOPED_WAVELENGTHS",
                               **{"observed_lower_cm-1": lo_energy, "observed_upper_cm-1": hi_energy})
                    if 3000 <= vac <= 25010:
                        air = float(vac_to_air(vac))
                        # Carbon level energies are printed at 1 cm^-1 precision.
                        wave_tol = max(.015, (1.0 if species == "C I" else .02) * 1e8 / gap**2)
                        ep = float(lo_energy) / 8065.544
                        ids = sorted({key for wave, excitation, key in canonical[species]
                                      if abs(wave - air) <= wave_tol and abs(ep - excitation) <= .02})
                        row.update(ritz_air_A=f"{air:.4f}", canonical_candidates=";".join(ids),
                                   state=("PROVISIONAL_CANONICAL_TUPLE" if len(ids) == 1 else
                                          "HOLD_CANONICAL_TUPLE"))
            writer.writerow(row)


if __name__ == "__main__":
    main()
