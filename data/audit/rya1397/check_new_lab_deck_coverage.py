"""Check synthesis-deck tuples for unbound laboratory candidates on Sirius."""
import csv
import hashlib
import os
from collections import defaultdict
from pathlib import Path

AUDIT = Path("data/audit/rya1397")
DECKS = {
    "near-UV": Path("data/linelists/ispec_nearuv_3000_3780/atomic_lines.tsv"),
    "VIS": Path(os.environ["ISPEC_DIR"]) /
    "input/linelists/transitions/GESv6_atom_hfs_iso.420_920nm/atomic_lines.tsv",
    "NIR": Path("data/linelists/ispec_ir_9200_13000/atomic_lines.tsv"),
    "K": Path("data/linelists/ispec_k_19452_24846/atomic_lines.tsv"),
}
ROMAN = {"1": "I", "2": "II", "3": "III"}


def similar(left, right, tolerance):
    try:
        return abs(float(left) - float(right)) <= tolerance
    except (TypeError, ValueError):
        return False


def main():
    with (AUDIT / "solar_line_review.csv").open() as stream:
        candidates = [
            r for r in csv.DictReader(stream)
            if r["source_kind"] == "LAB_TABLE" and r["decision"] == "REVIEW_NEW_CANDIDATE"
        ]
    deck_rows = defaultdict(list)
    hashes = {}
    for name, path in DECKS.items():
        hashes[name] = hashlib.sha256(path.read_bytes()).hexdigest()
        with path.open() as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                parts = row["element"].split()
                if len(parts) == 2 and parts[1] in ROMAN:
                    deck_rows[name, f"{parts[0]} {ROMAN[parts[1]]}"].append(row)

    fields = ["target", "species", "band", "source_doi", "source_row", "wavelength_air_A",
              "lower_energy_eV", "upper_energy_eV", "lower_J", "upper_J", "source_loggf",
              "deck_sha256", "deck_wavelengths_A", "deck_loggf", "disposition"]
    with (AUDIT / "new_lab_deck_coverage.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields, lineterminator="\n")
        writer.writeheader()
        for r in candidates:
            band = r["band"]
            deck = "VIS" if band == "red-optical" else band
            found = [line for line in deck_rows[deck, r["species"]]
                     if similar(line["wave_A"], r["wavelength"], .05)
                     and similar(line["lower_state_eV"], r["lower_energy_eV"], .02)
                     and similar(line["upper_state_eV"], r["upper_energy_eV"], .02)
                     and (not r["lower_J"] or similar(line["lower_j"], r["lower_J"], .01))
                     and (not r["upper_J"] or similar(line["upper_j"], r["upper_J"], .01))]
            if deck not in DECKS:
                state = "OUTSIDE_SCOPED_DECKS"
            elif not found:
                state = "NO_CONCORDANT_DECK_TUPLE"
            else:
                state = "DECK_TUPLE_REQUIRES_COMPONENT_AND_STRENGTH_REVIEW"
            match = dict(wavelength_air_A=r["wavelength"], disposition=state,
                         deck_sha256=hashes.get(deck, ""),
                         deck_wavelengths_A=";".join(line["wave_A"] for line in found),
                         deck_loggf=";".join(line["loggf"] for line in found))
            writer.writerow({k: match.get(k, r.get(k, "")) for k in fields})


if __name__ == "__main__":
    main()
