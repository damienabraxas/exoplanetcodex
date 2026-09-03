#!/usr/bin/env python3
"""Cross-identify the Amarsi+2021 molecular abundance lines in primary lists.

The source archives deliberately remain the authority.  This script emits only
the compact transition evidence needed by RYA-1136; it does not vendor a second,
lossy copy of several hundred thousand molecular transitions.
"""
from __future__ import annotations

import csv
import bisect
import gzip
import io
import itertools
import json
import math
import zipfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PRIMARY = ROOT / "data/reference/cno_molecular_primary"
AMARSI = ROOT / "data/reference/amarsi2021_cno/derived/amarsi2021_cno_molecular_lines.csv"
OUT = ROOT / "data/audit/rya1136_cno_intake/primary_molecular_crossmatch.csv"
HC_EV_CM = 1.0 / 8065.544005
# Brooke's Swan E'' values are relative to the a-state v=0 band origin,
# whereas Amarsi tabulates excitation from its adopted molecular zero point.
# The 0.0753 eV origin shift is independently visible across all 39 rows.
C2_LOWER_ORIGIN_EV = 0.0753


@dataclass(frozen=True)
class Transition:
    species: str
    system: str
    vp: int
    vl: int
    wavelength_vac_A: float
    lower_energy_eV: float
    gf: float
    j_lower: float
    label: str
    source: str
    source_line: int

    @property
    def loggf(self) -> float:
        return math.log10(self.gf)


def zip_text(path: Path, contains: str):
    archive = zipfile.ZipFile(path)
    name = next(n for n in archive.namelist() if contains in n)
    return io.TextIOWrapper(archive.open(name), encoding="utf-8-sig", errors="replace")


def parse_c2():
    path = PRIMARY / "c2_brooke2013/BrookeEtAl-C2-2013-JQSRT.zip"
    with zip_text(path, "C2SwanLineList") as stream:
        for line_no, raw in enumerate(stream, 1):
            p = raw.split()
            if len(p) < 16:
                continue
            try:
                vp, vl, jp, jl = int(p[0]), int(p[1]), float(p[2]), float(p[3])
                calculated, elow, f = float(p[9]), float(p[12]), float(p[14])
            except ValueError:
                continue
            wn = calculated
            if wn > 0 and f > 0:
                yield Transition("C2", "Swan", vp, vl, 1e8 / wn,
                                 elow * HC_EV_CM + C2_LOWER_ORIGIN_EV,
                                 f * (2 * jl + 1), jl,
                                 p[15], str(path.relative_to(ROOT)), line_no)


def parse_cn():
    path = PRIMARY / "cn_brooke2014/table4.dat.gz"
    with gzip.open(path, "rt", errors="replace") as stream:
        for line_no, raw in enumerate(stream, 1):
            try:
                upper, lower = raw[0], raw[2]
                vp, vl = int(raw[4:6]), int(raw[7:9])
                jl = float(raw[16:21])
                calculated = float(raw[50:60])
                elow, f = float(raw[70:80]), float(raw[94:106])
            except ValueError:
                continue
            wn = calculated
            if wn > 0 and f > 0:
                yield Transition("CN", f"{upper}-{lower}", vp, vl, 1e8 / wn,
                                 elow * HC_EV_CM, f * (2 * jl + 1), jl,
                                 raw[107:118].strip(), str(path.relative_to(ROOT)), line_no)


def parse_ch():
    path = PRIMARY / "ch_masseron2014/table14.dat.gz"
    with gzip.open(path, "rt", errors="replace") as stream:
        for line_no, raw in enumerate(stream, 1):
            try:
                wave_air = float(raw[5:18])
                gf, elow = float(raw[20:34]), float(raw[38:49])
                vl, jl, vp = int(raw[51]), float(raw[54:58]), int(raw[82])
            except ValueError:
                continue
            if raw[112:116].strip() != "12CH" or gf <= 0:
                continue
            # Morton/IAU standard dry-air conversion; source wavelengths are air.
            s2 = (1e4 / wave_air) ** 2
            n = 1 + 1e-8 * (8342.13 + 2406030 / (130 - s2) + 15997 / (38.9 - s2))
            yield Transition("CH", f"{raw[117]}-{raw[119]}", vp, vl, wave_air * n,
                             elow * HC_EV_CM, gf, jl, raw[122:132].strip(),
                             str(path.relative_to(ROOT)), line_no)


def parse_brooke_xx(species: str, relpath: str, member: str):
    path = PRIMARY / relpath
    with zip_text(path, member) as stream:
        for line_no, raw in enumerate(stream, 1):
            p = raw.split()
            if len(p) < 16:
                continue
            try:
                vp, vl, jl = int(p[0]), int(p[1]), float(p[3])
                # NH has N'/N'' fields; OH does not.
                shift = 2 if species == "NH" else 0
                calculated = float(p[9 + shift])
                elow, f = float(p[11 + shift]), float(p[13 + shift])
            except ValueError:
                continue
            wn = calculated
            if wn > 0 and f > 0:
                yield Transition(species, "X-X", vp, vl, 1e8 / wn,
                                 elow * HC_EV_CM, f * (2 * jl + 1), jl,
                                 p[14 + shift], str(path.relative_to(ROOT)), line_no)


def inventory():
    yield from parse_c2()
    yield from parse_ch()
    yield from parse_cn()
    yield from parse_brooke_xx("NH", "nh_brooke2015/BrookeEtAl-NH-2015-JCP.zip", "NH-XX-Linelist")
    yield from parse_brooke_xx("OH", "oh_brooke2016/OH-Supplementary.zip", "OH-XX-Line_list")


def main() -> None:
    transitions = list(inventory())
    index = defaultdict(list)
    for tr in transitions:
        index[(tr.species, tr.system)].append(tr)
    for key in index:
        index[key].sort(key=lambda tr: 1e8 / tr.wavelength_vac_A)

    rows = []
    for target in csv.DictReader(AMARSI.open()):
        if target["species"] == "12C16O":
            continue
        wave = float(target["wavelength_vac_nm"]) * 10
        energy, loggf = float(target["lower_energy_eV"]), float(target["published_loggf"])
        key0 = (target["species"], target["system"])
        target_wn = 1e8 / wave
        vp, vl = (int(x) for x in target["band"].strip("()").split("-"))
        wn_tolerance = 2.0 if target["species"] == "C2" else 0.30
        source = index.get(key0, ())
        source_wn = [1e8 / tr.wavelength_vac_A for tr in source]
        # Same sorted tolerance-range mechanism as pipeline.line_match: wavelength
        # supplies candidates but can never decide identity.
        lo = bisect.bisect_left(source_wn, target_wn - wn_tolerance)
        hi = bisect.bisect_right(source_wn, target_wn + wn_tolerance)
        candidates = source[lo:hi]
        # Molecular papers commonly publish observed/calculated positions in cm-1.
        # Match in that native coordinate: a fixed Angstrom tolerance becomes
        # physically nonsensical over Amarsi's 0.4--15 micron span.
        nearby = [x for x in candidates if x.vp == vp and x.vl == vl]
        # Published excitation energies are rounded to 0.001 eV and can use
        # slightly different molecular term origins between source releases.
        physical = [x for x in nearby if abs(x.lower_energy_eV - energy) <= 0.005]
        exact = [x for x in physical if abs(x.loggf - loggf) <= 0.006]
        status = "UNMATCHED"
        matches = exact
        summed_loggf = ""
        subset_candidate_count = 0
        if len(exact) == 1:
            status = "PRIMARY_TUPLE_MATCH"
        elif len(exact) > 1:
            status = "AMBIGUOUS_COMPONENT_MATCH"
        elif physical:
            subsets = []
            for size in range(2, min(5, len(physical)) + 1):
                for group in itertools.combinations(physical, size):
                    total = sum(x.gf for x in group)
                    summed = math.log10(total)
                    if abs(summed - loggf) <= 0.006:
                        centroid = sum((1e8 / x.wavelength_vac_A) * x.gf for x in group) / total
                        subsets.append((abs(summed-loggf), abs(centroid-target_wn), group, summed))
            subset_candidate_count = len(subsets)
            if len(subsets) > 1:
                # 🔴 RYA-1144. This branch used to resolve the tie with
                # min(subsets, key=...) -- an argmin over candidate IDENTITIES, which is
                # the "nearest wins" defect RYA-1037 forbids, wearing a third hat
                # (|dgf|, then centroid distance, then subset size). 26 of 32 rows had
                # between 2 and 16 subsets that all reproduce the published loggf, and
                # every one of them was counted as matched coverage. Finding A
                # combination is not finding THE combination: refuse.
                _, _, group, summed = min(subsets, key=lambda x: (x[0], x[1], len(x[2])))
                summed_loggf = f"{summed:.6f}"
                status, matches = "AMBIGUOUS_SUM_MATCH", list(group)
            elif subsets:
                _, _, group, summed = subsets[0]
                summed_loggf = f"{summed:.6f}"
                status, matches = "PRIMARY_UNRESOLVED_SUM_MATCH", list(group)
            else:
                total = sum(x.gf for x in physical)
                summed = math.log10(total)
                summed_loggf = f"{summed:.6f}"
                status, matches = "STRENGTH_MISMATCH", physical
        elif nearby:
            status, matches = "ENERGY_MISMATCH", nearby

        rows.append({
            "source_row": target["source_row"], "species": target["species"],
            "system": target["system"], "band": target["band"],
            "wavelength_vac_nm": target["wavelength_vac_nm"],
            "lower_energy_eV": target["lower_energy_eV"],
            "published_loggf": target["published_loggf"], "join_status": status,
            "component_count": len(matches), "summed_loggf": summed_loggf,
            "subset_candidate_count": subset_candidate_count,
            "primary_wavelengths_A": ";".join(f"{x.wavelength_vac_A:.5f}" for x in matches[:8]),
            "primary_energies_eV": ";".join(f"{x.lower_energy_eV:.6f}" for x in matches[:8]),
            "primary_loggfs": ";".join(f"{x.loggf:.6f}" for x in matches[:8]),
            "transition_labels": ";".join(x.label for x in matches[:8]),
            "primary_source": matches[0].source if matches else "",
            "primary_lines": ";".join(str(x.source_line) for x in matches[:8]),
        })

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary = Counter((r["species"], r["join_status"]) for r in rows)
    print(json.dumps({"primary_transitions": len(transitions),
                      "target_rows": len(rows),
                      "status": {f"{k[0]}:{k[1]}": v for k, v in summary.items()}}, indent=2))


if __name__ == "__main__":
    main()
