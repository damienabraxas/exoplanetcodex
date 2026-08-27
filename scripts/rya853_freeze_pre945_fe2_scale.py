#!/usr/bin/env python3
"""
RYA-853 scope 3 — freeze the Fe II gf scale AS IT WAS BEFORE WE ADOPTED THE REFEREE.
====================================================================================
RYA-945 (`b545dc6`) ingested Den Hartog 2019 Table 6 into `canonical_gf`. That was the
right thing to do for the line list, and it destroyed the scope-3 experiment: the referee
now sits on both sides of the comparison, so `ours - DH` is 0.000 by construction and the
verdict printed "LEGITIMATE - our scale IS the pure-lab scale ... solar-fitting REFUTED"
off a self-match.

The question scope 3 asks is about the scale that UNDERWROTE THE IONIZATION BALANCE, and
that scale exists only in git history now. This script lifts it out and commits it as a
frozen artifact, so the referee has a fixed input that no future ingest can quietly change.

    source: {PRE945_SHA} = b545dc6^, the commit immediately before the RYA-945 ingest

What it is NOT: a claim that the pre-945 values are better. They are simply the ones being
refereed. Post-945 the 22 overlap lines ARE Den Hartog's, by adoption, and asking whether
they agree with Den Hartog is not a question.

Usage:
    python3 scripts/rya853_freeze_pre945_fe2_scale.py
"""
from __future__ import annotations

import hashlib
import io
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

#: The commit immediately before [RYA-945] ingested DH19 into canonical_gf (b545dc6^).
PRE945_SHA = "a18b9faf44c074ace8cb1ee5c2765adfe251b764"
INGEST_SHA = "b545dc6c98f3bb21a86095cae73575a3cfbcb97b"
CANONICAL_REL = "data/linelists/canonical_gf.csv"

OUT_CSV = ROOT / "data" / "reference" / "fe_gf_lab" / "fe2_pre945_scale_snapshot.csv"
OUT_PROV = OUT_CSV.with_suffix(".prov.json")

KEEP = ["species", "wavelength_air_A", "excitation_potential_eV", "log_gf",
        "loggf_reference", "nist_grade", "gf_tier", "adjudication_status"]

#: The referee's own name, in every spelling it appears under in the reference column.
#: A snapshot that contains any of these is not a snapshot of a scale to be refereed.
REFEREE_MARKERS = ("denhartog2019", "den hartog 2019", "dh19", "2019apjs..243")


def _git_show(sha: str, rel: str) -> str:
    out = subprocess.run(["git", "show", f"{sha}:{rel}"], cwd=ROOT,
                         capture_output=True, text=True)
    if out.returncode != 0:
        raise SystemExit(
            f"cannot read {rel} at {sha[:9]}: {out.stderr.strip()}\n"
            "This needs the full history; a shallow clone will not do.")
    return out.stdout


def main() -> None:
    raw = _git_show(PRE945_SHA, CANONICAL_REL)
    d = pd.read_csv(io.StringIO(raw), comment="#", low_memory=False)
    fe2 = d[d.species.astype(str) == "Fe II"].copy()
    fe2 = fe2[[c for c in KEEP if c in fe2.columns]]
    fe2 = fe2.sort_values("wavelength_air_A").reset_index(drop=True)

    ref = fe2.loggf_reference.astype(str).str.lower()
    contaminated = fe2[ref.str.contains("|".join(REFEREE_MARKERS), regex=True, na=False)]

    # CONTROL, not a formality. The entire point of freezing at b545dc6^ is that the
    # referee is absent from it. If the SHA is ever edited to something after the ingest
    # this fires, instead of the referee silently grading its own homework again.
    if len(contaminated):
        raise SystemExit(
            f"REFUSING to write the snapshot: {len(contaminated)} rows at {PRE945_SHA[:9]} "
            f"already cite the referee (Den Hartog 2019). That is not a pre-adoption "
            f"scale.\n{contaminated.head(10).to_string()}")

    post = pd.read_csv(ROOT / CANONICAL_REL, comment="#", low_memory=False)
    post = post[post.species.astype(str) == "Fe II"]
    n_post_dh = int(post.loggf_reference.astype(str).str.lower()
                    .str.contains("|".join(REFEREE_MARKERS), regex=True, na=False).sum())

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    fe2.to_csv(OUT_CSV, index=False)
    sha256 = hashlib.sha256(OUT_CSV.read_bytes()).hexdigest()

    OUT_PROV.write_text(json.dumps({
        "ticket": "RYA-853 scope 3",
        "artifact": OUT_CSV.name,
        "what": ("the Fe II log gf scale as it stood BEFORE [RYA-945] ingested Den Hartog "
                 "2019 into canonical_gf — the scale that underwrote the solar Fe I - Fe II "
                 "ionization balance, and the only one the DH19 referee can test"),
        "source_file": CANONICAL_REL,
        "source_sha": PRE945_SHA,
        "source_sha_is": "b545dc6^ — the commit immediately before the RYA-945 ingest",
        "superseded_by_sha": INGEST_SHA,
        "why_frozen": ("after RYA-945 the overlap lines in canonical_gf ARE Den Hartog's, so "
                       "ours-minus-DH is 0.000 by construction and the referee scores its "
                       "own values. Re-running scope 3 against live canonical_gf on "
                       "2026-08-27 printed 'LEGITIMATE ... solar-fitting REFUTED' off a "
                       "zero-width comparison"),
        "n_rows": int(len(fe2)),
        "wavelength_span_A": [float(fe2.wavelength_air_A.min()),
                              float(fe2.wavelength_air_A.max())],
        "controls": [
            f"REFEREE-ABSENT: 0 of {len(fe2)} snapshot rows cite Den Hartog 2019 PASS",
            f"REFEREE-PRESENT-TODAY: {n_post_dh} Fe II rows in live canonical_gf DO cite it "
            f"— which is exactly why the snapshot is needed",
            f"UNIQUE: {int(fe2.wavelength_air_A.duplicated().sum())} duplicate wavelengths",
        ],
        "sha256": sha256,
        "regenerate": "python3 scripts/rya853_freeze_pre945_fe2_scale.py",
    }, indent=2))

    print(f"[freeze] Fe II rows at {PRE945_SHA[:9]} (b545dc6^): {len(fe2)}")
    print(f"[freeze] cite the referee: 0 (control PASS)")
    print(f"[freeze] cite the referee in LIVE canonical_gf today: {n_post_dh}")
    print(f"[out] {OUT_CSV}")
    print(f"[out] {OUT_PROV}")


if __name__ == "__main__":
    main()
