"""
scripts/decouple_gold_rya814.py — RYA-814
=========================================
Explode the monolithic gold reference into per-element frozen records, re-assemble
it, and PROVE the round trip is byte-identical.

    python3 scripts/decouple_gold_rya814.py --roundtrip          # prove it (no writes)
    python3 scripts/decouple_gold_rya814.py --explode --apply    # write the records
    python3 scripts/decouple_gold_rya814.py --assemble           # print the view
    python3 scripts/decouple_gold_rya814.py --status             # complete vs partial

THE ROUND-TRIP GATE IS THE POINT
--------------------------------
RYA-814 is a REPRESENTATION refactor: no frozen value may move. sha256 equality
between the original monolith and the re-assembled one is what makes that claim
unfalsifiable-if-wrong. If the hashes differ, the decomposition is wrong and the
work stops there — it does not get "reconciled".

--roundtrip runs entirely IN MEMORY and writes nothing, so the proof can be run
against the live frozen gold at any time without touching it.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pipeline.element_freeze import (  # noqa: E402
    ElementFreezeError, ELEMENTS_DIR, assemble, assembly_status, freeze_element,
    preamble_path, read_element, serialise_fields)

SOLAR_DIR = ROOT / "data" / "reference" / "solar"


def _sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def _current_version() -> str:
    return (SOLAR_DIR / "CURRENT").read_text().strip()


def _gold_path(version: str) -> Path:
    return SOLAR_DIR / f"solar_abundances_{version}.csv"


def _git_commit() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "--short=12", "HEAD"], cwd=ROOT,
                              capture_output=True, text=True).stdout.strip()
    except Exception:
        return ""


def parse_monolith(version: str):
    """(preamble, columns, [(element, fields, verbatim_line)]) from a frozen file."""
    raw = _gold_path(version).read_text(encoding="utf-8")
    trailing_newline = raw.endswith("\n")
    lines = raw.split("\n")
    if trailing_newline:
        lines = lines[:-1]

    comment_lines = [l for l in lines if l.startswith("#")]
    body = [l for l in lines if l and not l.startswith("#")]
    columns = next(csv.reader([body[0]]))
    rows = body[1:]

    parsed = []
    for i, line in enumerate(rows):
        vals = next(csv.reader([line]))
        if len(vals) != len(columns):
            raise SystemExit(f"row {i} has {len(vals)} fields, header has {len(columns)}")
        fields = dict(zip(columns, vals))
        fields["_row_index"] = i          # row ORDER is part of the bytes
        parsed.append((vals[0], fields, line))

    preamble = {
        "comment_lines": comment_lines,
        "columns": columns,
        "trailing_newline": trailing_newline,
        "exploded_from": version,
    }
    return preamble, columns, parsed


def cmd_roundtrip(version: str) -> int:
    """Explode -> re-serialise -> re-assemble IN MEMORY and compare sha256."""
    path = _gold_path(version)
    original = path.read_bytes()
    original_sha = _sha256_bytes(original)
    preamble, columns, parsed = parse_monolith(version)

    print(f"round-trip gate — gold {version}")
    print(f"  source            {path}")
    print(f"  original sha256   {original_sha}")
    print(f"  rows              {len(parsed)}   columns {len(columns)}")

    # (1) every parsed row must re-serialise to its own source line
    drift = [el for el, fields, line in parsed
             if serialise_fields(fields, columns) != line]
    if drift:
        print(f"\n  FIELD DRIFT on {len(drift)} row(s): {drift[:6]}")
        print("  Parsed fields do not reproduce the source line. The decomposition "
              "is lossy — STOP (RYA-814: no value may move).")
        return 1
    print(f"  field re-serialisation  OK ({len(parsed)}/{len(parsed)} rows)")

    # (2) reassemble from the pieces and compare bytes
    buf = io.StringIO()
    buf.write("\n".join(preamble["comment_lines"]) + "\n")
    buf.write(",".join(columns) + "\n")
    for _el, _f, line in sorted(parsed, key=lambda t: t[1]["_row_index"]):
        buf.write(line + "\n")
    rebuilt = buf.getvalue()
    if not preamble["trailing_newline"]:
        rebuilt = rebuilt.rstrip("\n")
    rebuilt_sha = _sha256_bytes(rebuilt.encode("utf-8"))

    print(f"  reassembled sha256 {rebuilt_sha}")
    same = rebuilt_sha == original_sha
    print(f"\n  {'PASS — byte-identical' if same else 'FAIL — bytes differ'}")
    if not same:
        ob, rb = original.decode("utf-8"), rebuilt
        for i, (a, b) in enumerate(zip(ob.split("\n"), rb.split("\n"))):
            if a != b:
                print(f"    first difference at line {i}:")
                print(f"      original : {a[:120]}")
                print(f"      rebuilt  : {b[:120]}")
                break
        print(f"    length {len(ob)} vs {len(rb)}")
        return 1
    return 0


def cmd_explode(version: str, apply: bool) -> int:
    preamble, columns, parsed = parse_monolith(version)
    commit = _git_commit()
    print(f"explode gold {version} -> per-element records "
          f"({'APPLY' if apply else 'dry-run'})")
    if not apply:
        for el, _f, line in parsed[:4]:
            print(f"  {el:3s}  {line[:88]}")
        print(f"  ... {len(parsed)} elements total")
        print("\n  [dry-run] re-run with --apply to write.")
        return 0

    ELEMENTS_DIR.mkdir(parents=True, exist_ok=True)
    import json
    preamble_path().write_text(
        json.dumps(preamble, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    written = 0
    for el, fields, line in parsed:
        existing = read_element(el)
        if existing is not None and existing.version == version:
            print(f"  skip {el} (already frozen at {version})")
            continue
        try:
            freeze_element(el, fields, line, version=version, source_commit=commit,
                           provenance=f"exploded from solar_abundances_{version}.csv "
                                      f"(RYA-814 decoupling; representation only, no "
                                      f"value moved)",
                           columns=columns, allow_new_version=True)
            written += 1
        except ElementFreezeError as exc:
            print(f"  REFUSED {el}: {exc}")
            return 1
    print(f"  wrote {written} record(s) -> {ELEMENTS_DIR}")
    return 0


def cmd_assemble(version: str) -> int:
    text = assemble(version_label=version)
    sha = _sha256_bytes(text.encode("utf-8"))
    original_sha = _sha256_bytes(_gold_path(version).read_bytes())
    sys.stdout.write(text)
    print(f"\n# assembled sha256 {sha}", file=sys.stderr)
    print(f"# {version} sha256   {original_sha}", file=sys.stderr)
    print(f"# {'MATCH' if sha == original_sha else 'DIFFER'}", file=sys.stderr)
    return 0 if sha == original_sha else 1


def cmd_status() -> int:
    st = assembly_status()
    print(f"  {st['label']}")
    print(f"  frozen: {', '.join(st['frozen']) if st['frozen'] else '(none)'}")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--version", default=None, help="gold version (default: CURRENT)")
    ap.add_argument("--roundtrip", action="store_true")
    ap.add_argument("--explode", action="store_true")
    ap.add_argument("--assemble", action="store_true")
    ap.add_argument("--status", action="store_true")
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args()
    v = a.version or _current_version()
    if a.roundtrip:
        sys.exit(cmd_roundtrip(v))
    if a.explode:
        sys.exit(cmd_explode(v, a.apply))
    if a.assemble:
        sys.exit(cmd_assemble(v))
    if a.status:
        sys.exit(cmd_status())
    ap.error("pick --roundtrip / --explode / --assemble / --status")
