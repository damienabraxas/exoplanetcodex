"""Replay the pinned upstream native Si E2 control in a new retained directory.

This is software interoperability evidence, not an O I experimental benchmark.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--atsp", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    if platform.node().lower().split(".")[0] != "sirius":
        raise RuntimeError("Run atomic solvers on Sirius")
    args.out.mkdir(parents=True, exist_ok=False)
    inputs = {}
    for target, relative in {
        "I.c": "examples/files_c/1De1.4.c", "F.c": "examples/files_c/1Se1.4.c",
        "I.w": "examples/e1/e1.14_4_1D.w", "F.w": "examples/e1/e1.14_4_1S.w",
        "I.l": "examples/e1/e1.14_4_1D.l", "F.l": "examples/e1/e1.14_4_1S.l",
        "expected.ls": "examples/tr_ls/1De1.1Se1.14_4.ls",
    }.items():
        data = (args.atsp / relative).read_bytes()
        (args.out / target).write_bytes(data)
        inputs[target] = {"source": relative, "sha256": hashlib.sha256(data).hexdigest()}
    stdin = "I\nF\nn\nn\nc\nE2\n*\n"
    (args.out / "stdin.txt").write_text(stdin)
    with (args.out / "stdout.txt").open("xb") as stdout, (args.out / "stderr.txt").open("xb") as stderr:
        result = subprocess.run([str(args.atsp / "bin/biotr")], input=stdin.encode(), cwd=args.out,
                                stdout=stdout, stderr=stderr, timeout=180)
    text = (args.out / "stdout.txt").read_text(errors="replace")
    native = args.out / "I.F.ls"
    same = native.is_file() and native.read_bytes() == (args.out / "expected.ls").read_bytes()
    record = {"schema": "rya1354.atsp-native-control.v1", "inputs": inputs,
              "scope": "UPSTREAM_SI_E2_SOFTWARE_CONTROL_NOT_OXYGEN_VALIDATION",
              "process_exit": result.returncode,
              "semantic_failure": "CONFIGURATION DATA WRONG" in text or "PROGRAM HALTED" in text,
              "native_equals_expected_bytes": same,
              "binary_sha256": hashlib.sha256((args.atsp / "bin/biotr").read_bytes()).hexdigest()}
    (args.out / "receipt.json").write_text(json.dumps(record, indent=2) + "\n")
    return 0 if result.returncode == 0 and same and not record["semantic_failure"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
