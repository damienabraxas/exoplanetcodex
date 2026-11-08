"""Isolated Sirius-only O I solver diagnostic, never an admitted gf product.

Uses a point static nucleus (Z=8), single configurations and default solver
numerics. No measured level substitution, experimental gf, Solar abundance,
isotope claim or uncertainty estimate. Every prompt, binary and output is hashed.
This is an exploratory diagnostic, not a blinded or converged atomic benchmark.
"""
from __future__ import annotations

import argparse
from decimal import Decimal
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def modern_coupling_layout(text: str) -> str:
    """Export book configuration/eigenvector records to native ATSP2K columns.

    Sources: native cfgin2.f FORMAT 7, gstate.f FORMAT 36, eigvec.f MCHF
    branch. Add explicit closed/total orbital inventory, open-orbital list,
    64-column CSF fields and coefficient field; preserve every scientific value.
    """
    lines = text.splitlines(keepends=True)
    opened = []
    for i in range(2, len(lines)-1, 2):
        if lines[i].lstrip().startswith("*"):
            break
        n = lines[i].count("(")
        if not n:
            raise ValueError("expected a book CSF orbital record")
        orbital_record = lines[i].rstrip("\r\n")
        orbitals = [orbital_record[8*j+1:8*j+4] for j in range(n)]
        for orbital in orbitals:
            if orbital not in opened:
                opened.append(orbital)
        coefficient = Decimal(orbital_record[40:].strip())
        encoded_coefficient = f"{coefficient:15.12f}"
        if Decimal(encoded_coefficient) != coefficient:
            raise ValueError("native field cannot preserve coefficient exactly")
        lines[i] = orbital_record[:8*n].ljust(64) + encoded_coefficient + "\n"
        line = lines[i+1].rstrip("\r\n")
        labels = [line[8*j+5:8*j+8] for j in range(2*n-1)]
        if any(len(label) != 3 or label == "   " for label in labels):
            raise ValueError("missing book coupling labels; no inference allowed")
        converted = "".join(" " + label for label in labels)
        if [converted[4*j+1:4*j+4] for j in range(len(labels))] != labels:
            raise ValueError("coupling layout conversion changed a label")
        lines[i+1] = converted + "\n"
    nclosed = len(lines[1].strip().split())
    lines[0] = lines[0].rstrip("\r\n").ljust(30) + f"{nclosed:3d}{nclosed+len(opened):4d}\n"
    lines.insert(2, "".join(" " + orbital for orbital in opened) + "\n")
    return "".join(lines)


def run(executable: Path, stdin: str, cwd: Path, timeout: int = 180) -> dict:
    if not executable.is_file():
        raise FileNotFoundError(executable)
    name = executable.name
    (cwd / f"{name}.stdin").write_text(stdin)
    with (cwd / f"{name}.stdout").open("xb") as out, (cwd / f"{name}.stderr").open("xb") as err:
        p = subprocess.run([str(executable)], input=stdin.encode(), cwd=cwd,
                           stdout=out, stderr=err, timeout=timeout)
    record = {"executable": str(executable), "sha256": sha(executable), "returncode": p.returncode}
    # Both legacy packages can STOP with exit 0 on rejected scientific input.
    # These observed refusals must survive into the receipt and CLI exit status.
    output = ((cwd / f"{name}.stdout").read_text(errors="replace") +
              (cwd / f"{name}.stderr").read_text(errors="replace")).upper()
    record["semantic_errors"] = [marker for marker in
        ("TOO MANY ELECTRONS FOR THIS CODE", "CONFIGURATION DATA WRONG", "PROGRAM HALTED",
         "MAXIMUM ITERATIONS IN SCF EXCEEDED")
        if marker in output]
    if 'THE ORDER OF' in output and 'IS INCONSISTENT' in output:
        record['semantic_errors'].append('INCONSISTENT_ORBITAL_ORDER')
    if p.returncode:
        raise RuntimeError(f"{name}: exit {p.returncode}; retained logs in {cwd}")
    return record


def grasp_case(bindir: Path, out: Path, orbital: str, two_j: str, state_input: str) -> list:
    out.mkdir()
    records = []
    steps = [
        ("rnucleus", "8\n0\n0\n0\n0\n0\n"),
        ("rcsfgenerate", f"*\n0\n1s(2,i)2s(2,i)2p(3,i){orbital}(1,i)\n\n3s,3p\n{two_j}\n0\nn\n"),
    ]
    for name, stdin in steps:
        records.append(run(bindir / name, stdin, out))
    (out / "rcsf.inp").write_bytes((out / "rcsf.out").read_bytes())
    for name, stdin in [("rangular", "y\n"), ("rwfnestimate", "y\n2\n*\n"),
                        ("rmcdhf", f"y\n{state_input}*\n*\n100\n")]:
        records.append(run(bindir / name, stdin, out))
    for required in ("rmcdhf.sum", "rmix.out", "rwfn.out"):
        if not (out / required).stat().st_size:
            raise ValueError(f"empty solver artifact: {out / required}")
    for source, target in (("rcsf.inp", "state.c"), ("rmix.out", "state.m"), ("rwfn.out", "state.w")):
        (out / target).write_bytes((out / source).read_bytes())
    # Upstream unique-label postprocessing allocates 25.6 GB even for this tiny
    # case on the pinned revision. Ordinary LS compositions avoid that path;
    # retain the failed unique-label attempt separately, never patch the solver.
    records.append(run(bindir / "jj2lsj", "state\nn\nn\ny\n", out))
    return records


def atsp_case(bindir: Path, out: Path, orbital: str, term: str) -> list:
    out.mkdir()
    records = [run(bindir / "Gencl", f"\nO_{term}\n 1s  2s\n2p(3){orbital}(1)\n\n2p,{orbital}\n0\n\n{term}\n\n", out),
               run(bindir / "Nonh", "n\ny\n", out),
               run(bindir / "HF", f"O,{term},8.\n  1s  2s\n2p(3){orbital}(1)\nall\ny\ny\nn\nn\n", out)]
    (out / "wfn.inp").write_bytes((out / "wfn.out").read_bytes())
    with (out / "cfg.inp").open("a") as fh:
        fh.write("*\n")  # upstream MCHF example's explicit orthogonality terminator
    records.append(run(bindir / "Mchf", f"O,{term},8.\n=4\ny\ny\ny\n", out))
    return records


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--grasp-bin", type=Path, required=True)
    ap.add_argument("--atsp-bin", type=Path, required=True)
    ap.add_argument("--atsp2k-bin", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    if platform.node().lower().split(".")[0] != "sirius":
        raise RuntimeError("Atomic solver computation must run on Sirius")
    args.out.mkdir(parents=True, exist_ok=False)
    report = {"schema": "rya1354.solver-diagnostic.v1", "started_utc": datetime.now(timezone.utc).isoformat(),
              "host": platform.node(), "python": sys.version, "script_sha256": sha(Path(__file__)),
              "scope": "EXPLORATORY_SINGLE_CONFIGURATION_NOT_BLIND_NOT_CONVERGENCE_VALIDATED",
              "status": "IN_PROGRESS", "steps": [], "scientific_admission": "HOLD"}
    report_file = args.out / "receipt.json"
    try:
        # Hund-high-spin candidate blocks only. LS assignment remains unverified;
        # no numerical output is promoted to a named physical transition.
        report["steps"] += grasp_case(args.grasp_bin, args.out / "grasp_3s_J2", "3s", "4,4", "1\n")
        report["steps"] += grasp_case(args.grasp_bin, args.out / "grasp_3p_J123", "3p", "2,6", "1\n1\n1\n5\n")
        transition = args.out / "grasp_E1"
        transition.mkdir()
        for stem, directory in (("lower", "grasp_3s_J2"), ("upper", "grasp_3p_J123")):
            for ext in ("c", "m", "w"):
                (transition / f"{stem}.{ext}").write_bytes((args.out / directory / f"state.{ext}").read_bytes())
        (transition / "isodata").write_bytes((args.out / "grasp_3s_J2/isodata").read_bytes())
        report["steps"].append(run(args.grasp_bin / "rbiotransform", "y\nn\nlower\nupper\ny\n", transition))
        report["steps"].append(run(args.grasp_bin / "rtransition", "y\nn\nlower\nupper\nE1\n", transition))
        # Independent nonrelativistic method sanity check on O I ground 3P.
        # It does not compare unlike ground and excited-state energies or claim
        # that a ground-state smoke calculation benchmarks excited-state gf.
        atsp = args.out / "atsp_ground_3P"
        atsp.mkdir()
        report["steps"].append(run(args.atsp_bin / "HF", "O,3P,8.\n  1s  2s\n2p(4)\nall\ny\ny\nn\nn\n", atsp))
        report["steps"] += atsp_case(args.atsp_bin, args.out / "atsp_3s_5S", "3s", "5S")
        report["steps"] += atsp_case(args.atsp_bin, args.out / "atsp_3p_5P", "3p", "5P")
        independent = args.out / "atsp_E1"
        independent.mkdir()
        for stem, directory in (("lower", "atsp_3s_5S"), ("upper", "atsp_3p_5P")):
            for src, ext in (("cfg.out", "c"), ("wfn.out", "w")):
                (independent / f"{stem}.{ext}").write_bytes((args.out / directory / src).read_bytes())
        report["steps"].append(run(args.atsp_bin / "Mltpol", "lower\nupper\nn\n\ne1\n*\n", independent))
        report["steps"].append(run(args.atsp_bin / "Lstr", "lower\nupper\nn\ny\n", independent))
        if args.atsp2k_bin:
            modern = args.out / "atsp2k_biorthogonal_E1"
            modern.mkdir()
            for stem in ("lower", "upper"):
                for ext in ("c", "w"):
                    source = independent / f"{stem}.{ext}"
                    if ext == "c":
                        (modern / f"{stem}.{ext}").write_text(modern_coupling_layout(source.read_text()))
                    else:
                        (modern / f"{stem}.{ext}").write_bytes(source.read_bytes())
            report["steps"].append(run(args.atsp2k_bin / "biotr", "lower\nupper\nn\nn\nm\nE1\n*\n", modern))
        report["semantic_failures"] = [r for r in report["steps"] if r["semantic_errors"]]
        report["status"] = ("COMPLETED_WITH_SCIENTIFIC_FAILURES" if report["semantic_failures"]
                            else "SOLVER_STEPS_COMPLETED_REQUIRES_NUMERICAL_AND_IDENTITY_REVIEW")
    except (OSError, RuntimeError, ValueError, subprocess.TimeoutExpired) as exc:
        report["status"] = "FAILED_RETAINED_FOR_DIAGNOSIS"
        report["error"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        report["files"] = {str(p.relative_to(args.out)): {"sha256": sha(p), "bytes": p.stat().st_size}
                           for p in sorted(args.out.rglob("*")) if p.is_file() and p != report_file}
        report_file.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(report["status"])
    return 1 if report.get("semantic_failures") else 0


if __name__ == "__main__":
    raise SystemExit(main())
