import hashlib
import json
import shlex
import shutil
import subprocess
from pathlib import Path

root = Path("/mnt/codex-data/codex/rya1354")
grasp = root / "toolchains/grasp"
out = root / "runs/mcdhf_biotransform_v81"
out.mkdir()
shutil.copyfile(__file__, out / "driver.py")
status = "FAILED_RETAINED"
stage = "build"
try:
    native = grasp / "build/src/appl/rbiotransform90"
    cmake = native / "CMakeFiles/rbiotransform.dir"
    original = grasp / "src/appl/rbiotransform90/mcpin.f90"
    changed = out / original.name
    source = original.read_text()
    assert source.count("EVEC(1:500)") == 1
    changed.write_text(source.replace("EVEC(1:500)", "EVEC"))
    binary = out / "rbiotransform"
    obj = out / "mcpin.f90.o"
    flags_make = dict(
        line.split(" = ", 1) for line in (cmake / "flags.make").read_text().splitlines()
        if " = " in line
    )
    command = ["/usr/bin/gfortran", *shlex.split(flags_make["Fortran_FLAGS"]),
               *shlex.split(flags_make["Fortran_INCLUDES"]), "-I", str(native),
               "-c", str(changed), "-o", str(obj)]
    subprocess.run(command, cwd=out, check=True)
    link = shlex.split((cmake / "link.txt").read_text())
    old_obj = "CMakeFiles/rbiotransform.dir/mcpin.f90.o"
    assert link.count(old_obj) == 1
    link[link.index(old_obj)] = str(obj)
    link[link.index("-o") + 1] = str(binary)
    link = [arg for arg in link if not arg.startswith("-Wl,--dependency-file=")]
    subprocess.run(link, cwd=native, check=True)

    transition = grasp / "build/bin/rtransition"
    report = root / "runs/mcdhf_e1_v79/report.py"
    compare = root / "runs/mcdhf_e1_v79/compare.py"
    refs = root / "data/reference/li2023_oi/raw"
    a1, a2 = refs / "tablea1.tsv", refs / "tablea2.tsv"
    (out / "build_commands.json").write_text(json.dumps([command, link], indent=2) + "\n")
    provenance = [out / "driver.py", original, changed, binary, transition, report, compare, a1, a2]

    for layer, source_dir in (
        ("n5", root / "runs/mcdhf_layers_v77/n5_E1"),
        ("n6", root / "runs/mcdhf_layers_v78/n6_E1"),
    ):
        stage = f"{layer}_inputs"
        work = out / layer
        work.mkdir()
        for name in ("isodata", "lower.c", "lower.m", "lower.w", "lower.lsj.lbl",
                     "upper.c", "upper.m", "upper.w", "upper.lsj.lbl",
                     "rbiotransform.stdin", "rtransition.stdin"):
            dest = work / name
            shutil.copyfile(source_dir / name, dest)
            provenance.append(dest)
        (out / "launch_sha.txt").write_text("".join(
            f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path}\n" for path in provenance
        ))
        for exe, name, marker in (
            (binary, "rbiotransform", "BIOTRANSFORM: Execution complete."),
            (transition, "rtransition", "RTRANSITION: Execution complete."),
        ):
            stage = f"{layer}_{name}"
            log = work / f"{name}.stdout"
            with (work / f"{name}.stdin").open() as stdin, log.open("w") as stdout, (work / f"{name}.stderr").open("w") as stderr:
                subprocess.run([exe], stdin=stdin, stdout=stdout, stderr=stderr,
                               cwd=work, timeout=7200, check=True)
            assert marker in log.read_text(), f"{stage} completion marker missing"
        stage = f"{layer}_report"
        table = work / "lower.upper.t"
        assert table.stat().st_size > 0
        if layer == "n5":
            assert table.read_bytes() == (source_dir / table.name).read_bytes(), "n5 native E1 table differs"
        subprocess.run(["python3", report, work, a1, a2, work / "loggf_vs_li2023.csv"], check=True)
    stage = "comparison"
    with (out / "comparison.md").open("w") as stream:
        subprocess.run(["python3", compare, out / "n5/loggf_vs_li2023.csv",
                        out / "n6/loggf_vs_li2023.csv"], stdout=stream, check=True)
    status = "COMPLETED_REVIEW_REQUIRED"
finally:
    receipt = dict(status=status, stage=stage, scientific_admission="HOLD",
                   files=sorted(str(path.relative_to(out)) for path in out.rglob("*") if path.is_file()))
    (out / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
