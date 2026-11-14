"""Build isolated read-only config diagnostic alongside existing raw-rate API."""
import difflib
import platform
from pathlib import Path
import shutil
import subprocess
from rya1354_fac_bootstrap import sha


def build(source,probe,evidence):
    if platform.node().lower().split('.')[0]!='sirius':raise RuntimeError('Atomic frontend builds belong on Sirius')
    probe.mkdir(parents=True,exist_ok=False);evidence.mkdir(parents=True,exist_ok=False)
    original=(source/'sfac/sfac.c').read_text();root=Path(__file__).parent
    snippet=(root/'rya1354_fac_raw_transition.c').read_text()+'\n'+(root/'rya1354_fac_config_diagnostics.c').read_text()
    marker='#include "global.h"\n';registration='  {"TransitionTable", PTransitionTable, METH_VARARGS},'
    if original.count(marker)!=1 or original.count(registration)!=1:raise ValueError('Probe insertion points differ')
    patched=original.replace(marker,marker+'\n'+snippet+'\n',1).replace(registration,
        '  {"RawTransition", PRawTransition, METH_VARARGS},\n  {"ConfigDiagnostics", PConfigDiagnostics, METH_VARARGS},\n'+registration,1)
    (probe/'sfac.c').write_text(patched);(probe/'upstream.sfac.c').write_text(original)
    (probe/'observability.diff').write_text(''.join(difflib.unified_diff(original.splitlines(True),patched.splitlines(True),fromfile='pinned/sfac.c',tofile='read-only-diagnostics/sfac.c')))
    for name in ('stoken.c','stoken.h'):shutil.copyfile(source/'sfac'/name,probe/name)
    shutil.copyfile(source/'LICENSE',probe/'LICENSE');shutil.copyfile(source/'sfac/Makefile',probe/'source.Makefile')
    before=sha(source/'libfac.a')
    with (probe/'make.stdout').open('xb') as out,(probe/'make.stderr').open('xb') as err:
        p=subprocess.run(['make','-f',str(source/'sfac/Makefile'),'sfac'],cwd=probe,stdout=out,stderr=err,timeout=180,check=False)
    if p.returncode or sha(source/'libfac.a')!=before:raise ValueError('Config diagnostic build failed/changed library')
    for p in probe.iterdir():
        if p.is_file() and p.name!='sfac':shutil.copyfile(p,evidence/p.name)
    binary=probe/'sfac'
    return dict(path=str(binary),sha256=sha(binary),bytes=binary.stat().st_size,libfac_sha256=before,upstream_binary_sha256=sha(source/'sfac/sfac'))
