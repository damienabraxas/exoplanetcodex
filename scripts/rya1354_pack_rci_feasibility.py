"""Pack v48-v51 RCI-feasibility evidence on Sirius. Files >= 4 MiB stay off-archive, declared
with the receipt's sha256/bytes and the archived stdin that regenerates them."""
import json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from rya1354_pack_native_stream import pack

ROOT = Path('/mnt/codex-data/codex/rya1354')
CAP = 4 * 1024**2


def producer(rel):
    p = Path(rel); d = p.parent; n = p.name
    if d.name in ('full', 'mr'):
        return d / 'rcsfgenerate.stdin'
    if d.name in ('DC', 'DCB'):
        return (d.parent / 'full' / 'rcsfgenerate.stdin') if n == 'rcsf.inp' else d / 'rcsfinteract.stdin'
    if d.name == 'run':
        if n.startswith('mcp.'):
            return d / 'rangular.stdin'
        if n in ('rcsf.inp', 'state.c'):
            return d.parent / 'DCB' / 'rcsfinteract.stdin'
        if n in ('rwfn.inp', 'state.w'):
            return d / 'rwfnestimate.stdin'
        if n == 'rci.res':
            return d / 'rci.stdin'
    return None   # logs and other files without a regenerating input are stored (fragmented), never dropped


def main(version, label, out, name='rci_feasibility'):
    run = f'runs/{name}_{version}'
    receipt = json.loads((ROOT / run / 'receipt.json').read_text())
    off = {}
    for rel, rec in receipt['files'].items():
        stdin = producer(rel) if rec['bytes'] >= CAP else None
        if stdin is not None:
            if stdin.as_posix() not in receipt['files']:
                raise ValueError('Regeneration input not in receipt: ' + str(stdin))
            off[f'{run}/{rel}'] = dict(sha256=rec['sha256'], bytes=rec['bytes'], regenerate_stdin=f'{run}/{stdin.as_posix()}')
    pack(ROOT, [Path(f'evidence/{name}_driver_{version}'), Path(run)], Path(out), off_archive_regenerable=off or None)
    print(version, label, len(off), sum(v['bytes'] for v in off.values()))


if __name__ == '__main__':
    main(*sys.argv[1:5])
