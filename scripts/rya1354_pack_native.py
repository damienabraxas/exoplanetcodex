"""Lossless, size-bounded native evidence parts with a complete byte manifest.

Only explicitly listed task outputs are read. Native mcp.* caches stay on the
solver host; all other file paths/bytes are retained. No filesystem extraction.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import tarfile


def pack(root,paths,out,limit=5*1024*1024,*,omit_rci_restart_cache=False):
    root=root.resolve()
    out.mkdir(parents=True,exist_ok=False)
    files={}
    omitted=[]
    external={}
    for selected in paths:
        path=root/selected
        if not path.exists():
            raise FileNotFoundError(path)
        for p in sorted(path.rglob('*')) if path.is_dir() else [path]:
            if not p.is_file():
                continue
            name=str(p.relative_to(root))
            if p.name.startswith('mcp.'):
                omitted.append(name)
                continue
            if omit_rci_restart_cache and p.name=='rci.res':
                external[name]=dict(sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size,
                                    disposition='REGENERABLE_HAMILTONIAN_RESTART_CACHE_RETAINED_ON_SIRIUS',
                                    source='Pinned rci90/lodres.f90 and matrix.f90: restart matrix cache; wavefunctions/mixing/CSFs/transition tables retained in archive')
                continue
            files[name]=p.read_bytes()
    parts=[]
    def write(names):
        candidate=out/f'.candidate_{len(parts):03}.tar.gz'
        with tarfile.open(candidate,'w:gz') as archive:
            seen={}
            for name in names:
                data=files[name]
                digest=hashlib.sha256(data).hexdigest()
                info=tarfile.TarInfo(name)
                info.mode=0o644
                info.mtime=0
                if digest in seen:
                    info.type=tarfile.LNKTYPE
                    info.linkname=seen[digest]
                    archive.addfile(info)
                else:
                    info.size=len(data)
                    archive.addfile(info,io.BytesIO(data))
                    seen[digest]=name
        if candidate.stat().st_size>limit:
            # Remove only this newly created temporary representation; original
            # source outputs and already accepted parts are never changed.
            candidate.unlink()
            if len(names)<2:
                raise ValueError('Single native artifact exceeds part ceiling: '+names[0])
            middle=len(names)//2
            write(names[:middle])
            write(names[middle:])
        else:
            target=out/f'part{len(parts)+1:03}.tar.gz'
            candidate.rename(target)
            parts.append(target)
    write(sorted(files))
    replay={}
    for part in parts:
        with tarfile.open(part) as archive:
            replay.update({m.name:archive.extractfile(m).read() for m in archive.getmembers() if m.isfile() or m.islnk()})
    if replay!=files:
        raise ValueError('Native archive partition did not preserve every byte')
    manifest=dict(schema='rya1354.native-evidence-parts.v1',
                  parts={p.name:dict(sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size) for p in parts},
                  files={n:dict(sha256=hashlib.sha256(b).hexdigest(),bytes=len(b)) for n,b in files.items()},
                  native_paths_preserved_byte_for_byte=len(files),omitted_angular_caches=omitted,
                  external_restart_caches=external,
                  verification='Require every declared part, check all part/member hashes, then validate original receipts after joining. Original source outputs are untouched.')
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(parts=len(parts),native_files=len(files),compressed_bytes=sum(p.stat().st_size for p in parts))))


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--root',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('paths',nargs='+',type=Path)
    ap.add_argument('--omit-rci-restart-cache',action='store_true',help='Explicitly inventory/hash regenerable rci.res matrix caches retained on solver host')
    args=ap.parse_args()
    pack(args.root,args.paths,args.out,omit_rci_restart_cache=args.omit_rci_restart_cache)
