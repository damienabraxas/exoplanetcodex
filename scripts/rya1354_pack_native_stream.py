"""Single-pass lossless v2 capsules; bounded raw chunks and gzip level1."""
import argparse
import hashlib
import io
import json
from pathlib import Path
import tarfile
from rya1354_native_parts import load_parts

CHUNK=4*1024*1024
LIMIT=5*1024*1024


def pack(root,paths,out,chunk=CHUNK,*,references=None,repo_root=None,receipt_inventories=None,off_archive_regenerable=None):
    """references: {path: {sha256,bytes,archive,archive_path}} for files byte-identical to an
    earlier committed archive (schema v3); they are not stored again but restored and hash-checked."""
    if not 1024<=chunk<=CHUNK:raise ValueError('Raw chunk must be within proven archive ceiling')
    root=root.resolve();out.mkdir(parents=True,exist_ok=False)
    selected={}
    for path in paths:
        path=root/path
        if not path.exists():raise FileNotFoundError(path)
        for p in sorted(path.rglob('*')) if path.is_dir() else [path]:
            if p.is_file():selected[str(p.relative_to(root))]=p
    for key,record in (off_archive_regenerable or {}).items():   # declared off-archive: verified on disk, not stored
        if key in selected:
            h=hashlib.sha256()
            with selected[key].open('rb') as f:
                for block in iter(lambda:f.read(1<<24),b''):h.update(block)
            if (h.hexdigest(),selected[key].stat().st_size)!=(record['sha256'],record['bytes']):
                raise ValueError('Declared off-archive file differs on disk: '+key)
            del selected[key]
    if not selected:raise ValueError('No native files selected')
    references=dict(references or {})
    if set(references)&set(selected):raise ValueError('A path is both stored and referenced')
    parts={};files={};stored={};fragments={};pending={};pending_size=0
    def flush():
        nonlocal pending,pending_size
        if not pending:return
        path=out/f'part{len(parts)+1:03}.tar.gz'
        with tarfile.open(path,'w:gz',compresslevel=1) as archive:
            for name,data in pending.items():
                info=tarfile.TarInfo(name);info.size=len(data);info.mode=0o644;info.mtime=0
                archive.addfile(info,io.BytesIO(data))
        data=path.read_bytes()
        if len(data)>LIMIT:raise ValueError('Compressed archive exceeded ceiling')
        parts[path.name]=dict(sha256=hashlib.sha256(data).hexdigest(),bytes=len(data))
        pending={};pending_size=0
    def add(name,data):
        nonlocal pending_size
        if pending_size+len(data)>chunk:flush()
        if name in stored:raise ValueError('Duplicate stored member')
        stored[name]=dict(sha256=hashlib.sha256(data).hexdigest(),bytes=len(data));pending[name]=data;pending_size+=len(data)
    for name,p in sorted(selected.items()):
        size=p.stat().st_size;digest=hashlib.sha256();offset=0
        with p.open('rb') as src:
            if size<=chunk:
                data=src.read();digest.update(data);add(name,data);offset=len(data)
            else:
                fragments[name]=[];key=hashlib.sha256(name.encode()).hexdigest()
                while data:=src.read(chunk):
                    member=f'.rya1354_fragments/{key}/{len(fragments[name]):06}'
                    digest.update(data);add(member,data)
                    fragments[name].append(dict(member=member,offset=offset,bytes=len(data),sha256=hashlib.sha256(data).hexdigest()));offset+=len(data)
        if offset!=size or p.stat().st_size!=size:raise ValueError('Native file changed during packing')
        files[name]=dict(sha256=digest.hexdigest(),bytes=size)
    flush()
    manifest=dict(schema='rya1354.native-evidence-parts.v2',parts=parts,files=files,fragments=fragments,stored_members=stored,native_paths_preserved_byte_for_byte=len(files),omitted_angular_caches=[],external_restart_caches={},verification='Single pass, <=4MiB raw groups, gzip level1; every original path and byte retained. Existing independent v2 decoder verifies all hashes, reconstruction and original receipts.')
    if references:
        archives={}
        for name,r in references.items():
            if set(r)!={'sha256','bytes','archive','archive_path'}:raise ValueError('Reference record differs: '+name)
            if r['archive']=='.':
                if r['archive_path'] not in files:raise ValueError('Self-reference target is not stored here: '+name)
                archives['.']='.'
            else:archives[r['archive']]=hashlib.sha256((Path(repo_root)/r['archive']/'manifest.json').read_bytes()).hexdigest()
        manifest.update(schema='rya1354.native-evidence-parts.v3',referenced_archives=archives,referenced_files=references,
            native_paths_preserved_byte_for_byte=len(files)+len(references),
            verification='As v2 for stored files. Referenced files are byte-identical to the named, manifest-pinned earlier archive; the same decoder restores them from it and checks every hash before receipts are verified.')
    if receipt_inventories or off_archive_regenerable:
        manifest.update(receipt_inventories=dict(receipt_inventories or {}),off_archive_regenerable=dict(off_archive_regenerable or {}))
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    # Do not declare completion until the existing independent decoder passes.
    replay,_=load_parts(out,repo_root=repo_root)
    if set(replay)!=set(files)|set(references):raise ValueError('Reconstructed native file set differs')
    print(json.dumps(dict(parts=len(parts),native_files=len(files),referenced_files=len(references),compressed_bytes=sum(v['bytes'] for v in parts.values()))))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('paths',nargs='+',type=Path);a=ap.parse_args();pack(a.root,a.paths,a.out)
