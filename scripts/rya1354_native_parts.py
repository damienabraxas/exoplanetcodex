"""Read v1/v2 native evidence parts, restoring fragmented files without extraction."""
import hashlib
import json
from pathlib import Path

from rya1354_two_method_summary import capsule_files, verify_artifacts


def checked(data,record,label):
    if len(data)!=record['bytes'] or hashlib.sha256(data).hexdigest()!=record['sha256']:
        raise ValueError('Native bytes/hash differ: '+label)


def reassemble(stored,fragments):
    result=dict(stored);used=set()
    for name,records in fragments.items():
        if name in result or name.startswith('.rya1354_fragments/') or not records:
            raise ValueError('Invalid or colliding fragmented original')
        chunks=[];offset=0
        for record in records:
            member=record['member']
            if member in used or not member.startswith('.rya1354_fragments/') or record['offset']!=offset:
                raise ValueError('Duplicate, unordered or noncontiguous fragment')
            if member not in result or record['bytes']<=0:
                raise ValueError('Missing/empty native fragment')
            data=result.pop(member);checked(data,record,member)
            used.add(member);chunks.append(data);offset+=len(data)
        result[name]=b''.join(chunks)
    if any(name.startswith('.rya1354_fragments/') for name in result):
        raise ValueError('Unclaimed native fragment')
    return result


V3='rya1354.native-evidence-parts.v3'
_REFERENCE_CACHE={}


SELF='.'


def referenced(manifest,repo_root,own):
    """v3 only: restore files byte-identical to an EARLIER committed archive (pinned by its
    manifest hash, decoded by this same reader) or, with archive '.', to a file STORED in this one."""
    external={a for a in manifest['referenced_archives'] if a!=SELF}
    if external and repo_root is None:raise ValueError('Referenced native files need the repository root')
    if manifest['referenced_archives'].get(SELF,SELF)!=SELF:raise ValueError('Self-reference pin differs')
    out={}
    for archive in sorted(external):
        pin=manifest['referenced_archives'][archive]
        path=(Path(repo_root)/archive).resolve()
        if Path(archive).is_absolute() or '..' in Path(archive).parts:raise ValueError('Referenced archive path is not repository-relative')
        if hashlib.sha256((path/'manifest.json').read_bytes()).hexdigest()!=pin:raise ValueError('Referenced archive manifest differs: '+archive)
        if path not in _REFERENCE_CACHE:_REFERENCE_CACHE[path]=load_parts(path,repo_root=repo_root)[0]
    for name,record in manifest['referenced_files'].items():
        if record['archive'] not in manifest['referenced_archives']:raise ValueError('Undeclared referenced archive')
        source=own if record['archive']==SELF else _REFERENCE_CACHE[(Path(repo_root)/record['archive']).resolve()]
        if record['archive_path'] not in source:raise ValueError('Referenced original missing: '+name)
        data=source[record['archive_path']];checked(data,record,name);out[name]=data
    return out


def load_parts(directory,*,repo_root=None):
    manifest=json.loads((directory/'manifest.json').read_text())
    schema=manifest.get('schema')
    if schema not in ('rya1354.native-evidence-parts.v1','rya1354.native-evidence-parts.v2',V3) or not manifest['parts']:
        raise ValueError('Unsupported/empty native bundle')
    stored={}
    for name,record in manifest['parts'].items():
        if Path(name).name!=name:
            raise ValueError('Part filename is not local')
        path=directory/name;checked(path.read_bytes(),record,name)
        for key,data in capsule_files(path,verify=False).items():
            if key in stored and stored[key]!=data:
                raise ValueError('Conflicting native member')
            stored[key]=data
    expected=manifest['files'] if schema=='rya1354.native-evidence-parts.v1' else manifest['stored_members']
    if set(stored)!=set(expected):
        raise ValueError('Missing/extra stored member')
    for key,data in stored.items():checked(data,expected[key],key)
    files=reassemble(stored,manifest['fragments']) if schema!='rya1354.native-evidence-parts.v1' else stored
    if set(files)!=set(manifest['files']):
        raise ValueError('Missing/extra reconstructed native file')
    for key,data in files.items():checked(data,manifest['files'][key],key)
    if schema==V3:
        extra=referenced(manifest,repo_root,files)
        if set(extra)&set(files) or len(files)+len(extra)!=manifest['native_paths_preserved_byte_for_byte']:
            raise ValueError('Referenced/stored native path accounting differs')
        files.update(extra)
    inventoried=sidecar_inventories(manifest,files)
    verify_artifacts(files,external_restart_caches=manifest.get('external_restart_caches'),inventoried=inventoried)
    return files,manifest


def sidecar_inventories(manifest,files):
    """Optional `receipt_inventories` {receipt: sidecar} for a native receipt written without a files
    list. The sidecar must be archived; every file it names under the receipt's directory must be
    archived with that hash or listed in `off_archive_regenerable` with the same hash and bytes, and
    every archived file under that directory must be named. Returns the receipts so covered."""
    declared=manifest.get('receipt_inventories') or {}
    off=manifest.get('off_archive_regenerable') or {}
    covered=set()
    for receipt,sidecar in declared.items():
        if receipt not in files or sidecar not in files:raise ValueError('Declared receipt inventory is not archived: '+receipt)
        if 'files' in json.loads(files[receipt]):raise ValueError('Receipt already carries a files list: '+receipt)
        parent=receipt.rsplit('/',1)[0]+'/';inv=json.loads(files[sidecar])
        if str(inv.get('run'))+'/'!=parent:raise ValueError('Sidecar names another run: '+receipt)
        named={parent+k:v for k,v in inv['files'].items()}
        for key,record in named.items():
            if key in files:checked(files[key],record,key)
            elif key not in off or {k:off[key][k] for k in ('sha256','bytes')}!=record:
                raise ValueError('Inventoried native file missing or changed: '+key)
        if {k for k in files if k.startswith(parent)}-set(named):raise ValueError('Archived file absent from sidecar: '+receipt)
        covered|=set(named)
    if set(off)-covered:raise ValueError('Off-archive file not covered by a declared sidecar')
    for key,record in off.items():
        if record.get('regenerate_stdin') not in files:raise ValueError('Off-archive file lacks archived regeneration input: '+key)
    return set(declared)
