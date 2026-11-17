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


def referenced(manifest,repo_root):
    """v3 only: restore files that are byte-identical to an EARLIER committed archive.
    Each archive is decoded by this same reader and pinned by its manifest hash."""
    if repo_root is None:raise ValueError('Referenced native files need the repository root')
    out={}
    for archive,pin in manifest['referenced_archives'].items():
        path=(Path(repo_root)/archive).resolve()
        if Path(archive).is_absolute() or '..' in Path(archive).parts:raise ValueError('Referenced archive path is not repository-relative')
        if hashlib.sha256((path/'manifest.json').read_bytes()).hexdigest()!=pin:raise ValueError('Referenced archive manifest differs: '+archive)
        if path not in _REFERENCE_CACHE:_REFERENCE_CACHE[path]=load_parts(path,repo_root=repo_root)[0]
    for name,record in manifest['referenced_files'].items():
        if record['archive'] not in manifest['referenced_archives']:raise ValueError('Undeclared referenced archive')
        source=_REFERENCE_CACHE[(Path(repo_root)/record['archive']).resolve()]
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
        extra=referenced(manifest,repo_root)
        if set(extra)&set(files) or len(files)+len(extra)!=manifest['native_paths_preserved_byte_for_byte']:
            raise ValueError('Referenced/stored native path accounting differs')
        files.update(extra)
    verify_artifacts(files,external_restart_caches=manifest.get('external_restart_caches'))
    return files,manifest
