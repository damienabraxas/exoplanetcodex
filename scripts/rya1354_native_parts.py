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


def load_parts(directory):
    manifest=json.loads((directory/'manifest.json').read_text())
    schema=manifest.get('schema')
    if schema not in ('rya1354.native-evidence-parts.v1','rya1354.native-evidence-parts.v2') or not manifest['parts']:
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
    expected=manifest['stored_members'] if schema.endswith('.v2') else manifest['files']
    if set(stored)!=set(expected):
        raise ValueError('Missing/extra stored member')
    for key,data in stored.items():checked(data,expected[key],key)
    files=reassemble(stored,manifest['fragments']) if schema.endswith('.v2') else stored
    if set(files)!=set(manifest['files']):
        raise ValueError('Missing/extra reconstructed native file')
    for key,data in files.items():checked(data,manifest['files'][key],key)
    verify_artifacts(files,external_restart_caches=manifest.get('external_restart_caches'))
    return files,manifest
