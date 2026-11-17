"""Plan a v3 evidence pack: which native paths are byte-identical to a committed archive
(external reference), which duplicate another path of the same pack (self reference),
and which must be stored. Input is `sha256sum` output from the native host."""
import argparse
import json
from pathlib import Path


def plan(sha_lines,archives,manifests):
    """archives: repository-relative archive dirs, searched in order; manifests: their manifest dicts."""
    by={}
    for archive,manifest in zip(archives,manifests):
        for path,record in sorted(manifest['files'].items()):
            by.setdefault(record['sha256'],(archive,path,record['bytes']))
    refs={};stored=[];first={}
    for line in sha_lines:
        if not line.strip():continue
        digest,path=line.split(None,1);path=path.strip()
        if path in refs or path in stored:raise ValueError('Duplicate native path: '+path)
        if digest in by:
            archive,source,size=by[digest]
            refs[path]=dict(sha256=digest,bytes=size,archive=archive,archive_path=source)
        else:stored.append(path);first.setdefault(digest,[]).append(path)
    for digest,paths in first.items():
        keep=sorted(paths)[0]
        for path in paths:
            if path==keep:continue
            stored.remove(path);refs[path]=dict(sha256=digest,bytes=None,archive='.',archive_path=keep)
    return refs,sorted(stored)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--sha',type=Path,required=True);ap.add_argument('--repo-root',type=Path,required=True)
    ap.add_argument('--archive',action='append',default=[]);ap.add_argument('--refs-out',type=Path,required=True)
    ap.add_argument('--stored-out',type=Path,required=True);a=ap.parse_args()
    manifests=[json.loads((a.repo_root/x/'manifest.json').read_text()) for x in a.archive]
    refs,stored=plan(a.sha.read_text().splitlines(),a.archive,manifests)
    a.refs_out.write_text(json.dumps(refs,indent=1)+'\n');a.stored_out.write_text('\n'.join(stored)+'\n')
    print(json.dumps(dict(external=sum(r['archive']!='.' for r in refs.values()),self=sum(r['archive']=='.' for r in refs.values()),stored=len(stored))))


if __name__=='__main__':main()
