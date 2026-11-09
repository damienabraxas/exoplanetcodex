"""Losslessly repack native outputs with standard tar hardlinks; never edit inputs."""
import argparse
import hashlib
from pathlib import Path
import tarfile
import copy
import io


def repack(source,destination):
    if destination.exists() or source.resolve()==destination.resolve():
        raise ValueError('Choose a new destination; source evidence is immutable')
    seen={}
    with tarfile.open(source) as original, tarfile.open(destination,'w:gz') as compact:
        for member in original.getmembers():
            if member.isfile():
                data=original.extractfile(member).read()
                digest=hashlib.sha256(data).hexdigest()
                if digest in seen:
                    member.type=tarfile.LNKTYPE
                    member.linkname=seen[digest]
                    member.size=0
                    compact.addfile(member)
                else:
                    seen[digest]=member.name
                    compact.addfile(member,io.BytesIO(data))
            else:
                compact.addfile(member)
    # Compare every native path byte-for-byte, including all aliased artifacts.
    with tarfile.open(source) as a,tarfile.open(destination) as b:
        before={m.name:a.extractfile(m).read() for m in a.getmembers() if m.isfile() or m.islnk()}
        after={m.name:b.extractfile(m).read() for m in b.getmembers() if m.isfile() or m.islnk()}
    if before!=after:
        raise ValueError('Repack did not preserve all native artifact bytes')
    print(f'{len(before)} native paths preserved byte-for-byte; {destination.stat().st_size} bytes')


def split_cases(source,prefix):
    """Partition native case/state tables; validate shared receipts after joining."""
    with tarfile.open(source) as archive:
        members=[m for m in archive.getmembers() if m.isfile() or m.islnk()]
        files={m.name:archive.extractfile(m).read() for m in members}
    restored={}
    for case in ('DC','Breit_state','Breit_multipoles'):
        destination=Path(str(prefix)+f'.{case}.tar.gz')
        if destination.exists():
            raise ValueError('Choose a new partition destination')
        seen={}
        with tarfile.open(destination,'w:gz') as output:
            for m in members:
                group=('Breit_state' if '/Breit_static_limit/' in m.name else
                       'Breit_multipoles' if '/Breit_static_limit_multipoles/' in m.name else 'DC')
                if group!=case:
                    continue
                member=copy.copy(m)
                data=files[m.name]
                digest=hashlib.sha256(data).hexdigest()
                if digest in seen:
                    member.type=tarfile.LNKTYPE
                    member.linkname=seen[digest]
                    member.size=0
                    output.addfile(member)
                else:
                    seen[digest]=m.name
                    member.type=tarfile.REGTYPE
                    member.linkname=''
                    member.size=len(data)
                    output.addfile(member,io.BytesIO(data))
        with tarfile.open(destination) as check:
            restored.update({m.name:check.extractfile(m).read() for m in check.getmembers() if m.isfile() or m.islnk()})
        print(destination,destination.stat().st_size)
    if restored!=files:
        raise ValueError('Partition lost native artifacts')
    print(f'All {len(files)} native paths preserved byte-for-byte across the three parts')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('source',type=Path)
    ap.add_argument('destination',type=Path)
    ap.add_argument('--split-cases',action='store_true')
    args=ap.parse_args()
    if args.split_cases:
        split_cases(args.source,args.destination)
    else:
        repack(args.source,args.destination)
