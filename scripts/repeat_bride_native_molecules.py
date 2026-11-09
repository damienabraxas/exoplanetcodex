"""Native all/off repeatability and measured Darwin child-program memory."""
import argparse,json,os,platform,re,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--build',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if platform.system()!='Darwin':raise RuntimeError('Darwin native RSS units required')
    from config.constants import require_codex_path
    from pipeline.bride_reference import fingerprint
    from scripts.run_bride_native_molecules import native_status
    run=json.loads((a.source/'report.json').read_text());build=json.loads((a.build/'report.json').read_text())
    plan_path=ROOT/'config/experiments/rya1281_native_molecules_repeat.json';plan=json.loads(plan_path.read_text())
    if fingerprint(a.build/'report.json')!=run['build']:raise ValueError('Changed native build')
    entry=next(x for x in build['builds']if x['target']=='bsyn_lu');program=a.build/'bsyn_lu/bsyn_lu'
    if fingerprint(program)!=entry['program']:raise ValueError('Changed native program')
    a.output.mkdir(parents=True,exist_ok=False);result=dict(schema='codex.bride_native_molecule_repeat/1',scientific_admission='HOLD',generator=fingerprint(__file__),
        source=fingerprint(a.source/'report.json'),plan=fingerprint(plan_path),program=fingerprint(program),scope=plan['scope'],cases=[])
    for name in plan['cases']:
        old=a.source/name;entry=next(x for x in run['cases']if x['name']==name)
        if fingerprint(old/'input.txt')!=entry['input']or fingerprint(old/'native_flux.txt')!=entry['native_flux']:raise ValueError('Changed repeat source')
        directory=a.output/name;directory.mkdir();(directory/'DATA').symlink_to(require_codex_path('engines.ispec')/'synthesizer/turbospectrum/DATA',target_is_directory=True)
        original=(old/'input.txt').read_text();old_output=str(old/'native_flux.txt');new_output=directory/'native_flux.txt'
        if original.count(old_output)!=1:raise ValueError('Unknown native result-file layout')
        text=original.replace(old_output,str(new_output));(directory/'input.txt').write_text(text)
        before=time.perf_counter();command=['/usr/bin/time','-l',str(program)]
        r=subprocess.run(command,input=text,cwd=directory,text=True,capture_output=True,timeout=plan['seconds_per_case'],env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
        elapsed=time.perf_counter()-before;(directory/'stdout.txt').write_text(r.stdout);(directory/'stderr.txt').write_text(r.stderr);native_status(r.returncode,r.stdout)
        matches=re.findall(r'^\s*(\d+)\s+maximum resident set size\s*$',r.stderr,re.M)
        if len(matches)!=1 or int(matches[0])<=0:raise ValueError('Missing Darwin native RSS measurement')
        if new_output.read_bytes()!=(old/'native_flux.txt').read_bytes():raise ValueError('Repeated native numeric spectrum differs')
        result['cases'].append(dict(name=name,seconds=elapsed,native_maximum_resident_set_size_bytes=int(matches[0]),
            raw_spectrum_byte_identical=True,only_result_path_changed=True,command=command,input=fingerprint(directory/'input.txt'),
            native_flux=fingerprint(new_output),stdout=fingerprint(directory/'stdout.txt'),stderr=fingerprint(directory/'stderr.txt')))
    (a.output/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result['cases']))


if __name__=='__main__':main()
