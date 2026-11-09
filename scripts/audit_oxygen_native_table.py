"""Measure consumed native atom formatting and inspect pinned MULTI damping semantics."""
import argparse,csv,json,subprocess,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

def main():
    from pipeline.reference_engines import fingerprint,load_engine
    from pipeline.multi_atom import read_multi_atom
    from pipeline.wavelength_util import vac_to_air
    from config.constants import require_codex_path
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--native-run',type=Path,required=True);p.add_argument('--multi-source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    pin='4d6aa3c39d634b38da772e469a5692435eba1354'
    if subprocess.check_output(['git','rev-parse','HEAD'],cwd=a.multi_source,text=True).strip()!=pin:raise ValueError('Changed MULTI source pin')
    _,app,_,_=load_engine('openmulti3d');atom_file=require_codex_path('experiments.bride_oxygen_atom');atom=read_multi_atom(atom_file);table=a.native_run/'native/atoms/O.txt';lines=table.read_text().splitlines();start=next(i for i,x in enumerate(lines)if x.split()[:2]==['i','j']and'lam_air'in x)+1;rows=[]
    for index,L in enumerate(atom.lines):
        fields=lines[start+index].split();idx,i,j=map(int,fields[:3]);lam=float(fields[3]);loggf=float(fields[5]);lower=atom.levels[L.i-1];upper=atom.levels[L.j-1]
        if (idx,i+1,j+1)!=(index,L.i,L.j):raise ValueError('Native table line identity differs')
        raw_vac=1e8/(upper.E_cm-lower.E_cm);raw_air=float(vac_to_air(raw_vac));raw_loggf=float(np.log10(L.f*lower.g));rows.append(dict(line_index=index,i_1based=L.i,j_1based=L.j,source_canonical_A=raw_air,native_table_coordinate_A=lam,native_table_minus_canonical_A=lam-raw_air,canonical_frame=('vacuum'if raw_vac<2000 else'air'),source_loggf=raw_loggf,native_loggf=loggf,native_minus_source_loggf=loggf-raw_loggf,native_f_over_source=10**loggf/lower.g/L.f))
    source=a.multi_source/'pymulti/files/multi/mul23/source/mul23_sub.f';text=source.read_text(encoding='latin1')
    if 'IF(GQ(KR).GE.0.0) THEN'not in text or 'GS=GQ(KR)*NE(K)'not in text or "IF(ATOMID.EQ.'H   ' .AND. GQ(KR).EQ.0.0) THEN"not in text:raise ValueError('Changed original MULTI damping branches')
    inputs={str(x):fingerprint(x)for x in [atom_file,table,source,app/'HOP/atom_convert_mod.f90',app/'HOP/bbtrans_mod.f90',app/'HOP/handyfuncs_mod.f90']}
    selected=[x for x in rows if any(abs(x['source_canonical_A']-c)<.7 for c in [6158.17,6300.304,7771.944,7774.166,7775.388,8446.4])]
    report=dict(generator=fingerprint(__file__),inputs=inputs,multi_source_pin=pin,scientific_admission='HOLD',nlines=len(rows),selected=selected,max_table_coordinate_difference_A=max(abs(x['native_table_minus_canonical_A'])for x in rows),max_absolute_loggf_shift=max(abs(x['native_minus_source_loggf'])for x in rows),conclusions=['Pinned original MULTI2.3 sets GS=GQ*ne when GQ>=0; automatic GQ0 branch applies to Hydrogen only. Oxygen gStark0 therefore means no Stark in this source, consistent with Bride and unlike native HOP automatic Cowley. This does not prove producer identity of held MULTI record.','Native dispatch converts and reads back table: writer emits AIR wavelength to0.001A and loggf to0.001dex. Consumed source atom is identical, but these rounded transfer parameters differ from Bride full precision.','Report source-determined shifts/strength changes, not a fitted wavelength or abundance adjustment. No source/atom/default adoption or exact line-profile parity asserted. Below2000 the canonical coordinate remains vacuum; external table uses its own air conversion. Table-coordinate differences there are not direct RT wavelength errors; external roundtrip must be qualified separately.'])
    a.output.mkdir(parents=True,exist_ok=False)
    with(a.output/'lines.csv').open('w')as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(dict(selected=selected,max_air_shift=report['max_table_coordinate_difference_A'])))
if __name__=='__main__':main()
