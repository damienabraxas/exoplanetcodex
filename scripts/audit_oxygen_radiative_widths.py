"""Audit held O radiative widths against its own A graph; no atom correction or RT."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))

def main():
    from config.constants import require_codex_path
    from pipeline.multi_atom import read_multi_atom,to_lightweaver
    from pipeline.reference_engines import fingerprint
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--width-comparison',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    atom_file=require_codex_path('experiments.bride_oxygen_atom')
    reference=require_codex_path('experiments.bride_oxygen_reference')
    expected=json.loads(reference.read_text())['atom']
    if fingerprint(atom_file)!=expected:raise ValueError('Changed held atom')
    atom=read_multi_atom(atom_file);model=to_lightweaver(atom)
    if atom.symbol!='O' or len(atom.lines)!=256 or len(atom.levels)!=41:raise ValueError('Unexpected source scope')
    width_file=a.width_comparison/'widths.npz'
    width_report=a.width_comparison/'report.json'
    wr=json.loads(width_report.read_text())
    if fingerprint(width_file)!=wr['output'] or wr['inputs'].get(str(atom_file))!=expected:raise ValueError('Wrong width oracle/atom')
    checks={str(x):fingerprint(x)for x in [atom_file,reference,width_file,width_report,ROOT/'pipeline/multi_atom.py']}
    with np.load(width_file)as d:widths=d['native_widths'].copy()
    if widths.shape!=(256,56,5):raise ValueError('Unexpected oracle dimensions')
    natural=np.array([x.gRad for x in atom.lines])
    np.testing.assert_array_equal(widths[:,:,1],np.broadcast_to(natural[:,None],(256,56)))
    decay=np.zeros(41)
    for line in model.lines:decay[line.j]+=line.Aji
    rows=[]
    for index,(source,line)in enumerate(zip(atom.lines,model.lines)):
        partial=decay[line.i]+decay[line.j]
        # Internal graph diagnostic only: omitted channels may increase real decay.
        replacement=max(source.gRad,partial)
        fraction=(replacement-source.gRad)/widths[index,:,0]
        rows.append(dict(index=index,i_1based=source.i,j_1based=source.j,
            vacuum_A=1e8/(atom.levels[source.j-1].E_cm-atom.levels[source.i-1].E_cm),
            source_gamma_s_inv=source.gRad,own_Aji_s_inv=float(line.Aji),
            partial_graph_gamma_s_inv=float(partial),
            source_below_own_Aji=bool(source.gRad<line.Aji*(1-1e-4)),
            VDW_gamma_range_s_inv=[float(widths[index,:,2].min()),float(widths[index,:,2].max())],
            internal_graph_floor_total_width_fraction_range=[float(fraction.min()),float(fraction.max())]))
    for path,fp in checks.items():
        if fingerprint(path)!=fp:raise ValueError('Changed audit source')
    report=dict(generator=fingerprint(__file__),inputs=checks,scientific_admission='HOLD',
        method='Same held atom Aji graph; lower+upper outgoing spontaneous sums are incomplete radiative lifetime bounds. No observed fit, atom modification, spectrum/SE response or physical replacement inferred.',
        lifetime_reference='https://physics.nist.gov/Pubs/AtSpec/node18.html',
        own_Aji_comparison_relative_tolerance=1e-4,
        below_own_Aji_indices=[r['index']for r in rows if r['source_below_own_Aji']],records=rows)
    a.output.mkdir(parents=True,exist_ok=False)
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps([r for r in rows if r['source_below_own_Aji']],indent=2))

if __name__=='__main__':main()
