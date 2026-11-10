"""Algebraic original-DAMP coefficient export; one source evaluation, no regression fit."""
import argparse,json,subprocess,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))

def main():
    from pipeline.reference_engines import fingerprint
    from pipeline.multi_atom import read_multi_atom
    from config.constants import require_codex_path
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--base',type=Path,required=True);p.add_argument('--laws',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();inputs={}
    def report(folder):
        path=folder/'report.json';inputs[str(path)]=fingerprint(path);r=json.loads(path.read_text());return r
    base=report(a.base);laws=report(a.laws);bf=a.base/'widths.npz';lf=a.laws/'laws.npz'
    for path,expected in [(bf,base['output']),(lf,laws['output'])]:
        if fingerprint(path)!=expected:raise ValueError('Changed source width output')
        inputs[str(path)]=expected
    atom_file=require_codex_path('experiments.bride_oxygen_atom');atom=read_multi_atom(atom_file);atom_receipt=fingerprint(atom_file)
    if base['inputs'][str(atom_file)]!=atom_receipt or atom.symbol!='O':raise ValueError('Changed source atom')
    with np.load(bf)as d:gamma=d['native_widths'].copy();T=d['temperature_K'].copy();HI=d['native_HI_m3'].copy();gw=d['source_gvdw'].copy()
    if gamma.shape!=(256,56,5)or any(x.gStark!=0 for x in atom.lines):raise ValueError('Unsupported source law support')
    # This is the original routine's default-REAL .3 literal, promoted to DP.
    beta=np.where(gw>=20,(1-(gw-np.floor(gw)))/2,float(np.float32(.3)))
    coefficient=gamma[:,0,2]/(HI[0]*T[0]**beta);bound=float(128*np.finfo(np.float64).eps)
    predicted=coefficient[:,None]*HI[None,:]*T[None,:]**beta[:,None];error=float(np.max(abs(predicted/gamma[:,:,2]-1)))
    if error>bound:raise ValueError('Source coefficient not constant across all56 states')
    checks=[];low=float(T.min());high=float(T.max())
    with np.load(lf)as d:
        for name,c in laws['cases'].items():
            temp=T+c.get('delta_T',0);density=HI*c.get('h_scale',1);actual=d[name][:,:,2];test=coefficient[:,None]*density[None,:]*temp[None,:]**beta[:,None];err=float(np.max(abs(test/actual-1)))
            if err>bound:raise ValueError('Source coefficient fails native independent control')
            low=min(low,float(temp.min()));high=max(high,float(temp.max()));checks.append(dict(case=name,max_relative_error=err))
    records=[]
    for i,line in enumerate(atom.lines):
        if line.gVdw!=1 and line.gVdw<20:raise ValueError('Unsupported source VDW family')
        records.append(dict(index=i,i_1based=line.i,j_1based=line.j,source_gvdw=line.gVdw,natural_gamma_s_inv=line.gRad,coefficient_m3_s_inv_K_minus_beta=float(coefficient[i]),beta=float(beta[i]),family='ABO'if line.gVdw>=20 else'legacy-Unsold'))
    packet=dict(schema='codex.original_multi_vdw_law/1',purpose='experimental-fixed-O-width-comparison',scientific_admission='HOLD',source_family='original-MULTI2.3-DAMP-O41',source_pin='4d6aa3c39d634b38da772e469a5692435eba1354',producer_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),atom=atom_receipt,source_atomic_mass_amu=atom.weight,temperature_numeric_interval_K=[low,high],coefficient_rule='Gamma_vdw_s^-1 = coefficient_m3_s^-1_K^-beta * neutral_H_total_m^-3 * T_K^beta',extra_helium_policy='Source prescription owns treatment; no added LW He multiplier. Historical Unsold header states fixedHe/H0.1; runtime density is total neutral H.',records=records,source_inputs=inputs,normalization_sample_index=0,normalization_method='Algebraic normalization of native sample0, no average/regression or observed spectrum fit',baseline_max_relative_error=error,control_checks=checks,roundoff_bound=bound,limits=['Numerical source-law interval only, not a physical stellar applicability hull. No RT, collision, EOS, observed solar, 3D or all-band admission.','No H level populations implied; density role is the original O DAMP sumNH1:5 aggregate.','Gamma law independent of microturbulence; Doppler normalization remains the caller responsibility.'])
    a.output.mkdir(parents=True,exist_ok=False);(a.output/'law.json').write_text(json.dumps(packet,indent=2)+'\n');(a.output/'report.json').write_text(json.dumps(dict(generator=fingerprint(__file__),inputs=inputs,output=fingerprint(a.output/'law.json'),baseline_max_relative_error=error,control_checks=checks,scientific_admission='HOLD'),indent=2)+'\n');print(json.dumps(dict(error=error,control_checks=checks,temperature_numeric_interval_K=[low,high])))
if __name__=='__main__':main()
