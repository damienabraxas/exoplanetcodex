"""Decompose all-level physical-rate differences and explicit frequency controls."""
import argparse,csv,json,sys,re
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
CASES={'default':'rya1371-hop-native-physical-rates-v1','default_midpoint':'rya1371-hop-physical-midpoint1-v1','default_local8446':'rya1371-hop-physical-midpoint8446-v1','native_support_only':'rya1371-hop-physical-native-support-v1','native_grid':'rya1371-hop-physical-native-grid-v1','native_midpoint':'rya1371-hop-physical-native-midpoint-v1','native_local8446':'rya1371-hop-physical-native-midpoint8446-v1'}


def flow(r,c,n):
    off=r+c;incoming=(off*n[None,:,:]).sum(axis=1);outgoing=off.sum(axis=0)*n
    return incoming-outgoing,incoming+outgoing


def main():
    from pipeline.reference_engines import fingerprint
    from pipeline.multi_atom import read_multi_atom
    from config.constants import require_codex_path
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--raw-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();checks={};summary={};level_rows=[];pair_rows=[];maps=[];base=None
    atom_path=require_codex_path('experiments.bride_oxygen_atom');atom=read_multi_atom(atom_path);checks[str(atom_path)]=fingerprint(atom_path)
    for case,folder in CASES.items():
        directory=a.raw_root/folder;doc=json.loads((directory/'report.json').read_text())
        if not doc['computational_completed']or 'native_physical_rate_residual'not in doc:raise ValueError('Require qualified physical-rate replay')
        for name in ['report.json','rates.npz']:checks[str(directory/name)]=fingerprint(directory/name)
        with np.load(directory/'rates.npz')as d:arrays={key:d[key].copy()for key in ['native_R_s_inv','replay_R_s_inv','native_C_s_inv','replay_C_s_inv','populations','nstar','wave_vac_nm']}
        if base is None:base=arrays
        else:
            for key in ['native_R_s_inv','native_C_s_inv','replay_C_s_inv','populations','nstar']:np.testing.assert_array_equal(arrays[key],base[key])
        R=arrays['native_R_s_inv'];S=arrays['replay_R_s_inv'];C=arrays['native_C_s_inv'];D=arrays['replay_C_s_inv'];n=arrays['populations'];native_net,native_gross=flow(R,C,n);own_net,own_gross=flow(S,D,n)
        rad_net,_=flow(S-R,np.zeros_like(C),n);col_net,_=flow(np.zeros_like(R),D-C,n);scale=np.maximum(own_gross,np.finfo(float).tiny)
        if np.max(abs(own_net-native_net-rad_net-col_net)/scale)>256*np.finfo(float).eps:raise ValueError('Physical net-flow decomposition fails')
        relative=abs(own_net)/scale;delta=abs(own_net-native_net)/scale;k,z=np.unravel_index(np.argmax(relative),relative.shape);matrix=S-R;contributions=matrix[k,:,z]*n[:,z]-matrix[:,k,z]*n[k,z]
        partners=[]
        for partner in np.argsort(abs(contributions))[::-1][:6]:
            if partner==k or contributions[partner]==0:continue
            lo=min(k,int(partner));hi=max(k,int(partner));partners.append(dict(partner_1based=int(partner+1),kind='BB'if atom.level(lo+1).stage==atom.level(hi+1).stage else'BF',energy_gap_vac_A=1e8/(atom.level(hi+1).E_cm-atom.level(lo+1).E_cm),delta_net_m3_s=float(contributions[partner]),signed_fraction_of_own_gross=float(contributions[partner]/scale[k,z])))
        summary[case]=dict(points=len(arrays['wave_vac_nm']),seconds=doc['seconds'],max_physical_residual=float(relative[k,z]),native_physical_residual=float(np.max(abs(native_net)/native_gross)),max_radiative_delta_normalized=float(np.max(abs(rad_net)/scale)),max_collision_quantization_delta_normalized=float(np.max(abs(col_net)/scale)),worst=dict(level_1based=int(k+1),depth_top_to_bottom=int(z),partners=partners),component_7_13_top_up_native_s_inv=float(R[12,6,0]),component_7_13_top_up_replay_s_inv=float(S[12,6,0]),component_7_13_top_relative_up=float(S[12,6,0]/R[12,6,0]-1),scope='Frozen identical n/nStar/physical fields; numerical grid controls only. Residual is evaluated at native radiation-population stage, not a new coupled solve or solar validation.')
        maps.append(relative)
        for level in range(41):
            for depth in range(56):level_rows.append(dict(case=case,level_1based=level+1,depth_top_to_bottom=depth,residual=float(relative[level,depth]),radiative_delta_normalized=float(rad_net[level,depth]/scale[level,depth]),collision_delta_normalized=float(col_net[level,depth]/scale[level,depth])))
        for lo in range(41):
            for hi in range(lo+1,41):
                if not np.any(R[hi,lo])and not np.any(R[lo,hi]):continue
                f=matrix[hi,lo]*n[lo]-matrix[lo,hi]*n[hi];effect=np.maximum(abs(f)/scale[lo],abs(f)/scale[hi]);depth=int(np.argmax(effect))
                pair_rows.append(dict(case=case,lower_1based=lo+1,upper_1based=hi+1,kind='BB'if atom.level(lo+1).stage==atom.level(hi+1).stage else'BF',energy_gap_vac_A=1e8/(atom.level(hi+1).E_cm-atom.level(lo+1).E_cm),max_normalized_endpoint_effect=float(effect[depth]),depth_top_to_bottom=depth,signed_delta_into_upper_m3_s=float(f[depth])))
    params=require_codex_path('experiments.openmulti3d_oxygen_accurate_run')/'native/params.nml';source=require_codex_path('engines.openmulti3d_reference')/'experiments/OpenMULTI3D/SRC/rt_solvers_mod.f90'
    # The live AMR declaration, not a guessed default, determines interface spacing.
    amr=re.search(r"\bAMR\s*=\s*'([^']*)'",params.read_text())
    if amr is None:raise ValueError('Missing native interface convention')
    checks[str(params)]=fingerprint(params);checks[str(source)]=fingerprint(source)
    if "amr_boundary = amr/='' .and. self%dir < 0"not in source.read_text():raise ValueError('Changed source AMR interface convention')
    a.output.mkdir(parents=True,exist_ok=False)
    for name,rows in [('level_depth.csv',level_rows),('pairs.csv',pair_rows)]:
        with(a.output/name).open('w',newline='')as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,len(summary),figsize=(20,5),sharey=True)
    for ax,(case,stats),values in zip(axes,summary.items(),maps):
        im=ax.imshow(np.log10(np.maximum(values,1e-12)),origin='lower',aspect='auto',vmin=-8,vmax=-2);ax.set_title(case.replace('_','\n')+'\n'+f"max {stats['max_physical_residual']:.3g}");ax.set_xlabel('depth (top to bottom)')
    axes[0].set_ylabel('level index (zero-based)');fig.colorbar(im,ax=axes,label='log10 |physical net| / gross; display floor 1e-12',shrink=.8);fig.savefig(a.output/'residual_maps.png',dpi=150);plt.close(fig)
    for path,fp in checks.items():
        if fingerprint(path)!=fp:raise ValueError('Changed physical-rate audit input')
    report=dict(generator=fingerprint(__file__),inputs=checks,cases=summary,all_level_depth_rows=len(level_rows),all_pair_rows=len(pair_rows),native_amr=amr[1].strip(),scientific_admission='HOLD',computational_completed=True,interpretation='Coarse native-grid agreement is not numerical convergence or physical validation. Endpoint-only and spacing/support controls are distinct; both grid families require resolution checks. Native AMR is enabled, so the simple missing-upstream-interface-spacing hypothesis is not applicable.')
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({case:stats['max_physical_residual']for case,stats in summary.items()}))
if __name__=='__main__':main()
