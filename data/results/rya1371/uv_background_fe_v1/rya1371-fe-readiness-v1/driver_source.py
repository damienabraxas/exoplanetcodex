"""Read-only full Fe source inventory and dry grid/resource handoff; no RT/SE."""
import argparse,json,sys
from collections import Counter
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def main():
    from pipeline.reference_engines import fingerprint,load_engine
    from pipeline.multi_atom import read_multi_atom,to_lightweaver
    from pipeline.hop_packet import MAX_PACKET_WAVELENGTHS
    from config.constants import codex_path
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    _,app,_,engine=load_engine('openmulti3d');source=app/'input/atoms/atom.fe607a';atom=read_multi_atom(source)
    checks={str(source):fingerprint(source),str(Path(__file__)):fingerprint(__file__),str(ROOT/'pipeline/multi_atom.py'):fingerprint(ROOT/'pipeline/multi_atom.py')}
    legacy=ROOT/'data/results/rya1055/atom_ion_reach.json';old=json.loads(legacy.read_text());checks[str(legacy)]=fingerprint(legacy)
    by_stage=Counter(x.stage for x in atom.levels);BB=Counter((atom.level(x.i).stage,atom.level(x.j).stage)for x in atom.lines);BF=Counter((atom.level(x.i).stage,atom.level(x.j).stage)for x in atom.continua)
    waves=np.array([1e7/(atom.level(x.j).E_cm-atom.level(x.i).E_cm)for x in atom.lines])
    bands={name:int(np.count_nonzero((waves>=lo)&(waves<hi)))for name,lo,hi in [('UV_below400nm',0,400),('VIS400_700nm',400,700),('NIR700_1000nm',700,1000),('IR1000_2000nm',1000,2000),('IR_at_least2000nm',2000,float('inf'))]}
    profile_points=sum(2*(max(x.nq,5)-1)+1 for x in atom.lines);matrix_bytes=len(atom.levels)**2*56*8
    conversion=dict(status='not attempted')
    try:
        model=to_lightweaver(atom);grids=[];supports=[]
        for line in model.lines:
            line.quadrature.Nlambda=2*(line.quadrature.Nlambda-1)+1
            wave=line.quadrature.wavelength(line);grids.append(wave);supports.append((wave[0],wave[-1]))
        grids.extend(np.asarray(c.wavelengthGrid)for c in model.continua);grids.append(np.array([500.,845.]))
        merged=np.unique(np.concatenate(grids))
        if not np.isfinite(merged).all()or np.any(merged<=0):raise ValueError('Invalid dry Fe wavelength support')
        local_nodes=sum(int(np.searchsorted(merged,hi,side='right')-np.searchsorted(merged,lo))for lo,hi in supports)
        conversion=dict(status='AtomicModel construction and dry quadrature only; not physics qualification',merged_candidate_points=len(merged),fits_current_HOP_packet_scope=bool(len(merged)<=MAX_PACKET_WAVELENGTHS),candidate_BB_profile_samples=local_nodes,projected_BB_phi_bytes_at56depth5rays2directions=local_nodes*56*5*2*8,scope='No Context, profile array, atmosphere, collision matrix, EOS packet or RT/SE solve allocated. Candidate union uses inherited O diagnostic numerical prescription; Fe grid selection still unqualified. Byte projections are arithmetic, not measured peak memory or runtime.')
    except ValueError as error:conversion=dict(status='conversion/dry-grid refused',reason=str(error))
    temperature=np.concatenate([x.temperature if x.temperature is not None else atom.temperature for x in atom.collisions])
    legacy_path=codex_path('grids.gerber_ts')
    report=dict(generator=fingerprint(__file__),inputs=checks,engine_receipt=engine,source_atom=checks[str(source)],source_abundance_input=atom.abundance,source_weight=atom.weight,levels=len(atom.levels),levels_per_ion=dict(by_stage),BB=len(atom.lines),BB_stage_pairs={f'{i}->{j}':n for (i,j),n in BB.items()},BF=len(atom.continua),BF_stage_pairs={f'{i}->{j}':n for (i,j),n in BF.items()},collision_records=len(atom.collisions),collision_kinds=dict(Counter(x.kind for x in atom.collisions)),source_temperature_span_K=[float(np.min(temperature)),float(np.max(temperature))],nonzero_Stark_coefficients=sum(x.gStark!=0 for x in atom.lines),bound_bound_centre_span_vac_nm=[float(np.min(waves)),float(np.max(waves))],centre_counts_by_inventory_bin=bands,source_q_counts=dict(Counter(x.nq for x in atom.lines)),raw_doubled_line_quadrature_nodes=profile_points,raw_profile_bytes_projection=profile_points*56*5*2*8,one_dense607_square_matrix_at56depth_bytes=matrix_bytes,dry_conversion=conversion,legacy_identity=dict(record=str(legacy),recorded_sha256=old['sha256'],byte_identical_to_supplied_atom=checks[str(source)]['sha256']==old['sha256'],registered_directory=str(legacy_path),registered_directory_present=legacy_path.is_dir(),scope='Different byte hashes do not prove changed physical data; complete comparison to the legacy atom is still owed. Missing registered directory is explicit, with no silent path substitute.'),gates=['Oxygen continuum numerical convergence and reviewed engine foundation','Fe atom identity/provenance/rights and assigned-line coverage, with FeII BB absence explicit','Fe-specific collision interpolation/direction and full atomic/background continuum ownership qualification; O41-specific bindings cannot be reused as-is','Dry resource gate before any full Fe column; no size/tolerance bypass','Matched-input FeI1D reference/departure/rate benchmark at declared source abundance, then current solar observable/line-by-line validation','FeII ionization diagnostic needs a separate eligible model or explicit LTE limitation','Parameter sensitivity and justified broader-band expansion; no tuning to a3D7.46 abundance target'],computational_completed=True,scientific_admission='HOLD',scope='Next-element readiness only. Source FeI/FeII/FeIII bookkeeping is not general FeII NLTE capability. Existing production Fe values/dispositions unaffected. No restricted/super-level atom, collision approximation, new gf, benchmark result, solar abundance or timing claim. Fe is the next expansion target per Ryan; C/N remain queued. Existing1137/1139/1372/1368 ownership and later corrections must be reconciled; no duplicate ticket or outreach.')
    for path,fp in checks.items():
        if fingerprint(path)!=fp:raise ValueError('Changed Fe inventory source')
    a.output.mkdir(parents=True,exist_ok=False);(a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({k:report[k]for k in ['levels','BB','BF','BB_stage_pairs','dry_conversion','legacy_identity']}))


if __name__=='__main__':main()
