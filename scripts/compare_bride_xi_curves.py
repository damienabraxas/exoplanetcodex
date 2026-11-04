"""Combine compatible declared xi slices and check coarse-curve reconstruction."""
import argparse,json,sys
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))


def check_shared_inputs(reports):
    """Reject different physics/thermal fields even when output wavelengths agree."""
    import copy
    keys=('generator','stark_source','stark_kernel','pressure_source','pressure_provider','provider','line_coefficients','lines','partition','source_manifest','source_populations','source_binding','wide_baseline','sweep_adapter')
    anchor=reports[0]
    for r in reports:
        if r['schema']!='codex.bride_xi_flux/1' or r['mode']!='cowley' or r['population_mode']!='coupled':raise ValueError('Wrong xi result schema/physics mode')
        for key in keys:
            if r[key]!=anchor[key]:raise ValueError('Incompatible xi physics: '+key)
        for key in ('source_spline','source_mn','controls'):
            if r['coupled'][key]!=anchor['coupled'][key]:raise ValueError('Incompatible coupled control: '+key)
    fields=[]
    for r in reports:
        b=copy.deepcopy(r['local_binding']['local_binding']);b['scalars'].pop('vturb_m_s');fields.append(b)
    if any(f!=fields[0]for f in fields):raise ValueError('Incompatible fixed atmosphere/composition')


def interpolate_profile(x,line,continuum,target,method):
    """Interpolate positive line/continuum flux separately; no extrapolation."""
    x=np.asarray(x);line=np.asarray(line);continuum=np.asarray(continuum)
    if len(x)<2 or np.any(np.diff(x)<=0)or not x[0]<=target<=x[-1]:raise ValueError('Invalid profile grid/domain')
    if line.ndim!=2 or line.shape!=continuum.shape or line.shape[0]!=len(x):raise ValueError('Invalid profile shape')
    if not all(np.isfinite(v).all()and np.all(v>0)for v in (line,continuum)):raise ValueError('Invalid physical flux')
    if method=='pchip':
        numerator=PchipInterpolator(x,line,axis=0,extrapolate=False)(target)
        denominator=PchipInterpolator(x,continuum,axis=0,extrapolate=False)(target)
    elif method=='linear':
        numerator=np.array([np.interp(target,x,v)for v in line.T])
        denominator=np.array([np.interp(target,x,v)for v in continuum.T])
    else:raise ValueError('Unknown curve operator')
    return numerator/denominator


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--local-plan',type=Path,required=True)
    p.add_argument('--wide-plan',type=Path,required=True)
    p.add_argument('--local-root',type=Path,required=True)
    p.add_argument('--wide-root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    from pipeline.bride_reference import fingerprint
    from pipeline.bride_observables import window_integral
    groups=[];all_reports=[];receipts=[]
    for plan_path,folder in ((args.local_plan,args.local_root),(args.wide_plan,args.wide_root)):
        plan=json.loads(plan_path.read_text());group=[]
        for case in plan['cases']:
            d=folder/case['name'];r=json.loads((d/'report.json').read_text())
            if r['plan']!=fingerprint(plan_path) or r['case']!=case or r['flux']!=fingerprint(d/'flux.csv'):raise ValueError('Changed xi curve input')
            f=np.genfromtxt(d/'flux.csv',delimiter=',',names=True)
            areas=np.array([1000*window_integral(f['wavelength_air_A'],f['normalized_flux'],c-.5,c+.5)for c in (7771.944,7774.166,7775.388)])
            group.append((case['xi'],areas,f));all_reports.append(r);receipts.append(fingerprint(d/'report.json'))
        groups.append(group)
    check_shared_inputs(all_reports)
    nominal=[v for g in groups for v in g if v[0]==1.]
    if len(nominal)!=2:raise ValueError('Both slices require nominal counterpart')
    np.testing.assert_array_equal(nominal[0][2],nominal[1][2])
    unique={}
    for group in groups:
        for xi,areas,data in group:
            if xi in unique:np.testing.assert_array_equal(areas,unique[xi])
            unique[xi]=areas
    wide=sorted(groups[1],key=lambda v:v[0]);wx=np.array([v[0]for v in wide]);wy=np.array([v[1]for v in wide])
    for v in wide:np.testing.assert_array_equal(v[2]['wavelength_air_A'],wide[0][2]['wavelength_air_A'])
    line_cube=np.array([v[2]['Fnu_W_m2_Hz']for v in wide]);continuum_cube=np.array([v[2]['continuum_Fnu_W_m2_Hz']for v in wide])
    diagnostic=[]
    for xi,area,data in groups[0]:
        if xi in wx:continue
        if not wx[0]<xi<wx[-1]:raise ValueError('No curve extrapolation permitted')
        linear=np.array([np.interp(xi,wx,wy[:,j])for j in range(3)])
        cubic=PchipInterpolator(wx,wy,axis=0,extrapolate=False)(xi)
        wave=data['wavelength_air_A'];np.testing.assert_array_equal(wave,wide[0][2]['wavelength_air_A'])
        profile_checks={}
        for method in ('linear','pchip'):
            prediction=interpolate_profile(wx,line_cube,continuum_cube,xi,method)
            error=prediction-data['normalized_flux']
            profile_checks[method]=dict(max_abs_normalized_flux_error=float(max(abs(error))),
                oxygen_aperture_max_error=[float(max(abs(error[abs(wave-c)<=.5])))for c in (7771.944,7774.166,7775.388)])
        diagnostic.append(dict(xi_km_s=xi,linear_minus_actual_mA=(linear-area).tolist(),pchip_minus_actual_mA=(cubic-area).tolist(),
            profile_reconstruction=profile_checks,profile_method='Interpolate emergent Fnu and continuum separately, then divide; no SE or interpolated-population admission'))
    args.output.mkdir(parents=True,exist_ok=False)
    path=args.output/'curve.csv';xs=np.array(sorted(unique));ys=np.array([unique[x]for x in xs]);delta=ys-nominal[0][1]
    np.savetxt(path,np.column_stack((xs,ys,delta)),delimiter=',',header='xi_km_s,area7771_mA,area7774_mA,area7775_mA,delta7771_mA,delta7774_mA,delta7775_mA',comments='')
    import matplotlib
    matplotlib.use('Agg');import matplotlib.pyplot as plt
    fig,axes=plt.subplots(3,1,figsize=(9,9),sharex=True)
    for j,(ax,c)in enumerate(zip(axes,(7771.944,7774.166,7775.388))):
        ax.plot(xs,delta[:,j],'o-',label='All fresh coupled solves');ax.axhline(0,color='black',lw=.5)
        ax.set_ylabel('Area change (mA)');ax.set_title(f'{c:.3f} AIR A');ax.grid(alpha=.2)
    axes[-1].set_xlabel('Microturbulence (km/s)');fig.suptitle('Combined solar xi response; fixed thermal structure and O8.70\nCoarse-curve checks are interpolation diagnostics, not physical accuracy bounds')
    fig.tight_layout();figure=args.output/'combined_xi_response.png';fig.savefig(figure,dpi=150);plt.close(fig)
    report=dict(schema='codex.bride_xi_curve_review/1',scientific_admission='HOLD',generator=fingerprint(__file__),inputs=receipts,
        nominal_flux_array_identity=True,unique_xi_km_s=xs.tolist(),curve=fingerprint(path),figure=fingerprint(figure),coarse_curve_checks=diagnostic,
        caveats='Known local points0.9/1.1 check reconstruction from wider0.5/0.75/1/1.25/1.5 slice. Not blind validation, new stellar reference, physical interpolation accuracy or an abundance fit.')
    (args.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps(diagnostic))

if __name__=='__main__':main()
