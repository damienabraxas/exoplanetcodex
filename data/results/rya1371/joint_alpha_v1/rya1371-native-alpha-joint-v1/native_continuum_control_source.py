"""Opt-in source-Hermite alpha for fixed-population opacity/rate diagnostics."""
import hashlib,json,os,subprocess,time
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from pipeline.multi_continuum import MultiExplicitContinuum


@dataclass(eq=False,repr=False)
class NativeInterpolatedContinuum(MultiExplicitContinuum):
    control: object = None

    def alpha(self,wavelength):
        if self.control is None:raise ValueError('Missing explicit native continuum control')
        return self.control.evaluate((self.i,self.j),wavelength)


class NativeContinuumControl:
    def __init__(self,tables,executable,identities,output,deadline):
        self.tables={key:np.array(value,copy=True)for key,value in tables.items()}
        for table in self.tables.values():table.flags.writeable=False
        self.executable=Path(executable);self.identities=identities;self.output=Path(output);self.deadline=deadline
        self.cache={};self.calls=[];self.compiled_checks=[]
        self.audit=dict(method='native-frequency-hermite-raw-source',table_kind='raw-source',identities=identities,calls=self.calls,compiled_alpha_checks=self.compiled_checks,
            scope='One shared compiled continuum alpha supplies active opacity/emissivity and rates. Exact source-table knots use their canonical Angstrom coordinates; arbitrary points are not snapped. Native float32 interpolation output converted cm2->m2. Original model table/support/grid, BB and held-LW thermal constants unchanged. Fixed-population diagnostic only; not coupled SE or physical adoption.')

    def verify_inputs(self):
        from pipeline.reference_engines import fingerprint
        for path,expected in self.identities.items():
            if fingerprint(path)!=expected:raise ValueError('Native alpha input changed: '+path)
        for row in self.calls:
            for key in ['input','output','stderr']:
                if fingerprint(self.output/row[key]['name'])!=row[key]:raise ValueError('Native alpha runtime evidence changed')

    def evaluate(self,pair,wavelength):
        from pipeline.reference_engines import fingerprint
        from pipeline.hop_packet import MAX_PACKET_WAVELENGTHS
        wave=np.asarray(wavelength,dtype=float)
        if pair not in self.tables or wave.ndim!=1 or len(wave)>MAX_PACKET_WAVELENGTHS or not np.isfinite(wave).all()or np.any(wave<=0):raise ValueError('Invalid native alpha query/support')
        key=(pair,hashlib.sha256(wave.tobytes()).hexdigest())
        if key in self.cache:return self.cache[key].copy()
        table=self.tables[pair];source_nm=table[:,0]*.1;inside=(wave>=source_nm[0])&(wave<=source_nm[-1]);result=np.zeros_like(wave)
        if np.any(inside):
            query=wave[inside]*10;local=wave[inside];index=np.searchsorted(source_nm,local);hit=index<len(source_nm);hit[hit]&=source_nm[index[hit]]==local[hit]
            query[hit]=table[index[hit],0]
            if len(query)==1:query=np.repeat(query,2)
            if time.perf_counter()>=self.deadline:raise TimeoutError('Original native alpha diagnostic deadline exceeded')
            if os.environ.get('OMP_NUM_THREADS')!='1'or os.environ.get('OPENBLAS_NUM_THREADS')!='1':raise ValueError('Require single-thread native alpha evaluation')
            self.output.mkdir(parents=True,exist_ok=True);name=f'alpha_{len(self.calls):04d}_{pair[0]+1}_{pair[1]+1}'
            request=f'{len(table)} {len(query)}\n'+'\n'.join(f'{x:.17e} {y:.17e}'for x,y in table)+'\n'+'\n'.join(f'{x:.17e}'for x in query)+'\n'
            input_file=self.output/(name+'.in');input_file.write_text(request)
            proc=subprocess.run([str(self.executable.resolve())],input=request,text=True,capture_output=True,timeout=self.deadline-time.perf_counter(),check=False)
            out=self.output/(name+'.out');err=self.output/(name+'.err');out.write_text(proc.stdout);err.write_text(proc.stderr)
            if proc.returncode:raise RuntimeError('Native continuum oracle refused; raw output retained')
            values=np.fromstring(proc.stdout,sep=' ')
            if values.shape!=query.shape or not np.isfinite(values).all()or np.any(values<0):raise ValueError('Invalid native continuum oracle output')
            result[inside]=values[:np.count_nonzero(inside)]*1e-4
            self.calls.append(dict(lower_1based=pair[0]+1,upper_1based=pair[1]+1,input_points=len(wave),inside_points=int(inside.sum()),canonical_source_knots=int(hit.sum()),input=fingerprint(input_file),output=fingerprint(out),stderr=fingerprint(err)))
        self.cache[key]=result.copy();return result

    def apply(self,model):
        if (len(model.levels),len(model.lines),len(model.continua))!=(41,256,22):raise ValueError('Native alpha control requires complete O41 source')
        if {(t.i,t.j)for t in model.continua}!=set(self.tables):raise ValueError('Native continuum graph mismatch')
        for pair,table in self.tables.items():
            values=self.evaluate(pair,table[:,0]*.1)
            np.testing.assert_array_equal(values,table[:,1].astype(np.float32).astype(float)*1e-4)
        self.audit['source_tables_reproduced_at_knots']=len(self.tables)
        changed=[]
        for t in model.continua:
            table=self.tables[t.i,t.j]
            np.testing.assert_array_equal(np.asarray(t.wavelengthGrid),table[:,0]*.1)
            np.testing.assert_array_equal(np.asarray(t.alphaGrid),table[:,1]*1e-4)
            replacement=NativeInterpolatedContinuum(i=t.i,j=t.j,wavelengthGrid=t.wavelengthGrid,alphaGrid=t.alphaGrid,control=self)
            replacement.setup(model);changed.append(replacement)
        model.continua=changed
        return self.audit

    def verify_compiled(self,atom):
        count=0
        for t in atom.trans:
            if t.type=='Line':continue
            actual=np.asarray(t.alpha);expected=self.evaluate((t.i,t.j),np.asarray(t.wavelength));np.testing.assert_array_equal(actual,expected);count+=1
            self.compiled_checks.append(dict(lower_1based=t.i+1,upper_1based=t.j+1,points=len(actual),alpha_sha256=hashlib.sha256(actual.tobytes()).hexdigest()))
        if count!=22:raise ValueError('Incomplete compiled continuum control')


def load_native_continuum_control(spec,atom_path,atom,*,resolver,output,deadline):
    from pipeline.reference_engines import ROOT,fingerprint,load_engine
    if set(spec)!={'method','resource','receipt','atom'}or spec['method']!='native-frequency-hermite-raw-source':raise ValueError('Unknown native continuum interpolation control')
    if fingerprint(atom_path)!=spec['atom']or(len(atom.levels),len(atom.lines),len(atom.continua))!=(41,256,22):raise ValueError('Changed/incomplete native alpha source atom')
    directory=Path(resolver(spec['resource']));receipt=directory/'receipt.json';doc=json.loads(receipt.read_text())
    if fingerprint(receipt)!=spec['receipt']or not doc['completed']or len(doc.get('analytic_controls',[]))!=2:raise ValueError('Unqualified native continuum interpolation oracle')
    exe=directory/'continuum_probe';checks={str(receipt):fingerprint(receipt),str(exe):fingerprint(exe),str(atom_path):fingerprint(atom_path),str(Path(__file__)):fingerprint(__file__)}
    if fingerprint(exe)!=doc['executable']:raise ValueError('Changed native continuum executable')
    checks.update(doc['inputs'])
    checks[str(ROOT/'pipeline/multi_continuum.py')]=fingerprint(ROOT/'pipeline/multi_continuum.py')
    for name in ['LwTransition.hpp','SimdFullIterationTemplates.hpp','Dispatch_compute_full_operator_rates.ipp']:
        path=Path(resolver('engines.bride_eval'))/'lw_src/Source'/name;checks[str(path)]=fingerprint(path)
    _,app,_,_=load_engine('openmulti3d')
    for name in ['trans_mod.f90','interpol_mod.f90','const_mod.f90']:
        original=app/'HOP'/name;linked=[Path(p)for p in doc['inputs']if p.endswith('/HOP/'+name)]
        if len(linked)!=1 or fingerprint(linked[0])!=fingerprint(original):raise ValueError('Oracle differs from native continuum source')
        checks[str(original)]=fingerprint(original)
    tables={(t.i-1,t.j-1):np.column_stack([t.wavelength[::-1],t.alpha[::-1]])for t in atom.continua}
    control=NativeContinuumControl(tables,exe,checks,output,deadline);control.verify_inputs()
    return control,checks
