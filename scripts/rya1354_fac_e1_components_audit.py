"""Portable reader/audit for v44 signed native E1 component records.

Reconstructs FAC E1 reduced amplitudes from every logged term and checks them
against the native RawTransition amplitude. Exact identities are required where
the native arithmetic is identical (A-order sum, products, radial cache values,
cut flags). Where only summation ORDER differs (FAC merges CSF-pair terms after
an unspecified qsort order), the admitted difference is the standard recursive
summation bound gamma_n*sum|x| (Higham, Accuracy and Stability, eq. 4.4), with
u=2**-53; it is a floating-point bound, not a physical tolerance.
Diagnostics only: no rate, lifetime, uncertainty or grade is produced here.
"""
from collections import defaultdict
import math

U=2.0**-53
TAG='FAC_E1_COMPONENTS_V44'


def gamma(n):
    if n*U>=1:raise ValueError('Summation bound undefined')
    return n*U/(1-n*U)


def parse(text):
    """Split the append-only record file into per-call contexts, losslessly."""
    lines=text.splitlines()
    if not lines or lines[0].split('\t')[:2]!=['#',TAG]:raise ValueError('Component header differs')
    header=dict(x.split('=',1) for x in lines[0].split('\t')[2:])
    contexts=[];cur=None
    for no,line in enumerate(lines[1:],2):
        f=line.split('\t');tag=f[0]
        if tag=='C':
            if cur is not None and 'end' not in cur:raise ValueError(f'Unterminated context before line {no}')
            if len(f)!=15:raise ValueError(f'Context columns differ line {no}')
            low,up,j1,j2,p1,p2,nele,g,mode,ih1,ih2=map(int,f[1:12])
            cur=dict(lower=low,upper=up,two_j=(j1,j2),parity=(p1,p2),nele=nele,gauge=g,mode=mode,iham=(ih1,ih2),
                     energy_au=float(f[12]),aw=float(f[13]),swapped=int(f[14]),grid=[],levels=([],[]),norm={},A=[],B=[])
            contexts.append(cur);continue
        if cur is None or 'end' in cur:raise ValueError(f'Record outside a context line {no}')
        if tag=='G':
            if len(f)!=3 or int(f[1])!=len(cur['grid']):raise ValueError(f'Grid record differs line {no}')
            cur['grid'].append(float(f[2]))
        elif tag=='L':
            if len(f)!=13:raise ValueError(f'Level-basis columns differ line {no}')
            side,pos=int(f[1]),int(f[2])
            if pos!=len(cur['levels'][side]):raise ValueError(f'Basis order differs line {no}')
            cur['levels'][side].append(dict(ibasis=int(f[3]),state=int(f[4]),mixing=float(f[5]),kgroup=int(f[6]),
                kcfg=int(f[7]),kstate=int(f[8]),group=f[9],complex=f[10],config=f[11],name=f[12]))
        elif tag=='N':
            if len(f)!=5:raise ValueError(f'Norm columns differ line {no}')
            cur['norm'][int(f[1])]=dict(n_basis=int(f[2]),ham_nbasis=int(f[3]),sum_mix2=float(f[4]))
        elif tag=='A':
            if len(f)!=10 or int(f[1])!=len(cur['A']):raise ValueError(f'Angular record differs line {no}')
            cur['A'].append(dict(pair=((int(f[2]),int(f[3])),(int(f[4]),int(f[5]))),coeff=float(f[6]),radial=float(f[7]),
                term=float(f[8]),running=float(f[9])))
        elif tag=='B':
            if len(f)!=14:raise ValueError(f'CSF-term columns differ line {no}')
            cur['B'].append(dict(lower=int(f[1]),upper=int(f[2]),sub=int(f[3]),angular=float(f[4]),k=int(f[5]),
                pair=((int(f[6]),int(f[7])),(int(f[8]),int(f[9]))),coeff=float(f[10]),flags=int(f[11]),
                radial=float(f[12]),term=float(f[13])))
        elif tag=='E':
            if len(f)!=4 or int(f[1])!=cur['lower'] or int(f[2])!=len(cur['B']):raise ValueError(f'Context end differs line {no}')
            cur['end']=float(f[3])
        else:raise ValueError(f'Unknown record line {no}')
    if not contexts or 'end' not in contexts[-1]:raise ValueError('Missing or unterminated context')
    return header,contexts


def verify(ctx,native_amplitude,angzcut):
    """Exact native identities plus a derived summation bound; raises on any breach."""
    if ctx['mode']!=0 or ctx['two_j'][0]+ctx['two_j'][1]<2 or ctx['parity'][0]==ctx['parity'][1]:
        raise ValueError('Not a full-relativistic E1 context')
    s=0.0;seen=set()
    for a in ctx['A']:
        if a['pair'] in seen:raise ValueError('Duplicate merged orbital pair')
        seen.add(a['pair'])
        if a['term']!=a['radial']*a['coeff']:raise ValueError('Native A product differs')
        s+=a['term']
        if s!=a['running']:raise ValueError('Native A running sum differs')
    if s!=ctx['end'] or s!=native_amplitude:raise ValueError('A-order amplitude does not equal native amplitude')
    lev=ctx['levels']
    for side in (0,1):
        n=ctx['norm'][side]
        if n['n_basis']!=len(lev[side]) or n['n_basis']>n['ham_nbasis']:raise ValueError('Stored basis count differs')
        acc=0.0
        for b in lev[side]:acc+=b['mixing']*b['mixing']
        if acc!=n['sum_mix2']:raise ValueError('Stored mixing norm differs')
    radial={a['pair']:a['radial'] for a in ctx['A']}
    merged=defaultdict(list)
    for b in ctx['B']:
        m1=lev[0][b['lower']]['mixing'];m2=lev[1][b['upper']]['mixing']
        first,second=(m2,m1) if ctx['swapped'] else (m1,m2)
        prod=first*second;r0=b['angular']*prod
        if abs(b['coeff'])!=abs(r0):raise ValueError('CSF term is not angular x mixing product')
        flags=(abs(first)<angzcut)|((abs(second)<angzcut)<<1)|((abs(prod)<angzcut)<<2)|((abs(r0)<angzcut)<<3)
        if flags!=b['flags']:raise ValueError('Native angular-cut flags differ')
        if b['term']!=b['radial']*b['coeff']:raise ValueError('CSF term product differs')
        if b['k']!=2:raise ValueError('Non-E1 rank in E1 context')
        if not flags:
            if b['pair'] not in radial or radial[b['pair']]!=b['radial']:raise ValueError('Kept term radial differs from native cache')
            merged[b['pair']].append(b['coeff'])
    if set(merged)!=set(radial):raise ValueError('Kept CSF pairs do not span native orbital pairs')
    worst=0.0
    for a in ctx['A']:
        x=merged[a['pair']];bound=gamma(len(x))*math.fsum(abs(v) for v in x)
        d=abs(math.fsum(x)-a['coeff'])
        if d>bound:raise ValueError('Merged angular coefficient exceeds summation bound')
        worst=max(worst,d/bound if bound else 0.0)
    kept=[b['term'] for b in ctx['B'] if not b['flags']]
    total=math.fsum(kept);abs_total=math.fsum(abs(t) for t in kept)
    bound=2*gamma(len(kept)+1)*abs_total
    if abs(total-s)>bound:raise ValueError('Kept CSF reconstruction exceeds summation bound')
    return dict(native_amplitude=s,kept_terms=len(kept),cut_terms=len(ctx['B'])-len(kept),
                csf_reconstruction=total,reconstruction_difference=total-s,reconstruction_bound=bound,
                worst_coefficient_fraction_of_bound=worst)


def grid_context(ctx):
    g=ctx['grid'];aw=ctx['aw']
    exact=[i for i,x in enumerate(g) if x==aw]
    near=min(range(len(g)),key=lambda i:abs(g[i]-aw))
    return dict(points=g,aw=aw,exact_grid_index=exact[0] if exact else None,nearest_index=near,
                nearest_relative_offset=(g[near]-aw)/aw,inside=g[0]<=aw<=g[-1])


def decompose(ctx,top=25):
    """Signed decomposition by CSF pair, configuration pair and orbital pair."""
    lev=ctx['levels'];kept=[b for b in ctx['B'] if not b['flags']];cut=[b for b in ctx['B'] if b['flags']]
    S=math.fsum(b['term'] for b in kept);A=math.fsum(abs(b['term']) for b in kept)
    csf=defaultdict(list);cfg=defaultdict(list);orb=defaultdict(list)
    for b in kept:
        lo,up=lev[0][b['lower']],lev[1][b['upper']]
        csf[(b['lower'],b['upper'])].append(b['term'])
        cfg[(lo['config'],up['config'])].append(b['term'])
        orb[b['pair']].append(b['term'])
    def rows(groups,label):
        out=[]
        for key,v in groups.items():
            t=math.fsum(v);out.append(dict(key=label(key),amplitude=t,abs_sum=math.fsum(abs(x) for x in v),
                                          fraction_of_total=t/S if S else None,terms=len(v)))
        out.sort(key=lambda r:-abs(r['amplitude']))
        return out
    csf_rows=rows(csf,lambda k:dict(lower=k[0],upper=k[1],lower_name=lev[0][k[0]]['name'],upper_name=lev[1][k[1]]['name'],
        lower_config=lev[0][k[0]]['config'],upper_config=lev[1][k[1]]['config'],
        lower_mixing=lev[0][k[0]]['mixing'],upper_mixing=lev[1][k[1]]['mixing']))
    pl=max(range(len(lev[0])),key=lambda i:abs(lev[0][i]['mixing']))
    pu=max(range(len(lev[1])),key=lambda i:abs(lev[1][i]['mixing']))
    principal=math.fsum(csf.get((pl,pu),[]))
    return dict(amplitude=S,abs_sum=A,cancellation_index=A/abs(S) if S else math.inf,
        csf_pairs=len(csf_rows),top_csf_pairs=csf_rows[:top],
        principal=dict(lower=pl,upper=pu,lower_weight=lev[0][pl]['mixing']**2,upper_weight=lev[1][pu]['mixing']**2,
                       amplitude=principal,fraction_of_total=principal/S if S else None),
        configuration_pairs=rows(cfg,lambda k:dict(lower_config=k[0],upper_config=k[1])),
        orbital_pairs=rows(orb,lambda k:dict(pair=[list(k[0]),list(k[1])])),
        cut=dict(terms=len(cut),amplitude=math.fsum(b['term'] for b in cut),abs_sum=math.fsum(abs(b['term']) for b in cut)),
        stored_basis={str(s):dict(ctx['norm'][s],dropped_by_mixing_cut=ctx['norm'][s]['ham_nbasis']-ctx['norm'][s]['n_basis'])
                      for s in (0,1)})


def configuration_shift(first,second):
    """Signed change of each configuration-pair amplitude between two contexts.
    Labels join by nonrelativistic configuration names; orbitals differ between
    models, so this is a labelled diagnostic, not an adiabatic overlap proof."""
    a={(r['key']['lower_config'],r['key']['upper_config']):r['amplitude'] for r in first['configuration_pairs']}
    b={(r['key']['lower_config'],r['key']['upper_config']):r['amplitude'] for r in second['configuration_pairs']}
    out=[dict(lower_config=k[0],upper_config=k[1],first=a.get(k,0.0),second=b.get(k,0.0),change=b.get(k,0.0)-a.get(k,0.0))
         for k in set(a)|set(b)]
    out.sort(key=lambda r:-abs(r['change']))
    return dict(total_change=second['amplitude']-first['amplitude'],
                accounted=math.fsum(r['change'] for r in out),rows=out)
