"""Exact trace and explicitly labelled diagnostics of the held continuum rule."""
import numpy as np


def selection(wave,flux,centre,reference):
    from pipeline import local_continuum as lc
    w,f,m=[np.asarray(x,float) for x in (wave,flux,reference)]
    if w.ndim!=1 or f.shape!=w.shape or m.shape!=w.shape or np.any(np.diff(w)<=0):raise ValueError('Wrong continuum grid')
    eligible=np.isfinite(w)&np.isfinite(f)&np.isfinite(m)&(m>0)&(abs(w-centre)<=lc.ENV_HALF_WIDTH_A)&(abs(w-centre)>.5)
    if eligible.sum()<lc.MIN_MODEL_PIX:raise ValueError('Insufficient eligible reference pixels')
    threshold=float(np.percentile(m[eligible],lc.MODEL_Q))
    return eligible&(m>=threshold),threshold


def support_geometry(wave,centre,indices):
    x=np.asarray(wave)[indices]-centre
    design=np.column_stack((x,np.ones_like(x)))
    if len(x)<2 or np.linalg.matrix_rank(design)!=2:raise ValueError('Rank-deficient continuum support')
    weights=np.linalg.lstsq(design.T,np.array([0.,1.]),rcond=None)[0]
    return dict(count=len(x),blue_count=int(sum(x<0)),red_count=int(sum(x>0)),
        min_offset_A=float(x.min()),max_offset_A=float(x.max()),
        centre_outside_hull_A=float(max(x.min(),-x.max(),0.)),design_condition_2=float(np.linalg.cond(design)),
        centre_weight_l1=float(sum(abs(weights))),centre_weight_l2=float(np.linalg.norm(weights)),
        centre_weights=weights.tolist(),
        scope='Fixed-support linear prediction geometry only; not a noise model or physical error budget')


def trace_fit(wave,flux,centre,reference,mask,*,clip=True):
    from pipeline import local_continuum as lc
    w,f,m=[np.asarray(x,float) for x in (wave,flux,reference)];mask=np.asarray(mask,bool)
    if w.shape!=f.shape or m.shape!=w.shape or mask.shape!=w.shape:raise ValueError('Wrong trace arrays')
    indices=np.flatnonzero(mask)
    if len(indices)<lc.MIN_MODEL_PIX:
        return dict(applied=False,reason='INSUFFICIENT_ORIGINAL_ANCHORS',count=len(indices),
            coefficient_support_count=0,postfit_survivor_count=len(indices),
            clipping_candidate_count=None,clipping_pass_count=0,
            clipping_stop_reason='INSUFFICIENT_ORIGINAL_ANCHORS',clipping_limit_reached=False),None
    if not np.isfinite(f[indices]).all() or not np.isfinite(m[indices]).all() or np.any(m[indices]<=0):raise ValueError('Invalid anchor ratio')
    passes=[];stop_reason=None
    for number in range(5 if clip else 1):
        used=indices.copy();x=w[used]-centre;y=f[used]/m[used]
        geometry=support_geometry(w,centre,used)
        c=np.polyfit(x,y,lc.DEGREE);residual=y-np.polyval(c,x)
        mad=1.4826*float(np.median(abs(residual-np.median(residual))))
        keep=abs(residual)<=lc.CLIP_MAD_MODEL*max(mad,1e-4)
        passes.append(dict(pass_index=number+1,fit_indices=used.tolist(),clip_keep_indices=used[keep].tolist(),coefficients=c.tolist()))
        if not clip:
            stop_reason='CLIPPING_DISABLED';break
        if keep.all():
            stop_reason='ALL_RETAINED';break
        if keep.sum()<lc.MIN_MODEL_PIX:
            stop_reason='CLIP_WOULD_UNDERRUN_MINIMUM';break
        indices=used[keep]
    if stop_reason is None:stop_reason='ITERATION_LIMIT_AFTER_CLIP'
    level=float(c[-1]);applied=abs(level-1)<=lc.MAX_SHIFT
    trace=dict(applied=applied,reason='applied' if applied else 'CONTINUUM_UNPLACEABLE',
        level=level,slope_per_A=float(c[0]),reported_n_bins=len(indices) if applied else int(mask.sum()),fit_indices=used.tolist(),
        surviving_indices=indices.tolist(),passes=passes,geometry=geometry,
        postfit_survivors_differ=not np.array_equal(indices,used),
        coefficient_support_count=len(used),postfit_survivor_count=len(indices),
        clipping_candidate_count=int(keep.sum()),clipping_pass_count=len(passes),
        clipping_stop_reason=stop_reason,clipping_limit_reached=stop_reason=='ITERATION_LIMIT_AFTER_CLIP')
    # Independent fixed-support reconstruction of the centre prediction.
    estimate=float(np.dot(geometry['centre_weights'],y))
    bound=64*np.finfo(float).eps*geometry['design_condition_2']*max(1.,float(np.max(abs(y))))
    if abs(estimate-level)>bound:raise ValueError('Centre influence reconstruction failed')
    return trace,np.polyval(c,w-centre) if applied else None


def two_factor_changes(corners,key):
    """Symmetric finite path accounting; effects include induced clipping."""
    a,b,c,d=[corners[k][key] for k in ('synthesis_synthesis','synthesis_iag','iag_synthesis','iag_iag')]
    selector=.5*((c-a)+(d-b));denominator=.5*((b-a)+(d-c))
    return dict(total=d-a,selector=selector,denominator_plus_induced_clipping=denominator,
        interaction=d-c-b+a,closure=(selector+denominator)-(d-a),
        scope='Finite2x2 path-average accounting, not independent uncertainty terms or physical causal attribution')
