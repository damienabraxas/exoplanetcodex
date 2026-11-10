"""Fixed-support centre prediction accounting, not a noise or opacity model."""
import numpy as np


def fixed_support_blocks(offsets,ratio0,ratio1,weights,bin_width_A):
    x,a,b,h=[np.asarray(v,float) for v in (offsets,ratio0,ratio1,weights)]
    if x.ndim!=1 or len(x)<2 or any(v.shape!=x.shape for v in (a,b,h)):
        raise ValueError('Wrong fixed-support arrays')
    if not all(np.isfinite(v).all() for v in (x,a,b,h)) or np.any(np.diff(x)<=0):
        raise ValueError('Nonfinite or unordered support')
    if not np.isfinite(bin_width_A) or bin_width_A<=0:raise ValueError('Invalid declared bin width')
    contributions=h*(b-a)
    groups=np.floor(x/bin_width_A).astype(int)
    blocks=[]
    for group in np.unique(groups):
        mask=groups==group
        blocks.append(dict(bin_index=int(group),lo_offset_A=float(group*bin_width_A),
            hi_offset_A=float((group+1)*bin_width_A),count=int(mask.sum()),
            centre_delta=float(np.sum(contributions[mask])),
            absolute_contribution_sum=float(np.sum(abs(contributions[mask]))),
            weight_sum=float(np.sum(h[mask])),weight_l1=float(np.sum(abs(h[mask])))))
    return dict(centre_delta=float(np.sum(contributions)),
        absolute_contribution_sum=float(np.sum(abs(contributions))),blocks=blocks,
        closure=float(sum(r['centre_delta'] for r in blocks)-np.sum(contributions)),
        scope='Fixed-support linear centre accounting only; signed contributions can cancel. Bins are diagnostic bookkeeping, not identified blends or independent errors.')
