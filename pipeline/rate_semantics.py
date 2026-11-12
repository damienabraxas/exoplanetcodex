"""Prevent effective SE coefficients being treated as gross physical rates."""
import numpy as np

def gross_rate_relative_difference(reference,replay,*,reference_semantics,replay_semantics):
    if reference_semantics!='physical-gross'or replay_semantics!='physical-gross':
        raise ValueError('Gross-rate comparison requires fresh unpreconditioned physical rates on both sides')
    reference=np.asarray(reference);replay=np.asarray(replay)
    if reference.shape!=replay.shape or not np.isfinite(reference).all()or not np.isfinite(replay).all()or np.any(reference<0)or np.any(replay<0):raise ValueError('Invalid physical-rate arrays')
    if np.any((reference==0)&(replay!=0)):raise ValueError('Physical rate support differs')
    mask=reference>0
    return float(np.max(abs(replay[mask]/reference[mask]-1)))if np.any(mask)else 0.
