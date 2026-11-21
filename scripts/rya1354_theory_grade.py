"""CANDIDATE theoretical accuracy class for computed E1 rates, following Li, Jönsson, Amarsi,
Li & Grumer 2023 (O I, arXiv 2304.02310) Sect. 2.2, Table 4-5 and the 'fourth procedure':

  dT  = |A_B - A_C| / max(A_B, A_C)                       (Eq. 3; Froese Fischer 2009, Ekman+2014)
  CF  = (|sum_jk c_j <j||T||k> c'_k| / sum_jk |c_j <j||T||k> c'_k|)^2   (Eq. 4; Cowan 1981),
        evaluated in the Babushkin gauge (Table 5 note), summed per CSF PAIR (j,k)
  dT~ = max(dT, <dT> of the transition's A-magnitude group)  (six groups, Sect. 3.2)
  class from (dT~, CF) by Table 5; A/B/C/D/E collect NIST {AA,A+,A}, {B+,B}, {C+,C}, {D+,D}, E.

dT is a STATISTICAL indicator for groups of lines (Li 2023, Li 2023 N I, Papoulia 2019), not a
per-line calibrated sigma. Status: HOLD candidate; thresholds and adoption need Ryan's ruling
(RYA-1352 THEORETICAL_CODEX governance). Nothing here writes a grade into a product.
"""
import math
from collections import defaultdict

LI2023_GROUP_EDGES=(1e-2,1e1,1e3,5e5,3.5e6)   # s^-1, six groups (Li 2023 Sect. 3.2, third procedure)
STATUS='CANDIDATE_HOLD_NOT_ADOPTED'


def dT(a_b,a_c):
    if a_b<0 or a_c<0:raise ValueError('Rates must be non-negative')
    m=max(a_b,a_c)
    return 0.0 if m==0 else abs(a_b-a_c)/m


def group_of(a,edges=LI2023_GROUP_EDGES):
    return sum(a>=e for e in edges)


def dT_tilde(rows,edges=LI2023_GROUP_EDGES):
    """rows: dicts with A_B, A_C. Returns per-row dict(dT, group, group_mean_dT, dT_tilde).
    Group by A_B (Babushkin, Li's preferred form) magnitude."""
    out=[dict(dT=dT(r['A_B'],r['A_C']),group=group_of(r['A_B'],edges)) for r in rows]
    by=defaultdict(list)
    for o in out:by[o['group']].append(o['dT'])
    mean={g:math.fsum(v)/len(v) for g,v in by.items()}
    for o in out:o.update(group_mean_dT=mean[o['group']],group_size=len(by[o['group']]),dT_tilde=max(o['dT'],mean[o['group']]))
    return out


def accuracy_class(d,cf):
    """Li 2023 Table 5. d = dT or dT~ as a FRACTION; cf = Babushkin CF."""
    if not (d>=0 and cf>=0):raise ValueError('dT and CF must be non-negative numbers')
    p=100*d
    band=0 if p<=3 else 1 if p<=10 else 2 if p<=25 else 3 if p<=50 else 4
    if cf<0.1:band=min(band+1,4)
    return 'ABCDE'[band]


def cancellation_factor(ctx):
    """Li Eq. 4 from a v44 FAC component context: terms summed per (lower CSF, upper CSF) pair
    before taking magnitudes. Kept (uncut) terms only, as in the native amplitude."""
    pairs=defaultdict(list)
    for b in ctx['B']:
        if not b['flags']:pairs[(b['lower'],b['upper'])].append(b['term'])
    sums=[math.fsum(v) for v in pairs.values()]
    total=abs(math.fsum(sums));den=math.fsum(abs(s) for s in sums)
    return 0.0 if den==0 else (total/den)**2
