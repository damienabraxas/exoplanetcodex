"""Exact consumed-table strength precision control; never a fitted scale."""
import numpy as np


def strength_factors(lines,weights,path):
    text=path.read_text().splitlines();starts=[i+1 for i,x in enumerate(text)if x.split()[:2]==['i','j']and'lam_air'in x]
    if len(starts)!=1:raise ValueError('Ambiguous consumed strength table')
    factors=[];rows=[]
    for index,line in enumerate(lines):
        row=text[starts[0]+index].split()
        if len(row)!=13 or tuple(map(int,row[:3]))!=(index,line.i,line.j):raise ValueError('Consumed strength transition identity mismatch')
        loggf=float(row[5]);value=10**loggf/weights[line.i]/line.f
        if not np.isfinite(value)or value<=0:raise ValueError('Invalid source strength precision factor')
        factors.append(value);rows.append(dict(index=index,lower_1based=line.i+1,upper_1based=line.j+1,source_f=line.f,consumed_loggf=loggf,factor=value))
    return np.asarray(factors),rows
