"""Guarded input/units boundary for an isolated official calcrhox diagnostic."""
from pathlib import Path
import subprocess
import numpy as np


def read_rosseland_opacity(path, structure):
    lines = Path(path).read_text().splitlines()
    index = next(i for i, line in enumerate(lines) if 'KappaRoss' in line and 'Density' in line)
    n = structure['ndep']
    a = np.asarray([[float(x) for x in line.split()] for line in lines[index+1:index+1+n]])
    if a.shape != (n, 8) or not np.isfinite(a).all() or not np.array_equal(a[:, 0], np.arange(1, n+1)) or not np.array_equal(a[:, 1], structure['logtau_ross']) or np.any(a[:, 2] <= 0):
        raise ValueError('Malformed/mismatched native Rosseland opacity table')
    return a[:, 2]*.1  # cm2/g -> m2/kg


def column_input(logtau, opacity_m2_kg):
    tau = np.asarray(logtau, dtype=float)
    kappa = np.asarray(opacity_m2_kg, dtype=float)
    if tau.ndim != 1 or not 4 <= len(tau) <= 400 or kappa.shape != tau.shape or not np.isfinite(tau).all() or not np.isfinite(kappa).all() or np.any(np.diff(tau) <= 0) or np.any(kappa <= 0):
        raise ValueError('Require 4..400 ordered tau rows and positive finite opacity')
    # Preserve the original routine's log10 cgs opacity and default REAL inputs.
    with np.errstate(over='ignore', under='ignore', invalid='ignore'):
        t32 = tau.astype(np.float32)
        k32 = np.log10(kappa/.1).astype(np.float32)
        physical_tau = np.power(np.float32(10), t32)
        physical_opacity = np.power(np.float32(10), k32)
    if not np.isfinite(physical_tau).all() or not np.isfinite(physical_opacity).all() or np.any(physical_tau <= 0) or np.any(physical_opacity <= 0) or np.any(np.diff(physical_tau) <= 0):
        raise ValueError('Inputs collapse or overflow in native REAL/log coordinates')
    return str(len(tau))+'\n'+''.join(f'{t:.17e} {k:.17e}\n' for t, k in zip(tau, np.log10(kappa/.1)))


def run_column_kernel(binary, logtau, opacity_m2_kg):
    text = column_input(logtau, opacity_m2_kg)
    run = subprocess.run([str(Path(binary).resolve())], input=text, capture_output=True, text=True, timeout=10)
    if run.returncode:
        raise RuntimeError('Official column kernel failed: '+run.stderr)
    output = np.asarray([float(x.replace('D', 'E')) for x in run.stdout.split()])
    if output.shape != np.asarray(logtau).shape or not np.isfinite(output).all() or np.any(output <= 0) or np.any(np.diff(output) <= 0):
        raise ValueError('Invalid official column output')
    return output*10.  # g/cm2 -> kg/m2


def log_opacity_integral(logtau, opacity_m2_kg, top_column_kg_m2):
    """Exact integral of a declared piecewise power-law opacity, not native truth."""
    lt = np.asarray(logtau, dtype=float)
    k = np.asarray(opacity_m2_kg, dtype=float)
    if lt.ndim != 1 or len(lt) < 2 or k.shape != lt.shape or not np.isfinite(lt).all() or not np.isfinite(k).all() or np.any(np.diff(lt) <= 0) or np.any(k <= 0) or not np.isfinite(top_column_kg_m2) or top_column_kg_m2 <= 0:
        raise ValueError('Invalid log-opacity integral inputs')
    with np.errstate(over='ignore', under='ignore', invalid='ignore'):
        tau = 10.**lt
        width = np.diff(lt)*np.log(10.)
        q = width-np.diff(np.log(k))
        relative_exponential = np.ones_like(q)
        np.divide(np.expm1(q), q, out=relative_exponential, where=q != 0)
        increments = tau[:-1]/k[:-1]*width*relative_exponential
        result = top_column_kg_m2+np.r_[0., np.cumsum(increments)]
    if not np.isfinite(result).all() or np.any(np.diff(result) <= 0):
        raise ValueError('Unrepresentable log-opacity integral')
    return result
