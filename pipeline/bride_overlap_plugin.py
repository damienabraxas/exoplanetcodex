"""Load an explicitly receipted process-local numerical experiment plugin."""
import ctypes,json
from pathlib import Path


def load_overlap_plugin(directory):
    from lightweaver import LwCompiled
    from pipeline.bride_reference import fingerprint
    directory=Path(directory);receipt=json.loads((directory/'receipt.json').read_text())
    if receipt['backend']!=fingerprint(LwCompiled.__file__):raise ValueError('Plugin backend identity changed')
    path=directory/'bride_overlap_guard.dylib'
    if receipt['plugin']!=fingerprint(path)or path.stat().st_size>5*1024**2:raise ValueError('Plugin identity or size guard failed')
    ctypes.CDLL(LwCompiled.__file__,mode=ctypes.RTLD_GLOBAL)
    name=receipt.get('scheme','bride_overlap_guard_scalar')
    if name not in ['bride_overlap_guard_scalar','bride_transition_local_scalar','bride_full_coupling_scalar','bride_full_coupling_cap_scalar']:raise ValueError('Unknown iteration plugin scheme')
    if name not in LwCompiled.FsIterationSchemes.names:LwCompiled.FsIterationSchemes.load_fns_from_path(str(path))
    properties=LwCompiled.FsIterationSchemes.scheme_properties(name)
    if not all(properties[k]for k in ['respectsFormalSolver','defaultPerAtomStorage','defaultWlaGijStorage']):raise ValueError('Unexpected plugin storage contract')
    return name,dict(receipt=fingerprint(directory/'receipt.json'),plugin=receipt['plugin'],backend=receipt['backend'],properties=properties)
